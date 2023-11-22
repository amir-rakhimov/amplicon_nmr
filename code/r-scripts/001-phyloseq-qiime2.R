# Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank 
library(qiime2R)
library(phyloseq)
library(tidyverse)

asvlevel=F
truncationlvl<-"234"

authorname<-"pooled"
directory<-paste0("./data/",authorname,"-data/")
qiimedir<-paste0(directory,"qiime/")

if(asvlevel==TRUE){
  agglom.rank<-"OTU"
}else{
  agglom.rank<-"Genus"
}
read.end.type<-"single"
source("./r-scripts/make_ps_pretty.R")

# Import qza files and convert them into a phyloseq object ####
ps.q<-qza_to_phyloseq(
  # metadata = paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"),
  features = paste0(qiimedir,authorname,"-",read.end.type,"-filtered-table-trimmed-dada2-",
                    truncationlvl,".qza"),
  taxonomy = paste0(qiimedir,authorname,"-",read.end.type,"-taxonomy-trimmed-dada2-",
                    truncationlvl,".qza"),
  tree = paste0(qiimedir,authorname,"-",read.end.type,"-rooted-tree-trimmed-dada2-",
                truncationlvl,".qza")
)

# add custom metadata cause previous command loses metadata for some reason
# custom.md<-read.table(paste0("./data/pooled-data/","filenames-pooled-final-supercomp.tsv"),
#                       header = T)
custom.md<-read.table(paste0(directory,"filenames-single-",authorname,"-raw-supercomp.tsv"),
                      header = T)
colnames(custom.md)[1]<-"Sample"
custom.md<-custom.md%>%column_to_rownames(var = "Sample")

sample_data(ps.q)<-custom.md

# you can exclude some samples based on class
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
#                                               'ntccontrol','rabbitcontrol',
#                                               'harecontrol'),]
custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
custom.md<-custom.md[!rownames(custom.md) %in%
                       intersect(names(which(colSums(ps.q@otu_table)<20000)),
                                 rownames(custom.md)),]


## Construct the phyloseq object directly from dada2 output ####
# we combine qza with new metadata
ps.foo <- phyloseq(otu_table(ps.q),
                   sample_data(custom.md),
                   tax_table(ps.q),
                   phy_tree(ps.q))
ps.q<-ps.foo
rm(ps.foo)

# Create a dataframe of absolute abundances ####
## extract absolute abundances ####
if (asvlevel==TRUE){
  ps.q.agg<-ps.q %>%
    subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
    psmelt()  # transform the phyloseq object into an R dataframe
}else{
  ps.q.agg<-ps.q %>%
    tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
    subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
    psmelt()  # transform the phyloseq object into an R dataframe
}
ps.q.agg<-ps.q.agg%>%
  filter(Abundance!=0)

ps.q.agg$Kingdom<-
  gsub("d__","",ps.q.agg$Kingdom)
## Replace empty taxa with Unclassified+previous taxon ####
if (asvlevel==TRUE){
  ps.q.agg.pretty<-make_ps_pretty(ps.q.agg,"OTU")
}else{
  ps.q.agg.pretty<-make_ps_pretty(ps.q.agg,agglom.rank)
}
ps.q.agg<-ps.q.agg.pretty
rm(ps.q.agg.pretty)

# add relative abundance
ps.q.agg<-ps.q.agg%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)

# sanity check: is total relative abundance of each sample 100%?
ps.q.agg %>%
  group_by(Sample) %>% # Group by sample id
  summarise(RelativeAbundance = sum(RelativeAbundance)) %>% # Sum all abundances
  pull(RelativeAbundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value

# Total counts for each sample ####
ps.q.total<-ps.q.agg%>%
  group_by(Sample)%>%
  summarise(TotalAbundance=sum(Abundance))

# Extract taxa with highest mean relative abundance ####
# dataset of four columns: class, two taxonomic ranks (e.g Genus, Family) 
# and abundance 

# group the dataframe by classes
classcol<-which(colnames(ps.q.agg) =="class")
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}


# for each class, we take a genus, sum its abundances from all samples,
# then take a mean. This will be our MeanRelativeAbundance
if(asvlevel==TRUE){
  ps.q.agg<-ps.q.agg%>%
    group_by(class,OTU,Genus,Family)%>% # group by class,
    # agglom.rank and the preceding column (based on index)
    mutate(MeanRelativeAbundance = mean(RelativeAbundance)) # convert RelativeAbundance into MeanRelativeAbundance
}else{ 
  ps.q.agg<-ps.q.agg%>%
    group_by_at(c(classcol,agglom.rank.col,agglom.rank.col-1))%>% # group by class,
    # agglom.rank and the preceding column (based on index)
    mutate(MeanRelativeAbundance = mean(RelativeAbundance))# convert RelativeAbundance into MeanRelativeAbundance
}

# get otus with mean abundance >1% ####
ps.q.1pc<-distinct(ps.q.agg[ps.q.agg$MeanRelativeAbundance>1,
                            c("class",agglom.rank,"MeanRelativeAbundance")])


# we need to decide whether a taxon in a specific host is present in 1%
if(agglom.rank=="OTU"){
  ps.q.agg<-ps.q.agg%>%
    ungroup()%>%
    mutate(class_agglom.rank=paste(class,OTU))
  ps.q.1pc<-ps.q.1pc%>%
    ungroup()%>%
    mutate(class_agglom.rank=paste(class,OTU))
}else{
  ps.q.agg<-ps.q.agg%>%
    ungroup()%>%
    mutate(class_agglom.rank=paste(class,Genus))
  ps.q.1pc<-ps.q.1pc%>%
    ungroup()%>%
    mutate(class_agglom.rank=paste(class,Genus))
}

# If otu has mean rel.ab <1%, set it as Remainder ####
# Taxon.bp is for barplot
ps.q.agg$Taxon.bp<-ps.q.agg$Taxon
ps.q.agg$Taxon.bp<-ifelse(ps.q.agg$class_agglom.rank %in% ps.q.1pc$class_agglom.rank,
                          ps.q.agg$Taxon.bp,
                     "Remainder (Mean abundance < 1%)")

# if our Taxon is in the 1pc dataset (mean abundance >1%), keep it as it is
# otherwise, set it as Remainder (Mean abundance < 1%)
save.image(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                      truncationlvl,agglom.rank,
                                      "phyloseq-workspace.RData",sep = "-")))
