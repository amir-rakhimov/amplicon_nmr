# Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank 
library(qiime2R)
library(phyloseq)
library(tidyverse)

asvlevel=F # specify if we're agglomerating at ASV level (not Species)
truncationlvl<-"234" # truncation level that we chose in QIIME2

authorname<-"pooled" # name of the folder with QIIME2 output
qiimedir<-paste0("./data/qiime/",authorname,"-qiime/") # directory with QZA files
metadatadir<-paste0("./data/metadata/",authorname,"-metadata/") # directory with metadata

# this is the taxonomic rank that will be used for agglomeration
if(asvlevel==TRUE){
  agglom.rank<-"OTU"
}else{
  agglom.rank<-"Genus"
}
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
source("./code/r-scripts/make_ps_pretty.R") # import a script to make phyloseq
# object pretty
# specify the name of your metadata file
metadata.filename<-paste0(metadatadir,
                          paste("filenames",read.end.type,
                                authorname,"raw-supercomp.tsv", sep = "-"))

# Import qza files and convert them into a phyloseq object ####
ps.q<-qza_to_phyloseq(
  features = paste0(qiimedir,authorname,"-",read.end.type,"-filtered-table-trimmed-dada2-",
                    truncationlvl,".qza"), # feature table
  taxonomy = paste0(qiimedir,authorname,"-",read.end.type,"-taxonomy-trimmed-dada2-",
                    truncationlvl,".qza"), # taxonomy
  tree = paste0(qiimedir,authorname,"-",read.end.type,"-rooted-tree-trimmed-dada2-",
                truncationlvl,".qza") # rooted tree
)

# add custom metadata cause previous command loses metadata for some reason
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
# convert the Sample column into row names
custom.md<-custom.md%>%column_to_rownames(var = "Sample") 
# assign the custom metadata as your phyloseq object's metadata
sample_data(ps.q)<-custom.md

# you can exclude some samples based on class
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
#                                               'ntccontrol','rabbitcontrol',
#                                               'harecontrol'),]
custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
# you can exclude samples based on their library size (total number of reads)
custom.md<-custom.md[!rownames(custom.md) %in%
                       intersect(names(which(colSums(ps.q@otu_table)<20000)),
                                 rownames(custom.md)),]

## Construct the phyloseq object directly from dada2 output ####
# we combine the phyloseq object with new metadata (if we excluded samples)
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
# Remove entries with zero Abundance
ps.q.agg<-ps.q.agg%>%
  filter(Abundance!=0)
# change the name d__Kingdom to Kingdom
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

# add relative abundance column: Abundance divided by total abundance in a sample
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

# Extract taxa with highest mean relative abundance (i.e. mean relative abundance >1%) ####
# We will get a dataset of four columns: class (animal host), 
# two taxonomic ranks (e.g Genus, Family)  and abundance 

# group the dataframe by classes (animal hosts)
classcol<-which(colnames(ps.q.agg) =="class")
# find a column by which we agglomerated the dataset
# for ASV, we actually use Species
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}

# for each class, we take a agglom.rank, sum its abundances from all samples,
# then take a mean. This will be our MeanRelativeAbundance
if(asvlevel==TRUE){
  ps.q.agg<-ps.q.agg%>%
    group_by(class,OTU,Genus,Family)%>% # group by class (animal host),
    # agglom.rank and the preceding column (based on index)
    # TODO: I still haven't figured out how to do it with indices only
    # then, compute MeanRelativeAbundance from RelativeAbundance 
    mutate(MeanRelativeAbundance = mean(RelativeAbundance)) 
}else{ 
  ps.q.agg<-ps.q.agg%>%
    group_by_at(c(classcol,agglom.rank.col,agglom.rank.col-1))%>% # group by class,
    # agglom.rank and the preceding column (based on index)
    # compute MeanRelativeAbundance from RelativeAbundance 
    mutate(MeanRelativeAbundance = mean(RelativeAbundance))
}

## get taxa with mean abundance >1% ####
# ps.q.agg[ps.q.agg$MeanRelativeAbundance>1,] will output TRUE/FALSE
# depending which MeanRelativeAbundance values are >1
# then, we keep only TRUE values and take columns "class", agglom.rank, and
# MeanRelativeAbundance
ps.q.1pc<-distinct(ps.q.agg[ps.q.agg$MeanRelativeAbundance>1,
                            c("class",agglom.rank,"MeanRelativeAbundance")])


# we need to decide whether a taxon in a specific host is present in 1%
# create a dummy column with "class agglom.rank" strings
if(agglom.rank=="OTU"){
  ps.q.agg<-ps.q.agg%>%
    ungroup()%>% # ungroup because the dataset was grouped previously
    mutate(class_agglom.rank=paste(class,OTU)) # in the original dataset
  ps.q.1pc<-ps.q.1pc%>%
    ungroup()%>% # ungroup because the dataset was grouped previously
    mutate(class_agglom.rank=paste(class,OTU)) # and in the 1% dataset
}else{
  ps.q.agg<-ps.q.agg%>%
    ungroup()%>%# ungroup because the dataset was grouped previously
    mutate(class_agglom.rank=paste(class,Genus)) # in the original dataset
  ps.q.1pc<-ps.q.1pc%>%
    ungroup()%>% # ungroup because the dataset was grouped previously
    mutate(class_agglom.rank=paste(class,Genus)) # and in the 1% dataset
}

## If otu has MeanRelativeAbundance <1%, set it as Remainder ####
# Taxon.bp is for barplot
ps.q.agg$Taxon.bp<-ps.q.agg$Taxon
ps.q.agg$Taxon.bp<-ifelse(ps.q.agg$class_agglom.rank %in% ps.q.1pc$class_agglom.rank,
                          ps.q.agg$Taxon.bp,
                     "Remainder (Mean abundance < 1%)")

# if our Taxon is in the 1pc dataset (mean abundance >1%), keep it as it is
# otherwise, set it as Remainder (Mean abundance < 1%)

# Save the workspace
save.image(paste0("./output/rdafiles/",paste(authorname,read.end.type,"qiime2",
                                      truncationlvl,agglom.rank,
                                      "phyloseq-workspace.RData",sep = "-")))
