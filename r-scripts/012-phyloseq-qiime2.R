library(qiime2R)
library(phyloseq)
library(ggplot2)
library(tidyverse)
# library(pals)
# library(patchwork)
# library(Polychrome) 

truncationlvl<-"234"
# 274-229

authorname<-"alldir"
directory<-paste0("./data/",authorname,"-data/")
qiimedir<-paste0(directory,"qiime/")
agglom.rank<-"Genus"
read.end.type<-"single"
source("./r-scripts/make_ps_pretty.R")

# Import qza files and convert them into a phyloseq object ####
ps.q<-qza_to_phyloseq(
  # metadata = paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"),
  features = paste0(qiimedir,"pooled-",read.end.type,"-filtered-table-trimmed-dada2-",
                    truncationlvl,".qza"),
  taxonomy = paste0(qiimedir,"pooled-",read.end.type,"-taxonomy-trimmed-dada2-",
                    truncationlvl,".qza"),
  tree = paste0(qiimedir,"pooled-",read.end.type,"-rooted-tree-trimmed-dada2-",
                truncationlvl,".qza")
)

# add custom metadata cause previous command loses metadata for some reason
# custom.md<-read.table(paste0("./data/alldir-data/","filenames-pooled-final-supercomp.tsv"),
#                       header = T)
custom.md<-read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
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

# Create a dataframe of relative abundances ####
## Agglomerate by agglom.rank, transform into rel.ab, then ####
# create a dataframe
ps.q.agg.rel<-ps.q %>%
  #subset_samples(class=="NMR") %>% # choose only naked mole-rat
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  transform_sample_counts(function(x)100* x / sum(x)) %>% # transform into percentages
  psmelt()#%>%  # transform the phyloseq object into an R dataframe
  #filter(Abundance!=0)

head(ps.q.agg.rel) 
# ps.q.agg.rel$Sample<-paste(ps.q.agg.rel$Sample,"final",sep="_")

## rename d__Bacteria ####
ps.q.agg.rel$Kingdom<-
  gsub("d__","",ps.q.agg.rel$Kingdom)

# Total counts for each sample ####
ps.q.total<-ps.q %>%
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  psmelt() %>% # transform the phyloseq object into an R dataframe
  group_by(Sample) %>% # group by Sample id to perform operations on samples
  summarise(Abundance=sum(Abundance)) # sum abundances by sample id
# ps.q.total$Sample<-paste(ps.q.total$Sample,"final",sep="_")

# Replace empty taxa with Unclassified+previous taxon ####
ps.q.agg.rel.pretty<-make_ps_pretty(ps.q.agg.rel,"Genus")
ps.q.agg.rel<-ps.q.agg.rel.pretty
rm(ps.q.agg.rel.pretty)

# sanity check: is total relative abundance of each sample 100%?
ps.q.agg.rel %>%
  group_by(Sample) %>% # Group by sample id
  mutate(Abundance = as.numeric(Abundance))%>% # transform Abundance column into numeric
  summarise(Abundance = sum(Abundance)) %>% # Sum all abundances
  pull(Abundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value

# Extract taxa with highest mean relative abundance ####
# dataset of four columns: class, two taxonomic ranks (e.g Genus, Family) 
# and abundance 

# group the dataframe by classes
classcol<-which(colnames(ps.q.agg.rel) =="class")
agglom.rank.col<-which(colnames(ps.q.agg.rel) ==agglom.rank)
ps.q.agg.rel%>%
  group_by_at(c(classcol,agglom.rank.col,agglom.rank.col-1))%>% 
  mutate(Abundance = as.numeric(Abundance))%>% 
  summarise(Abundance = mean(Abundance)) %>%
  ungroup()%>%
  group_by(class)%>%
  summarise(sss=sum(Abundance))

# for each class, we take a genus, sum its abundances from all samples,
# then take a mean. This will be our new Abundance
ps.q.agg.rel.pooled<-ps.q.agg.rel%>%
  group_by_at(c(classcol,agglom.rank.col,agglom.rank.col-1))%>% # group by class,
  # agglom.rank and the preceding column (based on index)
  mutate(Abundance = as.numeric(Abundance))%>% # convert Abundance column into numeric
  summarise(Abundance = mean(Abundance)) %>% # convert Abundance into mean abundance
  group_by(class,Genus,Family)%>%
  arrange(class,desc(Abundance))%>% # order descending
  filter(Abundance!=0)

# get otus with mean abundance >1% ####
ps.q.1pc<-ps.q.agg.rel.pooled[ps.q.agg.rel.pooled$Abundance>1,]

# If otu has mean rel.ab <1%, set it as Remainder ####
# Taxon.bp is for barplot
ps.q.agg.rel$Taxon.bp<-ps.q.agg.rel$Taxon
ps.q.agg.rel$Taxon.bp<-ifelse(ps.q.agg.rel[[agglom.rank.col]] %in% ps.q.1pc[[agglom.rank]],
                           ps.q.agg.rel$Taxon.bp,
                            "Remainder (Mean abundance < 1%)")
# if our Taxon is in the 1pc dataset (mean abundance >1%), keep it as it is
# otherwise, set it as Remainder (Mean abundance < 1%)


# Create a dataframe of absolute abundances ####
# remove taxa with 0 abundance
# ps.q.nonzero<-prune_taxa(taxa_sums(ps.q)>0,ps.q)

## extract absolute abundances ####
# ps.q.agg.abs<-ps.q.nonzero %>%
ps.q.agg.abs<-ps.q %>%
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  psmelt()#%>%  # transform the phyloseq object into an R dataframe
  # filter(Abundance!=0)

ps.q.agg.abs$Kingdom<-
  gsub("d__","",ps.q.agg.abs$Kingdom)
## Replace empty taxa with Unclassified+previous taxon ####
ps.q.agg.abs.pretty<-make_ps_pretty(ps.q.agg.abs,"Genus")
ps.q.agg.abs<-ps.q.agg.abs.pretty
rm(ps.q.agg.abs.pretty)

save.image(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
                  "-phyloseq-workspace.RData"))
