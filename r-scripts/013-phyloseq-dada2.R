library(phyloseq)
library(ggplot2)
library(tidyverse)


truncationlvl<-"284-203"
# 313-229
# 284-229
# 284-203
# 273-203
# 273-198
# 265-203
# 265-198
# trunclvls<-c("313-229",
#              "284-229",
#              "284-203",
#              "273-203",
#              "273-198",
#              "265-203",
#              "265-198")

authorname<-"alldir"
directory<-paste0("./data/",authorname,"-data/")
agglom.rank<-"Genus"

# Open dada2 output and convert into phyloseq ####
load("./rdafiles/yashliu-r-dada2-284-203.RData")

# Phyloseq ####
samples.out <- as.data.frame(rownames(seqtab.nochim))
colnames(samples.out)<-"Sample"
# Add custom metadata
custom.md<-read.table(paste0("./data/alldir-data/","filenames-yashliu-final-supercomp.tsv"),
                      header = T)
colnames(custom.md)[1]<-"Sample"
custom.md$Sample<-gsub("mf_","mf-",custom.md$Sample)
custom.md$Sample<-paste(custom.md$Sample,"final",sep="_")

# Merge 
custom.md<-inner_join(custom.md,samples.out,
                      by="Sample")
# Exclude classes
custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx'),]

custom.md<-custom.md%>%column_to_rownames(var = "Sample")

## Construct the phyloseq object directly from dada2 output ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(custom.md),
               tax_table(taxa))

# Create a dataframe of relative abundances ####
## Agglomerate by agglom.rank, transform into rel.ab, then ####
# create a dataframe
ps.agg.rel<-ps %>%
  #subset_samples(class=="NMR") %>% # choose only naked mole-rat
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="Bacteria")%>% # choose only bacteria
  transform_sample_counts(function(x)100* x / sum(x)) %>% # transform into percentages
  psmelt()  # transform the phyloseq object into an R dataframe
head(ps.agg.rel) 

# Total counts for each sample ####
ps.total<-ps %>%
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="Bacteria")%>% # choose only bacteria
  psmelt() %>% # transform the phyloseq object into an R dataframe
  group_by(Sample) %>% # group by Sample id to perform operations on samples
  summarise(Abundance=sum(Abundance)) # sum abundances by sample id

# Substitute NA values
ps.agg.rel[is.na(ps.agg.rel)]<-"NonNA"

## Create a vector with taxonomic ranks ####
phylum<-which(colnames(ps.agg.rel) =="Phylum") # identify the position of phylum column
agglom.rank.col<-which(colnames(ps.agg.rel) ==agglom.rank) # identify the 
# position of agglom.rank column
taxcols<-phylum:agglom.rank.col # vector of taxonomic ranks

# Replace empty taxa with Unclassified+previous taxon ####
# we go through each taxonomic rank from phylum to agglom.rank
# we substitute "uncultured" or "unclassified" and make it like 
# Unclassified (Previous taxonomic rank)
# foo<-ps.agg.rel
# ps.agg.rel<-foo
for (i in 1:nrow(ps.agg.rel)){
  for (j in taxcols){
    if(ps.agg.rel[i,j]=="NonNA"){
      ps.agg.rel[i,j]<-paste0("Unclassified"," (",ps.agg.rel[i,j-1],")") 
    }
    if (grepl("Unclassified",ps.agg.rel[i, j-1])){
      ps.agg.rel[i,j]<-ps.agg.rel[i,j-1]
    }
    if (ps.agg.rel[i,j]=="uncultured"){
      ps.agg.rel[i,j]<-paste0("Uncultured"," (",ps.agg.rel[i,j-1],")")
    }
    if (grepl("Uncultured",ps.agg.rel[i, j-1])){
      ps.agg.rel[i,j]<-ps.agg.rel[i,j-1]
    }
  }
}

# sanity check: is total relative abundance of each sample 100%?
ps.agg.rel %>%
  group_by(Sample) %>% # Group by sample id
  mutate(Abundance = as.numeric(Abundance))%>% # transform Abundance column into numeric
  summarise(Abundance = sum(Abundance)) %>% # Sum all abundances
  pull(Abundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value

# Extract taxa with highest mean relative abundance ####
# dataset of three columns: two taxonomic ranks (e.g Genus, Family) 
# and abundance 
# group the dataframe by classes (animal classes)
classcol<-which(colnames(ps.agg.rel) =="class")
# sanity check
ps.agg.rel%>%
  group_by_at(c(classcol,agglom.rank.col,agglom.rank.col-1))%>% # group by class,
  # agglom.rank and the preceding column (based on index)
  mutate(Abundance = as.numeric(Abundance))%>% # convert Abundance column into numeric
  summarise(Abundance = mean(Abundance)) %>% # convert Abundance into mean abundance
  ungroup()%>%
  group_by(class)%>%
  summarise(sss=sum(Abundance)) # get total relative abundance for each class

# the actual pooling
ps.agg.rel.pooled<-ps.agg.rel%>%
  group_by_at(c(classcol,agglom.rank.col,agglom.rank.col-1))%>% # group by class,
  # agglom.rank and the preceding column (based on index)
  mutate(Abundance = as.numeric(Abundance))%>% # convert Abundance column into numeric
  summarise(Abundance = mean(Abundance)) %>% # convert Abundance into mean abundance
  group_by(class,Genus,Family)%>%
  arrange(class,desc(Abundance)) # order descending

# get otus with mean abundance >1% ####
ps.1pc<-ps.agg.rel.pooled[ps.agg.rel.pooled$Abundance>1,]


# Create a Taxon column where items are either Unclassified (previous rank) ####
# or agglom.rank (previous taxonomic rank)
# e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
ps.agg.rel$Taxon<-
  ifelse(grepl("Unclassified|Uncultured",
               ps.agg.rel[[agglom.rank.col]]),
         ps.agg.rel[[agglom.rank.col]], 
         paste0(ps.agg.rel[[agglom.rank.col]]," (",ps.agg.rel[[agglom.rank.col-1]],")")) 
# grepl finds which agglom.rank column items are uncultured or unclassified
# if true, keeps Unclassified (previous rank)
# if false (so, it's not Unclassified, actually has taxonomic classification), 
# converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)


# If otu has mean rel.ab <1%, set it as Remainder ####
ps.agg.rel$Taxon<-ifelse(ps.agg.rel[[agglom.rank.col]] %in% ps.1pc[[agglom.rank]],
                           ps.agg.rel$Taxon,
                           "Remainder (Mean abundance < 1%)")
# if our Taxon is in the 1pc dataset (mean abundance >1%), keep it as it is
# otherwise, set it as Remainder (Mean abundance < 1%)


# Create a dataframe of absolute abundances ####
# remove taxa with 0 abundance
ps.nonzero<-prune_taxa(taxa_sums(ps)>0,ps)

## extract absolute abundances ####
ps.agg.abs<-ps.nonzero %>%
  #subset_samples(class=="NMR") %>% # choose only naked mole-rat
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="Bacteria")%>% # choose only bacteria
  psmelt()  # transform the phyloseq object into an R dataframe

ps.agg.abs[is.na(ps.agg.abs)]<-"NonNA"

## Replace empty taxa with Unclassified+previous taxon ####
# find the column with your rank of interest
# then extract columns between phylum and your rank
phylum<-which(colnames(ps.agg.abs) =="Phylum")
agglom.rank.col<-which(colnames(ps.agg.abs) ==agglom.rank)
taxcols<-phylum:agglom.rank.col

for (i in 1:nrow(ps.agg.abs)){
  for (j in taxcols){
    if(ps.agg.abs[i,j]=="NonNA"){
      ps.agg.abs[i,j]<-paste0("Unclassified"," (",ps.agg.abs[i,j-1],")") 
    }
    if (grepl("Unclassified",ps.agg.abs[i, j-1])){
      ps.agg.abs[i,j]<-ps.agg.abs[i,j-1]
    }
    if (ps.agg.abs[i,j]=="uncultured"){
      ps.agg.abs[i,j]<-paste0("Uncultured"," (",ps.agg.abs[i,j-1],")")
    }
    if (grepl("Uncultured",ps.agg.abs[i, j-1])){
      ps.agg.abs[i,j]<-ps.agg.abs[i,j-1]
    }
  }
}

## add a column with merged taxonomic ranks ####
ps.agg.abs$taxa.full<-
  apply(ps.agg.abs[,c((phylum-1):agglom.rank.col)],
        1,paste,collapse="_")

save.image(paste0("./rdafiles/yashliu-dada2-",truncationlvl,"-",agglom.rank,
                  "-138-1-phyloseq-workspace.RData"))
