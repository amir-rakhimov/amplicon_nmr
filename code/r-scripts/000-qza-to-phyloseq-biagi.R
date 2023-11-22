library(qiime2R)
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(pals)
library(patchwork)
library(Polychrome) 

# For Biagi dataset
authorname<-"biagi"
truncationlvl<-"no-trunc"
directory<-paste0("./data/",authorname,"-data/")
qiimedir<-paste0(directory,"qiime/")
agglom.rank<-"Genus"
agglom.rank<-"Family"

# Import qza files and convert them into a phyloseq object ####
ps.b<-qza_to_phyloseq(
  # metadata = paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"),
  features = paste0(qiimedir,authorname,"-table-trimmed-dada2-",
                    truncationlvl,".qza"),
  taxonomy = paste0(qiimedir,authorname,"-taxonomy-trimmed-dada2-",
                    truncationlvl,".qza"),
  tree = paste0(qiimedir,authorname,"-rooted-tree-trimmed-dada2-",
                truncationlvl,".qza")
)
# Add custom metadata
custom.md<-read.table(paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"),
                      header = T)
custom.md<-data.frame(colnames(ps.b@otu_table),
                      class="NMRwt")

colnames(custom.md)[1]<-"Sample"
custom.md<-custom.md%>%column_to_rownames(var = "Sample")

sample_data(ps.b)<-custom.md
# View(ps.b@sam_data)
# read_q2metadata(paste0(qiimedir,"filenames-",authorname,"-supercomp.tsv"))%>%View()

# Create a dataframe of relative abundances ####
## Agglomerate by agglom.rank, transform into rel.ab, then ####
# create a dataframe
ps.b.agg.rel<-ps.b %>%
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  transform_sample_counts(function(x)100* x / sum(x)) %>% # transform into percentages
  psmelt()  # transform the phyloseq object into an R dataframe
head(ps.b.agg.rel) 



## rename d__Bacteria ####
ps.b.agg.rel$Kingdom<-
  gsub("d__","",ps.b.agg.rel$Kingdom)

# Total counts for each sample ####
ps.b.total<-ps.b %>%
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  psmelt() %>% # transform the phyloseq object into an R dataframe
  group_by(Sample) %>% # group by Sample id to perform operations on samples
  summarise(Abundance=sum(Abundance)) # sum abundances by sample id

# Substitute NA values
ps.b.agg.rel[is.na(ps.b.agg.rel)]<-"NonNA"

## Create a vector with taxonomic ranks ####
phylum<-which(colnames(ps.b.agg.rel) =="Phylum") # identify the position of phylum column
agglom.rank.col<-which(colnames(ps.b.agg.rel) ==agglom.rank) # identify the 
# position of agglom.rank column
taxcols<-phylum:agglom.rank.col # vector of taxonomic ranks

# Replace empty taxa with Unclassified+previous taxon ####
# we go through each taxonomic rank from phylum to agglom.rank
# we substitute "uncultured" or "unclassified" and make it like 
# Unclassified (Previous taxonomic rank)
# foo<-ps.b.agg.rel
# ps.b.agg.rel<-foo
for (i in 1:nrow(ps.b.agg.rel)){
  for (j in taxcols){
    if(ps.b.agg.rel[i,j]=="NonNA"){
      ps.b.agg.rel[i,j]<-paste0("Unclassified"," (",ps.b.agg.rel[i,j-1],")") 
    }
    if (grepl("Unclassified",ps.b.agg.rel[i, j-1])){
      ps.b.agg.rel[i,j]<-ps.b.agg.rel[i,j-1]
    }
    if (ps.b.agg.rel[i,j]=="uncultured"){
      ps.b.agg.rel[i,j]<-paste0("Uncultured"," (",ps.b.agg.rel[i,j-1],")")
    }
    if (grepl("Uncultured",ps.b.agg.rel[i, j-1])){
      ps.b.agg.rel[i,j]<-ps.b.agg.rel[i,j-1]
    }
  }
}

# sanity check: is total relative abundance of each sample 100%?
ps.b.agg.rel %>%
  group_by(Sample) %>% # Group by sample id
  mutate(Abundance = as.numeric(Abundance))%>% # transform Abundance column into numeric
  summarise(Abundance = sum(Abundance)) %>% # Sum all abundances
  pull(Abundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value


ps.b.agg.rel.pooled<-ps.b.agg.rel %>%
  group_by_at(c(agglom.rank.col,agglom.rank.col-1)) %>% # group by agglom.rank 
  # and the preceding column (based on index)
  mutate(Abundance = as.numeric(Abundance))%>% # convert Abundance column into numeric
  summarise(Abundance = mean(Abundance)) %>% # convert Abundance into mean abundance
  arrange(-Abundance) # order descending
head(ps.b.agg.rel.pooled)

# get otus with mean abundance >1% ####
ps.b.1pc<-ps.b.agg.rel.pooled[ps.b.agg.rel.pooled$Abundance>1,]

# Create a Taxon column where items are either Unclassified (previous rank) ####
# or agglom.rank (previous taxonomic rank)
# e.g Unclassified (Lachnospiraceae) or Lactobacillus (Lactobacillaceae)
ps.b.agg.rel$Taxon<-ifelse(grepl("Unclassified|Uncultured",
                                ps.b.agg.rel[[agglom.rank.col]]),  
                          ps.b.agg.rel[[agglom.rank.col]], 
                          paste0(ps.b.agg.rel[[agglom.rank.col]]," (",ps.b.agg.rel[[agglom.rank.col-1]],")")) 
# grepl finds which agglom.rank column items are uncultured or unclassified
# if true, keeps Unclassified (previous rank)
# if false (so, it's not Unclassified, actually has taxonomic classification), 
# converts into agglom.rank taxon (previous rank taxon) :Lactobacillus (Lactobacillaceae)

# If otu has mean rel.ab <1%, set it as Remainder ####
ps.b.agg.rel$Taxon<-ifelse(ps.b.agg.rel[[agglom.rank.col]] %in% ps.b.1pc[[1]],
                          ps.b.agg.rel$Taxon,
                          "Remainder (Mean abundance < 1%)")
# if our Taxon is in the 1pc dataset (mean abundance >1%), keep it as it is
# otherwise, set it as Remainder (Mean abundance < 1%)


# Create a dataframe of absolute abundances ####
# remove taxa with 0 abundance
ps.b.nonzero<-prune_taxa(taxa_sums(ps.b)>0,ps.b)

## extract absolute abundances ####
ps.b.agg.abs<-ps.b.nonzero %>%
  tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")%>% # choose only bacteria
  psmelt()  # transform the phyloseq object into an R dataframe

ps.b.agg.abs$Kingdom<-
  gsub("d__","",ps.b.agg.abs$Kingdom)

ps.b.agg.abs[is.na(ps.b.agg.abs)]<-"NonNA"

## Replace empty taxa with Unclassified+previous taxon ####
# find the column with your rank of interest
# then extract columns between phylum and your rank
phylum<-which(colnames(ps.b.agg.abs) =="Phylum")
agglom.rank.col<-which(colnames(ps.b.agg.abs) ==agglom.rank)
taxcols<-phylum:agglom.rank.col

for (i in 1:nrow(ps.b.agg.abs)){
  for (j in taxcols){
    if(ps.b.agg.abs[i,j]=="NonNA"){
      ps.b.agg.abs[i,j]<-paste0("Unclassified"," (",ps.b.agg.abs[i,j-1],")") 
    }
    if (grepl("Unclassified",ps.b.agg.abs[i, j-1])){
      ps.b.agg.abs[i,j]<-ps.b.agg.abs[i,j-1]
    }
    if (ps.b.agg.abs[i,j]=="uncultured"){
      ps.b.agg.abs[i,j]<-paste0("Uncultured"," (",ps.b.agg.abs[i,j-1],")")
    }
    if (grepl("Uncultured",ps.b.agg.abs[i, j-1])){
      ps.b.agg.abs[i,j]<-ps.b.agg.abs[i,j-1]
    }
  }
}

## add a column with merged taxonomic ranks ####
ps.b.agg.abs$taxa.full<-
  apply(ps.b.agg.abs[,c((phylum-1):agglom.rank.col)],
        1,paste,collapse="_")

save.image(paste0("./rdafiles/biagi-",agglom.rank,"-","phyloseq-workspace.RData"))
