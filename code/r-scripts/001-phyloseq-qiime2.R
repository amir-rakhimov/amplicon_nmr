# Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank ####

# Once you produce the feature table, taxonomic classification, and the 
# phylogenetic tree in QIIME2, it's time to perform downstream processing
# in R. First, we need to import the QZA files using `qiime2R` package.
# We convert the QZA files directly into phyloseq objects.

# The final output of this script is an ASV table with X columns:  
# `Sample`: 
#TODO: add the rest
# `Sample`                
# `Abundance`             
# `class`                 
# `animal`                
# `sex`                   
# `birthday`              
# `Kingdom`               
# `Phylum`                
# `Class`                 
# `Order`
# ...
# Genus (if agglomerated by Genus)
# `RelativeAbundance`
# `MeanRelativeAbundance`


## Import libraries ####
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(microViz)
## Specifying parameters and directory/file names #### 
# This is the taxonomic rank that will be used for agglomeration
agglom.rank<-"Phylum"
# specify if we're agglomerating at ASV level (not Species)
if(agglom.rank=="OTU"){
  asvlevel<-TRUE
}else{
  asvlevel<-FALSE
}
truncationlvl<-"234" # truncation level that we chose in QIIME2
authorname<-"pooled" # name of the folder with QIIME2 output
date_time<-"20240425_02_57_13"
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime"),
                    paste(date_time,read.end.type,truncationlvl,sep="-")) # directory with QZA files

metadatadir<-file.path("./data/metadata",
                       paste(authorname,"metadata",sep = "-")) # directory with metadata

# Specify the name of your metadata file
metadata.filename<-file.path(metadatadir,
                          paste("filenames",read.end.type,
                                authorname,"raw-supercomp.tsv", sep = "-"))

## Import qza files and convert them into a phyloseq object ####
ps.q<-qza_to_phyloseq(
  features = file.path(qiimedir, paste0(paste(authorname,read.end.type,
                                              "filtered-table-trimmed-dada2",
                    truncationlvl,sep="-"),".qza")), # feature table
  taxonomy = file.path(qiimedir,paste0(paste(authorname,read.end.type,
                                             "filtered-taxonomy-trimmed-dada2",
                    truncationlvl,sep="-"),".qza")), # taxonomy
  tree = file.path(qiimedir,paste0(paste(authorname,read.end.type,
                                         "rooted-tree-trimmed-dada2",
                truncationlvl,sep="-"),".qza")) # rooted tree
)
# Change the name d__Kingdom to Kingdom
ps.q.taxtab<-as.data.frame(tax_table(ps.q))
ps.q.taxtab$Kingdom<-
  gsub("d__","",ps.q.taxtab$Kingdom)
tax_table(ps.q)<-as.matrix(ps.q.taxtab)
rm(ps.q.taxtab)

# Add custom metadata cause previous command loses metadata for some reason
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
# convert the Sample column into row names because phyloseq needs samples
# as rownames
# Remove absolute.filepath column
custom.md<-custom.md%>%
  column_to_rownames(var = "Sample")%>%
  select(-absolute.filepath)
custom.md$class<-as.factor(custom.md$class)
custom.md$animal<-as.factor(custom.md$animal)
custom.md$sex<-as.factor(custom.md$sex)
custom.md$birthday<-as.Date(custom.md$birthday)
# assign the custom metadata as your phyloseq object's metadata
sample_data(ps.q)<-custom.md

# you can exclude some samples based on class. Specify the excluded classes
# in a vector, then use the `%in%` operator. It will remove entries 
# of the `class` column (animal hosts) from the `custom.md` object (metadata)
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
#                                               'ntccontrol','rabbitcontrol',
#                                               'harecontrol'),]
custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
# you can exclude samples based on their library size (total number of reads)
custom.md<-custom.md[!rownames(custom.md) %in%
                       intersect(names(which(colSums(ps.q@otu_table)<20000)),
                                 rownames(custom.md)),]

### Construct the phyloseq object directly from dada2 output ####
# We combine the phyloseq object with new metadata (if we excluded samples)
ps.foo <- phyloseq(otu_table(ps.q),
                   sample_data(custom.md),
                   tax_table(ps.q),
                   phy_tree(ps.q))
ps.q<-ps.foo
rm(ps.foo)

# Number of features in the unfiltered dataset
length(rownames(ps.q@tax_table@.Data))

# Total frequency in the unfiltered dataset
sum(colSums(ps.q@otu_table@.Data))

# Summary statistics (min, median, max, quartiles) of the unfiltered dataset
ps.q@otu_table@.Data%>%
  colSums()%>%
  summary()

# Select only Bacteria
ps.q<-ps.q %>%
  subset_taxa(Kingdom%in%"Bacteria")

### Fix empty taxa with higher rank taxon ####
# Because we want to remove NA values and make ambiguous "uncultured" or 
# "unclassified" taxa more understandable.
ps.q<-tax_fix(ps.q,unknowns = c("NA","uncultured","Unassigned",
                                "uncultured_bacterium","uncultured_rumen",
                                "gut_metagenome","human_gut","mouse_gut",
                                "wallaby_gut","uncultured_soil", 
                                "uncultured_organism","uncultured_prokaryote"))
## Create a dataframe of absolute abundances ####
### Extract absolute abundances ####
if (asvlevel==TRUE){
  ps.q.agg<-ps.q %>%
    subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
    psmelt()  # transform the phyloseq object into an R dataframe
}else{
  ps.q.agg<-ps.q %>%
    tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
    subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
    psmelt()  # transform the phyloseq object into an R dataframe
}

# Remove entries with zero Abundance
ps.q.agg<-ps.q.agg%>%
  filter(Abundance!=0)

# Number of samples in the filtered dataset ####
ps.q.agg%>%
  distinct(Sample)%>%
  tally()

# Number of features in the filtered dataset ####
ps.q.agg%>%
  distinct(OTU)%>%
  tally()

# Total frequency in the filtered dataset ####
ps.q.agg%>%
  summarise(TotalAbundance=sum(Abundance))

# Summary statistics (min, median, max, quartiles) of the filtered dataset ####
ps.q.agg%>%
  select(Sample,Abundance)%>%
  group_by(Sample)%>%
  summarise(FrequencyPerSample=sum(Abundance))%>%
  select(FrequencyPerSample)%>%
  summary()

# Add relative abundance column: Abundance divided by total abundance in a sample ####
ps.q.agg<-ps.q.agg%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample",agglom.rank))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)

# Sanity check: is total relative abundance of each sample 100%?
ps.q.agg %>%
  group_by(Sample) %>% # Group by sample id
  summarise(sumRelativeAbundance = sum(RelativeAbundance)) %>% # Sum all abundances
  pull(sumRelativeAbundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value

ps.q.agg %>%
  mutate(sumRelativeAbundance = sum(RelativeAbundance)) %>%
  ggplot(aes(x=Sample,y=sumRelativeAbundance))+
  geom_bar(stat="identity")
## Add mean relative abundance data ####
# We will group the dataset by three columns: class (animal host), 
# two taxonomic ranks (e.g Genus, Family), and maybe OTU (actually ASV)
# if we agglomerate at ASV level.

# Group the dataframe by classes (animal hosts)
classcol<-which(colnames(ps.q.agg) =="class")
# find a column by which we agglomerated the dataset
# for ASV, we actually use Species
# We find the column index because we can derive the previous taxonomic 
# rank with a simple agglom.rank.col-1.
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}

# First, we calculate the library size per sample.
# Then, inside each class, we take a agglom.rank, sum its abundances from all samples,
# then take a mean. This will be our MeanRelativeAbundance
if(asvlevel==TRUE){
  ps.q.agg<-ps.q.agg%>%
    group_by(class)%>% # group by class (animal host),
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",agglom.rank))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
}else{ 
  ps.q.agg<-ps.q.agg%>%
    group_by(class)%>% # group by class,
    # compute MeanRelativeAbundance from Abundance 
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",agglom.rank))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
    select(-OTU)
}
ps.q.agg<-ps.q.agg%>%
  select(-TotalClass,-TotalSample,-TotalAgglomRank)

objects.to.keep<-c("agglom.rank","ps.q.agg","asvlevel","custom.md",
                   "authorname","truncationlvl","read.end.type")
objects.to.keep<-which(ls()%in%objects.to.keep)
rm(list = ls()[-objects.to.keep])
# Save the workspace
save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

sessionInfo()
