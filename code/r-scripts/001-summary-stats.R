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
truncationlvl<-"234" # truncation level that we chose in QIIME2
authorname<-"pooled" # name of the folder with QIIME2 output
read.end.type<-"single" # single reads or paired reads: decided in QIIME2

# If you already created a workspace, just load it
phyloseq.date_time<-"20240502_17_45_41"
load(file.path("./output/rdafiles",paste(
  phyloseq.date_time,
  authorname,read.end.type,"phyloseq-summary-stats",
  truncationlvl,
  "workspace.RData",sep = "-")))

# If not, then create it
# date_time<-"20240425_02_57_13"
# qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime"),
#                     paste(date_time,read.end.type,truncationlvl,sep="-")) # directory with QZA files
# 
# metadatadir<-file.path("./data/metadata",
#                        paste(authorname,"metadata",sep = "-")) # directory with metadata
# 
# # Specify the name of your metadata file
# metadata.filename<-file.path(metadatadir,
#                              paste("filenames",read.end.type,
#                                    authorname,"raw-supercomp.tsv", sep = "-"))
# 
# 
# ## Import qza files and convert them into a phyloseq object ####
# ps.q<-qza_to_phyloseq(
#   features = file.path(qiimedir, paste0(paste(authorname,read.end.type,
#                                               "filtered-table-trimmed-dada2",
#                                               truncationlvl,sep="-"),".qza")), # feature table
#   taxonomy = file.path(qiimedir,paste0(paste(authorname,read.end.type,
#                                              "filtered-taxonomy-trimmed-dada2",
#                                              truncationlvl,sep="-"),".qza")), # taxonomy
#   tree = file.path(qiimedir,paste0(paste(authorname,read.end.type,
#                                          "rooted-tree-trimmed-dada2",
#                                          truncationlvl,sep="-"),".qza")) # rooted tree
# )
# # Change the name d__Kingdom to Kingdom
# ps.q.taxtab<-as.data.frame(tax_table(ps.q))
# ps.q.taxtab$Kingdom<-
#   gsub("d__","",ps.q.taxtab$Kingdom)
# tax_table(ps.q)<-as.matrix(ps.q.taxtab)
# rm(ps.q.taxtab)
# 
# # Add custom metadata cause previous command loses metadata for some reason
# custom.md<-read.table(metadata.filename, header = T)
# colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
# # convert the Sample column into row names because phyloseq needs samples
# # as rownames
# # Remove absolute.filepath column
# custom.md<-custom.md%>%
#   column_to_rownames(var = "Sample")%>%
#   select(-absolute.filepath)
# custom.md$class<-as.factor(custom.md$class)
# custom.md$animal<-as.factor(custom.md$animal)
# custom.md$sex<-as.factor(custom.md$sex)
# custom.md$birthday<-as.Date(custom.md$birthday)
# # assign the custom metadata as your phyloseq object's metadata
# sample_data(ps.q)<-custom.md
# 
# # you can exclude some samples based on class. Specify the excluded classes
# # in a vector, then use the `%in%` operator. It will remove entries 
# # of the `class` column (animal hosts) from the `custom.md` object (metadata)
# # custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
# #                                               'ntccontrol','rabbitcontrol',
# #                                               'harecontrol'),]
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
# # you can exclude samples based on their library size (total number of reads)
# custom.md<-custom.md[!rownames(custom.md) %in%
#                        intersect(names(which(colSums(ps.q@otu_table)<20000)),
#                                  rownames(custom.md)),]
# 
# ### Construct the phyloseq object directly from dada2 output ####
# # We combine the phyloseq object with new metadata (if we excluded samples)
# ps.foo <- phyloseq(otu_table(ps.q),
#                    sample_data(custom.md),
#                    tax_table(ps.q),
#                    phy_tree(ps.q))
# ps.q<-ps.foo
# rm(ps.foo)
# 
# # Number of features in the unfiltered dataset
# length(rownames(ps.q@tax_table@.Data))
# 
# # Total frequency in the unfiltered dataset
# sum(colSums(ps.q@otu_table@.Data))
# 
# # Summary statistics (min, median, max, quartiles) of the unfiltered dataset
# ps.q@otu_table@.Data%>%
#   colSums()%>%
#   summary()
# 
# # Select only Bacteria
# ps.q<-ps.q %>%
#   subset_taxa(Kingdom%in%"Bacteria")
# 
# ### Fix empty taxa with higher rank taxon ####
# # Because we want to remove NA values and make ambiguous "uncultured" or 
# # "unclassified" taxa more understandable.
# ps.q<-tax_fix(ps.q,unknowns = c("NA","uncultured","Unassigned",
#                                 "uncultured_bacterium","uncultured_rumen",
#                                 "gut_metagenome","human_gut","mouse_gut",
#                                 "wallaby_gut","uncultured_soil", 
#                                 "uncultured_organism","uncultured_prokaryote"))
# ## Create a dataframe of absolute abundances ####
# ### Extract absolute abundances ####
# ps.q.agg<-ps.q %>%
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()  # transform the phyloseq object into an R dataframe
# 
# ps.q.agg.phylum<-ps.q %>%
#   tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by agglom.rank
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()
# ps.q.agg.family<-ps.q %>%
#   tax_glom("Family",NArm = FALSE) %>% # agglomerate by agglom.rank
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()
# ps.q.agg.genus<-ps.q %>%
#   tax_glom("Genus",NArm = FALSE) %>% # agglomerate by agglom.rank
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()
# 
# # Remove entries with zero Abundance
# ps.q.agg<-ps.q.agg%>%
#   filter(Abundance!=0)
# ps.q.agg.phylum<-ps.q.agg.phylum%>%
#   filter(Abundance!=0)
# ps.q.agg.family<-ps.q.agg.family%>%
#   filter(Abundance!=0)
# ps.q.agg.genus<-ps.q.agg.genus%>%
#   filter(Abundance!=0)
# 
# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   authorname,read.end.type,"phyloseq-summary-stats",
#   truncationlvl,
#   "workspace.RData",sep = "-")))

# Total ASV/phyla/families/genera per class ####
n.asv.per.host<-ps.q.agg%>%
  group_by(class,OTU)%>% # group by class (animal host),
  distinct(OTU)%>%
  group_by(class)%>%
  tally()
n.phylum.per.host<-ps.q.agg.phylum%>%
  group_by(class,Phylum)%>% # group by class (animal host),
  distinct(Phylum)%>%
  group_by(class)%>%
  tally()
n.family.per.host<-ps.q.agg.family%>%
  group_by(class,Family)%>% # group by class (animal host),
  distinct(Family)%>%
  group_by(class)%>%
  tally()
n.genus.per.host<-ps.q.agg.genus%>%
  group_by(class,Genus)%>% # group by class (animal host),
  distinct(Genus)%>%
  group_by(class)%>%
  tally()

# Summary table ####
# Number of samples per host, total reads per host, mean library size,
# sd library size, asv per host, phyla per host, families per host, 
# genera per host
summary.table<-ps.q.agg%>%
  group_by(class)%>%
  mutate(TotalSamplesPerHost=n_distinct(Sample))%>%
  mutate(TotalReadsPerHost=sum(Abundance))%>%
  group_by(Sample)%>%
  mutate(LibrarySize=sum(Abundance))%>%
  distinct(Sample,.keep_all = T)%>%
  group_by(class)%>%
  mutate(MeanLibrarySize =round(mean(LibrarySize)),
         SDLibrarySize=round(sd(LibrarySize)))%>%
  select(class,
         TotalSamplesPerHost,
         TotalReadsPerHost,
         MeanLibrarySize,
         SDLibrarySize)%>%
  distinct(class,.keep_all = T)%>%
  arrange(class)%>%
  left_join(n.asv.per.host)%>%
  rename(ASVPerHost=n)%>%
  left_join(n.phylum.per.host)%>%
  rename(PhylaPerHost=n)%>%
  left_join(n.family.per.host)%>%
  rename(FamiliesPerHost=n)%>%
  left_join(n.genus.per.host)%>%
  rename(GeneraPerHost=n)

write.table(summary.table,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "summary-table.tsv",sep="-")),
            row.names = F,sep = "\t")



### 8. Check how many taxa are unclassified in each NMR sample ####
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
agglom.rank<-"Genus"
# Do a grep to find unclassified taxa: we take the vector of all ranks from 
# above except the agglom.rank (remove from the vector), and look for the 
# remaining strings in the agglom.rank column (e.g if agglom.rank="Phylum", we
# remove "Phylum" from all.ranks and look for other ranks in the ps.q.agg.for_bp Phylum
# column, i.e Kingdom, Class, Order, Family, Genus)
ps.q.agg.genus<-ps.q.agg.genus%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample",agglom.rank))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)

# Unclassified in NMR
ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample)%>%
  filter(class=="NMR")%>%
  summarise(total=sum(RelativeAbundance))%>%
  summary()

### 8.1 Check which animals have the highest proportion of unclassified taxa ####
# In filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank], collapse='|'),get(agglom.rank))),
# we are searching agglom.rank column for taxa that have a higher rank
# in their name. For example, Bacteria Kingdom is unclassified, and if it's in
# the agglom.rank column, we find the word "Bacteria" which makes it an 
# unclassified taxon.
# In summarise(TotalUnclassified=sum(RelativeAbundance)), we check the total
# proportion of unclassified taxa in each sample, then sort to find the 
# most unclassified samples
# Summary stats: mean, SD, min, max
# Note: stats are the same for rarefied and non-rarefied data
unclassified.genus.summary.table<-ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample,class)%>%
  summarise(TotalUnclassified=sum(RelativeAbundance))%>%
  group_by(class)%>%
  mutate(MeanTotalUnclassified=round(mean(TotalUnclassified)),
         SDTotalUnclassified=round(sd(TotalUnclassified)),
         minTotalUnclassified=round(min(TotalUnclassified)),
         maxTotalUnclassified=round(max(TotalUnclassified)))%>%
  select(-Sample,-TotalUnclassified)%>%
  distinct(class,.keep_all = T)%>%
  arrange(-MeanTotalUnclassified)

write.table(unclassified.genus.summary.table,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "unclassified-genus-summary-table.tsv",sep="-")),
            row.names = F,sep = "\t")

# Top 30 samples with unclassified taxa
ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample,class)%>%
  summarise(TotalUnclassified=sum(RelativeAbundance))%>%
  arrange(-TotalUnclassified)%>%
  head(n = 30)%>%
  group_by(class)%>%
  summary()


# Check the number of distinct OTU/taxa
ps.q.agg.genus%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(get(agglom.rank))%>%
  summarise(Taxon_count=n())%>%
  arrange(-Taxon_count)

# sanity check
ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  filter(Sample=="2D10")%>%
  select(RelativeAbundance)%>%
  ungroup()%>%
  summarise(sumab=sum(RelativeAbundance))



# Most abundant phyla ####
ps.q.agg.phylum%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Phylum)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)

# Most abundant families ####
ps.q.agg.family%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Family)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Family,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)



# Biagi ####
truncationlvl<-"0" # truncation level that we chose in QIIME2
authorname<-"biagi" # name of the folder with QIIME2 output
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
# If you already created a workspace, just load it
phyloseq.biagi.date_time<-"20240502_23_43_39"
load(file.path("./output/rdafiles",paste(
  phyloseq.biagi.date_time,
  authorname,read.end.type,"phyloseq-summary-stats",
  truncationlvl,
  "workspace.RData",sep = "-")))
# qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime")) # directory with QZA files

# metadatadir<-file.path("./data/metadata",
#                        paste(authorname,"metadata",sep = "-")) # directory with metadata
# 
# # Specify the name of your metadata file
# metadata.filename<-file.path(metadatadir,
#                              paste("filenames",read.end.type,
#                                    authorname,"raw-supercomp.tsv", sep = "-"))
# 
# ## Import qza files and convert them into a phyloseq object ####
# ps.q.biagi<-qza_to_phyloseq(
#   features = file.path(qiimedir, paste0(paste(authorname,read.end.type,
#                                               "filtered-table-trimmed-dada2",
#                                               truncationlvl,sep="-"),".qza")), # feature table
#   taxonomy = file.path(qiimedir,paste0(paste(authorname,read.end.type,
#                                              "filtered-taxonomy-trimmed-dada2",
#                                              truncationlvl,sep="-"),".qza")), # taxonomy
#   tree = file.path(qiimedir,paste0(paste(authorname,read.end.type,
#                                          "rooted-tree-trimmed-dada2",
#                                          truncationlvl,sep="-"),".qza")) # rooted tree
# )
# # Change the name d__Kingdom to Kingdom
# ps.q.biagi.taxtab<-as.data.frame(tax_table(ps.q.biagi))
# ps.q.biagi.taxtab$Kingdom<-
#   gsub("d__","",ps.q.biagi.taxtab$Kingdom)
# tax_table(ps.q.biagi)<-as.matrix(ps.q.biagi.taxtab)
# rm(ps.q.biagi.taxtab)
# 
# # Add custom metadata cause previous command loses metadata for some reason
# custom.md<-read.table(metadata.filename, header = T)
# colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
# # convert the Sample column into row names because phyloseq needs samples
# # as rownames
# # Remove absolute.filepath column
# custom.md<-custom.md%>%
#   column_to_rownames(var = "Sample")%>%
#   select(-absolute.filepath)
# custom.md$class<-as.factor(custom.md$class)
# custom.md$animal<-as.factor(custom.md$animal)
# custom.md$sex<-as.factor(custom.md$sex)
# custom.md$birthday<-as.factor(custom.md$birthday)
# # assign the custom metadata as your phyloseq object's metadata
# sample_data(ps.q.biagi)<-custom.md
# 
# # you can exclude some samples based on class. Specify the excluded classes
# # in a vector, then use the `%in%` operator. It will remove entries 
# # of the `class` column (animal hosts) from the `custom.md` object (metadata)
# # custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
# #                                               'ntccontrol','rabbitcontrol',
# #                                               'harecontrol'),]
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
# # you can exclude samples based on their library size (total number of reads)
# custom.md<-custom.md[!rownames(custom.md) %in%
#                        intersect(names(which(colSums(ps.q.biagi@otu_table)<20000)),
#                                  rownames(custom.md)),]
# 
# ### Construct the phyloseq object directly from dada2 output ####
# # We combine the phyloseq object with new metadata (if we excluded samples)
# ps.foo <- phyloseq(otu_table(ps.q.biagi),
#                    sample_data(custom.md),
#                    tax_table(ps.q.biagi),
#                    phy_tree(ps.q.biagi))
# ps.q.biagi<-ps.foo
# rm(ps.foo)
# 
# # Number of features in the unfiltered dataset ####
# length(rownames(ps.q.biagi@tax_table@.Data))
# 
# # Total frequency in the unfiltered dataset ####
# sum(colSums(ps.q.biagi@otu_table@.Data))
# 
# # Summary statistics (min, median, max, quartiles) of the unfiltered dataset ####
# ps.q.biagi@otu_table@.Data%>%
#   colSums()%>%
#   summary()
# 
# # Select only Bacteria
# ps.q.biagi<-ps.q.biagi %>%
#   subset_taxa(Kingdom%in%"Bacteria")
# 
# ### Fix empty taxa with higher rank taxon ####
# # Because we want to remove NA values and make ambiguous "uncultured" or 
# # "unclassified" taxa more understandable.
# ps.q.biagi<-tax_fix(ps.q.biagi,unknowns = c("NA","uncultured","Unassigned",
#                                 "uncultured_bacterium","uncultured_rumen",
#                                 "gut_metagenome","human_gut","mouse_gut",
#                                 "wallaby_gut","uncultured_soil", 
#                                 "uncultured_organism","uncultured_prokaryote"))
# ## Create a dataframe of absolute abundances ####
# ### Extract absolute abundances ####
# ps.q.biagi.agg<-ps.q.biagi %>%
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()  # transform the phyloseq object into an R dataframe
# 
# ps.q.biagi.agg.phylum<-ps.q.biagi %>%
#   tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by agglom.rank
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()
# ps.q.biagi.agg.family<-ps.q.biagi %>%
#   tax_glom("Family",NArm = FALSE) %>% # agglomerate by agglom.rank
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()
# ps.q.biagi.agg.genus<-ps.q.biagi %>%
#   tax_glom("Genus",NArm = FALSE) %>% # agglomerate by agglom.rank
#   subset_taxa(Kingdom%in%"Bacteria")%>% # choose only bacteria
#   psmelt()
# 
# # Remove entries with zero Abundance
# ps.q.biagi.agg<-ps.q.biagi.agg%>%
#   filter(Abundance!=0)
# ps.q.biagi.agg.phylum<-ps.q.biagi.agg.phylum%>%
#   filter(Abundance!=0)
# ps.q.biagi.agg.family<-ps.q.biagi.agg.family%>%
#   filter(Abundance!=0)
# ps.q.biagi.agg.genus<-ps.q.biagi.agg.genus%>%
#   filter(Abundance!=0)


# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   authorname,read.end.type,"phyloseq-summary-stats",
#   truncationlvl,
#   "workspace.RData",sep = "-")))

# Number of samples per class ####
n.samples.per.host.biagi<-ps.q.biagi.agg%>%
  group_by(class)%>%
  distinct(Sample)%>%
  tally()

# Total reads per class ####
n.reads.per.host.biagi<-ps.q.biagi.agg%>%
  group_by(class)%>%
  summarise(TotalReadsPerHost=sum(Abundance))

# Total ASV/phyla/families/genera per class ####
n.asv.per.host.biagi<-ps.q.biagi.agg%>%
  group_by(class,OTU)%>% # group by class (animal host),
  distinct(OTU)%>%
  group_by(class)%>%
  tally()
n.phylum.per.host.biagi<-ps.q.biagi.agg.phylum%>%
  group_by(class,Phylum)%>% # group by class (animal host),
  distinct(Phylum)%>%
  group_by(class)%>%
  tally()
n.family.per.host.biagi<-ps.q.biagi.agg.family%>%
  group_by(class,Family)%>% # group by class (animal host),
  distinct(Family)%>%
  group_by(class)%>%
  tally()
n.genus.per.host.biagi<-ps.q.biagi.agg.genus%>%
  group_by(class,Genus)%>% # group by class (animal host),
  distinct(Genus)%>%
  group_by(class)%>%
  tally()

summary.table.biagi<-ps.q.biagi.agg.genus%>%
  group_by(class)%>%
  mutate(TotalSamplesPerHost=n_distinct(Sample))%>%
  mutate(TotalReadsPerHost=sum(Abundance))%>%
  group_by(Sample)%>%
  mutate(LibrarySize=sum(Abundance))%>%
  distinct(Sample,.keep_all = T)%>%
  group_by(class)%>%
  mutate(MeanLibrarySize =mean(LibrarySize),
         SDLibrarySize=sd(LibrarySize))%>%
  select(class,
         TotalSamplesPerHost,
         TotalReadsPerHost,
         MeanLibrarySize,
         SDLibrarySize)%>%
  distinct(class,.keep_all = T)%>%
  arrange(class)%>%
  left_join(n.asv.per.host.biagi)%>%
  rename(ASVPerHost=n)%>%
  left_join(n.phylum.per.host.biagi)%>%
  rename(PhylaPerHost=n)%>%
  left_join(n.family.per.host.biagi)%>%
  rename(FamiliesPerHost=n)%>%
  left_join(n.genus.per.host.biagi)%>%
  rename(GeneraPerHost=n)


write.table(summary.table.biagi,
            file=file.path("./output/rtables","biagi",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "summary-table.tsv",sep="-")),
            row.names = F,sep = "\t")

# Most abundant phyla ####
ps.q.biagi.agg.phylum.dominant<-ps.q.biagi.agg.phylum%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Phylum)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)

# Most abundant families ####
ps.q.biagi.agg.family.dominant<-ps.q.biagi.agg.family%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Family)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Family,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)

ps.q.biagi.agg.family.dominant<-ps.q.biagi.agg.family.dominant%>%
  head(n=10)
save.image(paste0("./output/rdafiles/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "pooled-and-biagi-workspace.Rda")))



# Explore phyla in NMR and mice ####
ps.q.agg.phylum.dominant<-ps.q.agg.phylum%>%
  # filter(class%in%c("NMR","B6mouse"))%>%
  filter(class%in%c("NMR"))%>%
  group_by(class)%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(class,Phylum)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)


# Explore families in NMR and mice ####
ps.q.agg.family.dominant<-ps.q.agg.family%>%
  # filter(class%in%c("NMR","B6mouse"))%>%
  filter(class%in%c("NMR"))%>%
  group_by(class)%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(class,Family)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Family,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
ps.q.agg.family.dominant<-ps.q.agg.family.dominant%>%
  head(n=10)

ps.q.agg.family.dominant%>%
  ungroup()%>%
  select(Family,MeanRelativeAbundance)%>%
  full_join(ps.q.biagi.agg.family.dominant[,c("Family","MeanRelativeAbundance")],
            by="Family")%>%
  rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
  pivot_longer(cols = c("NMR","NMRwt"),
               names_to = "class",
               values_to = "MeanRelativeAbundance")%>%
  arrange(-MeanRelativeAbundance)%>%
  head(n=30)%>%
  mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
  ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
             y=MeanRelativeAbundance,
             fill=animal))+
  geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+
  xlab("") +
  ylab("Average relative abundance (%)")+
  labs(fill="")+
  coord_cartesian(expand = FALSE)+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), #TODO: what does it do?
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.title = element_text(size = 25), # size of legend title
        legend.position = c(0.9,0.8),
        legend.text = element_text(size = 20)) # legend under the plot
ggsave(file.path("./images/barplots",
                 paste(paste(format(Sys.time(),format="%Y%m%d"),
                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                       "barplot-wild-vs-lab-nmr-my-results.png")),
       plot=last_plot(),
       width = 5000,height = 4000,
       units = "px",dpi=300,device = "png")  



wild.nmr.families<-data.frame(
  Family=c("Lachnospiraceae",
           "Prevotellaceae",
           "Paraprevotellaceae",
           "Bacteroidales Order",
           "Clostridiales Order",
           "Ruminococcaceae",
           "Veillonellacaea",
           "Clostridiacaea",
           "Muribaculaceae",
           "Porphyromonadaceae"),
  MeanRelativeAbundance=c(17.6,
                          11,
                          8.8,
                          6.2,
                          6.1,
                          5.7,
                          4.7,
                          4.1,
                          4,
                          3)
)

ps.q.agg.family.dominant%>%
  ungroup()%>%
  select(Family,MeanRelativeAbundance)%>%
  full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
            by="Family")%>%
  rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
  pivot_longer(cols = c("NMR","NMRwt"),
               names_to = "class",
               values_to = "MeanRelativeAbundance")%>%
  arrange(-MeanRelativeAbundance)%>%
  head(n=30)%>%
  mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
  replace_na(replace =list( MeanRelativeAbundance=0))%>%
  ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
             y=MeanRelativeAbundance,
             fill=animal))+
  geom_bar(stat = "identity",position = position_dodge2())+
  theme_bw()+
  xlab("") +
  ylab("Average relative abundance (%)")+
  labs(fill="")+
  coord_cartesian(expand = FALSE)+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), #TODO: what does it do?
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.title = element_text(size = 25), # size of legend title
        legend.position = c(0.9,0.8),
        legend.text = element_text(size = 20)) # legend under the plot
ggsave(file.path("./images/barplots",
                 paste(paste(format(Sys.time(),format="%Y%m%d"),
                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                       "barplot-wild-vs-lab-nmr-from-paper.png")),
       plot=last_plot(),
       width = 5000,height = 4000,
       units = "px",dpi=300,device = "png")  
