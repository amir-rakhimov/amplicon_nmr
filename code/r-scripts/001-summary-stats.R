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

# write.table(summary.table,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")



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
# In summarise(TotalUnclassifiedPercent=sum(RelativeAbundance)), we check the total
# proportion of unclassified taxa in each sample, then sort to find the 
# most unclassified samples
# Summary stats: mean, SD, min, max
# Note: stats are the same for rarefied and non-rarefied data
unclassified.genus.summary.table<-ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample,class)%>%
  summarise(TotalUnclassifiedPercent=sum(RelativeAbundance))%>%
  group_by(class)%>%
  mutate(MeanTotalUnclassifiedPercent=round(mean(TotalUnclassifiedPercent)),
         SDTotalUnclassifiedPercent=round(sd(TotalUnclassifiedPercent)),
         minTotalUnclassifiedPercent=round(min(TotalUnclassifiedPercent)),
         maxTotalUnclassifiedPercent=round(max(TotalUnclassifiedPercent)),
         MedianTotalUnclassifiedPercent=round(median(TotalUnclassifiedPercent)))%>%
  select(-Sample,-TotalUnclassifiedPercent)%>%
  distinct(class,.keep_all = T)%>%
  arrange(-MeanTotalUnclassifiedPercent)

# write.table(unclassified.genus.summary.table,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

pretty.level.names<-c("NMR" = "Naked mole-rat", # better labels for facets
                       "B6mouse" = "B6 mouse",
                       "MSMmouse" = "MSM/Ms mouse",
                       "FVBNmouse" = "FVB/N mouse",
                       "DMR" = "Damaraland mole-rat",
                       "hare" = "European rabbit",
                       "rabbit" = "European brown hare",
                       "spalax" = "Spalax (blind mole-rat)",
                       "pvo" = "Siberian flying squirrel")
# Use only the taxa that are present in the workspace
# (custom.md is metadata from the rdafile)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
barplot.directory<-"./images/barplots/" # set the path where barplots will
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")

summary.table%>%
  ggplot(aes(x=reorder(class,-GeneraPerHost),y=GeneraPerHost))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x="",
       y="Genera per host")+
  scale_x_discrete(labels =pretty.level.names) +
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), #TODO: what does it do?
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "bottom")

# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "barplot-total-genera",paste(custom.levels,collapse = '-'),
#                       truncationlvl,agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=last_plot(),
#          width = 4000,height = 3000,
#          units = "px",dpi=300,device = image.format)
# }
  


# unclassified.genus.summary.table%>%
#   ggplot(aes(x=reorder(class,-MeanTotalUnclassifiedPercent),y=MeanTotalUnclassifiedPercent))+
#   geom_bar(stat = "identity")+
#   theme_bw()+
#   labs(x="",
#        y="Unclassified genera per sample (%)")+
#   scale_x_discrete(labels =pretty.level.names) +
#   theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
#         axis.line = element_blank(), #TODO: what does it do?
#         panel.spacing = unit(0.8, "cm"), # increase distance between facets
#         axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
#         # the x-axis labels by 45 degrees and shift to the right
#         axis.text.y = element_text(size=20), # size of y axis ticks
#         axis.title = element_text(size = 20), # size of axis names
#         plot.title = element_text(size = 25), # size of plot title
#         plot.caption = element_text(size=23), # size of plot caption
#         legend.text = element_text(size = 20), # size of legend text
#         legend.title = element_text(size = 25), # size of legend title
#         legend.position = "bottom")
# 
# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "barplot-unclas-genera",paste(custom.levels,collapse = '-'),
#                       truncationlvl,agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=last_plot(),
#          width = 4000,height = 3000,
#          units = "px",dpi=300,device = image.format)
# }


# boxplot of genera
ps.q.agg.genus%>%
  group_by(Sample,class,Genus)%>% # group by class (animal host),
  distinct(Genus)%>%
  group_by(class,Sample)%>%
  tally%>%
  ggplot(aes(x=class,y=n))+
  geom_boxplot()

# Arrange boxplot by mean unclassified genera
mean.factors<-unclassified.genus.summary.table%>%
  arrange(-MeanTotalUnclassifiedPercent)%>%
  pull(class)%>%
  as.character()
ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample,class)%>%
  summarise(TotalUnclassifiedPercent=sum(RelativeAbundance))%>%
  ggplot(aes(x=factor(class,levels=mean.factors),y=TotalUnclassifiedPercent))+
  geom_boxplot()+
  stat_summary(fun=mean, geom='point', shape=23,fill="red",size=4)+
  theme_bw()+
  labs(x="",
       y="Unclassified genera per sample (%)")+
  scale_x_discrete(labels =pretty.level.names) +
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), #TODO: what does it do?
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "bottom")
# for(image.format in image.formats){
#   ggsave(paste0(boxplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "boxplot-unclassified-genera",paste(custom.levels,collapse = '-'),
#                       truncationlvl,agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=last_plot(),
#          width = 4000,height = 3000,
#          units = "px",dpi=300,device = image.format)
# }



# Top 30 samples with unclassified taxa
ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample,class)%>%
  summarise(TotalUnclassifiedPercent=sum(RelativeAbundance))%>%
  arrange(-TotalUnclassifiedPercent)%>%
  head(n = 30)%>%
  group_by(class)%>%
  summary()


# Check the number of distinct OTU/taxa
num.distinct.genera<-ps.q.agg.genus%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(get(agglom.rank))%>%
  summarise(Genus_count=n())%>%
  arrange(-Genus_count)
write.table(num.distinct.genera,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "num-distinct-genera-summary-table.tsv",sep="-")),
            row.names = F,sep = "\t")


# sanity check
ps.q.agg.genus%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  filter(Sample=="2D10")%>%
  select(RelativeAbundance)%>%
  ungroup()%>%
  summarise(sumab=sum(RelativeAbundance))



# Most abundant phyla in NMR ####
ps.q.agg.dominant.phyla<-ps.q.agg.phylum%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Phylum)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.phyla,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-phyla-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# Most abundant families in NMR  ####
ps.q.agg.dominant.families<-ps.q.agg.family%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Family)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Family,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.families,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-families-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# Most abundant genera in NMR  ####
ps.q.agg.dominant.genera<-ps.q.agg.genus%>%
  ungroup()%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Genus)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Family,Genus,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.genera,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-genera-table.tsv",sep="-")),
#             row.names = F,sep = "\t")



# Compare with Debebe et al. ####
ps.q.agg.family.relab<-ps.q.agg.family%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample","Family"))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Family"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)

ps.q.agg.genus.relab<-ps.q.agg.genus%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample","Genus"))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Genus"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  
# How much Bacteroidaceae are in NMR  ####
bacteroidaceae.nmr<-ps.q.agg.family.relab%>%
  filter(Family=="Bacteroidaceae",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()
write.table(bacteroidaceae.nmr,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "bacteroidaceae-nmr-table.tsv",sep="-")),
            row.names = F,sep = "\t")

# What are the most dominant Bacteroidota families in NMR ####
bacteroidota.nmr<-ps.q.agg.family.relab%>%
  filter(Phylum=="Bacteroidota",class=="NMR")%>%
  group_by(Family)%>%
  distinct(Family,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Family,MeanRelativeAbundance)
write.table(bacteroidota.nmr,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "bacteroidota-nmr-table.tsv",sep="-")),
            row.names = F,sep = "\t")

# Spirochaetaceae and Treponema ####
spirochaetaceae.nmr<-ps.q.agg.family.relab%>%
  filter(Family=="Spirochaetaceae",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
            max=max(RelativeAbundance),
            mean=TotalAgglomRank/TotalClass*100,
            sd=sd(RelativeAbundance),
            n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()
write.table(spirochaetaceae.nmr,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "spirochaetaceae-nmr-table.tsv",sep="-")),
            row.names = F,sep = "\t")

treponema.nmr<-ps.q.agg.genus.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Genus,min,max,mean,sd,n)%>%
  distinct()
write.table(treponema.nmr,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "treponema-nmr-table.tsv",sep="-")),
            row.names = F,sep = "\t")

# How many ASVs in Treponema ####
ps.q.agg%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

# Mogibacteriaceae is renamed to Anaerovoracaceae ####
mogibacteriaceae_anaerovoracaceae.all<-ps.q.agg.family.relab%>%
  filter(Family=="Anaerovoracaceae")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()%>%
  arrange(-mean)
write.table(mogibacteriaceae_anaerovoracaceae.all,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
            row.names = F,sep = "\t")

# sulfur metabolising bacteria ####
desulfobacterota.nmr<-ps.q.agg.genus.relab%>%
  filter(Phylum=="Desulfobacterota",class=="NMR")%>%
  group_by(Genus)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Genus,MeanRelativeAbundance,sd)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
write.table(desulfobacterota.nmr,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "desulfobacterota-nmr-table.tsv",sep="-")),
            row.names = F,sep = "\t")


# Plot Desulfobacterota
library(Polychrome)
library(ggtext)
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "B6mouse" = "B6 mouse",
                       "MSMmouse" = "MSM/Ms mouse",
                       "FVBNmouse" = "FVB/N mouse",
                       "DMR" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pvo" = "*Pteromys volans orii*",
                       "NMRwt"="Wild *Heterocephalus glaber*"
)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

if(exists("excluded.samples")){
  custom.levels<-custom.levels[!custom.levels%in%excluded.samples]
  pretty.level.names<-
    pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  pretty.level.names<-pretty.level.names[!names(pretty.level.names)%in%excluded.samples]
}else{
  pretty.level.names<-
    pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
}
set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900",
                                          "#FFA07A"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

# It's a flipped plot
ps.q.agg.genus.relab%>%
  filter(Phylum=="Desulfobacterota")%>%
  group_by(class,Genus)%>%
  ggplot(aes(x=factor(class,level=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE)+
  facet_wrap(~Genus,scales = "free_x",
             ncol = 2)+
  theme_bw()+
  coord_flip()+
  labs(x="",
       y="Relative abundance (%)")+
  scale_color_manual(breaks = rev(unname(pretty.level.names)),
                     labels=rev(unname(pretty.level.names)))+
  scale_x_discrete(labels=rev(pretty.level.names),
                   limits=rev(custom.levels))+ # rename boxplot labels (x axis)
  scale_fill_manual(values = rev(custom.fill))+
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "none")+
  ggtitle(paste0("Relative abundance of Desulfobacterota phylum members"))

for (image_format in c("png","tiff")){
  ggsave(paste0("./images/taxaboxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "Desulfobacterota phylum members-all-hosts",
                      sep = "-"),".",image_format),
         plot=last_plot(),
         width = 4000,height = 6000,
         units = "px",dpi=300,device = image_format)
}
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
