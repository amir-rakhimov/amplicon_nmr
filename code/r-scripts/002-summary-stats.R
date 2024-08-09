# In this script, we are exploring the imported dataset from QIIME2.
# We do not use the data from 001-phyloseq-qiime2.R script because we will 
# create multiple agglomerated tables at phylum, family, genus, and OTU level.

## 1. Import libraries ####
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
## 2. Specifying parameters and directory/file names #### 
truncationlvl<-"234" # truncation level that we chose in QIIME2
# truncationlvl<-"0" # truncation level that we chose in QIIME2

authorname<-"pooled" # name of the folder with QIIME2 output
# authorname<-"biagi" # name of the folder with QIIME2 output

read.end.type<-"single" # single reads or paired reads: decided in QIIME2

# If you already created a workspace, just load it
if(authorname=="pooled"){
  phyloseq.workspace.date_time<-"20240502_17_45_41"
}else if(authorname =="biagi"){
  phyloseq.workspace.date_time<-"20240502_23_43_39" # Debebe phyloseq workspace
}else{
  warning("Phyloseq workspace not found")
}

load(file.path("./output/rdafiles",paste(
  phyloseq.workspace.date_time,
  authorname,read.end.type,"phyloseq-summary-stats",
  truncationlvl,
  "workspace.RData",sep = "-")))

barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")

# # If not, then create it
# if(authorname=="pooled"){
#   qza_file_date_time<-"20240425_02_57_13"
#   qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime"),
#                       paste(qza_file_date_time,read.end.type,truncationlvl,sep="-")) # directory with QZA files
# }else if(authorname=="biagi"){
#   qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime")) # directory with QZA files for Debebe
# }else{
#   stop("QZA file directory not found")
# }
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
# ## 3. Import qza files and convert them into a phyloseq object ####
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
# # Add custom metadata
# custom.md<-read.table(metadata.filename, header = T)
# colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
# # convert the Sample column into row names because phyloseq needs samples
# # as rownames
# # Remove absolute.filepath column
# rownames(custom.md)<-custom.md$Sample
# custom.md<-custom.md%>%
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
# ### 3.1 Construct the phyloseq object directly from dada2 output ####
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
# ### 3.2 Fix empty taxa with higher rank taxon ####
# # Because we want to remove NA values and make ambiguous "uncultured" or
# # "unclassified" taxa more understandable.
# ps.q<-tax_fix(ps.q,unknowns = c("NA","uncultured","Unassigned",
#                                 "uncultured_bacterium","uncultured_rumen",
#                                 "gut_metagenome","human_gut","mouse_gut",
#                                 "wallaby_gut","uncultured_soil",
#                                 "uncultured_organism","uncultured_prokaryote"))
# ## 4. Create a dataframe of absolute abundances ####
# # Extract absolute abundances
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
#   filter(Abundance!=0)%>%
#   select(-sample_Sample) # remove the duplicate column
# ps.q.agg.phylum<-ps.q.agg.phylum%>%
#   filter(Abundance!=0)%>%
#   select(-sample_Sample) # remove the duplicate column
# ps.q.agg.family<-ps.q.agg.family%>%
#   filter(Abundance!=0)%>%
#   select(-sample_Sample) # remove the duplicate column
# ps.q.agg.genus<-ps.q.agg.genus%>%
#   filter(Abundance!=0)%>%
#   select(-sample_Sample) # remove the duplicate column
# # Rename objects if you're using data from other papers (not pooled data)
# if(authorname!="pooled"){
#   for(workspace_obj in ls()[grep("ps.q.",ls())]){
#     new_workspace_obj<-paste(workspace_obj,authorname,sep = ".")
#     assign(new_workspace_obj,get(workspace_obj))
#   }
# }
# # 5. Save the workspace ####
# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   authorname,read.end.type,"phyloseq-summary-stats",
#   truncationlvl,
#   "workspace.RData",sep = "-")))
# Data creation finished here ^^^^^ ####

## Setup plots ####
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


# 6. Check the total number of unique ASV/phyla/families/genera per class ####
# The function get.n.uniq.taxa.per.host groups the dataframe of abundances
# by host and taxonomic rank (e.g. Genus), retains unique rows (unique taxa in 
# each host), then groups by class and counts observations. It returns the
# number of unique taxa (e.g. genera) per host. Works for any number of hosts.
get.n.uniq.taxa.per.host<-function(tax.df,tax.rank){
  n.taxa.per.host<-tax.df%>%
    group_by_at(c("class",tax.rank))%>%
    distinct(get(tax.rank))%>%
    group_by(class)%>%
    tally
  return(n.taxa.per.host)
}
n.asv.per.host<-get.n.uniq.taxa.per.host(ps.q.agg,"OTU")
n.phylum.per.host<-get.n.uniq.taxa.per.host(ps.q.agg.phylum,"Phylum")
n.family.per.host<-get.n.uniq.taxa.per.host(ps.q.agg.family,"Family")
n.genus.per.host<-get.n.uniq.taxa.per.host(ps.q.agg.genus,"Genus")

# 7. Summary statistics table ####
# The function create_summary_stats_table takes a dataframe
# with ASV abundances (tax.df), groups by class (animal host), and adds a column with 
# number of samples per host (TotalSamplesPerHost). It uses a dplyr function n_distinct
# 2. Then, the function adds a column with number of total reads per host (TotalReadsPerHost).
# It uses a function sum() on Abundance column
# 3. Then, the function groups by Sample and adds a column with number of total reads 
# per sample (LibrarySize). It uses sum() on Abundance column
# 4. Then, the function keeps unique samples with distinct() function
# 5. Then, the function groups by class and adds a column with average number of reads
# per sample (MeanLibrarySize). It uses mean() function on LibrarySize.
# 6. The function also adds a Standard deviation of library size (SDLibrarySize). The 
# rationale is the same as when calculating MeanLibrarySize, except the function is
# sd()
# 7. The function selects columns class, TotalSamplesPerHost, TotalReadsPerHost, 
# MeanLibrarySize, and SDLibrarySize.
# 8. The function selects unique rows (distinct() function).
# 9. The function sorts by class column, then adds several dataframes:
# 9.1: Dataframe with the number of unique ASVs per host (n.asv.table).
# It is added as `n` column, but is renamed to ASVPerHost.
# 9.2: Dataframe with the number of unique phyla per host (n.phyllum.table).
# It is added as `n` column, but is renamed to PhylaPerHost.
# 9.3: Dataframe with the number of unique families per host (n.family.table).
# It is added as `n` column, but is renamed to FamiliesPerHost.
# 9.4: Dataframe with the number of unique genera per host (n.genus.table).
# It is added as `n` column, but is renamed to GeneraPerHost.
create_summary_stats_table<-function(tax.df,
                                     n.asv.table,
                                     n.phyllum.table,
                                     n.family.table,
                                     n.genus.table){
  final.summary.stats.table<-tax.df%>%
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
    left_join(n.asv.table)%>%
    rename(ASVPerHost=n)%>%
    left_join(n.phyllum.table)%>%
    rename(PhylaPerHost=n)%>%
    left_join(n.family.table)%>%
    rename(FamiliesPerHost=n)%>%
    left_join(n.genus.table)%>%
    rename(GeneraPerHost=n)
  return(final.summary.stats.table)
}
summary.stats.table<-create_summary_stats_table(ps.q.agg,
                           n.asv.per.host,
                           n.phylum.per.host,
                           n.family.per.host,
                           n.genus.per.host)

# write.table(summary.stats.table,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")


# 8. Add relative abundance and average relative abundance columns ####
# The function add_relab_to_tax_df: 1. Takes a table with absolute abundances
# 2. Groups rows by class and samples
# 3. Then adds a column with total abundance (sum of all reads) per
# sample (TotalSample, which is a sum of Abundance column per sample per host)
# 4. Then groups by class, sample, and a the lowest taxonomic rank (tax.rank)
# 5. Then adds a column with relative abundances, which are absolute abundances
# of each taxon in a sample divided by the total number of reads in a sample,
# multiplied by 100.
# 6. Then, the dataframe is grouped by class 
# 7. And we add a TotalClass column, which is a sum of reads from all samples
# belonging to a class.
# 8. Then, we group by class and the lowest taxonomic rank (tax.rank).
# 9. We add a TotalAgglomRank column, which is a sum of all reads of each 
# tax.rank in each sample.
# 8. Finally, the MeanRelativeAbundance column is the total number of reads
# belonging to a certain taxon in a host divided by the total number of reads
# from that host, multiplied by 100.
# (MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
add_relab_to_tax_df<-function(tax.df,tax.rank){
  tax.df<-tax.df%>%
    group_by(class,Sample)%>%
    mutate(TotalSample=sum(Abundance))%>%
    group_by_at(c("class","Sample",tax.rank))%>%
    mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
    group_by(class)%>% # group by class (animal host),
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",tax.rank))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  return(tax.df)
}

ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")

# 9. Add agegroup (Must run for plotting) ####
custom.md.ages<-readRDS("./output/rdafiles/custom.md.ages.rds")
# we create these new levels for plots
# The function add_agegroup_to_tax_df takes a dataframe of abundances,
# the lowest taxonomic rank in the dataframe (e.g. genus), and a metadata with
# age groups as input. It joins the dataframe with metadata, groups by age 
# groups from the metadata, then sums the abundances from each age group.
# These sums are added as a column TotalAgegroup. Then, the function 
# groups by age group and the lowest taxonomic rank, and sums abundances for each
# taxonomic rank in each age group (TotalAgglomRankAge column). Finally,
# the function calculates the average relative abundance of each taxonomic rank
# inside each age group.
add_agegroup_to_tax_df<-function(tax.df,tax.rank,metadata.df){
  tax.df<-tax.df%>%
    left_join(metadata.df)%>%
    group_by(agegroup)%>% # group by class (animal host),
    mutate(TotalAgegroup=sum(Abundance))%>%
    group_by_at(c("agegroup",tax.rank))%>%
    mutate(TotalAgglomRankAge=sum(Abundance))%>%
    mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)%>%
    ungroup()%>%
    select(-TotalAgegroup,-TotalAgglomRankAge)
}

ps.q.agg.genus.relab<-add_agegroup_to_tax_df(ps.q.agg.genus.relab,"Genus",
                                             custom.md.ages)
ps.q.agg.family.relab<-add_agegroup_to_tax_df(ps.q.agg.family.relab,"Family",
                                              custom.md.ages)

# 10. Check the percentage of unclassified taxa in each animal ####
# Here, we are interested in unclassified genera, but you can also try Families,
# Orders, Classes, etc.
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
agglom.rank<-"Genus"
# In filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank], collapse='|'),get(agglom.rank))),
# we are searching agglom.rank column for taxa that have a higher rank
# in their name. For example, Bacteria Kingdom is unclassified, and if it's in
# the agglom.rank column, we find the word "Bacteria" which makes it an 
# unclassified taxon.
# In summarise(TotalUnclassifiedPercent=sum(RelativeAbundance)), we sum the 
# relative abundance of all taxa that were unclassified in a given sample
# from a given host (because we group by sample and class).
# We calculate Mean, SD, min, max, and median of unclassified percentages.
# After that, there is no need for the TotalUnclassifiedPercent column because 
# the summary statistics are calculated for every host. So, we keep unique rows
# Note: stats are the same for rarefied and non-rarefied data
unclassified.genus.summary.stats.table<-ps.q.agg.genus.relab%>%
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

# write.table(unclassified.genus.summary.stats.table,
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
# 11. Plot the total number of unique genera per host ####
summary.stats.table%>%
  ggplot(aes(x=reorder(class,-GeneraPerHost),y=GeneraPerHost))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x="",
       y="Genera per host")+
  scale_x_discrete(labels =pretty.level.names) + # x axis labels become 
  # pretty.level.names because the names in the vector are the same as the 
  # x axis ticks in the original plot
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
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

# Boxplot that shows the distribution of genera in samples from each host ####
ps.q.agg.genus%>%
  group_by(Sample,class,Genus)%>% # group by class (animal host),
  distinct(Genus)%>%
  group_by(class,Sample)%>%
  tally%>%
  ggplot(aes(x=class,y=n))+
  geom_boxplot()

# Arrange boxplot by mean unclassified genera
# We take the table of unclassified genera stats and sort it by the average
# percentage of unclassified genera. Then, extract the vector of class names
mean.factors<-unclassified.genus.summary.stats.table%>%
  arrange(-MeanTotalUnclassifiedPercent)%>%
  pull(class)%>%
  as.character()
# Boxplot that shows the distribution of percentages of unclassified genera
# in each host
ps.q.agg.genus.relab%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(Sample,class)%>%
  summarise(TotalUnclassifiedPercent=sum(RelativeAbundance))%>%
  ggplot(aes(x=factor(class,levels=mean.factors),y=TotalUnclassifiedPercent))+
  geom_boxplot()+
  stat_summary(fun=mean, geom='point', shape=23,fill="red",size=4)+ # add a point
  # that shows the mean value in each box
  theme_bw()+
  labs(x="",
       y="Unclassified genera per sample (%)")+
  scale_x_discrete(labels =pretty.level.names) +
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23) # size of plot caption
        )
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



# sanity check
ps.q.agg.genus.relab%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  filter(Sample=="2D10")%>%
  select(RelativeAbundance)%>%
  ungroup()%>%
  summarise(sumab=sum(RelativeAbundance))

# 12. Check the most abundant phyla, families, genera in NMR ####
# The function get_dominant_taxa_in_host uses a dataframe with abundances, 
# the taxonomic rank for which you want to find dominant taxa (e.g. dominant genera), 
# and the host name/names.
# If there is on host, the function filters the dataframe, adds a column
# with the total abundance (sum of all reads in the host), then groups
# by taxonomic rank (e.g. genus) and calculates the total abundance of each taxon.
# After that, the function adds a column with average relative abundance 
# of each taxon.
# If there are two or more hosts, the procedure is the same, except we 
# group by class before counting total abundance and group by both class and 
# taxonomic rank before calculating total abundance of each taxon.
# The function selects the necessary columns like Phylum or Family, sorts
# by average relative abundance and returns a filtered dataframe.
get_dominant_taxa_in_host<-function(tax.df,tax.rank,host){
  if(length(host)==1){
    tax.df<-tax.df%>%
      ungroup()%>%
      filter(class==host)%>%
      mutate(TotalClass=sum(Abundance))%>%
      group_by_at(c(tax.rank))%>%
      mutate(TotalAgglomRank=sum(Abundance))%>%
      mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
    
  }else if (length(host)>1){
    tax.df<-tax.df%>%
      filter(class%in%host)%>%
      group_by(class)%>%
      mutate(TotalClass=sum(Abundance))%>%
      group_by_at(c("class",tax.rank))%>%
      mutate(TotalAgglomRank=sum(Abundance))%>%
      mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  }
  if(tax.rank=="Phylum"){
    tax.df<-tax.df%>%
      select(all_of(c(tax.rank,"MeanRelativeAbundance")))
  }else if(tax.rank=="Family"){
    tax.df<-tax.df%>%
      select(Phylum,Family,MeanRelativeAbundance)
  }else if (tax.rank=="Genus"){
    tax.df<-tax.df%>%
      select(Phylum,Family,Genus,MeanRelativeAbundance)
  }
  tax.df<-tax.df%>%
    distinct()%>%
    arrange(-MeanRelativeAbundance)
  return(tax.df)
}


ps.q.agg.dominant.phyla.nmr<-get_dominant_taxa_in_host(ps.q.agg.phylum,
                                                       "Phylum","NMR")
ps.q.agg.dominant.families.nmr<-get_dominant_taxa_in_host(ps.q.agg.family,
                                                          "Family","NMR")
ps.q.agg.dominant.genera.nmr<-get_dominant_taxa_in_host(ps.q.agg.genus,
                                                        "Genus","NMR")
# write.table(ps.q.agg.dominant.phyla,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-phyla-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# write.table(ps.q.agg.dominant.families,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-families-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# write.table(ps.q.agg.dominant.genera,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-genera-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")


# 13. Load data from Debebe ####
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

# 14. Total number of unique ASV/phyla/families/genera per class in Debebe ####
n.asv.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.biagi,"OTU")
n.phylum.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.phylum.biagi,"Phylum")
n.family.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.family.biagi,"Family")
n.genus.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.genus.biagi,"Genus")

# 15. Summary statistics table from Debebe ####
summary.stats.table.biagi<-create_summary_stats_table(ps.q.agg.genus.biagi,
                                                      n.asv.per.host.biagi,
                                                      n.phylum.per.host.biagi,
                                                      n.family.per.host.biagi,
                                                      n.genus.per.host.biagi)

# write.table(summary.stats.table.biagi,
#             file=file.path("./output/rtables","biagi",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 17. Most abundant phyla, families, and genera in Debebe ####
ps.q.agg.dominant.phyla.biagi<-get_dominant_taxa_in_host(ps.q.agg.phylum.biagi,
                                                       "Phylum","NMRwt")
ps.q.agg.dominant.families.biagi<-get_dominant_taxa_in_host(ps.q.agg.family.biagi,
                                                          "Family","NMRwt")
ps.q.agg.dominant.genera.nmr<-get_dominant_taxa_in_host(ps.q.agg.genus,
                                                        "Genus","NMR")
# save.image(paste0("./output/rdafiles/",
#                   paste(paste(format(Sys.time(),format="%Y%m%d"),
#                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                         "pooled-and-biagi-workspace.Rda")))




# 18. Plot dominant families in lab vs wild NMR from our data (wild data reanalysed) ####
ps.q.agg.dominant.families.nmr%>%
  ungroup()%>%
  select(Family,MeanRelativeAbundance)%>%
  full_join(ps.q.agg.dominant.families.biagi[,c("Family","MeanRelativeAbundance")],
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
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "inside",
        legend.position.inside =  c(0.9,0.8),
        legend.text = element_text(size = 20)) # legend under the plot
# ggsave(file.path("./images/barplots",
#                  paste(paste(format(Sys.time(),format="%Y%m%d"),
#                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                        "barplot-wild-vs-lab-nmr-my-results.png")),
#        plot=last_plot(),
#        width = 5000,height = 4000,
#        units = "px",dpi=300,device = "png")  



## 18.1 Plot dominant families reported in the paper ####
wild.nmr.families<-data.frame(
  Family=c("Lachnospiraceae",
           "Prevotellaceae",
           "Paraprevotellaceae",
           "Bacteroidales Order",
           "Clostridiales Order",
           "Oscillospiraceae",
           "Veillonellaceae",
           "Clostridiaceae",
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
# Add a column if the family was reported as dominant or not
ps.q.agg.dominant.families.nmr.with_domin<-ps.q.agg.dominant.families.nmr%>%
  ungroup()%>%
  mutate(row.index=as.numeric(rownames(ps.q.agg.dominant.families.nmr)),
         dominant=ifelse(row.index<=10,"dominant","not_dominant"))

ps.q.agg.dominant.families.nmr.with_domin%>%
  select(Family,MeanRelativeAbundance,dominant)%>%
  distinct(Family,.keep_all = T)%>%
  filter(Family%in%wild.nmr.families$Family)%>%
  full_join(ps.q.agg.dominant.families.nmr.with_domin[1:10,c("Family","MeanRelativeAbundance","dominant")])%>%
  full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
            by="Family")%>%
  rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
  pivot_longer(cols = c("NMR","NMRwt"),
               names_to = "class",
               values_to = "MeanRelativeAbundance")%>%
  arrange(-MeanRelativeAbundance)%>%
  # head(n=30)%>%
  mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
  replace_na(replace =list( MeanRelativeAbundance=0))%>%
  ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
             y=MeanRelativeAbundance,
             fill=animal))+
  geom_bar(stat = "identity",position = position_dodge2())+
  theme_bw()+
  labs(x="",
       y="Average relative abundance (%)")+
  labs(fill="")+
  coord_cartesian(expand = FALSE)+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        legend.text = element_text(size = 20))
# ggsave(file.path("./images/barplots",
#                  paste(paste(format(Sys.time(),format="%Y%m%d"),
#                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                        "barplot-wild-vs-lab-nmr-from-paper.png")),
#        plot=last_plot(),
#        width = 5000,height = 4000,
#        units = "px",dpi=300,device = "png")  


# 19. Compare dominant phyla, and genera in NMR vs mice ####
ps.q.agg.dominant.phyla.nmr_b6mouse<-get_dominant_taxa_in_host(ps.q.agg.phylum,
                                                               "Phylum",c("NMR","B6mouse"))
ps.q.agg.dominant.families.nmr_b6mouse<-get_dominant_taxa_in_host(ps.q.agg.family,
                                                                  "Family",c("NMR","B6mouse"))



# 21. Check how much Bacteroidaceae are in NMR  ####
bacteroidaceae.nmr<-ps.q.agg.family.relab%>%
  filter(Family=="Bacteroidaceae",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()
# write.table(bacteroidaceae.nmr,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "bacteroidaceae-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

## 21.1 Check the most dominant Bacteroidota families in NMR ####
bacteroidota.nmr<-ps.q.agg.family.relab%>%
  filter(Phylum=="Bacteroidota",class=="NMR")%>%
  group_by(Family)%>%
  distinct(Family,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Family,MeanRelativeAbundance)
# write.table(bacteroidota.nmr,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "bacteroidota-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 22. Check Spirochaetaceae, Spirochaetota, and Treponema in NMR ####
spirochaetaceae.nmr<-ps.q.agg.family.relab%>%
  filter(Family=="Spirochaetaceae",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
            max=max(RelativeAbundance),
            mean=TotalAgglomRank/TotalClass*100,
            sd=sd(RelativeAbundance),
            n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()
spirochaetota.nmr<-ps.q.agg.phylum.relab%>%
  filter(Phylum=="Spirochaetota",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Phylum,min,max,mean,sd,n)%>%
  distinct()
# write.table(spirochaetaceae.nmr,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "spirochaetaceae-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

treponema.nmr<-ps.q.agg.genus.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Genus,min,max,mean,sd,n)%>%
  distinct()
# write.table(treponema.nmr,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "treponema-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

## 22.1 Check the number of ASVs in Treponema from NMR ####
ps.q.agg%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

# 23. Check Mogibacteriaceae (renamed to Anaerovoracaceae) in NMR ####
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
# write.table(mogibacteriaceae_anaerovoracaceae.all,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 24. Analyse sulfur-metabolising bacteria in NMR ####
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
# 24.1 Total Desulfobacterota in NMR
ps.q.agg.phylum.relab%>%
  filter(Phylum=="Desulfobacterota",class=="NMR")%>%
  pull(MeanRelativeAbundance)%>%
  unique
# write.table(desulfobacterota.nmr,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "desulfobacterota-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")


## Plot Desulfobacterota ####
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

# for (image_format in c("png","tiff")){
#   ggsave(paste0("./images/taxaboxplots/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "Desulfobacterota phylum members-all-hosts",
#                       sep = "-"),".",image_format),
#          plot=last_plot(),
#          width = 4000,height = 6000,
#          units = "px",dpi=300,device = image_format)
# }




## Setup plots of age groups ####
pretty.agegroup.names<-names(table(custom.md.ages$old_agegroup))
names(pretty.agegroup.names)<-names(table(custom.md.ages$agegroup))
pretty.agegroup.names<-c("agegroup0_10"="Young group",
                         "agegroup10_16"="Old group")
agegroup.levels<-names(pretty.agegroup.names)
gg.labs.name<-"Age group"
gg.title.groups<-"age groups"

set.seed(1)
agegroup.fill<-createPalette(length(agegroup.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900",
                                          "#FFA07A"))
names(agegroup.fill)<-agegroup.levels
swatch(agegroup.fill)

sample.levels<-custom.md.ages%>%
  ungroup()%>%
  select(Sample,agegroup)%>%
  arrange(agegroup)%>%
  distinct()
sample.levels$Sample<-factor(sample.levels$Sample,
                             levels=unique(sample.levels$Sample))

### Groups 1, 2, and 3 ####
group1.genera<-c("Faecalibacterium",
                 "Roseburia",
                 "Coprococcus",
                 "Prevotella")
group2.increased.genera<-c("Eggerthella",
                           "Bilophila",
                           "Desulfovibrio",
                           "Fusobacterium",
                           "Anaerotruncus",
                           "Streptococcus",
                           "Escherichia")
group2.unhealthy.genera<-c("Eggerthella",
                           "Coprobacillus",
                           "Streptococcus",
                           "Bilophila",
                           "Actinomyces",
                           "Desulfovibrio",
                           "Campylobacter",
                           "Veillonella",
                           "Enterococcus")
group2.unhealthy.families<-c("Atopobiaceae",
                             "Enterobacteriaceae")
group3.genera<-c("Akkermansia",
                 "Odoribacter",
                 "Butyricimonas",
                 "Butyrivibrio",
                 "Barnesiella",
                 "Oscillospira")
group3.families<-c("Christensenellaceae")

group1.genera%in%ps.q.agg.genus.relab$Genus
group2.increased.genera%in%ps.q.agg.genus.relab$Genus
group2.unhealthy.genera%in%ps.q.agg.genus.relab$Genus
group2.unhealthy.families%in%ps.q.agg.family.relab$Family
group3.genera%in%ps.q.agg.genus.relab$Genus
group3.families%in%ps.q.agg.family.relab$Family

# Finished up to this line ################################################ 
# The function ggplot.species creates a faceted barplot of abundance for each
# taxon per sample. It uses 
# If we are plotting species from a vector of names, use filter(Species%in%species.to.plot)
# If we are plotting all species from a bigger group like entire genus, use
# filter(get(agglom.rank)==species.to.plot)
ggplot.species<-function(taxa.to.plot,
                         tax.df,
                         tax.rank,
                         pretty.agegroup.names,
                         agegroup.fill,
                         sample.levels){ #TODO: fix in other functions
  # if(tax.rank=="Species"){
  #   ggplot.object<-tax.df%>%
  #     filter(Species%in%taxa.to.plot)
  # }else{
  #   ggplot.object<-tax.df%>%
  #     filter(get(tax.rank)==taxa.to.plot)
  # }
  ggplot.object<-tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot)
  
  # order the plot by species vector
  taxa.to.plot<-gsub("_"," ",taxa.to.plot)
  taxa.to.plot<-paste0("<i>",taxa.to.plot,"</i>")
  
  ggplot.object<-ggplot.object%>%
    mutate(!!tax.rank:=gsub("_"," ",get(tax.rank)),
           !!tax.rank:=paste0("<i>",get(tax.rank),"</i>"),
           !!tax.rank:=factor(get(tax.rank),levels=taxa.to.plot))%>%
    mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
    group_by_at(c("class",tax.rank))%>%
    ggplot(aes(x=Sample,
               y=RelativeAbundance,
               fill=factor(agegroup)))+
    geom_bar(stat="identity")+
    facet_wrap(~get(tax.rank),
               scales = "free",
               ncol = 2)+
    theme_bw()+
    labs(x="",
         y="Relative abundance (%)",
         fill=gg.labs.name)+
    scale_color_manual(breaks = unname(pretty.agegroup.names),
                       labels=unname(pretty.agegroup.names))+
    scale_x_discrete(labels=pretty.agegroup.names,
                     limits=sample.levels$Sample)+
    scale_fill_manual(values = agegroup.fill,
                      labels=pretty.agegroup.names)+
    theme(axis.title.y = element_text(size = 25),
          axis.title = element_text(size = 20),
          axis.text.y = ggtext::element_markdown(size=18),
          axis.text.x = element_text(size=20),
          strip.text.x = ggtext::element_markdown(size=20),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  return(ggplot.object)
}


add_zero_rows<-function(taxa.vector,tax.df,tax.rank){
  tax.df<-tax.df%>%ungroup()
  for (sample.name in unique(tax.df$Sample)) {
    missing_taxa <- setdiff(taxa.vector, tax.df[[tax.rank]][tax.df$Sample == sample.name])
    if (length(missing_taxa) > 0) {
      for (taxa in missing_taxa) {
        if(tax.rank=="Species"){
          new_row <- tibble(Sample = sample.name, 
                            Species= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0,
                            Genus=unique(pull(ps.q.agg.relab[which(ps.q.agg.relab$Species==taxa),"Genus"])))
        }else{
          new_row <- tibble(Sample = sample.name, 
                            !!tax.rank:= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0)
        }
        
        tax.df <- tax.df %>% add_row(.before = nrow(df), !!!new_row)
      }
    }
  }
  # To fill the NA values in the empty columns based on non-empty rows in the Sample column 
  tax.df<- tax.df %>%
    group_by(Sample) %>%
    fill(age, .direction = "down")%>%
    fill(agegroup, .direction = "down")%>%
    fill(old_agegroup, .direction = "down")
  
}

# Mean relative abundance of group 1 genera by age group
# The function show.mean.abundances.and.plot filters the dataframe with 
# absolute abundances (tax.df) to retain ones specified in the taxa.to.plot vector.
# It groups the dataframe by taxonomic rank column (e.g. genera) and age groups,
# keeps columns with taxa, their average relative abundances in the dataframe
# and average relative abundances in each age group.
show.mean.abundances.and.plot<-function(taxa.to.plot,
                                        tax.df,
                                        tax.rank){
  tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    group_by_at(c(tax.rank,"agegroup"))%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance","MeanRelativeAbundanceAgegroup")))%>%
    ungroup()%>%
    distinct()%>%
    print()
  taxa.to.plot<-tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    ungroup()%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance")))%>%
    distinct()%>%
    arrange(-MeanRelativeAbundance)%>%
    pull(tax.rank)
  tax.df<-add_zero_rows(taxa.to.plot,tax.df,tax.rank)
  
  return(ggplot.species(taxa.to.plot,tax.df,tax.rank))
}





# ps.q.agg.genus.relab%>%
#   filter(Genus%in%group1.genera&!is.na(MeanRelativeAbundance))%>%
#   group_by(Genus,agegroup)%>%
#   select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
#   ungroup()%>%
#   distinct()
# ps.q.agg.genus.relab%>%
#   filter(Genus%in%group1.genera&!is.na(MeanRelativeAbundance))%>%
#   ungroup()%>%
#   select(Genus,MeanRelativeAbundance)%>%
#   distinct()%>%
#   arrange(-MeanRelativeAbundance)%>%
#   pull(Genus) -> group1.genera
# 
# ps.q.agg.genus.relab<-add_zero_rows(group1.genera,ps.q.agg.genus.relab,"Genus")
# 
# ggplot.species(group1.genera,ps.q.agg.genus.relab,"Genus")+
#   ggtitle(paste0("Relative abundance of Group 1 members"))


group1.genera
group1.genera.plot<-show.mean.abundances.and.plot(group1.genera,
                                                  ps.q.agg.genus.relab,
                                                  "Genus")
group1.genera.plot<-group1.genera.plot+
  ggtitle(paste0("Relative abundance of Group 1 members"))
group1.genera.plot # 4 plots



# Group 2 increased genera
group2.increased.genera
group2.increased.genera.plot<-show.mean.abundances.and.plot(group2.increased.genera,
                                                            ps.q.agg.genus.relab,
                                                            "Genus")
group2.increased.genera.plot<-group2.increased.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members increased with age"))
group2.increased.genera.plot # 4 plots

# Group 2 unhealthy genera
group2.unhealthy.genera
group2.unhealthy.genera.plot<-show.mean.abundances.and.plot(group2.unhealthy.genera,
                                                            ps.q.agg.genus.relab,
                                                            "Genus")
group2.unhealthy.genera.plot<-group2.unhealthy.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.genera.plot # 7 plots

# Group 2 unhealthy families
group2.unhealthy.families
group2.unhealthy.families.plot<-show.mean.abundances.and.plot(group2.unhealthy.families,
                                                              ps.q.agg.family.relab,
                                                              "Family")
group2.unhealthy.families.plot<-group2.unhealthy.families.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.families.plot # 2 plots



# Group 3 genera
group3.genera
group3.genera.plot<-show.mean.abundances.and.plot(group3.genera,
                                                  ps.q.agg.genus.relab,
                                                  "Genus")
group3.genera.plot<-group3.genera.plot+
  ggtitle(paste0("Relative abundance of Group 3 members"))
group3.genera.plot # 5 plots


# Group 3 genera
group3.families.plot<-show.mean.abundances.and.plot(group3.families,
                                                    ps.q.agg.family.relab,
                                                    "Family")
group3.families.plot<-group3.families.plot+
  ggtitle(paste0("Relative abundance of Group 3 members"))
group3.families.plot # 1 plots

for (image_format in c("png","tiff")){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group1.genera-nmr",
                      sep = "-"),".",image_format),
         plot=group1.genera.plot,
         width = 7000,height = 4000,
         units = "px",dpi=300,device = image_format) # 4 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.increased.genera-nmr",
                      sep = "-"),".",image_format),
         plot=group2.increased.genera.plot,
         width = 7000,height = 4000,
         units = "px",dpi=300,device = image_format) # 4 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.genera-nmr",
                      sep = "-"),".",image_format),
         plot=group2.unhealthy.genera.plot,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format) # 7 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.families-nmr",
                      sep = "-"),".",image_format),
         plot=group2.unhealthy.families.plot,
         width = 7000,height = 2000,
         units = "px",dpi=300,device = image_format)    # 2 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      " group3.genera-nmr",
                      sep = "-"),".",image_format),
         plot= group3.genera.plot ,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format)    # 5 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      " group3.families-nmr",
                      sep = "-"),".",image_format),
         plot= group3.families.plot ,
         width = 4000,height = 2000,
         units = "px",dpi=300,device = image_format)    # 1 plot
  
}


# Do rare genera become more abundant? ####
ps.q.agg.genus.relab.nmr<-ps.q.agg.genus.relab%>%
  filter(class=="NMR")
ps.q.agg.genus.relab.nmr%>%
  filter(Genus%in%c("Bacteroides","Roseburia","Faecalibacterium"))%>%
  # distinct(Genus,.keep_all = T)%>%
  # select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus,scales="free_y")

dominant.genera.agegroup0_10<-ps.q.agg.genus.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,agegroup,row.index)

dominant.genera.agegroup10_16<-ps.q.agg.genus.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,agegroup,row.index)

  
ps.q.agg.genus.nmr.meanrelab<-ps.q.agg.genus.relab.nmr%>%
  distinct(Genus,agegroup,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup,agegroup)%>%
  left_join(dominant.genera.agegroup0_10,by=join_by(Genus,agegroup))%>%
  left_join(dominant.genera.agegroup10_16,by=join_by(Genus,agegroup))%>%
  mutate(row.index.x=ifelse(is.na(row.index.x)&!is.na(row.index.y),row.index.y,row.index.x),
         row.index.y=ifelse(is.na(row.index.y)&!is.na(row.index.x),row.index.x,row.index.y))%>%
  select(-row.index.y)%>%
  rename(row.index=row.index.x)

ps.q.agg.genus.nmr.meanrelab.full<-ps.q.agg.genus.nmr.meanrelab%>%
  group_by(Genus)%>%
  mutate(num.obs=n())%>%
  arrange(num.obs)%>%
  filter(num.obs==1)%>%
  mutate(agegroup=ifelse(agegroup=="agegroup0_10","agegroup10_16","agegroup0_10"),
         MeanRelativeAbundanceAgegroup=0,
         row.index=NA)%>%
  full_join(ps.q.agg.genus.nmr.meanrelab,by=
              join_by(Genus,
                      agegroup,
                      MeanRelativeAbundance,
                      MeanRelativeAbundanceAgegroup,
                      row.index))%>%
  select(-num.obs)%>%
  mutate(row.index=ifelse(is.na(row.index),0,row.index))%>%
  # mutate(row.index = unique(row.index[!is.na(row.index)]))%>%
  ungroup()



ps.q.agg.genus.nmr.meanrelab.wide<-ps.q.agg.genus.nmr.meanrelab.full%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(-MeanRelativeAbundanceAgegroup)%>%
  pivot_wider(names_from = agegroup,
              values_from = row.index,
              values_fill = 0,
              names_prefix = "row.ind_")%>%
  left_join(ps.q.agg.genus.nmr.meanrelab.full)%>%
  select(-row.index)%>%
  pivot_wider(names_from = agegroup,
              values_from = MeanRelativeAbundanceAgegroup,
              values_fill = 0)%>%
  mutate(young_vs_old=agegroup10_16-agegroup0_10)%>%
  mutate(young_vs_old_fold=(agegroup10_16-agegroup0_10)/agegroup0_10)%>%
  # filter(agegroup0_10!=0,agegroup10_16!=0)%>%
  arrange(-young_vs_old_fold)%>%
  
  mutate(ind.dif=row.ind_agegroup0_10-row.ind_agegroup10_16)%>%
  filter(row.ind_agegroup10_16<=100,row.ind_agegroup10_16!=0)%>%
  mutate(agegroup0_10=round(agegroup0_10,digits=4),
         agegroup10_16=round(agegroup10_16,digits = 4),
         MeanRelativeAbundance=round(MeanRelativeAbundance,digits=4),
         young_vs_old=round(young_vs_old,digits=4))
  

ps.q.agg.genus.nmr.meanrelab.wide<-ps.q.agg.genus.relab.nmr%>%
  group_by(Genus)%>%
  mutate(num.samples=n())%>%
  select(Genus,num.samples)%>%
  ungroup()%>%
  distinct(Genus,.keep_all = T)%>%
  right_join(ps.q.agg.genus.nmr.meanrelab.wide)
  
nmr.genera.selected<-ps.q.agg.genus.nmr.meanrelab.wide%>%
  filter(abs(young_vs_old)>1)

ps.q.agg.genus.relab.nmr%>%
  filter(Genus%in%nmr.genera.selected$Genus)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus,scales="free_y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))

ps.q.agg.genus.relab.nmr%>%
  filter(Genus%in%c("[Clostridium]_innocuum_group",
                    "Coriobacteriales_Incertae_Sedis Family",
                    "[Eubacterium]_ventriosum_group",
                    "Paludicola",
                    "Synergistaceae Family",
                    "Clostridia_UCG-014",
                    "Blautia"))%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus,scales="free_y")+
  theme_bw()+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        strip.text = element_text(size=20),
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "bottom")

ggsave(file.path("./images/barplots",
                 paste(paste(format(Sys.time(),format="%Y%m%d"),
                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                       "barplot-outliers.png")),
       plot=last_plot(),
       width = 6500,height = 4500,
       units = "px",dpi=300,device = "png")  

outliers<-ps.q.agg.genus.relab.nmr%>%
  group_by(class)%>% # group by class,
  # compute MeanRelativeAbundance from Abundance 
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Genus"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  
  select(-sex,-birthday,-TotalAgglomRank,-TotalClass,-age)%>%
  mutate(MedianRelativeAbundance=median(RelativeAbundance),
         difmed=MeanRelativeAbundance-MedianRelativeAbundance,
         q25_RA=quantile(RelativeAbundance,prob=c(.25,.5,.75))[1],
         q75_RA=quantile(RelativeAbundance,prob=c(.25,.5,.75))[3])%>%
  filter(RelativeAbundance>q75_RA+(q75_RA-q25_RA)*2)%>%
  mutate(n_outliers=n())%>%
  filter(n_outliers<=3)%>%
  group_by(Genus)%>%
  arrange(-difmed)%>%
  ungroup()


# Fitler and recalculate everything ####
filtered.df<-ps.q.agg.genus.relab.nmr%>%
  anti_join(.,outliers)%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample","Genus"))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  
  group_by(class)%>% # group by class,
  # compute MeanRelativeAbundance from Abundance 
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Genus"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  
  group_by(agegroup)%>% # group by class (animal host),
  mutate(TotalAgegroup=sum(Abundance))%>%
  group_by_at(c("agegroup","Genus"))%>%
  mutate(TotalAgglomRankAge=sum(Abundance))%>%
  mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)%>%
  ungroup()%>%
  select(-TotalAgegroup,-TotalAgglomRankAge)
  

######### After filtering ############
filtered.dominant.genera.agegroup0_10<-filtered.df%>%
  filter(agegroup=="agegroup0_10")%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,agegroup,row.index)

filtered.dominant.genera.agegroup10_16<-filtered.df%>%
  filter(agegroup=="agegroup10_16")%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,agegroup,row.index)


filtered.df.meanrelab<-filtered.df%>%
  distinct(Genus,agegroup,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup,agegroup)%>%
  left_join(filtered.dominant.genera.agegroup0_10,by=join_by(Genus,agegroup))%>%
  left_join(filtered.dominant.genera.agegroup10_16,by=join_by(Genus,agegroup))%>%
  mutate(row.index.x=ifelse(is.na(row.index.x)&!is.na(row.index.y),row.index.y,row.index.x),
         row.index.y=ifelse(is.na(row.index.y)&!is.na(row.index.x),row.index.x,row.index.y))%>%
  select(-row.index.y)%>%
  rename(row.index=row.index.x)

filtered.df.meanrelab.full<-filtered.df.meanrelab%>%
  group_by(Genus)%>%
  mutate(num.obs=n())%>%
  arrange(num.obs)%>%
  filter(num.obs==1)%>%
  mutate(agegroup=ifelse(agegroup=="agegroup0_10","agegroup10_16","agegroup0_10"),
         MeanRelativeAbundanceAgegroup=0,
         row.index=NA)%>%
  full_join(filtered.df.meanrelab,by=
              join_by(Genus,
                      agegroup,
                      MeanRelativeAbundance,
                      MeanRelativeAbundanceAgegroup,
                      row.index))%>%
  select(-num.obs)%>%
  mutate(row.index=ifelse(is.na(row.index),0,row.index))%>%
  # mutate(row.index = unique(row.index[!is.na(row.index)]))%>%
  ungroup()



filtered.df.meanrelab.wide<-filtered.df.meanrelab.full%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(-MeanRelativeAbundanceAgegroup)%>%
  pivot_wider(names_from = agegroup,
              values_from = row.index,
              values_fill = 0,
              names_prefix = "row.ind_")%>%
  left_join(filtered.df.meanrelab.full)%>%
  select(-row.index)%>%
  pivot_wider(names_from = agegroup,
              values_from = MeanRelativeAbundanceAgegroup,
              values_fill = 0)%>%
  mutate(young_vs_old=agegroup10_16-agegroup0_10)%>%
  mutate(young_vs_old_fold=(agegroup10_16-agegroup0_10)/agegroup0_10)%>%
  # filter(agegroup0_10!=0,agegroup10_16!=0)%>%
  arrange(-young_vs_old_fold)%>%
  
  mutate(ind.dif=row.ind_agegroup0_10-row.ind_agegroup10_16)%>%
  filter(row.ind_agegroup10_16<=100,row.ind_agegroup10_16!=0)%>%
  mutate(agegroup0_10=round(agegroup0_10,digits=4),
         agegroup10_16=round(agegroup10_16,digits = 4),
         MeanRelativeAbundance=round(MeanRelativeAbundance,digits=4),
         young_vs_old=round(young_vs_old,digits=4))


filtered.df.meanrelab.wide<-filtered.df%>%
  group_by(Genus)%>%
  mutate(num.samples=n())%>%
  select(Genus,num.samples)%>%
  ungroup()%>%
  distinct(Genus,.keep_all = T)%>%
  right_join(filtered.df.meanrelab.wide)

filtered.genera.selected<-filtered.df.meanrelab.wide%>%
  filter(abs(young_vs_old)>1)

filtered.df%>%
  filter(Genus%in%filtered.genera.selected$Genus)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus,scales="free_y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))

filtered.df%>%
  filter(Genus%in%c("[Clostridium]_innocuum_group",
                    "Coriobacteriales_Incertae_Sedis Family",
                    "[Eubacterium]_ventriosum_group",
                    "Paludicola",
                    "Synergistaceae Family",
                    "Clostridia_UCG-014",
                    "Blautia"))%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus,scales="free_y")+
  theme_bw()



# filtered.df%>%
#   arrange(-MeanRelativeAbundanceAgegroup)%>%
#   select(-MeanRelativeAbundanceAgegroup)%>%
#   pivot_wider(names_from = agegroup,
#               values_from = row.index,
#               values_fill = 0,
#               names_prefix = "row.ind_")%>%
#   left_join(ps.q.agg.genus.nmr.meanrelab.full)%>%
#   select(-row.index)%>%
#   pivot_wider(names_from = agegroup,
#               values_from = MeanRelativeAbundanceAgegroup,
#               values_fill = 0)%>%
#   mutate(young_vs_old=agegroup10_16-agegroup0_10)%>%
#   mutate(young_vs_old_fold=(agegroup10_16-agegroup0_10)/agegroup0_10)%>%
#   # filter(agegroup0_10!=0,agegroup10_16!=0)%>%
#   arrange(-young_vs_old_fold)%>%
#   
#   mutate(ind.dif=row.ind_agegroup0_10-row.ind_agegroup10_16)%>%
#   filter(row.ind_agegroup10_16<=100,row.ind_agegroup10_16!=0)%>%
#   mutate(agegroup0_10=round(agegroup0_10,digits=4),
#          agegroup10_16=round(agegroup10_16,digits = 4),
#          MeanRelativeAbundance=round(MeanRelativeAbundance,digits=4),
#          young_vs_old=round(young_vs_old,digits=4))
#   
#   
#   
#   
#   
#   mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
#   ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
#   geom_bar(stat="identity")+
#   facet_wrap(~Genus,scales="free_y")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=45,hjust=1))
