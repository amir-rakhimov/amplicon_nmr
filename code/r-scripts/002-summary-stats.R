# In this script, we will explore the imported dataset from QIIME2 (using qiime2R
# and phyloseq).
# We will also rarefy the data for future analyses.
# We will use the data from 001-phyloseq-qiime2.R script (ps.q.agg
# agglomerated tables at phylum, family, genus, and OTU level).

## 1. Import libraries ####
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
# install.packages(c("tidyverse","vegan","ggtext","Polychrome","ggrepel"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(tidyverse)
library(vegan)
library(Polychrome)
library(ggtext)
library(ggrepel)

## 2. Specifying parameters and directory/file names #### 
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)

truncationlvl<-"234" # truncation level that we chose in QIIME2
# truncationlvl<-"0" # truncation level that we chose in QIIME2

authorname<-"pooled" # name of the folder with QIIME2 output
# authorname<-"biagi" # name of the folder with QIIME2 output

read.end.type<-"single" # single reads or paired reads: decided in QIIME2
# Import datasets as rds files
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))


barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")

## Setup plots ####
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "DMR" = "*Fukomys Damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*")
# Use only the taxa that are present in the workspace
# (custom.md is metadata from the rdafile)
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
                           seedcolors = c("#EE2C2C","#BF3EFF","#5CACEE","#00CD66",
                                          "#FF8C00","#00EE00","#EEC900", "#00FFFF",
                                          "#FF6EB4",
                                          "#FFA07A"))
names(custom.fill)<-custom.levels
swatch(custom.fill)


# Setup general ggplot theme
mytheme<-theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
               axis.text.y = element_text(size=20), # size of y axis ticks
               axis.title = element_text(size = 20), # size of axis names
               plot.title = element_text(size = 25), # size of plot title
               plot.caption = element_text(size=23), # size of plot caption
               legend.text = element_text(size = 20), # size of legend text
               legend.title = element_text(size = 25) # size of legend title
)
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
#             file=file.path(rtables.directory,authorname,
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
ps.q.agg.relab<-add_relab_to_tax_df(ps.q.agg,"OTU")

# 9. Add agegroup (Must run for plotting) ####
custom.md.ages<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))
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
    left_join(metadata.df[,c("Sample","agegroup","old_agegroup")],by="Sample")%>%
    group_by(agegroup)%>% # group by class (animal host),
    mutate(TotalAgegroup=sum(Abundance))%>%
    group_by_at(c("agegroup",tax.rank))%>%
    mutate(TotalAgglomRankAge=sum(Abundance))%>%
    mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)%>%
    ungroup()%>%
    select(-TotalAgegroup,-TotalAgglomRankAge)
}


ps.q.agg.family.relab.nmr<-ps.q.agg.family.relab%>%
  filter(class=="NMR")
ps.q.agg.genus.relab.nmr<-ps.q.agg.genus.relab%>%
  filter(class=="NMR")
ps.q.agg.relab.nmr<-ps.q.agg.relab%>%
  filter(class=="NMR")

ps.q.agg.family.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.family.relab.nmr,"Family",
                                              custom.md.ages)
ps.q.agg.genus.relab<-add_agegroup_to_tax_df(ps.q.agg.genus.relab,"Genus",
                                             custom.md.ages)
ps.q.agg.relab<-add_agegroup_to_tax_df(ps.q.agg.relab,"OTU",
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
### Get the number of unclassified genera in each host ####
unclassified.genus.summary.stats.table<-ps.q.agg.genus.relab%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(class)%>%
  distinct(Genus,.keep_all = T)%>%
  tally()%>%
  arrange(-n)%>%
  rename(NumUnclassifiedGenera=n)%>%
  left_join(unclassified.genus.summary.stats.table)
### Sanity check: get the number of classified genera in each host (unrarefied) ####
unclassified.genus.summary.stats.table<-ps.q.agg.genus.relab%>%
  # filter(class=="NMR")%>%
  filter(!grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(class)%>%
  distinct(Genus,.keep_all = T)%>%
  tally()%>%
  arrange(-n)%>%
  rename(NumClassifiedGenera=n)%>%
  left_join(unclassified.genus.summary.stats.table)

# write.table(unclassified.genus.summary.stats.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 11. Rarefy the table and check the percentage of unclassified taxa ####
# Convert the data frame into wide format: rows are samples and columns
# are taxa
get_rarefied_table<-function(tax.df,tax.rank,host.classes){
  tax.df.wide<-tax.df%>%
    filter(class %in% host.classes,Abundance!=0)%>%
    select(Sample,Abundance,class,all_of(tax.rank))%>%
    filter(Abundance!=0)%>%
    pivot_wider(names_from = all_of(tax.rank),
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("Sample")%>% # Set sample names as row names
    select(-class)
  # Find the smallest sample size
  min.n_seqs.all<-tax.df%>%
    filter(class %in% host.classes,Abundance!=0)%>%
    select(Sample,all_of(tax.rank),Abundance)%>%
    group_by(Sample)%>%
    summarize(n_seqs=sum(Abundance))%>%
    summarize(min=min(n_seqs))%>%
    pull(min)
  print(paste("Smallest sample size:", min.n_seqs.all))
  
  ### Rarefied asv table with vegan ####
  set.seed(1)
  tax.df.rare<-rrarefy(tax.df.wide,sample=min.n_seqs.all)
  tax.df.rare<-tax.df.rare%>%
    as_tibble(rownames="Sample")%>%
    pivot_longer(-Sample)%>%
    as.data.frame()%>%
    left_join(unique(tax.df[,c("Sample","class")]),
              by="Sample")
  if(tax.rank=="OTU"){
    tax.df.rare<-tax.df.rare%>%
      rename(OTU=name,Abundance=value)%>%
      filter(Abundance!=0)  
  }else{
    # rename the 'name' column corresponding to the tax.rank
    tax.df.rare[,paste(tax.rank)]<-tax.df.rare$name
    tax.df.rare<-tax.df.rare%>%
      select(-name)%>%
      rename(Abundance=value)%>%
      filter(Abundance!=0)
  }
  # write.table(tax.df.rare,
  #             file = file.path(rtables.directory,paste0(
  #               paste(
  #                 paste(format(Sys.time(),format="%Y%m%d"),
  #                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #                 "ps.q.df.rare-nonfiltered",tax.rank,
  #                 paste(host.classes,collapse = '-'),sep = "-"),
  #               ".tsv")),
  #             row.names = F,
  #             sep = "\t")
  # saveRDS(tax.df.rare,
  #         file = file.path(rdafiles.directory,paste0(
  #           paste(
  #             paste(format(Sys.time(),format="%Y%m%d"),
  #                   format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #             "ps.q.df.rare-nonfiltered",tax.rank,
  #             paste(host.classes,collapse = '-'),sep = "-"),
  #           ".rds")))
}

ps.q.agg.genus.rare<-get_rarefied_table(ps.q.agg.genus,"Genus",custom.levels)
ps.q.agg.genus.nmr.rare<-get_rarefied_table(ps.q.agg.genus.relab.nmr,"Genus","NMR")
ps.q.agg.nmr.rare<-get_rarefied_table(ps.q.agg.relab.nmr,"OTU","NMR")

### Add relative abundances and taxonomic information to the rarefied dataframe ####
# All hosts (genus)
ps.q.agg.genus.rare.relab<-add_relab_to_tax_df(ps.q.agg.genus.rare,"Genus")
# Add other taxonomic ranks to the dataframe
ps.q.agg.genus.rare.relab<-ps.q.agg.genus.rare.relab%>%
  left_join(unique(ps.q.agg.genus[,c("Kingdom","Phylum","Class","Order","Family","Genus")]))

# Add relative abundances to NMR dataframe (ASV level)
ps.q.agg.nmr.rare.relab<-add_relab_to_tax_df(ps.q.agg.nmr.rare,"OTU")
ps.q.agg.nmr.rare.relab<-ps.q.agg.nmr.rare.relab%>%
  left_join(unique(ps.q.agg[,c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")]))

### 11.1 Plot a rarefaction curve ####
# ps.q.mat<-as(t(otu_table(ps.q)),"matrix") # from phyloseq
# ps.q.genus.mat<-ps.q.agg.genus%>%
#   filter(class %in% custom.levels,Abundance!=0)%>%
#   select(Sample,Abundance,class,all_of(agglom.rank))%>%
#   filter(Abundance!=0)%>%
#   pivot_wider(names_from = all_of(agglom.rank),
#               values_from = "Abundance",
#               values_fill = 0)%>%
#   as.data.frame()%>%
#   column_to_rownames("Sample")%>% # Set sample names as row names
#   select(-class)%>%
#   as.matrix() # convert to matrix
# set.seed(1)
# rare.df<-rarecurve(ps.q.genus.mat,step = 100,sample=min(rowSums(ps.q.genus.mat)),tidy = TRUE)
# rare.df%>%
#   # filter(Sample<=100000)%>%
#   group_by(Site)%>% # site is sample name
#   mutate(label=if_else(Sample==max(Sample),as.character(Site),NA_character_))%>%
#   # filter(Site%in%rownames(custom.md[which(custom.md$class=="NMR"),]))%>%
#   filter(Site%in%unique(ps.q.agg.genus$Sample[ps.q.agg.genus$class%in%custom.levels]))%>%
#   ggplot(.,aes(x=Sample,y=Species,col=Site))+
#   geom_line()+
#   # coord_cartesian(xlim=c(0,100000))+
#   geom_vline(xintercept = min(rowSums(ps.q.genus.mat)))+
#   annotate("text",
#            x=min(rowSums(ps.q.genus.mat))+2000,
#            y=10,
#            label=min(rowSums(ps.q.genus.mat)))+
#   geom_label_repel(aes(label = label),
#                    nudge_x = 1,
#                    na.rm = TRUE) +
#   theme_bw()+
#   labs(x="Sample size",
#        y="ASV")+
#   theme(legend.position = "none")
# ggsave(paste0("./images/lineplots/",
#               paste(paste(format(Sys.time(),format="%Y%m%d"),
#                           format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                     "rarecurve",
#                     truncationlvl,agglom.rank,
#                     sep = "-"),".png"),
#        plot=last_plot(),
#        width = 4500,height = 3000,
#        units = "px",dpi=300,device = "png")

### 11.3 Calculate summary stats of unclassified taxa for rarefied data ####
unclassified.genus.summary.stats.table.rare<-ps.q.agg.genus.rare.relab%>%
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
### 11.4 Get the number of unclassified genera for rarefied data ####
unclassified.genus.summary.stats.table.rare<-ps.q.agg.genus.rare.relab%>%
  filter(grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                     collapse='|'),get(agglom.rank)))%>%
  group_by(class)%>%
  distinct(Genus,.keep_all = T)%>%
  tally()%>%
  arrange(-n)%>%
  rename(NumUnclassifiedGenera=n)%>%
  left_join(unclassified.genus.summary.stats.table.rare)

### Sanity check: get the number of classified genera in rarefied data ####
unclassified.genus.summary.stats.table.rare<-ps.q.agg.genus.rare.relab%>%
  # filter(class=="NMR")%>%
  filter(!grepl(paste(all.ranks[! all.ranks %in% agglom.rank],
                      collapse='|'),get(agglom.rank)))%>%
  group_by(class)%>%
  distinct(Genus,.keep_all = T)%>%
  tally()%>%
  arrange(-n)%>%
  rename(NumClassifiedGenera=n)%>%
  left_join(unclassified.genus.summary.stats.table.rare)

# write.table(unclassified.genus.summary.stats.table.rare,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table-rarefied.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 11. Bar plot of the total number of unique genera per host ####
summary.stats.table%>%
  ggplot(aes(x=reorder(class,-GeneraPerHost),y=GeneraPerHost))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x="",
       y="Genera per host")+
  scale_x_discrete(labels =pretty.level.names) + # x axis labels become 
  # pretty.level.names because the names in the vector are the same as the 
  # x axis ticks in the original plot
  mytheme + # general theme
  theme(panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
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
  mutate(class=factor(class,levels=custom.levels))%>%
  group_by(Sample,class,Genus)%>% # group by class (animal host),
  distinct(Genus)%>%
  group_by(class,Sample)%>%
  tally%>%
  ggplot(aes(x=class,y=n,fill=class))+
  geom_boxplot(alpha=0.1)+ # colorless boxplot
  geom_jitter(aes(color=class),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=2)+
  theme_bw()

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
  ggplot(aes(x=factor(class,levels=mean.factors),y=TotalUnclassifiedPercent,
             fill=class))+
  geom_boxplot(alpha=0.1,show.legend = F)+
  stat_summary(fun=mean, 
               geom='point', 
               shape=23,
               fill="red",
               size=4,
               alpha=0.5)+ # add a point that shows the mean value in each box
  geom_jitter(aes(color=class),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=2,
              show.legend = FALSE)+
  scale_color_manual(values=custom.fill)+
  scale_fill_manual(values=custom.fill)+
  theme_bw()+
  labs(x="",
       y="Unclassified genera per sample (%)")+
  scale_x_discrete(labels =pretty.level.names) +
  mytheme + 
  theme(axis.text.x = element_markdown(angle=45,size=20,hjust=1) # rotate 
        # the x-axis labels by 45 degrees and shift to the right
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
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-phyla-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# write.table(ps.q.agg.dominant.families,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-families-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# write.table(ps.q.agg.dominant.genera,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-genera-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# # 13. Analyse data from Debebe ####
# ### 13.1 Load data from Debebe ####
# truncationlvl<-"0" # truncation level that we chose in QIIME2
# authorname<-"biagi" # name of the folder with QIIME2 output
# read.end.type<-"single" # single reads or paired reads: decided in QIIME2
# # If you already created a workspace, just load it
# phyloseq.biagi.date_time<-"20240502_23_43_39"
# load(file.path(rdafiles.directory,paste(
#   phyloseq.biagi.date_time,
#   authorname,read.end.type,"phyloseq-summary-stats",
#   truncationlvl,
#   "workspace.RData",sep = "-")))
# 
# ### 13.1 Total number of unique ASV/phyla/families/genera per class in Debebe ####
# n.asv.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.biagi,"OTU")
# n.phylum.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.phylum.biagi,"Phylum")
# n.family.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.family.biagi,"Family")
# n.genus.per.host.biagi<-get.n.uniq.taxa.per.host(ps.q.agg.genus.biagi,"Genus")
# 
# ### 13.2 Summary statistics table from Debebe ####
# summary.stats.table.biagi<-create_summary_stats_table(ps.q.agg.genus.biagi,
#                                                       n.asv.per.host.biagi,
#                                                       n.phylum.per.host.biagi,
#                                                       n.family.per.host.biagi,
#                                                       n.genus.per.host.biagi)
# 
# # write.table(summary.stats.table.biagi,
# #             file=file.path("./output/rtables","biagi",
# #                            paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                                  "summary-table.tsv",sep="-")),
# #             row.names = F,sep = "\t")
# 
# ### 13.3 Most abundant phyla, families, and genera in Debebe ####
# ps.q.agg.dominant.phyla.biagi<-get_dominant_taxa_in_host(ps.q.agg.phylum.biagi,
#                                                        "Phylum","NMRwt")
# ps.q.agg.dominant.families.biagi<-get_dominant_taxa_in_host(ps.q.agg.family.biagi,
#                                                           "Family","NMRwt")
# ps.q.agg.dominant.genera.nmr<-get_dominant_taxa_in_host(ps.q.agg.genus,
#                                                         "Genus","NMR")
# # save.image(paste0(rdafiles.directory,
# #                   paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                         "pooled-and-biagi-workspace.Rda")))
# 
# 
# 
# 
# ### 13.4 Bar plot of dominant families in lab vs wild NMR from our data (wild data reanalysed) ####
# ps.q.agg.dominant.families.nmr%>%
#   ungroup()%>%
#   select(Family,MeanRelativeAbundance)%>%
#   full_join(ps.q.agg.dominant.families.biagi[,c("Family","MeanRelativeAbundance")],
#             by="Family")%>%
#   rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
#   pivot_longer(cols = c("NMR","NMRwt"),
#                names_to = "class",
#                values_to = "MeanRelativeAbundance")%>%
#   arrange(-MeanRelativeAbundance)%>%
#   head(n=30)%>%
#   mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
#   ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
#              y=MeanRelativeAbundance,
#              fill=animal))+
#   geom_bar(stat = "identity",position = "dodge")+
#   theme_bw()+
#   xlab("") +
#   ylab("Average relative abundance (%)")+
#   labs(fill="")+
#   coord_cartesian(expand = FALSE)+
#   mytheme + # general ggplot theme
#   theme(axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate
#         # the x-axis labels by 45 degrees and shift to the right
#         legend.position = "inside",
#         legend.position.inside =  c(0.9,0.8)) # legend under the plot
# # ggsave(file.path("./images/barplots",
# #                  paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                        "barplot-wild-vs-lab-nmr-my-results.png")),
# #        plot=last_plot(),
# #        width = 5000,height = 4000,
# #        units = "px",dpi=300,device = "png")  
# 
# 
# 
# #### 13.4.1 Bar plot of dominant families reported in the paper ####
# wild.nmr.families<-data.frame(
#   Family=c("Lachnospiraceae",
#            "Prevotellaceae",
#            "Paraprevotellaceae",
#            "Bacteroidales Order",
#            "Clostridiales Order",
#            "Oscillospiraceae",
#            "Veillonellaceae",
#            "Clostridiaceae",
#            "Muribaculaceae",
#            "Porphyromonadaceae"),
#   MeanRelativeAbundance=c(17.6,
#                           11,
#                           8.8,
#                           6.2,
#                           6.1,
#                           5.7,
#                           4.7,
#                           4.1,
#                           4,
#                           3)
# )
# # Add a column if the family was reported as dominant or not
# ps.q.agg.dominant.families.nmr.with_domin<-ps.q.agg.dominant.families.nmr%>%
#   ungroup()%>%
#   mutate(row.index=as.numeric(rownames(ps.q.agg.dominant.families.nmr)),
#          dominant=ifelse(row.index<=10,"dominant","not_dominant"))
# 
# ps.q.agg.dominant.families.nmr.with_domin%>%
#   select(Family,MeanRelativeAbundance,dominant)%>%
#   distinct(Family,.keep_all = T)%>%
#   filter(Family%in%wild.nmr.families$Family)%>%
#   full_join(ps.q.agg.dominant.families.nmr.with_domin[1:10,c("Family","MeanRelativeAbundance","dominant")])%>%
#   full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
#             by="Family")%>%
#   rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
#   pivot_longer(cols = c("NMR","NMRwt"),
#                names_to = "class",
#                values_to = "MeanRelativeAbundance")%>%
#   arrange(-MeanRelativeAbundance)%>%
#   # head(n=30)%>%
#   mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
#   replace_na(replace =list( MeanRelativeAbundance=0))%>%
#   ggplot(aes(x=reorder(Family,desc(MeanRelativeAbundance)),
#              y=MeanRelativeAbundance,
#              fill=animal))+
#   geom_bar(stat = "identity",position = position_dodge2())+
#   theme_bw()+
#   labs(x="",
#        y="Average relative abundance (%)")+
#   labs(fill="")+
#   coord_cartesian(expand = FALSE)+
#   mytheme + # general ggplot theme
#   theme(axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate
#         # the x-axis labels by 45 degrees and shift to the right
#         legend.position = "inside",
#         legend.position.inside = c(0.9,0.8))
# # ggsave(file.path("./images/barplots",
# #                  paste(paste(format(Sys.time(),format="%Y%m%d"),
# #                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
# #                        "barplot-wild-vs-lab-nmr-from-paper.png")),
# #        plot=last_plot(),
# #        width = 5000,height = 4000,
# #        units = "px",dpi=300,device = "png")  


# 14. Compare dominant phyla, and genera in NMR vs mice ####
ps.q.agg.dominant.phyla.nmr_b6mouse<-get_dominant_taxa_in_host(ps.q.agg.phylum,
                                                               "Phylum",c("NMR","B6mouse"))
ps.q.agg.dominant.families.nmr_b6mouse<-get_dominant_taxa_in_host(ps.q.agg.family,
                                                                  "Family",c("NMR","B6mouse"))

# 15. Check how much Bacteroidaceae are in NMR  ####
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
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "bacteroidaceae-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

### 15.1 Check the most dominant Bacteroidota families in NMR ####
bacteroidota.nmr<-ps.q.agg.family.relab%>%
  filter(Phylum=="Bacteroidota",class=="NMR")%>%
  group_by(Family)%>%
  distinct(Family,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Family,MeanRelativeAbundance)
# write.table(bacteroidota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "bacteroidota-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 16. Check Spirochaetaceae, Spirochaetota, and Treponema in NMR ####
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
#             file=file.path(rtables.directory,
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
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "treponema-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

### 16.1 Check the number of ASVs in Treponema from NMR ####
ps.q.agg%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

# 17. Check Mogibacteriaceae (renamed to Anaerovoracaceae) in NMR ####
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
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 18. Analyse sulfur-metabolising bacteria in NMR ####
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
# Desulfobacterota in other hosts
desulfobacterota.all.mean<-ps.q.agg.genus.relab%>%
  filter(Phylum=="Desulfobacterota")%>%
  group_by(Genus,class)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Genus,MeanRelativeAbundance,class)%>%
  distinct()%>%
  pivot_wider(names_from = Genus,
              values_from = MeanRelativeAbundance,
              values_fill = 0)
colnames(desulfobacterota.all.mean)[which(colnames(desulfobacterota.all.mean)!="class")]<-
  paste0("Mean",colnames(desulfobacterota.all.mean)[which(colnames(desulfobacterota.all.mean)!="class")])
desulfobacterota.all.sd<-ps.q.agg.genus.relab%>%
  filter(Phylum=="Desulfobacterota")%>%
  group_by(Genus,class)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Genus,sd,class)%>%
  distinct()%>%
  pivot_wider(names_from = Genus,
              values_from = sd,
              values_fill = 0)
colnames(desulfobacterota.all.sd)[which(colnames(desulfobacterota.all.sd)!="class")]<-
  paste0("SD",colnames(desulfobacterota.all.sd)[which(colnames(desulfobacterota.all.sd)!="class")])
desulfobacterota.all.mean%>%left_join(desulfobacterota.all.sd)%>%View
### 18.1 Total Desulfobacterota in NMR ####
ps.q.agg.phylum.relab%>%
  filter(Phylum=="Desulfobacterota",class=="NMR")%>%
  pull(MeanRelativeAbundance)%>%
  unique
# write.table(desulfobacterota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "desulfobacterota-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(desulfobacterota.all.mean,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "desulfobacterota-table-all-mean.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(desulfobacterota.all.mean,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "desulfobacterota-table-all-sd.tsv",sep="-")),
#             row.names = F,sep = "\t")

### 18.2 Plot Desulfobacterota ####
# It's a flipped plot
ps.q.agg.genus.relab%>%
  mutate(class=factor(class,levels=custom.levels))%>%
  filter(Phylum=="Desulfobacterota")%>%
  group_by(class,Genus)%>%
  ggplot(aes(x=factor(class,levels=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE,
               outliers = FALSE)+
  geom_jitter(width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=3,
              show.legend = FALSE)+
  facet_wrap(~Genus,scales = "free_x",
             ncol = 2)+
  theme_bw()+
  coord_flip()+
  labs(x="",
       y="Relative abundance (%)")+
  # VVV scale_x_discrete works on x axis but the final labels are on y because the coord_flip() flipped the plot.
  # So, we use rev() to address the flipping VVV
  scale_x_discrete(labels=rev(pretty.level.names), # new labels (named vector) on the axis
                   limits=rev(custom.levels) # limits adjust which levels (and in what order) are displayed
                   )+ # rename boxplot labels 
  scale_fill_manual(values = rev(custom.fill))+ # use custom color for boxplots (values matches data to the named vector)
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




## 19. Setup sample levels for NMR####
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))


## 20. Groups 1, 2, and 3 ####
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

group1.genera%in%ps.q.agg.genus.relab.nmr$Genus
group2.increased.genera%in%ps.q.agg.genus.relab.nmr$Genus
group2.unhealthy.genera%in%ps.q.agg.genus.relab.nmr$Genus
group2.unhealthy.families%in%ps.q.agg.family.relab.nmr$Family
group3.genera%in%ps.q.agg.genus.relab.nmr$Genus
group3.families%in%ps.q.agg.family.relab.nmr$Family

# The function ggplot.species creates a faceted barplot of abundance for each
# taxon per sample (only NMR). It uses a vector with names of taxa, a dataframe with
# abundances, and a taxonomic rank (e.g. Genus). It creates age groups for the 
# barplot and creates levels for the color palette.
# The function keeps only taxa that are found in the taxa.to.plot. vector using
# filter(get(tax.rank)%in%taxa.to.plot). It doesn't matter if we plot only
# one taxon or multiple. Taxonomic rank also doesn't matter.
# The barplot is faceted by taxon and bars are colored by age group. The
# palette is created with Polychrome library. The function returns a ggplot object.
ggplot.species<-function(taxa.to.plot,
                         tax.df,
                         tax.rank){ 
  gg.labs.name<-"Age group"
  pretty.agegroup.names<-c("Young","Old")
  # add names to age groups to match the dataframe
  names(pretty.agegroup.names)<-tax.df%>%
    ungroup()%>%
    select(agegroup)%>%
    distinct(agegroup)%>%
    arrange(agegroup)%>%
    pull
  agegroup.levels<-names(pretty.agegroup.names)
  set.seed(1)
  agegroup.fill<-createPalette(length(agegroup.levels),
                               seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                              "#FF8C00","#BF3EFF", "#00FFFF",
                                              "#FF6EB4","#00EE00","#EEC900",
                                              "#FFA07A"))
  names(agegroup.fill)<-agegroup.levels
  
  ggplot.object<-tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot)
  
  # order the plot by species vector
  taxa.to.plot<-gsub("_"," ",taxa.to.plot)
  taxa.to.plot<-paste0("<i>",taxa.to.plot,"</i>")
  # create a dataframe of samples in each age group
  sample.levels<-tax.df%>%
    ungroup()%>%
    filter(class=="NMR")%>%
    select(Sample,age)%>%
    arrange(age,Sample)%>%
    distinct()
  # convert samples to factors
  sample.levels$Sample<-factor(sample.levels$Sample,
                               levels=unique(sample.levels$Sample))
  # !!tax.rank:= taxa (bang bang) evaluates the tax.rank before the rest is evaluated
  # to substitute an environment-variable (created with <-) with a data-variable (inside a data frame).
  # https://rlang.r-lib.org/reference/topic-inject.html
  ggplot.object<-ggplot.object%>%
    mutate(!!tax.rank:=gsub("_"," ",get(tax.rank)),# remove underscores
           !!tax.rank:=paste0("<i>",get(tax.rank),"</i>"), # convert into italic
           !!tax.rank:=factor(get(tax.rank),levels=taxa.to.plot))%>% # taxa as factors
    mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>% # samples as factors
    group_by_at(c("class",tax.rank))%>%
    ggplot(aes(x=Sample,
               y=RelativeAbundance,
               fill=factor(agegroup)))+
    geom_bar(stat="identity")+
    facet_wrap(~get(tax.rank), # faceted barplot
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

# The function add_zero_rows adds rows with Abundance = 0 to a dataframe tax.df
# based on a vector taxa.vector in the taxonomic rank tax.rank. When we create
# a barplot using the ggplot.species function, we want to show all bars, even 
# those with 0 value. So, we use the function add_zero_rows for adding empty bars.
# 1. The function loops through each sample in the dataframe: for (sample.name in unique(tax.df$Sample)).
# 2. The function creates a missing_taxa vector: setdiff(taxa.vector, tax.df[[tax.rank]][tax.df$Sample == sample.name]).
# The rationale is this: find which rows belong to a sample sample.name: tax.df$Sample == sample.name.
# The result is a vector of TRUE/FALSE values.
# Then, extract all rows with taxa at the tax.rank column as a vector: tax.df[[tax.rank]].
# Then, keep only those taxa that are found in the sample.name sample (if
# tax.df$Sample == sample.name is TRUE, then we keep tax.df[[tax.rank]] ).
# Finally, the function checks if the taxa in the taxa.vector are found in the filtered vector.
# If some taxon isn't found, the function will add a row with 0 Abundance. This
# taxon will be in the missing_taxa vector
# 3. Once the function creates the missing_taxa vector, it checks its length.
# If the length >0, the function will proceed to creation of zero rows: if (length(missing_taxa) > 0) 
# 4. While looping through each taxon in the missing_taxa vector, the function
# creates a new row as a tibble with 5 columns: Sample is the sample.name that 
# it's looping through in step 1. 
# The Abundance and RelativeAbundance columns will have a 0 value.
# If the tax.rank is Species, the next column is "Species", and the value is 
# the missing taxon. Otherwise, the function creates a column according to tax.rank.
# If the tax.rank is Species, the Genus column will have a value from the tax.df
# (tax.df is filtered by Species column and Genus is kep, then pooled): 
# unique(pull(tax.df[which(tax.df$Species==taxa),"Genus"])))
# 5. The new_row is added to the tax.df using the add_row function at the end 
# of the dataframe. !!! is a splicing operator that injects a list of 
# arguments: tax.df %>% add_row(.before = nrow(df), !!!new_row)
# https://rlang.r-lib.org/reference/splice-operator.html
# https://rlang.r-lib.org/reference/topic-inject.html
# !!! is not the same as !!
# 6. The columns that have NA value after injection will be filled by the 
# corresponding value from the same sample. So, if "age" column is empty in the
# 2D10 sample, dplyr will take the value from the other row that is not empty.
add_zero_rows<-function(taxa.vector,tax.df,tax.rank){
  tax.df<-tax.df%>%ungroup() # just in case
  for (sample.name in unique(tax.df$Sample)) {
    missing_taxa <- setdiff(taxa.vector, tax.df[[tax.rank]][tax.df$Sample == sample.name])
    if (length(missing_taxa) > 0) {
      for (taxa in missing_taxa) {
        if(tax.rank=="Species"){
          new_row <- tibble(Sample = sample.name, 
                            Species= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0,
                            Genus=unique(pull(tax.df[which(tax.df$Species==taxa),"Genus"])))
        }else{
          new_row <- tibble(Sample = sample.name, 
                            !!tax.rank:= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0)
        }
        # !!! injects a list of arguments
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
# 1. The function show.mean.abundances.and.plot filters the dataframe with 
# absolute abundances (tax.df) to retain ones specified in the taxa.to.plot vector.
# It groups the dataframe by taxonomic rank column (e.g. genera) and age groups,
# keeps columns with taxa, their average relative abundances in the dataframe
# and average relative abundances in each age group. The output is printed into
# the console (agegroup, Genus, MeanRelativeAbundance, MeanRelativeAbundanceAgegroup). 
# 2. The function sorts the vector of taxa to be plotted according to their 
# average relative abundance: arrange(-MeanRelativeAbundance).
# 3. The function calls add_zero_rows function to add empty rows to the 
# tax.df before plotting.
# 4. The function calls ggplot.species function and returns the ggplot object
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
  sorted.taxa.to.plot<-tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    ungroup()%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance")))%>%
    distinct()%>%
    arrange(-MeanRelativeAbundance)%>%
    pull(tax.rank)
  tax.df<-add_zero_rows(sorted.taxa.to.plot,tax.df,tax.rank)
  
  return(ggplot.species(sorted.taxa.to.plot,tax.df,tax.rank))
}

# Mean relative abundance of group 1 genera by age group 
# ggplot.species(group1.genera,ps.q.agg.genus.relab.nmr,"Genus")+
#   ggtitle(paste0("Relative abundance of Group 1 members"))


group1.genera
group1.genera.plot<-show.mean.abundances.and.plot(group1.genera,
                                                  ps.q.agg.genus.relab.nmr,
                                                  "Genus")
group1.genera.plot<-group1.genera.plot+
  ggtitle(paste0("Relative abundance of Group 1 members"))
group1.genera.plot # 4 plots

# Group 2 increased genera
group2.increased.genera
group2.increased.genera.plot<-show.mean.abundances.and.plot(group2.increased.genera,
                                                            ps.q.agg.genus.relab.nmr,
                                                            "Genus")
group2.increased.genera.plot<-group2.increased.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members increased with age"))
group2.increased.genera.plot # 4 plots

# Group 2 unhealthy genera
group2.unhealthy.genera
group2.unhealthy.genera.plot<-show.mean.abundances.and.plot(group2.unhealthy.genera,
                                                            ps.q.agg.genus.relab.nmr,
                                                            "Genus")
group2.unhealthy.genera.plot<-group2.unhealthy.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.genera.plot # 7 plots

# Group 2 unhealthy families
group2.unhealthy.families
group2.unhealthy.families.plot<-show.mean.abundances.and.plot(group2.unhealthy.families,
                                                              ps.q.agg.family.relab.nmr,
                                                              "Family")
group2.unhealthy.families.plot<-group2.unhealthy.families.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.families.plot # 2 plots

# Group 3 genera
group3.genera
group3.genera.plot<-show.mean.abundances.and.plot(group3.genera,
                                                  ps.q.agg.genus.relab.nmr,
                                                  "Genus")
group3.genera.plot<-group3.genera.plot+
  ggtitle(paste0("Relative abundance of Group 3 members"))
group3.genera.plot # 5 plots


# Group 3 genera
group3.families.plot<-show.mean.abundances.and.plot(group3.families,
                                                    ps.q.agg.family.relab.nmr,
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

# Analyse ASVs shared between age group ####
### Give ASVs shorter names: Genus, "ASV", OTU, OTU number.
# For example, Allobaculum_ASV_22
nmr.asv.names<-ps.q.agg.relab.nmr%>%
  select(OTU,Genus)%>%
  group_by(Genus)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(Genus,OTU)%>%
  mutate(row.index=row_number())%>%
  mutate(ASV_name=paste(Genus,row.index,sep="_ASV_"))

# The new name becomes the OTU column. The old name becomes OTU_old_name column.
ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")])%>%
  rename("OTU_old_name"="OTU",
         "OTU"="ASV_name")%>%
  relocate(OTU,.before = Sample)%>%
  relocate(OTU_old_name,.after = Genus)

### How many ASVs are shared between two age group ####
otu.young<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(OTU,agegroup,MeanRelativeAbundanceAgegroup)
otu.old<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(OTU,agegroup,MeanRelativeAbundanceAgegroup)

shared.otu<-intersect(otu.young$OTU,otu.old$OTU)
length(shared.otu)

### How much % do shared ASVs take on average? ####
ps.q.agg.relab.nmr%>%
  # no separation by age
  filter(OTU%in%shared.otu)%>%
  group_by(Sample)%>%
  summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
  summarise(MeanRelAbSharedASVTotal=mean(SumRelAbSharedASV))

### How much % do shared ASVs take in each age group? ####
ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  # separation by age
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbSharedASVTotalAge=mean(SumRelAbSharedASV))

### Extract ASVs and the corresponding genera ####
shared.otu.genera<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  # keep unique rows
  distinct(OTU,.keep_all = T)%>%
  select(Genus,OTU)%>%
  group_by(Genus)%>%
  # count rows (ASVs) for each genus; add as a column for sorting
  mutate(num_otus=n())%>%
  # genera with the highest number of ASVs will be on top
  arrange(-num_otus)%>%
  # remove the column
  select(-num_otus)%>%
  ungroup()

#### Are the shared ASVs enriched in certain genera? ####
shared.otu.genera%>%
  group_by(Genus)%>%
  # count rows (ASVs) for each genus
  tally%>%
  arrange(-n)%>%
  # cumulative sum shows how many ASVs the top genera take
  mutate(cum_sum=cumsum(n))

# Are the asvs of the five most common genera all present in the dataframe of shared ASVs? ####
muribaculaceae.asv<-ps.q.agg.relab.nmr%>%
  filter(Genus=="Muribaculaceae")%>%
  distinct(OTU)%>%
  pull
# Subset the dataframe of shared ASVs to keep only Muribaculaceae ASVs, then
# extract ASVs as vector
table(muribaculaceae.asv%in%subset(shared.otu.genera,Genus=="Muribaculaceae")$OTU)

treponema.asv<-ps.q.agg.relab.nmr%>%
  filter(Genus=="Treponema")%>%
  distinct(OTU)%>%
  pull
table(treponema.asv%in%subset(shared.otu.genera,Genus=="Treponema")$OTU)

lachnospiraceae_family.asv<-ps.q.agg.relab.nmr%>%
  filter(Genus=="Lachnospiraceae Family")%>%
  distinct(OTU)%>%
  pull
table(lachnospiraceae_family.asv%in%subset(shared.otu.genera,Genus=="Lachnospiraceae Family")$OTU)

oscillospiraceae_family.asv<-ps.q.agg.relab.nmr%>%
  filter(Genus=="Oscillospiraceae Family")%>%
  distinct(OTU)%>%
  pull
table(oscillospiraceae_family.asv%in%subset(shared.otu.genera,Genus=="Oscillospiraceae Family")$OTU)

bacteria_kingdom.asv<-ps.q.agg.relab.nmr%>%
  filter(Genus=="Bacteria Kingdom")%>%
  distinct(OTU)%>%
  pull
table(bacteria_kingdom.asv%in%subset(shared.otu.genera,Genus=="Bacteria Kingdom")$OTU)

### What is the average relative abundance of each genus in the top 5 of shared.otu.genera? ####
ps.q.agg.relab.nmr%>%
  # the subset command relies on the fact that we sorted the shared.otu.genera
  # by the number of ASVs. So, the unique(shared.otu.genera$Genus)[1:5] has 
  # genera with the highest number of ASVs
  filter(OTU%in%subset(shared.otu.genera,Genus%in%unique(shared.otu.genera$Genus)[1:5])$OTU)%>%
  group_by(Genus)%>%
  # sum Abundance for each genus
  mutate(TotalGenus=sum(Abundance))%>%
  # Average relative abundance of each genus (we're in the dataframe of ASVs, not
  # genera!)
  mutate(MeanRelAbGenus=TotalGenus/TotalClass*100)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelAbGenus)%>%
  select(Genus,MeanRelAbGenus)

# Bar plot of abundances: high variability
ps.q.agg.relab.nmr%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  
  filter(OTU%in%subset(shared.otu.genera,Genus%in%unique(shared.otu.genera$Genus)[1:5])$OTU)%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=Genus))+
  geom_bar(stat="identity")
  # facet_grid(~agegroup,
  #            space = "free", # bars will have same widths
  #            scales="free")

# Finished up to this line ^^^^ ################################################ 
# Which shared ASVs are the most abundant? ####
### Most abundant ASVs on average ####
top10.asv.average<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Genus,OTU,MeanRelativeAbundance)%>%
  head(n=10)

### Bar plot shows the 10 most abundant ASVs on average account for 30-40% of samples ####
ps.q.agg.relab.nmr%>%
  filter(OTU%in%top10.asv.average$OTU)%>%
  mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=new_OTU))+
  geom_bar(stat="identity")+
  facet_grid(~agegroup,
             scales="free",
             space = "free")

### Most abundant ASVs in each age group ####
ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU,agegroup)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Genus,OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)

top10.asv.young<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu,agegroup=="agegroup0_10")%>%
  group_by(OTU,agegroup)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)

top10.asv.old<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu,agegroup=="agegroup10_16")%>%
  group_by(OTU,agegroup)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)

# Are top 10 most abundant ASVs same in two age groups?
setequal(top10.asv.young$OTU[1:10],top10.asv.old$OTU[1:10])
intersect(top10.asv.young$OTU[1:10],top10.asv.old$OTU[1:10])

# Intersect of the top 10 ASVs in each of the two age groups
ps.q.agg.relab.nmr%>%
  filter(OTU%in%intersect(top10.asv.young$OTU[1:10],top10.asv.old$OTU[1:10]))%>%
  mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=new_OTU))+
  geom_bar(stat="identity")+
  facet_grid(~agegroup,
             scales="free",
             space = "free")
# Union of the top 10 ASVs in each of the two age groups 
top10.asv.union<-sort(union(top10.asv.young$OTU[1:10],top10.asv.old$OTU[1:10]))
set.seed(1)
otu.fill<-createPalette(length(top10.asv.union),
                             seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(otu.fill)<-top10.asv.union

ps.q.agg.relab.nmr%>%
  filter(OTU%in%top10.asv.union)%>%
  mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU))+
  geom_bar(stat="identity")+
  scale_fill_manual(labels=names(otu.fill),
                    values=otu.fill)+
  theme_bw()+
  coord_cartesian(expand = FALSE)+
  labs(x="Sample",
       y="Relative abundance (%)",
       fill="ASV",
       title="Top 10 most abundant ASVs across age")+
  theme(axis.title.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        strip.text.x = ggtext::element_markdown(size=20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")
ggsave(file.path("./images/barplots",
                 paste(paste(format(Sys.time(),format="%Y%m%d"),
                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                       "top10-asv.png")),
       plot=last_plot(),
       width = 4500,height = 3500,
       units = "px",dpi=300,device = "png")

# M40 sample is very different #### 
m40.asvs<-ps.q.agg.relab.nmr%>% 
  filter(Sample=="M40")%>%
  select(OTU,Genus,RelativeAbundance,MeanRelativeAbundance,
         MeanRelativeAbundanceAgegroup)%>%
  arrange(-RelativeAbundance)%>%
  head(n=10)

set.seed(1)
m40.otu.fill<-createPalette(length(m40.asvs$OTU),
                        seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(m40.otu.fill)<-sort(m40.asvs$OTU)

ps.q.agg.relab.nmr%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  filter(OTU%in%m40.asvs$OTU)%>%
  # mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU))+
  geom_bar(stat="identity")+
  scale_fill_manual(labels=names(m40.otu.fill),
                    values=m40.otu.fill)+
  theme_bw()+
  coord_cartesian(expand = FALSE)+
  labs(x="Sample",
       y="Relative abundance (%)",
       fill="ASV",
       title="Top 10 most abundant ASVs in M40 sample")+
  theme(axis.title.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        strip.text.x = ggtext::element_markdown(size=20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")
ggsave(file.path("./images/barplots",
                 paste(paste(format(Sys.time(),format="%Y%m%d"),
                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                       "top10-asv-m40.png")),
       plot=last_plot(),
       width = 4500,height = 3500,
       units = "px",dpi=300,device = "png")






# Import the rarefied dataframe ####
ps.q.df.preprocessed.date_time<-"20240524_13_58_11"
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables/",authorname,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      paste0("ps.q.df.","rare"),"nonfiltered","OTU",
      paste("NMR",collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

# Or between species
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables/",authorname,paste0(
    paste(
      "20240426_22_00_04",
      paste0("ps.q.df.","rare"),"nonfiltered","Genus",
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

# Let's find which major ASVs are specific to one age group
young.ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

old.ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)
  


all.div.filtered<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%intersect(young.ps.q.agg.relab.nmr$OTU,
                   old.ps.q.agg.relab.nmr$OTU))%>%
  group_by(Sample)%>%
  reframe(sobs=specnumber(Abundance), # richness (num of species)
          shannon=diversity(Abundance,index = "shannon"),
          # simpson=diversity(Abundance, index="simpson"),
          invsimpson=diversity(Abundance, index="invsimpson"), # inverse simpson
          tot=sum(Abundance),
          agegroup=agegroup,
          sex=sex)%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()

# There are still differences
ggplot(all.div.filtered[all.div.filtered$metric %in%
                 plot.metrics,],
       # aes(x=reorder(class,-value),y=value,fill=class))+
       aes(x=factor(all.div.filtered$agegroup,
                    level=custom.levels),y=value,fill=factor(agegroup)))+
  geom_boxplot(show.legend = FALSE)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             ncol=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()

# How many genera are shared between two age groups ####
genera.young<-ps.q.agg.genus.relab.nmr%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,agegroup,MeanRelativeAbundanceAgegroup)
genera.old<-ps.q.agg.genus.relab.nmr%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,agegroup,MeanRelativeAbundanceAgegroup)

shared.genera<-intersect(genera.young$Genus,genera.old$Genus)
length(shared.genera)

# How much % do common genera take in each age group on average?
ps.q.agg.relab.nmr%>%
  filter(Genus%in%shared.genera)%>%
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbCommonASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbCommonASVTotalAge=mean(SumRelAbCommonASV))

#### TODO ####
# Do rare genera become more abundant? ####
# Add rank to ASVs/genera
genera.young.ranked<-ps.q.agg.genus.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,agegroup,row.index,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)

genera.old.ranked<-ps.q.agg.genus.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,agegroup,row.index,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)

otu.young.ranked<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,OTU,agegroup,row.index,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)


otu.old.ranked<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  mutate(row.index=as.numeric(rownames(.)))%>%
  select(Genus,OTU,agegroup,row.index,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)


# Outliers
ps.q.agg.genus.relab.nmr%>%
  filter(Genus%in%c("Bacteroides","Roseburia","Faecalibacterium"))%>%
  # distinct(Genus,.keep_all = T)%>%
  # select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=agegroup))+
  geom_bar(stat="identity")+
  facet_wrap(~Genus,scales="free_y")

ps.q.agg.genus.nmr.meanrelab<-ps.q.agg.genus.relab.nmr%>%
  distinct(Genus,agegroup,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup,agegroup)%>%
  left_join(genera.young.ranked[,c("Genus","agegroup","row.index")],by=join_by(Genus,agegroup))%>%
  left_join(genera.old.ranked[,c("Genus","agegroup","row.index")],by=join_by(Genus,agegroup))%>%
  mutate(row.index.x=ifelse(is.na(row.index.x)&!is.na(row.index.y),row.index.y,row.index.x),
         row.index.y=ifelse(is.na(row.index.y)&!is.na(row.index.x),row.index.x,row.index.y))%>%
  select(-row.index.y)%>%
  rename(row.index=row.index.x)

ps.q.agg.nmr.meanrelab<-ps.q.agg.relab.nmr%>%
  distinct(OTU,agegroup,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(Genus,OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup,agegroup)%>%
  left_join(otu.young.ranked[,c("OTU","agegroup","row.index")],by=join_by(OTU,agegroup))%>%
  left_join(otu.old.ranked[,c("OTU","agegroup","row.index")],by=join_by(OTU,agegroup))%>%
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

ps.q.agg.nmr.meanrelab.full<-ps.q.agg.nmr.meanrelab%>%
  group_by(OTU)%>%
  mutate(num.obs=n())%>%
  arrange(num.obs)%>%
  filter(num.obs==1)%>%
  mutate(agegroup=ifelse(agegroup=="agegroup0_10","agegroup10_16","agegroup0_10"),
         MeanRelativeAbundanceAgegroup=0,
         row.index=NA)%>%
  full_join(ps.q.agg.nmr.meanrelab,by=
              join_by(OTU,
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
  mutate(young_vs_old=old-young)%>%
  mutate(young_vs_old_fold=(old-young)/young)%>%
  # filter(agegroup0_10!=0,agegroup10_16!=0)%>%
  arrange(-young_vs_old_fold)%>%
  
  mutate(ind.dif=row.ind_young-row.ind_old)%>%
  filter(row.ind_old<=100,row.ind_old!=0)%>%
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


# Filter from outliers and recalculate everything ####
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




## Hierarchical clustering ####
library(cluster)
library(factoextra)
relations<-read.table("./data/metadata/pooled-metadata/nmr-relations.tsv",
                      header = T,
                      sep = "\t")
#load data
df <- USArrests
ps.q.df.wide<-ps.q.df.preprocessed%>%
  filter(Sample!="M40")%>%
  pivot_wider(names_from = "OTU", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()%>%
  select(-class)%>%
  column_to_rownames("Sample")
# Transform the data
ps.q.df.wide<-decostand(ps.q.df.wide,method="rclr",logbase = exp)

#remove rows with missing values
df <- na.omit(df)
ps.q.df.wide<-na.omit(ps.q.df.wide)

#scale each variable to have a mean of 0 and sd of 1
df <- scale(df)
# No scaling on compositional data

#view first six rows of dataset
head(df)
head(ps.q.df.wide[,1:10])

#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}

get_agglom_coef<-function(tax.df,link.method){
  agnes(tax.df,method=link.method)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)
for(mm in m){
  print(paste(mm,get_agglom_coef(ps.q.df.wide,mm)))
}
  
agnes(df,method="average")
agnes(ps.q.df.wide,method="ward")
# We can see that Ward’s minimum variance method produces the highest 
# agglomerative coefficient, thus we’ll use that as the method for our final hierarchical clustering:
  
#perform hierarchical clustering using Ward's minimum variance
clust <- agnes(df, method = "ward")
clust <- agnes(ps.q.df.wide,method="ward")


#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram")

#calculate gap statistic for each number of clusters (up to 10 clusters)
gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
gap_stat <- clusGap(ps.q.df.wide, FUN = hcut, nstart = 25, K.max = 10, B = 50)

#produce plot of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
# gap statistic is highest at k = 4 clusters. Thus, we’ll choose to group our observations into 4 distinct clusters.
# Apply Cluster Labels to Original Dataset
#compute distance matrix
d <- dist(df, method = "euclidean")
d <- dist(ps.q.df.wide, method = "euclidean")

#perform hierarchical clustering using Ward's method
final_clust <- hclust(d, method = "ward.D2" )

#cut the dendrogram into 4 clusters
groups <- cutree(final_clust, k=4)

#find number of observations in each cluster
table(groups)

#append cluster labels to original data
final_data <- cbind(USArrests, cluster = groups)
final_data <- as.data.frame(groups)%>%
  rename("cluster"="groups")%>%
  rownames_to_column("Sample")%>%
  left_join(custom.md[,c("Sample","class")],by="Sample")

#display first six rows of final data
head(final_data)

#find mean values for each cluster
aggregate(final_data, by=list(cluster=final_data$cluster), mean)


####
#perform hierarchical clustering using Ward's minimum variance
clust.average<-agnes(ps.q.df.wide,method="average")
clust.single<-agnes(ps.q.df.wide,method="single")
clust.complete<-agnes(ps.q.df.wide,method="complete")
clust.ward<-agnes(ps.q.df.wide,method="ward")

par(mfrow=c(2,2))
#produce dendrogram
pltree(clust.average, cex = 0.6, hang = -1, main = "Dendrogram")
pltree(clust.single, cex = 0.6, hang = -1, main = "Dendrogram")
pltree(clust.complete, cex = 0.6, hang = -1, main = "Dendrogram")
pltree(clust.ward, cex = 0.6, hang = -1, main = "Dendrogram")
par(mfrow=c(1,1))

# Kmeans clustering
# Silhouette method to choose optimal # of clusters
fviz_nbclust(ps.q.df.wide,kmeans,method="silhouette") +
  labs(subtitle = "Silhouette method for kmeans")
# Elbow method to choose optimal # of clusters
fviz_nbclust(ps.q.df.wide, kmeans, method = "wss") +
  geom_vline(xintercept = 1, linetype = 2) + # add line for better visualisation
  labs(subtitle = "Elbow method for kmeans") # add subtitle
# Gap statistic
set.seed(1)
fviz_nbclust(ps.q.df.wide, kmeans,
             nstart = 25,
             method = "gap_stat",
             nboot = 500 # reduce it for lower computation time (but less precise results)
) +
  labs(subtitle = "Gap statistic method for kmeans")

# Silhouette for Hierarchical cluster
fviz_nbclust(ps.q.df.wide,hcut,method="silhouette") +
  labs(subtitle = "Silhouette method for Hclust")
# Elbow method to choose optimal # of clusters (hierarchical)
fviz_nbclust(ps.q.df.wide, hcut, method = "wss") +
  geom_vline(xintercept = 1, linetype = 2) + # add line for better visualisation
  labs(subtitle = "Elbow method for Hclust") # add subtitle
# Gap statistic for hcut
set.seed(1)
fviz_nbclust(ps.q.df.wide, hcut,
             nstart = 25,
             method = "gap_stat",
             nboot = 500 # reduce it for lower computation time (but less precise results)
) +
  labs(subtitle = "Gap statistic method for Hcut")
