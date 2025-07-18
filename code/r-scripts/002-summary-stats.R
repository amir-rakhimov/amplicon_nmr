# In this script, we will explore the imported dataset from QIIME2 (using qiime2R
# and phyloseq).
# We will also rarefy the data for future analyses.
# We will use the data from 001-phyloseq-qiime2.R script (ps.q.agg
# agglomerated tables at phylum, family, genus, and OTU level).

## 1. Import libraries ####
# install.packages(c("tidyverse","vegan","ggtext","Polychrome","ggrepel"))
library(tidyverse)
library(vegan)
library(Polychrome)
library(ggtext)
library(ggrepel)
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/create_summary_stats_table.R")
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
source("./code/r-scripts/add_relab_to_tax_df.R")
source("./code/r-scripts/add_agegroup_to_tax_df.R")
source("./code/r-scripts/get_unclassified_summary_stats.R")
source("./code/r-scripts/ggplot_species.R")
source("./code/r-scripts/add_zero_rows.R")

## 2. Specifying parameters and directory/file names #### 
authorname<-"pooled" # name of the folder with QIIME2 output
# authorname<-"biagi" # name of the folder with QIIME2 output

rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)

truncationlvl<-"234" # truncation level that we chose in QIIME2
# truncationlvl<-"0" # truncation level that we chose in QIIME2

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

## 3. Setup plots ####
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
## 4. Summary statistics table (Table 1) ####
### 4.1 Check the total number of unique ASV/phyla/families/genera per class ####
# We will use it for the summary table
n.asv.per.host<-get_n_uniq_taxa_per_host(ps.q.agg,"OTU")
n.phylum.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.phylum,"Phylum")
n.family.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.family,"Family")
n.genus.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.genus,"Genus")

### 4.2 Get the summary stats: ####
# * Total reads
# * Library size (mean Abundance ± SD)
# * Num of ASV per host
# * Num of phyla per host
# * Num of families per host
# * Num of genera per host
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


# TODO: Not use? ------------------
# # Samples used in whole metagenome sequencing
# nmr.metagenome.samples<-c("2D10","2D14","G14","G18","H15",
#                           "H21","H3","H4","O15","Y51b",
#                           "Y66b")
# # Get the number of ASVs in all NMR samples from 16S (1834)
# get_n_uniq_taxa_per_host(subset(ps.q.agg,class=="NMR"),"OTU")
# get_n_uniq_taxa_per_host(subset(ps.q.agg,class=="NMR"),"Phylum")
# get_n_uniq_taxa_per_host(subset(ps.q.agg,class=="NMR"),"Family")
# get_n_uniq_taxa_per_host(subset(ps.q.agg,class=="NMR"),"Genus")
# # Get the number of ASVs in NMR samples from 16S that were also used in WMS (1175)
# get_n_uniq_taxa_per_host(subset(ps.q.agg,Sample%in%nmr.metagenome.samples),"OTU")
# get_n_uniq_taxa_per_host(subset(ps.q.agg,Sample%in%nmr.metagenome.samples),"Phylum")
# get_n_uniq_taxa_per_host(subset(ps.q.agg,Sample%in%nmr.metagenome.samples),"Family")
# get_n_uniq_taxa_per_host(subset(ps.q.agg,Sample%in%nmr.metagenome.samples),"Genus")
# 
# # Import kraken2 table for comparison
# kraken2.table<-readRDS(file="../metagenome/output/rdafiles/20241003_13_52_43-phyloseq-kraken2-Species-table.rds")
# kraken2.table<-kraken2.table%>%
#   filter(Kingdom=="Bacteria")
# kraken2.table%>%
#   distinct(Phylum)
# match(sort(pull(unique(subset(kraken2.table,select=Phylum)))),
#       sort(pull(unique(subset(ps.q.agg,Sample%in%nmr.metagenome.samples,select=Phylum)))))
# intersect(sort(pull(unique(subset(kraken2.table,select=Phylum)))),
#           sort(pull(unique(subset(ps.q.agg,Sample%in%nmr.metagenome.samples,select=Phylum)))))
# # Only bacterial taxa here
# get_n_uniq_taxa_per_host(kraken2.table,"Species")
# get_n_uniq_taxa_per_host(kraken2.table,"Phylum")
# get_n_uniq_taxa_per_host(kraken2.table,"Family")
# get_n_uniq_taxa_per_host(kraken2.table,"Genus")
# --------------

# 5. Add relative abundance and average relative abundance columns ####
ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")
ps.q.agg.relab<-add_relab_to_tax_df(ps.q.agg,"OTU")

# 6. Add agegroup (Must run for plotting) ####
custom.md.ages<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))
# we create these new levels for plots
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

# 7. Calculate summary stats of unclassified genera in each animal (unrarefied) (Table 1 and 2) ####
# Here, we are interested in unclassified genera, but you can also try Families,
# Orders, Classes, etc.
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
agglom.rank<-"Genus"
unclassified.genus.summary.stats.table<-
  get_unclassified_summary_stats(ps.q.agg.genus.relab,"Genus")

# write.table(unclassified.genus.summary.stats.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 8. Rarefy the table and check the percentage of unclassified taxa (Table 3) ####
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

### 8.1 Add relative abundances and taxonomic information to the rarefied dataframe ####
# All hosts (genus)
ps.q.agg.genus.rare.relab<-add_relab_to_tax_df(ps.q.agg.genus.rare,"Genus")
# Add other taxonomic ranks to the dataframe
ps.q.agg.genus.rare.relab<-ps.q.agg.genus.rare.relab%>%
  left_join(unique(ps.q.agg.genus[,c("Kingdom","Phylum","Class","Order","Family","Genus")]))

# Add relative abundances to NMR dataframe (ASV level)
ps.q.agg.nmr.rare.relab<-add_relab_to_tax_df(ps.q.agg.nmr.rare,"OTU")
ps.q.agg.nmr.rare.relab<-ps.q.agg.nmr.rare.relab%>%
  left_join(unique(ps.q.agg[,c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")]))

### 8.2 Plot a rarefaction curve ####
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

## 9. Calculate summary stats of unclassified genera in rarefied data (Table 2) ####
unclassified.genus.summary.stats.table.rare<-get_unclassified_summary_stats(ps.q.agg.genus.rare.relab,
                                                                            "Genus")
# write.table(unclassified.genus.summary.stats.table.rare,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table-rarefied.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 10. Check the most abundant phyla, families, genera in NMR (Results and Discussion) ####
ps.q.agg.dominant.phyla.nmr<-get_dominant_taxa_in_host(ps.q.agg.phylum,
                                                       "Phylum","NMR")
head(ps.q.agg.dominant.phyla.nmr)
ps.q.agg.dominant.families.nmr<-get_dominant_taxa_in_host(ps.q.agg.family,
                                                          "Family","NMR")
head(ps.q.agg.dominant.families.nmr,n=20)
ps.q.agg.dominant.genera.nmr<-get_dominant_taxa_in_host(ps.q.agg.genus,
                                                        "Genus","NMR")
head(ps.q.agg.dominant.genera.nmr)
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

# 11. Check how much Bacteroidaceae are in NMR  ####
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

### 11.1 Check the most dominant Bacteroidota families in NMR ####
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

# 12. Check Spirochaetaceae, Spirochaetota, and Treponema in NMR ####
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

### 12.1 Check the number of ASVs in Treponema from NMR ####
ps.q.agg%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

# 13. Check Mogibacteriaceae (renamed to Anaerovoracaceae) in NMR ####
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

# 14. Analyse sulfur-metabolising bacteria in NMR ####
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

### 14.1 Total Desulfobacterota in NMR ####
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

### 14.2 Plot Desulfobacterota ####
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


# Not in the paper ####
# Bar plot of the total number of unique genera per host ####
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
