#' ```{r, setup phyloseq-qiime2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("code/r-scripts/002-summary-stats-qiime2.R", knit = FALSE)
#' #file.rename("code/r-scripts/002-summary-stats-qiime2.Rmd", "markdown/002-summary-stats-qiime2.Rmd")
#' #rmarkdown::render('./markdown/002-summary-stats-qiime2.Rmd', 'html_document',
#' # knit_root_dir="/home/rakhimov/projects/amplicon_nmr/")
#' ```


#' In this script, we will explore the imported dataset from QIIME2 (using qiime2R
#' and phyloseq).
#' 
#' We will also rarefy the data for future analyses.
#' 
#' We will use the data from 001-phyloseq-qiime2.R script (ps.q.agg
#' agglomerated tables at phylum, family, genus, and OTU level).

## 1. Import libraries. ####
# install.packages(c("tidyverse","vegan","ggtext","Polychrome","ggrepel"))
library(tidyverse)
library(vegan)
library(Polychrome)
library(ggtext)
library(ggrepel)
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
source("./code/r-scripts/create_summary_stats_table.R")
source("./code/r-scripts/add_relab_to_tax_df.R")
source("./code/r-scripts/add_agegroup_to_tax_df.R")
source("./code/r-scripts/get_unclassified_summary_stats.R")
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/ggplot_species.R")
source("./code/r-scripts/add_zero_rows.R")

## 2. Specifying parameters and directory/file names #### 
authorname<-"pooled" # name of the folder with QIIME2 output

rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)

truncationlvl<-"234" # truncation level that we chose in QIIME2

read.end.type<-"single" # single reads or paired reads: decided in QIIME2
# Import datasets as rds files
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,# old
  paste("20260211_17_01_07","phyloseq-qiime",authorname,"OTU",read.end.type, # new 
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
  paste("20260211_17_01_10","phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
  paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
  paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))


barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")

## 3. Setup plots ####
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "DMR" = "*Fukomys damarensis*",
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

# Setup general ggplot theme
mytheme<-theme(axis.text.y = element_text(size=10), # size of y axis ticks
               axis.title = element_text(size = 10), # size of axis names
               legend.text = element_text(size = 10), # size of legend text
               legend.title = element_text(size = 15), # size of legend title
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank()
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
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

# 5. Add relative abundance and average relative abundance columns ####
ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")
ps.q.agg.relab<-add_relab_to_tax_df(ps.q.agg,"OTU")

# 6. Add agegroup (Must run for plotting) ####
custom.md.ages<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")

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
ps.q.agg.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.relab.nmr,"OTU",
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

# 8. Rarefy the table and check the percentage of unclassified taxa (Table 2) ####
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
  return(tax.df.rare)
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
ps.q.agg.dominant.phyla.all_hosts<-ps.q.agg.phylum.relab%>%
  group_by(class)%>%
  distinct(class,Phylum, MeanRelativeAbundance,sdRelativeAbundance)%>%
  group_by(class,Phylum)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()

ps.q.agg.dominant.phyla.nmr<-ps.q.agg.dominant.phyla.all_hosts%>%
  filter(class=="NMR")

ps.q.agg.dominant.families.all_hosts<-ps.q.agg.family.relab%>%
  group_by(class)%>%
  distinct(class,Phylum,Family, MeanRelativeAbundance,sdRelativeAbundance)%>%
  group_by(class,Phylum,Family)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()

ps.q.agg.dominant.families.nmr<-ps.q.agg.dominant.families.all_hosts%>%
  filter(class=="NMR")

ps.q.agg.dominant.genera.all_hosts<-ps.q.agg.genus.relab%>%
  group_by(class)%>%
  distinct(class,Phylum,Family,Genus, MeanRelativeAbundance,sdRelativeAbundance)%>%
  group_by(class,Phylum,Family,Genus)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()

ps.q.agg.dominant.genera.nmr<-ps.q.agg.dominant.genera.all_hosts%>%
  filter(class=="NMR")

# TODO: This function is not necessary? Maybe metagenome uses it
# ps.q.agg.dominant.phyla.nmr<-get_dominant_taxa_in_host(ps.q.agg.phylum,
#                                                        "Phylum","NMR")
# head(ps.q.agg.dominant.phyla.nmr)
# ps.q.agg.dominant.families.nmr<-get_dominant_taxa_in_host(ps.q.agg.family,
#                                                           "Family","NMR")
# head(ps.q.agg.dominant.families.nmr,n=20)
# ps.q.agg.dominant.genera.nmr<-get_dominant_taxa_in_host(ps.q.agg.genus,
#                                                         "Genus","NMR")
# head(ps.q.agg.dominant.genera.nmr)


write.table(ps.q.agg.dominant.phyla.nmr,
            file=file.path(rtables.directory,
                           # paste("20240523_12_02_33", 
                            paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "dominant-phyla-table-nmr.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.families.nmr,
            file=file.path(rtables.directory,
                           # paste("20240523_12_07_02",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "dominant-families-table-nmr.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.genera.nmr,
            file=file.path(rtables.directory,
                           # paste("20240523_12_31_19",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "dominant-genera-table-nmr.tsv",sep="-")),
            row.names = F,sep = "\t")
# For all hosts
write.table(ps.q.agg.dominant.phyla.all_hosts,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "dominant-phyla-table-all_hosts.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.families.all_hosts,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "dominant-families-table-all_hosts.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.genera.all_hosts,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "dominant-genera-table-all_hosts.tsv",sep="-")),
            row.names = F,sep = "\t")

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
desulfobacterota.all.mean%>%left_join(desulfobacterota.all.sd)

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
desulfobacterota.plot<-ps.q.agg.genus.relab%>%
  mutate(class=factor(class,levels=custom.levels))%>%
  filter(Phylum=="Desulfobacterota")%>%
  group_by(class,Genus)%>%
  ggplot(aes(x=factor(class,levels=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE,
               outliers = FALSE)+
  geom_jitter(aes(color=class),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=2,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
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
  scale_fill_viridis_d(option="C")+
  stat_summary(fun=median, geom="point", shape=23, size=2, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  # ggtitle(paste0("Relative abundance of Desulfobacterota phylum members"))+
  mytheme+
  theme(axis.text.y = ggtext::element_markdown(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(size=10),
        strip.text.x = element_text(size=10),
        plot.title = element_text(size = 8), # size of plot title
        plot.caption = element_text(size=8), # size of plot caption
        legend.position = "none")

for (image.format in c("png","tiff")){
  ggsave(paste0("./images/taxaboxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "Desulfobacterota phylum members-all-hosts",
                      sep = "-"),".",image.format),
         plot=desulfobacterota.plot,
         width=8, height=10,units="in",
         # width = 4000,height = 6000,
         # width = 1200,
         # units = "px",
         dpi=300,device = image.format)
}


# Plot Treponema ####
desulfobacterota.plot<-ps.q.agg.genus.relab%>%
  mutate(class=factor(class,levels=custom.levels))%>%
  filter(Phylum=="Spirochaetota")%>%
  group_by(class,Genus)%>%
  ggplot(aes(x=factor(class,levels=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE,
               outliers = FALSE)+
  geom_jitter(aes(color=class),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=2,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
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
  scale_fill_viridis_d(option="C")+
  stat_summary(fun=median, geom="point", shape=23, size=2, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  # ggtitle(paste0("Relative abundance of Spirochaetota phylum members"))+
  mytheme+
  theme(axis.text.y = ggtext::element_markdown(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(size=10),
        strip.text.x = element_text(size=10),
        plot.title = element_text(size = 8), # size of plot title
        plot.caption = element_text(size=8), # size of plot caption
        legend.position = "none")
for (image.format in c("png","tiff")){
  ggsave(paste0("./images/taxaboxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "Spirochaetota phylum members-all-hosts",
                      sep = "-"),".",image.format),
         plot=desulfobacterota.plot,
         width=8, height=10,units="in",
         # width = 4000,height = 6000,
         # width = 1200,
         # units = "px",
         dpi=300,device = image.format)
}



# Analyse ASVs shared between age groups ####
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
  # left_join(custom.md.ages[,c("Sample","agegroup")])%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(OTU,agegroup,MeanRelativeAbundanceAgegroup)
otu.old<-ps.q.agg.relab.nmr%>%
  # left_join(custom.md.ages[,c("Sample","agegroup")])%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  select(OTU,agegroup,MeanRelativeAbundanceAgegroup)

shared.otu<-intersect(otu.young$OTU,otu.old$OTU)
length(shared.otu) # 668 shared

length(otu.young$OTU) # 1745  in young individuals
length(otu.old$OTU) # 771  in old individuals

length(otu.old$OTU)-length(shared.otu)
# Check 103 ASV only in old individuals #####
ps.q.agg.relab.nmr%>%
  group_by(OTU,agegroup)%>%
  mutate(n_samples=n_distinct(Sample))%>% # Find the number of old samples the ASV was found in
  filter(OTU %in% otu.old$OTU,
         !OTU %in% shared.otu)%>%
  filter(n_samples >= 3)%>% # No ASVs with n_samples >3! 
  arrange(desc(MeanRelativeAbundanceAgegroup))%>%
  select(OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)


ps.q.agg.relab.nmr%>%
  filter(OTU %in% otu.old$OTU, # ASV in old but not shared vector
         !OTU %in% shared.otu)%>%
  select(OTU,Family,Genus, MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
  # filter(MeanRelativeAbundanceAgegroup>0.01)%>%
  mutate(Family_Genus=paste(Family, Genus,sep="_"))%>%
  group_by(Family_Genus)%>%
  summarise(n=n_distinct(OTU))%>%
  arrange(desc(n))

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


## 19. Setup sample levels for NMR####
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))


# Bar plot of abundances: high variability
ps.q.agg.relab.nmr%>%
  left_join(custom.md.ages)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  
  filter(OTU%in%subset(shared.otu.genera,Genus%in%unique(shared.otu.genera$Genus)[1:5])$OTU)%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=Genus))+
  geom_bar(stat="identity")
# facet_grid(~agegroup,
#            space = "free", # bars will have same widths
#            scales="free")

# Which shared ASVs are the most abundant? ####
### Most abundant ASVs on average ####
top10.asv.average<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Genus,OTU,MeanRelativeAbundance)%>%
  head(n=10)
top10.asv.average%>%
  arrange(Genus)

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
top10.asv.union<-top10.asv.union%>%
  as_tibble()%>%
  rename("OTU"="value")%>%
  mutate(new_OTU = OTU,
         new_OTU = gsub("Allobaculum_", "Allob. ",new_OTU),
         new_OTU = gsub("Erysipelotrichaceae Family_", "Erysip. F. ",new_OTU),
         new_OTU = gsub("Eubacteriaceae Family_", "Eubac. F. ",new_OTU),
         new_OTU = gsub("Fibrobacter_", "Fibrob. ",new_OTU),
         new_OTU = gsub("Muribaculaceae_", "Murib. ",new_OTU),
         new_OTU = gsub("Paludibacteraceae Family_", "Palud. F. ",new_OTU),
         new_OTU = gsub("o5_ASV", "o5 ASV",new_OTU),
         new_OTU = gsub("Prevotella_", "Prev. ",new_OTU),
         new_OTU = gsub("Prevotellaceae Family_", "Prevot. F. ",new_OTU),
         new_OTU = gsub("Prevotellaceae_UCG", "Prevot. UCG",new_OTU),
         new_OTU = gsub("UCG-001_", "UCG-001 ",new_OTU),
         new_OTU = gsub("UCG-003_", "UCG-003 ",new_OTU)
         )

set.seed(1)
otu.fill<-createPalette(nrow(top10.asv.union),
                        seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(otu.fill)<-top10.asv.union$new_OTU

top10.asv.plot<-ps.q.agg.relab.nmr%>%
  left_join(custom.md.ages)%>%
  filter(OTU%in%top10.asv.union$OTU) %>%
  # mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  left_join(top10.asv.union)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=new_OTU))+
  geom_bar(stat="identity")+
  scale_fill_manual(labels=names(otu.fill),
                    values=otu.fill)+
  theme_bw()+
  coord_cartesian(expand = FALSE)+
  labs(x="Sample",
       y="Relative abundance (%)",
       # title="Top 10 most abundant ASVs across age",
       fill="ASV"
       )+
  mytheme +
  theme(
    legend.title = element_text(size=10),
    legend.text = element_text(size=9),
    legend.key.size = unit(0.3, 'cm'), #change legend key size
    legend.key.spacing.y = unit(0.1, "lines"), # distant between key text
    legend.box.spacing = unit(0.1,"lines"), # space between the plot and the legend box
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle=45,size=10,hjust=1),# rotate 
    strip.text.x = ggtext::element_markdown(size=10),
    panel.spacing = unit(0.8, "cm"), # increase distance between facets
    plot.title = element_text(size = 8), # size of plot title
    plot.caption = element_text(size=8), # size of plot caption
    legend.position = "bottom")
for(image.format in image.formats){
  ggsave(file.path("./images/barplots",
                 paste0(paste(format(Sys.time(),format="%Y%m%d"),
                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                       "-top10-asv", ".",image.format)),
       plot=top10.asv.plot,
       width=8, height=6,units="in",
       # width = 5000,height = 3500,
       # units = "px",
       dpi=300,device = image.format)
}
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
  left_join(custom.md.ages)%>%
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
# Between NMR (ASV level)
# ps.q.df.preprocessed.date_time<-"20240524_13_58_11"
ps.q.df.preprocessed.date_time<-"20260211_17_14_21"
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables/",authorname,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      paste0("ps.q.df.","rare"),"nonfiltered","OTU",
      paste("NMR",collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

# Or between species (Genus level)
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables/",authorname,paste0(
    paste(
      # "20240426_22_00_04",
      "20260211_17_14_18",
      paste0("ps.q.df.","rare"),"nonfiltered","Genus",
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

# Let's find which major ASVs are specific to one age group ####
young.ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

old.ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

length(unique(young.ps.q.agg.relab.nmr$OTU))
length(unique(old.ps.q.agg.relab.nmr$OTU))


# How many genera are shared between two age groups ####
genera.young<-ps.q.agg.genus.relab.nmr%>%
  add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  left_join(custom.md.ages)%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)
genera.old<-ps.q.agg.genus.relab.nmr%>%
  add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  left_join(custom.md.ages)%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)

length(genera.young$Genus)
length(genera.old$Genus)
shared.genera<-intersect(genera.young$Genus,genera.old$Genus)
length(shared.genera)


# Show 5 genera found in old but not young NMR ####
unique.to_old.genera<-setdiff(genera.old$Genus,genera.young$Genus)
ps.q.agg.genus.relab.nmr%>%
  add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  filter(Genus %in% unique.to_old.genera)%>%
  select(Sample,Family,Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup )


# How much % do common genera take in each age group on average?
ps.q.agg.relab.nmr%>%
  filter(Genus%in%shared.genera)%>%
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbCommonASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbCommonASVTotalAge=mean(SumRelAbCommonASV))