#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---
#' 
#' ```{r, setup 002-summary-stats-qiime2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# Analysing QIIME2 data with phyloseq. ####
#' # Analysing QIIME2 data with phyloseq.
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' In this script, we will explore the imported dataset from QIIME2 using qiime2R
#' and phyloseq.
#' 
#' We will also rarefy the data for future analyses.
#' 
#' We will use the data from 001-phyloseq-qiime2.R script (ps.q.agg
#' agglomerated tables at phylum, family, genus, and OTU level).

#+ echo=FALSE
## 1. Load necessary libraries and scripts. ####
#'
#' ## Load necessary libraries and scripts.
# install.packages(c("tidyverse","vegan","ggtext","Polychrome","ggrepel"))
library(tidyverse)
library(vegan)
library(Polychrome)
library(ggtext)
library(ggrepel)

#+ echo=FALSE
## 2. Import the config file. ####
#' 
#' ## Import the config file.
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
print(paste("The analysis focus is:", active.analysis))

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
#+ echo=FALSE
## 3. Import datasets as rds files. ####
#'
#' Import datasets as rds files.
ps.q.agg.asv<-readRDS(file = ps.q.agg.asv.fname)
ps.q.agg.genus <- readRDS(file = ps.q.agg.genus.fname)
ps.q.agg.family <- readRDS(file = ps.q.agg.family.fname)
ps.q.agg.phylum<-readRDS(file = ps.q.agg.phylum.fname)
custom.md <- readRDS(custom.md.path)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
#' Load necessary scripts.
source(file.path(util.functions.r,"get_n_uniq_taxa_per_host.R"))
source(file.path(util.functions.r,"create_summary_stats_table.R"))
source(file.path(util.functions.r,"add_relab_to_tax_df.R"))
source(file.path(util.functions.r,"add_agegroup_to_tax_df.R"))
source(file.path(util.functions.r,"get_unclassified_summary_stats.R"))
source(file.path(util.functions.r,"get_dominant_taxa_in_host.R"))
source(file.path(util.functions.r,"ggplot_species.R"))
source(file.path(util.functions.r,"get_rarefied_table.R"))

#+ echo=FALSE
## 4. Calculating summary statistics. ####
#'
#' ## Calculating summary statistics.
#' 
#+ echo=FALSE
### 4.1 Check the total number of unique ASV/phyla/families/genera per class. ####
#' 
#' ### Check the total number of unique ASV/phyla/families/genera per class.
#' We will use it for the summary table in the next section.
n.asv.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.asv,"OTU")
n.phylum.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.phylum,"Phylum")
n.family.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.family,"Family")
n.genus.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.genus,"Genus")

#+ echo=FALSE
### 4.2 Create a summary table. ####
#'
#' ### Create a summary table.
#' The columns are:  
#' - Total reads  
#' - Library size (mean Abundance ± SD)  
#' - Number of ASV per host  
#' - Number of phyla per host  
#' - Number of families per host  
#' - Number of genera per host  
summary.stats.table<-create_summary_stats_table(ps.q.agg.asv,
                                                n.asv.per.host,
                                                n.phylum.per.host,
                                                n.family.per.host,
                                                n.genus.per.host)
knitr::kable(as.data.frame(summary.stats.table),format = "simple")
summary.stats.table.fname<-file.path(community.composition.tables,
  "summary-table.tsv")
if(!file.exists(summary.stats.table.fname)){
  write.table(summary.stats.table,
              file=summary.stats.table.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 5. Add relative abundance and average relative abundance columns. ####
#'
#' ## Add relative abundance and average relative abundance columns.
ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")
ps.q.agg.asv.relab<-add_relab_to_tax_df(ps.q.agg.asv,"OTU")
knitr::kable(head(ps.q.agg.phylum.relab),format = "simple")
knitr::kable(head(ps.q.agg.family.relab),format = "simple")
knitr::kable(head(ps.q.agg.genus.relab),format = "simple")
knitr::kable(head(ps.q.agg.asv.relab),format = "simple")

#+ echo=FALSE
## 6. Add agegroup variable to NMR data. ####
#'
#' ## Add agegroup variable to NMR data.
custom.md.ages<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")

ps.q.agg.genus.relab.nmr<-ps.q.agg.genus.relab%>%
  filter(class=="NMR")

ps.q.agg.asv.relab.nmr<-ps.q.agg.asv.relab%>%
  filter(class=="NMR")

#' Add the age groups.
ps.q.agg.genus.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.genus.relab.nmr,"Genus",
                                             custom.md.ages)
ps.q.agg.asv.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.asv.relab.nmr,"OTU",
                                       custom.md.ages)
knitr::kable(head(ps.q.agg.genus.relab.nmr),format = "simple")
knitr::kable(head(ps.q.agg.asv.relab.nmr),format = "simple")

#+ echo=FALSE
## 7. Calculate summary stats of unclassified genera in each animal (unrarefied). ####
#'
#' ## Calculate summary stats of unclassified genera in each animal (unrarefied). ####
#' Here, we are interested in unclassified genera, but you can also try Families,
#' Orders, Classes, etc.
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
agglom.rank<-"Genus"
unclassified.genus.summary.stats.table<-
  get_unclassified_summary_stats(ps.q.agg.genus.relab,"Genus")
knitr::kable(unclassified.genus.summary.stats.table, format ="simple")

unclassified.genus.summary.stats.table.fname <-
  file.path(community.composition.tables,
            "unclassified-genus-summary-table.tsv")
if (! file.exists(unclassified.genus.summary.stats.table.fname)){
  write.table(unclassified.genus.summary.stats.table,
              file = unclassified.genus.summary.stats.table.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 8. Rarefy the table and check the percentage of unclassified taxa. ####
#'
#' ## Rarefy the table and check the percentage of unclassified taxa. ####
ps.q.agg.asv.rare<-get_rarefied_table(ps.q.agg.asv,"OTU",custom.levels)
ps.q.agg.phylum.rare<-get_rarefied_table(ps.q.agg.phylum,"Phylum",custom.levels)
ps.q.agg.family.rare<-get_rarefied_table(ps.q.agg.family,"Family",custom.levels)
ps.q.agg.genus.rare<-get_rarefied_table(ps.q.agg.genus,"Genus",custom.levels)
ps.q.agg.genus.nmr.rare<-get_rarefied_table(ps.q.agg.genus.relab.nmr,"Genus","NMR")
ps.q.agg.asv.nmr.rare<-get_rarefied_table(ps.q.agg.asv.relab.nmr,"OTU","NMR")

knitr::kable(head(ps.q.agg.asv.relab.nmr),format = "simple")
knitr::kable(head(ps.q.agg.genus.nmr.rare), format = "simple")
knitr::kable(head(ps.q.agg.asv.nmr.rare), format = "simple")


#+ echo=FALSE
### 8.1 Add relative abundances and taxonomic information to the rarefied dataframe. ####
#' 
#' ### Add relative abundances and taxonomic information to the rarefied dataframe. 
#' All hosts (genus)
ps.q.agg.genus.rare.relab<-add_relab_to_tax_df(ps.q.agg.genus.rare,"Genus")
knitr::kable(head(ps.q.agg.genus.rare.relab),format = "simple")

#' Add other taxonomic ranks to the dataframe
ps.q.agg.genus.rare.relab<-ps.q.agg.genus.rare.relab%>%
  left_join(unique(ps.q.agg.genus[,c("Kingdom","Phylum","Class","Order","Family","Genus")]))
knitr::kable(head(ps.q.agg.genus.rare.relab),format = "simple")

#' Add relative abundances to NMR dataframe (ASV level)
ps.q.agg.asv.nmr.rare.relab<-add_relab_to_tax_df(ps.q.agg.asv.nmr.rare,"OTU")
ps.q.agg.asv.nmr.rare.relab<-ps.q.agg.asv.nmr.rare.relab%>%
  left_join(unique(ps.q.agg.asv[,c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")]))

### 8.2 Plot a rarefaction curve. ####
# ps.q.mat<-as(t(otu_table(ps.q)),"matrix") # from phyloseq
# ps.q.genus.mat<-ps.q.agg.genus%>%
#   filter(class %in% custom.levels,Abundance!=0)%>%
#   dplyr::select(Sample,Abundance,class,all_of(agglom.rank))%>%
#   filter(Abundance!=0)%>%
#   pivot_wider(names_from = all_of(agglom.rank),
#               values_from = "Abundance",
#               values_fill = 0)%>%
#   as.data.frame()%>%
#   column_to_rownames("Sample")%>% # Set sample names as row names
#   dplyr::select(-class)%>%
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
# ggsave(file.path(community.composition.figures,
#               paste0(paste("rarecurve",
#                     truncationlvl,agglom.rank,
#                     sep = "-"),".png")),
#        plot=last_plot(),
#        width = 4500,height = 3000,
#        units = "px",dpi=300,device = "png")


#+ echo=FALSE
### 8.3 Create a summary stats table for the rarefied dataframe. ####
#' 
#' ### Create a summary stats table for the rarefied dataframe.
n.asv.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.asv.rare,"OTU")
n.phylum.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.phylum.rare,"Phylum")
n.family.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.family.rare,"Family")
n.genus.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.genus.rare,"Genus")

summary.stats.table.rare<-create_summary_stats_table(ps.q.agg.asv.rare,
                                                n.asv.per.host.rare,
                                                n.phylum.per.host.rare,
                                                n.family.per.host.rare,
                                                n.genus.per.host.rare)
knitr::kable(summary.stats.table.rare, format ="simple")
summary.stats.table.rare.fname<-file.path(community.composition.tables,
                                     "summary-table-rarefied.tsv")
if(!file.exists(summary.stats.table.rare.fname)){
  write.table(summary.stats.table.rare,
              file=summary.stats.table.rare.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 9. Calculate summary stats of unclassified genera in rarefied data. ####
#' 
#' ## Calculate summary stats of unclassified genera in rarefied data.
unclassified.genus.summary.stats.table.rare<-get_unclassified_summary_stats(ps.q.agg.genus.rare.relab,
                                                                            "Genus")
knitr::kable(unclassified.genus.summary.stats.table.rare, format ="simple")
unclassified.genus.summary.stats.table.rare.fname<-
  file.path(community.composition.tables,
            paste("unclassified-genus-summary-table-rarefied.tsv",sep="-"))
if(!file.exists(unclassified.genus.summary.stats.table.rare.fname)){
  write.table(unclassified.genus.summary.stats.table.rare,
              file = unclassified.genus.summary.stats.table.rare.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 10. Check the most abundant phyla, families, genera in NMR and other hosts. ####
#'
#' ## Check the most abundant phyla, families, genera in NMR and other hosts.
#' Phyla:
ps.q.agg.dominant.phyla.all_hosts<-ps.q.agg.phylum.relab%>%
  group_by(class,Phylum)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class,Phylum, MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class,Phylum)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
knitr::kable(head(ps.q.agg.dominant.phyla.all_hosts), format = "simple")

ps.q.agg.dominant.phyla.nmr<-ps.q.agg.dominant.phyla.all_hosts%>%
  filter(class=="NMR")
knitr::kable(head(ps.q.agg.dominant.phyla.nmr), format ="simple")

#' Families:
ps.q.agg.dominant.families.all_hosts<-ps.q.agg.family.relab%>%
  group_by(class,Family)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class, Phylum, Family, MeanRelativeAbundance, sdRelativeAbundance, 
           min,max,n)%>%
  group_by(class,Phylum,Family)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
knitr::kable(head(ps.q.agg.dominant.families.all_hosts), format = "simple")

ps.q.agg.dominant.families.nmr<-ps.q.agg.dominant.families.all_hosts%>%
  filter(class=="NMR")
knitr::kable(head(ps.q.agg.dominant.families.nmr),format = "simple")

#' Genera:
ps.q.agg.dominant.genera.all_hosts<-ps.q.agg.genus.relab%>%
  group_by(class, Genus)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class,Phylum,Family,Genus, MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class,Phylum,Family,Genus)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
knitr::kable(head(ps.q.agg.dominant.genera.all_hosts), format = "simple")

ps.q.agg.dominant.genera.nmr<-ps.q.agg.dominant.genera.all_hosts%>%
  filter(class=="NMR")
knitr::kable(head(ps.q.agg.dominant.genera.nmr), format = "simple")

ps.q.agg.dominant.phyla.nmr.fname <- 
  file.path(community.composition.tables,# paste("20240523_12_02_33",
            "dominant-phyla-table-nmr.tsv")
if(!file.exists(ps.q.agg.dominant.phyla.nmr.fname)){
  write.table(ps.q.agg.dominant.phyla.nmr,
              file = ps.q.agg.dominant.phyla.nmr.fname,
              row.names = F,sep = "\t")
}

ps.q.agg.dominant.families.nmr.fname <-
  file.path(community.composition.tables,# paste("20240523_12_07_02",
            "dominant-families-table-nmr.tsv")
if(!file.exists(ps.q.agg.dominant.families.nmr.fname)){
  write.table(ps.q.agg.dominant.families.nmr,
              file = ps.q.agg.dominant.families.nmr.fname,
              row.names = F, sep = "\t")
}
ps.q.agg.dominant.genera.nmr.fname <-
  file.path(community.composition.tables,# paste("20240523_12_31_19",
            "dominant-genera-table-nmr.tsv")
if(!file.exists(ps.q.agg.dominant.genera.nmr.fname)){
  write.table(ps.q.agg.dominant.genera.nmr,
              file=ps.q.agg.dominant.genera.nmr.fname,
              row.names = F,sep = "\t")
}

# For all hosts
ps.q.agg.dominant.phyla.all_hosts.fname <- 
  file.path(community.composition.tables,
            "dominant-phyla-table-all_hosts.tsv")
if(!file.exists(ps.q.agg.dominant.phyla.all_hosts.fname)){
  write.table(ps.q.agg.dominant.phyla.all_hosts,
              file = ps.q.agg.dominant.phyla.all_hosts.fname,
              row.names = F,sep = "\t")
}

ps.q.agg.dominant.families.all_hosts.fname <-
  file.path(community.composition.tables,
            "dominant-families-table-all_hosts.tsv")
if(!file.exists (ps.q.agg.dominant.families.all_hosts.fname)){
  write.table(ps.q.agg.dominant.families.all_hosts,
              file = ps.q.agg.dominant.families.all_hosts.fname,
              row.names = F,sep = "\t")
}

ps.q.agg.dominant.genera.all_hosts.fname <-
  file.path(community.composition.tables,
            "dominant-genera-table-all_hosts.tsv")
if(! file.exists (ps.q.agg.dominant.genera.all_hosts.fname)){
  write.table(ps.q.agg.dominant.genera.all_hosts,
            file = ps.q.agg.dominant.genera.all_hosts.fname,
            row.names = F,sep = "\t")
}
#+ echo=FALSE
## 11. Check how much Bacteroidaceae are in NMR.  ####
#'
#' ## Check how much Bacteroidaceae are in NMR.
bacteroidaceae.nmr<-ps.q.agg.dominant.families.nmr%>%
  filter(Family=="Bacteroidaceae")
knitr::kable(bacteroidaceae.nmr, format = "simple")
bacteroidaceae.nmr.fname<- file.path(community.composition.tables,
                                     "bacteroidaceae-table-nmr.tsv")
if(! file.exists (bacteroidaceae.nmr.fname)){
  write.table(bacteroidaceae.nmr,
              file = bacteroidaceae.nmr.fname,
              row.names = F,sep = "\t")
}
#+ echo=FALSE
### 11.1 Check the most dominant Bacteroidota families in NMR. ####
#'
#'### Check the most dominant Bacteroidota families in NMR.
bacteroidota.nmr<-ps.q.agg.dominant.families.nmr%>%
  filter(Phylum=="Bacteroidota")
knitr::kable(head(bacteroidota.nmr), format = "simple")
bacteroidota.nmr.fname <-file.path(community.composition.tables,
                                   "bacteroidota-table-nmr.tsv")
if(! file.exists (bacteroidota.nmr.fname)){
  write.table(bacteroidota.nmr,
              file = bacteroidota.nmr.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 12. Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. ####
#'
#' ## Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. 
spirochaetaceae.nmr<-ps.q.agg.dominant.families.nmr%>%
  filter(Family=="Spirochaetaceae",class=="NMR")
knitr::kable(spirochaetaceae.nmr, format = "simple")

spirochaetota.nmr<-ps.q.agg.dominant.phyla.nmr%>%
  filter(Phylum=="Spirochaetota")
knitr::kable(spirochaetota.nmr, format = "simple")
spirochaetaceae.nmr.fname <-file.path(community.composition.tables,
                                      "spirochaetaceae-table-nmr.tsv")
if(! file.exists (spirochaetaceae.nmr.fname)){
  write.table(spirochaetaceae.nmr,
              file = spirochaetaceae.nmr.fname,
              row.names = F,sep = "\t")
  
}

treponema.nmr<-ps.q.agg.dominant.genera.nmr%>%
  filter(Genus=="Treponema")
knitr::kable(treponema.nmr, format = "simple")
treponema.nmr.fname <-file.path(community.composition.tables,
                                      "treponema-table-nmr.tsv")
if(! file.exists (treponema.nmr.fname)){
  write.table(treponema.nmr,
              file = treponema.nmr.fname,
              row.names = F,sep = "\t")
  
}

#+ echo=FALSE
### 12.1 Check the number of ASVs in Treponema from NMR. ####
#' 
#' ### Check the number of ASVs in Treponema from NMR.
ps.q.agg.asv%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

#+ echo=FALSE
## 13. Check Mogibacteriaceae (renamed to Anaerovoracaceae) in all hosts. ####
#' 
#' ## Check Mogibacteriaceae (renamed to Anaerovoracaceae) in all hosts.
mogibacteriaceae_anaerovoracaceae.all<-ps.q.agg.dominant.families.all_hosts%>%
  filter(Family=="Anaerovoracaceae")
knitr::kable(head(mogibacteriaceae_anaerovoracaceae.all), format = "simple")
mogibacteriaceae_anaerovoracaceae.all.fname <-file.path(community.composition.tables,
                                "mogibacteriaceae_anaerovoracaceae-all-table.tsv")
if(! file.exists (mogibacteriaceae_anaerovoracaceae.all.fname)){
  write.table(mogibacteriaceae_anaerovoracaceae.all,
              file = mogibacteriaceae_anaerovoracaceae.all.fname,
              row.names = F,sep = "\t")
  
}

#+ echo=FALSE
## 14. Analyse sulfur-metabolising bacteria in NMR. ####
#'
#' ## Analyse sulfur-metabolising bacteria in NMR. 
desulfobacterota.nmr<-ps.q.agg.dominant.genera.nmr%>%
  filter(Phylum=="Desulfobacterota",class=="NMR")
knitr::kable(head(desulfobacterota.nmr), format = "simple")

#' Desulfobacterota in other hosts
desulfobacterota.all<-ps.q.agg.dominant.genera.all_hosts%>%
  filter(Phylum=="Desulfobacterota")
knitr::kable(head(desulfobacterota.all), format = "simple")

#+ echo=FALSE
### 14.1 Total Desulfobacterota in NMR. ####
#' 
#' ### Total Desulfobacterota in NMR.
ps.q.agg.dominant.phyla.nmr%>%
  filter(Phylum=="Desulfobacterota")%>%
  knitr::kable(format = "simple")
desulfobacterota.nmr.fname <-file.path(community.composition.tables,
                                                        "desulfobacterota-table-nmr.tsv")
if(! file.exists (desulfobacterota.nmr.fname)){
  write.table(desulfobacterota.nmr,
              file = desulfobacterota.nmr.fname,
              row.names = F,sep = "\t")
  
}
desulfobacterota.all.fname <- file.path(community.composition.tables,
                                        "desulfobacterota-table-all.tsv")
if(! file.exists (desulfobacterota.all.fname)){
  write.table(desulfobacterota.all,
              file = desulfobacterota.all.fname,
              row.names = F,sep = "\t")
  
}

#+ echo=FALSE
### 14.2 Plot Desulfobacterota. ####
#'
#' ### Plot Desulfobacterota.
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
  taxa.plot.theme
#+ fig.height = 10, fig.width = 8
print(desulfobacterota.plot + 
  ggtitle(paste0("Relative abundance of Desulfobacterota phylum members"))+
  theme(plot.title = element_text(size = 14))
)

for (image.format in c("png","tiff")){
  ggsave(paste0("Desulfobacterota-phylum-members-all-hosts.",
                      image.format),
         plot = desulfobacterota.plot,
         path = community.composition.figures,
         width=8, height=10,units="in",
         dpi=300,device = image.format)
}


#+ echo=FALSE
## 15. Plot Treponema and other Spirochaetota. ####
#'
#' ## Plot Treponema and other Spirochaetota.
spirochaetota.plot<-ps.q.agg.genus.relab%>%
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
  taxa.plot.theme
#+ fig.height = 10, fig.width = 8
print(spirochaetota.plot + 
        ggtitle(paste0("Relative abundance of Spirochaetota phylum members"))+
        theme(plot.title = element_text(size = 14))
        )

for (image.format in c("png","tiff")){
  ggsave(paste0("Spirochaetota-phylum-members-all-hosts.",
                      image.format),
         plot=spirochaetota.plot,
         path = community.composition.figures,
         width=8, height=10,units="in",
         dpi=300,device = image.format)
}

#+ echo=FALSE
## 16. Analysis of naked mole-rat data ASVs. ####
#' 
#' ## Analysis of naked mole-rat data ASVs.
#' Give ASVs shorter names: Genus, "ASV", OTU, OTU number. For example, Allobaculum_ASV_22.
nmr.asv.names<-ps.q.agg.asv.relab.nmr%>%
  dplyr::select(OTU,Genus)%>%
  group_by(Genus)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(Genus,OTU)%>%
  mutate(row.index=row_number())%>%
  mutate(ASV_name=paste(Genus,row.index,sep="_ASV_"))
knitr::kable(head(nmr.asv.names), format = "simple")

#' The new name becomes the OTU column. The old name becomes OTU_old_name column.
ps.q.agg.asv.relab.nmr<-ps.q.agg.asv.relab.nmr%>%
  left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")])%>%
  rename("OTU_old_name"="OTU",
         "OTU"="ASV_name")%>%
  relocate(OTU,.before = Sample)%>%
  relocate(OTU_old_name,.after = Genus)

#+ echo=FALSE
### 16.1 How many ASVs are shared between two age groups? ####
#' ### How many ASVs are shared between two age groups?
#' First, find ASVs in young samples
otu.young<-ps.q.agg.asv.relab.nmr%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  dplyr::select(OTU,agegroup,MeanRelativeAbundanceAgegroup, sdRelativeAbundance)
#' Next, find ASVs in old samples
otu.old<-ps.q.agg.asv.relab.nmr%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  dplyr::select(OTU,agegroup,MeanRelativeAbundanceAgegroup)

#' 1745 ASVs in young individuals:
nrow(otu.young) 
#' 771 ASVs in old individuals:
nrow(otu.old) 

#' 668 shared ASVs:
shared.otu<-intersect(otu.young$OTU,otu.old$OTU)
length(shared.otu)

#' ASVs unique to old samples: 771 - 668 = 103
length(otu.old$OTU)-length(shared.otu)

#+ echo=FALSE
### 16.2 Check 103 ASV only in old individuals. ####
#'
#' ### Check 103 ASV only in old individuals.
ps.q.agg.asv.relab.nmr%>%
  group_by(OTU,agegroup)%>%
  mutate(n_samples=n_distinct(Sample))%>% # Find the number of old samples the ASV was found in
  filter(OTU %in% otu.old$OTU,
         !OTU %in% shared.otu)%>%
  filter(n_samples >= 3)%>% 
  arrange(desc(MeanRelativeAbundanceAgegroup))%>%
  dplyr::select(OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)
#' No ASVs with at least 3 samples! The 103 ASVs are individual-specific.

#' Which genera do the 103 old-specific ASVs belong to?
ps.q.agg.asv.relab.nmr%>%
  filter(OTU %in% otu.old$OTU, # ASV in old but not shared vector
         !OTU %in% shared.otu)%>%
  dplyr::select(OTU,Family,Genus, MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
  group_by(Family, Genus)%>%
  summarise(n=n_distinct(OTU))%>%
  arrange(desc(n))%>%
  head%>%
  knitr::kable(format = "simple")

#+ echo=FALSE
### 16.3 How much % do shared ASVs take on average? ####
#' 
#' ### How much % do shared ASVs take on average?
ps.q.agg.asv.relab.nmr%>%
  # no separation by age
  filter(OTU%in%shared.otu)%>%
  group_by(Sample)%>%
  summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
  summarise(MeanRelAbSharedASVTotal = mean(SumRelAbSharedASV),
            sdRelaAbSharedASVTotal = sd(SumRelAbSharedASV))
#' On average 94.2% ± 10.9%
#' 
#+ echo=FALSE
### 16.4 How much % do shared ASVs take in each age group? ####
#'
#'### How much % do shared ASVs take in each age group?
ps.q.agg.asv.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  # separation by age
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbSharedASVTotalAge=mean(SumRelAbSharedASV),
            sdRelaAbSharedASVTotal = sd(SumRelAbSharedASV))%>%
  knitr::kable(format = "simple")
#' Higher variation in young individuals
#' 
#+ echo=FALSE
### 16.5 Are the shared ASVs enriched in certain genera? ####
#'
#' ### Are the shared ASVs enriched in certain genera? 
#' First, get the table of shared ASVs and their genera.
shared.otu.genera<-ps.q.agg.asv.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  # keep unique rows
  distinct(Genus,OTU)%>%
  group_by(Genus)
knitr::kable(head(shared.otu.genera), format = "simple")

#' Calculate the ASVs in each genus with cumsum() and pull the most numerous genera.
shared.otu.genera.cumsum<-shared.otu.genera%>%
  group_by(Genus)%>%
  # count rows (ASVs) for each genus; add as a column for sorting
  summarise(num_asvs =n_distinct(OTU))%>%
  # genera with the highest number of ASVs will be on top
  arrange(-num_asvs,Genus)%>%
  ungroup()%>%
  # cumulative sum shows how many ASVs the top genera take
  mutate(cum_sum=cumsum(num_asvs))

knitr::kable(head(shared.otu.genera.cumsum), format = "simple")
#' Six genera account for 1/3 shared ASVs (225 out of 668 ASVS)

#+ echo=FALSE
### 16.6 How many ASVs of the top 5 genera in shared.otu.genera are actually found there? ####
#'
#' ### How many ASVs of the top 5 genera in shared.otu.genera are actually found there?
#' The five most common genera according to the shared.otu.genera dataframe are 
#' Lachnospiraceae Family (55 ASVs), Muribaculaceae (54 ASVS), 
#' Treponema (42 ASVs), Bacteria Kingdom (30 ASVs), and
#' Oscillospiraceae Family (22 ASVs).  
ps.q.agg.asv.relab.nmr%>%
  filter(Genus%in%c("Lachnospiraceae Family","Treponema",
                    "Muribaculaceae","Bacteria Kingdom",
                    "Oscillospiraceae Family"))%>%
  distinct(Genus,OTU)%>%
  mutate(is_shared=ifelse(OTU %in%shared.otu.genera$OTU, TRUE, FALSE))%>%
  group_by(Genus,is_shared)%>%
  summarise(n_asvs=n())%>%
  ungroup()%>%
  arrange(Genus,desc(is_shared))%>%
  knitr::kable(format = "simple")
#' Lachnospiraceae Family: 55 ASVs are shared, 108 are not.
#' Muribaculaceae: 54 are shared, 71 are not.
#' Treponema: 42 ASVs are shared, 67 are not.
#' Bacteria Kingdom: 30 ASVs are shared, 125 are not.
#' Oscillospiraceae Family: 22 ASVs are shared, 25 are not.

#+ echo=FALSE
### 16.7 What is the average relative abundance of each genus in the top 5 of shared.otu.genera? ####
#'
#' ### What is the average relative abundance of each genus in the top 5 of shared.otu.genera?
ps.q.agg.asv.relab.nmr%>%
  # the subset command relies on the fact that we sorted the shared.otu.genera
  # by the number of ASVs. So, the unique(shared.otu.genera$Genus)[1:5] has 
  # genera with the highest number of ASVs
  filter(OTU%in%subset(shared.otu.genera, 
                       Genus %in% unique(shared.otu.genera.cumsum$Genus)[1:5])$OTU)%>%
  group_by(Genus,Sample)%>%
  summarise(TotalGenus=sum(RelativeAbundance))%>%
  group_by(Genus)%>%
  summarise(MeanRelAbGenus = mean(TotalGenus),
            sdRelAbGenus = sd(TotalGenus))%>%
  arrange(-MeanRelAbGenus)%>%
  knitr::kable(format = "simple")

#+ echo=FALSE
## 17. Plot ASVs in NMR. ####
#' 
#' ## Plot ASVs in NMR. 
#' Setup sample levels for NMR for barplots.
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  dplyr::select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))

#' Bar plot of abundances: high variability
#+ fig.height = 6, fig.width = 11
ps.q.agg.asv.relab.nmr%>%
  left_join(custom.md.ages)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  filter(OTU%in%subset(shared.otu.genera,
                       Genus%in%unique(shared.otu.genera.cumsum$Genus)[1:5])$OTU)%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=Genus))+
  geom_bar(stat="identity")+
  labs(x = "Sample")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#+ echo=FALSE
### 17.1 10 most abundant shared ASVs on average account for 30-40% of samples (barplot). ####
#'
#' ### 10 most abundant shared ASVs on average account for 30-40% of samples (barplot).
#' Find the most abundant ASVs on average.
top10.asv.average<-ps.q.agg.asv.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,MeanRelativeAbundance)%>%
  head(n=10)
top10.asv.average%>%
  arrange(Genus)%>%
  knitr::kable(format = "simple")

#' Bar plot shows the 10 most abundant shared ASVs on average account for 30-40% of samples
#+ fig.height = 6, fig.width = 11
ps.q.agg.asv.relab.nmr%>%
  filter(OTU%in%top10.asv.average$OTU)%>%
  mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=new_OTU))+
  geom_bar(stat="identity")+
  facet_grid(~agegroup,
             scales="free",
             space = "free")+
  labs(x = "Sample")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


### 17.2 Most abundant shared ASVs in each age group. ####
#'
#' ### Most abundant shared ASVs in each age group.
top.shared.asvs.by_age<-ps.q.agg.asv.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU,agegroup)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,agegroup,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
  ungroup()
knitr::kable(head(top.shared.asvs.by_age), format = "simple")

top10.shared.asv.young<-top.shared.asvs.by_age%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  head(n=10)

top10.shared.asv.old<-top.shared.asvs.by_age%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  head(n=10)

#' Are top 10 most abundant ASVs same in two age groups?
setequal(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)
#' No. How many are common?
intersect(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)
#' Seven ASVs
#'
#' Union of the top 10 ASVs in each of the two age groups.
top10.shared.asv.union<-sort(union(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU))
#' Make the names shorter for the barplot
top10.shared.asv.union<-top10.shared.asv.union%>%
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
#' Prepare a custom fill with Polychrome package
set.seed(1)
otu.fill<-createPalette(nrow(top10.shared.asv.union),
                        seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(otu.fill)<-top10.shared.asv.union$new_OTU

#+ echo=FALSE
### 17.3 Barplot of the most abundant ASVs. ####
#' 
#' ### Barplot of the most abundant ASVs.
top10.shared.asv.plot<-ps.q.agg.asv.relab.nmr%>%
  left_join(custom.md.ages)%>%
  filter(OTU%in%top10.shared.asv.union$OTU) %>%
  # mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  left_join(top10.shared.asv.union)%>%
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
  project_theme +
  asv.barplot.theme +
  theme (plot.title = element_text(size = 8))

#+ fig.height = 6, fig.width = 8
print(top10.shared.asv.plot+
  ggtitle("Top 10 most abundant ASVs across age")+
  theme(plot.title = element_text(size = 14)))

for(image.format in image.formats){
  ggsave(paste0("top10-asv.", image.format),
       plot=top10.shared.asv.plot,
       path = community.composition.figures,
       width=8, height=6,units="in",
       dpi=300,device = image.format)
}

#+ echo=FALSE
### 17.4 M40 sample is very different. #### 
#' 
#' ### M40 sample is very different. 
m40.asvs<-ps.q.agg.asv.relab.nmr%>% 
  filter(Sample=="M40")%>%
  dplyr::select(OTU,Genus,RelativeAbundance,MeanRelativeAbundance,
         MeanRelativeAbundanceAgegroup)%>%
  arrange(-RelativeAbundance)%>%
  head(n=10)

set.seed(1)
m40.otu.fill<-createPalette(length(m40.asvs$OTU),
                            seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                            range=c(30, 80))
names(m40.otu.fill)<-sort(m40.asvs$OTU)
m40.asv.plot<-ps.q.agg.asv.relab.nmr%>%
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
  guides(fill = guide_legend(ncol = 3))+
  asv.barplot.theme

#+ fig.height = 6, fig.width = 10
print(m40.asv.plot)
ggsave(paste("top10-asv-m40.png"),
       plot = m40.asv.plot,
       path = community.composition.figures,
       width=8, height=6,units="in",
       dpi=300,device = "png")

#+ echo=FALSE
## 18. Analyse age groups in detail. ####
#'
#' ## Analyse age groups in detail. 

#+ echo=FALSE
### 18.1 Let's find which major ASVs are specific to one age group ####
#' 
#' ### Let's find which major ASVs are specific to one age group ####
young.ps.q.agg.asv.relab.nmr<-ps.q.agg.asv.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

old.ps.q.agg.asv.relab.nmr<-ps.q.agg.asv.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

length(unique(young.ps.q.agg.asv.relab.nmr$OTU))
length(unique(old.ps.q.agg.asv.relab.nmr$OTU))

#+ echo=FALSE
### 18.2 How many genera are shared between two age groups ####
#' 
#' ### How many genera are shared between two age groups ####
genera.young<-ps.q.agg.genus.relab.nmr%>%
  # add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  # left_join(custom.md.ages)%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)
genera.old<-ps.q.agg.genus.relab.nmr%>%
  # add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  # left_join(custom.md.ages)%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)

length(genera.young$Genus)
length(genera.old$Genus)
shared.genera<-intersect(genera.young$Genus,genera.old$Genus)
length(shared.genera)

#+ echo=FALSE
### 18.3 Show 5 genera found in old but not young NMR ####
#' 
#' ### Show 5 genera found in old but not young NMR ####
unique.to_old.genera<-setdiff(genera.old$Genus,genera.young$Genus)
ps.q.agg.genus.relab.nmr%>%
  # add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  filter(Genus %in% unique.to_old.genera)%>%
  dplyr::select(Sample,Family,Genus,MeanRelativeAbundance,
                MeanRelativeAbundanceAgegroup )%>%
  knitr::kable(format = "simple")

# How much % do common genera take in each age group on average?
ps.q.agg.asv.relab.nmr%>%
  filter(Genus%in%shared.genera)%>%
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbCommonASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbCommonASVTotalAge=mean(SumRelAbCommonASV))%>%
  knitr::kable(format = "simple")

## 19. Check non-bacterial data ####
ps.q.agg.genus%>%
  filter(Kingdom != "Bacteria")%>%
  distinct(Genus, .keep_all = T)%>%
  select(class, Kingdom:Genus, MeanRelativeAbundance)%>%
  knitr::kable(format = "simple")

ps.q.agg.asv%>%
  filter(Kingdom != "Bacteria")%>%
  group_by(class, Genus)%>%
  summarise(n = n())

ps.q.agg.asv%>%
  filter(Order == "Methanomassiliicoccales")
  

sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), "active.analysis"))
gc()