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
library(phyloseq)
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
#' Load necessary scripts.
source(file.path(util.functions.r,"get_mean_and_sd_relab_from_phyloseq.R"))
source(file.path(util.functions.r,"create_summary_stats.R"))
source(file.path(util.functions.r,"multi_rarefy_phyloseq.R"))
source(file.path(util.functions.r,"get_unclassified_summary_stats.R"))
source(file.path(util.functions.r,"get_average_summary_stats.R"))
source(file.path(util.functions.r,"get_average_unclassified_summary_stats_rare.R"))

# source(file.path(util.functions.r,"add_relab_to_tax_df.R"))
# source(file.path(util.functions.r,"add_agegroup_to_tax_df.R"))
# source(file.path(util.functions.r,"get_dominant_taxa_in_host.R"))
# source(file.path(util.functions.r,"ggplot_species.R"))

#+ echo=FALSE
## 3. Import raw datasets as rds files. ####
#'
#' Import raw datasets as rds files.
ps.q <- readRDS(file = ps.q.raw.fname)
ps.q.rel <- readRDS(file = ps.q.rel.raw.fname)
custom.md <- readRDS(custom.md.path)

#+ echo=FALSE
## 4. Analysis of naked mole-rat data ASVs. ####
#' 
#' ## Analysis of naked mole-rat data ASVs.
#' Load the NMR age metadata.
custom.md.ages<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
#' Filter metadata to keep only samples found in the ps.q object.
custom.md.ages <- custom.md.ages%>%
  filter(Sample %in% ps.q@sam_data$Sample)%>%
  rownames_to_column(var = "original_sample")%>%
  column_to_rownames("Sample")
#' Filter the ps.q object to keep only NMR samples.
ps.q.nmr <- subset_samples(ps.q, class == "NMR")
#' Change the sample data to the dataframe with NMR ages.
sample_data(ps.q.nmr) <- custom.md.ages
#' Create relative abundance data by transforming the ps.q.nmr object.
ps.q.nmr.rel <-transform_sample_counts(ps.q.nmr, function(x) 100*x/sum(x)) 
#' Calculate mean, SD, min, and max relative abundances of genera and ASVs by
#' age group.
mean_sd_relab.nmr_ages.asv <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.nmr.rel, "agegroup", "OTU")
#' Calculate mean, SD, min, and max relative abundances of NMR ASVs in 
#' total (no grouping).
mean_sd_relab.nmr_no_groups.asv <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.nmr.rel, "class", "OTU")

knitr::kable(head(mean_sd_relab.nmr_ages.asv),format = "simple")
knitr::kable(head(mean_sd_relab.nmr_no_groups.asv),format = "simple")

#' Setup sample levels for NMR for barplots.
sample.levels<- custom.md.ages%>%
  ungroup()%>%
  as_tibble(rownames = "Sample")%>%
  dplyr::select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))

analyse_nmr_asvs <-function(ps, # phyloseq object with raw abundances
                            ps.rel, # phyloseq object with relative abundances
                            mean_sd_rel.ages.asv, # mean, sd, min, max statistics for age groups
                            mean_sd_rel.total.asv, # mean, sd, min, max statistics without grouping
                            sample.levels # sample names with ages
                            ){
  #' Give ASVs shorter names: Genus, "ASV", OTU, OTU number. For example, Allobaculum_ASV_22.
  nmr.asv.names <- ps.rel %>%
    psmelt()%>%
    rename("RelativeAbundance" = "Abundance")%>%
    filter(RelativeAbundance !=0)%>%
    dplyr::select(OTU,Genus)%>%
    group_by(Genus)%>%
    distinct(OTU,.keep_all = T)%>%
    arrange(Genus,OTU)%>%
    mutate(row.index=row_number())%>%
    mutate(ASV_name=paste(Genus,row.index,sep="_ASV_"))%>%  
    ungroup()
  print( paste("Total number of ASVs:", nrow(nmr.asv.names)))
  
  #' The `ASV_name` column is renamed to `OTU_tidy`.
  ps.q.nmr.rel.asv<- ps.rel %>%
    psmelt()%>%
    rename("RelativeAbundance" = "Abundance")%>%
    as_tibble()%>%
    filter(RelativeAbundance !=0)%>%
    left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")],
              join_by(OTU, Genus))%>%
    rename("OTU"="OTU",
           "OTU_tidy"="ASV_name")%>%
    relocate(OTU_tidy,.before = Sample)%>%
    relocate(OTU,.after = Genus)
  
  
  #+ echo=FALSE
  ### 4.1 How many ASVs are shared between two age groups? ####
  #' ### How many ASVs are shared between two age groups?
  #' First, find ASVs in young samples
  otu.young <- mean_sd_rel.ages.asv%>%
    filter(agegroup=="agegroup0_10",MeanRelativeAbundance !=0)%>%
    distinct(OTU,.keep_all = T)%>%
    arrange(-MeanRelativeAbundance)%>%
    dplyr::select(OTU,agegroup,MeanRelativeAbundance, sdRelativeAbundance,
                  minRelativeAbundance, maxRelativeAbundance,
                  n_samples)
  
  #' Next, find ASVs in old samples
  otu.old <- mean_sd_rel.ages.asv%>%
    filter(agegroup=="agegroup10_16",MeanRelativeAbundance !=0)%>%
    distinct(OTU,.keep_all = T)%>%
    arrange(-MeanRelativeAbundance)%>%
    dplyr::select(OTU,agegroup,MeanRelativeAbundance, sdRelativeAbundance,
                  minRelativeAbundance, maxRelativeAbundance,
                  n_samples)
  #' `r nrow(otu.young)` ASVs in young individuals:
  print(paste( nrow(otu.young), "ASVs are found in young individuals" ))
  #' `r nrow(otu.old)` ASVs in old individuals:
  print(paste(nrow(otu.old), "ASVs are found in old individuals" ))
  
  #' Shared ASVs:
  shared.otu <- intersect(otu.young$OTU,otu.old$OTU)
  print(paste(length(shared.otu), "ASVs are shared between young and old individuals" ))
  
  #' ASVs unique to young samples: `r nrow(otu.young)` - `r length(shared.otu)` = 
  #' `r nrow(otu.young)-length(shared.otu)`
  print(paste (nrow(otu.young) - length(shared.otu), "ASVs are unique to young individuals") )
  
  #' ASVs unique to old samples: `r nrow(otu.old)` - `r length(shared.otu)` = 
  #' `r nrow(otu.old)-length(shared.otu)`
  print(paste (nrow(otu.old) - length(shared.otu), "ASVs are unique to old individuals") )
  
  #+ echo=FALSE
  ### 4.2 Find ASVs that are unique to young and old individuals. ####
  #'
  #' ### Find ASVs that are unique to young and old individuals.
  
  #' Which genera do the `r length(otu.young$OTU)-length(shared.otu)` young-specific ASVs belong to? 
  #' Showing only the first 10 rows because there's too many.
  print("Showing genera that possess young-specific ASVS (first 10 rows only)")
  mean_sd_rel.ages.asv%>%
    filter(OTU %in% otu.young$OTU, # ASV in old but not shared vector
           !OTU %in% shared.otu)%>%
    dplyr::select(OTU,Family,Genus)%>%
    group_by(Family, Genus)%>%
    summarise(n=n_distinct(OTU))%>%
    arrange(desc(n))%>%
    head(n = 10)%>%
    knitr::kable(format = "simple")%>%
    print()
  
  #' Which genera do the `r length(otu.old$OTU)-length(shared.otu)` old-specific ASVs belong to?
  print("Showing genera that possess old-specific ASVS (first 10 rows only)")
  mean_sd_rel.ages.asv%>%
    filter(OTU %in% otu.old$OTU, # ASV in old but not shared vector
           !OTU %in% shared.otu)%>%
    dplyr::select(OTU,Family,Genus)%>%
    group_by(Family, Genus)%>%
    summarise(n=n_distinct(OTU))%>%
    arrange(desc(n))%>%
    head(n = 10)%>%
    knitr::kable(format = "simple")%>%
    print()
  
  #' Are there any ASVs that are unique to young individuals AND present in at 
  #' least 3 samples? Showing only the first 10 rows because there's too many.
  print("Showing ASVs unique to young individuals AND present in at least 3 samples (first 10 rows only)")
  mean_sd_rel.ages.asv%>%
    filter(OTU %in% otu.young$OTU,
           !OTU %in% shared.otu)%>%
    filter(n_samples >= 3)%>% 
    arrange(desc(MeanRelativeAbundance))%>%
    dplyr::select(Phylum, Family, Genus, OTU, MeanRelativeAbundance, 
                  sdRelativeAbundance,n_samples)%>%
    head(n = 10)%>%
    knitr::kable(format = "simple")%>%
    print()
  
  print("Showing ASVs unique to old individuals AND present in at least 3 samples (first 10 rows only)")
  #' Are there any ASVs that are unique to old individuals AND present in at 
  #' least 3 samples? 
  mean_sd_rel.ages.asv%>%
    filter(OTU %in% otu.old$OTU,
           !OTU %in% shared.otu)%>%
    filter(n_samples >= 3)%>% 
    arrange(desc(MeanRelativeAbundance))%>%
    dplyr::select(Phylum, Family, Genus, OTU, MeanRelativeAbundance, 
                  sdRelativeAbundance,n_samples)%>%
    head(n = 10)%>%
    knitr::kable(format = "simple")%>%
    print()
  
  #+ echo=FALSE
  ### 4.3 How much % do `r length(shared.otu)` shared ASVs take on average? ####
  #' 
  #' ### How much % do `r length(shared.otu)` shared ASVs take on average?
  print("Showing how much % do shared ASVs take on average")
  ps.q.nmr.rel.asv%>%
    filter(OTU%in%shared.otu)%>%
    group_by(Sample)%>%
    summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
    summarise(MeanRelAbSharedASVTotal = mean(SumRelAbSharedASV),
              sdRelaAbSharedASVTotal = sd(SumRelAbSharedASV))%>%
    print()
  #' Most of the NMR samples are occupied by the shared ASVs.
  #' 
  #+ echo=FALSE
  ### 4.4 How much % do `r length(shared.otu)` shared ASVs take in each age group? ####
  #'
  #'### How much % do `r length(shared.otu)` shared ASVs take in each age group?
  print("Showing how much % do shared ASVs take in each age group")
  ps.q.nmr.rel.asv%>%
    filter(OTU %in% shared.otu)%>%
    # separation by age
    group_by(Sample,agegroup)%>%
    summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
    arrange(agegroup)%>%
    group_by(agegroup)%>%
    summarise(MeanRelAbSharedASVTotalAge=mean(SumRelAbSharedASV),
              sdRelaAbSharedASVTotal = sd(SumRelAbSharedASV))%>%
    knitr::kable(format = "simple")%>%
    print()
  #' Higher variation in young individuals
  #+ echo=FALSE
  ### 4.5 Are the shared ASVs enriched in certain genera? ####
  #'
  #' ### Are the shared ASVs enriched in certain genera? 
  
  #' Calculate the ASVs in each genus with cumsum() and pull the most numerous genera.
  # shared.otu.genera.cumsum<-shared.otu.genera%>%
  shared.otu.genera.cumsum<-mean_sd_rel.ages.asv %>%
    filter(OTU %in% shared.otu)%>%
    # keep unique rows
    distinct(Genus,OTU)%>%
    group_by(Genus)%>%
    # count rows (ASVs) for each genus; add as a column for sorting
    summarise(num_asvs =n_distinct(OTU))%>%
    # genera with the highest number of ASVs will be on top
    arrange(-num_asvs,Genus)%>%
    ungroup()%>%
    # cumulative sum shows how many ASVs the top genera take
    mutate(cum_sum=cumsum(num_asvs))
  
  knitr::kable(head(shared.otu.genera.cumsum), format = "simple")%>%
    print()
  print(paste0("Six genera account for ",
              round(shared.otu.genera.cumsum$cum_sum[6]/length(shared.otu) *100),
              "% shared ASVs (", shared.otu.genera.cumsum$cum_sum[6],
              " out of ",length(shared.otu), ")"))
  
  #+ echo=FALSE
  ### 4.6 How many ASVs of the top 5 genera in `shared.otu.genera.cumsum` are shared? ####
  #'
  #' ### How many ASVs of the top 5 genera in `shared.otu.genera.cumsum` are shared?
  top.shared.and.not_shared <- ps.q.nmr.rel.asv%>%
    filter(Genus%in% pull(shared.otu.genera.cumsum[1:5,1]) )%>%
    distinct(Genus,OTU, OTU_tidy)%>%
    mutate(is_shared=ifelse(OTU %in% shared.otu, "Shared", "Not_shared"))%>%
    group_by(Genus,is_shared)%>%
    summarise(n_asvs=n())%>%
    ungroup()%>%
    arrange(Genus,desc(is_shared))%>%
    pivot_wider(names_from = is_shared,
                values_from = n_asvs)
  print("Showing the number of ASVs shared and not shared between young and old samlples")
  top.shared.and.not_shared%>%
    knitr::kable(format = "simple")%>%
    print()
  
  #+ echo=FALSE
  ### 4.7 What is the relative abundance of shared ASVs at the genus level? ####
  #'
  #' ### What is the relative abundance of shared ASVs at the genus level?
  shared.otu.in.top.5.genera.by.asv.cumsum <- mean_sd_rel.ages.asv%>%
    filter(OTU %in% shared.otu)%>%
    # keep unique rows
    distinct(Genus,OTU)%>%
    filter(Genus %in% shared.otu.genera.cumsum$Genus[1:5])%>%
    pull(OTU)
  print(paste("Top 5 genera comprise",length(shared.otu.in.top.5.genera.by.asv.cumsum),"shared ASVs."))
  
  
  #' Bar plot of abundances: Even if the number of shared ASVs might be high, 
  #' the relative abundance in top 5 genera is up to 25%
  #+ fig.height = 6, fig.width = 11
  print(ps.q.nmr.rel.asv%>%
    left_join(sample.levels, by = "Sample")%>%
    mutate(Sample = factor (Sample, levels = sample.levels$Sample))%>%
    filter(OTU %in% shared.otu.in.top.5.genera.by.asv.cumsum )%>%
    ggplot(aes(x = NewSample,y = RelativeAbundance,fill = Genus))+
    geom_bar(stat = "identity")+
    labs(x = "Sample",
         title = "Relative abundance of top 5 genera by the number of shared ASVs")+
    theme(axis.text.x = element_text(angle = 45, hjust=1)))
  
  
  ### 4.8 Find 10 most abundant shared ASVs in each age group. ####
  #'
  #' ### Find 10 most abundant shared ASVs in each age group.
  top10.shared.asv.young<- mean_sd_rel.ages.asv%>%  
    filter(agegroup=="agegroup0_10")%>%
    left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")],
              by = c("Genus", "OTU"))%>%
    rename("OTU"="OTU",
           "OTU_tidy"="ASV_name")%>%
    arrange(-MeanRelativeAbundance)%>%
    dplyr::select(-Kingdom, -Phylum, -Family)%>%
    relocate(OTU_tidy, .before = OTU)%>%
    head(n=10)
  
  top10.shared.asv.old<-mean_sd_rel.ages.asv%>%  
    filter(agegroup=="agegroup10_16")%>%
    left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")],
              by = c("Genus", "OTU"))%>%
    rename("OTU"="OTU",
           "OTU_tidy"="ASV_name")%>%
    arrange(-MeanRelativeAbundance)%>%
    dplyr::select(-Kingdom, -Phylum, -Family)%>%
    relocate(OTU_tidy, .before = OTU)%>%
    head(n=10)
  
  #' How many top ASVS are common between two age groups ?
  intersect(top10.shared.asv.young$OTU_tidy,top10.shared.asv.old$OTU_tidy)
  print(paste(length(intersect(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)),
              "top ASVs are shared between young and old individuals"))

  #' Union of the top 10 ASVs in each of the two age groups.
  top10.shared.asv.union <- sort(union(top10.shared.asv.young$OTU_tidy,
                                       top10.shared.asv.old$OTU_tidy))
  print(paste(length(top10.shared.asv.union), "top ASVs were identified among two groups"))
  #' Make the names shorter for the barplot
  top10.shared.asv.union<-top10.shared.asv.union%>%
    as_tibble()%>%
    rename("OTU_tidy"="value")%>%
    mutate(OTU_tidy_short = OTU_tidy,
           OTU_tidy_short = gsub("Allobaculum_", "Allob. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Erysipelotrichaceae Family_", "Erysip. F. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Eubacteriaceae Family_", "Eubac. F. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Fibrobacter_", "Fibrob. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Muribaculaceae_", "Murib. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Paludibacteraceae Family_", "Palud. F. ",OTU_tidy_short),
           OTU_tidy_short = gsub("o5_ASV", "o5 ASV",OTU_tidy_short),
           OTU_tidy_short = gsub("Phascolarctobacterium_", "Phascol.",OTU_tidy_short),
           OTU_tidy_short = gsub("Prevotella_", "Prev. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Prevotellaceae Family_", "Prevot. F. ",OTU_tidy_short),
           OTU_tidy_short = gsub("Prevotellaceae_UCG", "Prevot. UCG",OTU_tidy_short),
           OTU_tidy_short = gsub("UCG-001_", "UCG-001 ",OTU_tidy_short),
           OTU_tidy_short = gsub("UCG-003_", "UCG-003 ",OTU_tidy_short)
    )
  
  
  return(list(otu.young = otu.young,
              otu.old = otu.old,
              shared.otu = shared.otu,
              shared.otu.genera.cumsum = shared.otu.genera.cumsum,
              nmr.asv.names = nmr.asv.names,
              ps.q.nmr.rel.asv = ps.q.nmr.rel.asv,
              top10.shared.asv.union = top10.shared.asv.union))
}

#' Run the script on raw data
nmr.raw <- analyse_nmr_asvs(ps = ps.q.nmr,
                            ps.rel = ps.q.nmr.rel,
                 mean_sd_rel.ages.asv = mean_sd_relab.nmr_ages.asv,
                 mean_sd_rel.total.asv = mean_sd_relab.nmr_no_groups.asv,
                 sample.levels = sample.levels)

#+ echo=FALSE
## 5. Plot ASVs in NMR to see anomalous samples. ####
#' 
#' ## Plot ASVs in NMR to see anomalous samples. 

#+ echo=FALSE
### 5.1 How much % do 10 most abundant shared ASVs on average account for (barplot)? ####
#'
#' ### How much % do 10 most abundant shared ASVs on average account for (barplot)?
#' Find the most abundant ASVs on average.
top10.asv.average <- mean_sd_relab.nmr_no_groups.asv%>%
  filter(OTU %in% nmr.raw$shared.otu)%>%
  distinct(OTU,.keep_all = T)%>%
  ungroup() %>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,MeanRelativeAbundance)%>%
  head(n=10)%>%
  left_join(nmr.raw$ps.q.nmr.rel.asv[,c("OTU", "OTU_tidy")]%>%distinct() ,
            by = "OTU")%>%
  relocate(MeanRelativeAbundance , .after = last_col())

top10.asv.average%>%
  arrange(Genus)%>%
  knitr::kable(format = "simple")

#' This union contains mostly the same ASVs as the list of 10 most abundant 
#' ASVs on average
print(paste(length(intersect(nmr.raw$top10.shared.asv.union$OTU_tidy, top10.asv.average$OTU_tidy) ),
            "out of", length(nmr.raw$top10.shared.asv.union$OTU_tidy), "top ASVs are identical among two groups"))

#' Bar plot shows the 10 most abundant shared ASVs on average account for 
#' 30-40% of samples
#+ fig.height = 6, fig.width = 11
nmr.raw$ps.q.nmr.rel.asv%>%
  filter(OTU%in%top10.asv.average$OTU)%>%
  mutate(OTU_tidy_Genus=paste0(OTU_tidy," (",Genus,")"))%>%
  left_join(sample.levels, by = "Sample")%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU_tidy_Genus))+
  geom_bar(stat="identity")+
  labs(x = "Sample",
       fill = "ASV")+
  coord_cartesian(expand = c(T,F,F,F))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#' Prepare a custom fill with Polychrome package
set.seed(1)
otu.fill<-createPalette(N = 100,
                        seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(otu.fill)<-nmr.raw$top10.shared.asv.union$OTU_tidy_short

#+ echo=FALSE
### 5.3 Barplot of the most abundant ASVs. ####
#' 
#' ### Barplot of the most abundant ASVs.
top10.shared.asv.plot<-nmr.raw$ps.q.nmr.rel.asv%>%
  # left_join(custom.md.ages)%>%
  filter(OTU_tidy%in%nmr.raw$top10.shared.asv.union$OTU_tidy) %>%
  # mutate(OTU_tidy_short=paste0(OTU," (",Genus,")"))%>%
  left_join(nmr.raw$top10.shared.asv.union, by = "OTU_tidy")%>%
  left_join(sample.levels, by = "Sample")%>%
  mutate(Sample = factor (Sample, levels = sample.levels$Sample))%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU_tidy_short))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=otu.fill)+
  theme_bw()+
  coord_cartesian(expand = c(T,F,F,F))+
  labs(x="Sample",
       y="Relative abundance (%)",
       fill="ASV"
  )+
  project_theme +
  asv.barplot.theme +
  theme (plot.title = element_text(size = 8))

#+ fig.height = 6, fig.width = 8
print(top10.shared.asv.plot+
        ggtitle("Top 10 most abundant ASVs across age in raw data")+
        theme(plot.title = element_text(size = 14)))

for(image.format in image.formats){
  ggsave(paste0("top10-asv-raw.", image.format),
         plot=top10.shared.asv.plot,
         path = community.composition.figures,
         width=8, height=6,units="in",
         dpi=300,device = image.format)
}

#+ echo=FALSE
### 5.4 M40 sample is very different. #### 
#' 
#' ### M40 sample is very different. 
m40.asvs <- nmr.raw$ps.q.nmr.rel.asv%>%
  filter(Sample=="M40")%>%
  dplyr::select(OTU,OTU_tidy,Genus,RelativeAbundance)%>%
  arrange(-RelativeAbundance)%>%
  head(n=10)%>%
  mutate(OTU_tidy_short = OTU_tidy,
         OTU_tidy_short = gsub("Akkermansia_", "Akker. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Bacteroides_", "Bacteroid. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Blautia_", "Blautia ",OTU_tidy_short),
         OTU_tidy_short = gsub("Christensenellaceae_R-7_group_", "Christensen. R-7 group ",OTU_tidy_short),
         OTU_tidy_short = gsub("Eubacteriaceae Family_", "Eubac. F. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Ileibacterium_", "Ileibac. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Fibrobacter_", "Fibrob. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Lachnospiraceae Family_", "Lachnosp. F. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Muribaculaceae_", "Murib. ",OTU_tidy_short),
         OTU_tidy_short = gsub("Paludicola_", "Palud. ",OTU_tidy_short),
  )%>%
  select(-RelativeAbundance)

m40.asvs.uniq <- setdiff(m40.asvs$OTU_tidy_short, nmr.raw$top10.shared.asv.union$OTU_tidy_short)
names(otu.fill) <- c(nmr.raw$top10.shared.asv.union$OTU_tidy_short, m40.asvs.uniq)

m40.asv.plot <- nmr.raw$ps.q.nmr.rel.asv%>%
  left_join(sample.levels, by = "Sample")%>%
  mutate(Sample = factor (Sample, levels = sample.levels$Sample))%>%
  left_join(m40.asvs, by =c("OTU_tidy", "OTU","Genus"))%>%
  filter(OTU%in%m40.asvs$OTU)%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU_tidy_short))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = otu.fill)+
  theme_bw()+
  coord_cartesian(expand = c(T,F,F,F))+
  labs(x="Sample",
       y="Relative abundance (%)",
       fill="ASV",
       title="Top 10 most abundant ASVs in M40 sample")+
  guides(fill = guide_legend(ncol = 3))+
  project_theme +
  asv.barplot.theme 

#+ fig.height = 6, fig.width = 10
print(m40.asv.plot)
ggsave(paste("top10-asv-m40.png"),
       plot = m40.asv.plot,
       path = community.composition.figures,
       width=8, height=6,units="in",
       dpi=300,device = "png")


#+ echo=FALSE
## 6. Remove M40 and Y29 samples. ####
#'
#' ## Remove M40 and Y29 samples.
#' M40 shows anomalous composition, while Y29 was treated with antibiotics
#' prior to sequencing. You can see in the barplot that 
#' *Erysipelotrichaceae Family ASV_5* has been expanded in Y29. 
#' Therefore, we can't use Y29 for downstream analyses
#' because its microbiota composition has been significantly changed.
excluded.samples <- c("M40","Y29")

sample.levels.filtered <- sample.levels%>%
  filter(!Sample %in% excluded.samples)
ps.q.filtered <- subset_samples(ps.q, ! sample_names(ps.q) %in% excluded.samples)
ps.q.rel.filtered <- subset_samples(ps.q.rel, ! sample_names(ps.q.rel) %in% excluded.samples)
ps.q.nmr.filtered <- subset_samples(ps.q.nmr, ! sample_names(ps.q.nmr) %in% excluded.samples)
ps.q.nmr.filtered.rel <- transform_sample_counts(ps.q.nmr.filtered, function(x) 100*x/sum(x)) 
mean_sd_relab.nmr_ages.asv.filtered <- get_mean_and_sd_relab_from_phyloseq(ps.q.nmr.filtered.rel,
                                                                           group.name = "agegroup",
                                                                           tax.rank = "OTU")
mean_sd_relab.nmr_no_groups.asv.filtered <- get_mean_and_sd_relab_from_phyloseq(ps.q.nmr.filtered.rel,
                                                                           group.name = "class",
                                                                           tax.rank = "OTU")
nmr.filtered <- analyse_nmr_asvs(ps = ps.q.nmr.filtered,
                                 ps.rel = ps.q.nmr.filtered.rel,
                                 mean_sd_rel.ages.asv = mean_sd_relab.nmr_ages.asv.filtered,
                                 mean_sd_rel.total.asv = mean_sd_relab.nmr_no_groups.asv.filtered,
                                 sample.levels = sample.levels.filtered)
#+ echo=FALSE
### 6.1 Are the top 10 most abundant shared ASVs different in the filtered data? If so, how much % samples do they occupy?  ####
#'
#' ### Are the top 10 most abundant shared ASVs different in the filtered data? If so, how much % samples do they occupy?
#+ echo=FALSE
#' Find the most abundant ASVs on average in the filtered data.
top10.asv.average.filtered <- mean_sd_relab.nmr_no_groups.asv.filtered%>%
  filter(OTU %in% nmr.filtered$shared.otu)%>%
  distinct(OTU,.keep_all = T)%>%
  ungroup() %>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,MeanRelativeAbundance)%>%
  head(n=10)%>%
  left_join(nmr.filtered$ps.q.nmr.rel.asv[,c("OTU", "OTU_tidy")]%>%distinct() ,
            by = "OTU")%>%
  relocate(MeanRelativeAbundance , .after = last_col())

#' There are no differences between top 10 most abundant shared ASVs.
setdiff(top10.asv.average$OTU_tidy,top10.asv.average.filtered$OTU_tidy)
setdiff(top10.asv.average.filtered$OTU_tidy,top10.asv.average$OTU_tidy)

top10.asv.average.filtered%>%
  arrange(Genus)%>%
  knitr::kable(format = "simple")

#' Bar plot shows the 10 most abundant shared ASVs on average account for 
#' 30-40% of samples
#+ fig.height = 6, fig.width = 11
nmr.filtered$ps.q.nmr.rel.asv%>%
  filter(OTU%in%top10.asv.average.filtered$OTU)%>%
  mutate(OTU_tidy_Genus=paste0(OTU_tidy," (",Genus,")"))%>%
  left_join(sample.levels, by = "Sample")%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU_tidy_Genus))+
  geom_bar(stat="identity")+
  labs(x = "Sample",
       fill = "ASV")+
  coord_cartesian(expand = c(T,F,F,F))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#+ echo=FALSE
### 6.2 Barplot of the most abundant ASVs in the filtered data. ####
#' 
#' ### Barplot of the most abundant ASVs in the filtered data. 
top10.shared.asv.plot.filtered<-nmr.filtered$ps.q.nmr.rel.asv%>%
  # left_join(custom.md.ages)%>%
  filter(OTU_tidy%in%nmr.filtered$top10.shared.asv.union$OTU_tidy) %>%
  # mutate(OTU_tidy_short=paste0(OTU," (",Genus,")"))%>%
  left_join(nmr.filtered$top10.shared.asv.union, by = "OTU_tidy")%>%
  left_join(sample.levels.filtered, by = "Sample")%>%
  mutate(Sample = factor (Sample, levels = sample.levels.filtered$Sample))%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU_tidy_short))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=otu.fill)+
  theme_bw()+
  coord_cartesian(expand = c(T,F,F,F))+
  labs(x="Sample",
       y="Relative abundance (%)",
       fill="ASV"
  )+
  project_theme +
  asv.barplot.theme +
  theme (plot.title = element_text(size = 8))

#+ fig.height = 6, fig.width = 8
print(top10.shared.asv.plot.filtered+
        ggtitle("Top 10 most abundant ASVs across age in filtered data")+
        theme(plot.title = element_text(size = 14)))

for(image.format in image.formats){
  ggsave(file = paste0("top10-asv-filtered.", image.format),
         plot = top10.shared.asv.plot.filtered,
         path = community.composition.figures,
         width = 8, height = 6,units = "in",
         dpi=300,device = image.format)
}

#+ echo=FALSE
## 7. Convert the filtered phyloseq object into a dataframe and save. ####
#'
#' ## Convert the filtered phyloseq object into a dataframe and save.
#' From now on, ps.q and ps.q.rel refer to the filtered dataset. We will delete
#' the raw data from the workspace.
ps.q <- ps.q.filtered
ps.q.rel <- ps.q.rel.filtered
ps.q.nmr <- ps.q.nmr.filtered
ps.q.nmr.rel <- ps.q.nmr.filtered.rel
sample.levels <- sample.levels.filtered
mean_sd_relab.nmr_ages.asv <- mean_sd_relab.nmr_ages.asv.filtered
mean_sd_relab.nmr_no_groups.asv <- mean_sd_relab.nmr_no_groups.asv.filtered

rm(ps.q.filtered, ps.q.rel.filtered)
rm(m40.asvs, m40.asv.plot, m40.asvs.uniq)
rm(ps.q.nmr.filtered, ps.q.nmr.filtered.rel)
rm(mean_sd_relab.nmr_ages.asv.filtered, mean_sd_relab.nmr_no_groups.asv.filtered)
rm(sample.levels.filtered)
rm(top10.asv.average, top10.asv.average.filtered,
   top10.shared.asv.plot, top10.shared.asv.plot.filtered,
   otu.fill)
#' Agglomerate the `ps.q.filtered` object and convert it into a dataframe with `psmelt()` 
#' function. Then, save it in tab-separated format and as an rds object. 
melt_phyloseq_and_save <-function (ps, ps.index){
  # ps.indexes <- c("OTU", "Phylum", "Family", "Genus")
  if(ps.index == "OTU"){
    ps.df <- ps.q %>%
      psmelt()
  }else{
    ps.df <- ps.q %>%
      tax_glom(ps.index, NArm=FALSE)%>%
      psmelt()
  }
  ps.df <- ps.df %>%
    dplyr::select(-sample_Sample)%>% # remove the duplicate column
    dplyr::select(-sex,-birthday,-animal)
  
  if(ps.index != "OTU"){
    ps.df<-ps.df%>%
      dplyr::select(-OTU)
  }
  # Save the tables in TSV format and as an RDS object
  ps.df.fname <- case_when(ps.index =="OTU" ~ ps.q.agg.asv.fname.no_ext,
                           ps.index =="Genus" ~ ps.q.agg.genus.fname.no_ext,
                           ps.index =="Family" ~ps.q.agg.family.fname.no_ext,
                           ps.index =="Phylum" ~ ps.q.agg.phylum.fname.no_ext
  ) 
  tsv.fname <- file.path(community.composition.tables, paste0(ps.df.fname,".tsv"))
  rda.fname <- file.path(community.composition.rdafiles, paste0(ps.df.fname,".rds"))
  if(!file.exists(tsv.fname)){
    write.table(ps.df,
                file = tsv.fname,
                row.names = F,sep = "\t")
  }
  if(!file.exists(rda.fname)){
    saveRDS(ps.df,
            file = rda.fname)
  }
  return(ps.df)
}

ps.q.agg.asv<-melt_phyloseq_and_save(ps.q, ps.index = "OTU")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.genus <- melt_phyloseq_and_save(ps.q, ps.index = "Genus")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.family <- melt_phyloseq_and_save(ps.q, ps.index = "Family")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.phylum<-melt_phyloseq_and_save(ps.q, ps.index = "Phylum")%>%
  filter(Abundance!=0)%>%
  as_tibble()
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-asv.tsv")
if(! file.exists (mean_sd_relab.nmr_ages.asv.fname)){
  write.table(mean_sd_relab.nmr_ages.asv,
              file = mean_sd_relab.nmr_ages.asv.fname,
              row.names = F,sep = "\t")
}

#' Save filtered phyloseq objects and the contents (OTU table, taxonomy, tree).
if(!file.exists(ps.q.filtered.fname)){
  saveRDS(ps.q,
          file = ps.q.filtered.fname)
}

if(!file.exists(ps.q.rel.filtered.fname)){
  saveRDS(ps.q.rel,
          file = ps.q.rel.filtered.fname)
}

if(!file.exists(ps.q.filtered.otu_table.fname)){
  saveRDS(ps.q@otu_table%>% as.matrix() %>%as.data.frame(),
          file = ps.q.filtered.otu_table.fname)
}
if(!file.exists(ps.q.filtered.tax_table.fname)){
  saveRDS(ps.q@tax_table%>% as.matrix() %>%as.data.frame(),
          file = ps.q.filtered.tax_table.fname)
}
if(!file.exists(ps.q.filtered.tree.fname)){
  saveRDS(ps.q@phy_tree,
          file = ps.q.filtered.tree.fname)
}

#+ echo=FALSE
## 8. Calculating summary statistics per host. ####
#'
#' ## Calculating summary statistics per host.
#' 
#+ echo=FALSE
#' The columns are:  
#' - Total reads  
#' - Library size (mean Abundance ± SD)  
#' - Number of ASV per host  
#' - Number of phyla per host  
#' - Number of families per host  
#' - Number of genera per host  
summary.stats<-create_summary_stats(ps.q.agg.asv,
                                                ps.q.agg.phylum,
                                                ps.q.agg.family,
                                                ps.q.agg.genus)
knitr::kable(as.data.frame(summary.stats),format = "simple")
summary.stats.fname<-file.path(community.composition.tables,
  "summary-table.tsv")
if(!file.exists(summary.stats.fname)){
  write.table(summary.stats,
              file=summary.stats.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 9. Calculate average, min, max, and SD relative abundances from phyloseq object. ####
#'
#' ## Calculate average, min, max, and SD relative abundances from phyloseq object.
mean_sd_relab.all_hosts.phyla <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Phylum")
mean_sd_relab.all_hosts.families <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Family")
mean_sd_relab.all_hosts.genera <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Genus")
mean_sd_relab.all_hosts.asv <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "OTU")

knitr::kable(head(mean_sd_relab.all_hosts.phyla),format = "simple")
knitr::kable(head(mean_sd_relab.all_hosts.families),format = "simple")
knitr::kable(head(mean_sd_relab.all_hosts.genera),format = "simple")
knitr::kable(head(mean_sd_relab.all_hosts.asv),format = "simple")

#' Save tables.
mean_sd_relab.all_hosts.phyla.fname <- 
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-all_hosts-Phylum.tsv")
if(!file.exists(mean_sd_relab.all_hosts.phyla.fname)){
  write.table(mean_sd_relab.all_hosts.phyla,
              file = mean_sd_relab.all_hosts.phyla.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.all_hosts.families.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-all_hosts-Family.tsv")
if(!file.exists (mean_sd_relab.all_hosts.families.fname)){
  write.table(mean_sd_relab.all_hosts.families,
              file = mean_sd_relab.all_hosts.families.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.all_hosts.genera.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-all_hosts-Genus.tsv")
if(! file.exists (mean_sd_relab.all_hosts.genera.fname)){
  write.table(mean_sd_relab.all_hosts.genera,
              file = mean_sd_relab.all_hosts.genera.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.all_hosts.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-all_hosts-OTU.tsv")
if(! file.exists (mean_sd_relab.all_hosts.asv.fname)){
  write.table(mean_sd_relab.all_hosts.asv,
              file = mean_sd_relab.all_hosts.asv.fname,
              row.names = F,sep = "\t")
}



#+ echo=FALSE
## 10. Calculate summary stats of unclassified genera in each animal (unrarefied). ####
#'
#' ## Calculate summary stats of unclassified genera in each animal (unrarefied). ####
#' Here, we are interested in unclassified genera, but you can also try Families,
#' Orders, Classes, etc.
unclassified.genus.summary.stats <- phyloseq_to_long(ps.q.rel,"Genus")%>%
  get_unclassified_summary_stats("Genus")%>%
  arrange(as.character(class))
knitr::kable(unclassified.genus.summary.stats, format ="simple")

unclassified.genus.summary.stats.fname <-
  file.path(community.composition.tables,
            "unclassified-genus-summary-table-raw.tsv")
if (! file.exists(unclassified.genus.summary.stats.fname)){
  write.table(unclassified.genus.summary.stats,
              file = unclassified.genus.summary.stats.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 11. Rarefy the table and check the percentage of unclassified taxa. ####
#'
#' ## Rarefy the table and check the percentage of unclassified taxa. ####
#' First, perform multiple rarefaction:
print(paste("Rarefying data with", rare.num_samples, "subsamples"))
ps.q.agg.asv.rare<-multi_rarefy_phyloseq(ps.q,
                                         tax.rank = "OTU",
                                         niter = rare.num_samples,
                                         output.filename.rds =  file.path(community.composition.rdafiles,paste0(
                                           paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
                                                 paste(custom.levels,collapse = '-'),sep = "-"),
                                           ".rds")))
ps.q.agg.phylum.rare<-multi_rarefy_phyloseq(ps.q,
                                            tax.rank = "Phylum",
                                            niter = rare.num_samples,
                                            output.filename.rds =  file.path(community.composition.rdafiles,paste0(
                                              paste("ps.q.rare","Phylum",paste0(rare.num_samples,"_iter"),
                                                    paste(custom.levels,collapse = '-'),sep = "-"),
                                              ".rds")))
ps.q.agg.family.rare<-multi_rarefy_phyloseq(ps.q,
                                            tax.rank = "Family",
                                            niter = rare.num_samples,
                                            output.filename.rds =  file.path(community.composition.rdafiles,paste0(
                                              paste("ps.q.rare","Family",paste0(rare.num_samples,"_iter"),
                                                    paste(custom.levels,collapse = '-'),sep = "-"),
                                              ".rds")))
ps.q.agg.genus.rare<-multi_rarefy_phyloseq(ps.q,
                                          tax.rank = "Genus",
                                          niter = rare.num_samples,
                                          output.filename.rds =  file.path(community.composition.rdafiles,paste0(
                                            paste("ps.q.rare","Genus",paste0(rare.num_samples,"_iter"),
                                                  paste(custom.levels,collapse = '-'),sep = "-"),
                                            ".rds")))
#' Rarefy naked mole-rat data specifically.
ps.q.agg.genus.nmr.rare<-multi_rarefy_phyloseq(subset_samples(ps.q, 
                                                              class =="NMR"),
                                               "Genus",
                                               niter = rare.num_samples,
                                               output.filename.rds =  file.path(community.composition.rdafiles,paste0(
                                                 paste("ps.q.rare","Genus",paste0(rare.num_samples,"_iter"),
                                                       paste("NMR",collapse = '-'),sep = "-"),
                                                 ".rds")))
ps.q.agg.asv.nmr.rare<-multi_rarefy_phyloseq(subset_samples(ps.q, class =="NMR"),
                                             "OTU",
                                             niter = rare.num_samples,
                                             output.filename.rds =  file.path(community.composition.rdafiles,paste0(
                                               paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
                                                     paste("NMR",collapse = '-'),sep = "-"),
                                               ".rds")))

#+ echo=FALSE
### 11.1 Add relative abundances and taxonomic information to the rarefied dataframe. ####
#' 
#' ### Add relative abundances and taxonomic information to the rarefied dataframe. 
#' All hosts (genus)
#' ps.q.agg.genus.rare.relab<-add_relab_to_tax_df(ps.q.agg.genus.rare,"Genus")
#' knitr::kable(head(ps.q.agg.genus.rare.relab),format = "simple")
#' 
#' #' Add other taxonomic ranks to the dataframe
#' ps.q.agg.genus.rare.relab<-ps.q.agg.genus.rare.relab%>%
#'   left_join(unique(ps.q.agg.genus[,c("Kingdom","Phylum","Class","Order","Family","Genus")]))
#' knitr::kable(head(ps.q.agg.genus.rare.relab),format = "simple")
#' 
#' #' Add relative abundances to NMR dataframe (ASV level)
#' ps.q.agg.asv.nmr.rare.relab<-add_relab_to_tax_df(ps.q.agg.asv.nmr.rare,"OTU")
#' ps.q.agg.asv.nmr.rare.relab<-ps.q.agg.asv.nmr.rare.relab%>%
#'   left_join(unique(ps.q.agg.asv[,c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")]))



#+ echo=FALSE
### 11.3 Create a summary stats table for the rarefied dataframe. ####
#' 
#' ### Create a summary stats table for the rarefied dataframe.
summary.stats.rare<-get_average_summary_stats(rare.asv = ps.q.agg.asv.rare,
                                                    rare.phylum = ps.q.agg.phylum.rare,
                                                    rare.family = ps.q.agg.family.rare,
                                                    rare.genus =  ps.q.agg.genus.rare,
                                   metadata = custom.md)
knitr::kable(summary.stats.rare, format ="simple")
summary.stats.rare.fname<-file.path(community.composition.tables,
                                     "summary-table-rarefied.tsv")
if(!file.exists(summary.stats.rare.fname)){
  write.table(summary.stats.rare,
              file=summary.stats.rare.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 12. Calculate summary stats of unclassified genera in rarefied data. ####
#' 
#' ## Calculate summary stats of unclassified genera in rarefied data.
unclassified.genus.summary.stats.rare<-
  get_average_unclassified_summary_stats_rare(ps.q.agg.genus.rare,"Genus",
                                              custom.md)%>%
  arrange(as.character(class))
knitr::kable(unclassified.genus.summary.stats.rare, format ="simple")
unclassified.genus.summary.stats.rare.fname<-
  file.path(community.composition.tables,
            paste("unclassified-genus-summary-table-rarefied.tsv",sep="-"))
if(!file.exists(unclassified.genus.summary.stats.rare.fname)){
  write.table(unclassified.genus.summary.stats.rare,
              file = unclassified.genus.summary.stats.rare.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 13. Analyse taxa in detail. ####
#'
#' ## 13. Analyse taxa in detail. ####

#+ echo=FALSE
### 13.1 Check how much Bacteroidaceae are in NMR.  ####
#'
#' ### Check how much Bacteroidaceae and Bacteroidota are in NMR.
bacteroidaceae.nmr<-mean_sd_relab.all_hosts.families%>%
  filter(class=="NMR")%>%
  filter(Family=="Bacteroidaceae")
knitr::kable(bacteroidaceae.nmr, format = "simple")
bacteroidaceae.nmr.fname<- file.path(community.composition.tables,
                                     "bacteroidaceae-table-nmr.tsv")
if(! file.exists (bacteroidaceae.nmr.fname)){
  write.table(bacteroidaceae.nmr,
              file = bacteroidaceae.nmr.fname,
              row.names = F,sep = "\t")
}

bacteroidota.nmr<-mean_sd_relab.all_hosts.families%>%
  filter(class=="NMR")%>%
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
### 13.2 Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. ####
#'
#' ### Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. 
spirochaetaceae.nmr<-mean_sd_relab.all_hosts.families%>%
  filter(class=="NMR")%>%
  filter(Family=="Spirochaetaceae")
knitr::kable(spirochaetaceae.nmr, format = "simple")

spirochaetota.nmr<-mean_sd_relab.all_hosts.phyla%>%
  filter(class=="NMR")%>%
  filter(Phylum=="Spirochaetota")
knitr::kable(spirochaetota.nmr, format = "simple")
spirochaetaceae.nmr.fname <-file.path(community.composition.tables,
                                      "spirochaetaceae-table-nmr.tsv")
if(! file.exists (spirochaetaceae.nmr.fname)){
  write.table(spirochaetaceae.nmr,
              file = spirochaetaceae.nmr.fname,
              row.names = F,sep = "\t")
  
}

treponema.nmr<-mean_sd_relab.all_hosts.genera%>%
  filter(class=="NMR")%>%
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
#### 13.2.1 Check the number of ASVs in Treponema from NMR. ####
#' 
#' #### Check the number of ASVs in Treponema from NMR.
ps.q.agg.asv%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

#+ echo=FALSE
#### 13.2.2 Plot Treponema and other Spirochaetota in all hosts. ####
#'
#' #### Plot Treponema and other Spirochaetota in all hosts.
spirochaetota.plot<-ps.q.agg.genus%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance = 100*Abundance/sum(Abundance))%>%
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
### 13.3 Analyse sulfur-metabolising bacteria. ####
#'
#' ### Analyse sulfur-metabolising bacteria. 
#' Desulfobacterota in NMR:
desulfobacterota.nmr<-mean_sd_relab.all_hosts.genera%>%
  filter(Phylum=="Desulfobacterota",class=="NMR")
knitr::kable(head(desulfobacterota.nmr), format = "simple")

#' Desulfobacterota in all hosts:
desulfobacterota.all<-mean_sd_relab.all_hosts.genera%>%
  filter(Phylum=="Desulfobacterota")%>%
  arrange(desc(MeanRelativeAbundance))
knitr::kable(head(desulfobacterota.all), format = "simple")

#+ echo=FALSE
#### 13.3.1 Total Desulfobacterota in NMR. ####
#' 
#' #### Total Desulfobacterota in NMR.
mean_sd_relab.all_hosts.genera%>%
  filter(Phylum=="Desulfobacterota", class=="NMR")%>%
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
#### 13.3.2 Plot Desulfobacterota in all hosts. ####
#'
#' #### Plot Desulfobacterota in all hosts.
# It's a flipped plot
desulfobacterota.plot<-ps.q.agg.genus%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance = 100*Abundance/sum(Abundance))%>%
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
## 15. Check non-bacterial data. ####
#' 
#' ## Check non-bacterial data. ####
#' Show the hosts and their non-bacterial taxa
mean_sd_relab.all_hosts.genera %>%
  filter(Kingdom != "Bacteria")%>%
  knitr::kable(format = "simple")

ps.q.agg.asv%>%
  filter(Order == "Methanomassiliicoccales")%>%
  knitr::kable(format = "simple")

#' Non-bacterial taxa plot:
ps.q.agg.asv%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance = 100*Abundance/sum(Abundance))%>%
  filter(Kingdom!="Bacteria")%>%
  mutate(class=factor(class,levels=custom.levels))%>%
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
#' Very few taxa were detected.

sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()