#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---
#' ```{r, setup 001-phyloseq-qiime2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank ####
#' # Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank
#'  
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' Once you produce the feature table, taxonomic classification, 
#' and the phylogenetic tree in QIIME2, it's time to perform 
#' downstream processing in R. First, we need to import the 
#' QZA files using `qiime2R` package.
#' We will convert the QZA files directly into phyloseq objects.
#'
#' The final output is:
#' 1) `phyloseq` object `ps.q`, which contains the QIIME2 
#' output (feature table, taxonomy, and tree). Saved as an RDS object.
#' 
#' 2) `phyloseq` object `ps.q.rel`, which is same as `ps.q`, but absolute 
#' abundance values are transformed into relative abundances.
#' Saved as an RDS object.
#' 
#' 3) Metadata `custom.md` and `custom.md.ages`: these are tables with sample
#' metadata. The first file contains information for all samples, while the 
#' second one only for naked mole-rat samples. `custom.md` will be added 
#' to `ps.q` as a `sample_data`. Saved as an RDS object and a TSV file.
#' 
#'
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
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

#+ echo=FALSE
## 2. Import data from QIIME2. #### 
#'
#' ## Import data from QIIME2.  
#' 
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
print(paste("The analysis focus is:", active.analysis))
source(here::here("config/R/config.R"))# config file with global variables

dir.create(community.composition.tables,recursive = TRUE)
dir.create(community.composition.rdafiles,recursive = TRUE)
dir.create(community.composition.figures,recursive = TRUE)

#+ echo=FALSE
## 3. Import qza files and convert them into a phyloseq object. ####
#'
#' ## Import qza files and convert them into a phyloseq object.
ps.q<-qza_to_phyloseq(
  features = file.path(qiime2.output.dir, "trimmed-dada2-table-merged-filtered.qza"), # feature table
  taxonomy = file.path(qiime2.output.dir, "trimmed-dada2-merged-taxonomy.qza"), # taxonomy
  tree = file.path(qiime2.output.dir, "trimmed-dada2-merged-filtered-rooted-tree.qza") # rooted tree
)

#' Change the name d__Kingdom to Kingdom.
ps.q.taxtab<-as.data.frame(tax_table(ps.q))
ps.q.taxtab$Kingdom<-
  gsub("d__","",ps.q.taxtab$Kingdom)
tax_table(ps.q)<-as.matrix(ps.q.taxtab)
rm(ps.q.taxtab)

#+ echo=FALSE
## 4. Add custom metadata. ####
#'
#' ## Add custom metadata. 
custom.md<-read.table(qiime2.metadata.filename, header = T,sep = "\t")
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
#' Convert the Sample column into row names because phyloseq 
#' needs samples as rownames.
#' 
#' Remove absolute.filepath column.
custom.md<-custom.md%>%
  dplyr::select(-forward.absolute.filepath,
                -reverse.absolute.filepath)%>%
  mutate(class= as.factor(class),
         sex = as.factor(sex),
         birthday = as.Date(birthday),
         animal = as.factor(animal))
# birthday=as.Date(ifelse(class=="B6mouse",sampling_date-weeks(12),birthday)))
rownames(custom.md)<-custom.md$Sample

#' Assign the custom metadata as your phyloseq object's metadata.
sample_data(ps.q)<-custom.md

#+ echo=FALSE
### 4.1 For NMR, we create a separate metadata object with age groups. ####
#'
#' ### For NMR, we create a separate metadata object with age groups. ####
custom.md.ages<-custom.md%>%
  dplyr::select(Sample,class)%>%
  filter(class=="NMR")%>%
  left_join(biosample.md,by = join_by("Sample" =="host_subject_id"))%>%
  rename("sex"= host_sex,
         "age" = host_age)%>%
  dplyr::select(-organism,-env_broad_scale,
         -env_local_scale, -env_medium,
         -geo_loc_name, -lat_lon,
         -host)%>%
  mutate(age = as.numeric(age),
         agegroup=cut(age, breaks =c(0,10,16),
                      right = FALSE))%>%
  ungroup()
#' We create these new levels for differential microbial abundance.
unique_levels <- custom.md.ages %>%
  distinct(agegroup)%>%
  arrange(agegroup) %>%
  mutate(new_agegroup = paste0("agegroup", agegroup))%>%
  mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
  mutate(new_agegroup = gsub("\\,","_",new_agegroup))
custom.md.ages <- custom.md.ages %>%
  left_join(unique_levels, by = "agegroup")%>%
  rename("old_agegroup" = "agegroup")%>%
  rename("agegroup" = "new_agegroup")%>%
  dplyr::select(-old_agegroup)%>%
  as.data.frame()%>%
  column_to_rownames("sample_name")

custom.md.ages.fname.no_ext <- file.path(new.metadata.dir,"custom.md.ages")
custom.md.ages.fname <- paste(custom.md.ages.fname.no_ext, "rds", sep = ".")
if(!file.exists(custom.md.ages.fname)){
  saveRDS(custom.md.ages,file = custom.md.ages.fname)
  write.table(custom.md.ages,file=paste(custom.md.ages.fname.no_ext, "tsv", sep = "."),
              row.names = F,sep = "\t")
}

#+ echo=FALSE
### 4.2 Filter metadata. ####
#'
#' ### Filter metadata. ####
#' You can exclude samples based on their library size (total number of reads).
low.abundance.samples <- names(which(sample_sums(ps.q)<20000))
print("Samples with library size < 20000:")
print(low.abundance.samples)

ps.q <- subset_samples(ps.q, ! sample_names(ps.q) %in% low.abundance.samples)

#' You can exclude samples based config file (they failed in QIIME2).
if( "excluded.classes" %in% ls(all.names = T) ){
  print("Excluded hosts:")
  print(excluded.classes)
  if( excluded.classes %in% sample_data(ps.q)$class ){
    ps.q <- subset_samples(ps.q, ! class %in% excluded.classes)
  }
}

#' Save the sample data of the ps.q object as custom.md
custom.md <- sample_data(ps.q)%>%
  as.matrix()%>%
  as.data.frame()%>%
  mutate(class = factor(class,  levels = unique(class)),
         sex = as.factor(sex),
         birthday = as.Date(birthday),
         animal = as.factor(animal))
custom.md.fname.no_ext <- file.path(new.metadata.dir,"custom.md")
custom.md.fname <- paste(custom.md.fname.no_ext, "rds", sep = ".")
if (!file.exists (custom.md.fname)){
  saveRDS(custom.md,file = custom.md.fname)
  write.table(custom.md,file = paste(custom.md.fname.no_ext, "tsv", sep = "."),
              row.names = F,sep = "\t")
}

#' Number of features in the unfiltered dataset:
ntaxa(ps.q)

#' Total frequency in the unfiltered dataset:
sum(sample_sums(ps.q))

#' Summary statistics (min, median, max, quartiles) of the unfiltered dataset:
summary(sample_sums(ps.q))

#' Select only Bacteria. Remove chloroplast and mitochondria
if(active.analysis == "bacteria_only"){
  ps.q<-ps.q %>%
    subset_taxa(Kingdom%in%"Bacteria")
}
ps.q <- ps.q%>%
  subset_taxa(!Order %in% "Chloroplast")%>%
  subset_taxa(!Family %in% "Mitochondria")
#+ echo=FALSE
### 4.3 Fix empty taxa with higher rank taxon. ####
#'
#' ### Fix empty taxa with higher rank taxon. ####
#' Because we want to remove NA values and make ambiguous "uncultured" or 
#' "unclassified" taxa more understandable.
ps.q<-tax_fix(ps.q,unknowns = c("NA","uncultured",#"Unassigned",
                                "uncultured_bacterium","uncultured_rumen",
                                "gut_metagenome","human_gut","mouse_gut",
                                "wallaby_gut","uncultured_soil", 
                                "uncultured_organism","uncultured_prokaryote",NA))
#+ echo=FALSE
## 5. Calculate relative abundance: Abundance divided by total abundance in a sample. ####
#'
#' ## 5. Calculate relative abundance: Abundance divided by total abundance in a sample. ####
ps.q.rel <- transform_sample_counts(ps.q, function(x) 100*x/sum(x)) 

#' Save phyloseq objects and the contents (OTU table, taxonomy, tree).
if(!file.exists(ps.q.raw.fname)){
  saveRDS(ps.q,
          file = ps.q.raw.fname)
}

if(!file.exists(ps.q.rel.raw.fname)){
  saveRDS(ps.q.rel,
          file = ps.q.rel.raw.fname)
}

if(!file.exists(ps.q.raw.otu_table.fname)){
  saveRDS(ps.q@otu_table%>% as.matrix() %>%as.data.frame(),
          file = ps.q.raw.otu_table.fname)
}
if(!file.exists(ps.q.raw.tax_table.fname)){
  saveRDS(ps.q@tax_table%>% as.matrix() %>%as.data.frame(),
          file = ps.q.raw.tax_table.fname)
}
if(!file.exists(ps.q.raw.tree.fname)){
  saveRDS(ps.q@phy_tree,
          file = ps.q.raw.tree.fname)
}


sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()