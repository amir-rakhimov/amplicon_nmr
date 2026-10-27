#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---

#' ```{r, setup 007-diffabund-tests.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# Differential microbial abundance tests with MaAsLin2, ALDEx2, and ANCOM-BC ####
#' # Differential microbial abundance tests with MaAsLin2, ALDEx2, and ANCOM-BC
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script performs differential microbial abundance tests on 
#' different hosts. We will use MaAsLin2, ALDEx2, and ANCOM-BC.

#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# Current: Bioconductor version 3.20 (BiocManager 1.30.27), R 4.4.3 (2025-02-28 ucrt)
#' ALDEx2: 1.38.0
#' Maaslin2: 1.20.0
#' ANCOM-BC: 2.8.1  
# BiocManager::install(c("Maaslin2","ALDEx2","ANCOMBC","phyloseq"), version = "3.20")
# BiocManager::install("ALDEx2", version = "3.17")
# BiocManager::install("ANCOMBC", version = "3.17")
# BiocManager::install("phyloseq", version = "3.17")
# install.packages(c("tidyverse","Polychrome"))
library(tidyverse)
library(Maaslin2)
library(ALDEx2)
library(ANCOMBC)
library(phyloseq)
library(gplots)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
rare.status<-"unrare"
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
message("The analysis focus is: ", active.analysis)

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

dir.create(diffabund.tables,recursive = TRUE)
dir.create(diffabund.rdafiles,recursive = TRUE)
dir.create(diffabund.figures,recursive = TRUE)

#' Import the phyloseq object
ps.q<-readRDS(file= ps.q.filtered.fname)
ps.q.rel<-readRDS(file = ps.q.rel.filtered.fname)

#+ echo=FALSE
## 3. Prepare necessary functions. ####
#'
#' ## Prepare necessary functions. ####
#' Function that performs differential abundance tests with different 
#' parameters. Instead of repeating the code for each comparison,
#' call the function and change the parameters.

### Create reference levels. ####
prepare_data <- function(ps,
                         agglom.rank, comparison, 
                         ref.level, metadata,
                         relations = NULL){
  message("Preparing data for differential microbial abundance tests on ", comparison, " variable at ", agglom.rank, " level. ",
          "Reference: ", ref.level)
  comparison_meta <- list(
    host = list(custom.levels = intersect(names(pretty.level.names),metadata$class),
                maaslin2.reference = paste0("class,",ref.level),
                maaslin2.fixed_effects = "class",
                maaslin2.random_effects = NULL,
                maaslin2.random_effects.file_suffix = "" ,
                aldex.reference = paste0("covariates",ref.level),
                ancombc.comparison = "class",
                ggplot.levels = intersect(names(pretty.level.names),metadata$class),
                ggplot.labels = pretty.level.names,
                gg.labs.name = "Host",
                gg.title.groups = "hosts",
                output.filename = paste("all_hosts",agglom.rank,comparison,
                                        "ref",ref.level,
                                        sep="-")
    ),
    age = list(custom.levels = "NMR",
               maaslin2.reference = paste("agegroup", ref.level, sep = ","),
               maaslin2.fixed_effects = "agegroup",
               maaslin2.random_effects = c("relation"),
               maaslin2.random_effects.file_suffix = "with-relations",
               aldex.reference = paste0("covariates",ref.level), 
               ancombc.comparison = "agegroup",
               ggplot.levels = c("agegroup0_10",
                                 "agegroup10_16" ),
               ggplot.labels = c("agegroup0_10" = "agegroup0_10",
                                 "agegroup10_16" = "agegroup10_16"),
               gg.labs.name = "Age group",
               gg.title.groups = "age groups",
               output.filename = paste("NMR",agglom.rank,comparison,
                                       "ref",ref.level,
                                       sep="-")
    ),
    sex = list(custom.levels = "NMR",
               maaslin2.reference = paste("sex", ref.level, sep = ","),
               maaslin2.fixed_effects = "sex" ,
               maaslin2.random_effects = c("relation"),
               maaslin2.random_effects.file_suffix = "with-relations",
               aldex.reference = paste0("covariates",ref.level),# M or F
               ancombc.comparison = "sex",
               ggplot.levels = c("female" ,
                                 "male" ),
               ggplot.labels = c("female" = "Females",
                                 "male" = "Males"),
               gg.labs.name = "Host sex",
               gg.title.groups = "groups",
               output.filename = paste("NMR",agglom.rank,comparison,
                                       "ref",ref.level,
                                       sep="-")
    )
  )
  current.comparison_meta <- comparison_meta[[comparison]]
  ### Add relations metadata ####
  if (setequal(current.comparison_meta$custom.levels,"NMR")){
    metadata<-metadata %>%
      right_join(relations,by="Sample")
    rownames(metadata)<-metadata$Sample
  }
  # add metadata to the list
  current.comparison_meta$metadata <- metadata
  
  # add sample groups variable
  if(setequal(current.comparison_meta$custom.levels,"NMR")){
    if(comparison=="age"){
      current.comparison_meta$sample.groups <-c("agegroup0_10",
                                                "agegroup10_16")
    }else if(comparison=="sex"){
      current.comparison_meta$sample.groups<-c("female",
                                               "male")
    }
  }else{
    current.comparison_meta$sample.groups <- current.comparison_meta$custom.levels
  }
  
  #' Relevel the comparison vector for ANCOM-BC. The first level will be the reference 
  #' for custom leveling
  current.comparison_meta$ancombc.levels <- 
    c(ref.level, 
      current.comparison_meta$sample.groups[current.comparison_meta$sample.groups!=ref.level])
  
  # Create a phyloseq with relative abundances
  ps.rel <- transform_sample_counts(ps, function(x) 100*x/sum(x))
  ### Melt the phyloseq object. ####
  if(agglom.rank == "OTU"){
    ps.agg <- ps%>%
      psmelt()
    
    ps.rel.agg <- ps.rel%>%
      psmelt()
  }else{
    ps.agg <- ps%>%
      tax_glom(taxrank = agglom.rank)%>%
      psmelt()%>%
      dplyr::select(-OTU)
    ps.rel.agg <- ps.rel%>%
      tax_glom(taxrank = agglom.rank)%>%
      psmelt()%>%
      dplyr::select(-OTU)
  }
  
  ps.agg.wide <- ps.agg %>%
    dplyr::select(all_of(c("Sample", "Abundance", agglom.rank)))%>%
    pivot_wider(names_from = agglom.rank, # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("Sample")
  
  current.comparison_meta$ps.agg <- ps.agg
  current.comparison_meta$ps.rel.agg <- ps.rel.agg
  current.comparison_meta$ps.agg.wide <- ps.agg.wide
  
  
  #' ### Create covariates for ALDEX2: take the samples in the ps.agg.wide, and find 
  #' the positions of these samples in the metadata. 
  #' Ref: match returns a vector of the positions of (first) matches of its first argument in its second. 
  #' We're searching the elements of the first vector in the second vector. 
  current.comparison_meta$covariates <-
    current.comparison_meta$metadata[[current.comparison_meta$ancombc.comparison]] [ # ancombc.comparison is same as for aldex
      match(rownames(current.comparison_meta$ps.agg.wide),
            rownames(current.comparison_meta$metadata))
    ]
  # if(setequal(current.comparison_meta$sample.groups,"NMR")& comparison=="age"){
  #   current.comparison_meta$covariates<-
  #     current.comparison_meta$metadata$agegroup[
  #       match(rownames(current.comparison_meta$ps.agg.wide),
  #             rownames(current.comparison_meta$metadata))]
  # }else if (comparison == "host"){
  #   current.comparison_meta$covariates<-
  #     current.comparison_meta$metadata$class[
  #       match(rownames(current.comparison_meta$ps.agg.wide),
  #             rownames(current.comparison_meta$metadata))]
  # }
  #' Model matrix is covariates and the reference group is the first covariate.
  mm <- model.matrix(~ current.comparison_meta$covariates-1)
  #' Reorder model.matrix columns to put aldex.reference as first column and
  # all the other columns after. 
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == current.comparison_meta$aldex.reference), 
                          which(colnames(mm) !=current.comparison_meta$aldex.reference))]]
  current.comparison_meta$mm <- mm
  
  return(current.comparison_meta)
}


#' Function that performs a test with MaAsLin2:
perform_maaslin2_test <- function(ps.agg.wide, 
                                  output.filename,
                                  comparison,
                                  metadata,
                                  ref.level,
                                  maaslin2.random_effects,
                                  maaslin2.random_effects.file_suffix,
                                  maaslin2.fixed_effects,
                                  maaslin2.reference
){
  #' Maaslin2 test
  maaslin2.main.dir <- file.path(maaslin2.out,
                                 paste0(authorname,"-output"),
                                 rare.status)
  if(!dir.exists(maaslin2.main.dir)){
    dir.create(maaslin2.main.dir,recursive = TRUE)
  }
  maaslin2.output.path <- file.path(maaslin2.out,
                                    paste0(authorname,"-output"),
                                    rare.status,paste(
                                      output.filename, 
                                      maaslin2.random_effects.file_suffix,
                                      sep = "-"))
  fit.rds.fname <- file.path(diffabund.rdafiles,
                             paste(
                               paste("maaslin.fit_data",comparison, sep = "."),
                               paste0(output.filename,".rds"),sep = "-"))
  
  if(dir.exists(maaslin2.output.path) & file.exists(fit.rds.fname)){
    message(maaslin2.output.path, " exists. Skipping")
    maaslin.fit_data <- readRDS(fit.rds.fname)
  }else{
    set.seed(1)
    maaslin.fit_data = 
      Maaslin2(input_data = ps.agg.wide, 
               input_metadata = metadata, 
               min_prevalence = 0,
               normalization = "TSS",
               transform = "LOG",
               analysis_method = "LM",
               random_effects = maaslin2.random_effects, 
               standardize = FALSE,
               output = maaslin2.output.path, 
               fixed_effects = c(maaslin2.fixed_effects),
               reference = maaslin2.reference,
               max_significance = 0.05,
               save_models = TRUE)
    #' Save the fit data object as an rds file.
    if(!file.exists(fit.rds.fname)){
      saveRDS(maaslin.fit_data,
              file = fit.rds.fname)
      message("Saved maaslin.fit_data as ", fit.rds.fname)
    }
    
  }
  return(maaslin.fit_data)
}


#' Function to process the MaAsLin2 output:
process_maaslin2_output<-function(maaslin.fit_data,
                                  ps.agg,
                                  sample.groups,
                                  agglom.rank,
                                  output.filename){
  #' Downstream processing of Maaslin2 output.
  source(file.path(util.functions.r,"make_features_maaslin.R"))
  #' Extract features with qvalue < 0.05.
  if(min(maaslin.fit_data$results$qval)<0.05){
    maaslin2.signif.features<-maaslin.fit_data$results%>%
      filter(qval<0.05) # should be qval
  }else{
    maaslin2.signif.features<-maaslin.fit_data$results%>%
      arrange(qval)%>%
      head(n = 15) # if no significant results found
  }
  
  #' Make features pretty. We have to make this exchange because maaslin 
  #' output treats space and hyphen as the same thing.
  if(agglom.rank=="OTU"){
    maaslin2.signif.features$feature<-gsub("^X","",
                                           maaslin2.signif.features$feature)
  }else{
    foo<-ps.agg
    foo<-foo%>%mutate("maaslin"=get(agglom.rank))
    foo<-make_features_maaslin(foo,"maaslin")
    foo<-unique(foo[,c("maaslin",agglom.rank)])
    maaslin2.signif.features<-maaslin2.signif.features%>%
      left_join(foo[,c("maaslin",agglom.rank)],by=c("feature"="maaslin"))%>%
      distinct()
    rm(foo)
    maaslin2.signif.features$feature<-maaslin2.signif.features[,agglom.rank]
    maaslin2.signif.features<-subset(maaslin2.signif.features, 
                                     select=-get(agglom.rank))
  }
  
  #' Extract features that are downregulated in all hosts simultaneously.
  # Downreglated features have coef<0.
  # We also add association strength column
  maaslin2.signif.decreased<-maaslin2.signif.features%>%
    as_tibble()%>%
    filter(coef<0)%>%
    group_by(feature)%>%
    filter(n()==length(sample.groups)-1)%>% 
    arrange(pval,feature)%>%
    mutate(n=n())%>%
    mutate(assoc.str=-log(qval)*sign(coef))%>%
    ungroup()
  
  
  #' Check if all features (OTU/taxa/etc) are also present in the 
  #' ps.agg dataframe.
  table(maaslin2.signif.decreased$feature%in%ps.agg[,agglom.rank])%>%
    print()
  
  signif.features.fname <-file.path(diffabund.tables,paste(
    "maaslin2",output.filename,
    "signif.tsv",sep="-"))
  signif.decreased.fname <- file.path(diffabund.tables,
                                      paste("maaslin2.signif.decreased",
                                            output.filename,
                                            "signif.tsv",sep="-"))
  
  if(file.exists(signif.features.fname)){
    message(signif.features.fname, " exists. Skipping")
  }else{
    message("Saving maaslin2.signif.features as ", signif.features.fname)
    write.table(maaslin2.signif.features,
                file = signif.features.fname,
                row.names = F,sep = "\t")
    
  }
  if(file.exists(signif.decreased.fname)){
    message(signif.decreased.fname, " exists. Skipping")
  }else{
    message("Saving maaslin2.signif.decreased as ", signif.decreased.fname)
    write.table(maaslin2.signif.decreased,
                file = signif.decreased.fname,
                row.names = F,sep = "\t")
    
  }
  return(list(maaslin2.signif.features = maaslin2.signif.features,
              maaslin2.signif.decreased = maaslin2.signif.decreased))
  
}


#' Function that performs a test with ANCOM-BC:
perform_ancombc_test<-function(ps.agg.wide,
                               metadata, 
                               ancombc.comparison,
                               ancombc.levels,
                               output.filename){
  
  metadata[[ancombc.comparison]] <- factor(metadata[[ancombc.comparison]],
                                           levels = ancombc.levels)
  # Main output
  ancombc.out.fname <- file.path(diffabund.rdafiles,
                                 paste("ancombc.out",output.filename,
                                       ".rds",sep="-"))
  # All significant features
  ancombc.signif.features.fname <- file.path(diffabund.tables,paste(
    "ancombc",output.filename,
    "signif.tsv",sep="-"))
  # Significantly decreased features.
  ancombc.signif.decreased.fname <-file.path(diffabund.tables,
                                             paste("ancombc.signif.decreased",
                                                   output.filename,
                                                   "signif.tsv",sep="-"))
  
  if(file.exists(ancombc.signif.features.fname) & 
     file.exists(ancombc.signif.decreased.fname) &
     file.exists(ancombc.out.fname)){
    message("ANCOM-BC output exists. skipping")
    ancombc.out <-readRDS(ancombc.out.fname)
    ancombc.signif.features <-read.table(ancombc.signif.features.fname,
                                         header = TRUE,
                                         sep = "\t")
    ancombc.signif.decreased <-read.table(ancombc.signif.decreased.fname,
                                          header = TRUE,
                                          sep = "\t")
  }else{
    #' Perform the differential abundance test.
    set.seed(1)
    ancombc.out<-ancombc(
      data = ps.agg.wide,
      meta_data = metadata,
      taxa_are_rows = FALSE,
      # tax_level = ,
      formula = ancombc.comparison,
      p_adj_method = "fdr", 
      prv_cut = 0, # by default prevalence filter of 10% is applied
      lib_cut = 0, 
      group = ancombc.comparison,
      struc_zero = TRUE, 
      neg_lb = TRUE, 
      tol = 1e-5, 
      max_iter = 100, 
      conserve = TRUE, 
      alpha = 0.05, 
      global = TRUE
    )
    
    if(!file.exists(ancombc.out.fname)){
      saveRDS(ancombc.out,
              file = ancombc.out.fname)
    }
    ancombc.res<-ancombc.out$res
    #' Find differentially abundant taxa by multiplying fold change with TRUE/FALSE.
    df_lfc = data.frame(ancombc.res$lfc[, -1] * ancombc.res$diff_abn[, -1], 
                        check.names = FALSE) %>%
      mutate(taxon_id = ancombc.res$diff_abn$taxon) %>%
      dplyr::select(taxon_id, everything())
    ancombc.signif.features<-subset(df_lfc,
                                    rowSums(df_lfc[,-c(1,2)] != 0) ==
                                      ncol(df_lfc[,-c(1,2)]))
    ancombc.signif.decreased<-subset(ancombc.signif.features,
                                     rowSums(ancombc.signif.features[,-c(1,2)] < 0)==
                                       ncol(ancombc.signif.features[,-c(1,2)]))
    head(ancombc.signif.features)%>%
      knitr::kable(format = "simple")%>%
      print()
    head(ancombc.signif.decreased)%>%
      knitr::kable(format = "simple")%>%
      print()
    
    #' Save the workspace.
    ancombc.workspace.fname <-file.path(diffabund.rdafiles,paste(
      "ancombc", output.filename, "workspace.RData",sep="-"))
    if(!file.exists(ancombc.workspace.fname)){
      save.image(ancombc.workspace.fname)
    }
    #' Save all significant features. 
    # "20260213_13_41_58"
    if(!file.exists(ancombc.signif.features.fname)){
      write.table(ancombc.signif.features,
                  file = ancombc.signif.features.fname,
                  row.names = F,sep = "\t")
    }
    #' Save decreased significant features.
    # "20260213_13_42_00"
    if(!file.exists(ancombc.signif.decreased.fname)){
      write.table(ancombc.signif.decreased,
                  file = ancombc.signif.decreased.fname,
                  row.names = F,sep = "\t")
    }
  }
  
  return(list(ancombc.out = ancombc.out,
              ancombc.signif.features = ancombc.signif.features,
              ancombc.signif.decreased = ancombc.signif.decreased))
}

perform_aldex2_test <-function(comparison,
                               ps.agg.wide,
                               metadata, 
                               mm,
                               sample.groups,
                               aldex.reference,
                               output.filename){
  #' Check if ALDEx2 output already exists.
  ps.q.glm.effect.fname <- file.path(diffabund.rdafiles,
                                     paste("ps.q.glm.effect",
                                           output.filename,
                                           ".rds",sep="-"))
  ps.q.aldex.clr.fname <- file.path(diffabund.rdafiles,
                                    paste("ps.q.aldex.clr",
                                          output.filename,
                                          ".rds",sep="-"))
  ps.q.glm.test.fname <- file.path(diffabund.rdafiles,
                                   paste("ps.q.glm.test",
                                         output.filename,
                                         ".rds",sep="-"))
  
  if (file.exists(ps.q.glm.effect.fname) & file.exists(ps.q.aldex.clr.fname) &
      file.exists(ps.q.glm.test.fname)){
    message("ALDEx2 output exists. skipping")
    ps.q.glm.effect <- readRDS(ps.q.glm.effect.fname)
    ps.q.aldex.clr <- readRDS(ps.q.aldex.clr.fname)
    ps.q.glm.test <- readRDS(ps.q.glm.test.fname)
  }else{
    #' Run the test: Aldex glm for a complex case.
    set.seed(1)
    ps.q.aldex.clr <- aldex.clr(t(ps.agg.wide), mm, 
                                mc.samples=1000, 
                                denom="all", verbose=F)
    ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
    #' Save the workspace just in case.
    glm.test.workspace.fname <-file.path(diffabund.rdafiles,
                                         paste("aldex2",output.filename,
                                               "workspace-test.RData",sep="-"))
    if(!file.exists(glm.test.workspace.fname)){
      save.image(glm.test.workspace.fname)
      
    }
    ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
    #' Save the workspace just in case.
    glm.effect.workspace.fname <-file.path(diffabund.rdafiles,
                                           paste("aldex2",output.filename,
                                                 "workspace-effect.RData",sep="-"))
    if(!file.exists(glm.effect.workspace.fname)){
      save.image(glm.effect.workspace.fname)
    }
    if(!file.exists(ps.q.aldex.clr.fname)){
      saveRDS(ps.q.aldex.clr,
              file = ps.q.aldex.clr.fname)
      message("Saved ps.q.aldex.clr as ", ps.q.aldex.clr.fname)
    }
    if(!file.exists(ps.q.glm.test.fname)){
      saveRDS(ps.q.glm.test,
              file = ps.q.glm.test.fname)
      message("Saved ps.q.glm.test as ", ps.q.glm.test.fname)
    }
    if(!file.exists(ps.q.glm.effect.fname)){
      saveRDS(ps.q.glm.effect,
              file = ps.q.glm.effect.fname)
      message("Saved ps.q.glm.effect as ", ps.q.glm.effect.fname)
    }
    
  }
  
  #' Check if processed ALDEx2 output already exists.
  aldex2.neg.effect.fname <- file.path(diffabund.tables,
                                       paste("aldex2.neg.effect",
                                             output.filename,
                                             "signif.tsv",sep="-"))
  aldex2.signif.features.fname <- file.path(diffabund.tables,paste(
    "aldex2",output.filename,"signif.tsv",sep="-")) 
  
  if (file.exists (aldex2.neg.effect.fname) & file.exists(aldex2.signif.features.fname)){
    message("Data with significant features from ALDEx2 exists. Skipping")
    aldex2.neg.effect <- read.table(aldex2.neg.effect.fname,header = TRUE)
    aldex2.signif.features <- read.table(aldex2.signif.features.fname, header = TRUE)
  }else{
    #' Extract significant features.
    aldex2.signif.features<-list()
    for (i in 1:length(ps.q.glm.effect)){
      # take all features that have good CI (not overlapping zero) and high effect size
      # identify features with significant effect size and good CI
      sig<-which((ps.q.glm.effect[[i]]$effect.low>0 & 
                    ps.q.glm.effect[[i]]$effect.high>0)|
                   (ps.q.glm.effect[[i]]$effect.low<0 & 
                      ps.q.glm.effect[[i]]$effect.high<0))
      sig<-intersect(which(abs(ps.q.glm.effect[[i]]$effect)>1), sig)
      
      if(length(sig)!=0){
        signif.df<-ps.q.glm.effect[[i]][sig,]
        
        if(comparison == "age" | comparison == "sex"){
          signif.df$OTU<-rownames(signif.df)
        }else{
          signif.df$Taxon<-rownames(signif.df)
        }
        
        # signif.df$class<-names(ps.q.glm.effect[i])
        rownames(signif.df)<-1:nrow(signif.df)
        aldex2.signif.features[[names(ps.q.glm.effect[i])]]<-signif.df
      } else{
        signif.df<-data.frame()
        aldex2.signif.features[[names(ps.q.glm.effect[i])]]<-signif.df
      }
      
    }
    rm(signif.df)
    aldex2.signif.features<-bind_rows(aldex2.signif.features,.id = "class")
    #' Check for duplicates.
    if(comparison == "age" | comparison == "sex"){
      aldex2.signif.features%>%
        dplyr::group_by(OTU)%>%
        summarise(n=n())%>%
        arrange(-n) 
    }else{
      aldex2.signif.features%>%
        dplyr::group_by(Taxon)%>%
        summarise(n=n())%>%
        arrange(-n)
    }
    aldex2.neg.effect<-aldex2.signif.features%>%
      filter(effect<0)
    head(aldex2.neg.effect)%>%
      knitr::kable(format = "simple")%>%
      print()
    
    #' Save the significant features with negative effect size as a 
    #'  tab-separated file.
    if(file.exists(aldex2.neg.effect.fname)){
      message(aldex2.neg.effect.fname, " exists. Skipping")
      
    }else{
      write.table(aldex2.neg.effect,
                  file = aldex2.neg.effect.fname,
                  row.names = F,sep = "\t")
    }
    #' Save all significant features.
    if(file.exists (aldex2.signif.features.fname)){
      message(aldex2.signif.features.fname, " exists. Skipping")
    }else{
      write.table(aldex2.signif.features,
                  file = aldex2.signif.features.fname,
                  row.names = F,sep = "\t")    
    }
  }
  
  #' Save the workspace just in case.
  final.workspace.fname<-file.path(diffabund.rdafiles,paste(
    "aldex2",output.filename,"workspace.RData",sep="-"))
  if(!file.exists(final.workspace.fname)){
    save.image(final.workspace.fname)
    
  }
  
  return(list(ps.q.glm.effect = ps.q.glm.effect,
              ps.q.aldex.clr = ps.q.aldex.clr,
              ps.q.glm.test = ps.q.glm.test,
              aldex2.neg.effect = aldex2.neg.effect,
              aldex2.signif.features = aldex2.signif.features))
}

#' Function that performs downstream processing of output from three tools 
#' (if you compare between hosts) or one tool (if you compare within host).
analyse_test_output <-function(agglom.rank,
                               comparison, 
                               comparison.col_name, # the column that corresponds to comparison
                               ref.level, 
                               output.filename,
                               mean_sd_relab,
                               maaslin2.signif.decreased,
                               sample.groups = NULL ,
                               aldex2.neg.effect = NULL,
                               ancombc.signif.decreased = NULL){
  
  # if( comparison =="host"){
  # function that filters a list to remove NULL vectors or empty vectors 
  drop_empty <- function(list_x ) list_x[sapply(list_x, function(x) !is.null(x) & !length(x) == 0)]
  # create a list
  if (!is.null(aldex2.neg.effect)){
    aldex2.neg.effect <- aldex2.neg.effect %>%
      group_by(Taxon)%>%
      add_count(Taxon)%>%
      filter(n ==length(sample.groups)-1)      
  }
  
  available.decreased <- list(maaslin2 = unique(maaslin2.signif.decreased$feature),
                              aldex2 = aldex2.neg.effect$Taxon,
                              ancombc = ancombc.signif.decreased$taxon_id)
  available.decreased <- drop_empty(available.decreased)
  
  #' Find common significantly decreased features between three tools.
  if(length(available.decreased) >0){
    # Find all intersections (venn diagram style, but in a vectorised form)
    venn.decreased <- venn(available.decreased, show.plot = FALSE)
    plot(venn.decreased)
    # Get all sets (both shared and unique to each tool)
    venn.decreased.all_combinations <- attr(venn.decreased, "intersections")
    # Merge the list into a dataframe 
    # The dataframe will have two columns (values, ind), where values are 
    # differentially abundant features, and ind are tool names
    venn.decreased.df <- stack(venn.decreased.all_combinations)
    names(venn.decreased.df) <- c(agglom.rank, "Tool")
    
    if(length(available.decreased) >=2){
      message("Tools that detected significantly decreased features: ", 
              paste(names(available.decreased), collapse = ", "))
      common.decreased<-Reduce(intersect, available.decreased)
      
      # If no common decreased features, unite them all:
      if(length(common.decreased) == 0){
        message("No shared features between all tools.")
        # Check pairwise combinations: look for colon (":") in the names of intersections
        # "intersections" can also have elements not shared between groups (single group)
        venn.decreased.pairs <- venn.decreased.all_combinations[grepl(":", names(venn.decreased.all_combinations))]
        message("Pairs of tools that have shared features: ", names(venn.decreased.pairs))
        print (venn.decreased.pairs)
        
        message("Combining all features into one vector:")
        common.decreased <- unlist(available.decreased, use.names = FALSE)
        common.decreased <- unique(common.decreased)
        print(common.decreased)
      }else{
        message("Features that are shared between all tools: ")
        print(common.decreased)
      }
    }else if(length(available.decreased) == 1){
      message("Only one tool detected significantly decreased features: ", names(available.decreased))
      common.decreased<- venn.decreased.df[[agglom.rank]]
      print(common.decreased)
    }
    
  }else if (length(available.decreased) ==0){
    message("No differentially abundant features!")
    return(NULL)
  }
  
  # If we identified differentially abundant features, save them
  if(nrow(venn.decreased.df) >=1){
    # Merge the list into a dataframe and filter features
    # The dataframe will have two columns (values, ind), where values are 
    # differentially abundant features, and ind are tool names
    # common.decreased.df <- stack(available.decreased)%>%
    #   filter(values %in% common.decreased )
    # names(common.decreased.df) <- c(agglom.rank, "Tool")
    # tool names for saving the df
    venn.decreased.df.tools <- unique (venn.decreased.df$Tool)
    
    venn.decreased.df<- mean_sd_relab%>% 
      filter(get(comparison.col_name) == ref.level, 
             get(agglom.rank)%in%venn.decreased.df[[agglom.rank]])%>%
      right_join(venn.decreased.df, by = agglom.rank)%>%
      relocate(Tool, .after = agglom.rank)%>%
      arrange(Tool, -MeanRelativeAbundance)%>%
      distinct()
    # Save the common significantly decreased features with their relative abundance
    message("Showing a few differentially abundant taxa")
    head(venn.decreased.df)%>%
      knitr::kable(format = "simple")%>%
      print()
    venn.decreased.df.across_tools.fname <- file.path(diffabund.tables,
                                                      paste("common.decreased",
                                                            paste(names(available.decreased),collapse = "-"),
                                                            output.filename,
                                                            "signif.tsv",sep="-"))
    if(file.exists(venn.decreased.df.across_tools.fname)){
      message(venn.decreased.df.across_tools.fname, " exists. Skipping")
    }else{
      write.table(venn.decreased.df,
                  file = venn.decreased.df.across_tools.fname,
                  row.names = F,sep = "\t")
      
    }
    return( list(venn.decreased = venn.decreased, # vectorised form that can be plotted
                 venn.decreased.df = venn.decreased.df)) # dataframe form with taxonomy
    
  }
  # }
}

#' Plot differentially abundant features.
plot_differentially_abundant_taxa<- function(agglom.rank,
                                             comparison,
                                             ref.level,
                                             output.filename,
                                             metadata,
                                             
                                             custom.levels,
                                             ps.rel.agg,
                                             ggplot.levels,
                                             ggplot.labels,
                                             gg.labs.name,
                                             gg.title.groups,
                                             features.to_plot,
                                             maaslin2.signif.features = NULL
){
  if(setequal(custom.levels,"NMR")){
    pretty.asv.names.df<-ps.rel.agg%>%
      ungroup()%>%
      filter(get(agglom.rank)%in%features.to_plot,class=="NMR")%>%
      distinct(get(agglom.rank),.keep_all = T)%>%
      dplyr::select(OTU,Genus)%>%
      left_join(maaslin2.signif.features, by = join_by(OTU == "feature") )%>%
      mutate(Taxon=paste0("ASV from ","<i>",Genus,"</i> (p = ",
                          round(qval,digits=4),")"))
    pretty.asv.names<-c(pretty.asv.names.df$Taxon)
    names(pretty.asv.names)<-pretty.asv.names.df$OTU
    
    if(comparison =="age"){
      feature.plot<-ps.rel.agg%>%
        filter(get(agglom.rank)%in%features.to_plot,class=="NMR")%>%
        left_join(metadata, by = "Sample")%>%
        group_by_at(c(agglom.rank,"agegroup"))%>%
        ggplot(aes(x=factor(agegroup,levels=names(ggplot.labels)),
                   y=Abundance,
                   fill=factor(agegroup)))
    }else if(comparison =="sex"){
      feature.plot<-ps.rel.agg%>%
        filter(get(agglom.rank)%in%features.to_plot,class=="NMR")%>%
        group_by_at(c(agglom.rank,"sex"))%>%
        ggplot(aes(x=factor(sex,levels=names(ggplot.labels)),
                   y=Abundance,
                   fill=factor(sex)))
    }
    feature.plot<-feature.plot +
      geom_boxplot(show.legend = FALSE)+
      facet_wrap(~OTU,scales = "free_y",
                 ncol = 2,
                 labeller = as_labeller(pretty.asv.names))+
      theme_bw()+
      coord_flip()+
      labs(x="",
           y="Relative abundance (%)")+
      scale_x_discrete(labels=ggplot.labels,
                       limits=ggplot.levels)+ # rename boxplot labels (x axis)
      # scale_fill_manual(values = custom.fill)+
      scale_fill_viridis_d(direction = (-1),
                           end = 0.9,
                           alpha = 0.5)+
      # ggtitle(paste0("Relative abundance of differentially abundant ASVs \nin different naked mole-rat groups"))+
      diffabund.plot.theme
    
    if(length(features.to_plot < 10 )){ # increase the size if we have few features
      feature.plot <- feature.plot +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = ggtext::element_markdown(size = 12),
              axis.title.x = element_text(size = 12),
              strip.text.x = ggtext::element_markdown(size = 12))
    }
    
    for(image.format in image.formats){
      ggsave(paste0(paste(ref.level,"specific-bacteria",output.filename,
                          sep = "-"),".",image.format),
             plot = feature.plot,
             path = diffabund.figures,
             width = 11, height = 8,units="in",
             dpi = 300,device = image.format)
    }
  }else{
    feature.plot<-ps.rel.agg%>%
      filter(get(agglom.rank)%in%features.to_plot)%>%
      group_by_at(c("class",agglom.rank))%>%
      ggplot(aes(x=factor(class,levels=rev(ggplot.levels)),
                 y=Abundance,
                 fill=factor(class,levels=rev(ggplot.levels))))+
      geom_boxplot(show.legend = FALSE)+
      facet_wrap(~Genus,scales = "free_x",
                 ncol = 4)+
      theme_bw(base_size = 5)+
      coord_flip()+
      labs(x="",
           y="Relative abundance (%)")+
      scale_color_viridis_d(breaks = rev(unname(ggplot.labels)),
                            labels=rev(unname(ggplot.labels)))+
      scale_x_discrete(labels=rev(ggplot.labels),
                       limits=rev(ggplot.levels))+ # rename boxplot labels (x axis)
      scale_fill_viridis_d(option="C")+
      # ggtitle(paste0("Relative abundance of naked mole-rat-specific taxa"))+
      diffabund.plot.theme
    
    
    for(image.format in image.formats){
      ggsave(paste0(paste(output.filename,"specific-bacteria",
                          sep = "-"),".",image.format),
             plot=feature.plot,
             path = diffabund.figures,
             width=8, height=11,units="in",
             dpi=300,device = image.format)
    }
    
  }
  print(feature.plot)
  return(feature.plot)
}



#+ echo=FALSE
## 4. Compare hosts (Genus level). ####
#' 
#' ## Compare hosts (Genus level).
#' Choose what to compare:
comparison<-"host"
#' Choose the reference level:
ref.level<-"NMR" 
agglom.rank<-"Genus"

#' Prepare data for tests.
custom.md <- readRDS(custom.md.path)
mean_sd_relab.all_hosts<- read.table(file.path(community.composition.tables,
                                               paste0("mean-max-min-sd-relab-all_hosts-",agglom.rank,".tsv")),
                                     sep = "\t", header = T)

host.data.for_test <- prepare_data(ps = ps.q,
                                   agglom.rank = agglom.rank, 
                                   comparison = comparison, 
                                   ref.level = ref.level, 
                                   metadata = custom.md)

#+ echo=FALSE
### 4.1 Run a test with MaAsLin2. ####
#'
#' ### Run a test with MaAsLin2.
maaslin.fit_data.host <- perform_maaslin2_test(ps.agg.wide = host.data.for_test$ps.agg.wide, 
                                               output.filename = host.data.for_test$output.filename,
                                               comparison = comparison,
                                               metadata = host.data.for_test$metadata,
                                               ref.level = ref.level,
                                               maaslin2.random_effects.file_suffix = host.data.for_test$maaslin2.random_effects.file_suffix,
                                               maaslin2.fixed_effects = host.data.for_test$maaslin2.fixed_effects,
                                               maaslin2.reference = host.data.for_test$maaslin2.reference
)
#' Process the output 
maaslin.processed_output.host <- 
  process_maaslin2_output(maaslin.fit_data = maaslin.fit_data.host,
                          ps.agg = host.data.for_test$ps.agg,
                          sample.groups = host.data.for_test$sample.groups,
                          agglom.rank = agglom.rank,
                          output.filename = host.data.for_test$output.filename)

#+ echo=FALSE
### 4.2 Run a test with ALDEx2. ####
#'
#' ### Run a test with ALDEx2.
aldex2.host <- perform_aldex2_test(comparison = comparison,
                                   ps.agg.wide = host.data.for_test$ps.agg.wide,
                                   metadata = host.data.for_test$metadata,
                                   mm = host.data.for_test$mm,
                                   sample.groups = host.data.for_test$sample.groups,
                                   aldex.reference = host.data.for_test$aldex.reference,
                                   output.filename = host.data.for_test$output.filename)


#+ echo=FALSE
### 4.3 Run a test with ANCOMBC. ####
#'
#' ### Run a test with ANCOMBC.
ancombc.host <-perform_ancombc_test(ps.agg.wide = host.data.for_test$ps.agg.wide,
                                    metadata = host.data.for_test$metadata,
                                    ancombc.comparison = host.data.for_test$ancombc.comparison,
                                    ancombc.levels = host.data.for_test$ancombc.levels,
                                    output.filename = host.data.for_test$output.filename)
#+ echo=FALSE
### 4.4 Downstream analysis and plotting of differentially abundant features. ####
#' 
#' ### Downstream analysis and plotting of differentially abundant features. 
#+ fig.width=11, fig.height=8
diffabund.test.summary.host <- analyse_test_output(agglom.rank = agglom.rank,
                                                   comparison = comparison,
                                                   comparison.col_name = host.data.for_test$maaslin2.fixed_effects,
                                                   ref.level = ref.level,
                                                   output.filename = host.data.for_test$output.filename,
                                                   mean_sd_relab = mean_sd_relab.all_hosts,
                                                   sample.groups = host.data.for_test$sample.groups,
                                                   maaslin2.signif.decreased = maaslin.processed_output.host$maaslin2.signif.decreased,
                                                   aldex2.neg.effect = aldex2.host$aldex2.neg.effect,
                                                   ancombc.signif.decreased = ancombc.host$ancombc.signif.decreased)
plot(diffabund.test.summary.host$venn.decreased)

plot_taxa.host <- plot_differentially_abundant_taxa(agglom.rank = agglom.rank,
                                                    comparison = comparison,
                                                    ref.level = ref.level,
                                                    output.filename = host.data.for_test$output.filename,
                                                    metadata = host.data.for_test$metadata,
                                                    custom.levels = host.data.for_test$custom.levels,
                                                    ps.rel.agg = host.data.for_test$ps.rel.agg,
                                                    ggplot.levels = host.data.for_test$ggplot.levels,
                                                    ggplot.labels = host.data.for_test$ggplot.labels,
                                                    gg.labs.name = host.data.for_test$gg.labs.name,
                                                    gg.title.groups = host.data.for_test$gg.title.groups,
                                                    features.to_plot = unique(diffabund.test.summary.host$venn.decreased.df[[agglom.rank]])
)


#+ echo=FALSE
## 5. Compare NMR age groups (ASV level). ####
#' 
#' ## Compare NMR age groups (ASV level).
#' Choose what to compare:
comparison<-"age"
#' Choose the reference level:
ref.level<-"agegroup0_10" 
agglom.rank<-"OTU"
ps.q.nmr <-readRDS(ps.q.filtered.fname) %>%
  subset_samples(class == "NMR")
nmr.relations<-read.table(file.path(new.metadata.dir,"nmr-relations.tsv"),
                          header = T,
                          sep = "\t")
#' Prepare the data for testing.
custom.md<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")

mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-asv.tsv")
mean_sd_relab.nmr_ages.asv <- read.table(mean_sd_relab.nmr_ages.asv.fname,
                                         header = TRUE,
                                         sep = "\t")

nmr.age.data.for_test<-prepare_data(ps = ps.q.nmr,
                                    agglom.rank = agglom.rank,
                                    comparison = comparison,
                                    ref.level = ref.level,
                                    metadata = custom.md,
                                    relations = nmr.relations)

#' Run MaAsLin2 test.
maaslin.fit_data.age<-
  perform_maaslin2_test(ps.agg.wide = nmr.age.data.for_test$ps.agg.wide,
                        output.filename = nmr.age.data.for_test$output.filename,
                        comparison = comparison,
                        metadata = nmr.age.data.for_test$metadata,
                        ref.level = ref.level,
                        maaslin2.random_effects = nmr.age.data.for_test$maaslin2.random_effects,
                        maaslin2.random_effects.file_suffix =nmr.age.data.for_test$maaslin2.random_effects.file_suffix,
                        maaslin2.fixed_effects = nmr.age.data.for_test$maaslin2.fixed_effects,
                        maaslin2.reference = nmr.age.data.for_test$maaslin2.reference)

maaslin.processed_output.age <- 
  process_maaslin2_output(maaslin.fit_data = maaslin.fit_data.age,
                          ps.agg = nmr.age.data.for_test$ps.agg,
                          sample.groups = nmr.age.data.for_test$sample.groups,
                          agglom.rank = agglom.rank,
                          output.filename = nmr.age.data.for_test$output.filename)

#' ### Downstream analysis and plotting of differentially abundant features. 
#+ fig.width=11, fig.height=8
diffabund.test.summary.age <- analyse_test_output(agglom.rank = agglom.rank,
                                                  comparison = comparison,
                                                  comparison.col_name = nmr.age.data.for_test$maaslin2.fixed_effects,
                                                  ref.level = ref.level,
                                                  output.filename = nmr.age.data.for_test$output.filename,
                                                  mean_sd_relab = mean_sd_relab.nmr_ages.asv,
                                                  maaslin2.signif.decreased = maaslin.processed_output.age$maaslin2.signif.decreased)


plot_taxa.age <- 
  plot_differentially_abundant_taxa(agglom.rank = agglom.rank,
                                    comparison = comparison,
                                    ref.level = ref.level,
                                    output.filename = nmr.age.data.for_test$output.filename,
                                    metadata = nmr.age.data.for_test$metadata,
                                    custom.levels = nmr.age.data.for_test$custom.levels,
                                    ps.rel.agg = nmr.age.data.for_test$ps.rel.agg,
                                    ggplot.levels = nmr.age.data.for_test$ggplot.levels,
                                    ggplot.labels = nmr.age.data.for_test$ggplot.labels,
                                    gg.labs.name = nmr.age.data.for_test$gg.labs.name,
                                    gg.title.groups = nmr.age.data.for_test$gg.title.groups,
                                    features.to_plot = unique(diffabund.test.summary.age$venn.decreased.df[[agglom.rank]]),
                                    maaslin2.signif.features = maaslin.processed_output.age$maaslin2.signif.features
  )

#+ echo=FALSE
## 6. Compare NMR sexes (ASV level). ####
#' 
#' ## Compare NMR sexes (ASV level).
#' Choose what to compare:
comparison<-"sex"
agglom.rank<-"OTU"
#' Choose the reference level:
ref.level<-"female"
#' Prepare the data for testing.
custom.md<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
nmr.sex.data.for_test<-prepare_data(ps = ps.q.nmr,
                                    agglom.rank = agglom.rank,
                                    comparison = comparison,
                                    ref.level = ref.level,
                                    metadata = custom.md,
                                    relations = nmr.relations)
mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-asv.tsv")
mean_sd_relab.nmr_ages.asv <- read.table(mean_sd_relab.nmr_ages.asv.fname,
                                         header = TRUE,
                                         sep = "\t")
#' Run MaAsLin2 test.
maaslin.fit_data.sex<-
  perform_maaslin2_test(ps.agg.wide = nmr.sex.data.for_test$ps.agg.wide,
                        output.filename = nmr.sex.data.for_test$output.filename,
                        comparison = comparison,
                        metadata = nmr.sex.data.for_test$metadata,
                        ref.level = ref.level,
                        maaslin2.random_effects = nmr.sex.data.for_test$maaslin2.random_effects,
                        maaslin2.random_effects.file_suffix = nmr.sex.data.for_test$maaslin2.random_effects.file_suffix,
                        maaslin2.fixed_effects = nmr.sex.data.for_test$maaslin2.fixed_effects,
                        maaslin2.reference = nmr.sex.data.for_test$maaslin2.reference
  )

maaslin.processed_output.sex <- 
  process_maaslin2_output(maaslin.fit_data = maaslin.fit_data.sex,
                          ps.agg = nmr.sex.data.for_test$ps.agg,
                          sample.groups = nmr.sex.data.for_test$sample.groups,
                          agglom.rank = agglom.rank,
                          output.filename = nmr.sex.data.for_test$output.filename)

ps.q.nmr.rel <- readRDS(ps.q.rel.filtered.fname) %>%
  subset_samples(class == "NMR")
source(file.path(util.functions.r,"get_mean_and_sd_relab_from_phyloseq.R"))
mean_sd_relab.nmr_sexes.asv <- get_mean_and_sd_relab_from_phyloseq(ps.q.nmr.rel,
                                                                   group.name = "sex",
                                                                   tax.rank = "OTU")
mean_sd_relab.nmr_sexes.asv <- mean_sd_relab.nmr_sexes.asv %>%
  mutate(sex = case_when(sex == "F" ~ "female",
                         sex == "M" ~ "male"))
nmr.sex.data.for_test$ps.rel.agg <- nmr.sex.data.for_test$ps.rel.agg%>%
  mutate(sex = case_when(sex == "F" ~ "female",
                         sex == "M" ~ "male"))

#' ### Downstream analysis and plotting of differentially abundant features. 
#+ fig.width=11, fig.height=8
diffabund.test.summary.sex <- 
  analyse_test_output(agglom.rank = agglom.rank,
                      comparison = comparison,
                      comparison.col_name = nmr.sex.data.for_test$maaslin2.fixed_effects,
                      ref.level = ref.level,
                      output.filename = nmr.sex.data.for_test$output.filename,
                      mean_sd_relab = mean_sd_relab.nmr_sexes.asv,
                      maaslin2.signif.decreased = maaslin.processed_output.sex$maaslin2.signif.decreased)
plot_taxa.sex <- 
  plot_differentially_abundant_taxa(agglom.rank = agglom.rank,
                                    comparison = comparison,
                                    ref.level = ref.level,
                                    output.filename = nmr.sex.data.for_test$output.filename,
                                    metadata = nmr.sex.data.for_test$metadata,
                                    custom.levels = nmr.sex.data.for_test$custom.levels,
                                    ps.rel.agg = nmr.sex.data.for_test$ps.rel.agg,
                                    ggplot.levels = nmr.sex.data.for_test$ggplot.levels,
                                    ggplot.labels = nmr.sex.data.for_test$ggplot.labels,
                                    gg.labs.name = nmr.sex.data.for_test$gg.labs.name,
                                    gg.title.groups = nmr.sex.data.for_test$gg.title.groups,
                                    features.to_plot = unique(diffabund.test.summary.sex$venn.decreased.df[[agglom.rank]]),
                                    maaslin2.signif.features = maaslin.processed_output.sex$maaslin2.signif.features
  )



sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()
