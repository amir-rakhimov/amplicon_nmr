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
library(Polychrome)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
rare.status<-"rare"
filter.status<-"nonfiltered"
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
print(paste("The analysis focus is:", active.analysis))

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

dir.create(diffabund.tables,recursive = TRUE)
dir.create(diffabund.rdafiles,recursive = TRUE)
dir.create(diffabund.figures,recursive = TRUE)

#+ echo=FALSE
## 3. Prepare necessary functions. ####
#'
#' ## Prepare necessary functions. ####
#' Function that performs differential abundance tests with different 
#' parameters. Instead of repeating the code for each comparison,
#' call the function and change the parameters.
prepare_data <- function(agglom.rank, comparison, 
                         ref.level, inside.host){
  print(paste("Performing differential microbial abundance tests on", comparison, "variable at", agglom.rank, "level.",
              "Reference:", ref.level))
  #' Import abundance table as an rds file (NOT rarefied):
  if(inside.host == TRUE){
    custom.levels<-"NMR"
    custom.md<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
      filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
  }else if(inside.host ==FALSE){
    if(agglom.rank == "Genus" ){
      custom.levels<-c("NMR",
                       "DMR",
                       "B6mouse",
                       "MSMmouse",
                       "FVBNmouse",
                       "spalax",
                       "pvo",
                       "hare",
                       "rabbit")
      custom.md <- readRDS(custom.md.path)
    } 
  }
  #' Import abundance table as an rds file (NOT rarefied):
  ps.q.agg<-readRDS(file=file.path(
    community.composition.rdafiles,
    paste("phyloseq-qiime",authorname,agglom.rank,read.end.type,
          truncationlvl,"table.rds",sep = "-")))
  ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
  
  #' Rarefied abundance table:
  ps.q.df.preprocessed<-readRDS(
    file.path(community.composition.rdafiles,
              paste0(
                paste(paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
                      paste(custom.levels,collapse = '-'),sep = "-"),
                ".rds")))%>%
    filter(class%in%custom.levels, Abundance!=0)%>%
    dplyr::select(all_of(c("Sample",agglom.rank,"Abundance","class")))
  
  # if(inside.host ==TRUE){
  #   ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
  #     filter(class=="NMR",Abundance!=0)
  # }
  
  #' Specify output file name:
  output.filename<-paste(paste(custom.levels,collapse = '-'),
                         agglom.rank,comparison,
                         "ref",ref.level,
                         sep="-")
  if(inside.host == TRUE){
    short.filename<-paste("NMR",agglom.rank,comparison,
                           "ref",ref.level,
                           sep="-")
  }else if(inside.host ==FALSE){
    short.filename<-paste("all_hosts",agglom.rank,comparison,
                           "ref",ref.level,
                           sep="-")
  }
  print(paste("Output filename prefix:", output.filename))
  # Preparing the wide dataset.
  ps.q.df.wide<-ps.q.df.preprocessed%>%
    
    # filter(Abundance!=0)%>%
    pivot_wider(names_from = agglom.rank, # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("Sample")%>%
    dplyr::select(-class)
  # Colnames are OTUs and rownames are sample IDs.
  if(setequal(custom.levels,"NMR")){
    if(comparison=="age"){
      sample.groups<-c("agegroup0_10",
                       "agegroup10_16")
    }else if(comparison=="sex"){
      sample.groups<-c("female",
                       "male")
    }
  }else{
    sample.groups<-custom.levels
  }
  print(sample.groups)
  print("Prepared the data for tests")
  
  return(list(dataset = ps.q.df.wide, 
              ps.q.agg = ps.q.agg,
              output.filename = output.filename, 
              short.filename = short.filename,
              metadata = custom.md, 
              custom.levels = custom.levels,
              sample.groups = sample.groups))
}

#' Function that performs a test with MaAsLin2:
perform_maaslin2_test <- function(ps.q.df.wide, 
                                 output.filename,
                                 comparison,
                                 custom.md,
                                 custom.levels,
                                 ref.level,
                                 inside.host){
  # Create reference levels.
  if (comparison=="age"){
    maaslin.reference<-paste("agegroup", ref.level, sep = ",")
    maaslin.comparison<-"agegroup"
  }else if(comparison=="sex"){
    maaslin.reference<-paste("sex", ref.level, sep = ",")
    maaslin.comparison<-"sex"
  }
  
  if (inside.host == TRUE){
    if(setequal(custom.levels,"NMR")){
      relations<-read.table("./data/metadata/pooled-metadata/nmr-relations.tsv",
                            header = T,
                            sep = "\t")
      custom.md<-custom.md%>%
        left_join(relations,by="Sample")
      rownames(custom.md)<-custom.md$Sample
    }
  }
  
  # Run Maaslin2.
  if(inside.host==TRUE){
    set.seed(1)
    maaslin.fit_data = 
      Maaslin2(input_data = ps.q.df.wide, 
               input_metadata = custom.md, 
               min_prevalence = 0,
               normalization = "TSS",
               transform = "LOG",
               analysis_method = "LM",
               random_effects = c("relation"), 
               standardize = FALSE,
               output = file.path(maaslin2.out,
                                  paste0(authorname,"-output"),
                                  rare.status,paste(
                                    output.filename, "with-relations",sep = "-")), 
               fixed_effects = maaslin.comparison,
               reference = maaslin.reference,
               max_significance = 0.05)
  }else{
    set.seed(1)
    maaslin.fit_data = 
      Maaslin2(input_data = ps.q.df.wide, 
               input_metadata = custom.md, 
               min_prevalence = 0,
               analysis_method = "LM",
               normalization = "TSS",
               transform = "LOG",
               random_effects = NULL,
               standardize = FALSE,
               output = file.path(maaslin2.out,
                                  paste0(authorname,"-output"),
                                  rare.status,paste(
                                    output.filename,sep = "-")), 
               fixed_effects = c("class"),
               reference = paste0("class,",ref.level),
               max_significance = 0.05)
    
  }
  return(maaslin.fit_data)
}

#' Function to process the MaAsLin2 output:
process_maaslin2_output<-function(agglom.rank,
                                  short.filename,
                                  maaslin.fit_data,
                                  ps.q.agg,
                                  sample.groups){
  #' Downstream processing of Maaslin2 output.
  source(file.path(util.functions.r,"make_features_maaslin.R"))
  #' Extract features with qvalue < 0.05.
  if(min(maaslin.fit_data$results$qval)<0.05){
    maaslin.signif.features<-maaslin.fit_data$results%>%
      filter(qval<0.05) # should be qval
  }else{
    maaslin.signif.features<-maaslin.fit_data$results%>%
      arrange(qval)%>%
      head(n = 15) # if no significant results found
  }
  
  #' Make features pretty. We have to make this exchange because maaslin 
  #' output treats space and hyphen as the same thing.
  if(agglom.rank=="OTU"){
    maaslin.signif.features$feature<-gsub("^X","",
                                          maaslin.signif.features$feature)
  }else{
    foo<-ps.q.agg
    foo<-foo%>%mutate("maaslin"=get(agglom.rank))
    foo<-make_features_maaslin(foo,"maaslin")
    foo<-unique(foo[,c("maaslin",agglom.rank)])
    maaslin.signif.features<-maaslin.signif.features%>%
      left_join(foo[,c("maaslin",agglom.rank)],by=c("feature"="maaslin"))%>%
      distinct()
    rm(foo)
    maaslin.signif.features$feature<-maaslin.signif.features[,agglom.rank]
    maaslin.signif.features<-subset(maaslin.signif.features, 
                                    select=-get(agglom.rank))
  }
  
  #' Extract features that are downregulated in all hosts simultaneously.
  # if(setequal(custom.levels,"NMR")){
  #   if(comparison=="age"){
  #     sample.groups<-c("agegroup0_10",
  #                      "agegroup10_16")
  #   }else if(comparison=="sex"){
  #     sample.groups<-c("F",
  #                      "M")
  #   }
  # }else{
  #   sample.groups<-custom.levels
  # }
  # Downreglated features have coef<0.
  # We also add association strength column
  maaslin.signif.decreased<-maaslin.signif.features%>%
    as_tibble()%>%
    filter(coef<0)%>%
    group_by(feature)%>%
    filter(n()==length(sample.groups)-1)%>% 
    arrange(pval,feature)%>%
    mutate(n=n())%>%
    mutate(assoc.str=-log(qval)*sign(coef))%>%
    ungroup()
  
  
  #' Check if all features (OTU/taxa/etc) are also present in the 
  #' ps.q.agg dataframe.
  if(agglom.rank=="OTU"){
    table(maaslin.signif.decreased$feature%in%ps.q.agg$OTU)
  }else{
    table(maaslin.signif.decreased$feature%in%ps.q.agg[,agglom.rank])
  }
  write.table(maaslin.signif.features,
              file=file.path(diffabund.tables,paste(
                "maaslin2",short.filename,
                "signif.tsv",sep="-")),
              row.names = F,sep = "\t")
  write.table(maaslin.signif.decreased,
              file=file.path(diffabund.tables,
                             paste("maaslin.signif.decreased",
                                   short.filename,
                                   "signif.tsv",sep="-")),
              row.names = F,sep = "\t")
  return(list(maaslin.signif.features = maaslin.signif.features,
         maaslin.signif.decreased = maaslin.signif.decreased))
  
}

#' Function that performs a test with ANCOM-BC:
perform_ancombc_test<-function(agglom.rank,
                               ps.q.df.wide,
                               ps.q.agg,
                               custom.md, 
                               sample.groups){
  #' Extract the taxonomic table.
  if(agglom.rank=="OTU"){
    taxmat<-ps.q.agg%>%
      dplyr::select(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>%
      distinct()%>%
      column_to_rownames(var = "OTU")%>%
      as.matrix()
  }else{
    all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
    agglom.rank.index<-match(agglom.rank,all.ranks)
    custom.ranks<-all.ranks[1:agglom.rank.index]
    
    taxmat<-ps.q.agg%>%
      ungroup()%>%
      dplyr::select(all_of(custom.ranks))%>%
      distinct()%>%
      column_to_rownames(var = all_of(agglom.rank))%>%
      as.matrix()
  }
  #' Create a phyloseq object for ANCOMBC.
  ps.q.OTU<-t(ps.q.df.wide)
  ps.q.OTU<-otu_table(ps.q.OTU,taxa_are_rows = T)
  ps.q.TAX<-tax_table(taxmat)
  ps.q.phyloseq.new<-phyloseq(otu_table(ps.q.OTU),
                              tax_table(ps.q.TAX),
                              sample_data(custom.md))
  
  #' Relevel the comparison vector. The first level will be the reference 
  #' for custom leveling
  ancombc.levels<-c(ref.level,
                    sample.groups[sample.groups!=ref.level])
  if(comparison=="host"){
    sample_data(ps.q.phyloseq.new)$class<-
      factor(sample_data(ps.q.phyloseq.new)$class,
             levels = ancombc.levels)
    ancombc.comparison<-"class"                                           
  }else if (comparison=="age"){
    sample_data(ps.q.phyloseq.new)$agegroup<-
      factor(sample_data(ps.q.phyloseq.new)$agegroup,
             levels = ancombc.levels)
    ancombc.comparison<-"agegroup"
  }else if(comparison=="sex"){
    sample_data(ps.q.phyloseq.new)$sex<-
      factor(sample_data(ps.q.phyloseq.new)$sex,
             levels = ancombc.levels)
    ancombc.comparison<-"sex"
  }
  
  #' Perform the differential abundance test.
  ancombc.out<-ancombc(
    data = ps.q.phyloseq.new,
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
  return(ancombc.out)
}


#' Function that performs downstream processing of output from three tools 
#' (if you compare between hosts) or one tool (if you compare within host).
analyse_test_output <-function(agglom.rank,
                               inside.host, 
                               comparison, 
                               custom.levels,
                               ref.level, 
                               sample.groups,
                               output.filename,
                               custom.md,
                               ps.q.agg,
                               maaslin.signif.features,
                               maaslin.signif.decreased,
                               aldex.signif.features = NULL,
                               aldex.neg.effect = NULL,
                               ancombc.signif.features = NULL,
                               ancombc.signif.decreased = NULL){

  if(inside.host==FALSE){
    #' Find common significant features between three tools.
    print(Reduce(intersect,list(maaslin.signif.features$feature,
                                aldex.signif.features$Taxon,
                                ancombc.signif.features$taxon_id)))
    
    #' Find common significantly decreased features between three tools.
    print(Reduce(intersect,list(maaslin.signif.decreased$feature,
                                aldex.neg.effect$Taxon,
                                ancombc.signif.decreased$taxon_id)))
    
    #' Common significant and decreased features between tools
    common.signif<-Reduce(intersect,list(maaslin.signif.features$feature,
                                         aldex.signif.features$Taxon,
                                         ancombc.signif.features$taxon_id))
    #' Only Maaslin2 and ANCOM-BC: they're decreased in other hosts, not in ref.
    common.decreased<-Reduce(intersect,list(maaslin.signif.decreased$feature,
                                            ancombc.signif.decreased$taxon_id))
    
    common.decreased.df<-ps.q.agg%>%
      filter(class==ref.level,Genus%in%common.decreased)%>%
      distinct(Genus,.keep_all = T)%>%
      ungroup()%>%
      arrange(-MeanRelativeAbundance)%>%
      dplyr::select(Genus,MeanRelativeAbundance)
    head(common.decreased.df)
    write.table(common.decreased.df,
                file=file.path(diffabund.tables,
                               paste("significant-features",
                                     output.filename,
                                     "signif.tsv",sep="-")),
                row.names = F,sep = "\t")
  }
  
  
  #' Plot differentially abundant features.
  if(comparison=="host"){
    pretty.level.names<-
      c("NMR" = "*Heterocephalus glaber*", # better labels for facets
        "DMR" = "*Fukomys damarensis*",
        "B6mouse" = "B6 mouse",
        "MSMmouse" = "MSM/Ms mouse",
        "FVBNmouse" = "FVB/N mouse",
        "spalax" = "*Nannospalax leucodon*",
        "pvo" = "*Pteromys volans orii*",
        "hare" = "*Lepus europaeus*",
        "rabbit" = "*Oryctolagus cuniculus*"
      )
    ggplot.levels<-names(pretty.level.names)
    gg.labs.name<-"Host"
    gg.title.groups<-"hosts"
  }else if (comparison=="age"){
    pretty.level.names<-names(table(custom.md$agegroup))
    
    names(pretty.level.names)<-custom.md%>%
      ungroup%>%
      dplyr::select(agegroup)%>%
      distinct(agegroup)%>%
      arrange(agegroup)%>%
      pull
    
    ggplot.levels<-names(pretty.level.names)
    gg.labs.name<-"Age group"
    gg.title.groups<-"age groups"
    
  }else if (comparison=="sex"){
    pretty.level.names<-
      c("female" = "Females",
        "male" = "Males")
    ggplot.levels<-names(pretty.level.names)
    pretty.level.names<-pretty.level.names[
      which(names(pretty.level.names)%in%ggplot.levels)]
    gg.labs.name<-"Host sex"
    gg.title.groups<-"groups"
  }
  
  if(inside.host==TRUE){
    pretty.asv.names.df<-ps.q.agg%>%
      ungroup()%>%
      filter(get(agglom.rank)%in%maaslin.signif.features$feature,class=="NMR")%>%
      distinct(get(agglom.rank),.keep_all = T)%>%
      dplyr::select(OTU,Genus)%>%
      mutate(Taxon=paste0("ASV from ","<i>",Genus,"</i> (p = ",
                          round(maaslin.signif.features$qval,digits=3),")"))
    pretty.asv.names<-c(pretty.asv.names.df$Taxon)
    names(pretty.asv.names)<-pretty.asv.names.df$OTU
    
    if(comparison =="age"){
      feature.plot<-ps.q.agg%>%
        filter(get(agglom.rank)%in%maaslin.signif.features$feature,class=="NMR")%>%
        left_join(custom.md)%>%
        group_by_at(c(agglom.rank,"agegroup"))%>%
        ggplot(aes(x=factor(agegroup,levels=names(pretty.level.names)),
                   y=RelativeAbundance,
                   fill=factor(agegroup)))
    }else if(comparison =="sex"){
      feature.plot<-ps.q.agg%>%
        filter(get(agglom.rank)%in%maaslin.signif.features$feature,class=="NMR")%>%
        left_join(custom.md)%>%
        group_by_at(c(agglom.rank,"sex"))%>%
        ggplot(aes(x=factor(sex,levels=names(pretty.level.names)),
                   y=RelativeAbundance,
                   fill=factor(sex)))
    }
    feature.plot<-feature.plot +
      geom_boxplot(show.legend = FALSE)+
      facet_wrap(~OTU,scales = "free_y",
                 ncol = 2,
                 labeller = as_labeller(pretty.asv.names))+
      theme_bw()+
      labs(x="",
           y="Relative abundance (%)")+
      scale_x_discrete(labels=pretty.level.names,
                       limits=ggplot.levels)+ # rename boxplot labels (x axis)
      # scale_fill_manual(values = custom.fill)+
      scale_fill_viridis_d(direction = (-1),
                           end = 0.9,
                           alpha = 0.5)+
      # ggtitle(paste0("Relative abundance of differentially abundant ASVs \nin different naked mole-rat groups"))+
      diffabund.plot.theme
    
    
    for(image.format in image.formats){
      ggsave(paste0(paste(ref.level,"specific-bacteria","NMR",comparison,
                          sep = "-"),".",image.format),
             plot = feature.plot,
             path = diffabund.figures,
             width = 11, height = 8,units="in",
             dpi = 300,device = image.format)
    }
  }else{
    feature.plot<-ps.q.agg%>%
      filter(get(agglom.rank)%in%common.decreased)%>%
      group_by_at(c("class",agglom.rank))%>%
      ggplot(aes(x=factor(class,levels=rev(ggplot.levels)),
                 y=RelativeAbundance,
                 fill=factor(class,levels=rev(ggplot.levels))))+
      geom_boxplot(show.legend = FALSE)+
      facet_wrap(~Genus,scales = "free_x",
                 ncol = 4)+
      theme_bw(base_size = 5)+
      coord_flip()+
      labs(x="",
           y="Relative abundance (%)")+
      scale_color_viridis_d(breaks = rev(unname(pretty.level.names)),
                            labels=rev(unname(pretty.level.names)))+
      scale_x_discrete(labels=rev(pretty.level.names),
                       limits=rev(ggplot.levels))+ # rename boxplot labels (x axis)
      scale_fill_viridis_d(option="C")+
      # ggtitle(paste0("Relative abundance of naked mole-rat-specific taxa"))+
      diffabund.plot.theme
    
    
    for(image.format in image.formats){
      ggsave(paste0(paste("NMR-specific-bacteria",
                          sep = "-"),".",image.format),
             plot=feature.plot,
             path = diffabund.figures,
             width=8, height=11,units="in",
             dpi=300,device = image.format)
    }
    
  }
  print(feature.plot)
  
  
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
inside.host<- FALSE
#' Prepare data for tests.
host.data.for_test<-prepare_data(agglom.rank = agglom.rank, 
             comparison = comparison, 
             ref.level = ref.level, 
             inside.host = inside.host)
#+ echo=FALSE
### 4.1 Run a test with MaAsLin2. ####
#'
#' ### Run a test with MaAsLin2.
maaslin.fit_data.host<-
  perform_maaslin2_test(ps.q.df.wide = host.data.for_test$dataset,
                        comparison = comparison,
                        output.filename = host.data.for_test$output.filename,
                        custom.md = host.data.for_test$metadata,
                        custom.levels = host.data.for_test$custom.levels,
                        ref.level = ref.level,
                        inside.host = inside.host
                     )

#' Save the fit data object as an rds file.
# 20260220_21_13_07
saveRDS(maaslin.fit_data.host,
        file = file.path(diffabund.rdafiles,
                         paste(
                           "maaslin.fit_data.host",
                           paste0(host.data.for_test$short.filename,".rds"),sep = "-")))

#' Process the output 
# "20260213_12_20_50" for signif.tsv
# "20260213_12_20_51" for maaslin.signif.decreased.tsv
maaslin.processed_output.host <- 
  process_maaslin2_output(agglom.rank = agglom.rank,
                          short.filename = host.data.for_test$short.filename,
                          maaslin.fit_data = maaslin.fit_data.host,
                          ps.q.agg = host.data.for_test$ps.q.agg,
                          sample.groups = host.data.for_test$sample.groups)

#' Save the workspace.
#' 20260213_12_20_28: genus workspace
save.image(file.path(diffabund.rdafiles,paste(
  "maaslin2",host.data.for_test$short.filename,"workspace.RData",sep="-")))

#+ echo=FALSE
### 4.2 Run a test with ALDEx2. ####
#'
#' ### Run a test with ALDEx2.
#' ### Create covariates: take the samples in the ps.q.df.wide, and find 
#' the positions of these samples in the custom.md. 
#' Ref: match returns a vector of the positions of (first) matches of its first argument in its second. 
#' We're searching the elements of the first vector in the second vector. 
if(setequal(host.data.for_test$custom.levels,"NMR")& comparison=="age"){
  covariates<-
    host.data.for_test$metadata$agegroup[
      match(rownames(host.data.for_test$dataset),
            rownames(host.data.for_test$metadata))]
}else if (inside.host==FALSE){
  covariates<-
    host.data.for_test$metadata$class[
    match(rownames(host.data.for_test$dataset),
          rownames(host.data.for_test$metadata))]
}
#' Model matrix is covariates and the reference group is the first covariate.
mm <- model.matrix(~ covariates-1)
#' Reorder model.matrix to put ref.level as first column.
if(inside.host==TRUE){
  if (comparison=="age"){
    aldex.reference<-paste0("covariates",ref.level)
  }else if(comparison=="sex"){
    aldex.reference<-"covaritatesF"
  }else if(comparison=="strain"){
    aldex.reference<-paste0("covariates",ref.level)
  }
  # Reorder columns of mm so that the first column is aldex.reference column and
  # all the other columns are after. 
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == aldex.reference), 
                          which(colnames(mm) !=aldex.reference))]]
  
}else{
  # Same reordering except this is the case when we have animal hosts. We 
  # directly put ref.level into the string.
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == paste0("covariates",ref.level)), 
                          which(colnames(mm) != paste0("covariates",ref.level)))]]
}


#' Run the test: Aldex glm for a complex case.
set.seed(1)
ps.q.aldex.clr <- aldex.clr(t(host.data.for_test$dataset), mm, 
                            mc.samples=1000, 
                            denom="all", verbose=F)
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
#' Save the workspace just in case.
save.image(file.path(diffabund.rdafiles,
                  paste("aldex2",host.data.for_test$short.filename,
                        "workspace-test.RData",sep="-")))
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
#' Save the workspace just in case.
save.image(file.path(diffabund.rdafiles,
                  paste("aldex2",host.data.for_test$short.filename,
                        "workspace-effect.RData",sep = "-")))


#' Extract significant features.
aldex.signif.features<-list()
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
    
    if(inside.host==TRUE){
      signif.df$OTU<-rownames(signif.df)
    }else{
      signif.df$Taxon<-rownames(signif.df)
    }
    
    # signif.df$class<-names(ps.q.glm.effect[i])
    rownames(signif.df)<-1:nrow(signif.df)
    aldex.signif.features[[names(ps.q.glm.effect[i])]]<-signif.df
  } else{
    signif.df<-data.frame()
    aldex.signif.features[[names(ps.q.glm.effect[i])]]<-signif.df
  }
  
}
rm(signif.df)
aldex.signif.features<-bind_rows(aldex.signif.features,.id = "class")
#' Check for duplicates.
if(inside.host==TRUE){
  aldex.signif.features%>%
    dplyr::group_by(OTU)%>%
    summarise(n=n())%>%
    arrange(-n) 
}else{
  aldex.signif.features%>%
    dplyr::group_by(Taxon)%>%
    summarise(n=n())%>%
    arrange(-n)
}
aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)
head(aldex.neg.effect)
#' Save the significant features with negative effect size as a 
#'  tab-separated file.
# "20260213_13_19_00"
write.table(aldex.neg.effect,
            file=file.path(diffabund.tables,
                           paste("aldex.neg.effect",
                                 host.data.for_test$short.filename,
                                 "signif.tsv",sep="-")),
            row.names = F,sep = "\t")
#' Save the workspace just in case.
save.image(file.path(diffabund.rdafiles,paste(
  "aldex2",host.data.for_test$short.filename,"workspace.RData",sep="-")))
#' Save all significant features.
# "20260213_13_20_05"
write.table(aldex.signif.features,
            file=file.path(diffabund.tables,paste(
              "aldex2",host.data.for_test$short.filename,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

#+ echo=FALSE
### 4.3 Run a test with ANCOMBC. ####
#'
#' ### Run a test with ANCOMBC.
ancombc.out<-perform_ancombc_test(agglom.rank = agglom.rank,
                                  ps.q.df.wide = host.data.for_test$dataset,
                                  custom.md = host.data.for_test$metadata,
                                  sample.groups = host.data.for_test$sample.groups,
                                  ps.q.agg = host.data.for_test$ps.q.agg)
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
head(ancombc.signif.features)
head(ancombc.signif.decreased)
#' Save the workspace.
save.image(file.path(diffabund.rdafiles,paste(
  "ancombc", host.data.for_test$short.filename, "workspace.RData",sep="-")))
#' Save all significant features. 
"20260213_13_41_58"
write.table(ancombc.signif.features,
            file=file.path(diffabund.tables,paste(
              "ancombc",host.data.for_test$short.filename,
              "signif.tsv",sep="-")),
            row.names = F,sep = "\t")
#' Save decreased significant features.
# "20260213_13_42_00"
write.table(ancombc.signif.decreased,
            file=file.path(diffabund.tables,
                           paste("ancombc.signif.decreased",
                                 host.data.for_test$short.filename,
                                 "signif.tsv",sep="-")),
            row.names = F,sep = "\t")

#+ echo=FALSE
### 4.4 Downstream analysis and plotting of differentially abundant features. ####
#' 
#' ### Downstream analysis and plotting of differentially abundant features. 
#+ fig.width=11, fig.height=8
analyse_test_output(agglom.rank = agglom.rank,
                    inside.host = inside.host, 
                    comparison = comparison, 
                    custom.levels = host.data.for_test$custom.levels,
                    ref.level = ref.level, 
                    ps.q.agg = host.data.for_test$ps.q.agg,
                    sample.groups = host.data.for_test$sample.groups,
                    output.filename = host.data.for_test$short.filename,
                    custom.md = host.data.for_test$metadata,
                    maaslin.signif.features = maaslin.processed_output.host$maaslin.signif.features,
                    maaslin.signif.decreased = maaslin.processed_output.host$maaslin.signif.decreased,
                    aldex.signif.features = aldex.signif.features,
                    aldex.neg.effect = aldex.neg.effect,
                    ancombc.signif.features = ancombc.signif.features,
                    ancombc.signif.decreased = ancombc.signif.decreased)
  
#+ echo=FALSE
## 5. Compare NMR age groups (ASV level). ####
#' 
#' ## Compare NMR age groups (ASV level).
#' Choose what to compare:
comparison<-"age"
#' Choose the reference level:
ref.level<-"agegroup0_10" 
agglom.rank<-"OTU"
inside.host=TRUE
#' Prepare the data for testing.
nmr.age.data.for_test<-prepare_data(agglom.rank = agglom.rank, 
             comparison = comparison, 
             ref.level = ref.level, 
             inside.host = TRUE)
#' Run MaAsLin2 test.
maaslin.fit_data.age<-
  perform_maaslin2_test(ps.q.df.wide = nmr.age.data.for_test$dataset,
                        comparison = comparison,
                        output.filename = nmr.age.data.for_test$output.filename,
                        custom.md = nmr.age.data.for_test$metadata,
                        custom.levels = nmr.age.data.for_test$custom.levels,
                        ref.level = ref.level,
                        inside.host = inside.host
  )

#' Save the fit data object as an rds file.
# 20260224_21_19_15
saveRDS(maaslin.fit_data.age,
        file = file.path(diffabund.rdafiles,
                         paste(
                           "maaslin.fit_data.age",
                           nmr.age.data.for_test$output.filename,
                           ".rds",sep = "-")))
# "20260213_14_39_57" for both signif and signif.decreased tsv
maaslin.processed_output.age <- 
  process_maaslin2_output(agglom.rank = agglom.rank,
                          short.filename = nmr.age.data.for_test$output.filename,
                          maaslin.fit_data = maaslin.fit_data.age,
                          ps.q.agg = nmr.age.data.for_test$ps.q.agg,
                          sample.groups = nmr.age.data.for_test$sample.groups)

#' Downstream analysis and plotting of differentially abundant features.
#+ fig.width=8, fig.height=6
analyse_test_output(agglom.rank = agglom.rank,
                    inside.host = inside.host, 
                    comparison = comparison, 
                    custom.levels = nmr.age.data.for_test$custom.levels,
                    ref.level = ref.level, 
                    ps.q.agg = nmr.age.data.for_test$ps.q.agg,
                    sample.groups = nmr.age.data.for_test$sample.groups,
                    output.filename = nmr.age.data.for_test$output.filename,
                    custom.md = nmr.age.data.for_test$metadata,
                    maaslin.signif.features = maaslin.processed_output.age$maaslin.signif.features,
                    maaslin.signif.decreased = maaslin.processed_output.age$maaslin.signif.decreased,
                    aldex.signif.features = NULL,
                    aldex.neg.effect = NULL,
                    ancombc.signif.features = NULL,
                    ancombc.signif.decreased = NULL)

#+ echo=FALSE
## 6. Compare NMR sexes (ASV level). ####
#' 
#' ## Compare NMR sexes (ASV level).
#' Choose what to compare:
comparison<-"sex"
agglom.rank<-"OTU"
#' Choose the reference level:
ref.level<-"female"
inside.host = TRUE
#' Prepare the data for testing.
nmr.sex.data.for_test<-prepare_data(agglom.rank = agglom.rank, 
             comparison = comparison, 
             ref.level = ref.level, 
             inside.host = TRUE)
#' Run MaAsLin2 test.
maaslin.fit_data.sex<-
  perform_maaslin2_test(ps.q.df.wide = nmr.sex.data.for_test$dataset,
                        comparison = comparison,
                        output.filename = nmr.sex.data.for_test$output.filename,
                        custom.md = nmr.sex.data.for_test$metadata,
                        custom.levels = nmr.sex.data.for_test$custom.levels,
                        ref.level = ref.level,
                        inside.host = inside.host
  )

#' Save the fit data object as an rds file.
# 20260224_21_21_02
saveRDS(maaslin.fit_data.sex,
        file = file.path(diffabund.rdafiles,
                         paste(
                           "maaslin.fit_data.sex",
                           nmr.sex.data.for_test$output.filename,".rds",sep = "-")))
# "20260215_19_30_37" for both signif and signif.decreased tsv
maaslin.processed_output.sex <-
  process_maaslin2_output(agglom.rank = agglom.rank,
                          maaslin.fit_data = maaslin.fit_data.sex,
                          short.filename = nmr.sex.data.for_test$output.filename,
                          ps.q.agg = nmr.sex.data.for_test$ps.q.agg,
                          sample.groups = nmr.sex.data.for_test$sample.groups)
#' Downstream analysis and plotting of differentially abundant features. 
#+ fig.width=8, fig.height=20
analyse_test_output(agglom.rank = agglom.rank,
                    inside.host = inside.host, 
                    comparison = comparison, 
                    custom.levels = nmr.sex.data.for_test$custom.levels,
                    ref.level = ref.level, 
                    ps.q.agg = nmr.sex.data.for_test$ps.q.agg,
                    sample.groups = nmr.sex.data.for_test$sample.groups,
                    output.filename = nmr.sex.data.for_test$output.filename,
                    custom.md = nmr.sex.data.for_test$metadata,
                    maaslin.signif.features = maaslin.processed_output.sex$maaslin.signif.features,
                    maaslin.signif.decreased = maaslin.processed_output.sex$maaslin.signif.decreased,
                    aldex.signif.features = NULL,
                    aldex.neg.effect = NULL,
                    ancombc.signif.features = NULL,
                    ancombc.signif.decreased = NULL)

ps.q.agg.genus.all_domains<-readRDS("analysis/1-qiime2/20260209_16_33_25-single-234/01-community_composition/all_domains/rdafiles/phyloseq-qiime-pooled-Genus-single-234-table.rds")

ps.q.agg.genus.bacteria<-readRDS("analysis/1-qiime2/20260209_16_33_25-single-234/01-community_composition/bacteria_only-old/rdafiles/phyloseq-qiime-pooled-Genus-single-234-table.rds")

rare.genus.all_domains<-read.table("analysis/1-qiime2/20260209_16_33_25-single-234/01-community_composition/all_domains/tables/ps.q.df.rare-nonfiltered-Genus-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit.tsv",header = T)%>%
  as_tibble()

rare.genus.bacteria<-read.table("analysis/1-qiime2/20260209_16_33_25-single-234/01-community_composition/bacteria_only-old/tables/renamed/ps.q.df.rare-nonfiltered-Genus-all_hosts.tsv",header = T)%>%
  as_tibble()

all.domains.maaslin.signif<-read.table(file = "analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/all_domains/tables/maaslin2-all_hosts-Genus-host-ref-NMR-signif.tsv",header = T)%>%
  as_tibble()
all.domains.maaslin.signif.decreased<-read.table(file = "analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/all_domains/tables/maaslin.signif.decreased-all_hosts-Genus-host-ref-NMR-signif.tsv",header = T)%>%
  as_tibble()

all.domains.signif.all<-read.table(file ="analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/all_domains/tables/significant-features-all_hosts-Genus-host-ref-NMR-signif.tsv", header= T)%>%
  as_tibble()

all.domains.ancom <-read.table(file ="analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/all_domains/tables/ancombc.signif.decreased-all_hosts-Genus-host-ref-NMR-signif.tsv", header= T)%>%
  as_tibble()

all.domains.maaslin.signif%>%arrange(feature)
all.domains.maaslin.signif.decreased%>%arrange(feature)
all.domains.signif.all

bacteria.only.maaslin.signif<-read.table(file = "analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/bacteria_only-old/tables/renamed/maaslin2-all_hosts-Genus-host-234-ref-NMR-signif.tsv",header = T)%>%
  as_tibble()

bacteria.only.maaslin.signif.decreased<-read.table(file = "analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/bacteria_only-old/tables/renamed/maaslin.signif.decreased-all_hosts-Genus-host-234-ref-NMR-signif.tsv",header = T)%>%
  as_tibble()

bacteria.only.signif.all<-read.table(file = "analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/bacteria_only-old/tables/renamed/significant-features-all_hosts-Genus-host-234-ref-NMR-signif.tsv",header = T)%>%
  as_tibble()
bacteria.only.ancom <-read.table(file = "analysis/1-qiime2/20260209_16_33_25-single-234/03-differential_abundance/bacteria_only-old/tables/renamed/ancombc.signif.decreased-all_hosts-Genus-host-234-ref-NMR-signif.tsv",header = T)%>%
  as_tibble()

bacteria.only.maaslin.signif.decreased%>%group_by(feature)%>%summarise(nf=n())%>%summarise(nf==8)%>%table()
all.domains.maaslin.signif.decreased%>%group_by(feature)%>%summarise(nf=n())%>%summarise(nf==8)%>%table()


# small discrepancy 
setdiff(unique(bacteria.only.maaslin.signif.decreased$feature),
        unique(all.domains.maaslin.signif.decreased$feature))
setdiff(unique(all.domains.maaslin.signif.decreased$feature),
        unique(bacteria.only.maaslin.signif.decreased$feature))

setdiff(unique(bacteria.only.ancom$taxon_id),
        unique(all.domains.ancom$taxon_id))
setdiff(unique(all.domains.ancom$taxon_id),
        unique(bacteria.only.ancom$taxon_id))

setdiff(all.domains.signif.all$Genus,
        bacteria.only.signif.all$Genus)
setdiff(bacteria.only.signif.all$Genus,
        all.domains.signif.all$Genus
  )


# Check
all.domains.maaslin.signif%>%
  filter(feature =="Acetanaerobacterium")

bacteria.only.maaslin.signif.decreased%>%
  filter(feature =="Acetanaerobacterium")

# No difference in the general table
ps.q.agg.genus.all_domains%>%
  filter(Genus =="Acetanaerobacterium")%>%
  distinct(class,Genus,MeanRelativeAbundance)%>%
  arrange(Genus,class)
ps.q.agg.genus.bacteria%>%
  filter(Genus =="Acetanaerobacterium")%>%
  distinct(class,Genus,MeanRelativeAbundance)%>%
  arrange(Genus,class)

ps.q.agg.genus.all_domains%>%
  filter(Genus %in%c("Methanobrevibacter", "Paludicola"  ))%>%
  distinct(class,Genus,MeanRelativeAbundance)%>%
  arrange(Genus,class)
ps.q.agg.genus.bacteria%>%
  filter(Genus %in%c("Methanobrevibacter", "Paludicola"  ))%>%
  distinct(class,Genus,MeanRelativeAbundance)%>%
  arrange(Genus,class)


# but in rarefied data, these bacteria are more abundant
rare.genus.all_domains%>%
  filter(Genus =="Acetanaerobacterium")
rare.genus.bacteria%>%
  filter(Genus =="Acetanaerobacterium")

rare.genus.all_domains%>%
  filter(Genus %in%c("Methanobrevibacter", "Paludicola"  ))%>%
  group_by(class, Genus)%>%
  summarise(MeanAbundance = mean(Abundance))%>%
  arrange(Genus,class)
rare.genus.bacteria%>%
  filter(Genus %in%c("Methanobrevibacter", "Paludicola"  ))%>%
  group_by(class, Genus)%>%
  summarise(MeanAbundance = mean(Abundance))%>%
  arrange(Genus,class)



sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), "active.analysis"))
gc()