---
output: 
  bookdown::html_document2:
     toc: true
---




# Differential microbial abundance tests with MaAsLin2, ALDEx2, and ANCOM-BC




## Introduction
This script performs differential microbial abundance tests on 
different hosts. We will use MaAsLin2, ALDEx2, and ANCOM-BC.




## Load necessary libraries.


``` r
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# Current: Bioconductor version 3.20 (BiocManager 1.30.27), R 4.4.3 (2025-02-28 ucrt)
```

ALDEx2: 1.38.0
Maaslin2: 1.20.0
ANCOM-BC: 2.8.1  


``` r
# BiocManager::install(c("Maaslin2","ALDEx2","ANCOMBC","phyloseq"), version = "3.20")
# BiocManager::install("ALDEx2", version = "3.17")
# BiocManager::install("ANCOMBC", version = "3.17")
# BiocManager::install("phyloseq", version = "3.17")
# install.packages(c("tidyverse","Polychrome"))
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   4.0.0     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.0.4     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ purrr::%||%()   masks base::%||%()
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(Maaslin2)
library(ALDEx2)
```

```
## Loading required package: zCompositions
## Loading required package: MASS
## 
## Attaching package: 'MASS'
## 
## The following object is masked from 'package:dplyr':
## 
##     select
## 
## Loading required package: truncnorm
## Loading required package: survival
## Loading required package: lattice
## Loading required package: latticeExtra
## 
## Attaching package: 'latticeExtra'
## 
## The following object is masked from 'package:ggplot2':
## 
##     layer
```

``` r
library(ANCOMBC)
library(phyloseq)
library(Polychrome)
```



## Specifying parameters and directory/file names. 
Name of the folder with QIIME2 output:


``` r
authorname<-"pooled"
```

Directories with input files:


``` r
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)
```

Truncation level that we chose in QIIME2:


``` r
truncationlvl<-"234" 
```

Single reads or paired reads (decided in QIIME2):


``` r
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
```



## Prepare necessary functions. ####
Function that performs differential abundance tests with different 
parameters. Instead of repeating the code for each comparison,
call the function and change the parameters.


``` r
prepare_data <- function(agglom.rank, comparison, 
                         ref.level, inside.host){
  rdafiles.directory<-"./output/rdafiles"
  rare.status<-"rare"
  filter.status<-"nonfiltered"
  print(paste("Performing differential microbial abundance tests on", comparison, "variable at", agglom.rank, "level.",
              "Reference:", ref.level))
  #' Import abundance table as an rds file (NOT rarefied):
  if(inside.host == TRUE){
    if(agglom.rank == "OTU"){
      ps.q.agg.date_time<-"20260211_17_01_07" # ASV NMR
      ps.q.df.preprocessed.date_time <- "20260211_17_14_21"
    }
    custom.levels<-"NMR"
    custom.md<-readRDS("./output/rdafiles/custom.md.ages.rds")%>%
      filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
  }else if(inside.host ==FALSE){
    if(agglom.rank == "Genus" ){
      ps.q.agg.date_time<-"20260211_17_01_10" # Genus all hosts
      custom.levels<-c("NMR",
                       "DMR",
                       "B6mouse",
                       "MSMmouse",
                       "FVBNmouse",
                       "spalax",
                       "pvo",
                       "hare",
                       "rabbit")
      custom.md<-readRDS("./output/rdafiles/custom.md.rds")
      ps.q.df.preprocessed.date_time <- "20260211_17_14_19"
      
    } 
  }
  #' Import abundance table as an rds file (NOT rarefied):
  ps.q.agg<-readRDS(file=file.path(
    rdafiles.directory,
    paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,read.end.type,
          truncationlvl,"table.rds",sep = "-")))
  ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
  
  #' Rarefied abundance table:
  ps.q.df.preprocessed<-readRDS(
    file.path(rdafiles.directory,paste0(
      paste(
        ps.q.df.preprocessed.date_time,
        paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
        paste(custom.levels,collapse = '-'),sep = "-"),
      ".rds")))
  if(inside.host ==TRUE){
    ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
      filter(class=="NMR",Abundance!=0)
  }
  
  #' Specify output file name:
  output.filename<-paste(paste(custom.levels,collapse = '-'),
                         agglom.rank,comparison,
                         truncationlvl,"ref",ref.level,
                         sep="-")
  print(paste("Output filename prefix:", output.filename))
  # Preparing the wide dataset.
  ps.q.df.wide<-ps.q.df.preprocessed%>%
    dplyr::select(all_of(c("Sample",agglom.rank,"Abundance","class")))%>%
    filter(Abundance!=0)%>%
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
              metadata = custom.md, 
              custom.levels = custom.levels,
              sample.groups = sample.groups))
}
```

Function that performs a test with MaAsLin2:


``` r
perform_maaslin2_test <- function(ps.q.df.wide, 
                                 output.filename,
                                 comparison,
                                 custom.md,
                                 custom.levels,
                                 ref.level,
                                 inside.host){
  authorname <-"pooled"
  rare.status<-"rare"
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
               output = file.path("./output/maaslin2",
                                  paste0(authorname,"-output"),
                                  rare.status,paste(
                                    paste(format(Sys.time(),format="%Y%m%d"),
                                          format(Sys.time(),
                                                 format = "%H_%M_%S"),sep = "_"),
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
               output = file.path("./output/maaslin2",
                                  paste0(authorname,"-output"),
                                  rare.status,paste(
                                    paste(format(Sys.time(),format="%Y%m%d"),
                                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                    output.filename,sep = "-")), 
               fixed_effects = c("class"),
               reference = paste0("class,",ref.level),
               max_significance = 0.05)
    
  }
  return(maaslin.fit_data)
}
```

Function to process the MaAsLin2 output:


``` r
process_maaslin2_output<-function(agglom.rank,
                                  maaslin.fit_data,
                                  ps.q.agg,
                                  sample.groups){
  #' Downstream processing of Maaslin2 output.
  source("./code/r-scripts/make_features_maaslin.R")
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
  # write.table(maaslin.signif.features,
  #             file=file.path("./output/rtables/pooled",paste(
  #               paste(format(Sys.time(),format="%Y%m%d"),
  #                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #               "maaslin2",output.filename,
  #               "signif.tsv",sep="-")),
  #             row.names = F,sep = "\t")
  # write.table(maaslin.signif.decreased,
  #             file=file.path("./output/rtables/pooled",
  #                            paste(paste(format(Sys.time(),format="%Y%m%d"),
  #                                        format(Sys.time(),format = "%H_%M_%S"),
  #                                        sep = "_"),
  #                                  "maaslin.signif.decreased",
  #                                  output.filename,
  #                                  "signif.tsv",sep="-")),
  #             row.names = F,sep = "\t")
  return(list(maaslin.signif.features = maaslin.signif.features,
         maaslin.signif.decreased = maaslin.signif.decreased))
  
}
```

Function that performs a test with ANCOM-BC:


``` r
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
```

Function that performs downstream processing of output from three tools 
(if you compare between hosts) or one tool (if you compare within host).


``` r
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
  image.formats<-c("png","tiff")
  rtables.directory<-file.path("./output/rtables",authorname)
  
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
    # write.table(common.decreased.df,
    #             file=file.path(rtables.directory,
    #                            paste(paste(format(Sys.time(),format="%Y%m%d"),
    #                                        format(Sys.time(),format = "%H_%M_%S"),
    #                                        sep = "_"),
    #                                  "significant-features",
    #                                  output.filename,
    #                                  "signif.tsv",sep="-")),
    #             row.names = F,sep = "\t")
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
  
  plot.theme <-theme(
    axis.title.y = element_blank(),
    axis.title = element_text(size = 5),
    axis.text.y = ggtext::element_markdown(size=5),
    # axis.text.y = element_text(size=5),
    axis.text.x = element_text(size=5),
    strip.text.x = ggtext::element_markdown(size=5),
    plot.title = element_text(size =5),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
  
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
      plot.theme
    
    
    # for(image.format in image.formats){
    #   ggsave(paste0("./images/taxaboxplots/",
    #                 paste(paste(format(Sys.time(),format="%Y%m%d"),
    #                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
    #                       ref.level,"specific-bacteria","NMR",comparison,
    #                       sep = "-"),".",image.format),
    #          plot=feature.plot,
    #          width=11, height=8,units="in",
    #          # width = 4500,height = 2500,units = "px",
    #          dpi=300,device = image.format)
    # }
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
      plot.theme
    
    
    # for(image.format in image.formats){
    #   ggsave(paste0("./images/taxaboxplots/",
    #                 paste(paste(format(Sys.time(),format="%Y%m%d"),
    #                             format(Sys.time(),format = "%H_%M_%S"),
    #                             sep = "_"),
    #                       "NMR-specific-bacteria",
    #                       sep = "-"),".",image.format),
    #          plot=feature.plot,
    #          width=8, height=11,units="in",
    #          # width = 4000,height = 12000,units = "px",
    #          dpi=300,device = image.format)
    # }
    
  }
  print(feature.plot)
  
  
}
```



## Compare hosts (Genus level).
Choose what to compare:


``` r
comparison<-"host"
```

Choose the reference level:


``` r
ref.level<-"NMR" 
agglom.rank<-"Genus"
inside.host<- FALSE
```

Prepare data for tests.


``` r
host.data.for_test<-prepare_data(agglom.rank = agglom.rank, 
             comparison = comparison, 
             ref.level = ref.level, 
             inside.host = inside.host)
```

```
## [1] "Performing differential microbial abundance tests on host variable at Genus level. Reference: NMR"
## [1] "Output filename prefix: NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR"
```

```
## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
## ℹ Please use `all_of()` or `any_of()` instead.
##   # Was:
##   data %>% select(agglom.rank)
## 
##   # Now:
##   data %>% select(all_of(agglom.rank))
## 
## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```
## [1] "NMR"       "DMR"       "B6mouse"   "MSMmouse"  "FVBNmouse" "spalax"   
## [7] "pvo"       "hare"      "rabbit"   
## [1] "Prepared the data for tests"
```



### Run a test with MaAsLin2.


``` r
maaslin.fit_data.host<-
  perform_maaslin2_test(ps.q.df.wide = host.data.for_test$dataset,
                        comparison = comparison,
                        output.filename = host.data.for_test$output.filename,
                        custom.md = host.data.for_test$metadata,
                        custom.levels = host.data.for_test$custom.levels,
                        ref.level = ref.level,
                        inside.host = inside.host
                     )
```

```
## [1] "Creating output folder"
## [1] "Creating output feature tables folder"
## [1] "Creating output fits folder"
## [1] "Creating output figures folder"
## 2026-03-04 18:15:55.35137 INFO::Writing function arguments to log file
## 2026-03-04 18:15:55.404847 INFO::Verifying options selected are valid
## 2026-03-04 18:15:55.433154 INFO::Determining format of input files
## 2026-03-04 18:15:55.437867 INFO::Input format is data samples as rows and metadata samples as rows
## 2026-03-04 18:15:55.450808 INFO::Formula for fixed effects: expr ~  class
## 2026-03-04 18:15:55.454911 INFO::Filter data based on min abundance and min prevalence
## 2026-03-04 18:15:55.460213 INFO::Total samples in data: 99
## 2026-03-04 18:15:55.463759 INFO::Min samples required with min abundance for a feature not to be filtered: 0.000000
## 2026-03-04 18:15:55.470459 INFO::Total filtered features: 0
## 2026-03-04 18:15:55.474859 INFO::Filtered feature names from abundance and prevalence filtering:
## 2026-03-04 18:15:55.483914 INFO::Total filtered features with variance filtering: 0
## 2026-03-04 18:15:55.487811 INFO::Filtered feature names from variance filtering:
## 2026-03-04 18:15:55.491301 INFO::Running selected normalization method: TSS
## 2026-03-04 18:15:55.8814 INFO::Bypass z-score application to metadata
## 2026-03-04 18:15:55.887475 INFO::Running selected transform method: LOG
## 2026-03-04 18:15:55.896752 INFO::Running selected analysis method: LM
## 2026-03-04 18:15:55.905246 INFO::Fitting model to feature number 1, Clostridium_sensu_stricto_1
## 2026-03-04 18:15:55.913145 INFO::Fitting model to feature number 2, Lachnospiraceae.Family
## 2026-03-04 18:15:55.917803 INFO::Fitting model to feature number 3, Ruminococcus
## 2026-03-04 18:15:55.923343 INFO::Fitting model to feature number 4, Akkermansia
## 2026-03-04 18:15:55.929355 INFO::Fitting model to feature number 5, Escherichia.Shigella
## 2026-03-04 18:15:55.934478 INFO::Fitting model to feature number 6, Bacteroidales.Order
## 2026-03-04 18:15:55.938964 INFO::Fitting model to feature number 7, X.Eubacterium._coprostanoligenes_group
## 2026-03-04 18:15:55.94367 INFO::Fitting model to feature number 8, Clostridia_UCG.014
## 2026-03-04 18:15:55.949131 INFO::Fitting model to feature number 9, Bacteroides
## 2026-03-04 18:15:55.954539 INFO::Fitting model to feature number 10, NK4A214_group
## 2026-03-04 18:15:55.95909 INFO::Fitting model to feature number 11, Christensenellaceae_R.7_group
## 2026-03-04 18:15:55.963774 INFO::Fitting model to feature number 12, Eubacteriaceae.Family
## 2026-03-04 18:15:55.968368 INFO::Fitting model to feature number 13, UCG.010
## 2026-03-04 18:15:55.97281 INFO::Fitting model to feature number 14, Alistipes
## 2026-03-04 18:15:55.978269 INFO::Fitting model to feature number 15, Rikenellaceae_RC9_gut_group
## 2026-03-04 18:15:55.983098 INFO::Fitting model to feature number 16, Monoglobus
## 2026-03-04 18:15:55.988081 INFO::Fitting model to feature number 17, Firmicutes.Phylum
## 2026-03-04 18:15:55.992735 INFO::Fitting model to feature number 18, dgA.11_gut_group
## 2026-03-04 18:15:55.998366 INFO::Fitting model to feature number 19, Clostridia_vadinBB60_group
## 2026-03-04 18:15:56.003471 INFO::Fitting model to feature number 20, Ruminococcaceae.Family
## 2026-03-04 18:15:56.008093 INFO::Fitting model to feature number 21, V9D2013_group
## 2026-03-04 18:15:56.012516 INFO::Fitting model to feature number 22, Bacteria.Kingdom
## 2026-03-04 18:15:56.017128 INFO::Fitting model to feature number 23, Oscillospiraceae.Family
## 2026-03-04 18:15:56.021847 INFO::Fitting model to feature number 24, Subdoligranulum
## 2026-03-04 18:15:56.026463 INFO::Fitting model to feature number 25, Lachnospiraceae_NK4A136_group
## 2026-03-04 18:15:56.03099 INFO::Fitting model to feature number 26, Gastranaerophilales
## 2026-03-04 18:15:56.035477 INFO::Fitting model to feature number 27, Candidatus_Saccharimonas
## 2026-03-04 18:15:56.039991 INFO::Fitting model to feature number 28, Barnesiellaceae.Family
## 2026-03-04 18:15:56.044434 INFO::Fitting model to feature number 29, Ruminococcaceae
## 2026-03-04 18:15:56.048878 INFO::Fitting model to feature number 30, Tyzzerella
## 2026-03-04 18:15:56.053429 INFO::Fitting model to feature number 31, RF39
## 2026-03-04 18:15:56.058008 INFO::Fitting model to feature number 32, Incertae_Sedis
## 2026-03-04 18:15:56.063721 INFO::Fitting model to feature number 33, UCG.005
## 2026-03-04 18:15:56.070076 INFO::Fitting model to feature number 34, X.Eubacterium._siraeum_group
## 2026-03-04 18:15:56.075221 INFO::Fitting model to feature number 35, Colidextribacter
## 2026-03-04 18:15:56.079944 INFO::Fitting model to feature number 36, Hungateiclostridiaceae.Family
## 2026-03-04 18:15:56.085547 INFO::Fitting model to feature number 37, Eggerthellaceae.Family
## 2026-03-04 18:15:56.090566 INFO::Fitting model to feature number 38, Coriobacteriales.Order
## 2026-03-04 18:15:56.095061 INFO::Fitting model to feature number 39, Enterorhabdus
## 2026-03-04 18:15:56.099488 INFO::Fitting model to feature number 40, Clostridia.Class
## 2026-03-04 18:15:56.104009 INFO::Fitting model to feature number 41, Rhodospirillales.Order
## 2026-03-04 18:15:56.108831 INFO::Fitting model to feature number 42, Lachnospiraceae_NK4B4_group
## 2026-03-04 18:15:56.113828 INFO::Fitting model to feature number 43, Eubacterium
## 2026-03-04 18:15:56.118294 INFO::Fitting model to feature number 44, Eggerthella
## 2026-03-04 18:15:56.122679 INFO::Fitting model to feature number 45, Coriobacteriales_Incertae_Sedis.Family
## 2026-03-04 18:15:56.12704 INFO::Fitting model to feature number 46, DNF00809
## 2026-03-04 18:15:56.132037 INFO::Fitting model to feature number 47, Butyricicoccus
## 2026-03-04 18:15:56.136623 INFO::Fitting model to feature number 48, Atopobiaceae.Family
## 2026-03-04 18:15:56.142707 INFO::Fitting model to feature number 49, Family_XIII_AD3011_group
## 2026-03-04 18:15:56.149632 INFO::Fitting model to feature number 50, UCG.004
## 2026-03-04 18:15:56.154938 INFO::Fitting model to feature number 51, Christensenellaceae.Family
## 2026-03-04 18:15:56.160057 INFO::Fitting model to feature number 52, Ruminiclostridium
## 2026-03-04 18:15:56.167461 INFO::Fitting model to feature number 53, Natranaerobiaceae.Family
## 2026-03-04 18:15:56.173856 INFO::Fitting model to feature number 54, Anaerovoracaceae.Family
## 2026-03-04 18:15:56.180327 INFO::Fitting model to feature number 55, Peptococcaceae.Family
## 2026-03-04 18:15:56.186584 INFO::Fitting model to feature number 56, UCG.007
## 2026-03-04 18:15:56.191398 INFO::Fitting model to feature number 57, Family_XIII_UCG.001
## 2026-03-04 18:15:56.198268 INFO::Fitting model to feature number 58, Eubacteriales.Order
## 2026-03-04 18:15:56.203831 INFO::Fitting model to feature number 59, Coriobacteriaceae.Family
## 2026-03-04 18:15:56.20855 INFO::Fitting model to feature number 60, Paenibacillus
## 2026-03-04 18:15:56.213787 INFO::Fitting model to feature number 61, X.Clostridium._methylpentosum_group
## 2026-03-04 18:15:56.2186 INFO::Fitting model to feature number 62, Flavonifractor
## 2026-03-04 18:15:56.223272 INFO::Fitting model to feature number 63, Clostridiaceae.Family
## 2026-03-04 18:15:56.228028 INFO::Fitting model to feature number 64, X.Eubacterium._ventriosum_group
## 2026-03-04 18:15:56.232708 INFO::Fitting model to feature number 65, Roseburia
## 2026-03-04 18:15:56.238796 INFO::Fitting model to feature number 66, Clostridium_sensu_stricto_13
## 2026-03-04 18:15:56.244715 INFO::Fitting model to feature number 67, Enterococcus
## 2026-03-04 18:15:56.250156 INFO::Fitting model to feature number 68, Erysipelotrichaceae.Family
## 2026-03-04 18:15:56.255106 INFO::Fitting model to feature number 69, Erysipelatoclostridium
## 2026-03-04 18:15:56.25983 INFO::Fitting model to feature number 70, Bacillus
## 2026-03-04 18:15:56.264615 INFO::Fitting model to feature number 71, CAG.352
## 2026-03-04 18:15:56.269646 INFO::Fitting model to feature number 72, Desulfovibrio
## 2026-03-04 18:15:56.275792 INFO::Fitting model to feature number 73, Papillibacter
## 2026-03-04 18:15:56.281131 INFO::Fitting model to feature number 74, Acinetobacter
## 2026-03-04 18:15:56.28581 INFO::Fitting model to feature number 75, Oscillospirales.Order
## 2026-03-04 18:15:56.290495 INFO::Fitting model to feature number 76, Izemoplasmatales
## 2026-03-04 18:15:56.295251 INFO::Fitting model to feature number 77, Lachnospiraceae_NC2004_group
## 2026-03-04 18:15:56.300192 INFO::Fitting model to feature number 78, Oxalobacter
## 2026-03-04 18:15:56.305713 INFO::Fitting model to feature number 79, Candidatus_Soleaferrea
## 2026-03-04 18:15:56.310581 INFO::Fitting model to feature number 80, Sporobacter
## 2026-03-04 18:15:56.315548 INFO::Fitting model to feature number 81, Defluviitaleaceae_UCG.011
## 2026-03-04 18:15:56.320507 INFO::Fitting model to feature number 82, Tuzzerella
## 2026-03-04 18:15:56.325688 INFO::Fitting model to feature number 83, Pygmaiobacter
## 2026-03-04 18:15:56.330488 INFO::Fitting model to feature number 84, Pelomonas
## 2026-03-04 18:15:56.335223 INFO::Fitting model to feature number 85, Methylobacterium.Methylorubrum
## 2026-03-04 18:15:56.339761 INFO::Fitting model to feature number 86, Sphingomonadaceae.Family
## 2026-03-04 18:15:56.344673 INFO::Fitting model to feature number 87, Mycobacterium
## 2026-03-04 18:15:56.349871 INFO::Fitting model to feature number 88, Tundrisphaera
## 2026-03-04 18:15:56.35465 INFO::Fitting model to feature number 89, Arthrobacter
## 2026-03-04 18:15:56.359305 INFO::Fitting model to feature number 90, Peptostreptococcaceae.Family
## 2026-03-04 18:15:56.363994 INFO::Fitting model to feature number 91, Enterobacteriaceae.Family
## 2026-03-04 18:15:56.368846 INFO::Fitting model to feature number 92, Clostridium_sensu_stricto_2
## 2026-03-04 18:15:56.374976 INFO::Fitting model to feature number 93, Paludicola
## 2026-03-04 18:15:56.380325 INFO::Fitting model to feature number 94, Helicobacter
## 2026-03-04 18:15:56.385164 INFO::Fitting model to feature number 95, X.Eubacterium._ruminantium_group
## 2026-03-04 18:15:56.390062 INFO::Fitting model to feature number 96, Erysipelatoclostridiaceae.Family
## 2026-03-04 18:15:56.394781 INFO::Fitting model to feature number 97, Anaerofustis
## 2026-03-04 18:15:56.400194 INFO::Fitting model to feature number 98, Enterobacterales.Order
## 2026-03-04 18:15:56.405135 INFO::Fitting model to feature number 99, Bacillales.Order
## 2026-03-04 18:15:56.409983 INFO::Fitting model to feature number 100, Aquisphaera
## 2026-03-04 18:15:56.414959 INFO::Fitting model to feature number 101, Muribaculaceae
## 2026-03-04 18:15:56.421194 INFO::Fitting model to feature number 102, Prevotella
## 2026-03-04 18:15:56.429046 INFO::Fitting model to feature number 103, Prevotellaceae.Family
## 2026-03-04 18:15:56.43466 INFO::Fitting model to feature number 104, Bacteroidia.Class
## 2026-03-04 18:15:56.439468 INFO::Fitting model to feature number 105, Rikenellaceae.Family
## 2026-03-04 18:15:56.444451 INFO::Fitting model to feature number 106, Blautia
## 2026-03-04 18:15:56.449267 INFO::Fitting model to feature number 107, Butyricimonas
## 2026-03-04 18:15:56.454063 INFO::Fitting model to feature number 108, p.251.o5
## 2026-03-04 18:15:56.459663 INFO::Fitting model to feature number 109, Selenomonadaceae.Family
## 2026-03-04 18:15:56.464674 INFO::Fitting model to feature number 110, Treponema
## 2026-03-04 18:15:56.469565 INFO::Fitting model to feature number 111, Prevotellaceae_UCG.003
## 2026-03-04 18:15:56.474542 INFO::Fitting model to feature number 112, Pyramidobacter
## 2026-03-04 18:15:56.47955 INFO::Fitting model to feature number 113, Victivallaceae
## 2026-03-04 18:15:56.484614 INFO::Fitting model to feature number 114, Olivibacter
## 2026-03-04 18:15:56.490209 INFO::Fitting model to feature number 115, Muribaculaceae.Family
## 2026-03-04 18:15:56.496038 INFO::Fitting model to feature number 116, Gordonibacter
## 2026-03-04 18:15:56.502424 INFO::Fitting model to feature number 117, Rs.K70_termite_group
## 2026-03-04 18:15:56.507889 INFO::Fitting model to feature number 118, Desulfovibrionaceae.Family
## 2026-03-04 18:15:56.514308 INFO::Fitting model to feature number 119, Oscillibacter
## 2026-03-04 18:15:56.521077 INFO::Fitting model to feature number 120, Weissella
## 2026-03-04 18:15:56.52829 INFO::Fitting model to feature number 121, Mucispirillum
## 2026-03-04 18:15:56.534861 INFO::Fitting model to feature number 122, Anaerovorax
## 2026-03-04 18:15:56.540527 INFO::Fitting model to feature number 123, Prevotellaceae_UCG.004
## 2026-03-04 18:15:56.54586 INFO::Fitting model to feature number 124, Candidatus_Arthromitus
## 2026-03-04 18:15:56.55226 INFO::Fitting model to feature number 125, Spirochaetaceae.Family
## 2026-03-04 18:15:56.558293 INFO::Fitting model to feature number 126, Marinilabiliaceae.Family
## 2026-03-04 18:15:56.563931 INFO::Fitting model to feature number 127, horsej.a03
## 2026-03-04 18:15:56.569043 INFO::Fitting model to feature number 128, Asteroleplasma
## 2026-03-04 18:15:56.575416 INFO::Fitting model to feature number 129, Fretibacterium
## 2026-03-04 18:15:56.582314 INFO::Fitting model to feature number 130, Mitsuokella
## 2026-03-04 18:15:56.589602 INFO::Fitting model to feature number 131, Anaerofilum
## 2026-03-04 18:15:56.595851 INFO::Fitting model to feature number 132, Staphylococcus
## 2026-03-04 18:15:56.602343 INFO::Fitting model to feature number 133, Rickettsiales.Order
## 2026-03-04 18:15:56.607814 INFO::Fitting model to feature number 134, Aerococcus
## 2026-03-04 18:15:56.614068 INFO::Fitting model to feature number 135, Streptococcus
## 2026-03-04 18:15:56.620149 INFO::Fitting model to feature number 136, Gemella
## 2026-03-04 18:15:56.626521 INFO::Fitting model to feature number 137, SP3.e08
## 2026-03-04 18:15:56.632184 INFO::Fitting model to feature number 138, Pasteurellaceae.Family
## 2026-03-04 18:15:56.638852 INFO::Fitting model to feature number 139, X67.14
## 2026-03-04 18:15:56.644142 INFO::Fitting model to feature number 140, Ammoniphilus
## 2026-03-04 18:15:56.650551 INFO::Fitting model to feature number 141, Saccharopolyspora
## 2026-03-04 18:15:56.656045 INFO::Fitting model to feature number 142, Victivallis
## 2026-03-04 18:15:56.662404 INFO::Fitting model to feature number 143, WPS.2
## 2026-03-04 18:15:56.667794 INFO::Fitting model to feature number 144, Veillonellales.Selenomonadales.Order
## 2026-03-04 18:15:56.674045 INFO::Fitting model to feature number 145, Enteroscipio
## 2026-03-04 18:15:56.680842 INFO::Fitting model to feature number 146, Flavobacteriaceae.Family
## 2026-03-04 18:15:56.688204 INFO::Fitting model to feature number 147, Parabacteroides
## 2026-03-04 18:15:56.69528 INFO::Fitting model to feature number 148, UCG.002
## 2026-03-04 18:15:56.70399 INFO::Fitting model to feature number 149, Synergistaceae.Family
## 2026-03-04 18:15:56.711346 INFO::Fitting model to feature number 150, Burkholderiales.Order
## 2026-03-04 18:15:56.722247 INFO::Fitting model to feature number 151, X.Eubacterium._brachy_group
## 2026-03-04 18:15:56.735833 INFO::Fitting model to feature number 152, Parasutterella
## 2026-03-04 18:15:56.747191 INFO::Fitting model to feature number 153, Lachnospiraceae_UCG.001
## 2026-03-04 18:15:56.756958 INFO::Fitting model to feature number 154, Negativicutes.Class
## 2026-03-04 18:15:56.76594 INFO::Fitting model to feature number 155, Acetitomaculum
## 2026-03-04 18:15:56.773511 INFO::Fitting model to feature number 156, Peptococcus
## 2026-03-04 18:15:56.780886 INFO::Fitting model to feature number 157, ASF356
## 2026-03-04 18:15:56.788046 INFO::Fitting model to feature number 158, Tannerellaceae.Family
## 2026-03-04 18:15:56.797272 INFO::Fitting model to feature number 159, Odoribacter
## 2026-03-04 18:15:56.803772 INFO::Fitting model to feature number 160, Rikenella
## 2026-03-04 18:15:56.810285 INFO::Fitting model to feature number 161, Lachnospiraceae_UCG.002
## 2026-03-04 18:15:56.819187 INFO::Fitting model to feature number 162, Coprobacter
## 2026-03-04 18:15:56.825542 INFO::Fitting model to feature number 163, Lachnospiraceae_FCS020_group
## 2026-03-04 18:15:56.831938 INFO::Fitting model to feature number 164, Lachnospiraceae_UCG.010
## 2026-03-04 18:15:56.840516 INFO::Fitting model to feature number 165, DTU089
## 2026-03-04 18:15:56.847531 INFO::Fitting model to feature number 166, Mailhella
## 2026-03-04 18:15:56.853524 INFO::Fitting model to feature number 167, X.Ruminococcus._gnavus_group
## 2026-03-04 18:15:56.860225 INFO::Fitting model to feature number 168, Peptostreptococcales.Tissierellales.Order
## 2026-03-04 18:15:56.867229 INFO::Fitting model to feature number 169, Sphingomonas
## 2026-03-04 18:15:56.874169 INFO::Fitting model to feature number 170, Parvibacter
## 2026-03-04 18:15:56.881038 INFO::Fitting model to feature number 171, Rhizobiaceae.Family
## 2026-03-04 18:15:56.888681 INFO::Fitting model to feature number 172, Caulobacter
## 2026-03-04 18:15:56.895494 INFO::Fitting model to feature number 173, Roseateles
## 2026-03-04 18:15:56.902286 INFO::Fitting model to feature number 174, Cloacibacillus
## 2026-03-04 18:15:56.908969 INFO::Fitting model to feature number 175, Leuconostoc
## 2026-03-04 18:15:56.914418 INFO::Fitting model to feature number 176, Lactonifactor
## 2026-03-04 18:15:56.920576 INFO::Fitting model to feature number 177, Oligosphaeraceae.Family
## 2026-03-04 18:15:56.926508 INFO::Fitting model to feature number 178, Anaerotruncus
## 2026-03-04 18:15:56.93192 INFO::Fitting model to feature number 179, Hydrogenoanaerobacterium
## 2026-03-04 18:15:56.936918 INFO::Fitting model to feature number 180, Spirochaeta
## 2026-03-04 18:15:56.941942 INFO::Fitting model to feature number 181, Lachnospiraceae_UCG.004
## 2026-03-04 18:15:56.947241 INFO::Fitting model to feature number 182, Paracaedibacteraceae.Family
## 2026-03-04 18:15:56.952359 INFO::Fitting model to feature number 183, Atopobium
## 2026-03-04 18:15:56.957594 INFO::Fitting model to feature number 184, Coriobacteriia.Class
## 2026-03-04 18:15:56.962425 INFO::Fitting model to feature number 185, Succinatimonas
## 2026-03-04 18:15:56.967799 INFO::Fitting model to feature number 186, Glutamicibacter
## 2026-03-04 18:15:56.972959 INFO::Fitting model to feature number 187, Actinomyces
## 2026-03-04 18:15:56.978504 INFO::Fitting model to feature number 188, UCG.001
## 2026-03-04 18:15:56.983379 INFO::Fitting model to feature number 189, Planococcaceae.Family
## 2026-03-04 18:15:56.989401 INFO::Fitting model to feature number 190, Olsenella
## 2026-03-04 18:15:56.995169 INFO::Fitting model to feature number 191, Thermoanaerobacter
## 2026-03-04 18:15:57.001705 INFO::Fitting model to feature number 192, Caldicoprobacter
## 2026-03-04 18:15:57.008281 INFO::Fitting model to feature number 193, Ezakiella
## 2026-03-04 18:15:57.01509 INFO::Fitting model to feature number 194, Sporosarcina
## 2026-03-04 18:15:57.023376 INFO::Fitting model to feature number 195, Aeromicrobium
## 2026-03-04 18:15:57.031143 INFO::Fitting model to feature number 196, Solirubrobacteraceae.Family
## 2026-03-04 18:15:57.036542 INFO::Fitting model to feature number 197, Limosilactobacillus
## 2026-03-04 18:15:57.043126 INFO::Fitting model to feature number 198, Bilophila
## 2026-03-04 18:15:57.049703 INFO::Fitting model to feature number 199, Rothia
## 2026-03-04 18:15:57.056212 INFO::Fitting model to feature number 200, Corynebacterium
## 2026-03-04 18:15:57.062788 INFO::Fitting model to feature number 201, Ochrobactrum
## 2026-03-04 18:15:57.070352 INFO::Fitting model to feature number 202, TM7x
## 2026-03-04 18:15:57.078793 INFO::Fitting model to feature number 203, Angelakisella
## 2026-03-04 18:15:57.087339 INFO::Fitting model to feature number 204, Puniceicoccaceae.Family
## 2026-03-04 18:15:57.097332 INFO::Fitting model to feature number 205, Neisseriaceae.Family
## 2026-03-04 18:15:57.108235 INFO::Fitting model to feature number 206, Gaiellales.Order
## 2026-03-04 18:15:57.11935 INFO::Fitting model to feature number 207, Nocardioides
## 2026-03-04 18:15:57.129908 INFO::Fitting model to feature number 208, Micrococcaceae.Family
## 2026-03-04 18:15:57.137278 INFO::Fitting model to feature number 209, Micrococcus
## 2026-03-04 18:15:57.143162 INFO::Fitting model to feature number 210, Porphyromonas
## 2026-03-04 18:15:57.158226 INFO::Fitting model to feature number 211, Thermomonosporaceae.Family
## 2026-03-04 18:15:57.164554 INFO::Fitting model to feature number 212, Moraxella
## 2026-03-04 18:15:57.170399 INFO::Fitting model to feature number 213, Arcanobacterium
## 2026-03-04 18:15:57.175771 INFO::Fitting model to feature number 214, Lactobacillales.Order
## 2026-03-04 18:15:57.181348 INFO::Fitting model to feature number 215, Gemmataceae.Family
## 2026-03-04 18:15:57.187293 INFO::Fitting model to feature number 216, Brevibacillus
## 2026-03-04 18:15:57.194152 INFO::Fitting model to feature number 217, TK10
## 2026-03-04 18:15:57.200126 INFO::Fitting model to feature number 218, Rokubacteriales
## 2026-03-04 18:15:57.207587 INFO::Fitting model to feature number 219, Prevotellaceae_UCG.001
## 2026-03-04 18:15:57.214239 INFO::Fitting model to feature number 220, UCG.009
## 2026-03-04 18:15:57.221468 INFO::Fitting model to feature number 221, X.Eubacterium._nodatum_group
## 2026-03-04 18:15:57.228521 INFO::Fitting model to feature number 222, Microbacteriaceae.Family
## 2026-03-04 18:15:57.235541 INFO::Fitting model to feature number 223, Microtrichales.Order
## 2026-03-04 18:15:57.24377 INFO::Fitting model to feature number 224, Allobaculum
## 2026-03-04 18:15:57.251718 INFO::Fitting model to feature number 225, GCA.900066575
## 2026-03-04 18:15:57.260129 INFO::Fitting model to feature number 226, Actinobacteriota.Phylum
## 2026-03-04 18:15:57.265469 INFO::Fitting model to feature number 227, Holdemania
## 2026-03-04 18:15:57.270918 INFO::Fitting model to feature number 228, UCG.011
## 2026-03-04 18:15:57.276489 INFO::Fitting model to feature number 229, Raoultibacter
## 2026-03-04 18:15:57.282353 INFO::Fitting model to feature number 230, Syntrophomonadaceae.Family
## 2026-03-04 18:15:57.287469 INFO::Fitting model to feature number 231, Oxalobacteraceae.Family
## 2026-03-04 18:15:57.292944 INFO::Fitting model to feature number 232, Flavobacterium
## 2026-03-04 18:15:57.297861 INFO::Fitting model to feature number 233, Citrobacter
## 2026-03-04 18:15:57.3022 INFO::Fitting model to feature number 234, RS62_marine_group
## 2026-03-04 18:15:57.306458 INFO::Fitting model to feature number 235, Fusobacterium
## 2026-03-04 18:15:57.31092 INFO::Fitting model to feature number 236, CL500.29_marine_group
## 2026-03-04 18:15:57.315367 INFO::Fitting model to feature number 237, Anaerostipes
## 2026-03-04 18:15:57.319711 INFO::Fitting model to feature number 238, Coprococcus
## 2026-03-04 18:15:57.324119 INFO::Fitting model to feature number 239, Marvinbryantia
## 2026-03-04 18:15:57.329447 INFO::Fitting model to feature number 240, Dorea
## 2026-03-04 18:15:57.335224 INFO::Fitting model to feature number 241, Oscillospira
## 2026-03-04 18:15:57.340555 INFO::Fitting model to feature number 242, Lachnospiraceae_ND3007_group
## 2026-03-04 18:15:57.344771 INFO::Fitting model to feature number 243, Elusimicrobium
## 2026-03-04 18:15:57.34908 INFO::Fitting model to feature number 244, Erysipelotrichaceae_UCG.003
## 2026-03-04 18:15:57.353408 INFO::Fitting model to feature number 245, Adlercreutzia
## 2026-03-04 18:15:57.358755 INFO::Fitting model to feature number 246, Lachnospiraceae_UCG.006
## 2026-03-04 18:15:57.363725 INFO::Fitting model to feature number 247, Anaerosporobacter
## 2026-03-04 18:15:57.367949 INFO::Fitting model to feature number 248, Agathobacter
## 2026-03-04 18:15:57.372478 INFO::Fitting model to feature number 249, Prevotella_9
## 2026-03-04 18:15:57.376664 INFO::Fitting model to feature number 250, Shuttleworthia
## 2026-03-04 18:15:57.381162 INFO::Fitting model to feature number 251, Tannerellaceae
## 2026-03-04 18:15:57.385394 INFO::Fitting model to feature number 252, Megasphaera
## 2026-03-04 18:15:57.389843 INFO::Fitting model to feature number 253, Paraclostridium
## 2026-03-04 18:15:57.394101 INFO::Fitting model to feature number 254, Clostridium_sensu_stricto_11
## 2026-03-04 18:15:57.398469 INFO::Fitting model to feature number 255, Comamonadaceae.Family
## 2026-03-04 18:15:57.402654 INFO::Fitting model to feature number 256, Selenomonas
## 2026-03-04 18:15:57.407789 INFO::Fitting model to feature number 257, Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium
## 2026-03-04 18:15:57.41301 INFO::Fitting model to feature number 258, Brachybacterium
## 2026-03-04 18:15:57.418755 INFO::Fitting model to feature number 259, X11.24
## 2026-03-04 18:15:57.423103 INFO::Fitting model to feature number 260, Bordetella
## 2026-03-04 18:15:57.428261 INFO::Fitting model to feature number 261, Stenotrophomonas
## 2026-03-04 18:15:57.43323 INFO::Fitting model to feature number 262, Chryseobacterium
## 2026-03-04 18:15:57.438472 INFO::Fitting model to feature number 263, X.Eubacterium._hallii_group
## 2026-03-04 18:15:57.442955 INFO::Fitting model to feature number 264, Streptomyces
## 2026-03-04 18:15:57.449191 INFO::Fitting model to feature number 265, Rhodanobacter
## 2026-03-04 18:15:57.454228 INFO::Fitting model to feature number 266, Brevundimonas
## 2026-03-04 18:15:57.458975 INFO::Fitting model to feature number 267, Novosphingobium
## 2026-03-04 18:15:57.463849 INFO::Fitting model to feature number 268, Pedobacter
## 2026-03-04 18:15:57.468274 INFO::Fitting model to feature number 269, Marmoricola
## 2026-03-04 18:15:57.472607 INFO::Fitting model to feature number 270, Luteipulveratus
## 2026-03-04 18:15:57.477019 INFO::Fitting model to feature number 271, Pseudoxanthomonas
## 2026-03-04 18:15:57.482342 INFO::Fitting model to feature number 272, Luteimonas
## 2026-03-04 18:15:57.487556 INFO::Fitting model to feature number 273, Arachidicoccus
## 2026-03-04 18:15:57.49198 INFO::Fitting model to feature number 274, Saccharimonadales
## 2026-03-04 18:15:57.497143 INFO::Fitting model to feature number 275, Asticcacaulis
## 2026-03-04 18:15:57.503588 INFO::Fitting model to feature number 276, LWQ8
## 2026-03-04 18:15:57.50882 INFO::Fitting model to feature number 277, Ferruginibacter
## 2026-03-04 18:15:57.513177 INFO::Fitting model to feature number 278, Dysgonomonadaceae.Family
## 2026-03-04 18:15:57.518504 INFO::Fitting model to feature number 279, Turicibacter
## 2026-03-04 18:15:57.523516 INFO::Fitting model to feature number 280, vadinBE97
## 2026-03-04 18:15:57.52779 INFO::Fitting model to feature number 281, Sarcina
## 2026-03-04 18:15:57.53302 INFO::Fitting model to feature number 282, Sanguibacteroides
## 2026-03-04 18:15:57.53857 INFO::Fitting model to feature number 283, Clostridiales.Order
## 2026-03-04 18:15:57.543123 INFO::Fitting model to feature number 284, Micromonosporaceae.Family
## 2026-03-04 18:15:57.547299 INFO::Fitting model to feature number 285, Xanthobacteraceae.Family
## 2026-03-04 18:15:57.551643 INFO::Fitting model to feature number 286, Psychroglaciecola
## 2026-03-04 18:15:57.556822 INFO::Fitting model to feature number 287, Lactobacillus
## 2026-03-04 18:15:57.56204 INFO::Fitting model to feature number 288, Acetatifactor
## 2026-03-04 18:15:57.567368 INFO::Fitting model to feature number 289, Faecalibacterium
## 2026-03-04 18:15:57.572276 INFO::Fitting model to feature number 290, Cellulosilyticum
## 2026-03-04 18:15:57.577338 INFO::Fitting model to feature number 291, Mycoplasma
## 2026-03-04 18:15:57.581886 INFO::Fitting model to feature number 292, Isosphaeraceae.Family
## 2026-03-04 18:15:57.58614 INFO::Fitting model to feature number 293, Bosea
## 2026-03-04 18:15:57.5903 INFO::Fitting model to feature number 294, Kingella
## 2026-03-04 18:15:57.594641 INFO::Fitting model to feature number 295, Pedomicrobium
## 2026-03-04 18:15:57.599969 INFO::Fitting model to feature number 296, Lysobacter
## 2026-03-04 18:15:57.60483 INFO::Fitting model to feature number 297, Anaeroplasma
## 2026-03-04 18:15:57.609384 INFO::Fitting model to feature number 298, Candidatus_Stoquefichus
## 2026-03-04 18:15:57.613689 INFO::Fitting model to feature number 299, Campylobacter
## 2026-03-04 18:15:57.617896 INFO::Fitting model to feature number 300, Bacilli.Class
## 2026-03-04 18:15:57.62235 INFO::Fitting model to feature number 301, Acholeplasmataceae.Family
## 2026-03-04 18:15:57.626605 INFO::Fitting model to feature number 302, Cellulomonas
## 2026-03-04 18:15:57.63087 INFO::Fitting model to feature number 303, Pseudonocardia
## 2026-03-04 18:15:57.635238 INFO::Fitting model to feature number 304, Faecalibaculum
## 2026-03-04 18:15:57.639467 INFO::Fitting model to feature number 305, X.Eubacterium._xylanophilum_group
## 2026-03-04 18:15:57.644871 INFO::Fitting model to feature number 306, Lactococcus
## 2026-03-04 18:15:57.649919 INFO::Fitting model to feature number 307, Butyrivibrio
## 2026-03-04 18:15:57.655197 INFO::Fitting model to feature number 308, UBA1819
## 2026-03-04 18:15:57.659969 INFO::Fitting model to feature number 309, Ileibacterium
## 2026-03-04 18:15:57.665218 INFO::Fitting model to feature number 310, p.2534.18B5_gut_group
## 2026-03-04 18:15:57.670215 INFO::Fitting model to feature number 311, Erwiniaceae.Family
## 2026-03-04 18:15:57.675561 INFO::Fitting model to feature number 312, hgcI_clade
## 2026-03-04 18:15:57.679952 INFO::Fitting model to feature number 313, Rhodoferax
## 2026-03-04 18:15:57.684439 INFO::Fitting model to feature number 314, Pantoea
## 2026-03-04 18:15:57.688736 INFO::Fitting model to feature number 315, Polaromonas
## 2026-03-04 18:15:57.694129 INFO::Fitting model to feature number 316, Pseudorhodobacter
## 2026-03-04 18:15:57.699089 INFO::Fitting model to feature number 317, Christensenellaceae
## 2026-03-04 18:15:57.703969 INFO::Fitting model to feature number 318, Clade_III
## 2026-03-04 18:15:57.709276 INFO::Fitting model to feature number 319, Dongia
## 2026-03-04 18:15:57.714187 INFO::Fitting model to feature number 320, Hungatella
## 2026-03-04 18:15:57.718492 INFO::Fitting model to feature number 321, Epulopiscium
## 2026-03-04 18:15:57.723093 INFO::Fitting model to feature number 322, Ligilactobacillus
## 2026-03-04 18:15:57.727771 INFO::Fitting model to feature number 323, Dubosiella
## 2026-03-04 18:15:57.731891 INFO::Fitting model to feature number 324, Alloprevotella
## 2026-03-04 18:15:57.737782 INFO::Fitting model to feature number 325, Muribaculum
## 2026-03-04 18:15:57.743875 INFO::Fitting model to feature number 326, Romboutsia
## 2026-03-04 18:15:57.750111 INFO::Fitting model to feature number 327, Jeotgalicoccus
## 2026-03-04 18:15:57.754865 INFO::Fitting model to feature number 328, Negativibacillus
## 2026-03-04 18:15:57.760002 INFO::Fitting model to feature number 329, Erysipelatoclostridiaceae
## 2026-03-04 18:15:57.766183 INFO::Fitting model to feature number 330, Atopostipes
## 2026-03-04 18:15:57.77219 INFO::Fitting model to feature number 331, Frisingicoccus
## 2026-03-04 18:15:57.778362 INFO::Fitting model to feature number 332, Veillonella
## 2026-03-04 18:15:57.783889 INFO::Fitting model to feature number 333, Intestinimonas
## 2026-03-04 18:15:57.789292 INFO::Fitting model to feature number 334, Fournierella
## 2026-03-04 18:15:57.794477 INFO::Fitting model to feature number 335, Lachnoclostridium
## 2026-03-04 18:15:57.801936 INFO::Fitting model to feature number 336, Erysipelotrichaceae
## 2026-03-04 18:15:57.808095 INFO::Fitting model to feature number 337, Fibrobacter
## 2026-03-04 18:15:57.813927 INFO::Fitting model to feature number 338, Bifidobacterium
## 2026-03-04 18:15:57.81935 INFO::Fitting model to feature number 339, Coriobacteriaceae_UCG.002
## 2026-03-04 18:15:57.824143 INFO::Fitting model to feature number 340, Harryflintia
## 2026-03-04 18:15:57.829003 INFO::Fitting model to feature number 341, Bacillaceae.Family
## 2026-03-04 18:15:57.834478 INFO::Fitting model to feature number 342, Ralstonia
## 2026-03-04 18:15:57.840191 INFO::Fitting model to feature number 343, X.Ruminococcus._torques_group
## 2026-03-04 18:15:57.846074 INFO::Fitting model to feature number 344, X.Eubacterium._fissicatena_group
## 2026-03-04 18:15:57.851694 INFO::Fitting model to feature number 345, X.Clostridium._innocuum_group
## 2026-03-04 18:15:57.85696 INFO::Fitting model to feature number 346, Ureaplasma
## 2026-03-04 18:15:57.861757 INFO::Fitting model to feature number 347, Rummeliibacillus
## 2026-03-04 18:15:57.867148 INFO::Fitting model to feature number 348, Fructilactobacillus
## 2026-03-04 18:15:57.87307 INFO::Fitting model to feature number 349, Family_XI.Family
## 2026-03-04 18:15:57.879478 INFO::Fitting model to feature number 350, Enterobacter
## 2026-03-04 18:15:57.885713 INFO::Fitting model to feature number 351, Actinomycetaceae.Family
## 2026-03-04 18:15:57.891167 INFO::Fitting model to feature number 352, Pseudomonas
## 2026-03-04 18:15:57.897048 INFO::Fitting model to feature number 353, Verruc.01
## 2026-03-04 18:15:57.902444 INFO::Fitting model to feature number 354, Eremococcus
## 2026-03-04 18:15:57.907074 INFO::Fitting model to feature number 355, Sutterella
## 2026-03-04 18:15:57.912655 INFO::Fitting model to feature number 356, Lactiplantibacillus
## 2026-03-04 18:15:57.918741 INFO::Fitting model to feature number 357, Erysipelotrichales.Order
## 2026-03-04 18:15:57.923664 INFO::Fitting model to feature number 358, PeM15
## 2026-03-04 18:15:57.929168 INFO::Fitting model to feature number 359, Phascolarctobacterium
## 2026-03-04 18:15:57.934157 INFO::Fitting model to feature number 360, Paludibacteraceae.Family
## 2026-03-04 18:15:57.938737 INFO::Fitting model to feature number 361, Sphaerochaeta
## 2026-03-04 18:15:57.943835 INFO::Fitting model to feature number 362, p.1088.a5_gut_group
## 2026-03-04 18:15:57.94863 INFO::Fitting model to feature number 363, Bacteroidales_RF16_group
## 2026-03-04 18:15:57.954134 INFO::Fitting model to feature number 364, UCG.008
## 2026-03-04 18:15:57.959987 INFO::Fitting model to feature number 365, WCHB1.41
## 2026-03-04 18:15:57.965266 INFO::Fitting model to feature number 366, Oribacterium
## 2026-03-04 18:15:57.971019 INFO::Fitting model to feature number 367, Eisenbergiella
## 2026-03-04 18:15:57.975969 INFO::Fitting model to feature number 368, Mogibacterium
## 2026-03-04 18:15:57.981751 INFO::Fitting model to feature number 369, Succinivibrionaceae.Family
## 2026-03-04 18:15:57.988257 INFO::Fitting model to feature number 370, Anaerobiospirillum
## 2026-03-04 18:15:57.994975 INFO::Fitting model to feature number 371, Acetanaerobacterium
## 2026-03-04 18:15:58.001738 INFO::Fitting model to feature number 372, Caproiciproducens
## 2026-03-04 18:15:58.008001 INFO::Fitting model to feature number 373, GWE2.42.42
## 2026-03-04 18:15:58.01374 INFO::Fitting model to feature number 374, Neisseria
## 2026-03-04 18:15:58.019376 INFO::Fitting model to feature number 375, Proteus
## 2026-03-04 18:15:58.025126 INFO::Fitting model to feature number 376, Intrasporangiaceae.Family
## 2026-03-04 18:15:58.030923 INFO::Fitting model to feature number 377, Undibacterium
## 2026-03-04 18:15:58.036665 INFO::Fitting model to feature number 378, Sphingobacteriales.Order
## 2026-03-04 18:15:58.042503 INFO::Fitting model to feature number 379, Limnohabitans
## 2026-03-04 18:15:58.048159 INFO::Fitting model to feature number 380, Saccharimonadales.Order
## 2026-03-04 18:15:58.054343 INFO::Fitting model to feature number 381, Shewanella
## 2026-03-04 18:15:58.059398 INFO::Fitting model to feature number 382, Actinomycetaceae
## 2026-03-04 18:15:58.064123 INFO::Fitting model to feature number 383, Conexibacter
## 2026-03-04 18:15:58.069789 INFO::Fitting model to feature number 384, OM43_clade
## 2026-03-04 18:15:58.075786 INFO::Fitting model to feature number 385, Lachnospiraceae_AC2044_group
## 2026-03-04 18:15:58.080675 INFO::Fitting model to feature number 386, Halarsenatibacter
## 2026-03-04 18:15:58.086178 INFO::Fitting model to feature number 387, Nesterenkonia
## 2026-03-04 18:15:58.09222 INFO::Fitting model to feature number 388, Staphylococcaceae.Family
## 2026-03-04 18:15:58.099866 INFO::Fitting model to feature number 389, Actinobacteria.Class
## 2026-03-04 18:15:58.108091 INFO::Fitting model to feature number 390, Megamonas
## 2026-03-04 18:15:58.118482 INFO::Fitting model to feature number 391, Prevotella_7
## 2026-03-04 18:15:58.129025 INFO::Fitting model to feature number 392, Providencia
## 2026-03-04 18:15:58.141685 INFO::Fitting model to feature number 393, Facklamia
## 2026-03-04 18:15:58.153458 INFO::Fitting model to feature number 394, Solibacillus
## 2026-03-04 18:15:58.165862 INFO::Fitting model to feature number 395, Virgibacillus
## 2026-03-04 18:15:58.176684 INFO::Fitting model to feature number 396, Gracilibacillus
## 2026-03-04 18:15:58.186533 INFO::Fitting model to feature number 397, Desulfovibrionales.Order
## 2026-03-04 18:15:58.193976 INFO::Fitting model to feature number 398, Opitutales.Order
## 2026-03-04 18:15:58.201866 INFO::Fitting model to feature number 399, Marinilactibacillus
## 2026-03-04 18:15:58.209326 INFO::Fitting model to feature number 400, GCA.900066755
## 2026-03-04 18:15:58.21756 INFO::Fitting model to feature number 401, Christensenella
## 2026-03-04 18:15:58.224399 INFO::Fitting model to feature number 402, Cuneatibacter
## 2026-03-04 18:15:58.233417 INFO::Fitting model to feature number 403, Collinsella
## 2026-03-04 18:15:58.241522 INFO::Fitting model to feature number 404, Lachnospiraceae_UCG.008
## 2026-03-04 18:15:58.247591 INFO::Fitting model to feature number 405, Salinicoccus
## 2026-03-04 18:15:58.255754 INFO::Fitting model to feature number 406, Bergeyella
## 2026-03-04 18:15:58.261855 INFO::Fitting model to feature number 407, Aerococcaceae.Family
## 2026-03-04 18:15:58.26652 INFO::Fitting model to feature number 408, Catenibacillus
## 2026-03-04 18:15:58.272155 INFO::Fitting model to feature number 409, Catabacter
## 2026-03-04 18:15:58.278343 INFO::Fitting model to feature number 410, Butyricicoccaceae.Family
## 2026-03-04 18:15:58.285881 INFO::Fitting model to feature number 411, Bradyrhizobium
## 2026-03-04 18:15:58.295531 INFO::Fitting model to feature number 412, Caulobacteraceae.Family
## 2026-03-04 18:15:58.304166 INFO::Fitting model to feature number 413, Rhodococcus
## 2026-03-04 18:15:58.311328 INFO::Fitting model to feature number 414, Micrococcales.Order
## 2026-03-04 18:15:58.316979 INFO::Fitting model to feature number 415, Altererythrobacter
## 2026-03-04 18:15:58.322738 INFO::Fitting model to feature number 416, Variovorax
## 2026-03-04 18:15:58.329197 INFO::Fitting model to feature number 417, Ilumatobacteraceae.Family
## 2026-03-04 18:15:58.33698 INFO::Fitting model to feature number 418, YC.ZSS.LKJ147
## 2026-03-04 18:15:58.34505 INFO::Fitting model to feature number 419, Terrimonas
## 2026-03-04 18:15:58.351002 INFO::Fitting model to feature number 420, WD2101_soil_group
## 2026-03-04 18:15:58.35643 INFO::Fitting model to feature number 421, Xanthomonadaceae.Family
## 2026-03-04 18:15:58.361619 INFO::Fitting model to feature number 422, Actinomycetospora
## 2026-03-04 18:15:58.366651 INFO::Fitting model to feature number 423, Luteibacter
## 2026-03-04 18:15:58.372292 INFO::Fitting model to feature number 424, IMCC26256
## 2026-03-04 18:15:58.378752 INFO::Fitting model to feature number 425, Micropepsis
## 2026-03-04 18:15:58.384624 INFO::Fitting model to feature number 426, Chryseolinea
## 2026-03-04 18:15:58.391207 INFO::Fitting model to feature number 427, BIyi10
## 2026-03-04 18:15:58.397386 INFO::Fitting model to feature number 428, Rhodopirellula
## 2026-03-04 18:15:58.403527 INFO::Fitting model to feature number 429, Lachnospira
## 2026-03-04 18:15:58.409212 INFO::Fitting model to feature number 430, Comamonas
## 2026-03-04 18:15:58.414442 INFO::Fitting model to feature number 431, Absconditabacteriales_.SR1.
## 2026-03-04 18:15:58.420472 INFO::Fitting model to feature number 432, Gracilibacteria
## 2026-03-04 18:15:58.42592 INFO::Fitting model to feature number 433, Fibrobacteraceae.Family
## 2026-03-04 18:15:58.431485 INFO::Fitting model to feature number 434, Brevibacterium
## 2026-03-04 18:15:58.437781 INFO::Fitting model to feature number 435, U29.B03
## 2026-03-04 18:15:58.445521 INFO::Fitting model to feature number 436, Microbacterium
## 2026-03-04 18:15:58.453723 INFO::Fitting model to feature number 437, Dielma
## 2026-03-04 18:15:58.460814 INFO::Fitting model to feature number 438, TM7a
## 2026-03-04 18:15:58.468391 INFO::Fitting model to feature number 439, Cytophaga
## 2026-03-04 18:15:58.47562 INFO::Fitting model to feature number 440, Coxiella
## 2026-03-04 18:15:58.486386 INFO::Fitting model to feature number 441, Natronincola
## 2026-03-04 18:15:58.497083 INFO::Fitting model to feature number 442, Sphingobacterium
## 2026-03-04 18:15:58.510172 INFO::Fitting model to feature number 443, Delftia
## 2026-03-04 18:15:58.523969 INFO::Fitting model to feature number 444, Achromobacter
## 2026-03-04 18:15:58.536416 INFO::Fitting model to feature number 445, Taibaiella
## 2026-03-04 18:15:58.547918 INFO::Fitting model to feature number 446, Massilia
## 2026-03-04 18:15:58.557499 INFO::Fitting model to feature number 447, Sphingobium
## 2026-03-04 18:15:58.564367 INFO::Fitting model to feature number 448, Dyadobacter
## 2026-03-04 18:15:58.571103 INFO::Fitting model to feature number 449, Microscillaceae.Family
## 2026-03-04 18:15:58.576228 INFO::Fitting model to feature number 450, Actinoallomurus
## 2026-03-04 18:15:58.582721 INFO::Fitting model to feature number 451, Cohnella
## 2026-03-04 18:15:58.587861 INFO::Fitting model to feature number 452, Roseomonas
## 2026-03-04 18:15:58.594452 INFO::Fitting model to feature number 453, Pseudosphingobacterium
## 2026-03-04 18:15:58.687385 INFO::Counting total values for each feature
## 2026-03-04 18:15:58.760738 INFO::Writing filtered data to file ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/features/filtered_data.tsv
## 2026-03-04 18:15:58.837134 INFO::Writing filtered, normalized data to file ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/features/filtered_data_norm.tsv
## 2026-03-04 18:15:58.943316 INFO::Writing filtered, normalized, transformed data to file ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/features/filtered_data_norm_transformed.tsv
## 2026-03-04 18:15:59.188275 INFO::Writing residuals to file ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/fits/residuals.rds
## 2026-03-04 18:15:59.218236 INFO::Writing fitted values to file ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/fits/fitted.rds
## 2026-03-04 18:15:59.235613 INFO::Writing all results to file (ordered by increasing q-values): ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/all_results.tsv
## 2026-03-04 18:15:59.381414 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.050000 ) to file (ordered by increasing q-values): ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/significant_results.tsv
## 2026-03-04 18:15:59.428238 INFO::Writing heatmap of significant results to file: ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR/heatmap.pdf
## 2026-03-04 18:15:59.568805 INFO::Writing association plots (one for each significant association) to output folder: ./output/maaslin2/pooled-output/rare/20260304_18_15_55-NMR-DMR-B6mouse-MSMmouse-FVBNmouse-spalax-pvo-hare-rabbit-Genus-host-234-ref-NMR
## 2026-03-04 18:15:59.601561 INFO::Plotting associations from most to least significant, grouped by metadata
## 2026-03-04 18:15:59.606066 INFO::Plotting data for metadata number 1, class
## 2026-03-04 18:15:59.611961 INFO::Creating boxplot for categorical data, class vs Veillonellales.Selenomonadales.Order
```

```
## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
## ℹ Please use the `linewidth` argument instead.
## ℹ The deprecated feature was likely used in the Maaslin2 package.
##   Please report the issue at <https://github.com/biobakery/maaslin2/issues>.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```
## 2026-03-04 18:15:59.892918 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:16:00.09212 INFO::Creating boxplot for categorical data, class vs Flavobacteriaceae.Family
```

```
## 2026-03-04 18:16:00.263879 INFO::Creating boxplot for categorical data, class vs Muribaculum
```

```
## 2026-03-04 18:16:00.442944 INFO::Creating boxplot for categorical data, class vs Victivallis
```

```
## 2026-03-04 18:16:00.642483 INFO::Creating boxplot for categorical data, class vs dgA.11_gut_group
```

```
## 2026-03-04 18:16:00.825046 INFO::Creating boxplot for categorical data, class vs UCG.002
```

```
## 2026-03-04 18:16:01.003256 INFO::Creating boxplot for categorical data, class vs Negativicutes.Class
```

```
## 2026-03-04 18:16:01.192835 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:01.381962 INFO::Creating boxplot for categorical data, class vs Barnesiellaceae.Family
```

```
## 2026-03-04 18:16:01.563697 INFO::Creating boxplot for categorical data, class vs Helicobacter
```

```
## 2026-03-04 18:16:01.779042 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:01.966828 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:02.150543 INFO::Creating boxplot for categorical data, class vs V9D2013_group
```

```
## 2026-03-04 18:16:02.345427 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:02.525773 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:02.704247 INFO::Creating boxplot for categorical data, class vs Muribaculum
```

```
## 2026-03-04 18:16:02.898654 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:03.084917 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:03.267833 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:03.460199 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:03.648726 INFO::Creating boxplot for categorical data, class vs Subdoligranulum
```

```
## 2026-03-04 18:16:03.828399 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4A136_group
```

```
## 2026-03-04 18:16:04.027503 INFO::Creating boxplot for categorical data, class vs Eggerthella
```

```
## 2026-03-04 18:16:04.218608 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:04.39842 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:04.592559 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:04.776853 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:04.963887 INFO::Creating boxplot for categorical data, class vs WPS.2
```

```
## 2026-03-04 18:16:05.157642 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:16:05.342632 INFO::Creating boxplot for categorical data, class vs Gordonibacter
```

```
## 2026-03-04 18:16:05.531826 INFO::Creating boxplot for categorical data, class vs Rikenellaceae.Family
```

```
## 2026-03-04 18:16:05.720181 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:05.908276 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:06.088327 INFO::Creating boxplot for categorical data, class vs Burkholderiales.Order
```

```
## 2026-03-04 18:16:06.279864 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:16:06.465911 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:16:06.655315 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:06.84536 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:07.028763 INFO::Creating boxplot for categorical data, class vs Enteroscipio
```

```
## 2026-03-04 18:16:07.214132 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae
```

```
## 2026-03-04 18:16:07.401734 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:16:07.582754 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:07.770613 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:07.961736 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:08.141832 INFO::Creating boxplot for categorical data, class vs dgA.11_gut_group
```

```
## 2026-03-04 18:16:08.328169 INFO::Creating boxplot for categorical data, class vs Parabacteroides
```

```
## 2026-03-04 18:16:08.519814 INFO::Creating boxplot for categorical data, class vs Butyricimonas
```

```
## 2026-03-04 18:16:08.706923 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:16:08.889557 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:09.116047 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:09.304514 INFO::Creating boxplot for categorical data, class vs Alloprevotella
```

```
## 2026-03-04 18:16:09.498295 INFO::Creating boxplot for categorical data, class vs Parabacteroides
```

```
## 2026-03-04 18:16:09.694966 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:09.878517 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:10.057724 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:10.251074 INFO::Creating boxplot for categorical data, class vs Muribaculum
```

```
## 2026-03-04 18:16:10.437055 INFO::Creating boxplot for categorical data, class vs Helicobacter
```

```
## 2026-03-04 18:16:10.619737 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:16:10.817182 INFO::Creating boxplot for categorical data, class vs Butyricicoccus
```

```
## 2026-03-04 18:16:11.004205 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:16:11.186286 INFO::Creating boxplot for categorical data, class vs Prevotellaceae.Family
```

```
## 2026-03-04 18:16:11.384295 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:16:11.576652 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:16:11.762059 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:11.958836 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:16:12.149177 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:12.335804 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:12.5376 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:12.740426 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:12.929394 INFO::Creating boxplot for categorical data, class vs Bifidobacterium
```

```
## 2026-03-04 18:16:13.125926 INFO::Creating boxplot for categorical data, class vs Butyricicoccus
```

```
## 2026-03-04 18:16:13.32008 INFO::Creating boxplot for categorical data, class vs Barnesiellaceae.Family
```

```
## 2026-03-04 18:16:13.503915 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4A136_group
```

```
## 2026-03-04 18:16:13.70645 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:16:13.90388 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:16:14.091829 INFO::Creating boxplot for categorical data, class vs Olivibacter
```

```
## 2026-03-04 18:16:14.71804 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:16:14.900297 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:15.086207 INFO::Creating boxplot for categorical data, class vs Alloprevotella
```

```
## 2026-03-04 18:16:15.286568 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:16:15.476175 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:15.658192 INFO::Creating boxplot for categorical data, class vs RF39
```

```
## 2026-03-04 18:16:15.852638 INFO::Creating boxplot for categorical data, class vs Eggerthellaceae.Family
```

```
## 2026-03-04 18:16:16.037581 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:16.221845 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:16:16.419729 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:16.60312 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:16.79453 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:16:16.991337 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:16:17.179142 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:17.368318 INFO::Creating boxplot for categorical data, class vs Peptococcus
```

```
## 2026-03-04 18:16:17.566453 INFO::Creating boxplot for categorical data, class vs Tyzzerella
```

```
## 2026-03-04 18:16:17.753831 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:16:17.947523 INFO::Creating boxplot for categorical data, class vs Turicibacter
```

```
## 2026-03-04 18:16:18.149184 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:16:18.337519 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:16:18.52595 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:16:18.719158 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._coprostanoligenes_group
```

```
## 2026-03-04 18:16:18.911114 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4A136_group
```

```
## 2026-03-04 18:16:19.097753 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:16:19.287254 INFO::Creating boxplot for categorical data, class vs UCG.007
```

```
## 2026-03-04 18:16:19.466192 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:16:19.645201 INFO::Creating boxplot for categorical data, class vs Bacteroidales.Order
```

```
## 2026-03-04 18:16:19.833594 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:20.015485 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:20.199142 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:16:20.391343 INFO::Creating boxplot for categorical data, class vs Butyricimonas
```

```
## 2026-03-04 18:16:20.570415 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:16:20.75573 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:16:20.944019 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:16:21.123132 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:16:21.316727 INFO::Creating boxplot for categorical data, class vs Sporobacter
```

```
## 2026-03-04 18:16:21.499218 INFO::Creating boxplot for categorical data, class vs Coprobacter
```

```
## 2026-03-04 18:16:21.683528 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:21.876018 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:22.058045 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:16:22.23668 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.010
```

```
## 2026-03-04 18:16:22.434743 INFO::Creating boxplot for categorical data, class vs Clostridia_vadinBB60_group
```

```
## 2026-03-04 18:16:22.619155 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:16:22.80311 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:16:22.998885 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:23.181814 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:16:23.369107 INFO::Creating boxplot for categorical data, class vs Parabacteroides
```

```
## 2026-03-04 18:16:23.565454 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:16:23.747249 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.001
```

```
## 2026-03-04 18:16:23.926031 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:24.116932 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:24.292775 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:24.478056 INFO::Creating boxplot for categorical data, class vs Parabacteroides
```

```
## 2026-03-04 18:16:24.665916 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.003
```

```
## 2026-03-04 18:16:24.844207 INFO::Creating boxplot for categorical data, class vs Prevotellaceae.Family
```

```
## 2026-03-04 18:16:25.029402 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:16:25.218026 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:16:25.399338 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:16:25.587372 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:25.76686 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:16:25.947233 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:16:26.19549 INFO::Creating boxplot for categorical data, class vs Bacteroidales.Order
```

```
## 2026-03-04 18:16:26.39141 INFO::Creating boxplot for categorical data, class vs Rhodospirillales.Order
```

```
## 2026-03-04 18:16:26.59858 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:16:26.809314 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:16:27.001746 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:16:27.186261 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:16:27.382053 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:16:27.574032 INFO::Creating boxplot for categorical data, class vs Negativibacillus
```

```
## 2026-03-04 18:16:27.754647 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:16:27.94949 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:16:28.13627 INFO::Creating boxplot for categorical data, class vs Butyricicoccus
```

```
## 2026-03-04 18:16:28.317078 INFO::Creating boxplot for categorical data, class vs Prevotellaceae.Family
```

```
## 2026-03-04 18:16:28.511385 INFO::Creating boxplot for categorical data, class vs Alloprevotella
```

```
## 2026-03-04 18:16:28.699788 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:16:28.88151 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:16:29.076488 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:29.263669 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:16:29.445763 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4A136_group
```

```
## 2026-03-04 18:16:29.639817 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:16:29.827883 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:16:30.010258 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:16:30.202736 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.001
```

```
## 2026-03-04 18:16:30.391216 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:16:30.574552 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:30.765388 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:30.955469 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:16:31.156298 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:16:31.35937 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridium
```

```
## 2026-03-04 18:16:31.549664 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4B4_group
```

```
## 2026-03-04 18:16:31.736315 INFO::Creating boxplot for categorical data, class vs Eggerthellaceae.Family
```

```
## 2026-03-04 18:16:31.942273 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:16:32.126474 INFO::Creating boxplot for categorical data, class vs Bacteroides
```

```
## 2026-03-04 18:16:32.309329 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:16:32.507866 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:16:32.691397 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:16:32.872173 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:16:33.077271 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:16:33.259416 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:33.446638 INFO::Creating boxplot for categorical data, class vs Rikenella
```

```
## 2026-03-04 18:16:33.646129 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:16:33.835204 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:16:34.016502 INFO::Creating boxplot for categorical data, class vs Butyricicoccus
```

```
## 2026-03-04 18:16:34.214531 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:16:34.397178 INFO::Creating boxplot for categorical data, class vs Prevotellaceae.Family
```

```
## 2026-03-04 18:16:34.593895 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:34.783879 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:16:34.966254 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:16:35.159996 INFO::Creating boxplot for categorical data, class vs Bifidobacterium
```

```
## 2026-03-04 18:16:35.34376 INFO::Creating boxplot for categorical data, class vs X.Ruminococcus._gnavus_group
```

```
## 2026-03-04 18:16:35.528371 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:16:35.726327 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4A136_group
```

```
## 2026-03-04 18:16:35.914143 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:16:36.100892 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:16:36.299679 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:16:36.482425 INFO::Creating boxplot for categorical data, class vs Mitsuokella
```

```
## 2026-03-04 18:16:36.665267 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:16:36.861327 INFO::Creating boxplot for categorical data, class vs Clostridia_vadinBB60_group
```

```
## 2026-03-04 18:16:37.048111 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:16:37.231088 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:16:37.43101 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:16:37.616971 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:16:37.84627 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:16:38.041735 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:16:38.229827 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:16:38.421681 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._siraeum_group
```

```
## 2026-03-04 18:16:38.6193 INFO::Creating boxplot for categorical data, class vs Akkermansia
```

```
## 2026-03-04 18:16:38.814606 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:16:39.009451 INFO::Creating boxplot for categorical data, class vs Fibrobacter
```

```
## 2026-03-04 18:16:39.205869 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:16:39.393586 INFO::Creating boxplot for categorical data, class vs Rikenella
```

```
## 2026-03-04 18:16:39.58656 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:16:39.784117 INFO::Creating boxplot for categorical data, class vs Tyzzerella
```

```
## 2026-03-04 18:16:39.977189 INFO::Creating boxplot for categorical data, class vs RF39
```

```
## 2026-03-04 18:16:40.171094 INFO::Creating boxplot for categorical data, class vs Eggerthellaceae.Family
```

```
## 2026-03-04 18:16:40.369589 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:16:40.558354 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:16:40.748258 INFO::Creating boxplot for categorical data, class vs Adlercreutzia
```

```
## 2026-03-04 18:16:40.948845 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:16:41.137441 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ventriosum_group
```

```
## 2026-03-04 18:16:41.333936 INFO::Creating boxplot for categorical data, class vs RF39
```

```
## 2026-03-04 18:16:41.534414 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:16:41.727876 INFO::Creating boxplot for categorical data, class vs Defluviitaleaceae_UCG.011
```

```
## 2026-03-04 18:16:41.923438 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:16:42.117539 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:16:42.306309 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:16:42.508216 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:16:42.707155 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:42.899354 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:16:43.11483 INFO::Creating boxplot for categorical data, class vs Marinilabiliaceae.Family
```

```
## 2026-03-04 18:16:43.309197 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:16:43.503529 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:16:43.715363 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridiaceae.Family
```

```
## 2026-03-04 18:16:43.908087 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:16:44.098463 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:16:44.30198 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:16:44.487621 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:16:44.672588 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:16:44.877597 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:16:45.06649 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:16:45.251171 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:16:45.462769 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:16:45.652075 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:16:45.838741 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:16:46.046086 INFO::Creating boxplot for categorical data, class vs UCG.009
```

```
## 2026-03-04 18:16:46.23473 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:16:46.425226 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:16:46.636786 INFO::Creating boxplot for categorical data, class vs Tuzzerella
```

```
## 2026-03-04 18:16:46.834334 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:16:47.038226 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:16:47.243902 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:16:47.434796 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae
```

```
## 2026-03-04 18:16:47.644409 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.002
```

```
## 2026-03-04 18:16:47.839166 INFO::Creating boxplot for categorical data, class vs Muribaculaceae.Family
```

```
## 2026-03-04 18:16:48.025842 INFO::Creating boxplot for categorical data, class vs Gemella
```

```
## 2026-03-04 18:16:48.232987 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:16:48.420152 INFO::Creating boxplot for categorical data, class vs Subdoligranulum
```

```
## 2026-03-04 18:16:48.604049 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:16:48.806399 INFO::Creating boxplot for categorical data, class vs Ruminiclostridium
```

```
## 2026-03-04 18:16:48.993491 INFO::Creating boxplot for categorical data, class vs Family_XIII_AD3011_group
```

```
## 2026-03-04 18:16:49.178172 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:16:49.382606 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:16:49.573716 INFO::Creating boxplot for categorical data, class vs Oxalobacter
```

```
## 2026-03-04 18:16:49.817841 INFO::Creating boxplot for categorical data, class vs Bacteroidales.Order
```

```
## 2026-03-04 18:16:50.015358 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:16:50.210976 INFO::Creating boxplot for categorical data, class vs Oxalobacter
```

```
## 2026-03-04 18:16:50.412005 INFO::Creating boxplot for categorical data, class vs Turicibacter
```

```
## 2026-03-04 18:16:50.620879 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:16:50.814601 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:16:50.999838 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:16:51.210054 INFO::Creating boxplot for categorical data, class vs Rikenella
```

```
## 2026-03-04 18:16:51.4065 INFO::Creating boxplot for categorical data, class vs Marvinbryantia
```

```
## 2026-03-04 18:16:51.592511 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ventriosum_group
```

```
## 2026-03-04 18:16:51.799334 INFO::Creating boxplot for categorical data, class vs Coprococcus
```

```
## 2026-03-04 18:16:51.992506 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:16:52.18433 INFO::Creating boxplot for categorical data, class vs Bacteria.Kingdom
```

```
## 2026-03-04 18:16:52.393758 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae_UCG.003
```

```
## 2026-03-04 18:16:52.586418 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:16:52.771253 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:16:52.978102 INFO::Creating boxplot for categorical data, class vs Oscillospirales.Order
```

```
## 2026-03-04 18:16:53.170002 INFO::Creating boxplot for categorical data, class vs Prevotellaceae.Family
```

```
## 2026-03-04 18:16:53.3651 INFO::Creating boxplot for categorical data, class vs Victivallaceae
```

```
## 2026-03-04 18:16:53.58217 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:16:53.771339 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:16:53.973675 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:16:54.170441 INFO::Creating boxplot for categorical data, class vs Colidextribacter
```

```
## 2026-03-04 18:16:54.360828 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._hallii_group
```

```
## 2026-03-04 18:16:54.563259 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:16:54.750489 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:16:54.941427 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:16:55.149529 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:16:55.340603 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:16:55.537184 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:16:55.741395 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:16:55.936038 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:16:56.126448 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.001
```

```
## 2026-03-04 18:16:56.333786 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:16:56.528292 INFO::Creating boxplot for categorical data, class vs Odoribacter
```

```
## 2026-03-04 18:16:56.718933 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:16:56.92236 INFO::Creating boxplot for categorical data, class vs Oscillospira
```

```
## 2026-03-04 18:16:57.111921 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:16:57.301397 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:16:57.500527 INFO::Creating boxplot for categorical data, class vs Dubosiella
```

```
## 2026-03-04 18:16:57.692501 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:16:57.887287 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:16:58.090934 INFO::Creating boxplot for categorical data, class vs Natranaerobiaceae.Family
```

```
## 2026-03-04 18:16:58.28811 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:16:58.480739 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:16:58.686678 INFO::Creating boxplot for categorical data, class vs Desulfovibrionaceae.Family
```

```
## 2026-03-04 18:16:58.880808 INFO::Creating boxplot for categorical data, class vs Desulfovibrio
```

```
## 2026-03-04 18:16:59.07697 INFO::Creating boxplot for categorical data, class vs Coriobacteriales_Incertae_Sedis.Family
```

```
## 2026-03-04 18:16:59.284251 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:16:59.477069 INFO::Creating boxplot for categorical data, class vs Anaerovoracaceae.Family
```

```
## 2026-03-04 18:16:59.675481 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:16:59.870681 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:17:00.072145 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:17:00.278342 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:17:00.489285 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:17:00.686687 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:17:00.886167 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:17:01.081004 INFO::Creating boxplot for categorical data, class vs UCG.010
```

```
## 2026-03-04 18:17:01.269411 INFO::Creating boxplot for categorical data, class vs Mogibacterium
```

```
## 2026-03-04 18:17:01.474725 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:17:01.673549 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:17:01.85905 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:17:02.104459 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:17:02.294133 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:17:02.481833 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:17:02.693743 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:17:02.887648 INFO::Creating boxplot for categorical data, class vs Desulfovibrio
```

```
## 2026-03-04 18:17:03.08613 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:17:03.301372 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:17:03.500435 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:17:03.689499 INFO::Creating boxplot for categorical data, class vs Lachnoclostridium
```

```
## 2026-03-04 18:17:03.908285 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:17:04.099625 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:17:04.287586 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:17:04.501119 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:17:04.690727 INFO::Creating boxplot for categorical data, class vs Succinatimonas
```

```
## 2026-03-04 18:17:04.884076 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:17:05.093796 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:17:05.287062 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:17:05.482615 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:17:05.697484 INFO::Creating boxplot for categorical data, class vs Negativibacillus
```

```
## 2026-03-04 18:17:05.896086 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._nodatum_group
```

```
## 2026-03-04 18:17:06.089513 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.001
```

```
## 2026-03-04 18:17:06.302891 INFO::Creating boxplot for categorical data, class vs Methylobacterium.Methylorubrum
```

```
## 2026-03-04 18:17:06.497239 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:17:06.685809 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:17:06.907055 INFO::Creating boxplot for categorical data, class vs Rhodospirillales.Order
```

```
## 2026-03-04 18:17:07.101818 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:17:07.30412 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:17:07.505971 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:17:07.70503 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:17:07.917257 INFO::Creating boxplot for categorical data, class vs RF39
```

```
## 2026-03-04 18:17:08.111751 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:17:08.302949 INFO::Creating boxplot for categorical data, class vs Oscillospirales.Order
```

```
## 2026-03-04 18:17:08.518569 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:17:08.708613 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae
```

```
## 2026-03-04 18:17:08.898861 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:17:09.113515 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:17:09.302759 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:17:09.492053 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:17:09.699061 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:17:09.892152 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:17:10.08788 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:17:10.302673 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:17:10.500581 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:17:10.690419 INFO::Creating boxplot for categorical data, class vs Streptococcus
```

```
## 2026-03-04 18:17:10.909471 INFO::Creating boxplot for categorical data, class vs Anaerosporobacter
```

```
## 2026-03-04 18:17:11.097847 INFO::Creating boxplot for categorical data, class vs horsej.a03
```

```
## 2026-03-04 18:17:11.286386 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:17:11.505026 INFO::Creating boxplot for categorical data, class vs Rhodospirillales.Order
```

```
## 2026-03-04 18:17:11.703431 INFO::Creating boxplot for categorical data, class vs Odoribacter
```

```
## 2026-03-04 18:17:11.9143 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:17:12.11478 INFO::Creating boxplot for categorical data, class vs UCG.009
```

```
## 2026-03-04 18:17:12.310841 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:17:12.525144 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:17:12.718423 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:17:12.908011 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._xylanophilum_group
```

```
## 2026-03-04 18:17:13.118557 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:17:13.311999 INFO::Creating boxplot for categorical data, class vs Bacteroidales.Order
```

```
## 2026-03-04 18:17:13.504425 INFO::Creating boxplot for categorical data, class vs Odoribacter
```

```
## 2026-03-04 18:17:13.724549 INFO::Creating boxplot for categorical data, class vs Selenomonas
```

```
## 2026-03-04 18:17:13.916163 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:17:14.104016 INFO::Creating boxplot for categorical data, class vs Intestinimonas
```

```
## 2026-03-04 18:17:14.369953 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:17:14.571415 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:17:14.77211 INFO::Creating boxplot for categorical data, class vs Defluviitaleaceae_UCG.011
```

```
## 2026-03-04 18:17:14.985292 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:17:15.181574 INFO::Creating boxplot for categorical data, class vs Atopostipes
```

```
## 2026-03-04 18:17:15.374811 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:17:15.589246 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:17:15.783441 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:17:15.976878 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:17:16.191371 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:17:16.382949 INFO::Creating boxplot for categorical data, class vs Blautia
```

```
## 2026-03-04 18:17:16.574939 INFO::Creating boxplot for categorical data, class vs Papillibacter
```

```
## 2026-03-04 18:17:16.789204 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:17:16.983794 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:17:17.174618 INFO::Creating boxplot for categorical data, class vs Parvibacter
```

```
## 2026-03-04 18:17:17.390921 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:17:17.58291 INFO::Creating boxplot for categorical data, class vs Gordonibacter
```

```
## 2026-03-04 18:17:17.772055 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:17:17.99005 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae
```

```
## 2026-03-04 18:17:18.186657 INFO::Creating boxplot for categorical data, class vs Incertae_Sedis
```

```
## 2026-03-04 18:17:18.393036 INFO::Creating boxplot for categorical data, class vs Eubacteriales.Order
```

```
## 2026-03-04 18:17:18.594889 INFO::Creating boxplot for categorical data, class vs Clostridia_vadinBB60_group
```

```
## 2026-03-04 18:17:18.7881 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:17:18.997384 INFO::Creating boxplot for categorical data, class vs Asteroleplasma
```

```
## 2026-03-04 18:17:19.198529 INFO::Creating boxplot for categorical data, class vs Clostridia_vadinBB60_group
```

```
## 2026-03-04 18:17:19.394451 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:17:19.604454 INFO::Creating boxplot for categorical data, class vs Bacteria.Kingdom
```

```
## 2026-03-04 18:17:19.80085 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:17:19.990887 INFO::Creating boxplot for categorical data, class vs Tuzzerella
```

```
## 2026-03-04 18:17:20.195729 INFO::Creating boxplot for categorical data, class vs Rothia
```

```
## 2026-03-04 18:17:20.388613 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:17:20.579223 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:17:20.787158 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:17:20.988092 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:17:21.179772 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridiaceae
```

```
## 2026-03-04 18:17:21.386695 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:17:21.579457 INFO::Creating boxplot for categorical data, class vs ASF356
```

```
## 2026-03-04 18:17:21.769404 INFO::Creating boxplot for categorical data, class vs Fretibacterium
```

```
## 2026-03-04 18:17:21.982808 INFO::Creating boxplot for categorical data, class vs Candidatus_Stoquefichus
```

```
## 2026-03-04 18:17:22.176808 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:17:22.369495 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:17:22.580984 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:17:22.780322 INFO::Creating boxplot for categorical data, class vs V9D2013_group
```

```
## 2026-03-04 18:17:22.97465 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:17:23.186524 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ruminantium_group
```

```
## 2026-03-04 18:17:23.384163 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:17:23.571824 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:17:23.785054 INFO::Creating boxplot for categorical data, class vs Marvinbryantia
```

```
## 2026-03-04 18:17:23.98091 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:17:24.188865 INFO::Creating boxplot for categorical data, class vs Succinivibrionaceae.Family
```

```
## 2026-03-04 18:17:24.383934 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:17:24.585819 INFO::Creating boxplot for categorical data, class vs Faecalibaculum
```

```
## 2026-03-04 18:17:24.799905 INFO::Creating boxplot for categorical data, class vs Oscillospirales.Order
```

```
## 2026-03-04 18:17:25.001599 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:17:25.19155 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:17:25.408122 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:17:25.604164 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:17:25.792775 INFO::Creating boxplot for categorical data, class vs Sphaerochaeta
```

```
## 2026-03-04 18:17:26.000409 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:17:26.192294 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:17:26.849577 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:17:27.020803 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:17:27.194864 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._brachy_group
```

```
## 2026-03-04 18:17:27.376253 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:17:27.563594 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:17:27.748922 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:17:27.92776 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:17:28.117363 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.001
```

```
## 2026-03-04 18:17:28.300483 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:17:28.484566 INFO::Creating boxplot for categorical data, class vs Caldicoprobacter
```

```
## 2026-03-04 18:17:28.675498 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:17:28.852919 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:17:29.033648 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:17:29.21985 INFO::Creating boxplot for categorical data, class vs Izemoplasmatales
```

```
## 2026-03-04 18:17:29.398586 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:17:29.583859 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:17:29.768762 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:17:29.951349 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:17:30.130425 INFO::Creating boxplot for categorical data, class vs Streptococcus
```

```
## 2026-03-04 18:17:30.320312 INFO::Creating boxplot for categorical data, class vs Lactococcus
```

```
## 2026-03-04 18:17:30.503991 INFO::Creating boxplot for categorical data, class vs Acetatifactor
```

```
## 2026-03-04 18:17:30.683435 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae.Family
```

```
## 2026-03-04 18:17:30.875134 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:17:31.05829 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:17:31.239966 INFO::Creating boxplot for categorical data, class vs Desulfovibrionaceae.Family
```

```
## 2026-03-04 18:17:31.430343 INFO::Creating boxplot for categorical data, class vs Tuzzerella
```

```
## 2026-03-04 18:17:31.611256 INFO::Creating boxplot for categorical data, class vs Oxalobacter
```

```
## 2026-03-04 18:17:31.793259 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:17:31.981405 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._xylanophilum_group
```

```
## 2026-03-04 18:17:32.166843 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae
```

```
## 2026-03-04 18:17:32.362427 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:17:32.541612 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:17:32.72253 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:17:32.91307 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:17:33.104909 INFO::Creating boxplot for categorical data, class vs Coriobacteriales_Incertae_Sedis.Family
```

```
## 2026-03-04 18:17:33.292077 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae.Family
```

```
## 2026-03-04 18:17:33.483623 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:17:33.660311 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:17:33.849585 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:17:34.038623 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:17:34.218026 INFO::Creating boxplot for categorical data, class vs NK4A214_group
```

```
## 2026-03-04 18:17:34.402587 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:17:34.596137 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:17:34.773956 INFO::Creating boxplot for categorical data, class vs Eggerthella
```

```
## 2026-03-04 18:17:34.95206 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.006
```

```
## 2026-03-04 18:17:35.138948 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:17:35.315455 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_ND3007_group
```

```
## 2026-03-04 18:17:35.505398 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:17:35.700806 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._xylanophilum_group
```

```
## 2026-03-04 18:17:35.881487 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4B4_group
```

```
## 2026-03-04 18:17:36.06443 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:17:36.252426 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:17:36.435061 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:17:36.618213 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:17:36.813403 INFO::Creating boxplot for categorical data, class vs p.251.o5
```

```
## 2026-03-04 18:17:36.991102 INFO::Creating boxplot for categorical data, class vs Acetitomaculum
```

```
## 2026-03-04 18:17:37.18162 INFO::Creating boxplot for categorical data, class vs Candidatus_Saccharimonas
```

```
## 2026-03-04 18:17:37.361544 INFO::Creating boxplot for categorical data, class vs Dorea
```

```
## 2026-03-04 18:17:37.546019 INFO::Creating boxplot for categorical data, class vs Holdemania
```

```
## 2026-03-04 18:17:37.734164 INFO::Creating boxplot for categorical data, class vs UCG.010
```

```
## 2026-03-04 18:17:37.913384 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridium
```

```
## 2026-03-04 18:17:38.098083 INFO::Creating boxplot for categorical data, class vs Agathobacter
```

```
## 2026-03-04 18:17:38.338384 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:17:38.525407 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridium
```

```
## 2026-03-04 18:17:38.723853 INFO::Creating boxplot for categorical data, class vs Spirochaeta
```

```
## 2026-03-04 18:17:38.925304 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:17:39.112618 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:17:39.29233 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:17:39.488323 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:17:39.670811 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ventriosum_group
```

```
## 2026-03-04 18:17:39.855985 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:17:40.052474 INFO::Creating boxplot for categorical data, class vs Akkermansia
```

```
## 2026-03-04 18:17:40.239561 INFO::Creating boxplot for categorical data, class vs Clostridia.Class
```

```
## 2026-03-04 18:17:40.426691 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:17:40.620193 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.001
```

```
## 2026-03-04 18:17:40.801815 INFO::Creating boxplot for categorical data, class vs Pyramidobacter
```

```
## 2026-03-04 18:17:40.984303 INFO::Creating boxplot for categorical data, class vs Christensenellaceae.Family
```

```
## 2026-03-04 18:17:41.181778 INFO::Creating boxplot for categorical data, class vs Defluviitaleaceae_UCG.011
```

```
## 2026-03-04 18:17:41.365211 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:17:41.556242 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:17:41.763096 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:17:41.950645 INFO::Creating boxplot for categorical data, class vs X.Clostridium._methylpentosum_group
```

```
## 2026-03-04 18:17:42.133233 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:17:42.332731 INFO::Creating boxplot for categorical data, class vs Ureaplasma
```

```
## 2026-03-04 18:17:42.521102 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._nodatum_group
```

```
## 2026-03-04 18:17:42.706653 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:17:42.904458 INFO::Creating boxplot for categorical data, class vs Microtrichales.Order
```

```
## 2026-03-04 18:17:43.088685 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:17:43.269839 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.004
```

```
## 2026-03-04 18:17:43.471546 INFO::Creating boxplot for categorical data, class vs Firmicutes.Phylum
```

```
## 2026-03-04 18:17:43.657517 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:17:43.844245 INFO::Creating boxplot for categorical data, class vs Streptococcus
```

```
## 2026-03-04 18:17:44.044825 INFO::Creating boxplot for categorical data, class vs Sutterella
```

```
## 2026-03-04 18:17:44.225954 INFO::Creating boxplot for categorical data, class vs Erysipelotrichales.Order
```

```
## 2026-03-04 18:17:44.409919 INFO::Creating boxplot for categorical data, class vs Faecalibaculum
```

```
## 2026-03-04 18:17:44.607868 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._siraeum_group
```

```
## 2026-03-04 18:17:44.7943 INFO::Creating boxplot for categorical data, class vs Gordonibacter
```

```
## 2026-03-04 18:17:44.980242 INFO::Creating boxplot for categorical data, class vs Oscillospirales.Order
```

```
## 2026-03-04 18:17:45.18756 INFO::Creating boxplot for categorical data, class vs Odoribacter
```

```
## 2026-03-04 18:17:45.377063 INFO::Creating boxplot for categorical data, class vs Treponema
```

```
## 2026-03-04 18:17:45.57187 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:17:45.760253 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:17:45.94519 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:17:46.138727 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:17:46.323309 INFO::Creating boxplot for categorical data, class vs Oribacterium
```

```
## 2026-03-04 18:17:46.50442 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:17:46.705818 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:17:46.892989 INFO::Creating boxplot for categorical data, class vs Firmicutes.Phylum
```

```
## 2026-03-04 18:17:47.079032 INFO::Creating boxplot for categorical data, class vs Sporosarcina
```

```
## 2026-03-04 18:17:47.282019 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:17:47.466424 INFO::Creating boxplot for categorical data, class vs ASF356
```

```
## 2026-03-04 18:17:47.650692 INFO::Creating boxplot for categorical data, class vs Hydrogenoanaerobacterium
```

```
## 2026-03-04 18:17:47.854081 INFO::Creating boxplot for categorical data, class vs Paludibacteraceae.Family
```

```
## 2026-03-04 18:17:48.041541 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:17:48.230638 INFO::Creating boxplot for categorical data, class vs Parasutterella
```

```
## 2026-03-04 18:17:48.430766 INFO::Creating boxplot for categorical data, class vs Anaerovoracaceae.Family
```

```
## 2026-03-04 18:17:48.611636 INFO::Creating boxplot for categorical data, class vs Hungateiclostridiaceae.Family
```

```
## 2026-03-04 18:17:48.810605 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:17:49.001247 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:17:49.187128 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:17:49.384266 INFO::Creating boxplot for categorical data, class vs Blautia
```

```
## 2026-03-04 18:17:49.577258 INFO::Creating boxplot for categorical data, class vs Dubosiella
```

```
## 2026-03-04 18:17:49.767106 INFO::Creating boxplot for categorical data, class vs Pygmaiobacter
```

```
## 2026-03-04 18:17:50.011562 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:17:50.209008 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:17:50.392299 INFO::Creating boxplot for categorical data, class vs Spirochaetaceae.Family
```

```
## 2026-03-04 18:17:50.602903 INFO::Creating boxplot for categorical data, class vs Incertae_Sedis
```

```
## 2026-03-04 18:17:50.791219 INFO::Creating boxplot for categorical data, class vs Streptococcus
```

```
## 2026-03-04 18:17:50.97464 INFO::Creating boxplot for categorical data, class vs Butyricicoccaceae.Family
```

```
## 2026-03-04 18:17:51.17617 INFO::Creating boxplot for categorical data, class vs Faecalibaculum
```

```
## 2026-03-04 18:17:51.363785 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:17:51.545297 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:17:51.749933 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:17:51.94232 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:17:52.124635 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:17:52.329327 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:17:52.518489 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:17:52.701498 INFO::Creating boxplot for categorical data, class vs CAG.352
```

```
## 2026-03-04 18:17:52.906178 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.001
```

```
## 2026-03-04 18:17:53.093726 INFO::Creating boxplot for categorical data, class vs Enterococcus
```

```
## 2026-03-04 18:17:53.281617 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:17:53.486401 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:17:53.670936 INFO::Creating boxplot for categorical data, class vs Tuzzerella
```

```
## 2026-03-04 18:17:53.856069 INFO::Creating boxplot for categorical data, class vs Faecalibaculum
```

```
## 2026-03-04 18:17:54.070604 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:17:54.259585 INFO::Creating boxplot for categorical data, class vs ASF356
```

```
## 2026-03-04 18:17:54.442341 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:17:54.654295 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:17:54.846706 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:17:55.029498 INFO::Creating boxplot for categorical data, class vs Campylobacter
```

```
## 2026-03-04 18:17:55.240014 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:17:55.430973 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:17:55.613547 INFO::Creating boxplot for categorical data, class vs UCG.009
```

```
## 2026-03-04 18:17:55.819048 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:17:56.004527 INFO::Creating boxplot for categorical data, class vs Lactiplantibacillus
```

```
## 2026-03-04 18:17:56.184876 INFO::Creating boxplot for categorical data, class vs Bacteroides
```

```
## 2026-03-04 18:17:56.393201 INFO::Creating boxplot for categorical data, class vs UCG.009
```

```
## 2026-03-04 18:17:56.581326 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:17:56.764256 INFO::Creating boxplot for categorical data, class vs Papillibacter
```

```
## 2026-03-04 18:17:56.97157 INFO::Creating boxplot for categorical data, class vs Blautia
```

```
## 2026-03-04 18:17:57.158629 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:17:57.358547 INFO::Creating boxplot for categorical data, class vs Blautia
```

```
## 2026-03-04 18:17:57.555029 INFO::Creating boxplot for categorical data, class vs Clostridia_vadinBB60_group
```

```
## 2026-03-04 18:17:57.741011 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ruminantium_group
```

```
## 2026-03-04 18:17:57.949469 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:17:58.136726 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:17:58.322849 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:17:58.537938 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:17:58.730636 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:17:58.916335 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:17:59.12609 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:17:59.323281 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:17:59.516 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:17:59.727105 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:17:59.914781 INFO::Creating boxplot for categorical data, class vs Bacteroides
```

```
## 2026-03-04 18:18:00.108602 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:18:00.320295 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:18:00.512073 INFO::Creating boxplot for categorical data, class vs Prevotella_9
```

```
## 2026-03-04 18:18:00.716945 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:18:00.907396 INFO::Creating boxplot for categorical data, class vs Oscillospirales.Order
```

```
## 2026-03-04 18:18:01.094964 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:18:01.300347 INFO::Creating boxplot for categorical data, class vs Peptococcus
```

```
## 2026-03-04 18:18:01.489053 INFO::Creating boxplot for categorical data, class vs Rs.K70_termite_group
```

```
## 2026-03-04 18:18:01.678413 INFO::Creating boxplot for categorical data, class vs Bacteroidales_RF16_group
```

```
## 2026-03-04 18:18:01.931037 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:18:02.120608 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:18:02.319243 INFO::Creating boxplot for categorical data, class vs Prevotella
```

```
## 2026-03-04 18:18:02.53118 INFO::Creating boxplot for categorical data, class vs Faecalibaculum
```

```
## 2026-03-04 18:18:02.729821 INFO::Creating boxplot for categorical data, class vs Coriobacteriales_Incertae_Sedis.Family
```

```
## 2026-03-04 18:18:02.9172 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae.Family
```

```
## 2026-03-04 18:18:03.136441 INFO::Creating boxplot for categorical data, class vs GCA.900066575
```

```
## 2026-03-04 18:18:03.331505 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:18:03.525297 INFO::Creating boxplot for categorical data, class vs Bacteroidia.Class
```

```
## 2026-03-04 18:18:03.730991 INFO::Creating boxplot for categorical data, class vs p.2534.18B5_gut_group
```

```
## 2026-03-04 18:18:03.920922 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:18:04.106771 INFO::Creating boxplot for categorical data, class vs Candidatus_Arthromitus
```

```
## 2026-03-04 18:18:04.312095 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._xylanophilum_group
```

```
## 2026-03-04 18:18:04.505682 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:18:04.696234 INFO::Creating boxplot for categorical data, class vs Peptostreptococcales.Tissierellales.Order
```

```
## 2026-03-04 18:18:04.898357 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:18:05.087323 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:18:05.273599 INFO::Creating boxplot for categorical data, class vs Papillibacter
```

```
## 2026-03-04 18:18:05.48139 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:18:05.67641 INFO::Creating boxplot for categorical data, class vs Anaerostipes
```

```
## 2026-03-04 18:18:05.86886 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:18:06.076901 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:18:06.27314 INFO::Creating boxplot for categorical data, class vs Dubosiella
```

```
## 2026-03-04 18:18:06.462071 INFO::Creating boxplot for categorical data, class vs Peptococcus
```

```
## 2026-03-04 18:18:06.675839 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:18:06.871351 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:18:07.064183 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:18:07.277892 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:18:07.472982 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:18:07.660646 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:18:07.8741 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:18:08.063407 INFO::Creating boxplot for categorical data, class vs Gastranaerophilales
```

```
## 2026-03-04 18:18:08.270855 INFO::Creating boxplot for categorical data, class vs Microbacteriaceae.Family
```

```
## 2026-03-04 18:18:08.46165 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:18:08.649631 INFO::Creating boxplot for categorical data, class vs Faecalibaculum
```

```
## 2026-03-04 18:18:08.853497 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:18:09.049009 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:18:09.247195 INFO::Creating boxplot for categorical data, class vs Olsenella
```

```
## 2026-03-04 18:18:09.457926 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:18:09.652473 INFO::Creating boxplot for categorical data, class vs Acetatifactor
```

```
## 2026-03-04 18:18:09.842332 INFO::Creating boxplot for categorical data, class vs Mycoplasma
```

```
## 2026-03-04 18:18:10.046945 INFO::Creating boxplot for categorical data, class vs Candidatus_Soleaferrea
```

```
## 2026-03-04 18:18:10.243393 INFO::Creating boxplot for categorical data, class vs Coriobacteriia.Class
```

```
## 2026-03-04 18:18:10.436366 INFO::Creating boxplot for categorical data, class vs DTU089
```

```
## 2026-03-04 18:18:10.656355 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4B4_group
```

```
## 2026-03-04 18:18:10.862541 INFO::Creating boxplot for categorical data, class vs Thermoanaerobacter
```

```
## 2026-03-04 18:18:11.056023 INFO::Creating boxplot for categorical data, class vs Butyrivibrio
```

```
## 2026-03-04 18:18:11.267758 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ventriosum_group
```

```
## 2026-03-04 18:18:11.464324 INFO::Creating boxplot for categorical data, class vs Oxalobacteraceae.Family
```

```
## 2026-03-04 18:18:11.657221 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:18:11.862797 INFO::Creating boxplot for categorical data, class vs Coriobacteriaceae_UCG.002
```

```
## 2026-03-04 18:18:12.050936 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:18:12.254683 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:18:12.443407 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:18:12.632592 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:18:12.838103 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:18:13.033776 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:18:13.227681 INFO::Creating boxplot for categorical data, class vs SP3.e08
```

```
## 2026-03-04 18:18:13.42931 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:18:13.618076 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:18:13.812315 INFO::Creating boxplot for categorical data, class vs SP3.e08
```

```
## 2026-03-04 18:18:14.057549 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:18:14.247232 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:18:14.436288 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:18:14.648723 INFO::Creating boxplot for categorical data, class vs UCG.001
```

```
## 2026-03-04 18:18:14.839534 INFO::Creating boxplot for categorical data, class vs Family_XIII_AD3011_group
```

```
## 2026-03-04 18:18:15.027315 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._coprostanoligenes_group
```

```
## 2026-03-04 18:18:15.238843 INFO::Creating boxplot for categorical data, class vs Bacillus
```

```
## 2026-03-04 18:18:15.437178 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.006
```

```
## 2026-03-04 18:18:15.626556 INFO::Creating boxplot for categorical data, class vs Eubacterium
```

```
## 2026-03-04 18:18:15.838862 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:18:16.02847 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:18:16.213808 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:18:16.423756 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:18:16.618836 INFO::Creating boxplot for categorical data, class vs Clostridiaceae.Family
```

```
## 2026-03-04 18:18:16.81323 INFO::Creating boxplot for categorical data, class vs Colidextribacter
```

```
## 2026-03-04 18:18:17.032791 INFO::Creating boxplot for categorical data, class vs CAG.352
```

```
## 2026-03-04 18:18:17.232012 INFO::Creating boxplot for categorical data, class vs Clostridium_sensu_stricto_1
```

```
## 2026-03-04 18:18:17.424178 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:18:17.636697 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:18:17.823682 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.001
```

```
## 2026-03-04 18:18:18.012041 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:18:18.226691 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:18:18.422273 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:18:18.627107 INFO::Creating boxplot for categorical data, class vs Oligosphaeraceae.Family
```

```
## 2026-03-04 18:18:18.832554 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:18:19.025049 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:18:19.227348 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:18:19.429474 INFO::Creating boxplot for categorical data, class vs Oscillibacter
```

```
## 2026-03-04 18:18:19.625074 INFO::Creating boxplot for categorical data, class vs Clostridiales.Order
```

```
## 2026-03-04 18:18:19.831368 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:18:20.024652 INFO::Creating boxplot for categorical data, class vs Coriobacteriales.Order
```

```
## 2026-03-04 18:18:20.217872 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:18:20.43547 INFO::Creating boxplot for categorical data, class vs NK4A214_group
```

```
## 2026-03-04 18:18:20.630577 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._siraeum_group
```

```
## 2026-03-04 18:18:20.821984 INFO::Creating boxplot for categorical data, class vs Pelomonas
```

```
## 2026-03-04 18:18:21.030174 INFO::Creating boxplot for categorical data, class vs Oscillibacter
```

```
## 2026-03-04 18:18:21.23223 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:18:21.478317 INFO::Creating boxplot for categorical data, class vs UCG.005
```

```
## 2026-03-04 18:18:21.752373 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:18:22.013669 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:18:22.247153 INFO::Creating boxplot for categorical data, class vs Marvinbryantia
```

```
## 2026-03-04 18:18:22.48729 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:18:22.721022 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:18:22.959163 INFO::Creating boxplot for categorical data, class vs Rikenellaceae_RC9_gut_group
```

```
## 2026-03-04 18:18:23.2158 INFO::Creating boxplot for categorical data, class vs DNF00809
```

```
## 2026-03-04 18:18:23.427621 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:18:23.644661 INFO::Creating boxplot for categorical data, class vs Shuttleworthia
```

```
## 2026-03-04 18:18:23.883663 INFO::Creating boxplot for categorical data, class vs Parasutterella
```

```
## 2026-03-04 18:18:24.092792 INFO::Creating boxplot for categorical data, class vs Caulobacteraceae.Family
```

```
## 2026-03-04 18:18:24.325443 INFO::Creating boxplot for categorical data, class vs Desulfovibrionaceae.Family
```

```
## 2026-03-04 18:18:24.543389 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:18:24.745226 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:18:24.965938 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:18:25.159512 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:18:25.358086 INFO::Creating boxplot for categorical data, class vs Tyzzerella
```

```
## 2026-03-04 18:18:25.569197 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:18:25.761818 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:18:25.955471 INFO::Creating boxplot for categorical data, class vs Selenomonadaceae.Family
```

```
## 2026-03-04 18:18:26.180119 INFO::Creating boxplot for categorical data, class vs Anaerovoracaceae.Family
```

```
## 2026-03-04 18:18:26.377711 INFO::Creating boxplot for categorical data, class vs UBA1819
```

```
## 2026-03-04 18:18:26.644369 INFO::Creating boxplot for categorical data, class vs Coriobacteriaceae.Family
```

```
## 2026-03-04 18:18:26.844245 INFO::Creating boxplot for categorical data, class vs Pygmaiobacter
```

```
## 2026-03-04 18:18:27.044898 INFO::Creating boxplot for categorical data, class vs Pygmaiobacter
```

```
## 2026-03-04 18:18:27.260854 INFO::Creating boxplot for categorical data, class vs Bradyrhizobium
```

```
## 2026-03-04 18:18:27.45677 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:18:27.652298 INFO::Creating boxplot for categorical data, class vs Candidatus_Stoquefichus
```

```
## 2026-03-04 18:18:27.858864 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:18:28.058621 INFO::Creating boxplot for categorical data, class vs Rhodospirillales.Order
```

```
## 2026-03-04 18:18:28.251187 INFO::Creating boxplot for categorical data, class vs Paracaedibacteraceae.Family
```

```
## 2026-03-04 18:18:28.463908 INFO::Creating boxplot for categorical data, class vs X.Ruminococcus._torques_group
```

```
## 2026-03-04 18:18:28.658691 INFO::Creating boxplot for categorical data, class vs Veillonella
```

```
## 2026-03-04 18:18:28.852495 INFO::Creating boxplot for categorical data, class vs Romboutsia
```

```
## 2026-03-04 18:18:29.064185 INFO::Creating boxplot for categorical data, class vs Proteus
```

```
## 2026-03-04 18:18:29.257679 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:18:29.449381 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:18:29.660108 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:18:29.851954 INFO::Creating boxplot for categorical data, class vs Blautia
```

```
## 2026-03-04 18:18:30.047719 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:18:30.259839 INFO::Creating boxplot for categorical data, class vs NK4A214_group
```

```
## 2026-03-04 18:18:30.456268 INFO::Creating boxplot for categorical data, class vs Actinobacteriota.Phylum
```

```
## 2026-03-04 18:18:30.655233 INFO::Creating boxplot for categorical data, class vs Monoglobus
```

```
## 2026-03-04 18:18:30.865058 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:18:31.061191 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:18:31.253795 INFO::Creating boxplot for categorical data, class vs Desulfovibrio
```

```
## 2026-03-04 18:18:31.465572 INFO::Creating boxplot for categorical data, class vs SP3.e08
```

```
## 2026-03-04 18:18:31.66506 INFO::Creating boxplot for categorical data, class vs Candidatus_Arthromitus
```

```
## 2026-03-04 18:18:31.858297 INFO::Creating boxplot for categorical data, class vs UCG.002
```

```
## 2026-03-04 18:18:32.068802 INFO::Creating boxplot for categorical data, class vs Parvibacter
```

```
## 2026-03-04 18:18:32.267723 INFO::Creating boxplot for categorical data, class vs Syntrophomonadaceae.Family
```

```
## 2026-03-04 18:18:32.456478 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:18:32.672098 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:18:32.865272 INFO::Creating boxplot for categorical data, class vs Bilophila
```

```
## 2026-03-04 18:18:33.05517 INFO::Creating boxplot for categorical data, class vs Eggerthellaceae.Family
```

```
## 2026-03-04 18:18:33.267853 INFO::Creating boxplot for categorical data, class vs Marvinbryantia
```

```
## 2026-03-04 18:18:33.465659 INFO::Creating boxplot for categorical data, class vs Colidextribacter
```

```
## 2026-03-04 18:18:33.658917 INFO::Creating boxplot for categorical data, class vs Anaerofilum
```

```
## 2026-03-04 18:18:33.874342 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._nodatum_group
```

```
## 2026-03-04 18:18:34.067821 INFO::Creating boxplot for categorical data, class vs Firmicutes.Phylum
```

```
## 2026-03-04 18:18:34.268302 INFO::Creating boxplot for categorical data, class vs Coriobacteriales_Incertae_Sedis.Family
```

```
## 2026-03-04 18:18:34.493506 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:18:34.694448 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:18:34.912251 INFO::Creating boxplot for categorical data, class vs Enterorhabdus
```

```
## 2026-03-04 18:18:35.112513 INFO::Creating boxplot for categorical data, class vs Rikenellaceae.Family
```

```
## 2026-03-04 18:18:35.303631 INFO::Creating boxplot for categorical data, class vs Bacteroidales.Order
```

```
## 2026-03-04 18:18:35.524432 INFO::Creating boxplot for categorical data, class vs Peptococcus
```

```
## 2026-03-04 18:18:35.725465 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:18:35.917804 INFO::Creating boxplot for categorical data, class vs Fournierella
```

```
## 2026-03-04 18:18:36.126177 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridiaceae.Family
```

```
## 2026-03-04 18:18:36.319815 INFO::Creating boxplot for categorical data, class vs Muribaculaceae.Family
```

```
## 2026-03-04 18:18:36.511142 INFO::Creating boxplot for categorical data, class vs Escherichia.Shigella
```

```
## 2026-03-04 18:18:36.728663 INFO::Creating boxplot for categorical data, class vs Weissella
```

```
## 2026-03-04 18:18:36.924838 INFO::Creating boxplot for categorical data, class vs Harryflintia
```

```
## 2026-03-04 18:18:37.115304 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:18:37.325219 INFO::Creating boxplot for categorical data, class vs Parasutterella
```

```
## 2026-03-04 18:18:37.521626 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:18:37.731601 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:18:37.927421 INFO::Creating boxplot for categorical data, class vs Angelakisella
```

```
## 2026-03-04 18:18:38.119074 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._nodatum_group
```

```
## 2026-03-04 18:18:38.33237 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:18:38.531047 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._coprostanoligenes_group
```

```
## 2026-03-04 18:18:38.73168 INFO::Creating boxplot for categorical data, class vs Parasutterella
```

```
## 2026-03-04 18:18:39.382292 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:18:39.556406 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:18:39.737887 INFO::Creating boxplot for categorical data, class vs GWE2.42.42
```

```
## 2026-03-04 18:18:39.929667 INFO::Creating boxplot for categorical data, class vs Negativibacillus
```

```
## 2026-03-04 18:18:40.103682 INFO::Creating boxplot for categorical data, class vs Acetitomaculum
```

```
## 2026-03-04 18:18:40.281806 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:18:40.467698 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:18:40.644325 INFO::Creating boxplot for categorical data, class vs Family_XIII_AD3011_group
```

```
## 2026-03-04 18:18:40.826584 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:18:41.015304 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:18:41.200226 INFO::Creating boxplot for categorical data, class vs Halarsenatibacter
```

```
## 2026-03-04 18:18:41.38129 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:18:41.576283 INFO::Creating boxplot for categorical data, class vs Shuttleworthia
```

```
## 2026-03-04 18:18:41.758932 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:18:41.945234 INFO::Creating boxplot for categorical data, class vs Caulobacter
```

```
## 2026-03-04 18:18:42.13687 INFO::Creating boxplot for categorical data, class vs Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium
```

```
## 2026-03-04 18:18:42.322823 INFO::Creating boxplot for categorical data, class vs Rhodanobacter
```

```
## 2026-03-04 18:18:42.499793 INFO::Creating boxplot for categorical data, class vs p.1088.a5_gut_group
```

```
## 2026-03-04 18:18:42.686811 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:18:42.865824 INFO::Creating boxplot for categorical data, class vs Acetitomaculum
```

```
## 2026-03-04 18:18:43.047546 INFO::Creating boxplot for categorical data, class vs Firmicutes.Phylum
```

```
## 2026-03-04 18:18:43.239352 INFO::Creating boxplot for categorical data, class vs Parvibacter
```

```
## 2026-03-04 18:18:43.419607 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.006
```

```
## 2026-03-04 18:18:43.605306 INFO::Creating boxplot for categorical data, class vs Prevotellaceae.Family
```

```
## 2026-03-04 18:18:43.799223 INFO::Creating boxplot for categorical data, class vs Hungatella
```

```
## 2026-03-04 18:18:43.980334 INFO::Creating boxplot for categorical data, class vs Sphingomonadaceae.Family
```

```
## 2026-03-04 18:18:44.163396 INFO::Creating boxplot for categorical data, class vs Jeotgalicoccus
```

```
## 2026-03-04 18:18:44.362046 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:18:44.546775 INFO::Creating boxplot for categorical data, class vs Aerococcus
```

```
## 2026-03-04 18:18:44.731483 INFO::Creating boxplot for categorical data, class vs UCG.002
```

```
## 2026-03-04 18:18:44.919841 INFO::Creating boxplot for categorical data, class vs UBA1819
```

```
## 2026-03-04 18:18:45.103219 INFO::Creating boxplot for categorical data, class vs Bacilli.Class
```

```
## 2026-03-04 18:18:45.285361 INFO::Creating boxplot for categorical data, class vs Tyzzerella
```

```
## 2026-03-04 18:18:45.477442 INFO::Creating boxplot for categorical data, class vs Parasutterella
```

```
## 2026-03-04 18:18:45.656675 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae.Family
```

```
## 2026-03-04 18:18:45.840363 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:18:46.025472 INFO::Creating boxplot for categorical data, class vs GCA.900066575
```

```
## 2026-03-04 18:18:46.204125 INFO::Creating boxplot for categorical data, class vs RF39
```

```
## 2026-03-04 18:18:46.383694 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:18:46.569647 INFO::Creating boxplot for categorical data, class vs Butyricicoccus
```

```
## 2026-03-04 18:18:46.750973 INFO::Creating boxplot for categorical data, class vs RS62_marine_group
```

```
## 2026-03-04 18:18:46.940969 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:18:47.117182 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:18:47.295258 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:18:47.496432 INFO::Creating boxplot for categorical data, class vs Pasteurellaceae.Family
```

```
## 2026-03-04 18:18:47.681193 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:18:47.865857 INFO::Creating boxplot for categorical data, class vs Oscillibacter
```

```
## 2026-03-04 18:18:48.060082 INFO::Creating boxplot for categorical data, class vs Family_XIII_UCG.001
```

```
## 2026-03-04 18:18:48.239117 INFO::Creating boxplot for categorical data, class vs X.Clostridium._methylpentosum_group
```

```
## 2026-03-04 18:18:48.420215 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._brachy_group
```

```
## 2026-03-04 18:18:48.628373 INFO::Creating boxplot for categorical data, class vs Nocardioides
```

```
## 2026-03-04 18:18:48.825917 INFO::Creating boxplot for categorical data, class vs Pygmaiobacter
```

```
## 2026-03-04 18:18:49.02534 INFO::Creating boxplot for categorical data, class vs Tyzzerella
```

```
## 2026-03-04 18:18:49.235583 INFO::Creating boxplot for categorical data, class vs Sanguibacteroides
```

```
## 2026-03-04 18:18:49.43643 INFO::Creating boxplot for categorical data, class vs Anaerofustis
```

```
## 2026-03-04 18:18:49.645792 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:18:49.856212 INFO::Creating boxplot for categorical data, class vs Bacteria.Kingdom
```

```
## 2026-03-04 18:18:50.053031 INFO::Creating boxplot for categorical data, class vs GWE2.42.42
```

```
## 2026-03-04 18:18:50.255952 INFO::Creating boxplot for categorical data, class vs Ruminococcus
```

```
## 2026-03-04 18:18:50.449708 INFO::Creating boxplot for categorical data, class vs Colidextribacter
```

```
## 2026-03-04 18:18:50.643328 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:18:50.9067 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:18:51.100722 INFO::Creating boxplot for categorical data, class vs Acetitomaculum
```

```
## 2026-03-04 18:18:51.299884 INFO::Creating boxplot for categorical data, class vs Coriobacteriales_Incertae_Sedis.Family
```

```
## 2026-03-04 18:18:51.509534 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:18:51.693624 INFO::Creating boxplot for categorical data, class vs Clostridium_sensu_stricto_1
```

```
## 2026-03-04 18:18:51.876478 INFO::Creating boxplot for categorical data, class vs GCA.900066575
```

```
## 2026-03-04 18:18:52.075387 INFO::Creating boxplot for categorical data, class vs Pelomonas
```

```
## 2026-03-04 18:18:52.265588 INFO::Creating boxplot for categorical data, class vs Burkholderiales.Order
```

```
## 2026-03-04 18:18:52.447559 INFO::Creating boxplot for categorical data, class vs Bilophila
```

```
## 2026-03-04 18:18:52.642748 INFO::Creating boxplot for categorical data, class vs Actinobacteria.Class
```

```
## 2026-03-04 18:18:52.826991 INFO::Creating boxplot for categorical data, class vs Bacillus
```

```
## 2026-03-04 18:18:53.00825 INFO::Creating boxplot for categorical data, class vs Peptostreptococcaceae.Family
```

```
## 2026-03-04 18:18:53.204936 INFO::Creating boxplot for categorical data, class vs Bacteroides
```

```
## 2026-03-04 18:18:53.392306 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:18:53.577949 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae.Family
```

```
## 2026-03-04 18:18:53.781461 INFO::Creating boxplot for categorical data, class vs Halarsenatibacter
```

```
## 2026-03-04 18:18:53.968715 INFO::Creating boxplot for categorical data, class vs UCG.010
```

```
## 2026-03-04 18:18:54.150022 INFO::Creating boxplot for categorical data, class vs Erysipelotrichaceae.Family
```

```
## 2026-03-04 18:18:54.352963 INFO::Creating boxplot for categorical data, class vs DTU089
```

```
## 2026-03-04 18:18:54.543432 INFO::Creating boxplot for categorical data, class vs Harryflintia
```

```
## 2026-03-04 18:18:54.728394 INFO::Creating boxplot for categorical data, class vs Rikenellaceae.Family
```

```
## 2026-03-04 18:18:54.926912 INFO::Creating boxplot for categorical data, class vs Intestinimonas
```

```
## 2026-03-04 18:18:55.119535 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:18:55.307752 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:18:55.504989 INFO::Creating boxplot for categorical data, class vs Allobaculum
```

```
## 2026-03-04 18:18:55.691094 INFO::Creating boxplot for categorical data, class vs Tannerellaceae
```

```
## 2026-03-04 18:18:55.87115 INFO::Creating boxplot for categorical data, class vs Weissella
```

```
## 2026-03-04 18:18:56.067096 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:18:56.254102 INFO::Creating boxplot for categorical data, class vs WCHB1.41
```

```
## 2026-03-04 18:18:56.4328 INFO::Creating boxplot for categorical data, class vs Prevotellaceae_UCG.001
```

```
## 2026-03-04 18:18:56.627973 INFO::Creating boxplot for categorical data, class vs UBA1819
```

```
## 2026-03-04 18:18:56.815932 INFO::Creating boxplot for categorical data, class vs Gastranaerophilales
```

```
## 2026-03-04 18:18:56.996706 INFO::Creating boxplot for categorical data, class vs Paracaedibacteraceae.Family
```

```
## 2026-03-04 18:18:57.194262 INFO::Creating boxplot for categorical data, class vs Paracaedibacteraceae.Family
```

```
## 2026-03-04 18:18:57.378361 INFO::Creating boxplot for categorical data, class vs Bacilli.Class
```

```
## 2026-03-04 18:18:57.566373 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:18:57.765379 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._brachy_group
```

```
## 2026-03-04 18:18:57.951201 INFO::Creating boxplot for categorical data, class vs Incertae_Sedis
```

```
## 2026-03-04 18:18:58.13347 INFO::Creating boxplot for categorical data, class vs Micrococcales.Order
```

```
## 2026-03-04 18:18:58.332289 INFO::Creating boxplot for categorical data, class vs CL500.29_marine_group
```

```
## 2026-03-04 18:18:58.514012 INFO::Creating boxplot for categorical data, class vs Micrococcaceae.Family
```

```
## 2026-03-04 18:18:58.707534 INFO::Creating boxplot for categorical data, class vs UCG.010
```

```
## 2026-03-04 18:18:58.90052 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae.Family
```

```
## 2026-03-04 18:18:59.084826 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:18:59.290774 INFO::Creating boxplot for categorical data, class vs Streptomyces
```

```
## 2026-03-04 18:18:59.477128 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:18:59.663628 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:18:59.866349 INFO::Creating boxplot for categorical data, class vs Polaromonas
```

```
## 2026-03-04 18:19:00.052159 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ruminantium_group
```

```
## 2026-03-04 18:19:00.237455 INFO::Creating boxplot for categorical data, class vs Tuzzerella
```

```
## 2026-03-04 18:19:00.436216 INFO::Creating boxplot for categorical data, class vs Clostridium_sensu_stricto_13
```

```
## 2026-03-04 18:19:00.623765 INFO::Creating boxplot for categorical data, class vs Tundrisphaera
```

```
## 2026-03-04 18:19:00.811951 INFO::Creating boxplot for categorical data, class vs Arthrobacter
```

```
## 2026-03-04 18:19:01.019415 INFO::Creating boxplot for categorical data, class vs Clostridium_sensu_stricto_2
```

```
## 2026-03-04 18:19:01.200412 INFO::Creating boxplot for categorical data, class vs Bacillales.Order
```

```
## 2026-03-04 18:19:01.379465 INFO::Creating boxplot for categorical data, class vs Aquisphaera
```

```
## 2026-03-04 18:19:01.578239 INFO::Creating boxplot for categorical data, class vs Planococcaceae.Family
```

```
## 2026-03-04 18:19:01.764193 INFO::Creating boxplot for categorical data, class vs Ezakiella
```

```
## 2026-03-04 18:19:01.965517 INFO::Creating boxplot for categorical data, class vs Aeromicrobium
```

```
## 2026-03-04 18:19:02.152824 INFO::Creating boxplot for categorical data, class vs Acholeplasmataceae.Family
```

```
## 2026-03-04 18:19:02.338419 INFO::Creating boxplot for categorical data, class vs Cellulomonas
```

```
## 2026-03-04 18:19:02.569548 INFO::Creating boxplot for categorical data, class vs Pseudonocardia
```

```
## 2026-03-04 18:19:02.760846 INFO::Creating boxplot for categorical data, class vs Fructilactobacillus
```

```
## 2026-03-04 18:19:02.944048 INFO::Creating boxplot for categorical data, class vs Family_XI.Family
```

```
## 2026-03-04 18:19:03.147261 INFO::Creating boxplot for categorical data, class vs Gastranaerophilales
```

```
## 2026-03-04 18:19:03.334653 INFO::Creating boxplot for categorical data, class vs Sanguibacteroides
```

```
## 2026-03-04 18:19:03.526963 INFO::Creating boxplot for categorical data, class vs Fusobacterium
```

```
## 2026-03-04 18:19:03.730323 INFO::Creating boxplot for categorical data, class vs Gastranaerophilales
```

```
## 2026-03-04 18:19:03.918292 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._coprostanoligenes_group
```

```
## 2026-03-04 18:19:04.108883 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:19:04.320398 INFO::Creating boxplot for categorical data, class vs UBA1819
```

```
## 2026-03-04 18:19:04.519922 INFO::Creating boxplot for categorical data, class vs X.Clostridium._methylpentosum_group
```

```
## 2026-03-04 18:19:04.71839 INFO::Creating boxplot for categorical data, class vs UCG.002
```

```
## 2026-03-04 18:19:04.929904 INFO::Creating boxplot for categorical data, class vs Paracaedibacteraceae.Family
```

```
## 2026-03-04 18:19:05.122839 INFO::Creating boxplot for categorical data, class vs GCA.900066755
```

```
## 2026-03-04 18:19:05.30969 INFO::Creating boxplot for categorical data, class vs Bilophila
```

```
## 2026-03-04 18:19:05.514879 INFO::Creating boxplot for categorical data, class vs Pygmaiobacter
```

```
## 2026-03-04 18:19:05.701227 INFO::Creating boxplot for categorical data, class vs Actinobacteria.Class
```

```
## 2026-03-04 18:19:05.88905 INFO::Creating boxplot for categorical data, class vs Frisingicoccus
```

```
## 2026-03-04 18:19:06.094191 INFO::Creating boxplot for categorical data, class vs Pseudoxanthomonas
```

```
## 2026-03-04 18:19:06.279156 INFO::Creating boxplot for categorical data, class vs Rikenellaceae.Family
```

```
## 2026-03-04 18:19:06.460998 INFO::Creating boxplot for categorical data, class vs Glutamicibacter
```

```
## 2026-03-04 18:19:06.671172 INFO::Creating boxplot for categorical data, class vs Pedobacter
```

```
## 2026-03-04 18:19:06.857717 INFO::Creating boxplot for categorical data, class vs UCG.011
```

```
## 2026-03-04 18:19:07.045073 INFO::Creating boxplot for categorical data, class vs Candidatus_Arthromitus
```

```
## 2026-03-04 18:19:07.255203 INFO::Creating boxplot for categorical data, class vs Christensenellaceae.Family
```

```
## 2026-03-04 18:19:07.471203 INFO::Creating boxplot for categorical data, class vs Desulfovibrionaceae.Family
```

```
## 2026-03-04 18:19:07.660253 INFO::Creating boxplot for categorical data, class vs Christensenellaceae.Family
```

```
## 2026-03-04 18:19:07.866417 INFO::Creating boxplot for categorical data, class vs Microbacteriaceae.Family
```

```
## 2026-03-04 18:19:08.055714 INFO::Creating boxplot for categorical data, class vs Incertae_Sedis
```

```
## 2026-03-04 18:19:08.24311 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridium
```

```
## 2026-03-04 18:19:08.449713 INFO::Creating boxplot for categorical data, class vs Weissella
```

```
## 2026-03-04 18:19:08.641316 INFO::Creating boxplot for categorical data, class vs Weissella
```

```
## 2026-03-04 18:19:08.822387 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:19:09.030952 INFO::Creating boxplot for categorical data, class vs Pasteurellaceae.Family
```

```
## 2026-03-04 18:19:09.218018 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae.Family
```

```
## 2026-03-04 18:19:09.407409 INFO::Creating boxplot for categorical data, class vs Variovorax
```

```
## 2026-03-04 18:19:09.619078 INFO::Creating boxplot for categorical data, class vs Acetanaerobacterium
```

```
## 2026-03-04 18:19:09.809486 INFO::Creating boxplot for categorical data, class vs Desulfovibrio
```

```
## 2026-03-04 18:19:10.006093 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae.Family
```

```
## 2026-03-04 18:19:10.211215 INFO::Creating boxplot for categorical data, class vs NK4A214_group
```

```
## 2026-03-04 18:19:10.401061 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._nodatum_group
```

```
## 2026-03-04 18:19:10.608855 INFO::Creating boxplot for categorical data, class vs GWE2.42.42
```

```
## 2026-03-04 18:19:10.804717 INFO::Creating boxplot for categorical data, class vs Lactobacillus
```

```
## 2026-03-04 18:19:10.993958 INFO::Creating boxplot for categorical data, class vs Caproiciproducens
```

```
## 2026-03-04 18:19:11.208967 INFO::Creating boxplot for categorical data, class vs Helicobacter
```

```
## 2026-03-04 18:19:11.404088 INFO::Creating boxplot for categorical data, class vs Bilophila
```

```
## 2026-03-04 18:19:11.599625 INFO::Creating boxplot for categorical data, class vs Bacteroidales.Order
```

```
## 2026-03-04 18:19:11.817455 INFO::Creating boxplot for categorical data, class vs Phascolarctobacterium
```

```
## 2026-03-04 18:19:12.018611 INFO::Creating boxplot for categorical data, class vs Halarsenatibacter
```

```
## 2026-03-04 18:19:12.214823 INFO::Creating boxplot for categorical data, class vs Rickettsiales.Order
```

```
## 2026-03-04 18:19:12.429921 INFO::Creating boxplot for categorical data, class vs UCG.010
```

```
## 2026-03-04 18:19:12.627823 INFO::Creating boxplot for categorical data, class vs Desulfovibrionaceae.Family
```

```
## 2026-03-04 18:19:12.818038 INFO::Creating boxplot for categorical data, class vs Rhodococcus
```

```
## 2026-03-04 18:19:13.031216 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:19:13.219768 INFO::Creating boxplot for categorical data, class vs Actinomyces
```

```
## 2026-03-04 18:19:13.416246 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._brachy_group
```

```
## 2026-03-04 18:19:13.617167 INFO::Creating boxplot for categorical data, class vs Weissella
```

```
## 2026-03-04 18:19:13.801964 INFO::Creating boxplot for categorical data, class vs Staphylococcus
```

```
## 2026-03-04 18:19:14.013784 INFO::Creating boxplot for categorical data, class vs X67.14
```

```
## 2026-03-04 18:19:14.204846 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NC2004_group
```

```
## 2026-03-04 18:19:14.395646 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae
```

```
## 2026-03-04 18:19:14.66196 INFO::Creating boxplot for categorical data, class vs UCG.002
```

```
## 2026-03-04 18:19:14.861725 INFO::Creating boxplot for categorical data, class vs Candidatus_Soleaferrea
```

```
## 2026-03-04 18:19:15.066521 INFO::Creating boxplot for categorical data, class vs Raoultibacter
```

```
## 2026-03-04 18:19:15.280785 INFO::Creating boxplot for categorical data, class vs Paraclostridium
```

```
## 2026-03-04 18:19:15.482564 INFO::Creating boxplot for categorical data, class vs Micromonosporaceae.Family
```

```
## 2026-03-04 18:19:15.67321 INFO::Creating boxplot for categorical data, class vs Xanthobacteraceae.Family
```

```
## 2026-03-04 18:19:15.880971 INFO::Creating boxplot for categorical data, class vs Psychroglaciecola
```

```
## 2026-03-04 18:19:16.07549 INFO::Creating boxplot for categorical data, class vs Neisseria
```

```
## 2026-03-04 18:19:16.264958 INFO::Creating boxplot for categorical data, class vs Megamonas
```

```
## 2026-03-04 18:19:16.467659 INFO::Creating boxplot for categorical data, class vs Prevotella_7
```

```
## 2026-03-04 18:19:16.663674 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.008
```

```
## 2026-03-04 18:19:16.848595 INFO::Creating boxplot for categorical data, class vs Lachnospira
```

```
## 2026-03-04 18:19:17.06202 INFO::Creating boxplot for categorical data, class vs Comamonas
```

```
## 2026-03-04 18:19:17.258241 INFO::Creating boxplot for categorical data, class vs Harryflintia
```

```
## 2026-03-04 18:19:17.447736 INFO::Creating boxplot for categorical data, class vs GCA.900066755
```

```
## 2026-03-04 18:19:17.657911 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:19:17.849646 INFO::Creating boxplot for categorical data, class vs Anaerobiospirillum
```

```
## 2026-03-04 18:19:18.041011 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:19:18.24782 INFO::Creating boxplot for categorical data, class vs Bacilli.Class
```

```
## 2026-03-04 18:19:18.440118 INFO::Creating boxplot for categorical data, class vs Frisingicoccus
```

```
## 2026-03-04 18:19:18.624084 INFO::Creating boxplot for categorical data, class vs Erysipelatoclostridium
```

```
## 2026-03-04 18:19:18.827352 INFO::Creating boxplot for categorical data, class vs Family_XIII_AD3011_group
```

```
## 2026-03-04 18:19:19.020702 INFO::Creating boxplot for categorical data, class vs Desulfovibrionaceae.Family
```

```
## 2026-03-04 18:19:19.209339 INFO::Creating boxplot for categorical data, class vs Rhodoferax
```

```
## 2026-03-04 18:19:19.413765 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:19:19.611285 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae
```

```
## 2026-03-04 18:19:19.799582 INFO::Creating boxplot for categorical data, class vs Synergistaceae.Family
```

```
## 2026-03-04 18:19:20.007004 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ventriosum_group
```

```
## 2026-03-04 18:19:20.205516 INFO::Creating boxplot for categorical data, class vs Salinicoccus
```

```
## 2026-03-04 18:19:20.397275 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:19:20.604564 INFO::Creating boxplot for categorical data, class vs Flavobacterium
```

```
## 2026-03-04 18:19:20.806724 INFO::Creating boxplot for categorical data, class vs Fournierella
```

```
## 2026-03-04 18:19:20.995276 INFO::Creating boxplot for categorical data, class vs Faecalibacterium
```

```
## 2026-03-04 18:19:21.208702 INFO::Creating boxplot for categorical data, class vs Clostridia_vadinBB60_group
```

```
## 2026-03-04 18:19:21.405716 INFO::Creating boxplot for categorical data, class vs Pasteurellaceae.Family
```

```
## 2026-03-04 18:19:21.610039 INFO::Creating boxplot for categorical data, class vs Pasteurellaceae.Family
```

```
## 2026-03-04 18:19:21.803192 INFO::Creating boxplot for categorical data, class vs Hungateiclostridiaceae.Family
```

```
## 2026-03-04 18:19:21.99524 INFO::Creating boxplot for categorical data, class vs CAG.352
```

```
## 2026-03-04 18:19:22.201447 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:19:22.395761 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:19:22.592062 INFO::Creating boxplot for categorical data, class vs Paludicola
```

```
## 2026-03-04 18:19:22.799662 INFO::Creating boxplot for categorical data, class vs Puniceicoccaceae.Family
```

```
## 2026-03-04 18:19:23.000653 INFO::Creating boxplot for categorical data, class vs Sanguibacteroides
```

```
## 2026-03-04 18:19:23.197331 INFO::Creating boxplot for categorical data, class vs Chryseobacterium
```

```
## 2026-03-04 18:19:23.40863 INFO::Creating boxplot for categorical data, class vs GWE2.42.42
```

```
## 2026-03-04 18:19:23.620665 INFO::Creating boxplot for categorical data, class vs Pygmaiobacter
```

```
## 2026-03-04 18:19:23.818391 INFO::Creating boxplot for categorical data, class vs Acetitomaculum
```

```
## 2026-03-04 18:19:24.035647 INFO::Creating boxplot for categorical data, class vs UCG.002
```

```
## 2026-03-04 18:19:24.235074 INFO::Creating boxplot for categorical data, class vs Fusobacterium
```

```
## 2026-03-04 18:19:24.432956 INFO::Creating boxplot for categorical data, class vs Gastranaerophilales
```

```
## 2026-03-04 18:19:24.664778 INFO::Creating boxplot for categorical data, class vs Clostridium_sensu_stricto_11
```

```
## 2026-03-04 18:19:24.863652 INFO::Creating boxplot for categorical data, class vs Actinobacteria.Class
```

```
## 2026-03-04 18:19:25.071447 INFO::Creating boxplot for categorical data, class vs Phascolarctobacterium
```

```
## 2026-03-04 18:19:25.273088 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._nodatum_group
```

```
## 2026-03-04 18:19:25.4775 INFO::Creating boxplot for categorical data, class vs Halarsenatibacter
```

```
## 2026-03-04 18:19:25.688038 INFO::Creating boxplot for categorical data, class vs Halarsenatibacter
```

```
## 2026-03-04 18:19:25.88502 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._ventriosum_group
```

```
## 2026-03-04 18:19:26.084551 INFO::Creating boxplot for categorical data, class vs Clostridia_UCG.014
```

```
## 2026-03-04 18:19:26.306795 INFO::Creating boxplot for categorical data, class vs Romboutsia
```

```
## 2026-03-04 18:19:26.510287 INFO::Creating boxplot for categorical data, class vs Stenotrophomonas
```

```
## 2026-03-04 18:19:26.702592 INFO::Creating boxplot for categorical data, class vs Elusimicrobium
```

```
## 2026-03-04 18:19:26.958558 INFO::Creating boxplot for categorical data, class vs Coriobacteriales_Incertae_Sedis.Family
```

```
## 2026-03-04 18:19:27.161502 INFO::Creating boxplot for categorical data, class vs Staphylococcus
```

```
## 2026-03-04 18:19:27.362624 INFO::Creating boxplot for categorical data, class vs X.Clostridium._innocuum_group
```

```
## 2026-03-04 18:19:27.584876 INFO::Creating boxplot for categorical data, class vs Ileibacterium
```

```
## 2026-03-04 18:19:27.786651 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_NK4A136_group
```

```
## 2026-03-04 18:19:27.98094 INFO::Creating boxplot for categorical data, class vs Pasteurellaceae.Family
```

```
## 2026-03-04 18:19:28.19838 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:19:28.395409 INFO::Creating boxplot for categorical data, class vs Catenibacillus
```

```
## 2026-03-04 18:19:28.589993 INFO::Creating boxplot for categorical data, class vs Eubacteriaceae.Family
```

```
## 2026-03-04 18:19:28.810525 INFO::Creating boxplot for categorical data, class vs V9D2013_group
```

```
## 2026-03-04 18:19:29.009023 INFO::Creating boxplot for categorical data, class vs Family_XIII_UCG.001
```

```
## 2026-03-04 18:19:29.201248 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:19:29.415275 INFO::Creating boxplot for categorical data, class vs Anaerobiospirillum
```

```
## 2026-03-04 18:19:29.608178 INFO::Creating boxplot for categorical data, class vs Staphylococcus
```

```
## 2026-03-04 18:19:29.795131 INFO::Creating boxplot for categorical data, class vs Bacilli.Class
```

```
## 2026-03-04 18:19:30.007846 INFO::Creating boxplot for categorical data, class vs Desulfovibrio
```

```
## 2026-03-04 18:19:30.20157 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:19:30.390219 INFO::Creating boxplot for categorical data, class vs Muribaculaceae
```

```
## 2026-03-04 18:19:30.608314 INFO::Creating boxplot for categorical data, class vs Veillonellales.Selenomonadales.Order
```

```
## 2026-03-04 18:19:30.802133 INFO::Creating boxplot for categorical data, class vs GWE2.42.42
```

```
## 2026-03-04 18:19:30.993442 INFO::Creating boxplot for categorical data, class vs Anaerovorax
```

```
## 2026-03-04 18:19:31.207076 INFO::Creating boxplot for categorical data, class vs Flavobacterium
```

```
## 2026-03-04 18:19:31.396002 INFO::Creating boxplot for categorical data, class vs GCA.900066575
```

```
## 2026-03-04 18:19:31.584426 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._brachy_group
```

```
## 2026-03-04 18:19:31.801509 INFO::Creating boxplot for categorical data, class vs Acetatifactor
```

```
## 2026-03-04 18:19:31.996849 INFO::Creating boxplot for categorical data, class vs Salinicoccus
```

```
## 2026-03-04 18:19:32.197278 INFO::Creating boxplot for categorical data, class vs Escherichia.Shigella
```

```
## 2026-03-04 18:19:32.401152 INFO::Creating boxplot for categorical data, class vs vadinBE97
```

```
## 2026-03-04 18:19:32.595485 INFO::Creating boxplot for categorical data, class vs Halarsenatibacter
```

```
## 2026-03-04 18:19:32.801467 INFO::Creating boxplot for categorical data, class vs Firmicutes.Phylum
```

```
## 2026-03-04 18:19:32.999428 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_FCS020_group
```

```
## 2026-03-04 18:19:33.192021 INFO::Creating boxplot for categorical data, class vs hgcI_clade
```

```
## 2026-03-04 18:19:33.401911 INFO::Creating boxplot for categorical data, class vs Clostridia.Class
```

```
## 2026-03-04 18:19:33.595369 INFO::Creating boxplot for categorical data, class vs Mailhella
```

```
## 2026-03-04 18:19:33.78747 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:19:34.000458 INFO::Creating boxplot for categorical data, class vs Oscillospiraceae.Family
```

```
## 2026-03-04 18:19:34.193804 INFO::Creating boxplot for categorical data, class vs Ileibacterium
```

```
## 2026-03-04 18:19:34.387448 INFO::Creating boxplot for categorical data, class vs Caldicoprobacter
```

```
## 2026-03-04 18:19:34.607012 INFO::Creating boxplot for categorical data, class vs Ruminococcaceae.Family
```

```
## 2026-03-04 18:19:34.803866 INFO::Creating boxplot for categorical data, class vs Rikenella
```

```
## 2026-03-04 18:19:34.992824 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:19:35.208385 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._siraeum_group
```

```
## 2026-03-04 18:19:35.402102 INFO::Creating boxplot for categorical data, class vs Sanguibacteroides
```

```
## 2026-03-04 18:19:35.594976 INFO::Creating boxplot for categorical data, class vs UBA1819
```

```
## 2026-03-04 18:19:35.81865 INFO::Creating boxplot for categorical data, class vs UCG.009
```

```
## 2026-03-04 18:19:36.020505 INFO::Creating boxplot for categorical data, class vs Fusobacterium
```

```
## 2026-03-04 18:19:36.219166 INFO::Creating boxplot for categorical data, class vs Fusobacterium
```

```
## 2026-03-04 18:19:36.427794 INFO::Creating boxplot for categorical data, class vs GCA.900066755
```

```
## 2026-03-04 18:19:36.620032 INFO::Creating boxplot for categorical data, class vs Actinobacteria.Class
```

```
## 2026-03-04 18:19:36.829332 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:19:37.024531 INFO::Creating boxplot for categorical data, class vs Fusobacterium
```

```
## 2026-03-04 18:19:37.218476 INFO::Creating boxplot for categorical data, class vs Christensenella
```

```
## 2026-03-04 18:19:37.433198 INFO::Creating boxplot for categorical data, class vs Alistipes
```

```
## 2026-03-04 18:19:37.627106 INFO::Creating boxplot for categorical data, class vs Eggerthellaceae.Family
```

```
## 2026-03-04 18:19:37.823321 INFO::Creating boxplot for categorical data, class vs Paracaedibacteraceae.Family
```

```
## 2026-03-04 18:19:38.038918 INFO::Creating boxplot for categorical data, class vs Oscillospirales.Order
```

```
## 2026-03-04 18:19:38.23276 INFO::Creating boxplot for categorical data, class vs Romboutsia
```

```
## 2026-03-04 18:19:38.416053 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:19:38.63416 INFO::Creating boxplot for categorical data, class vs Desulfovibrionales.Order
```

```
## 2026-03-04 18:19:38.827086 INFO::Creating boxplot for categorical data, class vs Family_XIII_AD3011_group
```

```
## 2026-03-04 18:19:39.086534 INFO::Creating boxplot for categorical data, class vs Brevundimonas
```

```
## 2026-03-04 18:19:39.278308 INFO::Creating boxplot for categorical data, class vs Novosphingobium
```

```
## 2026-03-04 18:19:39.478034 INFO::Creating boxplot for categorical data, class vs Marmoricola
```

```
## 2026-03-04 18:19:39.69301 INFO::Creating boxplot for categorical data, class vs Luteipulveratus
```

```
## 2026-03-04 18:19:39.890791 INFO::Creating boxplot for categorical data, class vs Luteimonas
```

```
## 2026-03-04 18:19:40.082961 INFO::Creating boxplot for categorical data, class vs Arachidicoccus
```

```
## 2026-03-04 18:19:40.287618 INFO::Creating boxplot for categorical data, class vs Asticcacaulis
```

```
## 2026-03-04 18:19:40.48615 INFO::Creating boxplot for categorical data, class vs LWQ8
```

```
## 2026-03-04 18:19:40.681095 INFO::Creating boxplot for categorical data, class vs Ferruginibacter
```

```
## 2026-03-04 18:19:40.885379 INFO::Creating boxplot for categorical data, class vs Eisenbergiella
```

```
## 2026-03-04 18:19:41.089167 INFO::Creating boxplot for categorical data, class vs Altererythrobacter
```

```
## 2026-03-04 18:19:41.285963 INFO::Creating boxplot for categorical data, class vs Ilumatobacteraceae.Family
```

```
## 2026-03-04 18:19:41.497591 INFO::Creating boxplot for categorical data, class vs YC.ZSS.LKJ147
```

```
## 2026-03-04 18:19:41.699236 INFO::Creating boxplot for categorical data, class vs Terrimonas
```

```
## 2026-03-04 18:19:41.891842 INFO::Creating boxplot for categorical data, class vs WD2101_soil_group
```

```
## 2026-03-04 18:19:42.105785 INFO::Creating boxplot for categorical data, class vs Xanthomonadaceae.Family
```

```
## 2026-03-04 18:19:42.303893 INFO::Creating boxplot for categorical data, class vs Actinomycetospora
```

```
## 2026-03-04 18:19:42.499554 INFO::Creating boxplot for categorical data, class vs Luteibacter
```

```
## 2026-03-04 18:19:42.713448 INFO::Creating boxplot for categorical data, class vs IMCC26256
```

```
## 2026-03-04 18:19:42.911919 INFO::Creating boxplot for categorical data, class vs Micropepsis
```

```
## 2026-03-04 18:19:43.102775 INFO::Creating boxplot for categorical data, class vs Chryseolinea
```

```
## 2026-03-04 18:19:43.316662 INFO::Creating boxplot for categorical data, class vs BIyi10
```

```
## 2026-03-04 18:19:43.514794 INFO::Creating boxplot for categorical data, class vs Rhodopirellula
```

```
## 2026-03-04 18:19:43.713375 INFO::Creating boxplot for categorical data, class vs TM7a
```

```
## 2026-03-04 18:19:43.919135 INFO::Creating boxplot for categorical data, class vs Cytophaga
```

```
## 2026-03-04 18:19:44.114962 INFO::Creating boxplot for categorical data, class vs Sphingobacterium
```

```
## 2026-03-04 18:19:44.304387 INFO::Creating boxplot for categorical data, class vs Delftia
```

```
## 2026-03-04 18:19:44.511057 INFO::Creating boxplot for categorical data, class vs Achromobacter
```

```
## 2026-03-04 18:19:44.703289 INFO::Creating boxplot for categorical data, class vs Taibaiella
```

```
## 2026-03-04 18:19:44.90301 INFO::Creating boxplot for categorical data, class vs Massilia
```

```
## 2026-03-04 18:19:45.113154 INFO::Creating boxplot for categorical data, class vs Sphingobium
```

```
## 2026-03-04 18:19:45.307124 INFO::Creating boxplot for categorical data, class vs Dyadobacter
```

```
## 2026-03-04 18:19:45.501223 INFO::Creating boxplot for categorical data, class vs Microscillaceae.Family
```

```
## 2026-03-04 18:19:45.7182 INFO::Creating boxplot for categorical data, class vs Actinoallomurus
```

```
## 2026-03-04 18:19:45.91591 INFO::Creating boxplot for categorical data, class vs Cohnella
```

```
## 2026-03-04 18:19:46.109881 INFO::Creating boxplot for categorical data, class vs Roseomonas
```

```
## 2026-03-04 18:19:46.327594 INFO::Creating boxplot for categorical data, class vs Pseudosphingobacterium
```

```
## 2026-03-04 18:19:46.528123 INFO::Creating boxplot for categorical data, class vs Mailhella
```

```
## 2026-03-04 18:19:46.727027 INFO::Creating boxplot for categorical data, class vs Ligilactobacillus
```

```
## 2026-03-04 18:19:46.954013 INFO::Creating boxplot for categorical data, class vs Pseudomonas
```

```
## 2026-03-04 18:19:47.158796 INFO::Creating boxplot for categorical data, class vs Enterobacter
```

```
## 2026-03-04 18:19:47.351303 INFO::Creating boxplot for categorical data, class vs Catenibacillus
```

```
## 2026-03-04 18:19:47.563486 INFO::Creating boxplot for categorical data, class vs UCG.004
```

```
## 2026-03-04 18:19:47.756906 INFO::Creating boxplot for categorical data, class vs Sanguibacteroides
```

```
## 2026-03-04 18:19:47.969404 INFO::Creating boxplot for categorical data, class vs Family_XIII_UCG.001
```

```
## 2026-03-04 18:19:48.171052 INFO::Creating boxplot for categorical data, class vs Roseateles
```

```
## 2026-03-04 18:19:48.36276 INFO::Creating boxplot for categorical data, class vs Fusobacterium
```

```
## 2026-03-04 18:19:48.571921 INFO::Creating boxplot for categorical data, class vs Intestinimonas
```

```
## 2026-03-04 18:19:48.766062 INFO::Creating boxplot for categorical data, class vs Actinobacteria.Class
```

```
## 2026-03-04 18:19:48.958075 INFO::Creating boxplot for categorical data, class vs Atopobiaceae.Family
```

```
## 2026-03-04 18:19:49.178816 INFO::Creating boxplot for categorical data, class vs Anaerosporobacter
```

```
## 2026-03-04 18:19:49.371762 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.001
```

```
## 2026-03-04 18:19:49.562694 INFO::Creating boxplot for categorical data, class vs Christensenellaceae_R.7_group
```

```
## 2026-03-04 18:19:49.778674 INFO::Creating boxplot for categorical data, class vs Phascolarctobacterium
```

```
## 2026-03-04 18:19:49.97964 INFO::Creating boxplot for categorical data, class vs Limosilactobacillus
```

```
## 2026-03-04 18:19:50.173258 INFO::Creating boxplot for categorical data, class vs Anaeroplasma
```

```
## 2026-03-04 18:19:50.396608 INFO::Creating boxplot for categorical data, class vs Opitutales.Order
```

```
## 2026-03-04 18:19:50.59279 INFO::Creating boxplot for categorical data, class vs Bacteroides
```

```
## 2026-03-04 18:19:50.805033 INFO::Creating boxplot for categorical data, class vs Eggerthellaceae.Family
```

```
## 2026-03-04 18:19:51.000688 INFO::Creating boxplot for categorical data, class vs Pantoea
```

```
## 2026-03-04 18:19:51.193955 INFO::Creating boxplot for categorical data, class vs Intestinimonas
```

```
## 2026-03-04 18:19:51.835246 INFO::Creating boxplot for categorical data, class vs Rikenellaceae.Family
```

```
## 2026-03-04 18:19:52.009746 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._coprostanoligenes_group
```

```
## 2026-03-04 18:19:52.193063 INFO::Creating boxplot for categorical data, class vs X.Eubacterium._fissicatena_group
```

```
## 2026-03-04 18:19:52.383162 INFO::Creating boxplot for categorical data, class vs GCA.900066755
```

```
## 2026-03-04 18:19:52.561454 INFO::Creating boxplot for categorical data, class vs GCA.900066755
```

```
## 2026-03-04 18:19:52.739101 INFO::Creating boxplot for categorical data, class vs Tyzzerella
```

```
## 2026-03-04 18:19:52.92744 INFO::Creating boxplot for categorical data, class vs Anaerotruncus
```

```
## 2026-03-04 18:19:53.104974 INFO::Creating boxplot for categorical data, class vs Peptococcaceae.Family
```

```
## 2026-03-04 18:19:53.285709 INFO::Creating boxplot for categorical data, class vs Anaerobiospirillum
```

```
## 2026-03-04 18:19:53.471666 INFO::Creating boxplot for categorical data, class vs Weissella
```

```
## 2026-03-04 18:19:53.654105 INFO::Creating boxplot for categorical data, class vs Lachnospiraceae_UCG.006
```

```
## 2026-03-04 18:19:53.836881 INFO::Creating boxplot for categorical data, class vs Christensenella
```

```
## 2026-03-04 18:19:54.02373 INFO::Creating boxplot for categorical data, class vs Paenibacillus
```

```
## 2026-03-04 18:19:54.205769 INFO::Creating boxplot for categorical data, class vs Frisingicoccus
```

```
## 2026-03-04 18:19:54.387736 INFO::Creating boxplot for categorical data, class vs Frisingicoccus
```

```
## 2026-03-04 18:19:54.576907 INFO::Creating boxplot for categorical data, class vs Mucispirillum
```

```
## 2026-03-04 18:19:54.757572 INFO::Creating boxplot for categorical data, class vs Tannerellaceae.Family
```

```
## 2026-03-04 18:19:54.939675 INFO::Creating boxplot for categorical data, class vs Brachybacterium
```

```
## 2026-03-04 18:19:55.140785 INFO::Creating boxplot for categorical data, class vs Christensenellaceae.Family
```

```
## 2026-03-04 18:19:55.318013 INFO::Creating boxplot for categorical data, class vs Roseburia
```

```
## 2026-03-04 18:19:55.501503 INFO::Creating boxplot for categorical data, class vs Faecalibacterium
```

```
## 2026-03-04 18:19:55.691203 INFO::Creating boxplot for categorical data, class vs Corynebacterium
```

```
## 2026-03-04 18:19:55.877724 INFO::Creating boxplot for categorical data, class vs Salinicoccus
```

Save the fit data object as an rds file.


``` r
# 20260220_21_13_07
# saveRDS(maaslin.fit_data.host,
#         file = file.path("output/rdafiles",
#                          paste(
#                            paste(format(Sys.time(),format="%Y%m%d"),
#                                  format(Sys.time(),format = "%H_%M_%S"),
#                                  sep = "_"),"maaslin.fit_data.host",
#                            host.data.for_test$output.filename,".rds",sep = "-")))
```

Process the output 


``` r
# "20260213_12_20_50" for signif.tsv
# "20260213_12_20_51" for maaslin.signif.decreased.tsv
maaslin.processed_output.host <- process_maaslin2_output(agglom.rank = agglom.rank,
                        maaslin.fit_data = maaslin.fit_data.host,
                        ps.q.agg = host.data.for_test$ps.q.agg,
                        sample.groups = host.data.for_test$sample.groups)
```

Save the workspace.
20260213_12_20_28: genus workspace


``` r
# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   "maaslin2",host.data.for_test$output.filename,"workspace.RData",sep="-")))
```



### Run a test with ALDEx2.
### Create covariates: take the samples in the ps.q.df.wide, and find 
the positions of these samples in the custom.md. 
Ref: match returns a vector of the positions of (first) matches of its first argument in its second. 
We're searching the elements of the first vector in the second vector. 


``` r
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
```

Model matrix is covariates and the reference group is the first covariate.


``` r
mm <- model.matrix(~ covariates-1)
```

Reorder model.matrix to put ref.level as first column.


``` r
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
```

Run the test: Aldex glm for a complex case.


``` r
set.seed(1)
ps.q.aldex.clr <- aldex.clr(t(host.data.for_test$dataset), mm, 
                            mc.samples=1000, 
                            denom="all", verbose=F)
```

```
## integer matrix provided
```

```
## using all features for denominator
```

```
## operating in serial mode
```

```
## computing center with all features
```

``` r
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
```

Save the workspace just in case.


``` r
# save.image(paste0("./output/rdafiles/",
#                   paste(paste(format(Sys.time(),format="%Y%m%d"),
#                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                         "aldex2",host.data.for_test$output.filename,
#                         "workspace-test.RData",sep="-")))
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

```
## operating in serial mode
```

```
## sanity check complete
```

```
## rab.all  complete
```

```
## rab.win  complete
```

```
## rab of samples complete
```

```
## within sample difference calculated
```

```
## between group difference calculated
```

```
## group summaries calculated
```

```
## unpaired effect size calculated
```

```
## summarizing output
```

Save the workspace just in case.


``` r
# save.image(paste0("./output/rdafiles/",
#                   paste(paste(format(Sys.time(),format="%Y%m%d"),
#                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                         "aldex2",output.filename,
#                         "workspace-effect.RData",sep = "-")))
```

Extract significant features.


``` r
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
```

Check for duplicates.


``` r
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
```

```
## # A tibble: 39 × 2
##    Taxon                             n
##    <chr>                         <int>
##  1 Muribaculaceae                    2
##  2 Bacteria Kingdom                  1
##  3 Bacteroidales Order               1
##  4 Bacteroides                       1
##  5 Bacteroidia Class                 1
##  6 Barnesiellaceae Family            1
##  7 Burkholderiales Order             1
##  8 Butyricimonas                     1
##  9 Christensenellaceae_R-7_group     1
## 10 Clostridia Class                  1
## # ℹ 29 more rows
```

``` r
aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)
head(aldex.neg.effect)
```

```
##                class  rab.all rab.win.0 rab.win.1  diff.btw diff.win    effect
## 1 covariatesMSMmouse 10.10326  10.27137 0.8207258  -9.18190 4.125691 -2.166613
## 2   covariatesrabbit 13.22247  13.37924 0.3962133 -13.19546 3.340937 -3.979079
##   effect.low effect.high      overlap                                 Taxon
## 1  -14.38756  -0.3588378 2.000325e-03 [Eubacterium]_coprostanoligenes_group
## 2  -25.74068  -0.9965875 2.005172e-05                        Muribaculaceae
```

Save the significant features with negative effect size as a 
 tab-separated file.


``` r
# "20260213_13_19_00"
# write.table(aldex.neg.effect,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),
#                                        sep = "_"),
#                                  "aldex.neg.effect",
#                                  output.filename,
#                                  "signif.tsv",sep="-")),
#             row.names = F,sep = "\t")
```

Save the workspace just in case.


``` r
# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   "aldex2",output.filename,"workspace.RData",sep="-")))
```

Save all significant features.


``` r
# "20260213_13_20_05"
# write.table(aldex.signif.features,
#             file=file.path("./output/rtables",authorname,paste(
#               paste(format(Sys.time(),format="%Y%m%d"),
#                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#               "aldex2",output.filename,"signif.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



### Run a test with ANCOMBC.


``` r
ancombc.out<-perform_ancombc_test(agglom.rank = agglom.rank,
                                  ps.q.df.wide = host.data.for_test$dataset,
                                  custom.md = host.data.for_test$metadata,
                                  sample.groups = host.data.for_test$sample.groups,
                                  ps.q.agg = host.data.for_test$ps.q.agg)
```

```
## Warning: Using `all_of()` outside of a selecting function was deprecated in tidyselect
## 1.2.0.
## ℹ See details at
##   <https://tidyselect.r-lib.org/reference/faq-selection-context.html>
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```
## 'ancombc' has been fully evolved to 'ancombc2'. 
## Explore the enhanced capabilities of our refined method!
```

```
## Checking the input data type ...
```

```
## The input data is of type: phyloseq
```

```
## PASS
```

```
## Checking the sample metadata ...
```

```
## The specified variables in the formula: class
```

```
## The available variables in the sample metadata: Sample, class, animal, sex, birthday
```

```
## PASS
```

```
## Checking other arguments ...
```

```
## The number of groups of interest is: 9
```

```
## The sample size per group is: NMR = 24, DMR = 20, B6mouse = 4, MSMmouse = 8, FVBNmouse = 3, spalax = 15, pvo = 10, hare = 8, rabbit = 7
```

```
## Warning: Small sample size detected for the following group(s): 
## B6mouse, FVBNmouse
## Variance estimation would be unstable when the sample size is < 5 per group
```

```
## PASS
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## Loading required package: foreach
```

```
## 
## Attaching package: 'foreach'
```

```
## The following objects are masked from 'package:purrr':
## 
##     accumulate, when
```

```
## Loading required package: rngtools
```

```
## ANCOM-BC primary results ...
```

```
## ANCOM-BC global test ...
```

```
## Merge the information of structural zeros ... 
## Note that taxa with structural zeros will have 0 p/q-values and SEs
```

``` r
ancombc.res<-ancombc.out$res
```

Find differentially abundant taxa by multiplying fold change with TRUE/FALSE.


``` r
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
```

```
##                       taxon_id (Intercept)   classDMR classB6mouse
## 14                   Alistipes    2.666286  2.0081199     3.233456
## 16                  Monoglobus    0.000000  3.3741910     2.376819
## 68  Erysipelotrichaceae Family    6.511852 -4.4226691    -4.220883
## 97                Anaerofustis    1.829179 -1.0532790    -1.408951
## 101             Muribaculaceae    7.701193  0.9109427     1.145490
## 102                 Prevotella    6.408340  1.3045003    -6.436052
##     classMSMmouse classFVBNmouse classspalax   classpvo  classhare classrabbit
## 14      1.8795204      2.7064949   0.9674499  2.1016352  3.6669313    2.207966
## 16      1.4576855      2.4729608   2.2073642  4.4374797  5.4301759    5.849992
## 68     -4.4673276     -3.0115540  -6.5440707 -6.5440121 -6.2345364   -4.078730
## 97     -1.2913571     -1.6260193  -1.8613977 -1.8613391 -1.8624767   -1.049186
## 101     0.5642364      0.9027479   1.9690254  0.4025274 -0.7299611   -7.334110
## 102    -6.4410614     -6.4362291  -6.0809833 -6.4404998 -6.4416374   -6.437341
```

``` r
head(ancombc.signif.decreased)
```

```
##                       taxon_id (Intercept)  classDMR classB6mouse classMSMmouse
## 68  Erysipelotrichaceae Family    6.511852 -4.422669    -4.220883     -4.467328
## 97                Anaerofustis    1.829179 -1.053279    -1.408951     -1.291357
## 108                   p-251-o5    6.185470 -3.420137    -5.866609     -6.218192
## 110                  Treponema    5.393237 -1.809481    -5.247662     -5.425959
## 111     Prevotellaceae_UCG-003    6.733609 -3.666130    -6.761321     -6.766331
## 123     Prevotellaceae_UCG-004    3.961322 -3.059358    -3.989034     -3.994043
##     classFVBNmouse classspalax  classpvo classhare classrabbit
## 68       -3.011554   -6.544071 -6.544012 -6.234536   -4.078730
## 97       -1.626019   -1.861398 -1.861339 -1.862477   -1.049186
## 108      -6.213360   -6.217689 -4.111775 -6.218768   -6.214471
## 110      -5.421126   -5.425456 -5.425397 -5.426535   -5.422238
## 111      -6.761498   -6.765828 -6.765769 -6.766907   -6.762610
## 123      -3.989211   -3.993540 -3.993482 -3.734689   -3.990323
```

Save the workspace.


``` r
# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   "ancombc",host.data.for_test$output.filename,"workspace.RData",sep="-")))
```

Save all significant features. 


``` r
# "20260213_13_41_58"
# write.table(ancombc.signif.features,
#             file=file.path("./output/rtables",authorname,paste(
#               paste(format(Sys.time(),format="%Y%m%d"),
#                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#               "ancombc",host.data.for_test$output.filename,
#               "signif.tsv",sep="-")), 
#             row.names = F,sep = "\t")
```

Save decreased significant features.


``` r
# "20260213_13_42_00"
# write.table(ancombc.signif.decreased,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),
#                                        sep = "_"),
#                                  "ancombc.signif.decreased",
#                                  host.data.for_test$output.filename,
#                                  "signif.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



### Downstream analysis and plotting of differentially abundant features. 


``` r
analyse_test_output(agglom.rank = agglom.rank,
                    inside.host = inside.host, 
                    comparison = comparison, 
                    custom.levels = host.data.for_test$custom.levels,
                    ref.level = ref.level, 
                    ps.q.agg = host.data.for_test$ps.q.agg,
                    sample.groups = host.data.for_test$sample.groups,
                    output.filename = host.data.for_test$output.filename,
                    custom.md = host.data.for_test$metadata,
                    maaslin.signif.features = maaslin.processed_output.host$maaslin.signif.features,
                    maaslin.signif.decreased = maaslin.processed_output.host$maaslin.signif.decreased,
                    aldex.signif.features = aldex.signif.features,
                    aldex.neg.effect = aldex.neg.effect,
                    ancombc.signif.features = ancombc.signif.features,
                    ancombc.signif.decreased = ancombc.signif.decreased)
```

```
## [1] "Muribaculaceae"    "Victivallis"       "UCG-002"          
## [4] "Prevotella"        "Gordonibacter"     "Turicibacter"     
## [7] "Ligilactobacillus" "Mycoplasma"        "Bacteroidia Class"
## character(0)
```

<img src="007-diffabund-tests_files/figure-html/unnamed-chunk-45-1.png" alt="" width="1056" />



## Compare NMR age groups (ASV level).
Choose what to compare:


``` r
comparison<-"age"
```

Choose the reference level:


``` r
ref.level<-"agegroup0_10" 
agglom.rank<-"OTU"
inside.host=TRUE
```

Prepare the data for testing.


``` r
nmr.age.data.for_test<-prepare_data(agglom.rank = agglom.rank, 
             comparison = comparison, 
             ref.level = ref.level, 
             inside.host = TRUE)
```

```
## [1] "Performing differential microbial abundance tests on age variable at OTU level. Reference: agegroup0_10"
## [1] "Output filename prefix: NMR-OTU-age-234-ref-agegroup0_10"
## [1] "agegroup0_10"  "agegroup10_16"
## [1] "Prepared the data for tests"
```

Run MaAsLin2 test.


``` r
maaslin.fit_data.age<-
  perform_maaslin2_test(ps.q.df.wide = nmr.age.data.for_test$dataset,
                        comparison = comparison,
                        output.filename = nmr.age.data.for_test$output.filename,
                        custom.md = nmr.age.data.for_test$metadata,
                        custom.levels = nmr.age.data.for_test$custom.levels,
                        ref.level = ref.level,
                        inside.host = inside.host
  )
```

```
## [1] "Creating output folder"
## [1] "Creating output feature tables folder"
## [1] "Creating output fits folder"
## [1] "Creating output figures folder"
## 2026-03-04 18:36:08.709761 INFO::Writing function arguments to log file
## 2026-03-04 18:36:08.757752 INFO::Verifying options selected are valid
## 2026-03-04 18:36:08.761578 INFO::Determining format of input files
## 2026-03-04 18:36:08.76495 INFO::Input format is data samples as rows and metadata samples as rows
## 2026-03-04 18:36:08.792263 INFO::Formula for random effects: expr ~ (1 | relation)
## 2026-03-04 18:36:08.796487 INFO::Formula for fixed effects: expr ~  agegroup
## 2026-03-04 18:36:08.802372 INFO::Filter data based on min abundance and min prevalence
## 2026-03-04 18:36:08.805966 INFO::Total samples in data: 24
## 2026-03-04 18:36:08.809735 INFO::Min samples required with min abundance for a feature not to be filtered: 0.000000
## 2026-03-04 18:36:08.83528 INFO::Total filtered features: 0
## 2026-03-04 18:36:08.840852 INFO::Filtered feature names from abundance and prevalence filtering:
## 2026-03-04 18:36:08.85969 INFO::Total filtered features with variance filtering: 0
## 2026-03-04 18:36:08.86524 INFO::Filtered feature names from variance filtering:
## 2026-03-04 18:36:08.870064 INFO::Running selected normalization method: TSS
## 2026-03-04 18:36:08.891325 INFO::Bypass z-score application to metadata
## 2026-03-04 18:36:08.896008 INFO::Running selected transform method: LOG
## 2026-03-04 18:36:08.916693 INFO::Running selected analysis method: LM
## 2026-03-04 18:36:08.92323 INFO::Fitting model to feature number 1, X775f1d61e616978239ee36c337b0403c
## 2026-03-04 18:36:08.984279 INFO::Fitting model to feature number 2, X5a818007a6452d5db9223ef1388e451c
## 2026-03-04 18:36:09.01923 INFO::Fitting model to feature number 3, X2a9b76ddc8fc1e5f0e10f1fd40de1c79
## 2026-03-04 18:36:09.051811 INFO::Fitting model to feature number 4, c513b4ef037af47125244161fb1eab50
## 2026-03-04 18:36:09.089317 INFO::Fitting model to feature number 5, X72aa784b6a3f05f45910fe33902caa81
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.12787 INFO::Fitting model to feature number 6, X37d6420f6843a99bbb892a3b37178db1
## 2026-03-04 18:36:09.160113 INFO::Fitting model to feature number 7, X5e6767d29a7bbc264a3edb4d4ef313e3
## 2026-03-04 18:36:09.190961 INFO::Fitting model to feature number 8, X087e97a98923593a6ca9d3292593baa8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.222014 INFO::Fitting model to feature number 9, a23076ded2ba0960eed8ce5e654dc804
## 2026-03-04 18:36:09.25692 INFO::Fitting model to feature number 10, X8607913ed1bbf9689f69dfccb73812bc
## 2026-03-04 18:36:09.286874 INFO::Fitting model to feature number 11, X6a9f0aa987f51dc6c85f04fb362f2257
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.318164 INFO::Fitting model to feature number 12, X050d316928bb2c2dcdc4b778c5d482af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.349569 INFO::Fitting model to feature number 13, d48f70a0a77c0821bd2dd0c91f355f30
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.382455 INFO::Fitting model to feature number 14, cbb7a2fe6431a3d874d38a8b6a104854
## 2026-03-04 18:36:09.41486 INFO::Fitting model to feature number 15, X3f058a90b6942f839656eb1dda7d9826
## 2026-03-04 18:36:09.453189 INFO::Fitting model to feature number 16, d3ba622809f71a4ea2725cffa1152398
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.484382 INFO::Fitting model to feature number 17, X4fd27ee04878ce466b07a8fdb6efb91e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.522944 INFO::Fitting model to feature number 18, X10023be92dba55239d6f4b32248afb74
## 2026-03-04 18:36:09.558026 INFO::Fitting model to feature number 19, X22b0f194732516c62929da18c03026fc
## 2026-03-04 18:36:09.589181 INFO::Fitting model to feature number 20, b096b6c93a5f390950dbf69ccedb31e5
## 2026-03-04 18:36:09.619772 INFO::Fitting model to feature number 21, X02e2fceb21147e8bf7ee6a805b54d187
## 2026-03-04 18:36:09.663719 INFO::Fitting model to feature number 22, X949c5cd1c3f0441a3332d7c892064517
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.702457 INFO::Fitting model to feature number 23, c41814121308544fe6e6dc6327604921
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.738728 INFO::Fitting model to feature number 24, d913c2cbde361858dec499be698ae593
## 2026-03-04 18:36:09.775302 INFO::Fitting model to feature number 25, X5d2a29b5344df06519f7d59f2041da54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.813236 INFO::Fitting model to feature number 26, X72d39f953d79e334e405cc903e866ef6
## 2026-03-04 18:36:09.844257 INFO::Fitting model to feature number 27, f6780cf12baa09329eeff07270ea2e08
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.884856 INFO::Fitting model to feature number 28, X63de5a8a3981bc0f00dbefe46ef47678
## 2026-03-04 18:36:09.91544 INFO::Fitting model to feature number 29, X9f2b8085820f945947d0ae63fb7ec339
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.946471 INFO::Fitting model to feature number 30, d886c00083ac9a3157fcb9dc989197b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:09.978513 INFO::Fitting model to feature number 31, X6034c86ff9a636506753fccb7cc73615
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.009693 INFO::Fitting model to feature number 32, X8ed8bcc26ac1f36c71117f0def9d0789
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.041984 INFO::Fitting model to feature number 33, X9426f864cbf3629f040d6c7c3ed64180
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.074256 INFO::Fitting model to feature number 34, X85ac797d67074ba62de3413d665354ba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.105237 INFO::Fitting model to feature number 35, X2e91b1556d45535a4b49b012917c5b19
## 2026-03-04 18:36:10.13622 INFO::Fitting model to feature number 36, dd4e02d8ebf6639915b74504fa6759a4
## 2026-03-04 18:36:10.167602 INFO::Fitting model to feature number 37, X0ca67ebe8a679e80efb2452016ef7125
## 2026-03-04 18:36:10.199247 INFO::Fitting model to feature number 38, e6289f64c0b74b4cfdf34cebc235b06f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.230686 INFO::Fitting model to feature number 39, dd87bb11986799f8870d596638f99f9e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.262329 INFO::Fitting model to feature number 40, X848bc2dee4d9f89037148afd929c3398
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.293561 INFO::Fitting model to feature number 41, X4940b5123681ac4251a9c017b731390e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.325394 INFO::Fitting model to feature number 42, X37677efe4c976fd9854262287add28ad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.355798 INFO::Fitting model to feature number 43, X398f3eae476d5c1bbb98af8d94ab3907
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.388766 INFO::Fitting model to feature number 44, X9d2b5cac004beb045c25df1c8a40b544
## 2026-03-04 18:36:10.419143 INFO::Fitting model to feature number 45, c38f378cecf99bde998dee1cd520fed4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.460939 INFO::Fitting model to feature number 46, X76a1f1eb3d45f64d33db506689665b59
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.492542 INFO::Fitting model to feature number 47, X390143a9a62dfa24c965f1c0843f9452
## 2026-03-04 18:36:10.523386 INFO::Fitting model to feature number 48, X35e23a83d3a6bba3e8af7ce15b881dc9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.554027 INFO::Fitting model to feature number 49, X22e069b421275afdc6b4e6cc0950ed4f
## 2026-03-04 18:36:10.583605 INFO::Fitting model to feature number 50, X501bbca65f0eef36751614dabbe08b36
## 2026-03-04 18:36:10.614147 INFO::Fitting model to feature number 51, X64f7efa6ab1e68ed7e0c98fd6ae1c915
## 2026-03-04 18:36:10.644231 INFO::Fitting model to feature number 52, X4bc9f240864b111cdf624408a87a6e1f
## 2026-03-04 18:36:10.674391 INFO::Fitting model to feature number 53, X665096a7e2a12378fcdaaacb3e58042e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.705335 INFO::Fitting model to feature number 54, X18e6ba724c5b286b604c2411d99de87a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.73621 INFO::Fitting model to feature number 55, X50d7799434925279751198338f77130c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.767306 INFO::Fitting model to feature number 56, cde69820d308b87b23fe6e239523e1a6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.79965 INFO::Fitting model to feature number 57, cbefcf8231ca919d37b62749a9881423
## 2026-03-04 18:36:10.830375 INFO::Fitting model to feature number 58, X2c1332de9241f56a6c7f8221758a9adc
## 2026-03-04 18:36:10.859896 INFO::Fitting model to feature number 59, X33b48213d11840fa28edd31fa038979b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.890067 INFO::Fitting model to feature number 60, b8ffd9b9fd2d561816c0c4ba72baaf85
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:10.920314 INFO::Fitting model to feature number 61, X6fba4056b223f50c3d51622fc1a84da0
## 2026-03-04 18:36:10.951912 INFO::Fitting model to feature number 62, dae1fc37c9776df8e228a443451de296
## 2026-03-04 18:36:10.989088 INFO::Fitting model to feature number 63, X6b9c74b7afbf96e071014516d8bc97c1
## 2026-03-04 18:36:11.01942 INFO::Fitting model to feature number 64, X9ee68e4ae7ff58ebfcfec4ee721842b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.050825 INFO::Fitting model to feature number 65, X5b956d9720139ff5385b0900d6a6bf28
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.08216 INFO::Fitting model to feature number 66, X311e9ec578b60769eaadf606a38a84cb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.11336 INFO::Fitting model to feature number 67, a58b408564f69ed083451dfacb552a90
## 2026-03-04 18:36:11.145116 INFO::Fitting model to feature number 68, a647fbcd189f56d189342f5ce294cd57
## 2026-03-04 18:36:11.176341 INFO::Fitting model to feature number 69, X3871a8798bd4cf074b1352e193c6f2e8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.216843 INFO::Fitting model to feature number 70, e7480012a7f751cd862916f0e6d468f7
## 2026-03-04 18:36:11.249182 INFO::Fitting model to feature number 71, a4d9c76f86b8ebf235f18ccbf6f67e86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.280366 INFO::Fitting model to feature number 72, abfe7842b6f39f0f8b0d7c1aebab8294
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.312073 INFO::Fitting model to feature number 73, X463ddfeb23025d3a52776c37581f0d73
## 2026-03-04 18:36:11.342447 INFO::Fitting model to feature number 74, X2ddfa5cd15c7edd114632d1ddbea899d
## 2026-03-04 18:36:11.374579 INFO::Fitting model to feature number 75, e55032c72052adb2f83ce0329ba984ce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.406231 INFO::Fitting model to feature number 76, b1433ade846d530f325945bb135f1bba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.436996 INFO::Fitting model to feature number 77, a6332efb97b7d387a65f5650829fb2e6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.469501 INFO::Fitting model to feature number 78, d1c2e6ccadfb83866f23fe93be5fada4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.500872 INFO::Fitting model to feature number 79, X50bc369ed002b3668d082f5965556786
## 2026-03-04 18:36:11.537472 INFO::Fitting model to feature number 80, b9c87bd9b886c6c2fc772e34369130b1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.572159 INFO::Fitting model to feature number 81, X8f0018f8ebf455bf3e078ccd7756f454
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.602754 INFO::Fitting model to feature number 82, X5bda8853404fe5cb2a95e1783ca59cd8
## 2026-03-04 18:36:11.633162 INFO::Fitting model to feature number 83, a7110f6b5ef234a8d21e3e5391faee2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.670531 INFO::Fitting model to feature number 84, X09a54c881aae0c9078c3cd72034eeddc
## 2026-03-04 18:36:11.700798 INFO::Fitting model to feature number 85, X867ea3db476d7118238f5b73d76c4841
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.734814 INFO::Fitting model to feature number 86, dd4b4f0ca9ebbbaab73466b57715192e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.767459 INFO::Fitting model to feature number 87, X05625778e3e287a8c0eb4e66e2e066f4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.798625 INFO::Fitting model to feature number 88, X08adbc5295e79ddbd27a1e49b05fc4dc
## 2026-03-04 18:36:11.829146 INFO::Fitting model to feature number 89, X4f71569f1364ca8c24d0971e533bf5f5
## 2026-03-04 18:36:11.861308 INFO::Fitting model to feature number 90, f83157419dd89248f996d2f278fbc713
## 2026-03-04 18:36:11.891667 INFO::Fitting model to feature number 91, X2341971441277c698da922ff13b6cf42
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.925363 INFO::Fitting model to feature number 92, X5cb934450030fdb4f808c8cca52961a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:11.955876 INFO::Fitting model to feature number 93, X679b1d8a517eb894be6da3eb3197ce90
## 2026-03-04 18:36:11.995103 INFO::Fitting model to feature number 94, X75e19d2211a08be3f623147f55e19aec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.03233 INFO::Fitting model to feature number 95, X10412032e5a4099296e96a115a57c7f6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.063053 INFO::Fitting model to feature number 96, af6b154cb33e86aa5495219241f05b22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.094069 INFO::Fitting model to feature number 97, X31947dbf83d9f2c7ae9c2d3d07ad669e
## 2026-03-04 18:36:12.125442 INFO::Fitting model to feature number 98, X85ca038ebdb000cf6701a43615349ea2
## 2026-03-04 18:36:12.159175 INFO::Fitting model to feature number 99, X233d661a3d8de3b4f9a28e926523d4c4
## 2026-03-04 18:36:12.189058 INFO::Fitting model to feature number 100, X8c6d9cefd98212bf81b189254d48b8f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.219358 INFO::Fitting model to feature number 101, X784386df30006e4485cc62dcf6743234
## 2026-03-04 18:36:12.250147 INFO::Fitting model to feature number 102, X503d066c2ba9bb8c8b2177cc6353501c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.281479 INFO::Fitting model to feature number 103, d4990c02b10ff6cfcb7e0f96c9196a49
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.311843 INFO::Fitting model to feature number 104, X5dff79a69be9a1b82db8a5638a1721a8
## 2026-03-04 18:36:12.342232 INFO::Fitting model to feature number 105, c1ae93de93f91093cffe77b5625360cb
## 2026-03-04 18:36:12.379801 INFO::Fitting model to feature number 106, bcfb05341d821b6242511abff6cf365e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.410849 INFO::Fitting model to feature number 107, X87a9a268b3ae6c7a3e20b3aa3ed33aa6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.443559 INFO::Fitting model to feature number 108, X36ad2071119d1950199ca71b1937b246
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.474153 INFO::Fitting model to feature number 109, X041f8052f648be811c7cc62cd2abf516
## 2026-03-04 18:36:12.51125 INFO::Fitting model to feature number 110, X164b12a7f94fed5ed43389908aab6c81
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.541146 INFO::Fitting model to feature number 111, X9c3d252c26b173da3d2afd11ffb176f8
## 2026-03-04 18:36:12.571471 INFO::Fitting model to feature number 112, X957b8fefd3f381946d85f43ee6086f74
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.601988 INFO::Fitting model to feature number 113, f13620b4f8abb891cc211b453c38244f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.632617 INFO::Fitting model to feature number 114, X38f6192f7df7849cbe5ce8e5b4aa3af0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.663208 INFO::Fitting model to feature number 115, c0c6de5f96a4bf636b60aff158542063
## 2026-03-04 18:36:12.699867 INFO::Fitting model to feature number 116, X98e5eedb473623066d9edab84e519c0f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.730999 INFO::Fitting model to feature number 117, X9b9f46c24a9200aa740f6c951155126d
## 2026-03-04 18:36:12.760744 INFO::Fitting model to feature number 118, X87531a4d94460fab40ad7b6c875690ca
## 2026-03-04 18:36:12.800266 INFO::Fitting model to feature number 119, e1e8a7389e298f5ed01997c426b8efbb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.830374 INFO::Fitting model to feature number 120, X9e260a0418dbb7f9ecf40402a3453429
## 2026-03-04 18:36:12.860441 INFO::Fitting model to feature number 121, X43d3e97de5be62fbc5bef7436779aef0
## 2026-03-04 18:36:12.890671 INFO::Fitting model to feature number 122, X4382fa89de91c2123c51b2c5a3526794
## 2026-03-04 18:36:12.921004 INFO::Fitting model to feature number 123, X7823dfa03b1068e544769d82d5108d3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.957826 INFO::Fitting model to feature number 124, X70f25a4988056da265638c4028e0e9b4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:12.988027 INFO::Fitting model to feature number 125, X8bdedac0d01bba7276354218f19a8703
## 2026-03-04 18:36:13.01867 INFO::Fitting model to feature number 126, f4a400e4fc836805ab7e7283605760d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.050097 INFO::Fitting model to feature number 127, b6b1ff24f46d79cb2b27464cddc4d695
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.080765 INFO::Fitting model to feature number 128, X9f4a338d0544b7ac1c7643304ac451bf
## 2026-03-04 18:36:13.111422 INFO::Fitting model to feature number 129, X546a95e7573470b5fdd2309be2684624
## 2026-03-04 18:36:13.142083 INFO::Fitting model to feature number 130, X68b62e646ef6c45a8d31a50e9b3f83e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.173116 INFO::Fitting model to feature number 131, c832b45fd6e005a61178908bb566ca82
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.203452 INFO::Fitting model to feature number 132, X2ba85532d09e46703511d81ccb3a4622
## 2026-03-04 18:36:13.233159 INFO::Fitting model to feature number 133, X7c992c13c3614197b1b331b6182bd3f5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.271044 INFO::Fitting model to feature number 134, X45cbf3ea056aff67a4a75d4e8da686c0
## 2026-03-04 18:36:13.302479 INFO::Fitting model to feature number 135, X6d60086ac416f97c8d3a774045920413
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.334145 INFO::Fitting model to feature number 136, X587d21cd26f3e0a14be44c7f15b7aaeb
## 2026-03-04 18:36:13.365413 INFO::Fitting model to feature number 137, adcbd92e759e2bd9c0bf798af532d8c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.397283 INFO::Fitting model to feature number 138, a30e54c24ebcf6bc5e5a1ea7b868934d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.429761 INFO::Fitting model to feature number 139, f06b506d343103b6969101eb88e864ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.467736 INFO::Fitting model to feature number 140, X69edbd69fd17306596e6f0765d8baaf5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.505016 INFO::Fitting model to feature number 141, e6d608aa96cffb876fdcb1c5e610454f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.540439 INFO::Fitting model to feature number 142, a1f9f5146e383650703e9ba1718c81ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.583814 INFO::Fitting model to feature number 143, X938b0c38f033d20f1b16c42c9fb10a36
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.622561 INFO::Fitting model to feature number 144, X0edc4ea81df8083fcd669e4d24f5c34c
## 2026-03-04 18:36:13.65275 INFO::Fitting model to feature number 145, bf4b19f3a06f13e371c47798b759f967
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.688889 INFO::Fitting model to feature number 146, X5a0cb8bb500ca6731dc4283cba2aaeaf
## 2026-03-04 18:36:13.72085 INFO::Fitting model to feature number 147, c77b5d1e3caa62de4030ae1460352956
## 2026-03-04 18:36:13.751611 INFO::Fitting model to feature number 148, f71741d8f1751a7511fbd4455a86593b
## 2026-03-04 18:36:13.782809 INFO::Fitting model to feature number 149, da3141bea3827e8b54dfc7403b6f3ab9
## 2026-03-04 18:36:13.813309 INFO::Fitting model to feature number 150, X3b86c743ad3ceff86644bc3ce0d66a16
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.848389 INFO::Fitting model to feature number 151, X7823c990b48e0d5fdbb67b0cf1511b0d
## 2026-03-04 18:36:13.879063 INFO::Fitting model to feature number 152, caae62e7d284909f0ba194e383c99196
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.911925 INFO::Fitting model to feature number 153, b2fcbcf8b5b2ee9ae672857fff1c4ae4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:13.944357 INFO::Fitting model to feature number 154, X002f3595b8a20738572f67dfcaf531a1
## 2026-03-04 18:36:13.975359 INFO::Fitting model to feature number 155, X02031267bf033fd8b5bc46043628acd1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.007671 INFO::Fitting model to feature number 156, X20ab8b9e32bd6d9bfc506b65fd696ef9
## 2026-03-04 18:36:14.038403 INFO::Fitting model to feature number 157, X587b64913a785e4d7cc810c80b2a9722
## 2026-03-04 18:36:14.075814 INFO::Fitting model to feature number 158, f41431a72bfd93b6f2c4a49ea7642ae8
## 2026-03-04 18:36:14.111598 INFO::Fitting model to feature number 159, X8ab93fef97076ce48913426901d23711
## 2026-03-04 18:36:14.143607 INFO::Fitting model to feature number 160, X8dca4bfdcbcbb43dcbbf5a6c7943b833
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.176336 INFO::Fitting model to feature number 161, X7fb1c8722f2ebb009c074460b08d5b5f
## 2026-03-04 18:36:14.208872 INFO::Fitting model to feature number 162, X27b102c4d0f321ccb73b4fdd3cf04efa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.240272 INFO::Fitting model to feature number 163, dd3c9baa3973a152a2bada71543454f5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.272085 INFO::Fitting model to feature number 164, X5ae982af003eaae773f2ac4ac8d1f518
## 2026-03-04 18:36:14.309708 INFO::Fitting model to feature number 165, X52f10c839bad102da0f9182fe2a65d63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.344357 INFO::Fitting model to feature number 166, X4cf488fc55937cceb0cf0d50c67f78c7
## 2026-03-04 18:36:14.387292 INFO::Fitting model to feature number 167, X0e7e09e2f26a49f7e392c9a979f935d0
## 2026-03-04 18:36:14.418736 INFO::Fitting model to feature number 168, b9ce57fad79d507e7f193ab62a1d590f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.450813 INFO::Fitting model to feature number 169, X1f1f739ae6bb988057eb2f9fd17f3131
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.483619 INFO::Fitting model to feature number 170, fc47c350e68145c3f6a5dc0a818fb683
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.515192 INFO::Fitting model to feature number 171, X3fb50d74ec763a9a895b6a69751056ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.546313 INFO::Fitting model to feature number 172, eee7a5d0e69d4129a7cd4fd8c7a6abea
## 2026-03-04 18:36:14.583703 INFO::Fitting model to feature number 173, c9e692a4f678ce3e9a9bf69290382934
## 2026-03-04 18:36:14.620231 INFO::Fitting model to feature number 174, a8d690ae47779a14d63923f302f23b44
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.657956 INFO::Fitting model to feature number 175, X4989b6aa56952d17a9728cb74ba4eccd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.688583 INFO::Fitting model to feature number 176, ac5128abccedf13cded0c2d74b428936
## 2026-03-04 18:36:14.719407 INFO::Fitting model to feature number 177, bf5cec3c7cd23cf88fedb590337edeb7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.750277 INFO::Fitting model to feature number 178, X196dd825ebed8d5e0372081332465e46
## 2026-03-04 18:36:14.780785 INFO::Fitting model to feature number 179, X1664dd935ff0a38d0bc5b804f2c5f6aa
## 2026-03-04 18:36:14.810272 INFO::Fitting model to feature number 180, a32ebbe366e9cd264ffb16a15d10dcb9
## 2026-03-04 18:36:14.839748 INFO::Fitting model to feature number 181, X30803d11bc6be93d144583cd061c1434
## 2026-03-04 18:36:14.869069 INFO::Fitting model to feature number 182, X534c9f3c39ff74b0751677a74ca2223f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.899412 INFO::Fitting model to feature number 183, fcda29c5385d6115d10b5e976e087ae1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:14.929606 INFO::Fitting model to feature number 184, X8a238b6279fbb51f2a43b9c11a1e7293
## 2026-03-04 18:36:14.960146 INFO::Fitting model to feature number 185, X5e489e67242303169c8781859aa2ea64
## 2026-03-04 18:36:14.990172 INFO::Fitting model to feature number 186, X7109c3732f21db3d8ca19a4be2231a80
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.020421 INFO::Fitting model to feature number 187, X7af0ca6b149c28ed1da4c04b112c4551
## 2026-03-04 18:36:15.051634 INFO::Fitting model to feature number 188, X4404bd2d87786e911fba13e90fd51988
## 2026-03-04 18:36:15.081792 INFO::Fitting model to feature number 189, X4ce796eca0eb0252c01451f0acfd7c86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.11388 INFO::Fitting model to feature number 190, ff162ce0992d8e762f0f09422814ed01
## 2026-03-04 18:36:15.16064 INFO::Fitting model to feature number 191, X0c0189fefd69f40756d340e48d5f7304
## 2026-03-04 18:36:15.193674 INFO::Fitting model to feature number 192, X5dfc13e1239763278b4edaf53ab409b5
## 2026-03-04 18:36:15.226556 INFO::Fitting model to feature number 193, X9c9d04df1b23dcd81f5bbdb09d422af4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.258911 INFO::Fitting model to feature number 194, X8585426a699ce2bd5bdd80954f4df256
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.291092 INFO::Fitting model to feature number 195, f1ca62a3e8e5e041350db2c975346cb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.323521 INFO::Fitting model to feature number 196, X59f4de3950a671dd935530636121f79a
## 2026-03-04 18:36:15.354244 INFO::Fitting model to feature number 197, X8dee14385c53c6adc0b6ed841329b318
## 2026-03-04 18:36:15.386082 INFO::Fitting model to feature number 198, fb6e884f98ccf551138d94557d49b6c4
## 2026-03-04 18:36:15.424057 INFO::Fitting model to feature number 199, d749b2d647bb124c92604b1bdb7c6511
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.465209 INFO::Fitting model to feature number 200, X318bbf7176eb37d514c95f727b120e3f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.504326 INFO::Fitting model to feature number 201, X53f63e6d0ded3831371b61c05e0d1df9
## 2026-03-04 18:36:15.536571 INFO::Fitting model to feature number 202, X810c56e7a79638b2fdf56ece16b9a32c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.567934 INFO::Fitting model to feature number 203, e3ca497d1ddf9a2efda5736924232efe
## 2026-03-04 18:36:15.599985 INFO::Fitting model to feature number 204, X5453c211a413666e9c72d8a12bd3bf2f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.637847 INFO::Fitting model to feature number 205, X6dbf7882c3e44b3f93cd770a5af49299
## 2026-03-04 18:36:15.672593 INFO::Fitting model to feature number 206, X62e6d4459d864709f6413926e3a88095
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.707643 INFO::Fitting model to feature number 207, b4a791cce24c1baaf3966a1d9efd18bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.742088 INFO::Fitting model to feature number 208, X4c944499b014eb72118a65700a136719
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.775793 INFO::Fitting model to feature number 209, X278a0e8506254a8e6bd3ad6e37df2e1d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.808033 INFO::Fitting model to feature number 210, X406caa0d474693f37020bc3b3e99467d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.838712 INFO::Fitting model to feature number 211, X580f23c911e277aeb396a2dd45af2675
## 2026-03-04 18:36:15.874545 INFO::Fitting model to feature number 212, X474ee1999580bb7f367858c97607d986
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.907049 INFO::Fitting model to feature number 213, X72ef153369491e40b5b0405e442bb1d4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.938348 INFO::Fitting model to feature number 214, X42e99751e005f348bbb3a852815e7dab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:15.983625 INFO::Fitting model to feature number 215, X796e01ba695bbc78aaa63196bcd06330
## 2026-03-04 18:36:16.014382 INFO::Fitting model to feature number 216, b1258e4ee5ca06eb4d5b62a466e413ff
## 2026-03-04 18:36:16.045803 INFO::Fitting model to feature number 217, c0248093940670814ffe413551c08592
## 2026-03-04 18:36:16.081797 INFO::Fitting model to feature number 218, cb5947bc03ed14060ebce14a528a30a9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.11238 INFO::Fitting model to feature number 219, X7c698019fe98da258752a851be4f0dbf
## 2026-03-04 18:36:16.143782 INFO::Fitting model to feature number 220, f819c2a6d74d9e73b54d2b2d782b42ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.181901 INFO::Fitting model to feature number 221, X91aaea5b58ac19ee8a8bde60f2368289
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.218715 INFO::Fitting model to feature number 222, X7d3ca27433dbfec1d46096d6f472c0e6
## 2026-03-04 18:36:16.256185 INFO::Fitting model to feature number 223, X46bb1d815c069f361d1832faf2e2e4dd
## 2026-03-04 18:36:16.293081 INFO::Fitting model to feature number 224, ae38f6e010625cb3b1b6cea1eea13a57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.330866 INFO::Fitting model to feature number 225, X6bca7aa2b33b0274c340d4facd92afaa
## 2026-03-04 18:36:16.367016 INFO::Fitting model to feature number 226, f430ec182a2fba6dffb1534634aba940
## 2026-03-04 18:36:16.39912 INFO::Fitting model to feature number 227, X9a67c9da63556d4863bab09cafd40940
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.430933 INFO::Fitting model to feature number 228, e92df5b83e69c4ed9bb9d268ecb358e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.462701 INFO::Fitting model to feature number 229, X1fa6739ed5fc2a89235fb28d33d3ccb3
## 2026-03-04 18:36:16.493935 INFO::Fitting model to feature number 230, X1d93a791c26759878a10132d5e6e33fa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.525353 INFO::Fitting model to feature number 231, c67ed17a9beed3c005cf97f0e0f20a18
## 2026-03-04 18:36:16.562571 INFO::Fitting model to feature number 232, db3aca72ce7e5948035457f7871945a2
## 2026-03-04 18:36:16.594291 INFO::Fitting model to feature number 233, a769521df2160b6d284dfd3e1d449c9d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.627237 INFO::Fitting model to feature number 234, X7a21b6c3fa63544dd1cde940adf52388
## 2026-03-04 18:36:16.660904 INFO::Fitting model to feature number 235, X882e40e7940f3b6e5c6614a047e97dc6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.693029 INFO::Fitting model to feature number 236, b9846f46f48c55df07601aaef88e17cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.727995 INFO::Fitting model to feature number 237, X920d9d9fa0d42f62df13a5f61eba6019
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.759152 INFO::Fitting model to feature number 238, bad8cf42670e052284875915a6e94c58
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.806696 INFO::Fitting model to feature number 239, X27b1308239803f1111e6871b0f7733d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.83913 INFO::Fitting model to feature number 240, fc4c1b4d08d17bc8bbbfeccfc267f51a
## 2026-03-04 18:36:16.873418 INFO::Fitting model to feature number 241, X54af66f5c264f0024789d99c4c7afe8d
## 2026-03-04 18:36:16.904043 INFO::Fitting model to feature number 242, X5e4b3686f1cefe645291fca5b936ceb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.938543 INFO::Fitting model to feature number 243, X82df3e7e293f65085e273db27eac0658
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:16.969931 INFO::Fitting model to feature number 244, X2419595d114a6b8b454c88b03531471e
## 2026-03-04 18:36:17.008889 INFO::Fitting model to feature number 245, X8bcc8165e984ed5f75e38224a26c4500
## 2026-03-04 18:36:17.046607 INFO::Fitting model to feature number 246, fed6e1a07e0a36780b2112ad59790baf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.080584 INFO::Fitting model to feature number 247, X3a81c2a8bcd8262593a6a28835b1d5d1
## 2026-03-04 18:36:17.112404 INFO::Fitting model to feature number 248, X07aa13125fef95029e4e578ff8d77194
## 2026-03-04 18:36:17.143323 INFO::Fitting model to feature number 249, X2ee55325dd6fff97f7e798daf0a1a173
## 2026-03-04 18:36:17.174218 INFO::Fitting model to feature number 250, cdda58bad36d42fffe349c500561263a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.21122 INFO::Fitting model to feature number 251, d8bc5ba9f36585e295c58310928b7610
## 2026-03-04 18:36:17.242597 INFO::Fitting model to feature number 252, X2443c7149b3c68f110a91c9464cd3a9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.275348 INFO::Fitting model to feature number 253, d3ce175f14e9bd4ce6fa566098070432
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.30711 INFO::Fitting model to feature number 254, f5852b6815968593b30bfc519b224b13
## 2026-03-04 18:36:17.340639 INFO::Fitting model to feature number 255, X96cec0b11138b0ddb91bc8fd3a827e63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.377021 INFO::Fitting model to feature number 256, faecb8ca1ca70b02713116054ba6a25e
## 2026-03-04 18:36:17.413508 INFO::Fitting model to feature number 257, X428241385b3872803e4e6f09856f1eaf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.445151 INFO::Fitting model to feature number 258, f7def96d9ec34bbc94788d5096e6b653
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.475842 INFO::Fitting model to feature number 259, X5bdcbe959f1cc3b608874a0899b2b431
## 2026-03-04 18:36:17.506917 INFO::Fitting model to feature number 260, X5df10ae748068240d0954e709aa85ef6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.538142 INFO::Fitting model to feature number 261, X773a6686911fcf580df6cc6ce27e25a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:17.569115 INFO::Fitting model to feature number 262, X962f048f3548d862a1fa9bfef48d3c8d
## 2026-03-04 18:36:18.073198 INFO::Fitting model to feature number 263, X10d3555fce8c1645f69075680f9644e6
## 2026-03-04 18:36:18.103974 INFO::Fitting model to feature number 264, X5eb1bd896aeef657980bb6319e0fa46c
## 2026-03-04 18:36:18.13469 INFO::Fitting model to feature number 265, a1163d29062c052bb3d59c49fdc0c44a
## 2026-03-04 18:36:18.172885 INFO::Fitting model to feature number 266, X56ee0e676dc712ca9636f01f082d12e1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.210527 INFO::Fitting model to feature number 267, e2773b5fd096352ad027595baeab01ef
## 2026-03-04 18:36:18.24767 INFO::Fitting model to feature number 268, X66eecd48e1caa9adf86f0ef816e22d87
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.283463 INFO::Fitting model to feature number 269, cfc06d9223155de233be7daf646b0cfe
## 2026-03-04 18:36:18.314162 INFO::Fitting model to feature number 270, X741671cead50e930b5ae2fdba6f09ac2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.348425 INFO::Fitting model to feature number 271, X85124c575561a28bdd0e3985153d33eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.38002 INFO::Fitting model to feature number 272, X8c0a48a6da32d89300557a08eb584eab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.417227 INFO::Fitting model to feature number 273, X3de3d0c798ac11c7eb2d3cd2c5eb6d22
## 2026-03-04 18:36:18.456707 INFO::Fitting model to feature number 274, fdc6ffd6aff43d49592d12868ef6edda
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.490248 INFO::Fitting model to feature number 275, X175e258fe4d2966115d9f5607a6dc868
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.520992 INFO::Fitting model to feature number 276, X2adacb50dbacd8fa223cd5107f430319
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.552645 INFO::Fitting model to feature number 277, X4c09b59637d648bb9d4cead257187763
## 2026-03-04 18:36:18.583313 INFO::Fitting model to feature number 278, X6314d0517e01bcc3d8907f875891bddd
## 2026-03-04 18:36:18.613702 INFO::Fitting model to feature number 279, X001ce79a8f8d035a30c8203a7d76b05e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.645164 INFO::Fitting model to feature number 280, X30b0cc0ca2a2da9533100ba3491175d6
## 2026-03-04 18:36:18.676387 INFO::Fitting model to feature number 281, d0750756e3311a835281cad227e6c4a0
## 2026-03-04 18:36:18.708125 INFO::Fitting model to feature number 282, e72dbc2cbbcaf8d32bf52077e98fb9d6
## 2026-03-04 18:36:18.739551 INFO::Fitting model to feature number 283, X0fdd81e04676039a9b2bcbcf6ac166ff
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.771924 INFO::Fitting model to feature number 284, X6d729250b3b458887e3362f17fb6b54b
## 2026-03-04 18:36:18.802275 INFO::Fitting model to feature number 285, X58dfc0bd94f7612268c6e40fa4819fc5
## 2026-03-04 18:36:18.834598 INFO::Fitting model to feature number 286, X62111f7adcf73e318470d3c36ac127a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.873087 INFO::Fitting model to feature number 287, a098856e088a414969785c312c7077dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.912965 INFO::Fitting model to feature number 288, f7db7195338b7afd38f1f529a1d6a8c2
## 2026-03-04 18:36:18.942435 INFO::Fitting model to feature number 289, X5a101d21f512bd402ba5dce9ff8853b3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:18.97676 INFO::Fitting model to feature number 290, X722c3b22b71c06ae914fdc65c8a93e1b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.007295 INFO::Fitting model to feature number 291, X28d37664dfaf5d0318674d3903759319
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.038481 INFO::Fitting model to feature number 292, d9d2694a2422852ddb9e381500af7cca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.070507 INFO::Fitting model to feature number 293, X3d9a5be20ced9a47b9cdba62499e9b67
## 2026-03-04 18:36:19.10164 INFO::Fitting model to feature number 294, X49cddebff224d6008280a502f630d0dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.133324 INFO::Fitting model to feature number 295, X8adece3b43b1d67576ec24a05663f526
## 2026-03-04 18:36:19.164431 INFO::Fitting model to feature number 296, X2f748f1a2371911911065af675ff0818
## 2026-03-04 18:36:19.194854 INFO::Fitting model to feature number 297, b985a8c66261968b9d8a329d29755c8c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.232229 INFO::Fitting model to feature number 298, c77e2df1b5831857b410fd099900e05a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.263502 INFO::Fitting model to feature number 299, b8c05b8dc05f4f544735add1f5f14ee9
## 2026-03-04 18:36:19.294454 INFO::Fitting model to feature number 300, e5aa8ee9cb7d24ef5df26f0bc2a72299
## 2026-03-04 18:36:19.325565 INFO::Fitting model to feature number 301, c47ab297ca5c2e467fb15ca968fbf667
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.356832 INFO::Fitting model to feature number 302, X48da3e64315f06b66b46c4a7caa33b26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.388 INFO::Fitting model to feature number 303, X58a318efa7e9171f7b768500c10b0a85
## 2026-03-04 18:36:19.418676 INFO::Fitting model to feature number 304, X680fcfeaf91f648d224c23cb0c49e41a
## 2026-03-04 18:36:19.448939 INFO::Fitting model to feature number 305, X2d3e451e9169d654b061d02d71d3b319
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.479858 INFO::Fitting model to feature number 306, X7df3c6b0b97f9d38b18f47462e2eb666
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.511115 INFO::Fitting model to feature number 307, fb9c18e59bc4a13d9c4b1830cc979ba8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.541534 INFO::Fitting model to feature number 308, X6bae82c66563cfba208caab642e21b9c
## 2026-03-04 18:36:19.573718 INFO::Fitting model to feature number 309, feaf3412296cd21fe71e288500f2909d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.606178 INFO::Fitting model to feature number 310, X2daaed7767de192f02b5ef6aaac2fd45
## 2026-03-04 18:36:19.637988 INFO::Fitting model to feature number 311, X88f2e5b2a5d102be5d9d22e56c5b5033
## 2026-03-04 18:36:19.670417 INFO::Fitting model to feature number 312, ab0a47c0d634005b50d204cc87f0300f
## 2026-03-04 18:36:19.710035 INFO::Fitting model to feature number 313, X2b8c21681ff428a1c834dd9b1fd35d4f
## 2026-03-04 18:36:19.740474 INFO::Fitting model to feature number 314, X590139931aa19c8ed5551ec58c8fc556
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.775986 INFO::Fitting model to feature number 315, X1958f1ea043fda09b0846074b8cf3387
## 2026-03-04 18:36:19.806332 INFO::Fitting model to feature number 316, bcfaaaae2ef79153b03ef0814d0d7abf
## 2026-03-04 18:36:19.843947 INFO::Fitting model to feature number 317, ada302e958702d87a162ac13e5225fd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.874791 INFO::Fitting model to feature number 318, d733da4c4eb11b498ccf79ced6da79d2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.906858 INFO::Fitting model to feature number 319, X1142cc5fdfb17bf3dda55341b1fc2ece
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.94464 INFO::Fitting model to feature number 320, X16a4be5bfb2132b4cf9b1ce6753abbf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:19.975452 INFO::Fitting model to feature number 321, X1daf20c298f03dfd4551813455d60a75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.005213 INFO::Fitting model to feature number 322, X5d982d4a4a7f3ce92e0e0dbfe0008be4
## 2026-03-04 18:36:20.035256 INFO::Fitting model to feature number 323, X6f3b3e805737133d5a15e6150096b7ad
## 2026-03-04 18:36:20.065698 INFO::Fitting model to feature number 324, a63f57972cf38df045ae50390d88758b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.103039 INFO::Fitting model to feature number 325, X21df1e53f11d76242c275d0844a43fdf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.138493 INFO::Fitting model to feature number 326, X7880d3e8bcca1feb2db8b9839741438f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.170709 INFO::Fitting model to feature number 327, ebd2eb0273b048077d97a8e5c87843be
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.201275 INFO::Fitting model to feature number 328, X17c62c1c91771a6d6a11af56744a4812
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.237401 INFO::Fitting model to feature number 329, X4e2b574c5c896013c886dcccfe212bb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.267279 INFO::Fitting model to feature number 330, c351ea5b351f3da50f771e157e96c96f
## 2026-03-04 18:36:20.297789 INFO::Fitting model to feature number 331, X59a9c9743d0de49303aa710f7f2dcf1a
## 2026-03-04 18:36:20.328298 INFO::Fitting model to feature number 332, bdb97410884f71c5be50c8881e2c98b0
## 2026-03-04 18:36:20.35833 INFO::Fitting model to feature number 333, a36840adf6ca830c32cba3a0bf05c24c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.395061 INFO::Fitting model to feature number 334, da20e0254113aed621b9f35dd3b17a74
## 2026-03-04 18:36:20.424628 INFO::Fitting model to feature number 335, X4cf353124acf76a83451339182554a2c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.455734 INFO::Fitting model to feature number 336, d699ad5c7d262166690df8e94610d11d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.487215 INFO::Fitting model to feature number 337, X0e95449d477a09d2ab60aad4e716e758
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.527096 INFO::Fitting model to feature number 338, c7a4a2780d280e743c16f5be29834ceb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.557399 INFO::Fitting model to feature number 339, d26a840ed062289b08114e43ed383a08
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.587198 INFO::Fitting model to feature number 340, d89f4af0cc9704172e5a75524368626a
## 2026-03-04 18:36:20.617093 INFO::Fitting model to feature number 341, cd007c70afbc5e207cbd6979ffc64409
## 2026-03-04 18:36:20.646634 INFO::Fitting model to feature number 342, d0050df77ec59e70e009b5af346e275f
## 2026-03-04 18:36:20.67638 INFO::Fitting model to feature number 343, f148106f9939021fa196057f3c53f008
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.706611 INFO::Fitting model to feature number 344, X6b8516936fd7f3084db293f299865b97
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.737118 INFO::Fitting model to feature number 345, a3cf3667e4778d77c70efcf821731124
## 2026-03-04 18:36:20.771244 INFO::Fitting model to feature number 346, X268786e064d71003f5ae04c6bc980978
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.80772 INFO::Fitting model to feature number 347, X3e398419b48fe7fd6882bfc636293558
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.838551 INFO::Fitting model to feature number 348, b5261b3c77cd990eef2c90209b8fd44f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.869216 INFO::Fitting model to feature number 349, cd892c02a21141c715fdef5c234c791c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.901287 INFO::Fitting model to feature number 350, d417a886163ff1310f4ff6be7122a500
## 2026-03-04 18:36:20.931948 INFO::Fitting model to feature number 351, faa27d3e3fc59a116b955ca7bad4567d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:20.963148 INFO::Fitting model to feature number 352, X30b59313f4ae8d27dd90ead569dd220c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.000557 INFO::Fitting model to feature number 353, X9fe431eb75ff0fb76d1b3597e7c6e018
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.030812 INFO::Fitting model to feature number 354, c5330203390f5f92c8631531cd1929dd
## 2026-03-04 18:36:21.068323 INFO::Fitting model to feature number 355, X37ab667cfa248f97195af5b3622bfc11
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.10252 INFO::Fitting model to feature number 356, X480c1786c18bdbaa003293d15e5aaf50
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.133098 INFO::Fitting model to feature number 357, ebfa013a17853303ce8dcb5b81825b1f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.163885 INFO::Fitting model to feature number 358, fee34b2a495eca65a0f5a95ba6f8348c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.195257 INFO::Fitting model to feature number 359, X44333eb795d84f9c6d46999cc2e9cbb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.22758 INFO::Fitting model to feature number 360, X61f6666d6fb9564cd5c8e249386bb508
## 2026-03-04 18:36:21.25928 INFO::Fitting model to feature number 361, ebfc68839dcbe864af2c9e957124a850
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.299278 INFO::Fitting model to feature number 362, X44f394b4abd8f4a0b77f5d918b2f1a3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.337071 INFO::Fitting model to feature number 363, e770697660b66ac58e6264e25fe09caf
## 2026-03-04 18:36:21.367162 INFO::Fitting model to feature number 364, X3a08076d3ca27f6a4e23b428c227f723
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.398434 INFO::Fitting model to feature number 365, X7d9f76c17978c3357bc6e5a39c3e5c15
## 2026-03-04 18:36:21.428972 INFO::Fitting model to feature number 366, c74b3b8f79a2710cd2b055e35178c301
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.462655 INFO::Fitting model to feature number 367, d67eb9feb3c3e913b0eda1a3ccfd6876
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.492927 INFO::Fitting model to feature number 368, e57ea9765333a478b6d4a150a0f4759c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.524787 INFO::Fitting model to feature number 369, X008d766cc2fe0e9a1ff0df6a107b02af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.562371 INFO::Fitting model to feature number 370, X17175a3de507f921278c078d71de6978
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.594663 INFO::Fitting model to feature number 371, X42db06a32d44f7e0b04833ac260cd01f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.631497 INFO::Fitting model to feature number 372, X53055f08d7d20cfe83c3e1ba05b092d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.662034 INFO::Fitting model to feature number 373, X9daace87b866febfba2396dea07f2e3b
## 2026-03-04 18:36:21.695433 INFO::Fitting model to feature number 374, X4e159cc1d60cb0bdf1d6504d2446411b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.727031 INFO::Fitting model to feature number 375, X6ceba268b6002d018f5be12e1d30ca48
## 2026-03-04 18:36:21.757873 INFO::Fitting model to feature number 376, X917d19bb0a01a3511a262d04c27211fd
## 2026-03-04 18:36:21.789823 INFO::Fitting model to feature number 377, c7ef7fb22153326e3ebb23932d0c8f8c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.820006 INFO::Fitting model to feature number 378, X5664cc65cd15a2a34cae8203f26b6a1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.85049 INFO::Fitting model to feature number 379, b11e31f8a98617d9b537f565aa9d5cd0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.881144 INFO::Fitting model to feature number 380, X0e64362e5d9e8e0393a2a25ef51344d8
## 2026-03-04 18:36:21.910849 INFO::Fitting model to feature number 381, X3f6a3469b81ee3684014d0fd4ca8ba22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.941237 INFO::Fitting model to feature number 382, X84ded3c82c02986d17ba89efc005ca1b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:21.971514 INFO::Fitting model to feature number 383, X7128c784b9f14437017f065ccd13d93e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.002667 INFO::Fitting model to feature number 384, a2f1c3dabb5720aa985c1c56c785a79e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.038557 INFO::Fitting model to feature number 385, c3a3131987eb1c89019993c61051a9b7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.069621 INFO::Fitting model to feature number 386, da69dbfd52386c4e8bec85dcf31fab91
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.116323 INFO::Fitting model to feature number 387, X00bec85d4cac8001d056f3641e4fe4d9
## 2026-03-04 18:36:22.147788 INFO::Fitting model to feature number 388, X2c756f3da2bf4aa3ba4cb108bfe02ef2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.18455 INFO::Fitting model to feature number 389, d564fadbbec41c3174402b6ed7c2b4e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.222254 INFO::Fitting model to feature number 390, eb4dffb13b2c5f8569c37ff99084fede
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.259671 INFO::Fitting model to feature number 391, ebe44d95e0c653b8424e7f4522c4bff3
## 2026-03-04 18:36:22.296269 INFO::Fitting model to feature number 392, fa81c053b819a9eeadc0f0306c07528b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.334278 INFO::Fitting model to feature number 393, X0538a2fd5f6431b8728a0e6770ffa2df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.365948 INFO::Fitting model to feature number 394, a32076c3bf3b662e045a8b42816f37bf
## 2026-03-04 18:36:22.398035 INFO::Fitting model to feature number 395, dba2b43ba071ef5f2ab8490e283e56e4
## 2026-03-04 18:36:22.429175 INFO::Fitting model to feature number 396, X8f7b886fa9032c9830c8a4232c454ab5
## 2026-03-04 18:36:22.460425 INFO::Fitting model to feature number 397, aa45112add34677769b6912d1749b794
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.493042 INFO::Fitting model to feature number 398, f070e6f921a4473b49d198509117a125
## 2026-03-04 18:36:22.523351 INFO::Fitting model to feature number 399, X145a90b7279b2987c8f0eb9c4edc9ea6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.553842 INFO::Fitting model to feature number 400, X4e04424d4099af6d663496c7a432bfc3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.58497 INFO::Fitting model to feature number 401, bcfc621313b4ba3bafd334ad6f2a10f5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.615367 INFO::Fitting model to feature number 402, X3d8b735fffa6bb77e92fe95a8042fa00
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.645597 INFO::Fitting model to feature number 403, X5ad134b39170fe286ea186ba2604897e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.677146 INFO::Fitting model to feature number 404, c41bf5e56e53e0e28470374de88e6ba7
## 2026-03-04 18:36:22.708594 INFO::Fitting model to feature number 405, X1ecd1c27ef8548472bc31b8d51d8c4dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.743696 INFO::Fitting model to feature number 406, X7120bc0c2afde34437f0970446441b39
## 2026-03-04 18:36:22.777154 INFO::Fitting model to feature number 407, X78c42a2db92f4e3f05634c101cc11259
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.809769 INFO::Fitting model to feature number 408, e2041d00358f2c29a5ef9809e8b1ba96
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.841672 INFO::Fitting model to feature number 409, X6d6787eb066d8db71d9e3631225a8904
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.871947 INFO::Fitting model to feature number 410, X6fa8ac18801d6dfe5d45d7336d9f5887
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.903114 INFO::Fitting model to feature number 411, X76e9b07ccc3729985a554bc91e58bbc0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:22.943857 INFO::Fitting model to feature number 412, e188b87dc228e850ec7535fba9be379f
## 2026-03-04 18:36:22.977342 INFO::Fitting model to feature number 413, X956a5b207bb050c7ce43853ce336d5e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.008179 INFO::Fitting model to feature number 414, X174c5b603cffda455b1a9969d47dbaae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.046144 INFO::Fitting model to feature number 415, X27a528bc7a1cdd2669ae30690d420170
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.077109 INFO::Fitting model to feature number 416, X27ef3adf14f05f2c77b6c64118cddb51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.107773 INFO::Fitting model to feature number 417, X7a92bc426ff0b6b2428d8b843bca7b6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.138859 INFO::Fitting model to feature number 418, efa9ff2bcd8d4e435d7982e3e317de43
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.169191 INFO::Fitting model to feature number 419, X39be2b3eb5948a1b4c0e9f823eba63d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.204736 INFO::Fitting model to feature number 420, X726ad819f9c01b960e30637dd15b78c0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.240534 INFO::Fitting model to feature number 421, X9286b1ca453940aa45d397ab3e2058f1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.274608 INFO::Fitting model to feature number 422, X9cc11a19e2926c7037eff26c869725fc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.305612 INFO::Fitting model to feature number 423, a29e744251cd4a2f29279c6aa377d8da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.336254 INFO::Fitting model to feature number 424, b511ce605bde5b3c893d53d043c00b2f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.367675 INFO::Fitting model to feature number 425, c16b00d5bab290920ae7d2d0affe4ca8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.402811 INFO::Fitting model to feature number 426, X1ab6690e22c490e00ad8c2ac9bf24596
## 2026-03-04 18:36:23.435435 INFO::Fitting model to feature number 427, cadb8b4f0b3bfaa180c90a78357042c9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.46855 INFO::Fitting model to feature number 428, d6c5071ba31869771d7462efa2058e02
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.503145 INFO::Fitting model to feature number 429, e9746dcc794f4362e6690657f6a4d2f4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.540279 INFO::Fitting model to feature number 430, X0836d28e24ff3d9df1d6dab03c659c54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.575042 INFO::Fitting model to feature number 431, X2b989c0ee8f9bb3ad5b1f572948e6ebc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.607383 INFO::Fitting model to feature number 432, X41813981fc56ece9178aba2271d500ed
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.64612 INFO::Fitting model to feature number 433, X29983e552e325387562ed544443fa784
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.678868 INFO::Fitting model to feature number 434, X453de2772d9c68cfaa47792539745455
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.711349 INFO::Fitting model to feature number 435, X7e947a204620ae22a12b10f2007e3bfd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.742738 INFO::Fitting model to feature number 436, X8cee7e94f2af5f8c3e4def902cb4a4e3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.785909 INFO::Fitting model to feature number 437, a45df882254eed21da6066c7ad4e7664
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.819309 INFO::Fitting model to feature number 438, fee2f4467e991748a4950fb60d0eccda
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.850072 INFO::Fitting model to feature number 439, X77f68b854b479c56f7cd35baa1bbb572
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.882041 INFO::Fitting model to feature number 440, fbcb14cd27216d317b57152441189a35
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.919341 INFO::Fitting model to feature number 441, X65441a076e1d645106982a3e2423e0b3
## 2026-03-04 18:36:23.95684 INFO::Fitting model to feature number 442, X08e77766d42e5c1a91eb3a03251301e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:23.995267 INFO::Fitting model to feature number 443, X1d7ac3928596c06f8fc853280fd38fb8
## 2026-03-04 18:36:24.026339 INFO::Fitting model to feature number 444, e59f0d7208f53430a8d7ac10438fd9cc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.059349 INFO::Fitting model to feature number 445, b0c3763d12e20e512b26e2db2e3d6746
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.091779 INFO::Fitting model to feature number 446, X830c92e2d730f1ac8f4e9435bf6a38a5
## 2026-03-04 18:36:24.125674 INFO::Fitting model to feature number 447, X185dbca6abbfd276dc552de656dcf054
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.159652 INFO::Fitting model to feature number 448, cf7fa17e32f98d8284e9100e9091ed87
## 2026-03-04 18:36:24.19198 INFO::Fitting model to feature number 449, X85ae192e589c308e221bea6cd2b908eb
## 2026-03-04 18:36:24.225425 INFO::Fitting model to feature number 450, cc9396a0baae4a61f1657f60aa953c58
## 2026-03-04 18:36:24.255712 INFO::Fitting model to feature number 451, X67eac2ae5c9d17d4559785aafb0c895e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.287476 INFO::Fitting model to feature number 452, X7f2b1c292a43ac11066e41099d33bf57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.319126 INFO::Fitting model to feature number 453, X3587846c90a1f5ac601a049f34c8ba6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.350708 INFO::Fitting model to feature number 454, b5256642df3fc7ccd4b8b4a0f9bb1fa2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.382993 INFO::Fitting model to feature number 455, f77feba135e606ed346677c4ba78fc71
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.414879 INFO::Fitting model to feature number 456, X8f1a91087f5fa731da48b8f2b92ad6c2
## 2026-03-04 18:36:24.446881 INFO::Fitting model to feature number 457, b73f90b9f991bfbdf47d82eaf2087ea8
## 2026-03-04 18:36:24.478022 INFO::Fitting model to feature number 458, X3e55881dab6323e1f88a1163c2a44acb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.514628 INFO::Fitting model to feature number 459, X290c59a265e798e44f2e3e519155f58c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.545333 INFO::Fitting model to feature number 460, aa75bb71f1660c6e77bcd7b8c6cae618
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.57827 INFO::Fitting model to feature number 461, b65e5f621c1b50c36382b4a7caf924d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.623913 INFO::Fitting model to feature number 462, d722d9be7b36fe09e126853ef8355c6c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.65488 INFO::Fitting model to feature number 463, f04f69734d1e61c54ad34db724229aa5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.688733 INFO::Fitting model to feature number 464, X466b76e8660a84049eea28a8c28b9585
## 2026-03-04 18:36:24.725731 INFO::Fitting model to feature number 465, e7f4b79576c6742e915e9c303be5078c
## 2026-03-04 18:36:24.757455 INFO::Fitting model to feature number 466, ad251f7dd0a17d8e8c36a3f4a5d28dee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.794752 INFO::Fitting model to feature number 467, c4bbccfc8721e556270f6371b42deb5c
## 2026-03-04 18:36:24.825556 INFO::Fitting model to feature number 468, X287d63e09f7a925d1a73d4355be6d435
## 2026-03-04 18:36:24.855494 INFO::Fitting model to feature number 469, X59ff6f3e5e1ea3243d69c0e3a574d33a
## 2026-03-04 18:36:24.885951 INFO::Fitting model to feature number 470, X0c6c3ad654b756b0d03f35141d790a1a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:24.923189 INFO::Fitting model to feature number 471, X0c51753914969dbd012afbd6bd03c9be
## 2026-03-04 18:36:24.960415 INFO::Fitting model to feature number 472, X7203fe878b61ef50e6220aa58dbe777e
## 2026-03-04 18:36:24.990757 INFO::Fitting model to feature number 473, X06f75349b62c026bd7ee4b1e83db846f
## 2026-03-04 18:36:25.0225 INFO::Fitting model to feature number 474, X85d08c22bb731d3f1be31c53307dd020
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.053819 INFO::Fitting model to feature number 475, c7c85584ed62586869cb651f2b2f7616
## 2026-03-04 18:36:25.08348 INFO::Fitting model to feature number 476, X09d64e76596ad78536416e251990587e
## 2026-03-04 18:36:25.115633 INFO::Fitting model to feature number 477, b59ac55dd4ff09c33b997121f6a818cf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.152697 INFO::Fitting model to feature number 478, b007b4058009181e2dbff0c93d0aa6f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.183172 INFO::Fitting model to feature number 479, X675e3bde986941b3b9117125b452649f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.220746 INFO::Fitting model to feature number 480, ba8218c95403215e7b6016df996ea0c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.255999 INFO::Fitting model to feature number 481, X07aa2d5d9a104b56b558e58e1cd206d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.289368 INFO::Fitting model to feature number 482, cefc53ad5dd34001d4c7115b3a40f868
## 2026-03-04 18:36:25.325216 INFO::Fitting model to feature number 483, bad88c92afdf5e15954b44dbf952aea0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.356368 INFO::Fitting model to feature number 484, X0af0efc520539e982c5d5d36ed1a2287
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.387897 INFO::Fitting model to feature number 485, e2eb429b518244dbbde8eb1d986782f6
## 2026-03-04 18:36:25.417997 INFO::Fitting model to feature number 486, X2d6e9a6ac74bc29426134bcb1308bba5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.458345 INFO::Fitting model to feature number 487, X5793519a2f18bf75476acdace2873d83
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.488912 INFO::Fitting model to feature number 488, e98c5cfe89c70d3a054bf7c0bb632f99
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.523871 INFO::Fitting model to feature number 489, X601cc059b0b76f341ba3599749a94e56
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.555091 INFO::Fitting model to feature number 490, X686e2a31809281a578955a44515b0986
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.58822 INFO::Fitting model to feature number 491, X152532a23717c40d3e3c4dbdfb5b7be9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.619334 INFO::Fitting model to feature number 492, X60442eb0bb422b1646950130e527ff43
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.650206 INFO::Fitting model to feature number 493, a3be4713e16df274dbd40b0e6637bb6d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.681308 INFO::Fitting model to feature number 494, X4be568d1f92c7f820a9aa816f1c837df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.715371 INFO::Fitting model to feature number 495, X6f18e8537683b2c00b5aa2d653299c54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.748306 INFO::Fitting model to feature number 496, be66db8135bee6134a930acfb7e05c8c
## 2026-03-04 18:36:25.780652 INFO::Fitting model to feature number 497, X3a95b5b429aaf1252730d3ce59aa1157
## 2026-03-04 18:36:25.818216 INFO::Fitting model to feature number 498, fc224cd4a1ffe36b74bdec848bbba352
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.857577 INFO::Fitting model to feature number 499, X88685671fbdb6afbee06252a97db86bb
## 2026-03-04 18:36:25.888215 INFO::Fitting model to feature number 500, fa0af2b39d42b4ba5ae59ed015cf909a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.919546 INFO::Fitting model to feature number 501, X2bd74498582c27631f2164b64afb353e
## 2026-03-04 18:36:25.953224 INFO::Fitting model to feature number 502, X08a62cc1a58e42059939f913aa50e77a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:25.985514 INFO::Fitting model to feature number 503, X82c0fba4e73ea914e1ef1236ab224aa1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.017481 INFO::Fitting model to feature number 504, d806628b6dbda1b2faf42d2620278fe1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.049974 INFO::Fitting model to feature number 505, X828d98340f164d377d2dda60bb06f4a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.08071 INFO::Fitting model to feature number 506, X8ca98afd451a98fa6cf01f63bed22a96
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.111581 INFO::Fitting model to feature number 507, X5a53bba7d54f6d0157de0bb670839036
## 2026-03-04 18:36:26.143053 INFO::Fitting model to feature number 508, cdae58c0f1656c07f87431029fb464b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.175053 INFO::Fitting model to feature number 509, ce1c05cf3572b9f9d78199a9ff5d19c2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.206907 INFO::Fitting model to feature number 510, ce70d8fb6e8b521ca9594a0134697e3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.23963 INFO::Fitting model to feature number 511, X8f70d1241904ec8e899ec6fab95d100f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.284153 INFO::Fitting model to feature number 512, X5c4542612816ac5896eb8c902e134056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.322566 INFO::Fitting model to feature number 513, X051517485aea932a161723ecf6c596ff
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.35293 INFO::Fitting model to feature number 514, a8a38efe4ef56bb85642c5526e0b0201
## 2026-03-04 18:36:26.387393 INFO::Fitting model to feature number 515, X33d599964dd4ead082d0504ab37f2e4e
## 2026-03-04 18:36:26.424535 INFO::Fitting model to feature number 516, X0845f16afd324e69265811e5f899b19e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.462524 INFO::Fitting model to feature number 517, e4ff331359aae89dbe43026ee8d9dff3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.493157 INFO::Fitting model to feature number 518, X5a47b6a3ee0587defec53d95e2c842a3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.523984 INFO::Fitting model to feature number 519, X0f828935964def127d6255cdef582020
## 2026-03-04 18:36:26.557365 INFO::Fitting model to feature number 520, X1095b8196dbd38d664d439fc8a3e7634
## 2026-03-04 18:36:26.589841 INFO::Fitting model to feature number 521, X9749d4ed81f1f57122e213e1f8de77ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.621433 INFO::Fitting model to feature number 522, e738ffccde00c7cf4f544a818c42172a
## 2026-03-04 18:36:26.653816 INFO::Fitting model to feature number 523, X0f0aab6e8fef47c40ea53cebf01f27f7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.68646 INFO::Fitting model to feature number 524, X0b759bf20b1f3dc4e3c29eaca43676aa
## 2026-03-04 18:36:26.718761 INFO::Fitting model to feature number 525, X4103c51e17a726ef97db183e51c5aabb
## 2026-03-04 18:36:26.757734 INFO::Fitting model to feature number 526, X29744e1859bbc3cd114636ac2f0e53dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.789003 INFO::Fitting model to feature number 527, X1d045853283b902a8ce4ccf636ee1c5d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.821886 INFO::Fitting model to feature number 528, X7256084d97c5df6b2aed179843a97804
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.853562 INFO::Fitting model to feature number 529, dfd8b8de8c4401fed140e29521b9b033
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.884155 INFO::Fitting model to feature number 530, X82c789134ae6f1ce1eb7a21113d59602
## 2026-03-04 18:36:26.924283 INFO::Fitting model to feature number 531, X3bbce5f40b435eb449659929e78a17f4
## 2026-03-04 18:36:26.955759 INFO::Fitting model to feature number 532, X46326463c37cc1a8c3236cf2ce812bf9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:26.992175 INFO::Fitting model to feature number 533, X55ce23b59a993eb2b7741fde91da583a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.022919 INFO::Fitting model to feature number 534, X463f1b971157c801fbf3528011b4546e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.054163 INFO::Fitting model to feature number 535, X382dd0303518af0483be8c4ea08c74c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.084872 INFO::Fitting model to feature number 536, fb1708d9521d280ffd129fc84389ab24
## 2026-03-04 18:36:27.127674 INFO::Fitting model to feature number 537, X959bc5f6827ce5a4150f8c5c2608141e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.159615 INFO::Fitting model to feature number 538, X3af0184defe88fe12dd776417e66637c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.190983 INFO::Fitting model to feature number 539, X124d64c0e562f0d45b2a443f8ad1fef9
## 2026-03-04 18:36:27.222444 INFO::Fitting model to feature number 540, b9d2f025be5b0d19c8d6b8249ed48287
## 2026-03-04 18:36:27.252372 INFO::Fitting model to feature number 541, X47289f4be627b4ebd6a107b9d767b38b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.28855 INFO::Fitting model to feature number 542, X1ec3481ae3ea83074d29ca1e64c26cba
## 2026-03-04 18:36:27.323347 INFO::Fitting model to feature number 543, X5c3085994e1bbbf662babdfaf0ccf8bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.354252 INFO::Fitting model to feature number 544, X084347470f75d9e030c777fab7bf5930
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.386432 INFO::Fitting model to feature number 545, X34a5cc3ff0c374acc93606948a77b109
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.417363 INFO::Fitting model to feature number 546, c57234d6515915a9ccd83265178c38a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.452725 INFO::Fitting model to feature number 547, X1630c97bc99db59912c92545f184a2c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.484281 INFO::Fitting model to feature number 548, aa7bae079433f6ed40a3d3016cb2b4e1
## 2026-03-04 18:36:27.521521 INFO::Fitting model to feature number 549, X0b9987df25b2cbf3e32983f0bb343d22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.552545 INFO::Fitting model to feature number 550, X4fbf84a3fdfe517ab4a03517aa50188d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.583452 INFO::Fitting model to feature number 551, cc628ddb3cb5ec7eca660a77860a0dfa
## 2026-03-04 18:36:27.628595 INFO::Fitting model to feature number 552, a5d3cef1b84bffa024ae0259ebbe7819
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.659609 INFO::Fitting model to feature number 553, X001d2e9ba87c0b0118ae91de96ff290f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.69084 INFO::Fitting model to feature number 554, X1f36cbd271a0a85fe9d8415377e611c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.728192 INFO::Fitting model to feature number 555, X421c87630f8796da746888ec4aa55f4f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.759055 INFO::Fitting model to feature number 556, X6eff287d92b730a0169e0b2142f31b8d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.791408 INFO::Fitting model to feature number 557, X039ec62846d8972fa2fa02f2bd2056f7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.826302 INFO::Fitting model to feature number 558, X33cec55e74e9b122f1e03586050d529e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.857119 INFO::Fitting model to feature number 559, a8982ccbc26415e739f9c1dea2471f4f
## 2026-03-04 18:36:27.89228 INFO::Fitting model to feature number 560, ae4d72c28426dd93660524df503c9d97
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.925584 INFO::Fitting model to feature number 561, X173f0a0e35c5a252dcaddc1c05e40cc6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:27.969304 INFO::Fitting model to feature number 562, c88ac9ec2ae50bb75443a61e09401ebd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.006784 INFO::Fitting model to feature number 563, X9c10353ae637b8f0c7fb09b60c0d18d3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.041289 INFO::Fitting model to feature number 564, ba4f8c14f1a6f625c90f05848475b43a
## 2026-03-04 18:36:28.075466 INFO::Fitting model to feature number 565, X42a703e9af69e2172af67521bbc356cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.107519 INFO::Fitting model to feature number 566, da4cb8754d5d4a0317d887304c60b9f8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.138982 INFO::Fitting model to feature number 567, a2f33842d6158bd3ff1162b1b575109c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.172866 INFO::Fitting model to feature number 568, X18e76539a844edacd926d64850b8aa36
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.204982 INFO::Fitting model to feature number 569, X218c575bc466bd7f3de4855d1cb785c9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.23582 INFO::Fitting model to feature number 570, X527c201e40ff4e8b5562ed88b67ead25
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.267377 INFO::Fitting model to feature number 571, b1f3d0e531c2a6e1fc1ca63b5c5506e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.300809 INFO::Fitting model to feature number 572, X252d977d5b531a6a265d557f412e2d44
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.339931 INFO::Fitting model to feature number 573, X3015c2163a624808acd0a7c62a3bf0c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.371706 INFO::Fitting model to feature number 574, X3e1bcd69c044a695703a1985de3bd920
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.405639 INFO::Fitting model to feature number 575, X5d42ffdab90cffce3e9c2cfb952cc2fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.44121 INFO::Fitting model to feature number 576, dbb13a4636bef149f210c12027d77f9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.472977 INFO::Fitting model to feature number 577, X5ad1570deadbd75f87675e54addf230d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.503101 INFO::Fitting model to feature number 578, X13562e00ec4db4d18cd5dd1bd73d2ef3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.538767 INFO::Fitting model to feature number 579, X69419db7063b243b238cc5bae949369a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.569789 INFO::Fitting model to feature number 580, X91f7efb223b32216b23bb9a02e9ec570
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.602055 INFO::Fitting model to feature number 581, dfc46dd0818f492b120349ace278d6d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.634788 INFO::Fitting model to feature number 582, fcecd1c29c2e37c8c3de46a1926cf892
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.673189 INFO::Fitting model to feature number 583, X4a396a487ee50aa361519dc2a2369a5f
## 2026-03-04 18:36:28.711679 INFO::Fitting model to feature number 584, eee8f061f4328896cbad20ccb6bf41f9
## 2026-03-04 18:36:28.746282 INFO::Fitting model to feature number 585, X5cedacc9704fb36ee1087cdc9bdf80ed
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.77797 INFO::Fitting model to feature number 586, b83acb9199413447ed91faec843cd908
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.826717 INFO::Fitting model to feature number 587, X6dc0a98dc381716943fa9393aa9cb416
## 2026-03-04 18:36:28.859457 INFO::Fitting model to feature number 588, aed586669a5eb4831e60da011eac8ca9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.891226 INFO::Fitting model to feature number 589, d2b971f10e7c69e1d4c89c63cfbe9b94
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.930512 INFO::Fitting model to feature number 590, X00c145ccfb95cec516f4fa326d896af4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:28.970936 INFO::Fitting model to feature number 591, X06b4246e0976c80b3b656da49f9599c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.008255 INFO::Fitting model to feature number 592, X320902985335212dfccb275337dccacc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.039161 INFO::Fitting model to feature number 593, X3765a7f3a9962418ea2f3e764d644898
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.071409 INFO::Fitting model to feature number 594, b1d5ee1846c492211f4fed9aeb829a84
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.102583 INFO::Fitting model to feature number 595, daa12bd954485797a7fd1a0f2b9024eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.134493 INFO::Fitting model to feature number 596, dfb91dda0f7b903c88c7082678e4ce4d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.16522 INFO::Fitting model to feature number 597, e9d0e82558715ca0ca438f41fbbe4ef6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.198453 INFO::Fitting model to feature number 598, d35bd24b097c811e40172a1fc6d15a5e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.233889 INFO::Fitting model to feature number 599, e737489b71d52736ee1a0e2a50549d4c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.264239 INFO::Fitting model to feature number 600, cd5308d7409829a5057c4029560fac2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.301813 INFO::Fitting model to feature number 601, fe72c667e0eb30b792b6b59fb292c54c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.333513 INFO::Fitting model to feature number 602, X1d4285e838c09afb01b627dae637beb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.364755 INFO::Fitting model to feature number 603, X612fbc07b28756a2e743392fde08e155
## 2026-03-04 18:36:29.39535 INFO::Fitting model to feature number 604, X76a776e3b51070fc91491000f4030ae1
## 2026-03-04 18:36:29.425347 INFO::Fitting model to feature number 605, X2e246d891a50a35543494663b056b11b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.456374 INFO::Fitting model to feature number 606, X07f2b771a3355046b9ac7e90670a48b6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.487167 INFO::Fitting model to feature number 607, X6cecdfc4f05450d653299bdbdc9f26db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.518945 INFO::Fitting model to feature number 608, c7f6f5b0e601f6374c521ec454e23a90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.549715 INFO::Fitting model to feature number 609, ddd96397751cb6c487e71dca3fb8eedb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.580616 INFO::Fitting model to feature number 610, ee391c5cbd032373f7fa130bbdfa78b2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.622198 INFO::Fitting model to feature number 611, X331d529d33e29b0d43456b0ea53b88fe
## 2026-03-04 18:36:29.66115 INFO::Fitting model to feature number 612, b483981d9b053abefc0935cb20c58b66
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.695703 INFO::Fitting model to feature number 613, e82710e82def354d4332ee5b055234c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.728447 INFO::Fitting model to feature number 614, bacc8fba734fc18196d62879e73b30c2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.762484 INFO::Fitting model to feature number 615, X00e9012a608d346ef966588502aff2df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.794901 INFO::Fitting model to feature number 616, X171c685c03b82e7f5fbe12bd0ac2e720
## 2026-03-04 18:36:29.826833 INFO::Fitting model to feature number 617, X387aece6e4ad6511151d09055a8d2d5f
## 2026-03-04 18:36:29.860155 INFO::Fitting model to feature number 618, X510b61e30ef7cba2b46623f71a09d2bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.892279 INFO::Fitting model to feature number 619, X925c84e4221675297ab469c3758bd13e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.923653 INFO::Fitting model to feature number 620, a4be22de4b5390e64bfe2b6011523f3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.962869 INFO::Fitting model to feature number 621, X1764219dbd7ca6405b38a4f1122b5270
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:29.998288 INFO::Fitting model to feature number 622, cec7bcba266e98babcbadffedf294db7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.03125 INFO::Fitting model to feature number 623, e7929ae3c3f9c7c26009da7f71875850
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.069937 INFO::Fitting model to feature number 624, f5cd996896ec42be1a2d4513c5d663af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.101057 INFO::Fitting model to feature number 625, X1e0bd8313bb83eed0c2ca563b106ea19
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.134173 INFO::Fitting model to feature number 626, X12a7a583a9d85fc8ca1bad3a26192c2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.173012 INFO::Fitting model to feature number 627, X1ecbf29ff5cb9d6bf1eb427bcb0290e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.210122 INFO::Fitting model to feature number 628, d5b669bb5d3bff1c87f4c50299e68a13
## 2026-03-04 18:36:30.247909 INFO::Fitting model to feature number 629, de809f907097c58da02f0c46a3f761d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.284896 INFO::Fitting model to feature number 630, X1208c5d5a63919dfdce05b8942fca118
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.317129 INFO::Fitting model to feature number 631, X263e8a3f16a6c22a7991e8ee33351391
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.351615 INFO::Fitting model to feature number 632, X68f8d7429d5fe7dc889f92705924aa70
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.383141 INFO::Fitting model to feature number 633, fda375197ec72a134ce3c008754dffae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.414456 INFO::Fitting model to feature number 634, X3f25c26eb5e2ad00677e1c3341010d70
## 2026-03-04 18:36:30.445479 INFO::Fitting model to feature number 635, f96dbf5f5deb0558ba732d9e264601a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.487199 INFO::Fitting model to feature number 636, daaf8b158006c74d73fee01ca65fb5e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.525748 INFO::Fitting model to feature number 637, X04314b4b349b7be1b1adca749349e9d5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.565784 INFO::Fitting model to feature number 638, X24cb510a92ff13206c2c50e208ea26da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.596524 INFO::Fitting model to feature number 639, X0e372f45908a23ea01a2cfc4db466bd5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.627062 INFO::Fitting model to feature number 640, f5375d7a11214ac9f11f523325b2348b
## 2026-03-04 18:36:30.659223 INFO::Fitting model to feature number 641, X82d1fcf8762f0303c6a4c50d61e5ff8d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.691667 INFO::Fitting model to feature number 642, e903cac37add63e8c94650b15bfefd86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.723807 INFO::Fitting model to feature number 643, X6fe94b71e60dfd02b295e38a5f3e8da2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.756674 INFO::Fitting model to feature number 644, X17d13995f3dc073bc150a3ce1f936667
## 2026-03-04 18:36:30.788483 INFO::Fitting model to feature number 645, X9c3ba0cc2b8bcadd5776b16bf7947e5a
## 2026-03-04 18:36:30.819452 INFO::Fitting model to feature number 646, X73a596d4c84e2cc607f447b8d86e410d
## 2026-03-04 18:36:30.856843 INFO::Fitting model to feature number 647, X6936efdd0d4dc2c78a52637bbca6698e
## 2026-03-04 18:36:30.894979 INFO::Fitting model to feature number 648, X1accfe46cbb5fa4b43f42f12a4b09a46
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.928163 INFO::Fitting model to feature number 649, X372ae668689d07c1776f9e39c19b7469
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:30.966231 INFO::Fitting model to feature number 650, X3586e5ed26b0fe9d2b88c41818b48636
## 2026-03-04 18:36:30.996056 INFO::Fitting model to feature number 651, X4ad6f2f290bec0999affbdae99e9d315
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.026745 INFO::Fitting model to feature number 652, X43722f6ec1f24a34bdef18e320dd41a4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.058319 INFO::Fitting model to feature number 653, X6205f70998dd9ef2c6c8769ed201ca0e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.089019 INFO::Fitting model to feature number 654, de7b4217cc9e4c4ccabf081ab8157af8
## 2026-03-04 18:36:31.119103 INFO::Fitting model to feature number 655, X8d25b30e61a1542d118f0a1d96420b8b
## 2026-03-04 18:36:31.157257 INFO::Fitting model to feature number 656, X16508d3f8de41b678f64b67a13087618
## 2026-03-04 18:36:31.187708 INFO::Fitting model to feature number 657, X7bf58f84d1c8c8d1d013725ba2a6d42c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.218508 INFO::Fitting model to feature number 658, c567686addbcf1ec7c0326a8cb050ee3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.256834 INFO::Fitting model to feature number 659, X8ed345e1519449d98c9ae839e788ed75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.294667 INFO::Fitting model to feature number 660, feea58b2edca226cce18140b7a96f212
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.336727 INFO::Fitting model to feature number 661, ffaf1f6f3aa2f47240308a9757599ff3
## 2026-03-04 18:36:31.371833 INFO::Fitting model to feature number 662, ee7a112d25ef18a5ac190015e789ab7c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.406444 INFO::Fitting model to feature number 663, X8b79ecbf1789202e32c2134d323f3fe5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.436933 INFO::Fitting model to feature number 664, c56a8b428e298d443b0046bd7da3b033
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.47152 INFO::Fitting model to feature number 665, e08bffe7c47dc649f45e6605eaf77c00
## 2026-03-04 18:36:31.508189 INFO::Fitting model to feature number 666, X390bd5af0dfc36880157790449a6934a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.545295 INFO::Fitting model to feature number 667, abe04980ffc876dd5ddefee10c71e8da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.575869 INFO::Fitting model to feature number 668, X03b7c9224944f7a08e43db638d43e1e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.607635 INFO::Fitting model to feature number 669, X6713a01fce8ed800fdea27cfef912766
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.642962 INFO::Fitting model to feature number 670, X8a7b1595d4369644e07e85b86c70e3db
## 2026-03-04 18:36:31.675522 INFO::Fitting model to feature number 671, e4a9ae0bba905d4ac7db962abab2e295
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.710696 INFO::Fitting model to feature number 672, X4855f68229e7d81a2e9da0111a5ee938
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.742694 INFO::Fitting model to feature number 673, X9700a2382c70667f3daa3cb2a55d0a7a
## 2026-03-04 18:36:31.773702 INFO::Fitting model to feature number 674, X24226920f12955b61a66789e6ad6f215
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.804982 INFO::Fitting model to feature number 675, X4c0ab4869eec1c6667ca4fc83f477e39
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.842101 INFO::Fitting model to feature number 676, X02f8bdb4409411a2601dc87c551942a1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.874831 INFO::Fitting model to feature number 677, edce391df3040a66ce3f198029a9873c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.912696 INFO::Fitting model to feature number 678, X3b3e296ef6b3c79cdc410d0d662a6237
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.943802 INFO::Fitting model to feature number 679, eaea5405910f416f1be4d4d0b71ac90b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:31.980613 INFO::Fitting model to feature number 680, X0b78143d11a6ca0fa5137a6b875c3a51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.012076 INFO::Fitting model to feature number 681, a832ff8e0373f88c777d1e717dc8cf21
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.049177 INFO::Fitting model to feature number 682, b9f2ffae17e0a92098de47a3b3ea81c8
## 2026-03-04 18:36:32.081173 INFO::Fitting model to feature number 683, bad508016651bf74cf3cb269911ec12b
## 2026-03-04 18:36:32.117855 INFO::Fitting model to feature number 684, X213c0a6d6aaa77c7237d1c27dd8f3a78
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.148214 INFO::Fitting model to feature number 685, X43c87171c86c75f5fd3ff3da1f251b0f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.189012 INFO::Fitting model to feature number 686, X49a98f8f4510c8c171fc2ba1461073e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.221263 INFO::Fitting model to feature number 687, X72388b9f07ab344f60eb50a9fb4ff6b0
## 2026-03-04 18:36:32.252712 INFO::Fitting model to feature number 688, X50577f758ed98fa0971590d1f3818602
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.28337 INFO::Fitting model to feature number 689, c74a32f53cfb587c14d6d0326b225a5e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.315373 INFO::Fitting model to feature number 690, e7eeb1fe634e3182ff6d596bcb0674c1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.346293 INFO::Fitting model to feature number 691, X4c1637f4177c0b2d8f35befab7c70508
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.379051 INFO::Fitting model to feature number 692, X13420d21573c308fc7d77a3f77b01f6f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.410497 INFO::Fitting model to feature number 693, X9a0c6caa7312d271153e400a3e15903a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.448295 INFO::Fitting model to feature number 694, X9e8b037cd3a3c0d0c1abb98279a78e5b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.480256 INFO::Fitting model to feature number 695, b311631596b18e1f7ae54d3b969844c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.511448 INFO::Fitting model to feature number 696, ddc7f083639aa23de81829ed19d1faa2
## 2026-03-04 18:36:32.541878 INFO::Fitting model to feature number 697, X92fdc7ec8d89c654aeb51aa18a540bf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.572691 INFO::Fitting model to feature number 698, X537d57d1ee4850a807f6c6bb9f477897
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.604352 INFO::Fitting model to feature number 699, edf3da53246c455f9b500ac2db192e1e
## 2026-03-04 18:36:32.636539 INFO::Fitting model to feature number 700, bae5efc0e9e141ca3aa86b1f58ee238e
## 2026-03-04 18:36:32.667344 INFO::Fitting model to feature number 701, d8ba154ce25da930b156278285fa1afd
## 2026-03-04 18:36:32.697904 INFO::Fitting model to feature number 702, X0c101f614e28abcff7de67c17ca51a55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.729483 INFO::Fitting model to feature number 703, X7ea0a5c2e540405f53d0c30fa63c82a9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.75989 INFO::Fitting model to feature number 704, ddce072baa0bb0342e2fd9a8cd7cf489
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.789951 INFO::Fitting model to feature number 705, ea88fdd7f58f1c24015d406597798885
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.819803 INFO::Fitting model to feature number 706, d294afa38c6e4bc2e50c1bc296583d22
## 2026-03-04 18:36:32.85371 INFO::Fitting model to feature number 707, X1d7bbdef9955168131018ddbaf14993f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.883749 INFO::Fitting model to feature number 708, X1de141b5e2a83b1fbef3d71352ba38fc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.913759 INFO::Fitting model to feature number 709, X3797616e5f13944c402ba8e5c97de710
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.94398 INFO::Fitting model to feature number 710, X82bc96ed61841410472a980eadcc456a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:32.985229 INFO::Fitting model to feature number 711, ac9d52c7b33c42e07fc95741ebdccfb4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.016501 INFO::Fitting model to feature number 712, ba1ae9ef946c2c8619f46f0a85e39876
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.048125 INFO::Fitting model to feature number 713, c61c0fc771ac44540c315cf1e48285b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.085693 INFO::Fitting model to feature number 714, cb97cafd5cd9403ff6238b886b39a5cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.121925 INFO::Fitting model to feature number 715, cf5e11d1b9fdca94498f62c9f6a9163e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.153909 INFO::Fitting model to feature number 716, X36135ac55ed4d0ff6e31d6ff7aa4499e
## 2026-03-04 18:36:33.186996 INFO::Fitting model to feature number 717, X8371b60e37ffedd9c1a9a37fef3e3286
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.225767 INFO::Fitting model to feature number 718, a1efd505557860248435d7bea1bf3913
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.264336 INFO::Fitting model to feature number 719, d3c7a110bef8e513fdc270909827b2ea
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.300605 INFO::Fitting model to feature number 720, X3a5942b1a3e4aab99de89f5bac550156
## 2026-03-04 18:36:33.331799 INFO::Fitting model to feature number 721, X43e163b10a20212c4a092bdafae4f956
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.363512 INFO::Fitting model to feature number 722, X5b3835118edb11aeb9d4bab295d72cb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.394219 INFO::Fitting model to feature number 723, X6617d00c29734b4e9411fa53e40f16c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.424415 INFO::Fitting model to feature number 724, dbc6c4998c0721a95236e008182051c4
## 2026-03-04 18:36:33.456703 INFO::Fitting model to feature number 725, X6c92c5cc06059722f26df4dbf9eda96e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.487063 INFO::Fitting model to feature number 726, X82c546ae87e738cdf074853e15bed022
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.517592 INFO::Fitting model to feature number 727, b20e3d817c8a6fc306663500024b2479
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.54782 INFO::Fitting model to feature number 728, c5530dfee5b10d94b0fe01f9dbcb002f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.584552 INFO::Fitting model to feature number 729, X48a8fe019dc6cc1cfb047a123527b5c0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.614488 INFO::Fitting model to feature number 730, X77d5f38002bb3780f14d3d9880aa9008
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.645267 INFO::Fitting model to feature number 731, X8a7640fd3fdee3d1d7dc786a9ddd6673
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.675945 INFO::Fitting model to feature number 732, a66208cbb348d255e7357fae6c83e793
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.708094 INFO::Fitting model to feature number 733, b153c75d7dae701723a3e55a0e7610a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.745233 INFO::Fitting model to feature number 734, X1230aedcda2b9095c4b718d403c11b65
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.777978 INFO::Fitting model to feature number 735, X1d785c4e2c7f07957b49d819a3c473d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.819393 INFO::Fitting model to feature number 736, X2630cd266dd2f91e8e73919d7ec61b3c
## 2026-03-04 18:36:33.854289 INFO::Fitting model to feature number 737, X3adbb59e7c520b6a4948e315e22aa48b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.893252 INFO::Fitting model to feature number 738, X4547cdaa4d234e3ca4e24d18dbb1a7cb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.927402 INFO::Fitting model to feature number 739, b41174e1a26fcc77c1d1fecf6a39e3eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.960273 INFO::Fitting model to feature number 740, bc31699d255fca6a180a44290780aa5d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:33.992017 INFO::Fitting model to feature number 741, X3b918737053456dd3b8f2c656205dc3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.031223 INFO::Fitting model to feature number 742, X3e9a95e009b5e5728ca7b5861a9f6886
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.064341 INFO::Fitting model to feature number 743, X723bb2a4f2424433622e58cb4c184f68
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.098651 INFO::Fitting model to feature number 744, X832b2920c7d4c0c22e2d9b24f716e64a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.137438 INFO::Fitting model to feature number 745, X936c188690c04d1eec4e44cd1f966362
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.16889 INFO::Fitting model to feature number 746, bc7431fe3e8d871068b783cb153d98ea
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.205406 INFO::Fitting model to feature number 747, de1af36e20ecfa838185c13b9e66c677
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.240984 INFO::Fitting model to feature number 748, df6b94f16e20f444ca41aa31497124eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.277393 INFO::Fitting model to feature number 749, X1a66cd529f8773b8626aac7adc240a1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.316667 INFO::Fitting model to feature number 750, X2fd2e9098fb21a08c0855b9f03580f79
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.353433 INFO::Fitting model to feature number 751, X91d9ae565c05f931c2ceb57abeffbce8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.389096 INFO::Fitting model to feature number 752, d2dcbb92bf9cc3ab8d247629e9e07526
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.426888 INFO::Fitting model to feature number 753, X075dc0a8cf360e04255c741bde0699d0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.458177 INFO::Fitting model to feature number 754, X25ee65871f64cda654d5f4db4e7aea31
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.490572 INFO::Fitting model to feature number 755, X8bb457c94ddbae8fd75111d1598017aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.522223 INFO::Fitting model to feature number 756, d89d081df8bbf354930bda36c31f8d00
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.553726 INFO::Fitting model to feature number 757, e91ff13a5dab0760ae86ed2b229a578d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.587738 INFO::Fitting model to feature number 758, X0a62aa5adaaa05fc8876ed148f03fee9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.618317 INFO::Fitting model to feature number 759, X5f54914fe6d2951bfab121606df7ab31
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.66489 INFO::Fitting model to feature number 760, X8fb569330957254554e77948ce449cac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.7023 INFO::Fitting model to feature number 761, a85443cff26b1389b3ecc414e37c7da6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.739351 INFO::Fitting model to feature number 762, b4273d59573cc45dc778e776f176b9aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.769712 INFO::Fitting model to feature number 763, de57ef35f117467f9b193d96e53ae163
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.805672 INFO::Fitting model to feature number 764, X2886d1db4c3b6787f79e50bef8f69c55
## 2026-03-04 18:36:34.841507 INFO::Fitting model to feature number 765, X743535eee04707e3424f814948e07b47
## 2026-03-04 18:36:34.878377 INFO::Fitting model to feature number 766, X70a5b7339c7527388c360d247e8df171
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.909017 INFO::Fitting model to feature number 767, X4879cd5a48e851e25a3ca727ce334b32
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.945026 INFO::Fitting model to feature number 768, X1be31d76bb8cf9406332b7282354d758
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:34.981253 INFO::Fitting model to feature number 769, X163e8909e8fcf7a3b229cdbb8f26275f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.018086 INFO::Fitting model to feature number 770, cc81bf9bf985caa468eaa6b7fdfd1bf8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.055065 INFO::Fitting model to feature number 771, X04ceef6f70d00472bdd1af80c2f30983
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.092691 INFO::Fitting model to feature number 772, ae1394fa5a7795914411b0d3ea6546f1
## 2026-03-04 18:36:35.125037 INFO::Fitting model to feature number 773, da9409cba5c8301a2b42580176f14265
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.157944 INFO::Fitting model to feature number 774, X22d3ed64f87a280651657ccc0b581145
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.191056 INFO::Fitting model to feature number 775, X4b17e66aa631354e9f3f33314d089589
## 2026-03-04 18:36:35.222233 INFO::Fitting model to feature number 776, X83d8438154115532d2f599781802ff7a
## 2026-03-04 18:36:35.253068 INFO::Fitting model to feature number 777, X36ebd620b98f5ae6205c417f489e5207
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.285123 INFO::Fitting model to feature number 778, X5b550a3ccef0538fc08e87246a6a35b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.316893 INFO::Fitting model to feature number 779, X168b9624445c90028b89ca20afea6e4f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.35436 INFO::Fitting model to feature number 780, f436778965fbeae0dea00fc287ba37da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.39279 INFO::Fitting model to feature number 781, f54c94989621c5be35e72ff3e29893d2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.425944 INFO::Fitting model to feature number 782, X323474086152e9c4666a43f5eb515d05
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.458283 INFO::Fitting model to feature number 783, X96d483b6964c2f2ee2b237efb9f01baa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.496545 INFO::Fitting model to feature number 784, X3e12212b21a98a5799ebdf835788a4a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.543982 INFO::Fitting model to feature number 785, bef4a1751e1ee0ba5ba0366c85965eca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.575659 INFO::Fitting model to feature number 786, X6c69bfaaa63d25c2bf7e055ee84e3a90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.606194 INFO::Fitting model to feature number 787, X39f844378114cbe9256d291fff454c1f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.636386 INFO::Fitting model to feature number 788, X6455b26d66a9a71fe4411b6ba3b3ec4b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.666745 INFO::Fitting model to feature number 789, X5bdd8b2852cfb534356a12a5e39d96bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.697313 INFO::Fitting model to feature number 790, X41d4dbb4027e632037882ead150e113f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.727532 INFO::Fitting model to feature number 791, X25db1ff4778c4f47f626c8929563e2c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.758722 INFO::Fitting model to feature number 792, X79c05c36d834303543f0231cb21a2eac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.789921 INFO::Fitting model to feature number 793, fe54fff680b4cf948d8d943c2cc4f091
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.821372 INFO::Fitting model to feature number 794, X37d384450b5ee65ed1cc2bab56af5e1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.853773 INFO::Fitting model to feature number 795, X9b399f779ade5e5b6500eec0de6d13b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.891885 INFO::Fitting model to feature number 796, X6c04b85ba8beb48dfb0adfb7db7fbe7b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.924238 INFO::Fitting model to feature number 797, X8e0dc3caf705d8984bea5e9c501296d4
## 2026-03-04 18:36:35.962388 INFO::Fitting model to feature number 798, X18d2efc0bdd68758bba83ab330bd15bd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:35.999467 INFO::Fitting model to feature number 799, X0f3a73a457224cf1b94c2830fb8c6e86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.030021 INFO::Fitting model to feature number 800, cfcdf20e1255691d28fe0c56ce2016dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.060779 INFO::Fitting model to feature number 801, X21057ab2a1913644b88384f538082f65
## 2026-03-04 18:36:36.092318 INFO::Fitting model to feature number 802, X4889a0edb0de0f30c5f53c9aa5417068
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.123914 INFO::Fitting model to feature number 803, aa09d6b56ace46671bf60d93a5aab529
## 2026-03-04 18:36:36.154334 INFO::Fitting model to feature number 804, X0af82abc9d7b12a57264ed2f7d7c635c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.189568 INFO::Fitting model to feature number 805, X89c6d2a235c32ea0290509ae65221057
## 2026-03-04 18:36:36.22036 INFO::Fitting model to feature number 806, X34f8f51f3aa37533d2fc7d4cba1aefb9
## 2026-03-04 18:36:36.25277 INFO::Fitting model to feature number 807, X738286c5d38996bf4df0acd3eabbcc70
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.291352 INFO::Fitting model to feature number 808, ecc5e663311d028a84832f7d6db3d64b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.323355 INFO::Fitting model to feature number 809, c68dad1066d5ef8221fcc1093790b7ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.372828 INFO::Fitting model to feature number 810, X1e07bf1cff0e9e56f62458ab13710672
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.403857 INFO::Fitting model to feature number 811, X5c08d93436793f5fc51b871db4d85870
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.441746 INFO::Fitting model to feature number 812, X8070db26c5b69b8f7a174f40e71b085f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.47558 INFO::Fitting model to feature number 813, eef930d6e152847be082212f21d6fc69
## 2026-03-04 18:36:36.507536 INFO::Fitting model to feature number 814, X8dfbe52cb8feca704bc7d52c371d052d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.537865 INFO::Fitting model to feature number 815, X2958d05407011ced17c694c365dee47c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.568243 INFO::Fitting model to feature number 816, X5ad50b5f8b90a5fd0a409661a772acb2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.599795 INFO::Fitting model to feature number 817, X6de335ece7f713b9a8c71b18aaea1d4c
## 2026-03-04 18:36:36.630925 INFO::Fitting model to feature number 818, a0347f3b93155f49f127482da58b5d25
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.664506 INFO::Fitting model to feature number 819, X1d1fee5b3f5f2f03a6d43b68c228c61e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.696874 INFO::Fitting model to feature number 820, X3e232c33e572aeba36d4ac20b90b2822
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.733085 INFO::Fitting model to feature number 821, X59687c67da77385562ad90818c85468a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.77244 INFO::Fitting model to feature number 822, d3f8153113111d55a5a906a131a161bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.808304 INFO::Fitting model to feature number 823, X45819013a025366b17b161a0eacfe02c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.840144 INFO::Fitting model to feature number 824, X4c52e18b01f156b0b4030db975dd2c90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.8731 INFO::Fitting model to feature number 825, b432978e0c1aea4c678fef2f27d03ff1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.907186 INFO::Fitting model to feature number 826, eb0c2a0ec0c24342ab7639a869020215
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.938809 INFO::Fitting model to feature number 827, X1633280dca1ce6b99f323fb357a788f4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:36.97606 INFO::Fitting model to feature number 828, X19900408557ac9f0e79bf455cd982c2b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.006279 INFO::Fitting model to feature number 829, X0b2e8985a71922126f4f44935b75e37f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.043236 INFO::Fitting model to feature number 830, X5598141e7b6d1966331ae03cbd17c8ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.074048 INFO::Fitting model to feature number 831, X82e9bb8861cef7612bb7ae065798491e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.104153 INFO::Fitting model to feature number 832, X9c9fdbf346a77e3a795b6d9a7480b9bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.13815 INFO::Fitting model to feature number 833, e0c60e831deae33ffe603459de3fceff
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.175214 INFO::Fitting model to feature number 834, e92ba24ee64bcf3f4cb5f9a3682e4413
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.216736 INFO::Fitting model to feature number 835, X4d9a5071591ffc708722c047fe488fdf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.248135 INFO::Fitting model to feature number 836, X63c0016bfe6d090f37f597cd4df623a9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.278618 INFO::Fitting model to feature number 837, X977936e75d88dab18ff37a8e776040fe
## 2026-03-04 18:36:37.315317 INFO::Fitting model to feature number 838, b0c41adf5d284fbabdce3574ce198155
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.345827 INFO::Fitting model to feature number 839, be27327880cec05e98bf577491a9f35c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.376861 INFO::Fitting model to feature number 840, ee5a2aed5e233e505a43fd5e5d77f7b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.409101 INFO::Fitting model to feature number 841, X144a95b2a5248802ffb64716fbcee74d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.440089 INFO::Fitting model to feature number 842, X23d44c8c40fb72845b2f99c3f746be86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.473982 INFO::Fitting model to feature number 843, X245af114f65c6751ab4082702eb5717b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.504696 INFO::Fitting model to feature number 844, X52b05bd5be5f8a56cdac5e80ac038690
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.534707 INFO::Fitting model to feature number 845, X9fea56b0e2353c17ac863a7431b495fa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.572321 INFO::Fitting model to feature number 846, fd5bf8c9e268f6c5f48131c0615125fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.602912 INFO::Fitting model to feature number 847, X162a0a70867bfcb444146c5462dd083f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.639753 INFO::Fitting model to feature number 848, X2d8c6fa19ed17f298d900b8b94ca3577
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.670893 INFO::Fitting model to feature number 849, X9708c6393469d007561fa62858323ef3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.70081 INFO::Fitting model to feature number 850, X95312c6f0ad3ebf16fc431a571042056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.73107 INFO::Fitting model to feature number 851, X9cf77d8c73a399bf0ce345130e8448a6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.761842 INFO::Fitting model to feature number 852, e8f211866f017fe1a3cb035bbbf0f487
## 2026-03-04 18:36:37.791867 INFO::Fitting model to feature number 853, X28c021635cbdb3e420e03898e6736980
## 2026-03-04 18:36:37.821727 INFO::Fitting model to feature number 854, X06d2ce585d904c386410802b9ca7f6dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.852802 INFO::Fitting model to feature number 855, X96c5ae4d64ae76c5f2e1dedfb7ee89b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.886408 INFO::Fitting model to feature number 856, X263f6b9f13db9a695ba7215622ccf3e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.920897 INFO::Fitting model to feature number 857, e6c6d6c4f97b71d8b0747ea36185da7a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.952998 INFO::Fitting model to feature number 858, X8bb590cfc8fc74aad9d4f8a958a850f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:37.986541 INFO::Fitting model to feature number 859, X5d7946238ac7676f2e6c32ce5554a326
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.026962 INFO::Fitting model to feature number 860, X4dd5e98c725ec5f2fcc4a0c284b585e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.064976 INFO::Fitting model to feature number 861, X4aeeaaabe391e0fbeba72e87f6468d6d
## 2026-03-04 18:36:38.095893 INFO::Fitting model to feature number 862, X300433f8069de9c7d93f326c2bcb5e71
## 2026-03-04 18:36:38.125601 INFO::Fitting model to feature number 863, X036a46c77db16f0cab39f25993f32f5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.157834 INFO::Fitting model to feature number 864, X3bd6c8564a3f23aac3247f77acfd8adb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.189851 INFO::Fitting model to feature number 865, f8d3badaf1aa1b3e2e544de4310e176f
## 2026-03-04 18:36:38.220346 INFO::Fitting model to feature number 866, X28d6d19e890138b314ae244a6511662e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.252374 INFO::Fitting model to feature number 867, e54123bfdac1a2a0b1468c72d3a72056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.283947 INFO::Fitting model to feature number 868, X10507e5ecf030ee88fa122218bfa039e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.314356 INFO::Fitting model to feature number 869, X6fde2683f8248a303b9d6880e4a4e604
## 2026-03-04 18:36:38.344952 INFO::Fitting model to feature number 870, ba90f0b9771f6de4992c42c2004c1f75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.380445 INFO::Fitting model to feature number 871, a107e80c66b63e13999c5ca78894559c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.411633 INFO::Fitting model to feature number 872, e1f2d98fc9b94c8e80eb3124ab697199
## 2026-03-04 18:36:38.44344 INFO::Fitting model to feature number 873, X04d946612707d3fc0a92f32dcc8d1c30
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.47822 INFO::Fitting model to feature number 874, X33dab83ed982b34105927f635d5ce4d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.514908 INFO::Fitting model to feature number 875, X7feea0d71ff9eb9f5ea3d95e10d7b7b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.553164 INFO::Fitting model to feature number 876, X9eb591d492e323d974da25b0f092f73c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.584309 INFO::Fitting model to feature number 877, X0a8faf4be9b0ab6ca854de7a5bcad02e
## 2026-03-04 18:36:38.622817 INFO::Fitting model to feature number 878, X417159e7c56b61db9020926381287e6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.660983 INFO::Fitting model to feature number 879, X62871a4f2e689b3d2a8e24b4b6e90704
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.69937 INFO::Fitting model to feature number 880, X98345182f8f896dcba7b94aa00dd7eb7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.730473 INFO::Fitting model to feature number 881, X117ce9780a6794585b2b195cd40c66e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.761662 INFO::Fitting model to feature number 882, X1bf12663366670a5853f8ef1f662349b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.796476 INFO::Fitting model to feature number 883, e445d1b02329b8d17b6c33e47b2a1d35
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.828734 INFO::Fitting model to feature number 884, f2a3db16c28d4e5de91808d1c325f84e
## 2026-03-04 18:36:38.869895 INFO::Fitting model to feature number 885, X04646829704cf62ef942ccf8674fa3ce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.909473 INFO::Fitting model to feature number 886, X262b03eb315b00d94056179bb2f4ee53
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.949923 INFO::Fitting model to feature number 887, X446609f70ba4143bd817f44c210a0a4c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:38.981165 INFO::Fitting model to feature number 888, X794e9f08539a1c21267d2a49997fdddd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.020915 INFO::Fitting model to feature number 889, ac08cf8b04484d0a38afe5c6bc57420f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.05665 INFO::Fitting model to feature number 890, bba6e5f189ab17100885684b1c62d3e1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.090581 INFO::Fitting model to feature number 891, X76b202a36a7aa7bfa3b4f63e52e516cf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.126282 INFO::Fitting model to feature number 892, f0bc69d6d2e088f0e7895bbbe6de0d23
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.158122 INFO::Fitting model to feature number 893, X7b1effb26b50937a183a7a2f2bea3ae5
## 2026-03-04 18:36:39.190178 INFO::Fitting model to feature number 894, X9bea3a65d55fbec7b1de0fd8856c9bdb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.222745 INFO::Fitting model to feature number 895, d977c29f2689d12a8a72c44b02574e1c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.253537 INFO::Fitting model to feature number 896, a1ebedf636d67f3f4b9a40f02c5e702a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.284872 INFO::Fitting model to feature number 897, X4f6e37f2dfb5bcea2498189b425e0496
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.315286 INFO::Fitting model to feature number 898, X5d1bbb79e15c275d052188c91db4efce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.346751 INFO::Fitting model to feature number 899, a135c59c1c9deb0701d1f2eef52a8d34
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.37956 INFO::Fitting model to feature number 900, ba012602394a70abf571be8722eb7b33
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.41141 INFO::Fitting model to feature number 901, X9c8c8a3cfb2bab0fbc82c97e05bd537b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.442871 INFO::Fitting model to feature number 902, c843b591d20dcb1e33898346f82c79fa
## 2026-03-04 18:36:39.475636 INFO::Fitting model to feature number 903, d53cbb2969f0ef212c8cd9a67d1e0e1b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.507378 INFO::Fitting model to feature number 904, X143a46b132a920a6f068edbee257d9c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.538399 INFO::Fitting model to feature number 905, X1c401fa1144c19f7bb293293b446b27b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.568621 INFO::Fitting model to feature number 906, X4c880f5126065d48868d8720a7aeddd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.607348 INFO::Fitting model to feature number 907, X32359fd2ac112d65e1f7e9ae5ecd8b4a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.637712 INFO::Fitting model to feature number 908, fdcf22ebf60cd046b27e706f488c6c47
## 2026-03-04 18:36:39.674412 INFO::Fitting model to feature number 909, X55c9e708227432632f53f62dc5c8e7fb
## 2026-03-04 18:36:39.721809 INFO::Fitting model to feature number 910, a0bea49cb183c27b663831c30239c2bb
## 2026-03-04 18:36:39.758585 INFO::Fitting model to feature number 911, de7dd3152ac71f84ae29298ecd828f06
## 2026-03-04 18:36:39.789581 INFO::Fitting model to feature number 912, X528f95a03adde9b5d7b6cb2805b140d7
## 2026-03-04 18:36:39.820352 INFO::Fitting model to feature number 913, X361d196fabf4a7e0f2b0dcfcebe5320d
## 2026-03-04 18:36:39.857639 INFO::Fitting model to feature number 914, X79f55d42a25d89fdbe2f7c4665679358
## 2026-03-04 18:36:39.888298 INFO::Fitting model to feature number 915, X555423d48ff522cb887cbef0a2370ff8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:39.920193 INFO::Fitting model to feature number 916, X3dcbcb2abd5f2838e7f226e4d09a1f2f
## 2026-03-04 18:36:39.951011 INFO::Fitting model to feature number 917, X8312da8265d49aac51fdd0d54b32d4fd
## 2026-03-04 18:36:39.982056 INFO::Fitting model to feature number 918, a54c597c99769042fb4a11f69b70ac67
## 2026-03-04 18:36:40.013031 INFO::Fitting model to feature number 919, X5bc23fba732fd04d4986fd231269d6f3
## 2026-03-04 18:36:40.043574 INFO::Fitting model to feature number 920, X6c3a59bb3a352b8b87bd7a1109ec5548
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.076941 INFO::Fitting model to feature number 921, X53512fd8672ee91eef4eea370182e655
## 2026-03-04 18:36:40.107693 INFO::Fitting model to feature number 922, X5f4079304db67863bd1a742918ef8182
## 2026-03-04 18:36:40.137649 INFO::Fitting model to feature number 923, X3133ff9657a0f895b8cdd9b72fd49b05
## 2026-03-04 18:36:40.171357 INFO::Fitting model to feature number 924, X49bfc506d5b3743431bb0eeae6069917
## 2026-03-04 18:36:40.200938 INFO::Fitting model to feature number 925, c3daa62f8c41fb7c3a0d2ea92fc6b4af
## 2026-03-04 18:36:40.230513 INFO::Fitting model to feature number 926, cb7c601e1ff33d0a8e7ec6a8a7cb276b
## 2026-03-04 18:36:40.260125 INFO::Fitting model to feature number 927, X167e93c6b19b4454ef559b9fbfb409cb
## 2026-03-04 18:36:40.297431 INFO::Fitting model to feature number 928, X0b842f865339edc07fdf5ff3bd3f8379
## 2026-03-04 18:36:40.327624 INFO::Fitting model to feature number 929, X3dc628a32ca4dd564074799ece230f67
## 2026-03-04 18:36:40.358526 INFO::Fitting model to feature number 930, daabbbd5db10cf6cc8312d919bbcb924
## 2026-03-04 18:36:40.388416 INFO::Fitting model to feature number 931, f4f07e406c1e85a13f899ebf9b95e29d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.418568 INFO::Fitting model to feature number 932, X146d9408efd14ffff6e63556d93b5108
## 2026-03-04 18:36:40.4508 INFO::Fitting model to feature number 933, X6a5f53fae6f84847d447098ae4e8f98f
## 2026-03-04 18:36:40.4906 INFO::Fitting model to feature number 934, X1efbe8e6fc06fbbd97b5f4bfcf4a4ed0
## 2026-03-04 18:36:40.523655 INFO::Fitting model to feature number 935, a90a928ed331d7e3c701c511619e214d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.554679 INFO::Fitting model to feature number 936, df23731588f3884c39b70b5d6e7811bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.585474 INFO::Fitting model to feature number 937, fa602366e449f46a0d8cfaf64a48e353
## 2026-03-04 18:36:40.619614 INFO::Fitting model to feature number 938, b5832c657c49cbf026cd29bacca46fff
## 2026-03-04 18:36:40.650215 INFO::Fitting model to feature number 939, f55431826646e45d3a5a05104fda4c2c
## 2026-03-04 18:36:40.687508 INFO::Fitting model to feature number 940, a4c443d483ee55325f29b06c5d45e011
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.724799 INFO::Fitting model to feature number 941, a833bd321294a6cebf81e23cc647560b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.755877 INFO::Fitting model to feature number 942, ebb6d47df55f900c7264fde90ecc0085
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.793125 INFO::Fitting model to feature number 943, aa42a80bfe760ce5ef13e31a63ea9b2e
## 2026-03-04 18:36:40.823032 INFO::Fitting model to feature number 944, X2e3b65e7f787eb3a955b046b3a399821
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.859709 INFO::Fitting model to feature number 945, X976ebe9b8c37410b4afb4acb24205cb9
## 2026-03-04 18:36:40.890488 INFO::Fitting model to feature number 946, X2470670cc5f1cc0d346e708bd238260a
## 2026-03-04 18:36:40.926236 INFO::Fitting model to feature number 947, X90e925964999f0db59e87db25b2efa8f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:40.963568 INFO::Fitting model to feature number 948, X942a41f0147db28c830cb529bd03add1
## 2026-03-04 18:36:41.00061 INFO::Fitting model to feature number 949, X9a2da91ca6be4cdeed6cb6edb5f6b4ac
## 2026-03-04 18:36:41.035367 INFO::Fitting model to feature number 950, a9fa05aa77b1ee436a3f27c2856db7a7
## 2026-03-04 18:36:41.06981 INFO::Fitting model to feature number 951, dc126ebf3a1c182a034f5acbf4e22c14
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.102577 INFO::Fitting model to feature number 952, df76282db54adacc08160fc0eaeb6e7d
## 2026-03-04 18:36:41.139046 INFO::Fitting model to feature number 953, X271b453412999e72adabfd71eb789300
## 2026-03-04 18:36:41.174479 INFO::Fitting model to feature number 954, d36d14ae5dbdbd3d91e8d5fa3d5ec78d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.206751 INFO::Fitting model to feature number 955, e1b6cc549ae04f1c81727b330eb8bdb4
## 2026-03-04 18:36:41.238504 INFO::Fitting model to feature number 956, X6bc05eda0c85a2a0767635570b137a45
## 2026-03-04 18:36:41.270388 INFO::Fitting model to feature number 957, e09fe686537ebbc38c983c79018ffccf
## 2026-03-04 18:36:41.301449 INFO::Fitting model to feature number 958, X5d70be0acfe412990838d37001eb24b2
## 2026-03-04 18:36:41.344437 INFO::Fitting model to feature number 959, e4efa7d4be218d347e0d930295e832a3
## 2026-03-04 18:36:41.37613 INFO::Fitting model to feature number 960, X0dded72ceaf4dca73efed63cfb19812b
## 2026-03-04 18:36:41.407693 INFO::Fitting model to feature number 961, e1d4f791ba1642d0da77adc5de295325
## 2026-03-04 18:36:41.444933 INFO::Fitting model to feature number 962, f244f35600a551dfb3eeeaf86d473de8
## 2026-03-04 18:36:41.480209 INFO::Fitting model to feature number 963, X498291ee6f8d735d5f0353c70be3166c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.510928 INFO::Fitting model to feature number 964, X4e661aa07dec1515abb2615f59b388af
## 2026-03-04 18:36:41.541954 INFO::Fitting model to feature number 965, X50ec0d0bf1a352d385b6338057ac9bb9
## 2026-03-04 18:36:41.578711 INFO::Fitting model to feature number 966, X5a1b8ebc20cd26618173eca7a37685f1
## 2026-03-04 18:36:41.614716 INFO::Fitting model to feature number 967, f3beee66d95cd06a16b1919dde141bd0
## 2026-03-04 18:36:41.644736 INFO::Fitting model to feature number 968, X8c266bccda7a4c3fb93cad5f8e79fc08
## 2026-03-04 18:36:41.675255 INFO::Fitting model to feature number 969, f52ce8860f150e26d9ed45e91770ee38
## 2026-03-04 18:36:41.705933 INFO::Fitting model to feature number 970, d8e67bf70feea21e2daeddf4410172a1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.736856 INFO::Fitting model to feature number 971, X50fa0bc6d88a3894aead1fc528d41058
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.768083 INFO::Fitting model to feature number 972, c555bb81f49b80858ee015a9eea42aef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.798481 INFO::Fitting model to feature number 973, X91c734b3c610cf63571b5a22f5d4dc2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.827949 INFO::Fitting model to feature number 974, a14e904a7d2d1daace916fdf9856d31d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.858682 INFO::Fitting model to feature number 975, X94d7d95a02e3395be3ab28a0fed2872c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.889525 INFO::Fitting model to feature number 976, fb7f32b378722a5ca83e91cb41927eef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.919598 INFO::Fitting model to feature number 977, X9ebb4f9911373eb720c9bff4af833ddd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.950791 INFO::Fitting model to feature number 978, X667b35160a1f273c609d48a98ceb4d55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:41.981811 INFO::Fitting model to feature number 979, X6d2da75b05c499e6c1acb7351f0cf660
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.012882 INFO::Fitting model to feature number 980, ca9f3b34c45626302fa1fec4e105314e
## 2026-03-04 18:36:42.049932 INFO::Fitting model to feature number 981, X260de73cfc2af53494272275576eaa65
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.081489 INFO::Fitting model to feature number 982, e97000fb903a8e88351eb8e8d1c9d42c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.113327 INFO::Fitting model to feature number 983, X363e9dd3a114361345cb1ce36e8dd520
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.155926 INFO::Fitting model to feature number 984, X1e0f313eacdce85ff324dd5823aedff5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.188295 INFO::Fitting model to feature number 985, X8b147b77e5de5f3d3c1f6c478ffbab2a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.225447 INFO::Fitting model to feature number 986, X38d3539bc11b4977d1d3f32e1de0dbf4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.256832 INFO::Fitting model to feature number 987, a5f83ed4dc472953510501a5f899c2a1
## 2026-03-04 18:36:42.287216 INFO::Fitting model to feature number 988, X5cc9033bb1a702b19b9a7b2e44ecd851
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.324836 INFO::Fitting model to feature number 989, X3b28c382abaa17c7f1f339967619ddd4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.355833 INFO::Fitting model to feature number 990, bca782ea0c57a7716ed3aea3f5376da9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.385926 INFO::Fitting model to feature number 991, X2a57beb2f11af3ae774c0df31aac1275
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.416744 INFO::Fitting model to feature number 992, X5a3930b7b374d9d2019a9b26531ffb6d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.448471 INFO::Fitting model to feature number 993, X7fc2f360dc73536bb0281761c75399d5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.478759 INFO::Fitting model to feature number 994, X5d2b45d606b5eafbcc1e9f87452d3135
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.516374 INFO::Fitting model to feature number 995, X66deedccb8ee4f471c0029a94b18ec49
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.547345 INFO::Fitting model to feature number 996, X4aa6ba9c7da5f33af213c02d65ecb391
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.578466 INFO::Fitting model to feature number 997, X8da71749ef2ac23a056fc49fdda9131e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.612908 INFO::Fitting model to feature number 998, c90f38304baf5fc52b3332e63a8016c2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.643665 INFO::Fitting model to feature number 999, X396daf5d94b484aba7fcc4d8c0223759
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.674355 INFO::Fitting model to feature number 1000, X7bb6e29f06620ad6f4f0e7478a1fea5f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.704718 INFO::Fitting model to feature number 1001, X8463aae5183f78a44bcd2369cdf5ed93
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.736485 INFO::Fitting model to feature number 1002, X8aa5eb1b9d4e710d35223d1ddaaccc22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.768543 INFO::Fitting model to feature number 1003, X8fb79654f3000a52e5b1be221f49106e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.806097 INFO::Fitting model to feature number 1004, ca76a5d3cef024cc5c9c1ea5b19d77fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.838094 INFO::Fitting model to feature number 1005, X1aae042eddcd396bbe739d785fb9df7e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.869029 INFO::Fitting model to feature number 1006, X2b52a41587c8ed22703a2bc2d978f584
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.899991 INFO::Fitting model to feature number 1007, X389f789d4f3fae4c7040028704fbbc03
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.934631 INFO::Fitting model to feature number 1008, X41d3c7fe72f0ac46769a543a2322c72b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:42.978191 INFO::Fitting model to feature number 1009, ba4cf501f96340f482a9525b4ab1f69d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.009824 INFO::Fitting model to feature number 1010, X2a4d95ff9ae281d8c13c486f3772d1bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.0397 INFO::Fitting model to feature number 1011, X5f05f4bb091a1d67b7c724a6aaffdcc4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.074722 INFO::Fitting model to feature number 1012, X7705f0fc350c0d3189688a4858948d1c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.10675 INFO::Fitting model to feature number 1013, a79885cf6c4f4e70eb18d351e69b8b78
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.141555 INFO::Fitting model to feature number 1014, a907f48fc8e04ed9ae516022117a446d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.17534 INFO::Fitting model to feature number 1015, ad901777a802617307dca87d69e1f57d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.207366 INFO::Fitting model to feature number 1016, X751ccf1e1783168abe381e454cd7f817
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.238686 INFO::Fitting model to feature number 1017, f4216c4a78111523f17cbc7de20025b6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.269829 INFO::Fitting model to feature number 1018, X3528210d91f62915dfd614d5f8c5246d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.307175 INFO::Fitting model to feature number 1019, c1e2db64acc21a102c6a7d833d730a67
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.338507 INFO::Fitting model to feature number 1020, X48e84a9b22d27e8f7b2854e2d0e6e151
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.376342 INFO::Fitting model to feature number 1021, ff2a02b6f268f3d8c34998b06e396f46
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.413411 INFO::Fitting model to feature number 1022, cc55eefd97b5534a346ae449625be5dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.445195 INFO::Fitting model to feature number 1023, X50e17521cf480fe34b2801f90cee6178
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.476948 INFO::Fitting model to feature number 1024, X6042f5b4bcfc4ca01bf9dba3a94f82de
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.512505 INFO::Fitting model to feature number 1025, X6079f6c9472d90f3fe6b93c2786229df
## 2026-03-04 18:36:43.542842 INFO::Fitting model to feature number 1026, a11644f2e991a5a53bc43924a8e0f75f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.573425 INFO::Fitting model to feature number 1027, cb5f94b1638b84a247685c32724e4b10
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.611304 INFO::Fitting model to feature number 1028, X2faf5487b3d0d36915ccbac8006e96cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.64793 INFO::Fitting model to feature number 1029, X961327649d0d8b1a3f1d8506f7abf9a3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.679947 INFO::Fitting model to feature number 1030, X8e867981ffca8750b7a43cb387f9f135
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.711475 INFO::Fitting model to feature number 1031, X35b607864efdc9e272585bca1d40016e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.742972 INFO::Fitting model to feature number 1032, X1a074fc9ba7138aca27f372e1b4a000e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.77467 INFO::Fitting model to feature number 1033, X92e8b9e50e7601db88a32a3b69963bb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.817789 INFO::Fitting model to feature number 1034, X45a2ac2f0e2a4c246a56beecda0e2dbd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.848935 INFO::Fitting model to feature number 1035, fab3fc1d161c45c577a91e58006417b1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.879548 INFO::Fitting model to feature number 1036, X1cffb93453b8e2a2a884422848606aa6
## 2026-03-04 18:36:43.909382 INFO::Fitting model to feature number 1037, X84540ed2eebf1efb8e0d8d1e40dabffd
## 2026-03-04 18:36:43.939373 INFO::Fitting model to feature number 1038, d012209e84dc07f3d3aeccff20da3d51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:43.973533 INFO::Fitting model to feature number 1039, X50a43b2ec96d6b8bb0f9a8680a05f379
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.004472 INFO::Fitting model to feature number 1040, X567a407c520a5cf74f4cd9bd5e13b199
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.037335 INFO::Fitting model to feature number 1041, X66c3e26fde9a995062499289924775a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.07059 INFO::Fitting model to feature number 1042, X7194882e587298e6c178faf29b1572bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.10869 INFO::Fitting model to feature number 1043, d75f983da6423ae3900c33e44a566605
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.140962 INFO::Fitting model to feature number 1044, X89651bf4f451a1f07f6ab1983f484fce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.179025 INFO::Fitting model to feature number 1045, X3d9449debdcefa1ef7caea17e9aee6d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.216694 INFO::Fitting model to feature number 1046, X58a5ed3a034ea45c3da6ccdfae4f3659
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.250611 INFO::Fitting model to feature number 1047, X467bd9bb0da5518cccb8ef64d3427bd9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.284506 INFO::Fitting model to feature number 1048, X4ba356cfaa02c5d2668304cf8a75d98e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.316695 INFO::Fitting model to feature number 1049, X59f2c893c976b2bed0ce5f3168b51b5e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.354101 INFO::Fitting model to feature number 1050, d49581d80b688f6ea99b76cae93334ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.385084 INFO::Fitting model to feature number 1051, e0abea0159042cd283f2560997bd6dd7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.421955 INFO::Fitting model to feature number 1052, fc77df84a1ad92454b15a3a729ca9fab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.456466 INFO::Fitting model to feature number 1053, X5de6b44b38194149036db56946e96b38
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.488456 INFO::Fitting model to feature number 1054, X883c6c1b596bf1c5cf45a00427cebbde
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.519902 INFO::Fitting model to feature number 1055, be1fe3655bebd478e3fe789cc433ec14
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.55172 INFO::Fitting model to feature number 1056, c7aa42e4a69d17cdb7a5374e8c510136
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.582439 INFO::Fitting model to feature number 1057, f3f9d3323e4d2f15298aaed5439af56a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.614018 INFO::Fitting model to feature number 1058, X765c815c5453e9ed5dc420c8119c5d01
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.657029 INFO::Fitting model to feature number 1059, d36986d02dbdc4507e68cd32d674ac05
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.689158 INFO::Fitting model to feature number 1060, d4c443b446021feb34aad85d999b11cb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.727453 INFO::Fitting model to feature number 1061, X08b4991c00fa39be14bc8f7327e213a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.75929 INFO::Fitting model to feature number 1062, X4d076914382017f211091f07aa8db632
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.791484 INFO::Fitting model to feature number 1063, X944c1fe867c513e7c1e38241d24352c1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.824502 INFO::Fitting model to feature number 1064, b7af501666cfd4d95b15a259a0d33534
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.857266 INFO::Fitting model to feature number 1065, d49f0e0482144b9ac84c3a0dad81f3dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.895381 INFO::Fitting model to feature number 1066, ea17d1fba365a308f925eb11b72264bf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.926808 INFO::Fitting model to feature number 1067, X53f66ce1db04a14261e921e5f8dd43db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:44.965015 INFO::Fitting model to feature number 1068, X54ea53f5747082b754171686efe7ae35
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.005058 INFO::Fitting model to feature number 1069, X93628c71a8cacf127ce7fae55d40c1d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.038189 INFO::Fitting model to feature number 1070, a0d658715dbad152018ebec95f5b2a13
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.071343 INFO::Fitting model to feature number 1071, faeb2051952acacd1ee1deff131c3bb0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.102848 INFO::Fitting model to feature number 1072, X58d7505e7eb30d2e9a47244a921e3fae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.134618 INFO::Fitting model to feature number 1073, X85b0ab39a82b695768422d01e88a25d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.173006 INFO::Fitting model to feature number 1074, cbb718ece292de4ff64d526aa5604372
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.205878 INFO::Fitting model to feature number 1075, X5ab675c9e714d085b6b113d8892cd685
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.237291 INFO::Fitting model to feature number 1076, X46b863013a51fc2faa3e8cbf8c77271a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.267727 INFO::Fitting model to feature number 1077, cd2bb0b70d041ae5bb64067b2ab990ec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.298158 INFO::Fitting model to feature number 1078, b1c849febb33944267c7047a39244974
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.33439 INFO::Fitting model to feature number 1079, X5749de9d92f48616e003d42550b8f79d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.365167 INFO::Fitting model to feature number 1080, X51ad6d6b27b7f1551bd818f9f253fb61
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.396144 INFO::Fitting model to feature number 1081, e78437188624ce623ae9bcd812e72d2b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.427118 INFO::Fitting model to feature number 1082, X0baa3c3e2e6042eede193ac9f25d53a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.467573 INFO::Fitting model to feature number 1083, X9c3e3ed907988deb778c5f0cdec6de0a
## 2026-03-04 18:36:45.499997 INFO::Fitting model to feature number 1084, X0f5e024743e654ca7872273094fbed5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.537108 INFO::Fitting model to feature number 1085, X784e1de605fa06ef3726907163c5027e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.568201 INFO::Fitting model to feature number 1086, X9e5b25a18c0bf19a5b288654dc4ea008
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.599472 INFO::Fitting model to feature number 1087, X1004ca7e65c7a7d359941799d33a67b3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.630583 INFO::Fitting model to feature number 1088, X8f57a9d96db41075fa1ba4c4328e26d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.661053 INFO::Fitting model to feature number 1089, X04ede807f45e5ac5b121a53a237f830b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.691668 INFO::Fitting model to feature number 1090, X3ead1fa5a380796e8a1f7c74c5212bc0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.724115 INFO::Fitting model to feature number 1091, d8edc720feb6c43730bcfe12a789b5e8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.761636 INFO::Fitting model to feature number 1092, X5c10c0918fb25e77e97086722698f1ba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.797764 INFO::Fitting model to feature number 1093, X2d84007ca1c2c3f513156f6f520716fc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.828058 INFO::Fitting model to feature number 1094, X9f60b24425f186b088fe6280c81a8388
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.858727 INFO::Fitting model to feature number 1095, cfc746d1cf45fc692c0e53c4c7f671d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.889939 INFO::Fitting model to feature number 1096, d874902a8d6182ca9798bc486ebb4a99
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.924102 INFO::Fitting model to feature number 1097, f392db5a807da4e8a2e708bd84e04ca0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.958309 INFO::Fitting model to feature number 1098, X6d3c164b889d28be2ec50d08fb62322e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:45.988493 INFO::Fitting model to feature number 1099, c18ca53cea68a2a7144fe36bce89340b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.024587 INFO::Fitting model to feature number 1100, ce68fc44b24d1212ae41371700862b72
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.054772 INFO::Fitting model to feature number 1101, X455814e0df838e53217a329758d86744
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.092107 INFO::Fitting model to feature number 1102, X68612badfa6c422975b70cf092a1b873
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.122408 INFO::Fitting model to feature number 1103, X7d5e0a18434e6d972b5c2dc28b8bf51d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.153005 INFO::Fitting model to feature number 1104, X86a4904762e51b23c5eac255935f0534
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.190094 INFO::Fitting model to feature number 1105, X9060dd97d853d275ec247f2db98c265f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.219904 INFO::Fitting model to feature number 1106, X49b65a3ca234740512a8dc5b92ccf5a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.250133 INFO::Fitting model to feature number 1107, X4f9ed6ab77896830974fd4df52525f76
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.291461 INFO::Fitting model to feature number 1108, X674f8b7dba9f080323fcec4292eb16e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.329195 INFO::Fitting model to feature number 1109, X1dbdd4ba3c7ae9b85cea8f0f4cfa8a10
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.35993 INFO::Fitting model to feature number 1110, X9590c032d57e4f5568b0eafa67708b3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.390949 INFO::Fitting model to feature number 1111, b39e9d4cb76e6b56dc9999a92c5283d4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.421931 INFO::Fitting model to feature number 1112, e978d95f90a6f389da9b57cf0868b847
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.452502 INFO::Fitting model to feature number 1113, X198065e0411746d3c31cbc810a1c8b6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.483859 INFO::Fitting model to feature number 1114, bc25945dde85e8f2f55c642669695015
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.514507 INFO::Fitting model to feature number 1115, ce3868416d546ecf20e2c204da9a1828
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.551078 INFO::Fitting model to feature number 1116, X8910dfaae22f25dfe7f0679ca7c94618
## 2026-03-04 18:36:46.583306 INFO::Fitting model to feature number 1117, e39c0f2f2bf6163e4f2da67826b09132
## 2026-03-04 18:36:46.618909 INFO::Fitting model to feature number 1118, X7ad2b0d2a4533eb2ce608b7852890f98
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.651386 INFO::Fitting model to feature number 1119, X597728f24a67cc2770bedadf83591fdd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.691279 INFO::Fitting model to feature number 1120, X8cfb2b6f68d3ac591ad8d67f325aecc5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.73215 INFO::Fitting model to feature number 1121, a48784576958e537b25f9fa5d4e3f98e
## 2026-03-04 18:36:46.767482 INFO::Fitting model to feature number 1122, X267edef74cf84d7548b0d3f88cb5ff89
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.799449 INFO::Fitting model to feature number 1123, cc9703e82e512181f882f6bf229cb138
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.831021 INFO::Fitting model to feature number 1124, d91a7f73b26a7d13a51d9f3e4fe3cc06
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.86313 INFO::Fitting model to feature number 1125, X6331ea1bb0c3c736beaf34fd6284fd31
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.900903 INFO::Fitting model to feature number 1126, X87eb2526a517ff1a76336ae23315c926
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.934965 INFO::Fitting model to feature number 1127, X2e27f1f31ebbd4b4c0ff20a7feb121a8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:46.972001 INFO::Fitting model to feature number 1128, bcfe7158c79ee0946dc47793ffbf0549
## 2026-03-04 18:36:47.003472 INFO::Fitting model to feature number 1129, X27a395b590a278e9085d7157656a37ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.035329 INFO::Fitting model to feature number 1130, f6fd2a76558c3266f11fe4c856f929e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.069027 INFO::Fitting model to feature number 1131, afde775f8b52d66edb513773b710a8e0
## 2026-03-04 18:36:47.102032 INFO::Fitting model to feature number 1132, X9939ab0bed6fd6fa3c489cfa467cf9c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.144112 INFO::Fitting model to feature number 1133, X11f98106eb45262aa2b864ed085f3369
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.175786 INFO::Fitting model to feature number 1134, X19609c541e304d9414eeed6a276f038a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.206588 INFO::Fitting model to feature number 1135, X33151a1928fbfcff924af5d627bc553b
## 2026-03-04 18:36:47.237727 INFO::Fitting model to feature number 1136, a2224b2e402597ee34ee3fdcd9fb7298
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.274511 INFO::Fitting model to feature number 1137, X29087b19e0622aa58bc91d8e521126e4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.310278 INFO::Fitting model to feature number 1138, fcde5c46c10eec800e6b38a39ea9e0b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.344647 INFO::Fitting model to feature number 1139, X219c33e8cda1089efcbe919a3694dbf6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.378125 INFO::Fitting model to feature number 1140, X29d5f5579a251818c2a9c49b911782ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.409486 INFO::Fitting model to feature number 1141, X6e72215f2c39c515db5fa57a00ba473b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.441682 INFO::Fitting model to feature number 1142, a6583feb07439f47775cf6216ebe3e9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.473177 INFO::Fitting model to feature number 1143, X2d5028d34b33f496c8ea9d664e3d9f24
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.505272 INFO::Fitting model to feature number 1144, X9c0ac0fed41ad08dcd91f61b0440f4bf
## 2026-03-04 18:36:47.540661 INFO::Fitting model to feature number 1145, ff5ab5918e8b9c763a2f1c3cf33d06f1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.572504 INFO::Fitting model to feature number 1146, X67e28f844816795eb9d42efcff09c29b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.603226 INFO::Fitting model to feature number 1147, X79fdee6143db1981342333ae4a977269
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.633748 INFO::Fitting model to feature number 1148, f5aa4ea6105b38599df1b4c9f2c526ab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.664912 INFO::Fitting model to feature number 1149, d315b9630a2ca429ed87f6ac6856f505
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.699861 INFO::Fitting model to feature number 1150, X0186f5ab59b4d03b51ca5bb5d69fec54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.73082 INFO::Fitting model to feature number 1151, X334079baf11f432d414d30bd1e0f09e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.764255 INFO::Fitting model to feature number 1152, X960304d4b9cb07256cceb3a66af3234f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.797048 INFO::Fitting model to feature number 1153, X9385f532468f13f7aeef2e694329d709
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.829476 INFO::Fitting model to feature number 1154, f5a9b9639b6ba9307b5166734c9a69fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.861332 INFO::Fitting model to feature number 1155, X26c053114523080cbc3b791b07e8b937
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.891855 INFO::Fitting model to feature number 1156, X54d7f7e9bf7bd4e00d059e3a9f80b832
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.922947 INFO::Fitting model to feature number 1157, X91237a8dd758dacb643569fe824f48b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:47.972859 INFO::Fitting model to feature number 1158, ca7f6db958f83d9a7daac017d413bcc8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.004766 INFO::Fitting model to feature number 1159, X586c90b4d2a884eeb3bbdfb2c129e908
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.038597 INFO::Fitting model to feature number 1160, X71ed20fd789eccdde477d5b95c58ce98
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.084607 INFO::Fitting model to feature number 1161, X76f64e5e9cbda0476bbb42ad34e4cdf9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.120519 INFO::Fitting model to feature number 1162, c46904bcacd273e0a9a37384d6e45a26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.154005 INFO::Fitting model to feature number 1163, cfcd2110898f3d3051894f249b744767
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.184202 INFO::Fitting model to feature number 1164, X10096fb55ea863485a884e2bf9542163
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.21787 INFO::Fitting model to feature number 1165, c7e8689997df61fe02151acdef28ef63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.254623 INFO::Fitting model to feature number 1166, dac90bcd762785b630843d62a7da6c9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.286617 INFO::Fitting model to feature number 1167, eed08fcc06d959756c2dd08f4ae201c9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.323666 INFO::Fitting model to feature number 1168, f5d15fd366b9bb809bbf1a53a126b694
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.355826 INFO::Fitting model to feature number 1169, X045f2c5846fb21d9f2dfbc9ec3bf64e1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.388173 INFO::Fitting model to feature number 1170, X1cb08a81bec75cf1cb6a5347637f32e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.421686 INFO::Fitting model to feature number 1171, X393a57d6c6b02ca59aa23b7a88a4b773
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.455778 INFO::Fitting model to feature number 1172, X478a6c05d04db6a58eaae7bb27eb04ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.488052 INFO::Fitting model to feature number 1173, X522fc4ccc38089e8f0d8d4a57d4e82f8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.519351 INFO::Fitting model to feature number 1174, b909bd96ae3a423bbef71d7add0c86df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.554506 INFO::Fitting model to feature number 1175, c937f4adda1f009333dd511ae285fff7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.585697 INFO::Fitting model to feature number 1176, d11d08c6a579ddd3cc910338ebbe6547
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.622553 INFO::Fitting model to feature number 1177, ec054373da44610fd5b9a258e5f793a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.656082 INFO::Fitting model to feature number 1178, f0af9b7703a95bfa52c7ebcfc5a94b82
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.689117 INFO::Fitting model to feature number 1179, f1f7e45805f74287c1564249d80cafa2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.724379 INFO::Fitting model to feature number 1180, f70c8600851d3d162e26ad9cf59675fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.757795 INFO::Fitting model to feature number 1181, f81c4bee45cfd34ed1fe726d0c844fb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.789701 INFO::Fitting model to feature number 1182, X3a323cf260c90c3306639d009f11b7fc
## 2026-03-04 18:36:48.831927 INFO::Fitting model to feature number 1183, b7513d5b69299f329da60ff370d4acf9
## 2026-03-04 18:36:48.865344 INFO::Fitting model to feature number 1184, X086c111d7310df390ca2ce0ac8144156
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.89691 INFO::Fitting model to feature number 1185, X462983a97faff888d34ac78f5f92825f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:48.927494 INFO::Fitting model to feature number 1186, X3a50580ca187917dbd5d8dd58f5c6a48
## 2026-03-04 18:36:48.96416 INFO::Fitting model to feature number 1187, X62ebbf17f633335b2fc737ebb581df16
## 2026-03-04 18:36:48.997968 INFO::Fitting model to feature number 1188, X727da14bd94f15b7dee15dbe146c0da9
## 2026-03-04 18:36:49.029815 INFO::Fitting model to feature number 1189, X3e737e299b0235bc5df1bc95b9d551cc
## 2026-03-04 18:36:49.061027 INFO::Fitting model to feature number 1190, X6680fc094b1ae654169e3f6596065387
## 2026-03-04 18:36:49.09803 INFO::Fitting model to feature number 1191, X54339f8b1100dc127f86478b017555bb
## 2026-03-04 18:36:49.137512 INFO::Fitting model to feature number 1192, X609000883c374a04e356bda18b419104
## 2026-03-04 18:36:49.169367 INFO::Fitting model to feature number 1193, X95b88fdb11b48228c744237ea2277420
## 2026-03-04 18:36:49.199614 INFO::Fitting model to feature number 1194, X9e2a90893c63429f7974a93724f58473
## 2026-03-04 18:36:49.235108 INFO::Fitting model to feature number 1195, X3646994e928de61adcc5460ca0df9cd0
## 2026-03-04 18:36:49.266756 INFO::Fitting model to feature number 1196, X4dbb7bbc94b14db1ffad5c6e3de8341c
## 2026-03-04 18:36:49.298136 INFO::Fitting model to feature number 1197, f521618f99a90bb8c936cb315d5b33cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:49.329872 INFO::Fitting model to feature number 1198, a473147383e13716f9365ff1926aa6a6
## 2026-03-04 18:36:49.361721 INFO::Fitting model to feature number 1199, f2c1249d1644237604788e8198b48d31
## 2026-03-04 18:36:49.393085 INFO::Fitting model to feature number 1200, fe3bfbd4d173653203a9d6d256d92240
## 2026-03-04 18:36:49.424055 INFO::Fitting model to feature number 1201, X077b2ad66183e1ab78bd0a27b5170970
## 2026-03-04 18:36:49.456155 INFO::Fitting model to feature number 1202, X39a893e34f9c30e78354a90670c464ff
## 2026-03-04 18:36:49.493633 INFO::Fitting model to feature number 1203, X89c032289bb7d50431817cd8ad8fb3e3
## 2026-03-04 18:36:49.525938 INFO::Fitting model to feature number 1204, X65970e707cca45512e17863477cbf9fd
## 2026-03-04 18:36:49.558385 INFO::Fitting model to feature number 1205, X9e7f26a7544795db8dd3f0dcad58d1d9
## 2026-03-04 18:36:49.591744 INFO::Fitting model to feature number 1206, X738e89f33bb2d9a90cc65b6e38395094
## 2026-03-04 18:36:49.636856 INFO::Fitting model to feature number 1207, e8d8f54456fddc89b8d3d072dfd0dae2
## 2026-03-04 18:36:49.668966 INFO::Fitting model to feature number 1208, X8e068f286e4459a2f92d83bf8fee4c75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:49.705016 INFO::Fitting model to feature number 1209, ff3ef656c9c5c1e8d9c808185a095229
## 2026-03-04 18:36:49.738843 INFO::Fitting model to feature number 1210, X44816a4453af39d938da26fe98668bf6
## 2026-03-04 18:36:49.769548 INFO::Fitting model to feature number 1211, dced475a8c58af383bb71156ede21a28
## 2026-03-04 18:36:49.799158 INFO::Fitting model to feature number 1212, X1fb31804045eed5be6772950d92bbea1
## 2026-03-04 18:36:49.82895 INFO::Fitting model to feature number 1213, X2d4eff8bb4b08b2ad545641114a1f95d
## 2026-03-04 18:36:49.858691 INFO::Fitting model to feature number 1214, X9b1e1af81c490f7167ba1f0f3a5c27f0
## 2026-03-04 18:36:49.888424 INFO::Fitting model to feature number 1215, X64c9855e3fb0a940f7160dc982b493aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:49.920467 INFO::Fitting model to feature number 1216, X8375b729fbec2ab3e94c4e34a6b89601
## 2026-03-04 18:36:49.951704 INFO::Fitting model to feature number 1217, adfe3978e202f9a097a1f8cce2fd35b7
## 2026-03-04 18:36:49.990769 INFO::Fitting model to feature number 1218, c52e54c93b3cc50a746340ebdd80aa55
## 2026-03-04 18:36:50.023987 INFO::Fitting model to feature number 1219, X5b1b6becbd52352e6a9618ddadf0d8d9
## 2026-03-04 18:36:50.056181 INFO::Fitting model to feature number 1220, a8fda41c9782acfa04a7d7df17761e59
## 2026-03-04 18:36:50.087485 INFO::Fitting model to feature number 1221, X4e4884420b31400d35db2a1417fa4ec7
## 2026-03-04 18:36:50.124577 INFO::Fitting model to feature number 1222, f236ac69ec7755b0247018f9b997b504
## 2026-03-04 18:36:50.157792 INFO::Fitting model to feature number 1223, X5bcdc43e09669d72979a94199b86a102
## 2026-03-04 18:36:50.191112 INFO::Fitting model to feature number 1224, dd0a2ab3c468bf12f8a8de815ab5dd63
## 2026-03-04 18:36:50.223876 INFO::Fitting model to feature number 1225, X9de219990e1666e0ba94f08161b49e4c
## 2026-03-04 18:36:50.255462 INFO::Fitting model to feature number 1226, ac6cc9a3a10683c803042c539a9caea1
## 2026-03-04 18:36:50.286756 INFO::Fitting model to feature number 1227, b683bb4bfec4e83e90550a48ad78ef9a
## 2026-03-04 18:36:50.318205 INFO::Fitting model to feature number 1228, e720001f2f9616203ce0cd35245a6218
## 2026-03-04 18:36:50.348063 INFO::Fitting model to feature number 1229, f50f30b5af2b1813e816f5aeedb1f22f
## 2026-03-04 18:36:50.384862 INFO::Fitting model to feature number 1230, X596e2f7be3a79400444a18242148fe31
## 2026-03-04 18:36:50.421008 INFO::Fitting model to feature number 1231, X8ccc17b289af3286eacf8da751dd6ffc
## 2026-03-04 18:36:50.46751 INFO::Fitting model to feature number 1232, fd7a5206567863e6bbc88ec1431ce83a
## 2026-03-04 18:36:50.499418 INFO::Fitting model to feature number 1233, X58448ec65533e295cd04dcfce9c81137
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:50.537444 INFO::Fitting model to feature number 1234, d6960b8def5f8b8cd6e93a90054cd2cb
## 2026-03-04 18:36:50.568097 INFO::Fitting model to feature number 1235, X083129e51191b2bf89f87bfafec40b48
## 2026-03-04 18:36:50.599233 INFO::Fitting model to feature number 1236, d923157d249ba936676cb5b0e460ba21
## 2026-03-04 18:36:50.630199 INFO::Fitting model to feature number 1237, X8b00ab07a4a5072a42e9ab7ab540313a
## 2026-03-04 18:36:50.664493 INFO::Fitting model to feature number 1238, f2480634f9b255a885f0577d0c84663d
## 2026-03-04 18:36:50.70195 INFO::Fitting model to feature number 1239, X2ce08b3a482c56941e6838df1a1d0db3
## 2026-03-04 18:36:50.739146 INFO::Fitting model to feature number 1240, X0c6a4054470a42e000bb0b69e5aa1a4e
## 2026-03-04 18:36:50.770437 INFO::Fitting model to feature number 1241, f1d223031e540e85dc5b60c84bbf80c6
## 2026-03-04 18:36:50.801911 INFO::Fitting model to feature number 1242, d52068433df4203a90b8e9da9e6e4dd2
## 2026-03-04 18:36:50.832412 INFO::Fitting model to feature number 1243, X1244fb68a2ab9cacab3a01719181ef35
## 2026-03-04 18:36:50.869986 INFO::Fitting model to feature number 1244, X25410db3bc0d1fa7ca52ce54385cb0af
## 2026-03-04 18:36:50.900387 INFO::Fitting model to feature number 1245, d11910302f2ea2a2dbca76d22dc4bcae
## 2026-03-04 18:36:50.936981 INFO::Fitting model to feature number 1246, X838980cc7edfa38dbcd8a24eeb97f9f0
## 2026-03-04 18:36:50.967012 INFO::Fitting model to feature number 1247, X19855f95f67ddb3a43733ce20d0e6d62
## 2026-03-04 18:36:50.997342 INFO::Fitting model to feature number 1248, X1dee71f32a863cdbaac0e10e0c89d5ec
## 2026-03-04 18:36:51.028839 INFO::Fitting model to feature number 1249, X618fd6e0f9e5467a143b628465cd3c8b
## 2026-03-04 18:36:51.06049 INFO::Fitting model to feature number 1250, X45a02d8c560f708f60e2eb952d45974d
## 2026-03-04 18:36:51.090093 INFO::Fitting model to feature number 1251, X820869c7c18af78f1bfe8f498e88e952
## 2026-03-04 18:36:51.119881 INFO::Fitting model to feature number 1252, dfaeccc139705834e879974f5b88c539
## 2026-03-04 18:36:51.155929 INFO::Fitting model to feature number 1253, X21cd87ba0ed75295c17dd698ef342af5
## 2026-03-04 18:36:51.186429 INFO::Fitting model to feature number 1254, X3c45cc52efc8e822effcb8d8fbec21d9
## 2026-03-04 18:36:51.218685 INFO::Fitting model to feature number 1255, X56e95bebfea3a564631df46966c413f8
## 2026-03-04 18:36:51.248693 INFO::Fitting model to feature number 1256, X1f12408f6bee875e10fcacbb5869bad0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.296964 INFO::Fitting model to feature number 1257, X5586fc56392c9262889f6ed54a2be971
## 2026-03-04 18:36:51.327495 INFO::Fitting model to feature number 1258, X7aaec5dc635a1d42d497edfe99e32c73
## 2026-03-04 18:36:51.363107 INFO::Fitting model to feature number 1259, X9b77571525020f2fb7a3d30e0161ce49
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.399066 INFO::Fitting model to feature number 1260, X9e3e8224c26b1c7c08dfa6e01de4ca5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.428721 INFO::Fitting model to feature number 1261, c886b40eb27052d8337856e81521ffbe
## 2026-03-04 18:36:51.459066 INFO::Fitting model to feature number 1262, e005c3a476aa46f8a7f918f704ef9cb6
## 2026-03-04 18:36:51.490487 INFO::Fitting model to feature number 1263, X13ca585ef9477f5bca479ffd2848097c
## 2026-03-04 18:36:51.519655 INFO::Fitting model to feature number 1264, X2242fe61584902366620d2a2096534f2
## 2026-03-04 18:36:51.550121 INFO::Fitting model to feature number 1265, X8f6c96e5604e4f77d84378b3d03782d0
## 2026-03-04 18:36:51.586628 INFO::Fitting model to feature number 1266, b71116fa177f4b07f53bf5d7996356a4
## 2026-03-04 18:36:51.623603 INFO::Fitting model to feature number 1267, e069f67ca919b8c627480de636066f4c
## 2026-03-04 18:36:51.653893 INFO::Fitting model to feature number 1268, X450cceb1cd3f129094c17ac29feb8d05
## 2026-03-04 18:36:51.688564 INFO::Fitting model to feature number 1269, X8ab70daaf7ddf3ee4baacee5f8c12693
## 2026-03-04 18:36:51.719604 INFO::Fitting model to feature number 1270, a0b24ff1cf7c2fc50822ba074982bfdd
## 2026-03-04 18:36:51.757505 INFO::Fitting model to feature number 1271, f11d3a89c1dae0dcdf1ff1e444155592
## 2026-03-04 18:36:51.787498 INFO::Fitting model to feature number 1272, e14d6eef3ea1efe7f99f0d1cfec2af99
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.81773 INFO::Fitting model to feature number 1273, X08030887b139d81f555f8f8352faf277
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.848225 INFO::Fitting model to feature number 1274, X1125d77abfbe71bd6698e9d969084307
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.87912 INFO::Fitting model to feature number 1275, X5862b017b7c3a0d098e0212f786ab5f0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.910039 INFO::Fitting model to feature number 1276, X7e3607fbbecab3ea931b4089a87ef58e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.941451 INFO::Fitting model to feature number 1277, c3027c24f295e7851649e1cd6e8b5013
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:51.979337 INFO::Fitting model to feature number 1278, e21e8275b16b5eb4781c76c5711b664c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.017871 INFO::Fitting model to feature number 1279, f222d3ac2aa36e3b93a61fcb909dcb28
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.051887 INFO::Fitting model to feature number 1280, X360335c2e90ac4b8274ba79192d02952
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.090113 INFO::Fitting model to feature number 1281, d0cf180c4dbb461db1db9e31f98e9017
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.138764 INFO::Fitting model to feature number 1282, b448327fe9304f5a67c5e8d7e63b47af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.173685 INFO::Fitting model to feature number 1283, X265d49e10b1beb36a822e75f8499e466
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.208161 INFO::Fitting model to feature number 1284, c86f8eb8159950eceaddf95a4ab5b57c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.238833 INFO::Fitting model to feature number 1285, ad7793e14bb253aa0e97e853e8b015f2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.271366 INFO::Fitting model to feature number 1286, bcc404be7315530aea16b5a5f3d1c976
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.303262 INFO::Fitting model to feature number 1287, X5f5ee11493e3bb7a92450ba053af4591
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.334204 INFO::Fitting model to feature number 1288, X12b74af5dfbb936b04b4e8852792a241
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.367202 INFO::Fitting model to feature number 1289, X71e1a15d4596752eaeba8282f55eda63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.40026 INFO::Fitting model to feature number 1290, a5638847d1da6845322b4c096586adf3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.438191 INFO::Fitting model to feature number 1291, X97e1d178e0060f5631143e0b6997db73
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.476596 INFO::Fitting model to feature number 1292, a32f6818da54b1c56a3f44425cbbccfd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.513953 INFO::Fitting model to feature number 1293, dcac11c3e60c3e677bb956a2e668718c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.545077 INFO::Fitting model to feature number 1294, f02787626d812f40886dff5a3b0e0fad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.576601 INFO::Fitting model to feature number 1295, X2554ec50a1739d46c1932d4b6d83bc24
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.607795 INFO::Fitting model to feature number 1296, X538b7a57f1012017b54549f162276e5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.63792 INFO::Fitting model to feature number 1297, X7fe642a71440319d58bcb8d4b05ebe59
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.672534 INFO::Fitting model to feature number 1298, d6f3ce9ec26d83ee4d0c9615049f8056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.702625 INFO::Fitting model to feature number 1299, f4cafac796fc9f62ae36a2571742ac56
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.73458 INFO::Fitting model to feature number 1300, X370b434407b901f334f8f70fcde243b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.766218 INFO::Fitting model to feature number 1301, d74a8b7d8302fd484fce300c7f7cb90b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.804654 INFO::Fitting model to feature number 1302, X6138c8580647a54710e38072430b97b6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.83594 INFO::Fitting model to feature number 1303, X71dbe01e632b84071394a25599cc1b41
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.871194 INFO::Fitting model to feature number 1304, X86b78e8356db22087c189f4fb82f3cf5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.903742 INFO::Fitting model to feature number 1305, b6bf63a9790b48aa1c1dcf66baae1882
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:52.938171 INFO::Fitting model to feature number 1306, ae999de82bf29bf047a871241ccfbe50
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.008768 INFO::Fitting model to feature number 1307, c460b8c9a6aeac17873f4e562d3c631e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.043075 INFO::Fitting model to feature number 1308, d60fd4f4bf69a870d2f8db25d29f6486
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.074413 INFO::Fitting model to feature number 1309, dbadde65054417dedcb3c441117232e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.112449 INFO::Fitting model to feature number 1310, f44cf3ad8a5561421c5c86750b907568
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.143065 INFO::Fitting model to feature number 1311, X61ecc4bd7227a33c05851dfaef146197
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.176397 INFO::Fitting model to feature number 1312, f4c831a96e731d1e318dc1f4043c2f0b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.206879 INFO::Fitting model to feature number 1313, X067720f804b765e567ec36a0554015c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.237392 INFO::Fitting model to feature number 1314, X2559f65fb87e56ab705a049961c3b596
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.26969 INFO::Fitting model to feature number 1315, e770005755e733e70992d93622110bac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.300286 INFO::Fitting model to feature number 1316, faa8d1a49c539195a242b95a3e78959a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.331955 INFO::Fitting model to feature number 1317, X0aba666ff3f58e07402855543abb0d77
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.369293 INFO::Fitting model to feature number 1318, X2ab2ed82aa11a2e91348a6323587cad0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.402339 INFO::Fitting model to feature number 1319, X7fce3acaf39b24f7c3eceaf8f4a509f1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.440757 INFO::Fitting model to feature number 1320, fc12d8847c55802fd660c17f96a4d239
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.47555 INFO::Fitting model to feature number 1321, X1237e3d469a445359f6344c4d7134df1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.510245 INFO::Fitting model to feature number 1322, be6eee56cb6557d3fde8a15a30c5f522
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.543648 INFO::Fitting model to feature number 1323, d978f3ca8b4c0a612f4a89fac040f433
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.5761 INFO::Fitting model to feature number 1324, X2d190136e77dd7bf9e1b271f6417bbbf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.607948 INFO::Fitting model to feature number 1325, X2eb915427713a21123065d4ac7d77481
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.639723 INFO::Fitting model to feature number 1326, df1d25de03331abaee085336549fccbf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.672476 INFO::Fitting model to feature number 1327, X2b33504d2a090ab8b7f93846f188fc84
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.704581 INFO::Fitting model to feature number 1328, X4514fcefdb82718d97b412184a893dc2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.742829 INFO::Fitting model to feature number 1329, X8605ed355fda390c4dbe1d350506b9bd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.775406 INFO::Fitting model to feature number 1330, e716a797ad8788f8c02701addcb9d49a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.818027 INFO::Fitting model to feature number 1331, f1c055f26db635e0551344b8e1c405ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.849869 INFO::Fitting model to feature number 1332, bf02fe466b9048dea2d3515ac8765402
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.883209 INFO::Fitting model to feature number 1333, d02e995ace556a85c3f1b34687268a29
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.91465 INFO::Fitting model to feature number 1334, X4289d6cbcecaffdc69d6ec33e18c1fcd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.945204 INFO::Fitting model to feature number 1335, X6dfdaa590e7618e0eaa870da258d29c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:53.976011 INFO::Fitting model to feature number 1336, X71faecc40533f09ce827e9ccd43f4d39
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.006608 INFO::Fitting model to feature number 1337, c5d7cdc0cd5fc338196e70876edc6839
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.038465 INFO::Fitting model to feature number 1338, X779d376dde67553423bf61868a8153ec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.069306 INFO::Fitting model to feature number 1339, c05488ee8eec8ca5ec58e73157e44cee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.101951 INFO::Fitting model to feature number 1340, X1477502529efa4a2c43b1146a94e4a1f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.13446 INFO::Fitting model to feature number 1341, X50ec7bb8c7892c4a2d8d687da6beefef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.167618 INFO::Fitting model to feature number 1342, X8a9cb38167b3fcd87381c8985bc6819a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.206495 INFO::Fitting model to feature number 1343, X9856ac7a656cb6596ba6928448ded890
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.239368 INFO::Fitting model to feature number 1344, fbc2f9405994e3c4ca9c528d362faba1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.278057 INFO::Fitting model to feature number 1345, X2d4fb4e16e2ec5372a78961fcc33e48c
## 2026-03-04 18:36:54.3111 INFO::Fitting model to feature number 1346, X5d2b612c3a59a5919ab0c68e518a0921
## 2026-03-04 18:36:54.345599 INFO::Fitting model to feature number 1347, cac45687d42dfeaa1ad41cada5a86692
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.376509 INFO::Fitting model to feature number 1348, X44d1c2442ea542df4a75d9eda9ffc130
## 2026-03-04 18:36:54.409489 INFO::Fitting model to feature number 1349, X9a9a3e95155b49f25bd224150719c489
## 2026-03-04 18:36:54.442076 INFO::Fitting model to feature number 1350, f7f661315c061fabf9c0d780d805a534
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.473378 INFO::Fitting model to feature number 1351, d2dfb8391439482e10d58726a5e4fd8a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.510584 INFO::Fitting model to feature number 1352, X20e8adbd5fcb3bd31b9a6367c226e7db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.541224 INFO::Fitting model to feature number 1353, X9ecff22f4870ac57e0897f807ee34989
## 2026-03-04 18:36:54.571317 INFO::Fitting model to feature number 1354, d87e16aa9a2cd32ddf77082b4b182168
## 2026-03-04 18:36:54.603587 INFO::Fitting model to feature number 1355, X7aa1908721c25de6fbb0de8dde870b06
## 2026-03-04 18:36:54.652447 INFO::Fitting model to feature number 1356, X5ab5afd02358060e35e6b47024d2487f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.683639 INFO::Fitting model to feature number 1357, X8b955d9fba9b023835776d9b85cf1152
## 2026-03-04 18:36:54.715107 INFO::Fitting model to feature number 1358, ab718d064dc90341d0d9feab5795cede
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.74755 INFO::Fitting model to feature number 1359, X16cab7aa6e20994e74feb61b82b156a8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.778226 INFO::Fitting model to feature number 1360, d3e5fcec01188b4f3131fd3783b39618
## 2026-03-04 18:36:54.809574 INFO::Fitting model to feature number 1361, X8299fbf1188f45f3816d1b83d91a293c
## 2026-03-04 18:36:54.840777 INFO::Fitting model to feature number 1362, d52967973c245a507023695bc1681673
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.872731 INFO::Fitting model to feature number 1363, ff1991d802d8c09b2afedb3900b6677e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:54.903853 INFO::Fitting model to feature number 1364, X8831ebf2293100dea886fc93618e58ca
## 2026-03-04 18:36:54.937057 INFO::Fitting model to feature number 1365, X4bea2aea04a20f8b9ff601799f76f4a7
## 2026-03-04 18:36:54.968937 INFO::Fitting model to feature number 1366, b2a039758e7ec54d4dc52cb1385b3086
## 2026-03-04 18:36:55.000116 INFO::Fitting model to feature number 1367, X5b025d73c766db00d0bfc9384a36c699
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.033506 INFO::Fitting model to feature number 1368, d78f0ae1704d75a9a51b5d3203b46a18
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.065066 INFO::Fitting model to feature number 1369, X12711049f5f0b0fe67416dba7d01910b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.095591 INFO::Fitting model to feature number 1370, X2d04f5a36847a364b32939c0c7cc27de
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.127509 INFO::Fitting model to feature number 1371, X5da0eb7fa58b1b776bfe7a7c5f62bcdd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.158201 INFO::Fitting model to feature number 1372, X04ce44e19521cecbaa8d4f042e0f0154
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.195035 INFO::Fitting model to feature number 1373, X806807af7f9c8e247c7a5246ab1073ba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.240248 INFO::Fitting model to feature number 1374, f2df8d8cc211effb39df7e76235b0acb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.270928 INFO::Fitting model to feature number 1375, ddb238e065dcd145239865183f17b8a1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.302029 INFO::Fitting model to feature number 1376, X23b44691db0da7bcc193c4abd5e2994f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.334614 INFO::Fitting model to feature number 1377, X78bd9cd77859e67bf1cee19881c72185
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.365975 INFO::Fitting model to feature number 1378, fa53e063ccb39a71e3528a3246ab05e0
## 2026-03-04 18:36:55.396607 INFO::Fitting model to feature number 1379, X4e13eb363f4ee6b7744b0f7af637898e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.428362 INFO::Fitting model to feature number 1380, b0022a49e8ef1e2c089251ff260e88a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.472992 INFO::Fitting model to feature number 1381, X0761fbdedb674ca859a2e6464cd54682
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.505084 INFO::Fitting model to feature number 1382, f6d259ec62674092083cb3537a21ad15
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.54091 INFO::Fitting model to feature number 1383, dd233ad31d8d9dac73d102d3e32b503e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.572356 INFO::Fitting model to feature number 1384, X375132e5c4ba0f3c06a8b7a28c0a11d5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.605207 INFO::Fitting model to feature number 1385, X481f2e633138bf123a336f9c6809e5a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.636211 INFO::Fitting model to feature number 1386, ae38dc32794bf8dbb622878b34e90fd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.667416 INFO::Fitting model to feature number 1387, b845c8e9fe940562b0a61fd239b9f6a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.705769 INFO::Fitting model to feature number 1388, X76d335ac41d3f074d0c27f5cfced7bd5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.736682 INFO::Fitting model to feature number 1389, X4d2f31f9e9a66db5837edbc35d6212ac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.774738 INFO::Fitting model to feature number 1390, X6b644754f9e17055aa0fb066a1f39660
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.805347 INFO::Fitting model to feature number 1391, eab94e28ac20551638369c32c1890678
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.843172 INFO::Fitting model to feature number 1392, e4c6d550da4a82c067f9fe2d23ad7d87
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.882592 INFO::Fitting model to feature number 1393, X6d643d243d7154690d17a90f34214a15
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.913949 INFO::Fitting model to feature number 1394, X858ef2ff37b4306d93fdc9491470e18b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.945813 INFO::Fitting model to feature number 1395, a4bf35ce0ea6e33b1ca5fce93f2c72b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:55.983817 INFO::Fitting model to feature number 1396, a61c2f470b263001d328b36c16ee6a8f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.019953 INFO::Fitting model to feature number 1397, X947d4cca69b99f31480944eb79193ce7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.052107 INFO::Fitting model to feature number 1398, X271cb325b611cad7db3939e8e11c5541
## 2026-03-04 18:36:56.083937 INFO::Fitting model to feature number 1399, X823d7476979a0b8fc9479ae9113b7bf7
## 2026-03-04 18:36:56.119757 INFO::Fitting model to feature number 1400, X5ddfcafe5718c9c9f989fe6d946e741c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.150431 INFO::Fitting model to feature number 1401, X6cdde78777913ddf95c98995c3b4f62a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.18414 INFO::Fitting model to feature number 1402, X1bd40086437793e1ba33f5cbd83d6fe5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.235027 INFO::Fitting model to feature number 1403, X9abe6885f490b0b13978b8e993cf1e95
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.268938 INFO::Fitting model to feature number 1404, X316ab514dc5c338e49f89b0732fc5fca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.305056 INFO::Fitting model to feature number 1405, X648bee2fa5f6c2bee7662df9d276a906
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.35046 INFO::Fitting model to feature number 1406, X4b2e82cadba9d973d121de999f3199f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.384761 INFO::Fitting model to feature number 1407, X21cd78cabd2ea1730146030a8ee5bb3a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.420566 INFO::Fitting model to feature number 1408, X5a2d5c1a794e467db1cf2588418be645
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.459382 INFO::Fitting model to feature number 1409, X169969cb9914e40725adc23728995292
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.491417 INFO::Fitting model to feature number 1410, c6a0cea6e8edca58bc207bfd93291981
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.523266 INFO::Fitting model to feature number 1411, fc1f59f65cb9dbb055d1b5e6c36ee5fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.560819 INFO::Fitting model to feature number 1412, X34c29a3124e99390fe16cb7e3cf1542b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.592158 INFO::Fitting model to feature number 1413, X43690cf299862046d6c4f3c14cc26fbf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.624281 INFO::Fitting model to feature number 1414, X6e151e6f762a3ddd31df62acda1fa509
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.655886 INFO::Fitting model to feature number 1415, X776f5ce8233c756ba7ae51611d4021c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.68743 INFO::Fitting model to feature number 1416, X8123f26005a06511ef36459cd14c4be7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.719296 INFO::Fitting model to feature number 1417, X9498050fda592707a74ae42cb43cf376
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.749828 INFO::Fitting model to feature number 1418, b7824ad5efd327ec9aa70463ffe41313
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.782816 INFO::Fitting model to feature number 1419, be122441d4bfbf1b799a85d8765fb831
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:56.821152 INFO::Fitting model to feature number 1420, f522ba24c2cefcb9eaecc3bd442fc926
## 2026-03-04 18:36:56.852309 INFO::Fitting model to feature number 1421, X7d5650596629e727155c043cefe80839
## 2026-03-04 18:36:56.890795 INFO::Fitting model to feature number 1422, a57c1c3550d34aabb91ec1aaeccc3da0
## 2026-03-04 18:36:56.921935 INFO::Fitting model to feature number 1423, a6c825f96532401aa96ebbb3ca8cc4d5
## 2026-03-04 18:36:56.962491 INFO::Fitting model to feature number 1424, e6100b80a38222d7f4b53640a2813c92
## 2026-03-04 18:36:56.999968 INFO::Fitting model to feature number 1425, X5e298cb51b192a2586c265bd13ecddec
## 2026-03-04 18:36:57.036887 INFO::Fitting model to feature number 1426, d6bbba5be34ced0b2ddcf6d419a290fd
## 2026-03-04 18:36:57.072566 INFO::Fitting model to feature number 1427, X1e5afa7c69a7d5bc728c0b1c178825a4
## 2026-03-04 18:36:57.103998 INFO::Fitting model to feature number 1428, a72f4146638b9975a2f70a429c7888aa
## 2026-03-04 18:36:57.133743 INFO::Fitting model to feature number 1429, a6dd7f11774533a463aa01d6c837dc00
## 2026-03-04 18:36:57.164278 INFO::Fitting model to feature number 1430, X394106ad4bcc01e28e654f8b904a43f0
## 2026-03-04 18:36:57.205636 INFO::Fitting model to feature number 1431, X3323196417bad3bec4e7a302d5dfdd61
## 2026-03-04 18:36:57.242793 INFO::Fitting model to feature number 1432, X568794df14ca80a20df326e942094f60
## 2026-03-04 18:36:57.280153 INFO::Fitting model to feature number 1433, dc70f2dde3a4e16635013a39d5bc4738
## 2026-03-04 18:36:57.311031 INFO::Fitting model to feature number 1434, X598d0ff1000e5b59358afb0c03dd8f91
## 2026-03-04 18:36:57.347343 INFO::Fitting model to feature number 1435, X48d757c9fc57146957b1f61c94e1cc19
## 2026-03-04 18:36:57.378059 INFO::Fitting model to feature number 1436, c7c955897cbfe6837c99e7cd44d0887e
## 2026-03-04 18:36:57.413579 INFO::Fitting model to feature number 1437, e369521c3e10197580daa61efd80c436
## 2026-03-04 18:36:57.444376 INFO::Fitting model to feature number 1438, f7e47a7f64694d3180b261250205ce4a
## 2026-03-04 18:36:57.481023 INFO::Fitting model to feature number 1439, X3481a13605bd74eb2f7cc589b20486b5
## 2026-03-04 18:36:57.517636 INFO::Fitting model to feature number 1440, X7e75bd65b31ecb2641e7317dac525229
## 2026-03-04 18:36:57.550116 INFO::Fitting model to feature number 1441, b70927f642e4add511f7bd662c667a58
## 2026-03-04 18:36:57.583263 INFO::Fitting model to feature number 1442, X224b3eccb629d9e4a636f3af73fb0268
## 2026-03-04 18:36:57.620645 INFO::Fitting model to feature number 1443, eac0c8b1f85bfc7dba195c75ea8f3561
## 2026-03-04 18:36:57.658314 INFO::Fitting model to feature number 1444, X56c30e736c54be6e83da0343c057e762
## 2026-03-04 18:36:57.689777 INFO::Fitting model to feature number 1445, X9d1c840bc95d7a1d599e26ec32d6be01
## 2026-03-04 18:36:57.721829 INFO::Fitting model to feature number 1446, X49f668308b55c43b37f983fc8bc7ce4b
## 2026-03-04 18:36:57.752371 INFO::Fitting model to feature number 1447, X8fa3bb6a74294c2f87ffa923221aaef1
## 2026-03-04 18:36:57.789558 INFO::Fitting model to feature number 1448, X7c54d19b47ce29e4e80f24ec6e8abcfc
## 2026-03-04 18:36:57.820656 INFO::Fitting model to feature number 1449, X4202fe0cef4d4c2bf106ba2027a23bb8
## 2026-03-04 18:36:57.852991 INFO::Fitting model to feature number 1450, X22bf5bfb7975b9c37934c7334e313408
## 2026-03-04 18:36:57.883406 INFO::Fitting model to feature number 1451, X388c2c1d2035cce97331699e6dcf29b4
## 2026-03-04 18:36:57.914468 INFO::Fitting model to feature number 1452, X439a7765de6fa0674816e837c3810090
## 2026-03-04 18:36:57.945001 INFO::Fitting model to feature number 1453, bd536adea9f0edde6e5eafb1da9be283
## 2026-03-04 18:36:57.977393 INFO::Fitting model to feature number 1454, e3690c2df119bc765abb284d23301e5c
## 2026-03-04 18:36:58.007957 INFO::Fitting model to feature number 1455, X86e184a98df6da6c11e688e0e794d017
## 2026-03-04 18:36:58.054369 INFO::Fitting model to feature number 1456, X20c0706bcdf2adf7708e01585bf70288
## 2026-03-04 18:36:58.088634 INFO::Fitting model to feature number 1457, X2414414f9af0314343395f48bca12ec7
## 2026-03-04 18:36:58.124431 INFO::Fitting model to feature number 1458, X0616a2d8a3a5c92df8a0acd7c90be22a
## 2026-03-04 18:36:58.156684 INFO::Fitting model to feature number 1459, X03d7aab491a1a55630aa6c485fafe72c
## 2026-03-04 18:36:58.188117 INFO::Fitting model to feature number 1460, b0031e259b262ad8ac0ff64393f013fa
## 2026-03-04 18:36:58.219931 INFO::Fitting model to feature number 1461, X46ef44d3a2f2d113352a55c8d8ffad04
## 2026-03-04 18:36:58.25163 INFO::Fitting model to feature number 1462, b713d99eb533b299629c1394f75b6926
## 2026-03-04 18:36:58.283392 INFO::Fitting model to feature number 1463, X14dbca841e6bb9568be475a3360d5a58
## 2026-03-04 18:36:58.320702 INFO::Fitting model to feature number 1464, X4e4c5f7c2c776156be82a7e1d24cdffb
## 2026-03-04 18:36:58.351947 INFO::Fitting model to feature number 1465, X34c3988ed3a5e3630b8488505b7f9f98
## 2026-03-04 18:36:58.388882 INFO::Fitting model to feature number 1466, X72ee97f042fd5f6508c8f55419c8862e
## 2026-03-04 18:36:58.421932 INFO::Fitting model to feature number 1467, X7b3fe7dea742fd319815e05d8c5f9fe1
## 2026-03-04 18:36:58.458585 INFO::Fitting model to feature number 1468, X5006a0d2b1a915d327889a08c9e8885f
## 2026-03-04 18:36:58.48984 INFO::Fitting model to feature number 1469, a4bd3b3c5909749d8f08456f744defa0
## 2026-03-04 18:36:58.523237 INFO::Fitting model to feature number 1470, bf2eac04920ee7f7f1d299a292613df2
## 2026-03-04 18:36:58.559727 INFO::Fitting model to feature number 1471, e79c3000d572833b5371f40ba4a529ea
## 2026-03-04 18:36:58.58996 INFO::Fitting model to feature number 1472, X539fab4811fbd687bbeca6e06a832e8e
## 2026-03-04 18:36:58.628171 INFO::Fitting model to feature number 1473, ca1c86768a0a1d127988348bb1a84a4e
## 2026-03-04 18:36:58.659235 INFO::Fitting model to feature number 1474, X5f9c85a0e24285e0d751487138070fcb
## 2026-03-04 18:36:58.692836 INFO::Fitting model to feature number 1475, X0432cc28a194dcd588a5a751fac80832
## 2026-03-04 18:36:58.724278 INFO::Fitting model to feature number 1476, a75f41af38c6efe6ff870067af0eecd5
## 2026-03-04 18:36:58.756995 INFO::Fitting model to feature number 1477, X5712fb6a5e9c90125f5029a8cc18f538
## 2026-03-04 18:36:58.789285 INFO::Fitting model to feature number 1478, e35b85bcf133e5c1406ade391f27f35e
## 2026-03-04 18:36:58.825478 INFO::Fitting model to feature number 1479, X377d5ee705688586a17936fb3804994a
## 2026-03-04 18:36:58.857429 INFO::Fitting model to feature number 1480, X5b73b77b0b395c24b331446c95c75533
## 2026-03-04 18:36:58.900882 INFO::Fitting model to feature number 1481, de4b3ada2cc38daf9fcdc0741de5dca7
## 2026-03-04 18:36:58.934179 INFO::Fitting model to feature number 1482, X9d9c864dccdebdc7a243f28d6510f2ef
## 2026-03-04 18:36:58.966812 INFO::Fitting model to feature number 1483, X47cf08ac4ef1b9bbfb862292206a9f7c
## 2026-03-04 18:36:58.997232 INFO::Fitting model to feature number 1484, X531af51c5105e6170ef66e92168256c1
## 2026-03-04 18:36:59.029909 INFO::Fitting model to feature number 1485, acd2936ee50c99be98115e15414cb216
## 2026-03-04 18:36:59.065425 INFO::Fitting model to feature number 1486, c0484a4bafa44a0dd28c9d0880a3d591
## 2026-03-04 18:36:59.096387 INFO::Fitting model to feature number 1487, X211046170b10ef39b5c78eae5a7249d3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:36:59.127647 INFO::Fitting model to feature number 1488, X75c94a2e7f736b8dd4894c788c15ead7
## 2026-03-04 18:36:59.158159 INFO::Fitting model to feature number 1489, b2c82de8806b29b0f4f22aff7ee128c6
## 2026-03-04 18:36:59.188874 INFO::Fitting model to feature number 1490, X44cbcbce5883cdcee679699310064855
## 2026-03-04 18:36:59.226516 INFO::Fitting model to feature number 1491, X28d4648bb079932a475e8b4e2cb7b203
## 2026-03-04 18:36:59.258265 INFO::Fitting model to feature number 1492, X73ed71d087d17b73df11c5ef4d509a48
## 2026-03-04 18:36:59.289383 INFO::Fitting model to feature number 1493, aee63dacdb65a4e383d739b9a823c6d5
## 2026-03-04 18:36:59.322412 INFO::Fitting model to feature number 1494, X02d8636ae9a2d97764ca155d429a3573
## 2026-03-04 18:36:59.352321 INFO::Fitting model to feature number 1495, X7c2c90b4dd80a6c5e232c200e12bb135
## 2026-03-04 18:36:59.385166 INFO::Fitting model to feature number 1496, X0de79de542d618ec0a3640507d175b58
## 2026-03-04 18:36:59.419614 INFO::Fitting model to feature number 1497, b6d0a5d1db20bb06207235f7ca1c80e6
## 2026-03-04 18:36:59.450311 INFO::Fitting model to feature number 1498, X61bedb706ed0ebf66136c2de5e894ae0
## 2026-03-04 18:36:59.480599 INFO::Fitting model to feature number 1499, X85e51deb7d526a41feb684d724a249f7
## 2026-03-04 18:36:59.511076 INFO::Fitting model to feature number 1500, X50e3f160c939685dd9fe2516743de428
## 2026-03-04 18:36:59.542288 INFO::Fitting model to feature number 1501, X9be98a384b63040c28b7ba80f5bcfaa6
## 2026-03-04 18:36:59.57233 INFO::Fitting model to feature number 1502, bb4981f12957c9153404475002597fcf
## 2026-03-04 18:36:59.602009 INFO::Fitting model to feature number 1503, X55b6147730ff0b48500535882c1c90c8
## 2026-03-04 18:36:59.63449 INFO::Fitting model to feature number 1504, X70cec6470fc234fda2298dd42df999c8
## 2026-03-04 18:36:59.675991 INFO::Fitting model to feature number 1505, ab7683701ff9cd914c709a432024c95a
## 2026-03-04 18:36:59.71525 INFO::Fitting model to feature number 1506, ebae2a99fbe85d6f7192a78fccdb8fdb
## 2026-03-04 18:36:59.746361 INFO::Fitting model to feature number 1507, X82949f67dacd7431521ab84ca247a5d3
## 2026-03-04 18:36:59.777454 INFO::Fitting model to feature number 1508, a2a7e83f01568ca2ffc7aaeb75be4515
## 2026-03-04 18:36:59.808957 INFO::Fitting model to feature number 1509, X7cb5a8b160490076fb6549561c3423f1
## 2026-03-04 18:36:59.840361 INFO::Fitting model to feature number 1510, X003efbd0f8008eb1a10db7323ae5b9cd
## 2026-03-04 18:36:59.872878 INFO::Fitting model to feature number 1511, X658685f5b6770693b0e2a9758ffee903
## 2026-03-04 18:36:59.910198 INFO::Fitting model to feature number 1512, X124e851f1317df1c44f69576ded77eea
## 2026-03-04 18:36:59.948448 INFO::Fitting model to feature number 1513, X247c6962f97e1327653c49a55baaf0ef
## 2026-03-04 18:36:59.984809 INFO::Fitting model to feature number 1514, X95acc8353370e4c6473aabd257f89145
## 2026-03-04 18:37:00.017206 INFO::Fitting model to feature number 1515, a6be673d83361f0db6b04f4f6368aac8
## 2026-03-04 18:37:00.051482 INFO::Fitting model to feature number 1516, e78522b1e525cf73de0a832161386cfd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:00.083758 INFO::Fitting model to feature number 1517, X40ead826bf0cabf57669891968fb8028
## 2026-03-04 18:37:00.11974 INFO::Fitting model to feature number 1518, X4eafbc59f43c8a0357c8c6c01196f799
## 2026-03-04 18:37:00.155608 INFO::Fitting model to feature number 1519, b69359c339f22ecb5db0d3f40004badc
## 2026-03-04 18:37:00.189536 INFO::Fitting model to feature number 1520, X462a59f090f13da057d0e2bcf7113ffc
## 2026-03-04 18:37:00.223971 INFO::Fitting model to feature number 1521, df8cea54149f075a0cf136ebd86633f1
## 2026-03-04 18:37:00.256468 INFO::Fitting model to feature number 1522, fb647843e3c65da0b4fcddbc5afa2900
## 2026-03-04 18:37:00.288927 INFO::Fitting model to feature number 1523, ff1b4e1105cab28d1b9274d9aad00686
## 2026-03-04 18:37:00.323469 INFO::Fitting model to feature number 1524, X350fce626781822ad4f661f8ab50950e
## 2026-03-04 18:37:00.354873 INFO::Fitting model to feature number 1525, X516edf87db3e79ac6cb028e4866128c4
## 2026-03-04 18:37:00.387944 INFO::Fitting model to feature number 1526, c6fba20dabb5158b66a8c275f1972c9c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:00.420624 INFO::Fitting model to feature number 1527, cea3551436db23c969c3846aaca74006
## 2026-03-04 18:37:00.453165 INFO::Fitting model to feature number 1528, X40ee6ad12f5a84a745007123217b2f14
## 2026-03-04 18:37:00.485608 INFO::Fitting model to feature number 1529, X41bfe90c919674be43e40cc1ed15f397
## 2026-03-04 18:37:00.531386 INFO::Fitting model to feature number 1530, X5888c34571979386d88d40e153020962
## 2026-03-04 18:37:00.56356 INFO::Fitting model to feature number 1531, X9d1c414e4e43f9035786f41fe51a2c37
## 2026-03-04 18:37:00.601152 INFO::Fitting model to feature number 1532, X2ed61e783e248ef47514084aba18efe5
## 2026-03-04 18:37:00.632278 INFO::Fitting model to feature number 1533, X85eb6d9d8178170857155c18c4747e83
## 2026-03-04 18:37:00.663279 INFO::Fitting model to feature number 1534, X192171db9d9d7eef504462e639aff14e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:00.695023 INFO::Fitting model to feature number 1535, X382cb4846b04b8a16502cbccaa8fd78b
## 2026-03-04 18:37:00.725824 INFO::Fitting model to feature number 1536, X3bc397e5ed09b0d3152b2715837be559
## 2026-03-04 18:37:00.757129 INFO::Fitting model to feature number 1537, X45d5c0bb1e1bb6300d395413db2070c8
## 2026-03-04 18:37:00.788408 INFO::Fitting model to feature number 1538, X712381dc322c52191a6601dfa1ea3d5b
## 2026-03-04 18:37:00.818657 INFO::Fitting model to feature number 1539, X8a5fcf814785cd96fc7be2582c438c7c
## 2026-03-04 18:37:00.856211 INFO::Fitting model to feature number 1540, ae2ef6affcdf43d6f0f516c8212da640
## 2026-03-04 18:37:00.890118 INFO::Fitting model to feature number 1541, X17122257dd2d67ac0e4167bd03651e51
## 2026-03-04 18:37:00.921023 INFO::Fitting model to feature number 1542, X30d593c68137f5c78b3750e68af4b9f7
## 2026-03-04 18:37:00.952038 INFO::Fitting model to feature number 1543, X4540e89c3c1314208895f374fe7c9a6a
## 2026-03-04 18:37:00.985555 INFO::Fitting model to feature number 1544, X6277a075fc7865b1e27d287c43a7e57c
## 2026-03-04 18:37:01.01718 INFO::Fitting model to feature number 1545, X7a7377a0454382bbad229ff7c5d80546
## 2026-03-04 18:37:01.048162 INFO::Fitting model to feature number 1546, X7c3929432548e1d0644f5b531b396ad0
## 2026-03-04 18:37:01.07983 INFO::Fitting model to feature number 1547, b40189cb15dd5644b0adea1ff4afb315
## 2026-03-04 18:37:01.117671 INFO::Fitting model to feature number 1548, X95d765a3b3dfd08fdb458a3c6a7054b9
## 2026-03-04 18:37:01.152525 INFO::Fitting model to feature number 1549, e8019eaaf310c14bc4b7be16e7d3bfa1
## 2026-03-04 18:37:01.185345 INFO::Fitting model to feature number 1550, X7fbb704ff0cbacc0253188c453d69e36
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.216563 INFO::Fitting model to feature number 1551, X9be5e759a3fdd983eedea94bfca182f9
## 2026-03-04 18:37:01.247971 INFO::Fitting model to feature number 1552, X13bc51f05823bc48301a5d1c4bde9f5a
## 2026-03-04 18:37:01.279055 INFO::Fitting model to feature number 1553, X3e82bcc53ad83cc493c60475fbcc727d
## 2026-03-04 18:37:01.309469 INFO::Fitting model to feature number 1554, X632aa2401824dfbff3dccab0c136de74
## 2026-03-04 18:37:01.34974 INFO::Fitting model to feature number 1555, X6fe547e6fe340a40c701e84a3d8eb1cb
## 2026-03-04 18:37:01.388998 INFO::Fitting model to feature number 1556, X8387bba230a7f9e754d83c9171c28b31
## 2026-03-04 18:37:01.427078 INFO::Fitting model to feature number 1557, a0a9b9c3660fa6ae48f562a494208215
## 2026-03-04 18:37:01.458891 INFO::Fitting model to feature number 1558, X56afd101747172bfbeae0ac912dfc469
## 2026-03-04 18:37:01.489828 INFO::Fitting model to feature number 1559, X84bfbbf1dcab3c3f73b0089d3bb32c7c
## 2026-03-04 18:37:01.520495 INFO::Fitting model to feature number 1560, X42d802e1f056d33a6fc7e8db5a5bbff7
## 2026-03-04 18:37:01.551115 INFO::Fitting model to feature number 1561, X6e330d70c19eb8bf817a65230895f026
## 2026-03-04 18:37:01.582518 INFO::Fitting model to feature number 1562, a1b0ba01f038cc733b4bfe022a6525dc
## 2026-03-04 18:37:01.613956 INFO::Fitting model to feature number 1563, bc8fd2265d620b8bb9eea8d62bf623c6
## 2026-03-04 18:37:01.644252 INFO::Fitting model to feature number 1564, X5091bd9b1943a71c946ae11b4e07b7e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.675541 INFO::Fitting model to feature number 1565, X5d3a339376909059f6ad187fdc28d6df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.706355 INFO::Fitting model to feature number 1566, X5deed8fedf1ec0ed7efc98743a7dbb25
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.736439 INFO::Fitting model to feature number 1567, d77a370884869d02025c707a1059ebf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.76672 INFO::Fitting model to feature number 1568, X96535061b17238b675b7c8397d18b1db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.799143 INFO::Fitting model to feature number 1569, c94467ba9b103e3de9d85eaca888a027
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.830125 INFO::Fitting model to feature number 1570, X0b7c20e95c7d503306b44bc853c0f948
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.86107 INFO::Fitting model to feature number 1571, X589a69391f6c119b48b9884552436536
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.892342 INFO::Fitting model to feature number 1572, X922ebffab737973fc23feaefe8df04e8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.924307 INFO::Fitting model to feature number 1573, X183006dbbfa96d0df78d911109617c04
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.954541 INFO::Fitting model to feature number 1574, X47add6d359e281d9307ef7b818b1132c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:01.984458 INFO::Fitting model to feature number 1575, accdc77f3b67a008ff30692229aa28ef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.014737 INFO::Fitting model to feature number 1576, X76215b5328d849f5bb5aa09fee01861f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.045478 INFO::Fitting model to feature number 1577, fe585ab842d898d8e35d20f9005df2d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.07615 INFO::Fitting model to feature number 1578, X52c94f53dc5b8610cb0a96ec22b853eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.107619 INFO::Fitting model to feature number 1579, X7c5b049887df626ca87a19fbbe1609d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.150052 INFO::Fitting model to feature number 1580, e8bb3d4f052a5add7027c46cbaead30b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.18069 INFO::Fitting model to feature number 1581, X5708f680f66bc7b8047b515f2b4c4bc3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.211764 INFO::Fitting model to feature number 1582, X1948dc0f10daf853102b0731a6186af4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.241988 INFO::Fitting model to feature number 1583, X43286dda280419b7f392dd6ce424759a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.274642 INFO::Fitting model to feature number 1584, ef67f9a0ad8f89623991c3f3d2ec7d16
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.306039 INFO::Fitting model to feature number 1585, bd17fc6e6fee81cc3e9ec0905eb366b3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.337191 INFO::Fitting model to feature number 1586, X5c0401ad6b1ee6bb92b542bd10210f51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.368537 INFO::Fitting model to feature number 1587, X4807347adae205610568f73ccc96c6a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.400319 INFO::Fitting model to feature number 1588, f4bc6194719dec9e61ae79c9b52fcbf3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.432692 INFO::Fitting model to feature number 1589, X1b76fcb2cecd29cdfff9a0af888d4169
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.465342 INFO::Fitting model to feature number 1590, X6158e79837bbf56a578b6d8b44ef7adb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.499618 INFO::Fitting model to feature number 1591, d71613d583910ae8267587f2b140ceb7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.529628 INFO::Fitting model to feature number 1592, X8ceade32cdc05781325d8ad50f31f2b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.560537 INFO::Fitting model to feature number 1593, X58c9c50389bbfe8ad20f1d089a45172d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.590879 INFO::Fitting model to feature number 1594, X9f4c8ff7db4a28f9d8a9d36c2c773eba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.621474 INFO::Fitting model to feature number 1595, eb740349dd1e6ad53ae89a46ec50ecb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.652404 INFO::Fitting model to feature number 1596, X915e1c30b64ef9a338d6c10600e4a617
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.683796 INFO::Fitting model to feature number 1597, f027bfc4feea257cddbf035993ba545a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.714594 INFO::Fitting model to feature number 1598, X6079966a4d2ff4c7cdc74e1741c30fd6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.745201 INFO::Fitting model to feature number 1599, c6ce5406ec57d4e6377b8c641e8975da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.775846 INFO::Fitting model to feature number 1600, X56c5746a26ab97c041720b362ad8661e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.808527 INFO::Fitting model to feature number 1601, X6baf37ed88a786881d547a08161907e4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.840204 INFO::Fitting model to feature number 1602, a34ec679cd800d55496935d18635358d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.871221 INFO::Fitting model to feature number 1603, cc01d85478bd68c5a0abdb6b92a6bb57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.910064 INFO::Fitting model to feature number 1604, d10cef946e0c66aff6e1b5b7976792ad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.952757 INFO::Fitting model to feature number 1605, X1596a42ede5440896d482a43e80cf4be
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:02.983542 INFO::Fitting model to feature number 1606, X6174b3382fc068d14dd5ff3e1b5f7657
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.013705 INFO::Fitting model to feature number 1607, X80ab0c108d89e0bea0795b22fdd78c45
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.04504 INFO::Fitting model to feature number 1608, X84ad86e97445e7787a94091968391d2d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.075282 INFO::Fitting model to feature number 1609, bc148193d702e654dd5beafc48ed3b3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.10582 INFO::Fitting model to feature number 1610, X5819fc035f12f49f114cdc42ebd7f872
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.138739 INFO::Fitting model to feature number 1611, c1ad43ca73c671841937ae7a7f12bd00
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.168058 INFO::Fitting model to feature number 1612, X3857b4ad66f5288bdfaee996a22c2dd5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.199732 INFO::Fitting model to feature number 1613, X5ec337477fdc4a209b4d08ba19627c12
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.231081 INFO::Fitting model to feature number 1614, X88735c01a26d318d22ab194b70d82670
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.265591 INFO::Fitting model to feature number 1615, X993b533f1f2a19c778a99ee2d73d6401
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.297331 INFO::Fitting model to feature number 1616, c50ece6688d9455f2089d10c593dc5e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.32828 INFO::Fitting model to feature number 1617, X1f07789f2a30d8406f4e1e6c4e59d512
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.365445 INFO::Fitting model to feature number 1618, X12957d1b9a2b84263f7c2577a0be5ffe
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.395777 INFO::Fitting model to feature number 1619, X8ee9e03c0b6df43238a2bac06e8157c6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.426833 INFO::Fitting model to feature number 1620, X7c1a951cbb503325aa9c776b49db8e62
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.458075 INFO::Fitting model to feature number 1621, X37af96a311f2b6d6d9e6f1aeac8adcbb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.488796 INFO::Fitting model to feature number 1622, aa0d9d7fa21622e919bc657a1c27acea
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.522416 INFO::Fitting model to feature number 1623, a8d44a386b74d5efeb2ab88dd7b14d81
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.552978 INFO::Fitting model to feature number 1624, X6e2c98b75514f3cec57b18d85e847d56
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.583443 INFO::Fitting model to feature number 1625, b6e3cadd45df3d2e2e465e4b2e49fbf6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.613932 INFO::Fitting model to feature number 1626, X8e714ede8d3de600e6e1e85b93ffff3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.644395 INFO::Fitting model to feature number 1627, e36960ffa2d5b1880b964b90d5e01ad1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.678747 INFO::Fitting model to feature number 1628, X1012dfedbd7005d166522a23a70e17e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.720818 INFO::Fitting model to feature number 1629, ee8bdf48063e7dee314d7106122de356
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.752832 INFO::Fitting model to feature number 1630, X73741e0f0fc53215dfc000d6bcd03bf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.78355 INFO::Fitting model to feature number 1631, b4f013a316dd9db71487186d4e5f8c48
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.814053 INFO::Fitting model to feature number 1632, d80f4c95bba89e950a462cde932b4a19
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.84446 INFO::Fitting model to feature number 1633, ec281d5415aa63e552133483f7b7f263
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.874479 INFO::Fitting model to feature number 1634, X3889da5feeb2f015f4e43df76d1ba1c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.90581 INFO::Fitting model to feature number 1635, X0fdfc03c69957857ea0aa99afb974be6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.938565 INFO::Fitting model to feature number 1636, bb71589a895e6aac9ac57ef82e01389f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:03.970336 INFO::Fitting model to feature number 1637, X6f8c98e3a54279e65fb5b3d6c4596b06
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.001049 INFO::Fitting model to feature number 1638, ef7d9398b53a4c43e226717ea5256aa8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.032276 INFO::Fitting model to feature number 1639, X0744af08f13a84f475cacdfe0a7eb563
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.064643 INFO::Fitting model to feature number 1640, b933d7cfe0a4a8671012717ad8f3d467
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.096071 INFO::Fitting model to feature number 1641, X8e898fac58bba1bbf6b1d86cda3d468a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.128235 INFO::Fitting model to feature number 1642, f90ae11968e525ce0b184a7d096010b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.159186 INFO::Fitting model to feature number 1643, X6149b30955bdf96ec6892f0d385ddeb4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.190308 INFO::Fitting model to feature number 1644, X878ab77c15d89ce76f14512b57b66554
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.221317 INFO::Fitting model to feature number 1645, X97fb8735d375018eef0bc292ced80c3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.253128 INFO::Fitting model to feature number 1646, X122ae66707c8158bbd1a25e11191df1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.283481 INFO::Fitting model to feature number 1647, X08d50335ad948012223e0bcb0e02a98c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.314337 INFO::Fitting model to feature number 1648, X65450912952dfde7e9482ff90cb96c14
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.345302 INFO::Fitting model to feature number 1649, X9a11808482e84a5875187726140985df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.376089 INFO::Fitting model to feature number 1650, bcab38bdf56dea2067767829a85b6081
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.407293 INFO::Fitting model to feature number 1651, X70aa06947222e61a0342014baad5b4fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.438908 INFO::Fitting model to feature number 1652, X75df890f90e4d0da4c0bad5bfe8344c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.470982 INFO::Fitting model to feature number 1653, X36a5741f495aca0b66ad97f97dc7a194
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.513138 INFO::Fitting model to feature number 1654, X10415149c4af2e9597ec437422950ed1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.543922 INFO::Fitting model to feature number 1655, X57f59e0ae6b1c8ce77629f2b8cb9b6ec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.574791 INFO::Fitting model to feature number 1656, X839baf767eec4eed3854e3ac4ff072c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.605776 INFO::Fitting model to feature number 1657, X2b799a142091659df01cda80e96b81ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.636313 INFO::Fitting model to feature number 1658, X947eb89bd4aa5e20d8bf953ea6c0d081
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.666736 INFO::Fitting model to feature number 1659, cc61a14c8cb9f6513a9103c1a21d0b55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.697734 INFO::Fitting model to feature number 1660, c65985845943a827744fa6d3de55535d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.728241 INFO::Fitting model to feature number 1661, X88117a1cd850c7fafeeef4fa07f31d41
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.759226 INFO::Fitting model to feature number 1662, bedc4e430c931cfaab616b853069f814
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.789956 INFO::Fitting model to feature number 1663, X54559e1264fb8ce1cc649277c8326fad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.820679 INFO::Fitting model to feature number 1664, X10c2b175dcb8d7348ee54da806d75dde
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.85189 INFO::Fitting model to feature number 1665, X5774e92983913d64432d9cee33e6fd3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.882532 INFO::Fitting model to feature number 1666, X8f6915518fcb472c94ac31408b61a0da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.913613 INFO::Fitting model to feature number 1667, d74a6273b1480e0d17fe444efd8c855d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.944312 INFO::Fitting model to feature number 1668, ebb85b3d0073e8567c7187e4d074f323
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:04.975719 INFO::Fitting model to feature number 1669, ff0b4d81a135df7e6b677c1e0e371d26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.008117 INFO::Fitting model to feature number 1670, c7f5deefec0090691af7e3564a82429d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.039831 INFO::Fitting model to feature number 1671, d624453d3dcd59833f711d70d8be5e2a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.071778 INFO::Fitting model to feature number 1672, X6a7065ab25b547f46bbb229e5ed5a8d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.103219 INFO::Fitting model to feature number 1673, X7b7686f9f714d532db91900814a6ffe9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.134758 INFO::Fitting model to feature number 1674, X858cadaa6bed41000fb485d316ba52a4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.165716 INFO::Fitting model to feature number 1675, eb56c11fff2b5d82b6918720792b3902
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.197564 INFO::Fitting model to feature number 1676, a693465c86c54e2943f05ad39652a938
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.229231 INFO::Fitting model to feature number 1677, f197afa977ca78606df5de118bbe0bb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.26061 INFO::Fitting model to feature number 1678, X9e1f8fa2e00d442b56a9262b77fdd34c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.304094 INFO::Fitting model to feature number 1679, X466af77afd608fc85dfc61d3a08010fa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.334763 INFO::Fitting model to feature number 1680, X615e8c85484fde9ca7fc10a9ca8d2aa5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.36707 INFO::Fitting model to feature number 1681, b34dcecfc1eb9380c919f8566479a79e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.397769 INFO::Fitting model to feature number 1682, X459037d1dd08629797b49e758ec86da7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.429565 INFO::Fitting model to feature number 1683, X30d07efaf39113af71b777723b986086
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.461065 INFO::Fitting model to feature number 1684, cf22fbcb43140202c9d69acac1d2d7d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.492162 INFO::Fitting model to feature number 1685, X708f08717de48853cf9afe00b1198a27
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.524081 INFO::Fitting model to feature number 1686, a1a62a39c30efb6136daa1c51cba279a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.554728 INFO::Fitting model to feature number 1687, a48690290fc79a4f0d211eb587077fa8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.58562 INFO::Fitting model to feature number 1688, f67744afabc94b4f4c0a3b49ac8d862d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.615962 INFO::Fitting model to feature number 1689, X05269a738f261937bd7740c6501454e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.647036 INFO::Fitting model to feature number 1690, X2eb6747a7aebca6666a096a37071d7e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.676984 INFO::Fitting model to feature number 1691, X3db4777463dd9b345036bcdb868ef547
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.706813 INFO::Fitting model to feature number 1692, X53b62dd9d227d3bc7297f0dbd48a4989
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.73752 INFO::Fitting model to feature number 1693, X628547bdac555ee7526d0216152bcf26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.768355 INFO::Fitting model to feature number 1694, X6f1fbf83f07a443433b7f1e88a3dd392
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.799425 INFO::Fitting model to feature number 1695, X9bb4fc212f3b982613792c538ea52694
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.830455 INFO::Fitting model to feature number 1696, a0d57e53d8bc50f10b7857f69a8c5fa0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:05.952964 INFO::Counting total values for each feature
## 2026-03-04 18:37:06.01542 INFO::Writing filtered data to file ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/features/filtered_data.tsv
## 2026-03-04 18:37:06.087282 INFO::Writing filtered, normalized data to file ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/features/filtered_data_norm.tsv
## 2026-03-04 18:37:06.179532 INFO::Writing filtered, normalized, transformed data to file ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/features/filtered_data_norm_transformed.tsv
## 2026-03-04 18:37:06.333168 INFO::Writing residuals to file ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/fits/residuals.rds
## 2026-03-04 18:37:06.356835 INFO::Writing fitted values to file ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/fits/fitted.rds
## 2026-03-04 18:37:06.37377 INFO::Writing extracted random effects to file ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/fits/ranef.rds
## 2026-03-04 18:37:06.38537 INFO::Writing all results to file (ordered by increasing q-values): ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/all_results.tsv
## 2026-03-04 18:37:06.427678 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.050000 ) to file (ordered by increasing q-values): ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/significant_results.tsv
## 2026-03-04 18:37:06.434258 INFO::Writing heatmap of significant results to file: ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations/heatmap.pdf
## [1] "There is not enough metadata in the associations to create a heatmap plot. Please review the associations in text output file."
## 2026-03-04 18:37:06.447141 INFO::Writing association plots (one for each significant association) to output folder: ./output/maaslin2/pooled-output/rare/20260304_18_36_08-NMR-OTU-age-234-ref-agegroup0_10-with-relations
## 2026-03-04 18:37:06.476166 INFO::Plotting associations from most to least significant, grouped by metadata
## 2026-03-04 18:37:06.481327 INFO::Plotting data for metadata number 1, agegroup
## 2026-03-04 18:37:06.487719 INFO::Creating boxplot for categorical data, agegroup vs d9d2694a2422852ddb9e381500af7cca
```

```
## 2026-03-04 18:37:06.627495 INFO::Creating boxplot for categorical data, agegroup vs X17c62c1c91771a6d6a11af56744a4812
```

```
## 2026-03-04 18:37:06.765947 INFO::Creating boxplot for categorical data, agegroup vs c7ef7fb22153326e3ebb23932d0c8f8c
```

Save the fit data object as an rds file.


``` r
# 20260224_21_19_15
# saveRDS(maaslin.fit_data.age,
#         file = file.path("output/rdafiles",
#                          paste(
#                            paste(format(Sys.time(),format="%Y%m%d"),
#                                  format(Sys.time(),format = "%H_%M_%S"),
#                                  sep = "_"),"maaslin.fit_data.age",
#                            nmr.age.data.for_test$output.filename,
#                            ".rds",sep = "-")))
# "20260213_14_39_57" for both signif and signif.decreased tsv
maaslin.processed_output.age <- 
  process_maaslin2_output(agglom.rank = agglom.rank,
                          maaslin.fit_data = maaslin.fit_data.age,
                          ps.q.agg = nmr.age.data.for_test$ps.q.agg,
                          sample.groups = nmr.age.data.for_test$sample.groups)
```

Downstream analysis and plotting of differentially abundant features.


``` r
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
```

```
## Joining with `by = join_by(Sample, class)`
```

<img src="007-diffabund-tests_files/figure-html/unnamed-chunk-52-1.png" alt="" width="768" />



## Compare NMR sexes (ASV level).
Choose what to compare:


``` r
comparison<-"sex"
agglom.rank<-"OTU"
```

Choose the reference level:


``` r
ref.level<-"female"
inside.host = TRUE
```

Prepare the data for testing.


``` r
nmr.sex.data.for_test<-prepare_data(agglom.rank = agglom.rank, 
             comparison = comparison, 
             ref.level = ref.level, 
             inside.host = TRUE)
```

```
## [1] "Performing differential microbial abundance tests on sex variable at OTU level. Reference: female"
## [1] "Output filename prefix: NMR-OTU-sex-234-ref-female"
## [1] "female" "male"  
## [1] "Prepared the data for tests"
```

Run MaAsLin2 test.


``` r
maaslin.fit_data.sex<-
  perform_maaslin2_test(ps.q.df.wide = nmr.sex.data.for_test$dataset,
                        comparison = comparison,
                        output.filename = nmr.sex.data.for_test$output.filename,
                        custom.md = nmr.sex.data.for_test$metadata,
                        custom.levels = nmr.sex.data.for_test$custom.levels,
                        ref.level = ref.level,
                        inside.host = inside.host
  )
```

```
## [1] "Creating output folder"
## [1] "Creating output feature tables folder"
## [1] "Creating output fits folder"
## [1] "Creating output figures folder"
## 2026-03-04 18:37:08.638548 INFO::Writing function arguments to log file
## 2026-03-04 18:37:08.679841 INFO::Verifying options selected are valid
## 2026-03-04 18:37:08.683098 INFO::Determining format of input files
## 2026-03-04 18:37:08.686405 INFO::Input format is data samples as rows and metadata samples as rows
## 2026-03-04 18:37:08.714502 INFO::Formula for random effects: expr ~ (1 | relation)
## 2026-03-04 18:37:08.718153 INFO::Formula for fixed effects: expr ~  sex
## 2026-03-04 18:37:08.722598 INFO::Filter data based on min abundance and min prevalence
## 2026-03-04 18:37:08.728134 INFO::Total samples in data: 24
## 2026-03-04 18:37:08.731389 INFO::Min samples required with min abundance for a feature not to be filtered: 0.000000
## 2026-03-04 18:37:08.743764 INFO::Total filtered features: 0
## 2026-03-04 18:37:08.748067 INFO::Filtered feature names from abundance and prevalence filtering:
## 2026-03-04 18:37:08.766231 INFO::Total filtered features with variance filtering: 0
## 2026-03-04 18:37:08.771494 INFO::Filtered feature names from variance filtering:
## 2026-03-04 18:37:08.775028 INFO::Running selected normalization method: TSS
## 2026-03-04 18:37:08.787869 INFO::Bypass z-score application to metadata
## 2026-03-04 18:37:08.792384 INFO::Running selected transform method: LOG
## 2026-03-04 18:37:08.81106 INFO::Running selected analysis method: LM
## 2026-03-04 18:37:08.817074 INFO::Fitting model to feature number 1, X775f1d61e616978239ee36c337b0403c
## 2026-03-04 18:37:08.84811 INFO::Fitting model to feature number 2, X5a818007a6452d5db9223ef1388e451c
## 2026-03-04 18:37:08.881434 INFO::Fitting model to feature number 3, X2a9b76ddc8fc1e5f0e10f1fd40de1c79
## 2026-03-04 18:37:08.913048 INFO::Fitting model to feature number 4, c513b4ef037af47125244161fb1eab50
## 2026-03-04 18:37:08.944032 INFO::Fitting model to feature number 5, X72aa784b6a3f05f45910fe33902caa81
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:08.975342 INFO::Fitting model to feature number 6, X37d6420f6843a99bbb892a3b37178db1
## 2026-03-04 18:37:09.005715 INFO::Fitting model to feature number 7, X5e6767d29a7bbc264a3edb4d4ef313e3
## 2026-03-04 18:37:09.036058 INFO::Fitting model to feature number 8, X087e97a98923593a6ca9d3292593baa8
## 2026-03-04 18:37:09.067677 INFO::Fitting model to feature number 9, a23076ded2ba0960eed8ce5e654dc804
## 2026-03-04 18:37:09.098255 INFO::Fitting model to feature number 10, X8607913ed1bbf9689f69dfccb73812bc
## 2026-03-04 18:37:09.128702 INFO::Fitting model to feature number 11, X6a9f0aa987f51dc6c85f04fb362f2257
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.175119 INFO::Fitting model to feature number 12, X050d316928bb2c2dcdc4b778c5d482af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.205211 INFO::Fitting model to feature number 13, d48f70a0a77c0821bd2dd0c91f355f30
## 2026-03-04 18:37:09.235924 INFO::Fitting model to feature number 14, cbb7a2fe6431a3d874d38a8b6a104854
## 2026-03-04 18:37:09.2663 INFO::Fitting model to feature number 15, X3f058a90b6942f839656eb1dda7d9826
## 2026-03-04 18:37:09.297022 INFO::Fitting model to feature number 16, d3ba622809f71a4ea2725cffa1152398
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.327686 INFO::Fitting model to feature number 17, X4fd27ee04878ce466b07a8fdb6efb91e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.36219 INFO::Fitting model to feature number 18, X10023be92dba55239d6f4b32248afb74
## 2026-03-04 18:37:09.39488 INFO::Fitting model to feature number 19, X22b0f194732516c62929da18c03026fc
## 2026-03-04 18:37:09.425966 INFO::Fitting model to feature number 20, b096b6c93a5f390950dbf69ccedb31e5
## 2026-03-04 18:37:09.457575 INFO::Fitting model to feature number 21, X02e2fceb21147e8bf7ee6a805b54d187
## 2026-03-04 18:37:09.488681 INFO::Fitting model to feature number 22, X949c5cd1c3f0441a3332d7c892064517
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.523693 INFO::Fitting model to feature number 23, c41814121308544fe6e6dc6327604921
## 2026-03-04 18:37:09.553874 INFO::Fitting model to feature number 24, d913c2cbde361858dec499be698ae593
## 2026-03-04 18:37:09.584842 INFO::Fitting model to feature number 25, X5d2a29b5344df06519f7d59f2041da54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.622633 INFO::Fitting model to feature number 26, X72d39f953d79e334e405cc903e866ef6
## 2026-03-04 18:37:09.653107 INFO::Fitting model to feature number 27, f6780cf12baa09329eeff07270ea2e08
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.683365 INFO::Fitting model to feature number 28, X63de5a8a3981bc0f00dbefe46ef47678
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.714359 INFO::Fitting model to feature number 29, X9f2b8085820f945947d0ae63fb7ec339
## 2026-03-04 18:37:09.744498 INFO::Fitting model to feature number 30, d886c00083ac9a3157fcb9dc989197b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.77588 INFO::Fitting model to feature number 31, X6034c86ff9a636506753fccb7cc73615
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.806484 INFO::Fitting model to feature number 32, X8ed8bcc26ac1f36c71117f0def9d0789
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.836813 INFO::Fitting model to feature number 33, X9426f864cbf3629f040d6c7c3ed64180
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.876484 INFO::Fitting model to feature number 34, X85ac797d67074ba62de3413d665354ba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:09.907088 INFO::Fitting model to feature number 35, X2e91b1556d45535a4b49b012917c5b19
## 2026-03-04 18:37:09.947753 INFO::Fitting model to feature number 36, dd4e02d8ebf6639915b74504fa6759a4
## 2026-03-04 18:37:09.980661 INFO::Fitting model to feature number 37, X0ca67ebe8a679e80efb2452016ef7125
## 2026-03-04 18:37:10.011811 INFO::Fitting model to feature number 38, e6289f64c0b74b4cfdf34cebc235b06f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.043509 INFO::Fitting model to feature number 39, dd87bb11986799f8870d596638f99f9e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.074056 INFO::Fitting model to feature number 40, X848bc2dee4d9f89037148afd929c3398
## 2026-03-04 18:37:10.104532 INFO::Fitting model to feature number 41, X4940b5123681ac4251a9c017b731390e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.135561 INFO::Fitting model to feature number 42, X37677efe4c976fd9854262287add28ad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.166607 INFO::Fitting model to feature number 43, X398f3eae476d5c1bbb98af8d94ab3907
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.199597 INFO::Fitting model to feature number 44, X9d2b5cac004beb045c25df1c8a40b544
## 2026-03-04 18:37:10.230275 INFO::Fitting model to feature number 45, c38f378cecf99bde998dee1cd520fed4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.26565 INFO::Fitting model to feature number 46, X76a1f1eb3d45f64d33db506689665b59
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.296784 INFO::Fitting model to feature number 47, X390143a9a62dfa24c965f1c0843f9452
## 2026-03-04 18:37:10.327609 INFO::Fitting model to feature number 48, X35e23a83d3a6bba3e8af7ce15b881dc9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.35873 INFO::Fitting model to feature number 49, X22e069b421275afdc6b4e6cc0950ed4f
## 2026-03-04 18:37:10.388942 INFO::Fitting model to feature number 50, X501bbca65f0eef36751614dabbe08b36
## 2026-03-04 18:37:10.419172 INFO::Fitting model to feature number 51, X64f7efa6ab1e68ed7e0c98fd6ae1c915
## 2026-03-04 18:37:10.449464 INFO::Fitting model to feature number 52, X4bc9f240864b111cdf624408a87a6e1f
## 2026-03-04 18:37:10.480845 INFO::Fitting model to feature number 53, X665096a7e2a12378fcdaaacb3e58042e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.512649 INFO::Fitting model to feature number 54, X18e6ba724c5b286b604c2411d99de87a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.544378 INFO::Fitting model to feature number 55, X50d7799434925279751198338f77130c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.576175 INFO::Fitting model to feature number 56, cde69820d308b87b23fe6e239523e1a6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.608202 INFO::Fitting model to feature number 57, cbefcf8231ca919d37b62749a9881423
## 2026-03-04 18:37:10.64002 INFO::Fitting model to feature number 58, X2c1332de9241f56a6c7f8221758a9adc
## 2026-03-04 18:37:10.670727 INFO::Fitting model to feature number 59, X33b48213d11840fa28edd31fa038979b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.701718 INFO::Fitting model to feature number 60, b8ffd9b9fd2d561816c0c4ba72baaf85
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.759936 INFO::Fitting model to feature number 61, X6fba4056b223f50c3d51622fc1a84da0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.79332 INFO::Fitting model to feature number 62, dae1fc37c9776df8e228a443451de296
## 2026-03-04 18:37:10.82556 INFO::Fitting model to feature number 63, X6b9c74b7afbf96e071014516d8bc97c1
## 2026-03-04 18:37:10.856937 INFO::Fitting model to feature number 64, X9ee68e4ae7ff58ebfcfec4ee721842b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.891641 INFO::Fitting model to feature number 65, X5b956d9720139ff5385b0900d6a6bf28
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.923789 INFO::Fitting model to feature number 66, X311e9ec578b60769eaadf606a38a84cb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.956283 INFO::Fitting model to feature number 67, a58b408564f69ed083451dfacb552a90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:10.98797 INFO::Fitting model to feature number 68, a647fbcd189f56d189342f5ce294cd57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.020042 INFO::Fitting model to feature number 69, X3871a8798bd4cf074b1352e193c6f2e8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.052278 INFO::Fitting model to feature number 70, e7480012a7f751cd862916f0e6d468f7
## 2026-03-04 18:37:11.082882 INFO::Fitting model to feature number 71, a4d9c76f86b8ebf235f18ccbf6f67e86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.114373 INFO::Fitting model to feature number 72, abfe7842b6f39f0f8b0d7c1aebab8294
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.144793 INFO::Fitting model to feature number 73, X463ddfeb23025d3a52776c37581f0d73
## 2026-03-04 18:37:11.174326 INFO::Fitting model to feature number 74, X2ddfa5cd15c7edd114632d1ddbea899d
## 2026-03-04 18:37:11.204942 INFO::Fitting model to feature number 75, e55032c72052adb2f83ce0329ba984ce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.238707 INFO::Fitting model to feature number 76, b1433ade846d530f325945bb135f1bba
## 2026-03-04 18:37:11.275792 INFO::Fitting model to feature number 77, a6332efb97b7d387a65f5650829fb2e6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.306224 INFO::Fitting model to feature number 78, d1c2e6ccadfb83866f23fe93be5fada4
## 2026-03-04 18:37:11.33591 INFO::Fitting model to feature number 79, X50bc369ed002b3668d082f5965556786
## 2026-03-04 18:37:11.370047 INFO::Fitting model to feature number 80, b9c87bd9b886c6c2fc772e34369130b1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.400974 INFO::Fitting model to feature number 81, X8f0018f8ebf455bf3e078ccd7756f454
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.432571 INFO::Fitting model to feature number 82, X5bda8853404fe5cb2a95e1783ca59cd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.471553 INFO::Fitting model to feature number 83, a7110f6b5ef234a8d21e3e5391faee2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.501997 INFO::Fitting model to feature number 84, X09a54c881aae0c9078c3cd72034eeddc
## 2026-03-04 18:37:11.54426 INFO::Fitting model to feature number 85, X867ea3db476d7118238f5b73d76c4841
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.577758 INFO::Fitting model to feature number 86, dd4b4f0ca9ebbbaab73466b57715192e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.609414 INFO::Fitting model to feature number 87, X05625778e3e287a8c0eb4e66e2e066f4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.641187 INFO::Fitting model to feature number 88, X08adbc5295e79ddbd27a1e49b05fc4dc
## 2026-03-04 18:37:11.671986 INFO::Fitting model to feature number 89, X4f71569f1364ca8c24d0971e533bf5f5
## 2026-03-04 18:37:11.70279 INFO::Fitting model to feature number 90, f83157419dd89248f996d2f278fbc713
## 2026-03-04 18:37:11.733175 INFO::Fitting model to feature number 91, X2341971441277c698da922ff13b6cf42
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.763811 INFO::Fitting model to feature number 92, X5cb934450030fdb4f808c8cca52961a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.801306 INFO::Fitting model to feature number 93, X679b1d8a517eb894be6da3eb3197ce90
## 2026-03-04 18:37:11.832183 INFO::Fitting model to feature number 94, X75e19d2211a08be3f623147f55e19aec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.863361 INFO::Fitting model to feature number 95, X10412032e5a4099296e96a115a57c7f6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.895254 INFO::Fitting model to feature number 96, af6b154cb33e86aa5495219241f05b22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:11.925798 INFO::Fitting model to feature number 97, X31947dbf83d9f2c7ae9c2d3d07ad669e
## 2026-03-04 18:37:11.963265 INFO::Fitting model to feature number 98, X85ca038ebdb000cf6701a43615349ea2
## 2026-03-04 18:37:11.993919 INFO::Fitting model to feature number 99, X233d661a3d8de3b4f9a28e926523d4c4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.024622 INFO::Fitting model to feature number 100, X8c6d9cefd98212bf81b189254d48b8f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.061994 INFO::Fitting model to feature number 101, X784386df30006e4485cc62dcf6743234
## 2026-03-04 18:37:12.09325 INFO::Fitting model to feature number 102, X503d066c2ba9bb8c8b2177cc6353501c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.125726 INFO::Fitting model to feature number 103, d4990c02b10ff6cfcb7e0f96c9196a49
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.15604 INFO::Fitting model to feature number 104, X5dff79a69be9a1b82db8a5638a1721a8
## 2026-03-04 18:37:12.187199 INFO::Fitting model to feature number 105, c1ae93de93f91093cffe77b5625360cb
## 2026-03-04 18:37:12.219807 INFO::Fitting model to feature number 106, bcfb05341d821b6242511abff6cf365e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.251946 INFO::Fitting model to feature number 107, X87a9a268b3ae6c7a3e20b3aa3ed33aa6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.283306 INFO::Fitting model to feature number 108, X36ad2071119d1950199ca71b1937b246
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.321829 INFO::Fitting model to feature number 109, X041f8052f648be811c7cc62cd2abf516
## 2026-03-04 18:37:12.371537 INFO::Fitting model to feature number 110, X164b12a7f94fed5ed43389908aab6c81
## 2026-03-04 18:37:12.410344 INFO::Fitting model to feature number 111, X9c3d252c26b173da3d2afd11ffb176f8
## 2026-03-04 18:37:12.450043 INFO::Fitting model to feature number 112, X957b8fefd3f381946d85f43ee6086f74
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.482239 INFO::Fitting model to feature number 113, f13620b4f8abb891cc211b453c38244f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.513782 INFO::Fitting model to feature number 114, X38f6192f7df7849cbe5ce8e5b4aa3af0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.545903 INFO::Fitting model to feature number 115, c0c6de5f96a4bf636b60aff158542063
## 2026-03-04 18:37:12.57705 INFO::Fitting model to feature number 116, X98e5eedb473623066d9edab84e519c0f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.614939 INFO::Fitting model to feature number 117, X9b9f46c24a9200aa740f6c951155126d
## 2026-03-04 18:37:12.651777 INFO::Fitting model to feature number 118, X87531a4d94460fab40ad7b6c875690ca
## 2026-03-04 18:37:12.689131 INFO::Fitting model to feature number 119, e1e8a7389e298f5ed01997c426b8efbb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.719986 INFO::Fitting model to feature number 120, X9e260a0418dbb7f9ecf40402a3453429
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.756646 INFO::Fitting model to feature number 121, X43d3e97de5be62fbc5bef7436779aef0
## 2026-03-04 18:37:12.787572 INFO::Fitting model to feature number 122, X4382fa89de91c2123c51b2c5a3526794
## 2026-03-04 18:37:12.819341 INFO::Fitting model to feature number 123, X7823dfa03b1068e544769d82d5108d3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.856426 INFO::Fitting model to feature number 124, X70f25a4988056da265638c4028e0e9b4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.887236 INFO::Fitting model to feature number 125, X8bdedac0d01bba7276354218f19a8703
## 2026-03-04 18:37:12.917033 INFO::Fitting model to feature number 126, f4a400e4fc836805ab7e7283605760d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.947957 INFO::Fitting model to feature number 127, b6b1ff24f46d79cb2b27464cddc4d695
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:12.978361 INFO::Fitting model to feature number 128, X9f4a338d0544b7ac1c7643304ac451bf
## 2026-03-04 18:37:13.009643 INFO::Fitting model to feature number 129, X546a95e7573470b5fdd2309be2684624
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.040332 INFO::Fitting model to feature number 130, X68b62e646ef6c45a8d31a50e9b3f83e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.070877 INFO::Fitting model to feature number 131, c832b45fd6e005a61178908bb566ca82
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.10071 INFO::Fitting model to feature number 132, X2ba85532d09e46703511d81ccb3a4622
## 2026-03-04 18:37:13.132228 INFO::Fitting model to feature number 133, X7c992c13c3614197b1b331b6182bd3f5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.163763 INFO::Fitting model to feature number 134, X45cbf3ea056aff67a4a75d4e8da686c0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.208577 INFO::Fitting model to feature number 135, X6d60086ac416f97c8d3a774045920413
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.241658 INFO::Fitting model to feature number 136, X587d21cd26f3e0a14be44c7f15b7aaeb
## 2026-03-04 18:37:13.272651 INFO::Fitting model to feature number 137, adcbd92e759e2bd9c0bf798af532d8c3
## 2026-03-04 18:37:13.30325 INFO::Fitting model to feature number 138, a30e54c24ebcf6bc5e5a1ea7b868934d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.33387 INFO::Fitting model to feature number 139, f06b506d343103b6969101eb88e864ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.365557 INFO::Fitting model to feature number 140, X69edbd69fd17306596e6f0765d8baaf5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.396704 INFO::Fitting model to feature number 141, e6d608aa96cffb876fdcb1c5e610454f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.429181 INFO::Fitting model to feature number 142, a1f9f5146e383650703e9ba1718c81ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.461222 INFO::Fitting model to feature number 143, X938b0c38f033d20f1b16c42c9fb10a36
## 2026-03-04 18:37:13.492188 INFO::Fitting model to feature number 144, X0edc4ea81df8083fcd669e4d24f5c34c
## 2026-03-04 18:37:13.522464 INFO::Fitting model to feature number 145, bf4b19f3a06f13e371c47798b759f967
## 2026-03-04 18:37:13.55527 INFO::Fitting model to feature number 146, X5a0cb8bb500ca6731dc4283cba2aaeaf
## 2026-03-04 18:37:13.585749 INFO::Fitting model to feature number 147, c77b5d1e3caa62de4030ae1460352956
## 2026-03-04 18:37:13.616064 INFO::Fitting model to feature number 148, f71741d8f1751a7511fbd4455a86593b
## 2026-03-04 18:37:13.646846 INFO::Fitting model to feature number 149, da3141bea3827e8b54dfc7403b6f3ab9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.67783 INFO::Fitting model to feature number 150, X3b86c743ad3ceff86644bc3ce0d66a16
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.708456 INFO::Fitting model to feature number 151, X7823c990b48e0d5fdbb67b0cf1511b0d
## 2026-03-04 18:37:13.742578 INFO::Fitting model to feature number 152, caae62e7d284909f0ba194e383c99196
## 2026-03-04 18:37:13.772615 INFO::Fitting model to feature number 153, b2fcbcf8b5b2ee9ae672857fff1c4ae4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.802928 INFO::Fitting model to feature number 154, X002f3595b8a20738572f67dfcaf531a1
## 2026-03-04 18:37:13.832813 INFO::Fitting model to feature number 155, X02031267bf033fd8b5bc46043628acd1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.863439 INFO::Fitting model to feature number 156, X20ab8b9e32bd6d9bfc506b65fd696ef9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.894424 INFO::Fitting model to feature number 157, X587b64913a785e4d7cc810c80b2a9722
## 2026-03-04 18:37:13.926094 INFO::Fitting model to feature number 158, f41431a72bfd93b6f2c4a49ea7642ae8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:13.985032 INFO::Fitting model to feature number 159, X8ab93fef97076ce48913426901d23711
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.018585 INFO::Fitting model to feature number 160, X8dca4bfdcbcbb43dcbbf5a6c7943b833
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.051431 INFO::Fitting model to feature number 161, X7fb1c8722f2ebb009c074460b08d5b5f
## 2026-03-04 18:37:14.088775 INFO::Fitting model to feature number 162, X27b102c4d0f321ccb73b4fdd3cf04efa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.126843 INFO::Fitting model to feature number 163, dd3c9baa3973a152a2bada71543454f5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.163958 INFO::Fitting model to feature number 164, X5ae982af003eaae773f2ac4ac8d1f518
## 2026-03-04 18:37:14.1954 INFO::Fitting model to feature number 165, X52f10c839bad102da0f9182fe2a65d63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.230592 INFO::Fitting model to feature number 166, X4cf488fc55937cceb0cf0d50c67f78c7
## 2026-03-04 18:37:14.262406 INFO::Fitting model to feature number 167, X0e7e09e2f26a49f7e392c9a979f935d0
## 2026-03-04 18:37:14.292428 INFO::Fitting model to feature number 168, b9ce57fad79d507e7f193ab62a1d590f
## 2026-03-04 18:37:14.323549 INFO::Fitting model to feature number 169, X1f1f739ae6bb988057eb2f9fd17f3131
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.35719 INFO::Fitting model to feature number 170, fc47c350e68145c3f6a5dc0a818fb683
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.394399 INFO::Fitting model to feature number 171, X3fb50d74ec763a9a895b6a69751056ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.43123 INFO::Fitting model to feature number 172, eee7a5d0e69d4129a7cd4fd8c7a6abea
## 2026-03-04 18:37:14.462549 INFO::Fitting model to feature number 173, c9e692a4f678ce3e9a9bf69290382934
## 2026-03-04 18:37:14.492854 INFO::Fitting model to feature number 174, a8d690ae47779a14d63923f302f23b44
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.529369 INFO::Fitting model to feature number 175, X4989b6aa56952d17a9728cb74ba4eccd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.559435 INFO::Fitting model to feature number 176, ac5128abccedf13cded0c2d74b428936
## 2026-03-04 18:37:14.591462 INFO::Fitting model to feature number 177, bf5cec3c7cd23cf88fedb590337edeb7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.628949 INFO::Fitting model to feature number 178, X196dd825ebed8d5e0372081332465e46
## 2026-03-04 18:37:14.66383 INFO::Fitting model to feature number 179, X1664dd935ff0a38d0bc5b804f2c5f6aa
## 2026-03-04 18:37:14.697383 INFO::Fitting model to feature number 180, a32ebbe366e9cd264ffb16a15d10dcb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.728096 INFO::Fitting model to feature number 181, X30803d11bc6be93d144583cd061c1434
## 2026-03-04 18:37:14.760934 INFO::Fitting model to feature number 182, X534c9f3c39ff74b0751677a74ca2223f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.791164 INFO::Fitting model to feature number 183, fcda29c5385d6115d10b5e976e087ae1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.835154 INFO::Fitting model to feature number 184, X8a238b6279fbb51f2a43b9c11a1e7293
## 2026-03-04 18:37:14.867111 INFO::Fitting model to feature number 185, X5e489e67242303169c8781859aa2ea64
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:14.89891 INFO::Fitting model to feature number 186, X7109c3732f21db3d8ca19a4be2231a80
## 2026-03-04 18:37:14.936963 INFO::Fitting model to feature number 187, X7af0ca6b149c28ed1da4c04b112c4551
## 2026-03-04 18:37:14.97438 INFO::Fitting model to feature number 188, X4404bd2d87786e911fba13e90fd51988
## 2026-03-04 18:37:15.008871 INFO::Fitting model to feature number 189, X4ce796eca0eb0252c01451f0acfd7c86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.042904 INFO::Fitting model to feature number 190, ff162ce0992d8e762f0f09422814ed01
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.075072 INFO::Fitting model to feature number 191, X0c0189fefd69f40756d340e48d5f7304
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.107707 INFO::Fitting model to feature number 192, X5dfc13e1239763278b4edaf53ab409b5
## 2026-03-04 18:37:15.155407 INFO::Fitting model to feature number 193, X9c9d04df1b23dcd81f5bbdb09d422af4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.201086 INFO::Fitting model to feature number 194, X8585426a699ce2bd5bdd80954f4df256
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.250354 INFO::Fitting model to feature number 195, f1ca62a3e8e5e041350db2c975346cb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.296752 INFO::Fitting model to feature number 196, X59f4de3950a671dd935530636121f79a
## 2026-03-04 18:37:15.346292 INFO::Fitting model to feature number 197, X8dee14385c53c6adc0b6ed841329b318
## 2026-03-04 18:37:15.391516 INFO::Fitting model to feature number 198, fb6e884f98ccf551138d94557d49b6c4
## 2026-03-04 18:37:15.428498 INFO::Fitting model to feature number 199, d749b2d647bb124c92604b1bdb7c6511
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.473399 INFO::Fitting model to feature number 200, X318bbf7176eb37d514c95f727b120e3f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.513449 INFO::Fitting model to feature number 201, X53f63e6d0ded3831371b61c05e0d1df9
## 2026-03-04 18:37:15.554206 INFO::Fitting model to feature number 202, X810c56e7a79638b2fdf56ece16b9a32c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.594981 INFO::Fitting model to feature number 203, e3ca497d1ddf9a2efda5736924232efe
## 2026-03-04 18:37:15.637592 INFO::Fitting model to feature number 204, X5453c211a413666e9c72d8a12bd3bf2f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.678545 INFO::Fitting model to feature number 205, X6dbf7882c3e44b3f93cd770a5af49299
## 2026-03-04 18:37:15.722275 INFO::Fitting model to feature number 206, X62e6d4459d864709f6413926e3a88095
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.758477 INFO::Fitting model to feature number 207, b4a791cce24c1baaf3966a1d9efd18bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.793449 INFO::Fitting model to feature number 208, X4c944499b014eb72118a65700a136719
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.843248 INFO::Fitting model to feature number 209, X278a0e8506254a8e6bd3ad6e37df2e1d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.878877 INFO::Fitting model to feature number 210, X406caa0d474693f37020bc3b3e99467d
## 2026-03-04 18:37:15.915928 INFO::Fitting model to feature number 211, X580f23c911e277aeb396a2dd45af2675
## 2026-03-04 18:37:15.951794 INFO::Fitting model to feature number 212, X474ee1999580bb7f367858c97607d986
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:15.986451 INFO::Fitting model to feature number 213, X72ef153369491e40b5b0405e442bb1d4
## 2026-03-04 18:37:16.019852 INFO::Fitting model to feature number 214, X42e99751e005f348bbb3a852815e7dab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.05293 INFO::Fitting model to feature number 215, X796e01ba695bbc78aaa63196bcd06330
## 2026-03-04 18:37:16.08638 INFO::Fitting model to feature number 216, b1258e4ee5ca06eb4d5b62a466e413ff
## 2026-03-04 18:37:16.118853 INFO::Fitting model to feature number 217, c0248093940670814ffe413551c08592
## 2026-03-04 18:37:16.150645 INFO::Fitting model to feature number 218, cb5947bc03ed14060ebce14a528a30a9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.182731 INFO::Fitting model to feature number 219, X7c698019fe98da258752a851be4f0dbf
## 2026-03-04 18:37:16.215201 INFO::Fitting model to feature number 220, f819c2a6d74d9e73b54d2b2d782b42ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.248263 INFO::Fitting model to feature number 221, X91aaea5b58ac19ee8a8bde60f2368289
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.279638 INFO::Fitting model to feature number 222, X7d3ca27433dbfec1d46096d6f472c0e6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.310655 INFO::Fitting model to feature number 223, X46bb1d815c069f361d1832faf2e2e4dd
## 2026-03-04 18:37:16.341459 INFO::Fitting model to feature number 224, ae38f6e010625cb3b1b6cea1eea13a57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.372478 INFO::Fitting model to feature number 225, X6bca7aa2b33b0274c340d4facd92afaa
## 2026-03-04 18:37:16.40962 INFO::Fitting model to feature number 226, f430ec182a2fba6dffb1534634aba940
## 2026-03-04 18:37:16.45339 INFO::Fitting model to feature number 227, X9a67c9da63556d4863bab09cafd40940
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.486375 INFO::Fitting model to feature number 228, e92df5b83e69c4ed9bb9d268ecb358e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.518722 INFO::Fitting model to feature number 229, X1fa6739ed5fc2a89235fb28d33d3ccb3
## 2026-03-04 18:37:16.55023 INFO::Fitting model to feature number 230, X1d93a791c26759878a10132d5e6e33fa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.582088 INFO::Fitting model to feature number 231, c67ed17a9beed3c005cf97f0e0f20a18
## 2026-03-04 18:37:16.613613 INFO::Fitting model to feature number 232, db3aca72ce7e5948035457f7871945a2
## 2026-03-04 18:37:16.659013 INFO::Fitting model to feature number 233, a769521df2160b6d284dfd3e1d449c9d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.692114 INFO::Fitting model to feature number 234, X7a21b6c3fa63544dd1cde940adf52388
## 2026-03-04 18:37:16.72319 INFO::Fitting model to feature number 235, X882e40e7940f3b6e5c6614a047e97dc6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.754745 INFO::Fitting model to feature number 236, b9846f46f48c55df07601aaef88e17cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.78569 INFO::Fitting model to feature number 237, X920d9d9fa0d42f62df13a5f61eba6019
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.816088 INFO::Fitting model to feature number 238, bad8cf42670e052284875915a6e94c58
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.847447 INFO::Fitting model to feature number 239, X27b1308239803f1111e6871b0f7733d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.878537 INFO::Fitting model to feature number 240, fc4c1b4d08d17bc8bbbfeccfc267f51a
## 2026-03-04 18:37:16.916085 INFO::Fitting model to feature number 241, X54af66f5c264f0024789d99c4c7afe8d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.947664 INFO::Fitting model to feature number 242, X5e4b3686f1cefe645291fca5b936ceb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:16.979354 INFO::Fitting model to feature number 243, X82df3e7e293f65085e273db27eac0658
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.011637 INFO::Fitting model to feature number 244, X2419595d114a6b8b454c88b03531471e
## 2026-03-04 18:37:17.042574 INFO::Fitting model to feature number 245, X8bcc8165e984ed5f75e38224a26c4500
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.073232 INFO::Fitting model to feature number 246, fed6e1a07e0a36780b2112ad59790baf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.104591 INFO::Fitting model to feature number 247, X3a81c2a8bcd8262593a6a28835b1d5d1
## 2026-03-04 18:37:17.139401 INFO::Fitting model to feature number 248, X07aa13125fef95029e4e578ff8d77194
## 2026-03-04 18:37:17.170209 INFO::Fitting model to feature number 249, X2ee55325dd6fff97f7e798daf0a1a173
## 2026-03-04 18:37:17.200681 INFO::Fitting model to feature number 250, cdda58bad36d42fffe349c500561263a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.231409 INFO::Fitting model to feature number 251, d8bc5ba9f36585e295c58310928b7610
## 2026-03-04 18:37:17.262709 INFO::Fitting model to feature number 252, X2443c7149b3c68f110a91c9464cd3a9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.293222 INFO::Fitting model to feature number 253, d3ce175f14e9bd4ce6fa566098070432
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.323934 INFO::Fitting model to feature number 254, f5852b6815968593b30bfc519b224b13
## 2026-03-04 18:37:17.354023 INFO::Fitting model to feature number 255, X96cec0b11138b0ddb91bc8fd3a827e63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.384332 INFO::Fitting model to feature number 256, faecb8ca1ca70b02713116054ba6a25e
## 2026-03-04 18:37:17.41546 INFO::Fitting model to feature number 257, X428241385b3872803e4e6f09856f1eaf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.460938 INFO::Fitting model to feature number 258, f7def96d9ec34bbc94788d5096e6b653
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.492484 INFO::Fitting model to feature number 259, X5bdcbe959f1cc3b608874a0899b2b431
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.523286 INFO::Fitting model to feature number 260, X5df10ae748068240d0954e709aa85ef6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.554074 INFO::Fitting model to feature number 261, X773a6686911fcf580df6cc6ce27e25a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.584655 INFO::Fitting model to feature number 262, X962f048f3548d862a1fa9bfef48d3c8d
## 2026-03-04 18:37:17.615633 INFO::Fitting model to feature number 263, X10d3555fce8c1645f69075680f9644e6
## 2026-03-04 18:37:17.646166 INFO::Fitting model to feature number 264, X5eb1bd896aeef657980bb6319e0fa46c
## 2026-03-04 18:37:17.67607 INFO::Fitting model to feature number 265, a1163d29062c052bb3d59c49fdc0c44a
## 2026-03-04 18:37:17.706092 INFO::Fitting model to feature number 266, X56ee0e676dc712ca9636f01f082d12e1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.736583 INFO::Fitting model to feature number 267, e2773b5fd096352ad027595baeab01ef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.767076 INFO::Fitting model to feature number 268, X66eecd48e1caa9adf86f0ef816e22d87
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.797756 INFO::Fitting model to feature number 269, cfc06d9223155de233be7daf646b0cfe
## 2026-03-04 18:37:17.827901 INFO::Fitting model to feature number 270, X741671cead50e930b5ae2fdba6f09ac2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.858236 INFO::Fitting model to feature number 271, X85124c575561a28bdd0e3985153d33eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.889096 INFO::Fitting model to feature number 272, X8c0a48a6da32d89300557a08eb584eab
## 2026-03-04 18:37:17.919462 INFO::Fitting model to feature number 273, X3de3d0c798ac11c7eb2d3cd2c5eb6d22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.950553 INFO::Fitting model to feature number 274, fdc6ffd6aff43d49592d12868ef6edda
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:17.981465 INFO::Fitting model to feature number 275, X175e258fe4d2966115d9f5607a6dc868
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.013265 INFO::Fitting model to feature number 276, X2adacb50dbacd8fa223cd5107f430319
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.045367 INFO::Fitting model to feature number 277, X4c09b59637d648bb9d4cead257187763
## 2026-03-04 18:37:18.083169 INFO::Fitting model to feature number 278, X6314d0517e01bcc3d8907f875891bddd
## 2026-03-04 18:37:18.113225 INFO::Fitting model to feature number 279, X001ce79a8f8d035a30c8203a7d76b05e
## 2026-03-04 18:37:18.143761 INFO::Fitting model to feature number 280, X30b0cc0ca2a2da9533100ba3491175d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.17454 INFO::Fitting model to feature number 281, d0750756e3311a835281cad227e6c4a0
## 2026-03-04 18:37:18.216458 INFO::Fitting model to feature number 282, e72dbc2cbbcaf8d32bf52077e98fb9d6
## 2026-03-04 18:37:18.252075 INFO::Fitting model to feature number 283, X0fdd81e04676039a9b2bcbcf6ac166ff
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.283326 INFO::Fitting model to feature number 284, X6d729250b3b458887e3362f17fb6b54b
## 2026-03-04 18:37:18.314593 INFO::Fitting model to feature number 285, X58dfc0bd94f7612268c6e40fa4819fc5
## 2026-03-04 18:37:18.344929 INFO::Fitting model to feature number 286, X62111f7adcf73e318470d3c36ac127a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.375622 INFO::Fitting model to feature number 287, a098856e088a414969785c312c7077dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.406172 INFO::Fitting model to feature number 288, f7db7195338b7afd38f1f529a1d6a8c2
## 2026-03-04 18:37:18.438282 INFO::Fitting model to feature number 289, X5a101d21f512bd402ba5dce9ff8853b3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.469685 INFO::Fitting model to feature number 290, X722c3b22b71c06ae914fdc65c8a93e1b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.500823 INFO::Fitting model to feature number 291, X28d37664dfaf5d0318674d3903759319
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.53203 INFO::Fitting model to feature number 292, d9d2694a2422852ddb9e381500af7cca
## 2026-03-04 18:37:18.563215 INFO::Fitting model to feature number 293, X3d9a5be20ced9a47b9cdba62499e9b67
## 2026-03-04 18:37:18.593404 INFO::Fitting model to feature number 294, X49cddebff224d6008280a502f630d0dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.624265 INFO::Fitting model to feature number 295, X8adece3b43b1d67576ec24a05663f526
## 2026-03-04 18:37:18.655478 INFO::Fitting model to feature number 296, X2f748f1a2371911911065af675ff0818
## 2026-03-04 18:37:18.685793 INFO::Fitting model to feature number 297, b985a8c66261968b9d8a329d29755c8c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.718072 INFO::Fitting model to feature number 298, c77e2df1b5831857b410fd099900e05a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.749365 INFO::Fitting model to feature number 299, b8c05b8dc05f4f544735add1f5f14ee9
## 2026-03-04 18:37:18.788183 INFO::Fitting model to feature number 300, e5aa8ee9cb7d24ef5df26f0bc2a72299
## 2026-03-04 18:37:18.818854 INFO::Fitting model to feature number 301, c47ab297ca5c2e467fb15ca968fbf667
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.849624 INFO::Fitting model to feature number 302, X48da3e64315f06b66b46c4a7caa33b26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.880817 INFO::Fitting model to feature number 303, X58a318efa7e9171f7b768500c10b0a85
## 2026-03-04 18:37:18.911475 INFO::Fitting model to feature number 304, X680fcfeaf91f648d224c23cb0c49e41a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.943987 INFO::Fitting model to feature number 305, X2d3e451e9169d654b061d02d71d3b319
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:18.9829 INFO::Fitting model to feature number 306, X7df3c6b0b97f9d38b18f47462e2eb666
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.028718 INFO::Fitting model to feature number 307, fb9c18e59bc4a13d9c4b1830cc979ba8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.063143 INFO::Fitting model to feature number 308, X6bae82c66563cfba208caab642e21b9c
## 2026-03-04 18:37:19.094239 INFO::Fitting model to feature number 309, feaf3412296cd21fe71e288500f2909d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.125359 INFO::Fitting model to feature number 310, X2daaed7767de192f02b5ef6aaac2fd45
## 2026-03-04 18:37:19.156456 INFO::Fitting model to feature number 311, X88f2e5b2a5d102be5d9d22e56c5b5033
## 2026-03-04 18:37:19.187546 INFO::Fitting model to feature number 312, ab0a47c0d634005b50d204cc87f0300f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.218548 INFO::Fitting model to feature number 313, X2b8c21681ff428a1c834dd9b1fd35d4f
## 2026-03-04 18:37:19.249944 INFO::Fitting model to feature number 314, X590139931aa19c8ed5551ec58c8fc556
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.280989 INFO::Fitting model to feature number 315, X1958f1ea043fda09b0846074b8cf3387
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.312195 INFO::Fitting model to feature number 316, bcfaaaae2ef79153b03ef0814d0d7abf
## 2026-03-04 18:37:19.34281 INFO::Fitting model to feature number 317, ada302e958702d87a162ac13e5225fd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.376527 INFO::Fitting model to feature number 318, d733da4c4eb11b498ccf79ced6da79d2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.408894 INFO::Fitting model to feature number 319, X1142cc5fdfb17bf3dda55341b1fc2ece
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.440647 INFO::Fitting model to feature number 320, X16a4be5bfb2132b4cf9b1ce6753abbf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.470996 INFO::Fitting model to feature number 321, X1daf20c298f03dfd4551813455d60a75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.501117 INFO::Fitting model to feature number 322, X5d982d4a4a7f3ce92e0e0dbfe0008be4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.531631 INFO::Fitting model to feature number 323, X6f3b3e805737133d5a15e6150096b7ad
## 2026-03-04 18:37:19.561483 INFO::Fitting model to feature number 324, a63f57972cf38df045ae50390d88758b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.59157 INFO::Fitting model to feature number 325, X21df1e53f11d76242c275d0844a43fdf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.622562 INFO::Fitting model to feature number 326, X7880d3e8bcca1feb2db8b9839741438f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.653841 INFO::Fitting model to feature number 327, ebd2eb0273b048077d97a8e5c87843be
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.685161 INFO::Fitting model to feature number 328, X17c62c1c91771a6d6a11af56744a4812
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.719797 INFO::Fitting model to feature number 329, X4e2b574c5c896013c886dcccfe212bb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.751377 INFO::Fitting model to feature number 330, c351ea5b351f3da50f771e157e96c96f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.799471 INFO::Fitting model to feature number 331, X59a9c9743d0de49303aa710f7f2dcf1a
## 2026-03-04 18:37:19.836335 INFO::Fitting model to feature number 332, bdb97410884f71c5be50c8881e2c98b0
## 2026-03-04 18:37:19.867929 INFO::Fitting model to feature number 333, a36840adf6ca830c32cba3a0bf05c24c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.899349 INFO::Fitting model to feature number 334, da20e0254113aed621b9f35dd3b17a74
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.931016 INFO::Fitting model to feature number 335, X4cf353124acf76a83451339182554a2c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.962579 INFO::Fitting model to feature number 336, d699ad5c7d262166690df8e94610d11d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:19.994304 INFO::Fitting model to feature number 337, X0e95449d477a09d2ab60aad4e716e758
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.025257 INFO::Fitting model to feature number 338, c7a4a2780d280e743c16f5be29834ceb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.057855 INFO::Fitting model to feature number 339, d26a840ed062289b08114e43ed383a08
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.089472 INFO::Fitting model to feature number 340, d89f4af0cc9704172e5a75524368626a
## 2026-03-04 18:37:20.122083 INFO::Fitting model to feature number 341, cd007c70afbc5e207cbd6979ffc64409
## 2026-03-04 18:37:20.154369 INFO::Fitting model to feature number 342, d0050df77ec59e70e009b5af346e275f
## 2026-03-04 18:37:20.186249 INFO::Fitting model to feature number 343, f148106f9939021fa196057f3c53f008
## 2026-03-04 18:37:20.219352 INFO::Fitting model to feature number 344, X6b8516936fd7f3084db293f299865b97
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.251576 INFO::Fitting model to feature number 345, a3cf3667e4778d77c70efcf821731124
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.283921 INFO::Fitting model to feature number 346, X268786e064d71003f5ae04c6bc980978
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.315342 INFO::Fitting model to feature number 347, X3e398419b48fe7fd6882bfc636293558
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.346257 INFO::Fitting model to feature number 348, b5261b3c77cd990eef2c90209b8fd44f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.377178 INFO::Fitting model to feature number 349, cd892c02a21141c715fdef5c234c791c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.408467 INFO::Fitting model to feature number 350, d417a886163ff1310f4ff6be7122a500
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.440641 INFO::Fitting model to feature number 351, faa27d3e3fc59a116b955ca7bad4567d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.472507 INFO::Fitting model to feature number 352, X30b59313f4ae8d27dd90ead569dd220c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.505067 INFO::Fitting model to feature number 353, X9fe431eb75ff0fb76d1b3597e7c6e018
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.537367 INFO::Fitting model to feature number 354, c5330203390f5f92c8631531cd1929dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.57591 INFO::Fitting model to feature number 355, X37ab667cfa248f97195af5b3622bfc11
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.62217 INFO::Fitting model to feature number 356, X480c1786c18bdbaa003293d15e5aaf50
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.655959 INFO::Fitting model to feature number 357, ebfa013a17853303ce8dcb5b81825b1f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.687411 INFO::Fitting model to feature number 358, fee34b2a495eca65a0f5a95ba6f8348c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.71985 INFO::Fitting model to feature number 359, X44333eb795d84f9c6d46999cc2e9cbb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.752606 INFO::Fitting model to feature number 360, X61f6666d6fb9564cd5c8e249386bb508
## 2026-03-04 18:37:20.785265 INFO::Fitting model to feature number 361, ebfc68839dcbe864af2c9e957124a850
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.817707 INFO::Fitting model to feature number 362, X44f394b4abd8f4a0b77f5d918b2f1a3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.849877 INFO::Fitting model to feature number 363, e770697660b66ac58e6264e25fe09caf
## 2026-03-04 18:37:20.889459 INFO::Fitting model to feature number 364, X3a08076d3ca27f6a4e23b428c227f723
## 2026-03-04 18:37:20.921471 INFO::Fitting model to feature number 365, X7d9f76c17978c3357bc6e5a39c3e5c15
## 2026-03-04 18:37:20.953796 INFO::Fitting model to feature number 366, c74b3b8f79a2710cd2b055e35178c301
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:20.986824 INFO::Fitting model to feature number 367, d67eb9feb3c3e913b0eda1a3ccfd6876
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.021005 INFO::Fitting model to feature number 368, e57ea9765333a478b6d4a150a0f4759c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.05401 INFO::Fitting model to feature number 369, X008d766cc2fe0e9a1ff0df6a107b02af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.088069 INFO::Fitting model to feature number 370, X17175a3de507f921278c078d71de6978
## 2026-03-04 18:37:21.118995 INFO::Fitting model to feature number 371, X42db06a32d44f7e0b04833ac260cd01f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.150914 INFO::Fitting model to feature number 372, X53055f08d7d20cfe83c3e1ba05b092d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.187179 INFO::Fitting model to feature number 373, X9daace87b866febfba2396dea07f2e3b
## 2026-03-04 18:37:21.219849 INFO::Fitting model to feature number 374, X4e159cc1d60cb0bdf1d6504d2446411b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.251659 INFO::Fitting model to feature number 375, X6ceba268b6002d018f5be12e1d30ca48
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.28445 INFO::Fitting model to feature number 376, X917d19bb0a01a3511a262d04c27211fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.316527 INFO::Fitting model to feature number 377, c7ef7fb22153326e3ebb23932d0c8f8c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.349522 INFO::Fitting model to feature number 378, X5664cc65cd15a2a34cae8203f26b6a1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.381731 INFO::Fitting model to feature number 379, b11e31f8a98617d9b537f565aa9d5cd0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.41437 INFO::Fitting model to feature number 380, X0e64362e5d9e8e0393a2a25ef51344d8
## 2026-03-04 18:37:21.461904 INFO::Fitting model to feature number 381, X3f6a3469b81ee3684014d0fd4ca8ba22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.496312 INFO::Fitting model to feature number 382, X84ded3c82c02986d17ba89efc005ca1b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.529394 INFO::Fitting model to feature number 383, X7128c784b9f14437017f065ccd13d93e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.562757 INFO::Fitting model to feature number 384, a2f1c3dabb5720aa985c1c56c785a79e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.595509 INFO::Fitting model to feature number 385, c3a3131987eb1c89019993c61051a9b7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.627862 INFO::Fitting model to feature number 386, da69dbfd52386c4e8bec85dcf31fab91
## 2026-03-04 18:37:21.659685 INFO::Fitting model to feature number 387, X00bec85d4cac8001d056f3641e4fe4d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.692005 INFO::Fitting model to feature number 388, X2c756f3da2bf4aa3ba4cb108bfe02ef2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.723368 INFO::Fitting model to feature number 389, d564fadbbec41c3174402b6ed7c2b4e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.756987 INFO::Fitting model to feature number 390, eb4dffb13b2c5f8569c37ff99084fede
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.789018 INFO::Fitting model to feature number 391, ebe44d95e0c653b8424e7f4522c4bff3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.821821 INFO::Fitting model to feature number 392, fa81c053b819a9eeadc0f0306c07528b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.855077 INFO::Fitting model to feature number 393, X0538a2fd5f6431b8728a0e6770ffa2df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:21.887241 INFO::Fitting model to feature number 394, a32076c3bf3b662e045a8b42816f37bf
## 2026-03-04 18:37:21.917734 INFO::Fitting model to feature number 395, dba2b43ba071ef5f2ab8490e283e56e4
## 2026-03-04 18:37:21.94995 INFO::Fitting model to feature number 396, X8f7b886fa9032c9830c8a4232c454ab5
## 2026-03-04 18:37:21.982082 INFO::Fitting model to feature number 397, aa45112add34677769b6912d1749b794
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.014279 INFO::Fitting model to feature number 398, f070e6f921a4473b49d198509117a125
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.046731 INFO::Fitting model to feature number 399, X145a90b7279b2987c8f0eb9c4edc9ea6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.078199 INFO::Fitting model to feature number 400, X4e04424d4099af6d663496c7a432bfc3
## 2026-03-04 18:37:22.109076 INFO::Fitting model to feature number 401, bcfc621313b4ba3bafd334ad6f2a10f5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.144474 INFO::Fitting model to feature number 402, X3d8b735fffa6bb77e92fe95a8042fa00
## 2026-03-04 18:37:22.175328 INFO::Fitting model to feature number 403, X5ad134b39170fe286ea186ba2604897e
## 2026-03-04 18:37:22.206116 INFO::Fitting model to feature number 404, c41bf5e56e53e0e28470374de88e6ba7
## 2026-03-04 18:37:22.248473 INFO::Fitting model to feature number 405, X1ecd1c27ef8548472bc31b8d51d8c4dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.282237 INFO::Fitting model to feature number 406, X7120bc0c2afde34437f0970446441b39
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.313781 INFO::Fitting model to feature number 407, X78c42a2db92f4e3f05634c101cc11259
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.343662 INFO::Fitting model to feature number 408, e2041d00358f2c29a5ef9809e8b1ba96
## 2026-03-04 18:37:22.373843 INFO::Fitting model to feature number 409, X6d6787eb066d8db71d9e3631225a8904
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.40512 INFO::Fitting model to feature number 410, X6fa8ac18801d6dfe5d45d7336d9f5887
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.436979 INFO::Fitting model to feature number 411, X76e9b07ccc3729985a554bc91e58bbc0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.469385 INFO::Fitting model to feature number 412, e188b87dc228e850ec7535fba9be379f
## 2026-03-04 18:37:22.501707 INFO::Fitting model to feature number 413, X956a5b207bb050c7ce43853ce336d5e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.532845 INFO::Fitting model to feature number 414, X174c5b603cffda455b1a9969d47dbaae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.564702 INFO::Fitting model to feature number 415, X27a528bc7a1cdd2669ae30690d420170
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.597968 INFO::Fitting model to feature number 416, X27ef3adf14f05f2c77b6c64118cddb51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.629178 INFO::Fitting model to feature number 417, X7a92bc426ff0b6b2428d8b843bca7b6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.660202 INFO::Fitting model to feature number 418, efa9ff2bcd8d4e435d7982e3e317de43
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.690941 INFO::Fitting model to feature number 419, X39be2b3eb5948a1b4c0e9f823eba63d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.721872 INFO::Fitting model to feature number 420, X726ad819f9c01b960e30637dd15b78c0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.760465 INFO::Fitting model to feature number 421, X9286b1ca453940aa45d397ab3e2058f1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.790652 INFO::Fitting model to feature number 422, X9cc11a19e2926c7037eff26c869725fc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.820896 INFO::Fitting model to feature number 423, a29e744251cd4a2f29279c6aa377d8da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.85238 INFO::Fitting model to feature number 424, b511ce605bde5b3c893d53d043c00b2f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.884468 INFO::Fitting model to feature number 425, c16b00d5bab290920ae7d2d0affe4ca8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.920856 INFO::Fitting model to feature number 426, X1ab6690e22c490e00ad8c2ac9bf24596
## 2026-03-04 18:37:22.954331 INFO::Fitting model to feature number 427, cadb8b4f0b3bfaa180c90a78357042c9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:22.989965 INFO::Fitting model to feature number 428, d6c5071ba31869771d7462efa2058e02
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.028097 INFO::Fitting model to feature number 429, e9746dcc794f4362e6690657f6a4d2f4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.079775 INFO::Fitting model to feature number 430, X0836d28e24ff3d9df1d6dab03c659c54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.11501 INFO::Fitting model to feature number 431, X2b989c0ee8f9bb3ad5b1f572948e6ebc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.153324 INFO::Fitting model to feature number 432, X41813981fc56ece9178aba2271d500ed
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.18734 INFO::Fitting model to feature number 433, X29983e552e325387562ed544443fa784
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.219827 INFO::Fitting model to feature number 434, X453de2772d9c68cfaa47792539745455
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.250941 INFO::Fitting model to feature number 435, X7e947a204620ae22a12b10f2007e3bfd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.282843 INFO::Fitting model to feature number 436, X8cee7e94f2af5f8c3e4def902cb4a4e3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.314319 INFO::Fitting model to feature number 437, a45df882254eed21da6066c7ad4e7664
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.345391 INFO::Fitting model to feature number 438, fee2f4467e991748a4950fb60d0eccda
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.376994 INFO::Fitting model to feature number 439, X77f68b854b479c56f7cd35baa1bbb572
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.408015 INFO::Fitting model to feature number 440, fbcb14cd27216d317b57152441189a35
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.439502 INFO::Fitting model to feature number 441, X65441a076e1d645106982a3e2423e0b3
## 2026-03-04 18:37:23.471193 INFO::Fitting model to feature number 442, X08e77766d42e5c1a91eb3a03251301e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.504822 INFO::Fitting model to feature number 443, X1d7ac3928596c06f8fc853280fd38fb8
## 2026-03-04 18:37:23.536688 INFO::Fitting model to feature number 444, e59f0d7208f53430a8d7ac10438fd9cc
## 2026-03-04 18:37:23.573368 INFO::Fitting model to feature number 445, b0c3763d12e20e512b26e2db2e3d6746
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.606813 INFO::Fitting model to feature number 446, X830c92e2d730f1ac8f4e9435bf6a38a5
## 2026-03-04 18:37:23.639927 INFO::Fitting model to feature number 447, X185dbca6abbfd276dc552de656dcf054
## 2026-03-04 18:37:23.672397 INFO::Fitting model to feature number 448, cf7fa17e32f98d8284e9100e9091ed87
## 2026-03-04 18:37:23.703823 INFO::Fitting model to feature number 449, X85ae192e589c308e221bea6cd2b908eb
## 2026-03-04 18:37:23.736713 INFO::Fitting model to feature number 450, cc9396a0baae4a61f1657f60aa953c58
## 2026-03-04 18:37:23.781359 INFO::Fitting model to feature number 451, X67eac2ae5c9d17d4559785aafb0c895e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.817619 INFO::Fitting model to feature number 452, X7f2b1c292a43ac11066e41099d33bf57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:23.852531 INFO::Fitting model to feature number 453, X3587846c90a1f5ac601a049f34c8ba6b
## 2026-03-04 18:37:23.902512 INFO::Fitting model to feature number 454, b5256642df3fc7ccd4b8b4a0f9bb1fa2
## 2026-03-04 18:37:23.939807 INFO::Fitting model to feature number 455, f77feba135e606ed346677c4ba78fc71
## 2026-03-04 18:37:23.974009 INFO::Fitting model to feature number 456, X8f1a91087f5fa731da48b8f2b92ad6c2
## 2026-03-04 18:37:24.007052 INFO::Fitting model to feature number 457, b73f90b9f991bfbdf47d82eaf2087ea8
## 2026-03-04 18:37:24.04 INFO::Fitting model to feature number 458, X3e55881dab6323e1f88a1163c2a44acb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.072199 INFO::Fitting model to feature number 459, X290c59a265e798e44f2e3e519155f58c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.104407 INFO::Fitting model to feature number 460, aa75bb71f1660c6e77bcd7b8c6cae618
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.135843 INFO::Fitting model to feature number 461, b65e5f621c1b50c36382b4a7caf924d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.167385 INFO::Fitting model to feature number 462, d722d9be7b36fe09e126853ef8355c6c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.19946 INFO::Fitting model to feature number 463, f04f69734d1e61c54ad34db724229aa5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.231365 INFO::Fitting model to feature number 464, X466b76e8660a84049eea28a8c28b9585
## 2026-03-04 18:37:24.264384 INFO::Fitting model to feature number 465, e7f4b79576c6742e915e9c303be5078c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.297359 INFO::Fitting model to feature number 466, ad251f7dd0a17d8e8c36a3f4a5d28dee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.328138 INFO::Fitting model to feature number 467, c4bbccfc8721e556270f6371b42deb5c
## 2026-03-04 18:37:24.358993 INFO::Fitting model to feature number 468, X287d63e09f7a925d1a73d4355be6d435
## 2026-03-04 18:37:24.389572 INFO::Fitting model to feature number 469, X59ff6f3e5e1ea3243d69c0e3a574d33a
## 2026-03-04 18:37:24.420084 INFO::Fitting model to feature number 470, X0c6c3ad654b756b0d03f35141d790a1a
## 2026-03-04 18:37:24.451228 INFO::Fitting model to feature number 471, X0c51753914969dbd012afbd6bd03c9be
## 2026-03-04 18:37:24.482539 INFO::Fitting model to feature number 472, X7203fe878b61ef50e6220aa58dbe777e
## 2026-03-04 18:37:24.513963 INFO::Fitting model to feature number 473, X06f75349b62c026bd7ee4b1e83db846f
## 2026-03-04 18:37:24.544604 INFO::Fitting model to feature number 474, X85d08c22bb731d3f1be31c53307dd020
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.575231 INFO::Fitting model to feature number 475, c7c85584ed62586869cb651f2b2f7616
## 2026-03-04 18:37:24.607662 INFO::Fitting model to feature number 476, X09d64e76596ad78536416e251990587e
## 2026-03-04 18:37:24.638694 INFO::Fitting model to feature number 477, b59ac55dd4ff09c33b997121f6a818cf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.668886 INFO::Fitting model to feature number 478, b007b4058009181e2dbff0c93d0aa6f9
## 2026-03-04 18:37:24.713105 INFO::Fitting model to feature number 479, X675e3bde986941b3b9117125b452649f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.746105 INFO::Fitting model to feature number 480, ba8218c95403215e7b6016df996ea0c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.77769 INFO::Fitting model to feature number 481, X07aa2d5d9a104b56b558e58e1cd206d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.809119 INFO::Fitting model to feature number 482, cefc53ad5dd34001d4c7115b3a40f868
## 2026-03-04 18:37:24.840171 INFO::Fitting model to feature number 483, bad88c92afdf5e15954b44dbf952aea0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.870264 INFO::Fitting model to feature number 484, X0af0efc520539e982c5d5d36ed1a2287
## 2026-03-04 18:37:24.90345 INFO::Fitting model to feature number 485, e2eb429b518244dbbde8eb1d986782f6
## 2026-03-04 18:37:24.936116 INFO::Fitting model to feature number 486, X2d6e9a6ac74bc29426134bcb1308bba5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.967971 INFO::Fitting model to feature number 487, X5793519a2f18bf75476acdace2873d83
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:24.99914 INFO::Fitting model to feature number 488, e98c5cfe89c70d3a054bf7c0bb632f99
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.0304 INFO::Fitting model to feature number 489, X601cc059b0b76f341ba3599749a94e56
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.062987 INFO::Fitting model to feature number 490, X686e2a31809281a578955a44515b0986
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.09428 INFO::Fitting model to feature number 491, X152532a23717c40d3e3c4dbdfb5b7be9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.125481 INFO::Fitting model to feature number 492, X60442eb0bb422b1646950130e527ff43
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.156244 INFO::Fitting model to feature number 493, a3be4713e16df274dbd40b0e6637bb6d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.187722 INFO::Fitting model to feature number 494, X4be568d1f92c7f820a9aa816f1c837df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.218156 INFO::Fitting model to feature number 495, X6f18e8537683b2c00b5aa2d653299c54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.248891 INFO::Fitting model to feature number 496, be66db8135bee6134a930acfb7e05c8c
## 2026-03-04 18:37:25.280148 INFO::Fitting model to feature number 497, X3a95b5b429aaf1252730d3ce59aa1157
## 2026-03-04 18:37:25.311722 INFO::Fitting model to feature number 498, fc224cd4a1ffe36b74bdec848bbba352
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.342998 INFO::Fitting model to feature number 499, X88685671fbdb6afbee06252a97db86bb
## 2026-03-04 18:37:25.373551 INFO::Fitting model to feature number 500, fa0af2b39d42b4ba5ae59ed015cf909a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.403662 INFO::Fitting model to feature number 501, X2bd74498582c27631f2164b64afb353e
## 2026-03-04 18:37:25.434143 INFO::Fitting model to feature number 502, X08a62cc1a58e42059939f913aa50e77a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.479839 INFO::Fitting model to feature number 503, X82c0fba4e73ea914e1ef1236ab224aa1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.51165 INFO::Fitting model to feature number 504, d806628b6dbda1b2faf42d2620278fe1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.542116 INFO::Fitting model to feature number 505, X828d98340f164d377d2dda60bb06f4a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.572568 INFO::Fitting model to feature number 506, X8ca98afd451a98fa6cf01f63bed22a96
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.603056 INFO::Fitting model to feature number 507, X5a53bba7d54f6d0157de0bb670839036
## 2026-03-04 18:37:25.633707 INFO::Fitting model to feature number 508, cdae58c0f1656c07f87431029fb464b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.667358 INFO::Fitting model to feature number 509, ce1c05cf3572b9f9d78199a9ff5d19c2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.698962 INFO::Fitting model to feature number 510, ce70d8fb6e8b521ca9594a0134697e3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.731425 INFO::Fitting model to feature number 511, X8f70d1241904ec8e899ec6fab95d100f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.763459 INFO::Fitting model to feature number 512, X5c4542612816ac5896eb8c902e134056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.794738 INFO::Fitting model to feature number 513, X051517485aea932a161723ecf6c596ff
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.826168 INFO::Fitting model to feature number 514, a8a38efe4ef56bb85642c5526e0b0201
## 2026-03-04 18:37:25.856475 INFO::Fitting model to feature number 515, X33d599964dd4ead082d0504ab37f2e4e
## 2026-03-04 18:37:25.887132 INFO::Fitting model to feature number 516, X0845f16afd324e69265811e5f899b19e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.917611 INFO::Fitting model to feature number 517, e4ff331359aae89dbe43026ee8d9dff3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.94785 INFO::Fitting model to feature number 518, X5a47b6a3ee0587defec53d95e2c842a3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:25.978581 INFO::Fitting model to feature number 519, X0f828935964def127d6255cdef582020
## 2026-03-04 18:37:26.008477 INFO::Fitting model to feature number 520, X1095b8196dbd38d664d439fc8a3e7634
## 2026-03-04 18:37:26.039267 INFO::Fitting model to feature number 521, X9749d4ed81f1f57122e213e1f8de77ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.070536 INFO::Fitting model to feature number 522, e738ffccde00c7cf4f544a818c42172a
## 2026-03-04 18:37:26.100969 INFO::Fitting model to feature number 523, X0f0aab6e8fef47c40ea53cebf01f27f7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.132653 INFO::Fitting model to feature number 524, X0b759bf20b1f3dc4e3c29eaca43676aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.162676 INFO::Fitting model to feature number 525, X4103c51e17a726ef97db183e51c5aabb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.192755 INFO::Fitting model to feature number 526, X29744e1859bbc3cd114636ac2f0e53dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.227714 INFO::Fitting model to feature number 527, X1d045853283b902a8ce4ccf636ee1c5d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.272993 INFO::Fitting model to feature number 528, X7256084d97c5df6b2aed179843a97804
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.306129 INFO::Fitting model to feature number 529, dfd8b8de8c4401fed140e29521b9b033
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.33696 INFO::Fitting model to feature number 530, X82c789134ae6f1ce1eb7a21113d59602
## 2026-03-04 18:37:26.367512 INFO::Fitting model to feature number 531, X3bbce5f40b435eb449659929e78a17f4
## 2026-03-04 18:37:26.397729 INFO::Fitting model to feature number 532, X46326463c37cc1a8c3236cf2ce812bf9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.429239 INFO::Fitting model to feature number 533, X55ce23b59a993eb2b7741fde91da583a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.460972 INFO::Fitting model to feature number 534, X463f1b971157c801fbf3528011b4546e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.493747 INFO::Fitting model to feature number 535, X382dd0303518af0483be8c4ea08c74c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.526396 INFO::Fitting model to feature number 536, fb1708d9521d280ffd129fc84389ab24
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.558682 INFO::Fitting model to feature number 537, X959bc5f6827ce5a4150f8c5c2608141e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.589832 INFO::Fitting model to feature number 538, X3af0184defe88fe12dd776417e66637c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.629169 INFO::Fitting model to feature number 539, X124d64c0e562f0d45b2a443f8ad1fef9
## 2026-03-04 18:37:26.659512 INFO::Fitting model to feature number 540, b9d2f025be5b0d19c8d6b8249ed48287
## 2026-03-04 18:37:26.689676 INFO::Fitting model to feature number 541, X47289f4be627b4ebd6a107b9d767b38b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.721062 INFO::Fitting model to feature number 542, X1ec3481ae3ea83074d29ca1e64c26cba
## 2026-03-04 18:37:26.751871 INFO::Fitting model to feature number 543, X5c3085994e1bbbf662babdfaf0ccf8bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.782764 INFO::Fitting model to feature number 544, X084347470f75d9e030c777fab7bf5930
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.820927 INFO::Fitting model to feature number 545, X34a5cc3ff0c374acc93606948a77b109
## 2026-03-04 18:37:26.851295 INFO::Fitting model to feature number 546, c57234d6515915a9ccd83265178c38a0
## 2026-03-04 18:37:26.881867 INFO::Fitting model to feature number 547, X1630c97bc99db59912c92545f184a2c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.913365 INFO::Fitting model to feature number 548, aa7bae079433f6ed40a3d3016cb2b4e1
## 2026-03-04 18:37:26.952092 INFO::Fitting model to feature number 549, X0b9987df25b2cbf3e32983f0bb343d22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:26.983949 INFO::Fitting model to feature number 550, X4fbf84a3fdfe517ab4a03517aa50188d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.015173 INFO::Fitting model to feature number 551, cc628ddb3cb5ec7eca660a77860a0dfa
## 2026-03-04 18:37:27.060243 INFO::Fitting model to feature number 552, a5d3cef1b84bffa024ae0259ebbe7819
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.092937 INFO::Fitting model to feature number 553, X001d2e9ba87c0b0118ae91de96ff290f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.125416 INFO::Fitting model to feature number 554, X1f36cbd271a0a85fe9d8415377e611c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.158169 INFO::Fitting model to feature number 555, X421c87630f8796da746888ec4aa55f4f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.1892 INFO::Fitting model to feature number 556, X6eff287d92b730a0169e0b2142f31b8d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.221324 INFO::Fitting model to feature number 557, X039ec62846d8972fa2fa02f2bd2056f7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.253098 INFO::Fitting model to feature number 558, X33cec55e74e9b122f1e03586050d529e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.284032 INFO::Fitting model to feature number 559, a8982ccbc26415e739f9c1dea2471f4f
## 2026-03-04 18:37:27.313448 INFO::Fitting model to feature number 560, ae4d72c28426dd93660524df503c9d97
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.345066 INFO::Fitting model to feature number 561, X173f0a0e35c5a252dcaddc1c05e40cc6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.375815 INFO::Fitting model to feature number 562, c88ac9ec2ae50bb75443a61e09401ebd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.407565 INFO::Fitting model to feature number 563, X9c10353ae637b8f0c7fb09b60c0d18d3
## 2026-03-04 18:37:27.438375 INFO::Fitting model to feature number 564, ba4f8c14f1a6f625c90f05848475b43a
## 2026-03-04 18:37:27.469219 INFO::Fitting model to feature number 565, X42a703e9af69e2172af67521bbc356cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.499382 INFO::Fitting model to feature number 566, da4cb8754d5d4a0317d887304c60b9f8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.529321 INFO::Fitting model to feature number 567, a2f33842d6158bd3ff1162b1b575109c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.560239 INFO::Fitting model to feature number 568, X18e76539a844edacd926d64850b8aa36
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.592481 INFO::Fitting model to feature number 569, X218c575bc466bd7f3de4855d1cb785c9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.622716 INFO::Fitting model to feature number 570, X527c201e40ff4e8b5562ed88b67ead25
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.653602 INFO::Fitting model to feature number 571, b1f3d0e531c2a6e1fc1ca63b5c5506e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.684191 INFO::Fitting model to feature number 572, X252d977d5b531a6a265d557f412e2d44
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.715075 INFO::Fitting model to feature number 573, X3015c2163a624808acd0a7c62a3bf0c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.746308 INFO::Fitting model to feature number 574, X3e1bcd69c044a695703a1985de3bd920
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.779106 INFO::Fitting model to feature number 575, X5d42ffdab90cffce3e9c2cfb952cc2fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.810625 INFO::Fitting model to feature number 576, dbb13a4636bef149f210c12027d77f9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.891876 INFO::Fitting model to feature number 577, X5ad1570deadbd75f87675e54addf230d
## 2026-03-04 18:37:27.923612 INFO::Fitting model to feature number 578, X13562e00ec4db4d18cd5dd1bd73d2ef3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.954462 INFO::Fitting model to feature number 579, X69419db7063b243b238cc5bae949369a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:27.98494 INFO::Fitting model to feature number 580, X91f7efb223b32216b23bb9a02e9ec570
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.015585 INFO::Fitting model to feature number 581, dfc46dd0818f492b120349ace278d6d8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.04841 INFO::Fitting model to feature number 582, fcecd1c29c2e37c8c3de46a1926cf892
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.080459 INFO::Fitting model to feature number 583, X4a396a487ee50aa361519dc2a2369a5f
## 2026-03-04 18:37:28.111204 INFO::Fitting model to feature number 584, eee8f061f4328896cbad20ccb6bf41f9
## 2026-03-04 18:37:28.141524 INFO::Fitting model to feature number 585, X5cedacc9704fb36ee1087cdc9bdf80ed
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.172155 INFO::Fitting model to feature number 586, b83acb9199413447ed91faec843cd908
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.203537 INFO::Fitting model to feature number 587, X6dc0a98dc381716943fa9393aa9cb416
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.233956 INFO::Fitting model to feature number 588, aed586669a5eb4831e60da011eac8ca9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.264363 INFO::Fitting model to feature number 589, d2b971f10e7c69e1d4c89c63cfbe9b94
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.295276 INFO::Fitting model to feature number 590, X00c145ccfb95cec516f4fa326d896af4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.327743 INFO::Fitting model to feature number 591, X06b4246e0976c80b3b656da49f9599c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.359462 INFO::Fitting model to feature number 592, X320902985335212dfccb275337dccacc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.391133 INFO::Fitting model to feature number 593, X3765a7f3a9962418ea2f3e764d644898
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.423306 INFO::Fitting model to feature number 594, b1d5ee1846c492211f4fed9aeb829a84
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.455527 INFO::Fitting model to feature number 595, daa12bd954485797a7fd1a0f2b9024eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.487214 INFO::Fitting model to feature number 596, dfb91dda0f7b903c88c7082678e4ce4d
## 2026-03-04 18:37:28.517999 INFO::Fitting model to feature number 597, e9d0e82558715ca0ca438f41fbbe4ef6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.549114 INFO::Fitting model to feature number 598, d35bd24b097c811e40172a1fc6d15a5e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.583561 INFO::Fitting model to feature number 599, e737489b71d52736ee1a0e2a50549d4c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.615255 INFO::Fitting model to feature number 600, cd5308d7409829a5057c4029560fac2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.659106 INFO::Fitting model to feature number 601, fe72c667e0eb30b792b6b59fb292c54c
## 2026-03-04 18:37:28.691515 INFO::Fitting model to feature number 602, X1d4285e838c09afb01b627dae637beb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.723653 INFO::Fitting model to feature number 603, X612fbc07b28756a2e743392fde08e155
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.755657 INFO::Fitting model to feature number 604, X76a776e3b51070fc91491000f4030ae1
## 2026-03-04 18:37:28.786879 INFO::Fitting model to feature number 605, X2e246d891a50a35543494663b056b11b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.819797 INFO::Fitting model to feature number 606, X07f2b771a3355046b9ac7e90670a48b6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.852103 INFO::Fitting model to feature number 607, X6cecdfc4f05450d653299bdbdc9f26db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.884014 INFO::Fitting model to feature number 608, c7f6f5b0e601f6374c521ec454e23a90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.91462 INFO::Fitting model to feature number 609, ddd96397751cb6c487e71dca3fb8eedb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.946123 INFO::Fitting model to feature number 610, ee391c5cbd032373f7fa130bbdfa78b2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:28.976943 INFO::Fitting model to feature number 611, X331d529d33e29b0d43456b0ea53b88fe
## 2026-03-04 18:37:29.007781 INFO::Fitting model to feature number 612, b483981d9b053abefc0935cb20c58b66
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.03846 INFO::Fitting model to feature number 613, e82710e82def354d4332ee5b055234c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.0703 INFO::Fitting model to feature number 614, bacc8fba734fc18196d62879e73b30c2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.100846 INFO::Fitting model to feature number 615, X00e9012a608d346ef966588502aff2df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.131174 INFO::Fitting model to feature number 616, X171c685c03b82e7f5fbe12bd0ac2e720
## 2026-03-04 18:37:29.162105 INFO::Fitting model to feature number 617, X387aece6e4ad6511151d09055a8d2d5f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.193091 INFO::Fitting model to feature number 618, X510b61e30ef7cba2b46623f71a09d2bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.223204 INFO::Fitting model to feature number 619, X925c84e4221675297ab469c3758bd13e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.253203 INFO::Fitting model to feature number 620, a4be22de4b5390e64bfe2b6011523f3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.283283 INFO::Fitting model to feature number 621, X1764219dbd7ca6405b38a4f1122b5270
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.314062 INFO::Fitting model to feature number 622, cec7bcba266e98babcbadffedf294db7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.345468 INFO::Fitting model to feature number 623, e7929ae3c3f9c7c26009da7f71875850
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.375954 INFO::Fitting model to feature number 624, f5cd996896ec42be1a2d4513c5d663af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.406145 INFO::Fitting model to feature number 625, X1e0bd8313bb83eed0c2ca563b106ea19
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.450554 INFO::Fitting model to feature number 626, X12a7a583a9d85fc8ca1bad3a26192c2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.482634 INFO::Fitting model to feature number 627, X1ecbf29ff5cb9d6bf1eb427bcb0290e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.513179 INFO::Fitting model to feature number 628, d5b669bb5d3bff1c87f4c50299e68a13
## 2026-03-04 18:37:29.543686 INFO::Fitting model to feature number 629, de809f907097c58da02f0c46a3f761d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.574149 INFO::Fitting model to feature number 630, X1208c5d5a63919dfdce05b8942fca118
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.605153 INFO::Fitting model to feature number 631, X263e8a3f16a6c22a7991e8ee33351391
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.635963 INFO::Fitting model to feature number 632, X68f8d7429d5fe7dc889f92705924aa70
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.666838 INFO::Fitting model to feature number 633, fda375197ec72a134ce3c008754dffae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.697718 INFO::Fitting model to feature number 634, X3f25c26eb5e2ad00677e1c3341010d70
## 2026-03-04 18:37:29.727999 INFO::Fitting model to feature number 635, f96dbf5f5deb0558ba732d9e264601a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.75927 INFO::Fitting model to feature number 636, daaf8b158006c74d73fee01ca65fb5e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.78982 INFO::Fitting model to feature number 637, X04314b4b349b7be1b1adca749349e9d5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.820096 INFO::Fitting model to feature number 638, X24cb510a92ff13206c2c50e208ea26da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.850782 INFO::Fitting model to feature number 639, X0e372f45908a23ea01a2cfc4db466bd5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.881769 INFO::Fitting model to feature number 640, f5375d7a11214ac9f11f523325b2348b
## 2026-03-04 18:37:29.912834 INFO::Fitting model to feature number 641, X82d1fcf8762f0303c6a4c50d61e5ff8d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.943078 INFO::Fitting model to feature number 642, e903cac37add63e8c94650b15bfefd86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:29.973119 INFO::Fitting model to feature number 643, X6fe94b71e60dfd02b295e38a5f3e8da2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.003094 INFO::Fitting model to feature number 644, X17d13995f3dc073bc150a3ce1f936667
## 2026-03-04 18:37:30.033324 INFO::Fitting model to feature number 645, X9c3ba0cc2b8bcadd5776b16bf7947e5a
## 2026-03-04 18:37:30.069332 INFO::Fitting model to feature number 646, X73a596d4c84e2cc607f447b8d86e410d
## 2026-03-04 18:37:30.099475 INFO::Fitting model to feature number 647, X6936efdd0d4dc2c78a52637bbca6698e
## 2026-03-04 18:37:30.130019 INFO::Fitting model to feature number 648, X1accfe46cbb5fa4b43f42f12a4b09a46
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.168278 INFO::Fitting model to feature number 649, X372ae668689d07c1776f9e39c19b7469
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.211259 INFO::Fitting model to feature number 650, X3586e5ed26b0fe9d2b88c41818b48636
## 2026-03-04 18:37:30.243221 INFO::Fitting model to feature number 651, X4ad6f2f290bec0999affbdae99e9d315
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.276096 INFO::Fitting model to feature number 652, X43722f6ec1f24a34bdef18e320dd41a4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.308468 INFO::Fitting model to feature number 653, X6205f70998dd9ef2c6c8769ed201ca0e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.34755 INFO::Fitting model to feature number 654, de7b4217cc9e4c4ccabf081ab8157af8
## 2026-03-04 18:37:30.385664 INFO::Fitting model to feature number 655, X8d25b30e61a1542d118f0a1d96420b8b
## 2026-03-04 18:37:30.416582 INFO::Fitting model to feature number 656, X16508d3f8de41b678f64b67a13087618
## 2026-03-04 18:37:30.450533 INFO::Fitting model to feature number 657, X7bf58f84d1c8c8d1d013725ba2a6d42c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.489638 INFO::Fitting model to feature number 658, c567686addbcf1ec7c0326a8cb050ee3
## 2026-03-04 18:37:30.520792 INFO::Fitting model to feature number 659, X8ed345e1519449d98c9ae839e788ed75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.55549 INFO::Fitting model to feature number 660, feea58b2edca226cce18140b7a96f212
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.588362 INFO::Fitting model to feature number 661, ffaf1f6f3aa2f47240308a9757599ff3
## 2026-03-04 18:37:30.619852 INFO::Fitting model to feature number 662, ee7a112d25ef18a5ac190015e789ab7c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.65217 INFO::Fitting model to feature number 663, X8b79ecbf1789202e32c2134d323f3fe5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.683636 INFO::Fitting model to feature number 664, c56a8b428e298d443b0046bd7da3b033
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.715026 INFO::Fitting model to feature number 665, e08bffe7c47dc649f45e6605eaf77c00
## 2026-03-04 18:37:30.74702 INFO::Fitting model to feature number 666, X390bd5af0dfc36880157790449a6934a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.784595 INFO::Fitting model to feature number 667, abe04980ffc876dd5ddefee10c71e8da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.820113 INFO::Fitting model to feature number 668, X03b7c9224944f7a08e43db638d43e1e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.852713 INFO::Fitting model to feature number 669, X6713a01fce8ed800fdea27cfef912766
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.884051 INFO::Fitting model to feature number 670, X8a7b1595d4369644e07e85b86c70e3db
## 2026-03-04 18:37:30.915001 INFO::Fitting model to feature number 671, e4a9ae0bba905d4ac7db962abab2e295
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.946016 INFO::Fitting model to feature number 672, X4855f68229e7d81a2e9da0111a5ee938
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:30.978538 INFO::Fitting model to feature number 673, X9700a2382c70667f3daa3cb2a55d0a7a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.010146 INFO::Fitting model to feature number 674, X24226920f12955b61a66789e6ad6f215
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.054887 INFO::Fitting model to feature number 675, X4c0ab4869eec1c6667ca4fc83f477e39
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.086564 INFO::Fitting model to feature number 676, X02f8bdb4409411a2601dc87c551942a1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.116614 INFO::Fitting model to feature number 677, edce391df3040a66ce3f198029a9873c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.147469 INFO::Fitting model to feature number 678, X3b3e296ef6b3c79cdc410d0d662a6237
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.179 INFO::Fitting model to feature number 679, eaea5405910f416f1be4d4d0b71ac90b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.210261 INFO::Fitting model to feature number 680, X0b78143d11a6ca0fa5137a6b875c3a51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.240416 INFO::Fitting model to feature number 681, a832ff8e0373f88c777d1e717dc8cf21
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.270745 INFO::Fitting model to feature number 682, b9f2ffae17e0a92098de47a3b3ea81c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.301226 INFO::Fitting model to feature number 683, bad508016651bf74cf3cb269911ec12b
## 2026-03-04 18:37:31.33433 INFO::Fitting model to feature number 684, X213c0a6d6aaa77c7237d1c27dd8f3a78
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.365489 INFO::Fitting model to feature number 685, X43c87171c86c75f5fd3ff3da1f251b0f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.395391 INFO::Fitting model to feature number 686, X49a98f8f4510c8c171fc2ba1461073e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.426124 INFO::Fitting model to feature number 687, X72388b9f07ab344f60eb50a9fb4ff6b0
## 2026-03-04 18:37:31.457054 INFO::Fitting model to feature number 688, X50577f758ed98fa0971590d1f3818602
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.488293 INFO::Fitting model to feature number 689, c74a32f53cfb587c14d6d0326b225a5e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.519931 INFO::Fitting model to feature number 690, e7eeb1fe634e3182ff6d596bcb0674c1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.550891 INFO::Fitting model to feature number 691, X4c1637f4177c0b2d8f35befab7c70508
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.582531 INFO::Fitting model to feature number 692, X13420d21573c308fc7d77a3f77b01f6f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.616695 INFO::Fitting model to feature number 693, X9a0c6caa7312d271153e400a3e15903a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.647409 INFO::Fitting model to feature number 694, X9e8b037cd3a3c0d0c1abb98279a78e5b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.679152 INFO::Fitting model to feature number 695, b311631596b18e1f7ae54d3b969844c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.714287 INFO::Fitting model to feature number 696, ddc7f083639aa23de81829ed19d1faa2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.744477 INFO::Fitting model to feature number 697, X92fdc7ec8d89c654aeb51aa18a540bf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.774467 INFO::Fitting model to feature number 698, X537d57d1ee4850a807f6c6bb9f477897
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.815936 INFO::Fitting model to feature number 699, edf3da53246c455f9b500ac2db192e1e
## 2026-03-04 18:37:31.848598 INFO::Fitting model to feature number 700, bae5efc0e9e141ca3aa86b1f58ee238e
## 2026-03-04 18:37:31.880925 INFO::Fitting model to feature number 701, d8ba154ce25da930b156278285fa1afd
## 2026-03-04 18:37:31.911313 INFO::Fitting model to feature number 702, X0c101f614e28abcff7de67c17ca51a55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.942127 INFO::Fitting model to feature number 703, X7ea0a5c2e540405f53d0c30fa63c82a9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:31.972472 INFO::Fitting model to feature number 704, ddce072baa0bb0342e2fd9a8cd7cf489
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.002698 INFO::Fitting model to feature number 705, ea88fdd7f58f1c24015d406597798885
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.03328 INFO::Fitting model to feature number 706, d294afa38c6e4bc2e50c1bc296583d22
## 2026-03-04 18:37:32.065189 INFO::Fitting model to feature number 707, X1d7bbdef9955168131018ddbaf14993f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.097154 INFO::Fitting model to feature number 708, X1de141b5e2a83b1fbef3d71352ba38fc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.127649 INFO::Fitting model to feature number 709, X3797616e5f13944c402ba8e5c97de710
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.15895 INFO::Fitting model to feature number 710, X82bc96ed61841410472a980eadcc456a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.18986 INFO::Fitting model to feature number 711, ac9d52c7b33c42e07fc95741ebdccfb4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.2208 INFO::Fitting model to feature number 712, ba1ae9ef946c2c8619f46f0a85e39876
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.250834 INFO::Fitting model to feature number 713, c61c0fc771ac44540c315cf1e48285b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.281143 INFO::Fitting model to feature number 714, cb97cafd5cd9403ff6238b886b39a5cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.311628 INFO::Fitting model to feature number 715, cf5e11d1b9fdca94498f62c9f6a9163e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.343074 INFO::Fitting model to feature number 716, X36135ac55ed4d0ff6e31d6ff7aa4499e
## 2026-03-04 18:37:32.372986 INFO::Fitting model to feature number 717, X8371b60e37ffedd9c1a9a37fef3e3286
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.402898 INFO::Fitting model to feature number 718, a1efd505557860248435d7bea1bf3913
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.433584 INFO::Fitting model to feature number 719, d3c7a110bef8e513fdc270909827b2ea
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.466166 INFO::Fitting model to feature number 720, X3a5942b1a3e4aab99de89f5bac550156
## 2026-03-04 18:37:32.497935 INFO::Fitting model to feature number 721, X43e163b10a20212c4a092bdafae4f956
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.528658 INFO::Fitting model to feature number 722, X5b3835118edb11aeb9d4bab295d72cb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.559895 INFO::Fitting model to feature number 723, X6617d00c29734b4e9411fa53e40f16c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.604925 INFO::Fitting model to feature number 724, dbc6c4998c0721a95236e008182051c4
## 2026-03-04 18:37:32.637008 INFO::Fitting model to feature number 725, X6c92c5cc06059722f26df4dbf9eda96e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.668184 INFO::Fitting model to feature number 726, X82c546ae87e738cdf074853e15bed022
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.69881 INFO::Fitting model to feature number 727, b20e3d817c8a6fc306663500024b2479
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.729918 INFO::Fitting model to feature number 728, c5530dfee5b10d94b0fe01f9dbcb002f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.760754 INFO::Fitting model to feature number 729, X48a8fe019dc6cc1cfb047a123527b5c0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.792924 INFO::Fitting model to feature number 730, X77d5f38002bb3780f14d3d9880aa9008
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.824207 INFO::Fitting model to feature number 731, X8a7640fd3fdee3d1d7dc786a9ddd6673
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.855717 INFO::Fitting model to feature number 732, a66208cbb348d255e7357fae6c83e793
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.888455 INFO::Fitting model to feature number 733, b153c75d7dae701723a3e55a0e7610a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.920794 INFO::Fitting model to feature number 734, X1230aedcda2b9095c4b718d403c11b65
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.95934 INFO::Fitting model to feature number 735, X1d785c4e2c7f07957b49d819a3c473d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:32.990308 INFO::Fitting model to feature number 736, X2630cd266dd2f91e8e73919d7ec61b3c
## 2026-03-04 18:37:33.021036 INFO::Fitting model to feature number 737, X3adbb59e7c520b6a4948e315e22aa48b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.0534 INFO::Fitting model to feature number 738, X4547cdaa4d234e3ca4e24d18dbb1a7cb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.091012 INFO::Fitting model to feature number 739, b41174e1a26fcc77c1d1fecf6a39e3eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.122035 INFO::Fitting model to feature number 740, bc31699d255fca6a180a44290780aa5d
## 2026-03-04 18:37:33.152312 INFO::Fitting model to feature number 741, X3b918737053456dd3b8f2c656205dc3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.183124 INFO::Fitting model to feature number 742, X3e9a95e009b5e5728ca7b5861a9f6886
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.213967 INFO::Fitting model to feature number 743, X723bb2a4f2424433622e58cb4c184f68
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.244919 INFO::Fitting model to feature number 744, X832b2920c7d4c0c22e2d9b24f716e64a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.275667 INFO::Fitting model to feature number 745, X936c188690c04d1eec4e44cd1f966362
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.306841 INFO::Fitting model to feature number 746, bc7431fe3e8d871068b783cb153d98ea
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.338418 INFO::Fitting model to feature number 747, de1af36e20ecfa838185c13b9e66c677
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.379394 INFO::Fitting model to feature number 748, df6b94f16e20f444ca41aa31497124eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.414014 INFO::Fitting model to feature number 749, X1a66cd529f8773b8626aac7adc240a1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.446863 INFO::Fitting model to feature number 750, X2fd2e9098fb21a08c0855b9f03580f79
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.478828 INFO::Fitting model to feature number 751, X91d9ae565c05f931c2ceb57abeffbce8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.509631 INFO::Fitting model to feature number 752, d2dcbb92bf9cc3ab8d247629e9e07526
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.540705 INFO::Fitting model to feature number 753, X075dc0a8cf360e04255c741bde0699d0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.574004 INFO::Fitting model to feature number 754, X25ee65871f64cda654d5f4db4e7aea31
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.606026 INFO::Fitting model to feature number 755, X8bb457c94ddbae8fd75111d1598017aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.637283 INFO::Fitting model to feature number 756, d89d081df8bbf354930bda36c31f8d00
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.668396 INFO::Fitting model to feature number 757, e91ff13a5dab0760ae86ed2b229a578d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.699453 INFO::Fitting model to feature number 758, X0a62aa5adaaa05fc8876ed148f03fee9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.730644 INFO::Fitting model to feature number 759, X5f54914fe6d2951bfab121606df7ab31
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.761568 INFO::Fitting model to feature number 760, X8fb569330957254554e77948ce449cac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.800081 INFO::Fitting model to feature number 761, a85443cff26b1389b3ecc414e37c7da6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.831716 INFO::Fitting model to feature number 762, b4273d59573cc45dc778e776f176b9aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.864243 INFO::Fitting model to feature number 763, de57ef35f117467f9b193d96e53ae163
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:33.895104 INFO::Fitting model to feature number 764, X2886d1db4c3b6787f79e50bef8f69c55
## 2026-03-04 18:37:33.925403 INFO::Fitting model to feature number 765, X743535eee04707e3424f814948e07b47
## 2026-03-04 18:37:33.95635 INFO::Fitting model to feature number 766, X70a5b7339c7527388c360d247e8df171
## 2026-03-04 18:37:33.986936 INFO::Fitting model to feature number 767, X4879cd5a48e851e25a3ca727ce334b32
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.017357 INFO::Fitting model to feature number 768, X1be31d76bb8cf9406332b7282354d758
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.048523 INFO::Fitting model to feature number 769, X163e8909e8fcf7a3b229cdbb8f26275f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.080196 INFO::Fitting model to feature number 770, cc81bf9bf985caa468eaa6b7fdfd1bf8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.11076 INFO::Fitting model to feature number 771, X04ceef6f70d00472bdd1af80c2f30983
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.148495 INFO::Fitting model to feature number 772, ae1394fa5a7795914411b0d3ea6546f1
## 2026-03-04 18:37:34.192364 INFO::Fitting model to feature number 773, da9409cba5c8301a2b42580176f14265
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.232204 INFO::Fitting model to feature number 774, X22d3ed64f87a280651657ccc0b581145
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.26588 INFO::Fitting model to feature number 775, X4b17e66aa631354e9f3f33314d089589
## 2026-03-04 18:37:34.302983 INFO::Fitting model to feature number 776, X83d8438154115532d2f599781802ff7a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.333248 INFO::Fitting model to feature number 777, X36ebd620b98f5ae6205c417f489e5207
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.365696 INFO::Fitting model to feature number 778, X5b550a3ccef0538fc08e87246a6a35b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.400097 INFO::Fitting model to feature number 779, X168b9624445c90028b89ca20afea6e4f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.438927 INFO::Fitting model to feature number 780, f436778965fbeae0dea00fc287ba37da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.469814 INFO::Fitting model to feature number 781, f54c94989621c5be35e72ff3e29893d2
## 2026-03-04 18:37:34.500605 INFO::Fitting model to feature number 782, X323474086152e9c4666a43f5eb515d05
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.532421 INFO::Fitting model to feature number 783, X96d483b6964c2f2ee2b237efb9f01baa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.563021 INFO::Fitting model to feature number 784, X3e12212b21a98a5799ebdf835788a4a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.593028 INFO::Fitting model to feature number 785, bef4a1751e1ee0ba5ba0366c85965eca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.623793 INFO::Fitting model to feature number 786, X6c69bfaaa63d25c2bf7e055ee84e3a90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.655015 INFO::Fitting model to feature number 787, X39f844378114cbe9256d291fff454c1f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.685963 INFO::Fitting model to feature number 788, X6455b26d66a9a71fe4411b6ba3b3ec4b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.716074 INFO::Fitting model to feature number 789, X5bdd8b2852cfb534356a12a5e39d96bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.746184 INFO::Fitting model to feature number 790, X41d4dbb4027e632037882ead150e113f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.776422 INFO::Fitting model to feature number 791, X25db1ff4778c4f47f626c8929563e2c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.80634 INFO::Fitting model to feature number 792, X79c05c36d834303543f0231cb21a2eac
## 2026-03-04 18:37:34.835718 INFO::Fitting model to feature number 793, fe54fff680b4cf948d8d943c2cc4f091
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.86648 INFO::Fitting model to feature number 794, X37d384450b5ee65ed1cc2bab56af5e1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.896867 INFO::Fitting model to feature number 795, X9b399f779ade5e5b6500eec0de6d13b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.927839 INFO::Fitting model to feature number 796, X6c04b85ba8beb48dfb0adfb7db7fbe7b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:34.968056 INFO::Fitting model to feature number 797, X8e0dc3caf705d8984bea5e9c501296d4
## 2026-03-04 18:37:35.001854 INFO::Fitting model to feature number 798, X18d2efc0bdd68758bba83ab330bd15bd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.041254 INFO::Fitting model to feature number 799, X0f3a73a457224cf1b94c2830fb8c6e86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.073047 INFO::Fitting model to feature number 800, cfcdf20e1255691d28fe0c56ce2016dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.104603 INFO::Fitting model to feature number 801, X21057ab2a1913644b88384f538082f65
## 2026-03-04 18:37:35.135548 INFO::Fitting model to feature number 802, X4889a0edb0de0f30c5f53c9aa5417068
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.166912 INFO::Fitting model to feature number 803, aa09d6b56ace46671bf60d93a5aab529
## 2026-03-04 18:37:35.198644 INFO::Fitting model to feature number 804, X0af82abc9d7b12a57264ed2f7d7c635c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.228872 INFO::Fitting model to feature number 805, X89c6d2a235c32ea0290509ae65221057
## 2026-03-04 18:37:35.25929 INFO::Fitting model to feature number 806, X34f8f51f3aa37533d2fc7d4cba1aefb9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.289687 INFO::Fitting model to feature number 807, X738286c5d38996bf4df0acd3eabbcc70
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.320939 INFO::Fitting model to feature number 808, ecc5e663311d028a84832f7d6db3d64b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.352452 INFO::Fitting model to feature number 809, c68dad1066d5ef8221fcc1093790b7ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.387115 INFO::Fitting model to feature number 810, X1e07bf1cff0e9e56f62458ab13710672
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.418061 INFO::Fitting model to feature number 811, X5c08d93436793f5fc51b871db4d85870
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.450438 INFO::Fitting model to feature number 812, X8070db26c5b69b8f7a174f40e71b085f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.480717 INFO::Fitting model to feature number 813, eef930d6e152847be082212f21d6fc69
## 2026-03-04 18:37:35.510967 INFO::Fitting model to feature number 814, X8dfbe52cb8feca704bc7d52c371d052d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.541422 INFO::Fitting model to feature number 815, X2958d05407011ced17c694c365dee47c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.571078 INFO::Fitting model to feature number 816, X5ad50b5f8b90a5fd0a409661a772acb2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.601386 INFO::Fitting model to feature number 817, X6de335ece7f713b9a8c71b18aaea1d4c
## 2026-03-04 18:37:35.631671 INFO::Fitting model to feature number 818, a0347f3b93155f49f127482da58b5d25
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.662725 INFO::Fitting model to feature number 819, X1d1fee5b3f5f2f03a6d43b68c228c61e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.693344 INFO::Fitting model to feature number 820, X3e232c33e572aeba36d4ac20b90b2822
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.724312 INFO::Fitting model to feature number 821, X59687c67da77385562ad90818c85468a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.767677 INFO::Fitting model to feature number 822, d3f8153113111d55a5a906a131a161bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.801401 INFO::Fitting model to feature number 823, X45819013a025366b17b161a0eacfe02c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.833045 INFO::Fitting model to feature number 824, X4c52e18b01f156b0b4030db975dd2c90
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.86441 INFO::Fitting model to feature number 825, b432978e0c1aea4c678fef2f27d03ff1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.895802 INFO::Fitting model to feature number 826, eb0c2a0ec0c24342ab7639a869020215
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.927729 INFO::Fitting model to feature number 827, X1633280dca1ce6b99f323fb357a788f4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.959311 INFO::Fitting model to feature number 828, X19900408557ac9f0e79bf455cd982c2b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:35.991003 INFO::Fitting model to feature number 829, X0b2e8985a71922126f4f44935b75e37f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.022419 INFO::Fitting model to feature number 830, X5598141e7b6d1966331ae03cbd17c8ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.054391 INFO::Fitting model to feature number 831, X82e9bb8861cef7612bb7ae065798491e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.0862 INFO::Fitting model to feature number 832, X9c9fdbf346a77e3a795b6d9a7480b9bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.116441 INFO::Fitting model to feature number 833, e0c60e831deae33ffe603459de3fceff
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.147036 INFO::Fitting model to feature number 834, e92ba24ee64bcf3f4cb5f9a3682e4413
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.177986 INFO::Fitting model to feature number 835, X4d9a5071591ffc708722c047fe488fdf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.20909 INFO::Fitting model to feature number 836, X63c0016bfe6d090f37f597cd4df623a9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.239857 INFO::Fitting model to feature number 837, X977936e75d88dab18ff37a8e776040fe
## 2026-03-04 18:37:36.273487 INFO::Fitting model to feature number 838, b0c41adf5d284fbabdce3574ce198155
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.307062 INFO::Fitting model to feature number 839, be27327880cec05e98bf577491a9f35c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.340462 INFO::Fitting model to feature number 840, ee5a2aed5e233e505a43fd5e5d77f7b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.373762 INFO::Fitting model to feature number 841, X144a95b2a5248802ffb64716fbcee74d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.404393 INFO::Fitting model to feature number 842, X23d44c8c40fb72845b2f99c3f746be86
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.435883 INFO::Fitting model to feature number 843, X245af114f65c6751ab4082702eb5717b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.468058 INFO::Fitting model to feature number 844, X52b05bd5be5f8a56cdac5e80ac038690
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.498989 INFO::Fitting model to feature number 845, X9fea56b0e2353c17ac863a7431b495fa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.542517 INFO::Fitting model to feature number 846, fd5bf8c9e268f6c5f48131c0615125fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.576922 INFO::Fitting model to feature number 847, X162a0a70867bfcb444146c5462dd083f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.609466 INFO::Fitting model to feature number 848, X2d8c6fa19ed17f298d900b8b94ca3577
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.640597 INFO::Fitting model to feature number 849, X9708c6393469d007561fa62858323ef3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.670906 INFO::Fitting model to feature number 850, X95312c6f0ad3ebf16fc431a571042056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.702966 INFO::Fitting model to feature number 851, X9cf77d8c73a399bf0ce345130e8448a6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.73297 INFO::Fitting model to feature number 852, e8f211866f017fe1a3cb035bbbf0f487
## 2026-03-04 18:37:36.765391 INFO::Fitting model to feature number 853, X28c021635cbdb3e420e03898e6736980
## 2026-03-04 18:37:36.802241 INFO::Fitting model to feature number 854, X06d2ce585d904c386410802b9ca7f6dd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.833003 INFO::Fitting model to feature number 855, X96c5ae4d64ae76c5f2e1dedfb7ee89b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.863948 INFO::Fitting model to feature number 856, X263f6b9f13db9a695ba7215622ccf3e7
## 2026-03-04 18:37:36.895594 INFO::Fitting model to feature number 857, e6c6d6c4f97b71d8b0747ea36185da7a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.925561 INFO::Fitting model to feature number 858, X8bb590cfc8fc74aad9d4f8a958a850f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.956269 INFO::Fitting model to feature number 859, X5d7946238ac7676f2e6c32ce5554a326
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:36.986211 INFO::Fitting model to feature number 860, X4dd5e98c725ec5f2fcc4a0c284b585e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.016987 INFO::Fitting model to feature number 861, X4aeeaaabe391e0fbeba72e87f6468d6d
## 2026-03-04 18:37:37.048554 INFO::Fitting model to feature number 862, X300433f8069de9c7d93f326c2bcb5e71
## 2026-03-04 18:37:37.079014 INFO::Fitting model to feature number 863, X036a46c77db16f0cab39f25993f32f5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.109147 INFO::Fitting model to feature number 864, X3bd6c8564a3f23aac3247f77acfd8adb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.139642 INFO::Fitting model to feature number 865, f8d3badaf1aa1b3e2e544de4310e176f
## 2026-03-04 18:37:37.172591 INFO::Fitting model to feature number 866, X28d6d19e890138b314ae244a6511662e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.203513 INFO::Fitting model to feature number 867, e54123bfdac1a2a0b1468c72d3a72056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.235737 INFO::Fitting model to feature number 868, X10507e5ecf030ee88fa122218bfa039e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.267007 INFO::Fitting model to feature number 869, X6fde2683f8248a303b9d6880e4a4e604
## 2026-03-04 18:37:37.297243 INFO::Fitting model to feature number 870, ba90f0b9771f6de4992c42c2004c1f75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.34772 INFO::Fitting model to feature number 871, a107e80c66b63e13999c5ca78894559c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.38226 INFO::Fitting model to feature number 872, e1f2d98fc9b94c8e80eb3124ab697199
## 2026-03-04 18:37:37.413011 INFO::Fitting model to feature number 873, X04d946612707d3fc0a92f32dcc8d1c30
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.447866 INFO::Fitting model to feature number 874, X33dab83ed982b34105927f635d5ce4d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.479011 INFO::Fitting model to feature number 875, X7feea0d71ff9eb9f5ea3d95e10d7b7b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.509646 INFO::Fitting model to feature number 876, X9eb591d492e323d974da25b0f092f73c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.540372 INFO::Fitting model to feature number 877, X0a8faf4be9b0ab6ca854de7a5bcad02e
## 2026-03-04 18:37:37.570469 INFO::Fitting model to feature number 878, X417159e7c56b61db9020926381287e6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.600716 INFO::Fitting model to feature number 879, X62871a4f2e689b3d2a8e24b4b6e90704
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.630989 INFO::Fitting model to feature number 880, X98345182f8f896dcba7b94aa00dd7eb7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.661889 INFO::Fitting model to feature number 881, X117ce9780a6794585b2b195cd40c66e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.693405 INFO::Fitting model to feature number 882, X1bf12663366670a5853f8ef1f662349b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.723407 INFO::Fitting model to feature number 883, e445d1b02329b8d17b6c33e47b2a1d35
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.75315 INFO::Fitting model to feature number 884, f2a3db16c28d4e5de91808d1c325f84e
## 2026-03-04 18:37:37.782408 INFO::Fitting model to feature number 885, X04646829704cf62ef942ccf8674fa3ce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.812869 INFO::Fitting model to feature number 886, X262b03eb315b00d94056179bb2f4ee53
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.844064 INFO::Fitting model to feature number 887, X446609f70ba4143bd817f44c210a0a4c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.874499 INFO::Fitting model to feature number 888, X794e9f08539a1c21267d2a49997fdddd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.904863 INFO::Fitting model to feature number 889, ac08cf8b04484d0a38afe5c6bc57420f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.935885 INFO::Fitting model to feature number 890, bba6e5f189ab17100885684b1c62d3e1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.966959 INFO::Fitting model to feature number 891, X76b202a36a7aa7bfa3b4f63e52e516cf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:37.997907 INFO::Fitting model to feature number 892, f0bc69d6d2e088f0e7895bbbe6de0d23
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.029418 INFO::Fitting model to feature number 893, X7b1effb26b50937a183a7a2f2bea3ae5
## 2026-03-04 18:37:38.060431 INFO::Fitting model to feature number 894, X9bea3a65d55fbec7b1de0fd8856c9bdb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.103567 INFO::Fitting model to feature number 895, d977c29f2689d12a8a72c44b02574e1c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.136822 INFO::Fitting model to feature number 896, a1ebedf636d67f3f4b9a40f02c5e702a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.168297 INFO::Fitting model to feature number 897, X4f6e37f2dfb5bcea2498189b425e0496
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.198973 INFO::Fitting model to feature number 898, X5d1bbb79e15c275d052188c91db4efce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.229552 INFO::Fitting model to feature number 899, a135c59c1c9deb0701d1f2eef52a8d34
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.260666 INFO::Fitting model to feature number 900, ba012602394a70abf571be8722eb7b33
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.291248 INFO::Fitting model to feature number 901, X9c8c8a3cfb2bab0fbc82c97e05bd537b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.322239 INFO::Fitting model to feature number 902, c843b591d20dcb1e33898346f82c79fa
## 2026-03-04 18:37:38.352596 INFO::Fitting model to feature number 903, d53cbb2969f0ef212c8cd9a67d1e0e1b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.384131 INFO::Fitting model to feature number 904, X143a46b132a920a6f068edbee257d9c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.415187 INFO::Fitting model to feature number 905, X1c401fa1144c19f7bb293293b446b27b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.44732 INFO::Fitting model to feature number 906, X4c880f5126065d48868d8720a7aeddd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.478366 INFO::Fitting model to feature number 907, X32359fd2ac112d65e1f7e9ae5ecd8b4a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.508603 INFO::Fitting model to feature number 908, fdcf22ebf60cd046b27e706f488c6c47
## 2026-03-04 18:37:38.539405 INFO::Fitting model to feature number 909, X55c9e708227432632f53f62dc5c8e7fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.569403 INFO::Fitting model to feature number 910, a0bea49cb183c27b663831c30239c2bb
## 2026-03-04 18:37:38.599406 INFO::Fitting model to feature number 911, de7dd3152ac71f84ae29298ecd828f06
## 2026-03-04 18:37:38.630105 INFO::Fitting model to feature number 912, X528f95a03adde9b5d7b6cb2805b140d7
## 2026-03-04 18:37:38.660979 INFO::Fitting model to feature number 913, X361d196fabf4a7e0f2b0dcfcebe5320d
## 2026-03-04 18:37:38.690492 INFO::Fitting model to feature number 914, X79f55d42a25d89fdbe2f7c4665679358
## 2026-03-04 18:37:38.720221 INFO::Fitting model to feature number 915, X555423d48ff522cb887cbef0a2370ff8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.750705 INFO::Fitting model to feature number 916, X3dcbcb2abd5f2838e7f226e4d09a1f2f
## 2026-03-04 18:37:38.780794 INFO::Fitting model to feature number 917, X8312da8265d49aac51fdd0d54b32d4fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:38.811934 INFO::Fitting model to feature number 918, a54c597c99769042fb4a11f69b70ac67
## 2026-03-04 18:37:38.842366 INFO::Fitting model to feature number 919, X5bc23fba732fd04d4986fd231269d6f3
## 2026-03-04 18:37:38.886279 INFO::Fitting model to feature number 920, X6c3a59bb3a352b8b87bd7a1109ec5548
## 2026-03-04 18:37:38.919794 INFO::Fitting model to feature number 921, X53512fd8672ee91eef4eea370182e655
## 2026-03-04 18:37:38.951001 INFO::Fitting model to feature number 922, X5f4079304db67863bd1a742918ef8182
## 2026-03-04 18:37:38.981989 INFO::Fitting model to feature number 923, X3133ff9657a0f895b8cdd9b72fd49b05
## 2026-03-04 18:37:39.01313 INFO::Fitting model to feature number 924, X49bfc506d5b3743431bb0eeae6069917
## 2026-03-04 18:37:39.04365 INFO::Fitting model to feature number 925, c3daa62f8c41fb7c3a0d2ea92fc6b4af
## 2026-03-04 18:37:39.074363 INFO::Fitting model to feature number 926, cb7c601e1ff33d0a8e7ec6a8a7cb276b
## 2026-03-04 18:37:39.105563 INFO::Fitting model to feature number 927, X167e93c6b19b4454ef559b9fbfb409cb
## 2026-03-04 18:37:39.136357 INFO::Fitting model to feature number 928, X0b842f865339edc07fdf5ff3bd3f8379
## 2026-03-04 18:37:39.16729 INFO::Fitting model to feature number 929, X3dc628a32ca4dd564074799ece230f67
## 2026-03-04 18:37:39.198067 INFO::Fitting model to feature number 930, daabbbd5db10cf6cc8312d919bbcb924
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.229808 INFO::Fitting model to feature number 931, f4f07e406c1e85a13f899ebf9b95e29d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.26082 INFO::Fitting model to feature number 932, X146d9408efd14ffff6e63556d93b5108
## 2026-03-04 18:37:39.290929 INFO::Fitting model to feature number 933, X6a5f53fae6f84847d447098ae4e8f98f
## 2026-03-04 18:37:39.320975 INFO::Fitting model to feature number 934, X1efbe8e6fc06fbbd97b5f4bfcf4a4ed0
## 2026-03-04 18:37:39.350845 INFO::Fitting model to feature number 935, a90a928ed331d7e3c701c511619e214d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.382269 INFO::Fitting model to feature number 936, df23731588f3884c39b70b5d6e7811bb
## 2026-03-04 18:37:39.412559 INFO::Fitting model to feature number 937, fa602366e449f46a0d8cfaf64a48e353
## 2026-03-04 18:37:39.443242 INFO::Fitting model to feature number 938, b5832c657c49cbf026cd29bacca46fff
## 2026-03-04 18:37:39.473854 INFO::Fitting model to feature number 939, f55431826646e45d3a5a05104fda4c2c
## 2026-03-04 18:37:39.504255 INFO::Fitting model to feature number 940, a4c443d483ee55325f29b06c5d45e011
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.53453 INFO::Fitting model to feature number 941, a833bd321294a6cebf81e23cc647560b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.565371 INFO::Fitting model to feature number 942, ebb6d47df55f900c7264fde90ecc0085
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.595929 INFO::Fitting model to feature number 943, aa42a80bfe760ce5ef13e31a63ea9b2e
## 2026-03-04 18:37:39.637709 INFO::Fitting model to feature number 944, X2e3b65e7f787eb3a955b046b3a399821
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.670561 INFO::Fitting model to feature number 945, X976ebe9b8c37410b4afb4acb24205cb9
## 2026-03-04 18:37:39.701854 INFO::Fitting model to feature number 946, X2470670cc5f1cc0d346e708bd238260a
## 2026-03-04 18:37:39.731869 INFO::Fitting model to feature number 947, X90e925964999f0db59e87db25b2efa8f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.762762 INFO::Fitting model to feature number 948, X942a41f0147db28c830cb529bd03add1
## 2026-03-04 18:37:39.793362 INFO::Fitting model to feature number 949, X9a2da91ca6be4cdeed6cb6edb5f6b4ac
## 2026-03-04 18:37:39.823796 INFO::Fitting model to feature number 950, a9fa05aa77b1ee436a3f27c2856db7a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.854889 INFO::Fitting model to feature number 951, dc126ebf3a1c182a034f5acbf4e22c14
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.892128 INFO::Fitting model to feature number 952, df76282db54adacc08160fc0eaeb6e7d
## 2026-03-04 18:37:39.923644 INFO::Fitting model to feature number 953, X271b453412999e72adabfd71eb789300
## 2026-03-04 18:37:39.955514 INFO::Fitting model to feature number 954, d36d14ae5dbdbd3d91e8d5fa3d5ec78d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:39.987575 INFO::Fitting model to feature number 955, e1b6cc549ae04f1c81727b330eb8bdb4
## 2026-03-04 18:37:40.01867 INFO::Fitting model to feature number 956, X6bc05eda0c85a2a0767635570b137a45
## 2026-03-04 18:37:40.049956 INFO::Fitting model to feature number 957, e09fe686537ebbc38c983c79018ffccf
## 2026-03-04 18:37:40.081027 INFO::Fitting model to feature number 958, X5d70be0acfe412990838d37001eb24b2
## 2026-03-04 18:37:40.111362 INFO::Fitting model to feature number 959, e4efa7d4be218d347e0d930295e832a3
## 2026-03-04 18:37:40.142483 INFO::Fitting model to feature number 960, X0dded72ceaf4dca73efed63cfb19812b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.173854 INFO::Fitting model to feature number 961, e1d4f791ba1642d0da77adc5de295325
## 2026-03-04 18:37:40.206803 INFO::Fitting model to feature number 962, f244f35600a551dfb3eeeaf86d473de8
## 2026-03-04 18:37:40.238685 INFO::Fitting model to feature number 963, X498291ee6f8d735d5f0353c70be3166c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.27067 INFO::Fitting model to feature number 964, X4e661aa07dec1515abb2615f59b388af
## 2026-03-04 18:37:40.301944 INFO::Fitting model to feature number 965, X50ec0d0bf1a352d385b6338057ac9bb9
## 2026-03-04 18:37:40.332565 INFO::Fitting model to feature number 966, X5a1b8ebc20cd26618173eca7a37685f1
## 2026-03-04 18:37:40.368255 INFO::Fitting model to feature number 967, f3beee66d95cd06a16b1919dde141bd0
## 2026-03-04 18:37:40.400286 INFO::Fitting model to feature number 968, X8c266bccda7a4c3fb93cad5f8e79fc08
## 2026-03-04 18:37:40.445867 INFO::Fitting model to feature number 969, f52ce8860f150e26d9ed45e91770ee38
## 2026-03-04 18:37:40.479205 INFO::Fitting model to feature number 970, d8e67bf70feea21e2daeddf4410172a1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.512903 INFO::Fitting model to feature number 971, X50fa0bc6d88a3894aead1fc528d41058
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.544482 INFO::Fitting model to feature number 972, c555bb81f49b80858ee015a9eea42aef
## 2026-03-04 18:37:40.575424 INFO::Fitting model to feature number 973, X91c734b3c610cf63571b5a22f5d4dc2e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.607048 INFO::Fitting model to feature number 974, a14e904a7d2d1daace916fdf9856d31d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.63838 INFO::Fitting model to feature number 975, X94d7d95a02e3395be3ab28a0fed2872c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.669419 INFO::Fitting model to feature number 976, fb7f32b378722a5ca83e91cb41927eef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.701383 INFO::Fitting model to feature number 977, X9ebb4f9911373eb720c9bff4af833ddd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.73315 INFO::Fitting model to feature number 978, X667b35160a1f273c609d48a98ceb4d55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.765564 INFO::Fitting model to feature number 979, X6d2da75b05c499e6c1acb7351f0cf660
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.797461 INFO::Fitting model to feature number 980, ca9f3b34c45626302fa1fec4e105314e
## 2026-03-04 18:37:40.82797 INFO::Fitting model to feature number 981, X260de73cfc2af53494272275576eaa65
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.859763 INFO::Fitting model to feature number 982, e97000fb903a8e88351eb8e8d1c9d42c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.893996 INFO::Fitting model to feature number 983, X363e9dd3a114361345cb1ce36e8dd520
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.924897 INFO::Fitting model to feature number 984, X1e0f313eacdce85ff324dd5823aedff5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.961703 INFO::Fitting model to feature number 985, X8b147b77e5de5f3d3c1f6c478ffbab2a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:40.998761 INFO::Fitting model to feature number 986, X38d3539bc11b4977d1d3f32e1de0dbf4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.030843 INFO::Fitting model to feature number 987, a5f83ed4dc472953510501a5f899c2a1
## 2026-03-04 18:37:41.064144 INFO::Fitting model to feature number 988, X5cc9033bb1a702b19b9a7b2e44ecd851
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.098103 INFO::Fitting model to feature number 989, X3b28c382abaa17c7f1f339967619ddd4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.136027 INFO::Fitting model to feature number 990, bca782ea0c57a7716ed3aea3f5376da9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.168835 INFO::Fitting model to feature number 991, X2a57beb2f11af3ae774c0df31aac1275
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.207615 INFO::Fitting model to feature number 992, X5a3930b7b374d9d2019a9b26531ffb6d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.257715 INFO::Fitting model to feature number 993, X7fc2f360dc73536bb0281761c75399d5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.29577 INFO::Fitting model to feature number 994, X5d2b45d606b5eafbcc1e9f87452d3135
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.329525 INFO::Fitting model to feature number 995, X66deedccb8ee4f471c0029a94b18ec49
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.368466 INFO::Fitting model to feature number 996, X4aa6ba9c7da5f33af213c02d65ecb391
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.401466 INFO::Fitting model to feature number 997, X8da71749ef2ac23a056fc49fdda9131e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.43875 INFO::Fitting model to feature number 998, c90f38304baf5fc52b3332e63a8016c2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.472408 INFO::Fitting model to feature number 999, X396daf5d94b484aba7fcc4d8c0223759
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.502772 INFO::Fitting model to feature number 1000, X7bb6e29f06620ad6f4f0e7478a1fea5f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.53583 INFO::Fitting model to feature number 1001, X8463aae5183f78a44bcd2369cdf5ed93
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.567225 INFO::Fitting model to feature number 1002, X8aa5eb1b9d4e710d35223d1ddaaccc22
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.605184 INFO::Fitting model to feature number 1003, X8fb79654f3000a52e5b1be221f49106e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.636593 INFO::Fitting model to feature number 1004, ca76a5d3cef024cc5c9c1ea5b19d77fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.667752 INFO::Fitting model to feature number 1005, X1aae042eddcd396bbe739d785fb9df7e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.69895 INFO::Fitting model to feature number 1006, X2b52a41587c8ed22703a2bc2d978f584
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.731955 INFO::Fitting model to feature number 1007, X389f789d4f3fae4c7040028704fbbc03
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.762403 INFO::Fitting model to feature number 1008, X41d3c7fe72f0ac46769a543a2322c72b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.794168 INFO::Fitting model to feature number 1009, ba4cf501f96340f482a9525b4ab1f69d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.827006 INFO::Fitting model to feature number 1010, X2a4d95ff9ae281d8c13c486f3772d1bb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.858793 INFO::Fitting model to feature number 1011, X5f05f4bb091a1d67b7c724a6aaffdcc4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.89162 INFO::Fitting model to feature number 1012, X7705f0fc350c0d3189688a4858948d1c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.923063 INFO::Fitting model to feature number 1013, a79885cf6c4f4e70eb18d351e69b8b78
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.954967 INFO::Fitting model to feature number 1014, a907f48fc8e04ed9ae516022117a446d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:41.986795 INFO::Fitting model to feature number 1015, ad901777a802617307dca87d69e1f57d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.017302 INFO::Fitting model to feature number 1016, X751ccf1e1783168abe381e454cd7f817
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.049926 INFO::Fitting model to feature number 1017, f4216c4a78111523f17cbc7de20025b6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.100459 INFO::Fitting model to feature number 1018, X3528210d91f62915dfd614d5f8c5246d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.136534 INFO::Fitting model to feature number 1019, c1e2db64acc21a102c6a7d833d730a67
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.169606 INFO::Fitting model to feature number 1020, X48e84a9b22d27e8f7b2854e2d0e6e151
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.201914 INFO::Fitting model to feature number 1021, ff2a02b6f268f3d8c34998b06e396f46
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.232897 INFO::Fitting model to feature number 1022, cc55eefd97b5534a346ae449625be5dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.264739 INFO::Fitting model to feature number 1023, X50e17521cf480fe34b2801f90cee6178
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.297413 INFO::Fitting model to feature number 1024, X6042f5b4bcfc4ca01bf9dba3a94f82de
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.329664 INFO::Fitting model to feature number 1025, X6079f6c9472d90f3fe6b93c2786229df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.361299 INFO::Fitting model to feature number 1026, a11644f2e991a5a53bc43924a8e0f75f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.392856 INFO::Fitting model to feature number 1027, cb5f94b1638b84a247685c32724e4b10
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.430413 INFO::Fitting model to feature number 1028, X2faf5487b3d0d36915ccbac8006e96cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.463622 INFO::Fitting model to feature number 1029, X961327649d0d8b1a3f1d8506f7abf9a3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.495025 INFO::Fitting model to feature number 1030, X8e867981ffca8750b7a43cb387f9f135
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.525499 INFO::Fitting model to feature number 1031, X35b607864efdc9e272585bca1d40016e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.55619 INFO::Fitting model to feature number 1032, X1a074fc9ba7138aca27f372e1b4a000e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.587585 INFO::Fitting model to feature number 1033, X92e8b9e50e7601db88a32a3b69963bb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.619871 INFO::Fitting model to feature number 1034, X45a2ac2f0e2a4c246a56beecda0e2dbd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.651259 INFO::Fitting model to feature number 1035, fab3fc1d161c45c577a91e58006417b1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.689012 INFO::Fitting model to feature number 1036, X1cffb93453b8e2a2a884422848606aa6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.726126 INFO::Fitting model to feature number 1037, X84540ed2eebf1efb8e0d8d1e40dabffd
## 2026-03-04 18:37:42.762223 INFO::Fitting model to feature number 1038, d012209e84dc07f3d3aeccff20da3d51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.79222 INFO::Fitting model to feature number 1039, X50a43b2ec96d6b8bb0f9a8680a05f379
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.826983 INFO::Fitting model to feature number 1040, X567a407c520a5cf74f4cd9bd5e13b199
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.857849 INFO::Fitting model to feature number 1041, X66c3e26fde9a995062499289924775a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.901551 INFO::Fitting model to feature number 1042, X7194882e587298e6c178faf29b1572bc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.939502 INFO::Fitting model to feature number 1043, d75f983da6423ae3900c33e44a566605
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:42.972595 INFO::Fitting model to feature number 1044, X89651bf4f451a1f07f6ab1983f484fce
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.004782 INFO::Fitting model to feature number 1045, X3d9449debdcefa1ef7caea17e9aee6d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.036424 INFO::Fitting model to feature number 1046, X58a5ed3a034ea45c3da6ccdfae4f3659
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.069276 INFO::Fitting model to feature number 1047, X467bd9bb0da5518cccb8ef64d3427bd9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.106845 INFO::Fitting model to feature number 1048, X4ba356cfaa02c5d2668304cf8a75d98e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.137718 INFO::Fitting model to feature number 1049, X59f2c893c976b2bed0ce5f3168b51b5e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.168583 INFO::Fitting model to feature number 1050, d49581d80b688f6ea99b76cae93334ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.199559 INFO::Fitting model to feature number 1051, e0abea0159042cd283f2560997bd6dd7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.238506 INFO::Fitting model to feature number 1052, fc77df84a1ad92454b15a3a729ca9fab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.268685 INFO::Fitting model to feature number 1053, X5de6b44b38194149036db56946e96b38
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.299183 INFO::Fitting model to feature number 1054, X883c6c1b596bf1c5cf45a00427cebbde
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.329916 INFO::Fitting model to feature number 1055, be1fe3655bebd478e3fe789cc433ec14
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.360814 INFO::Fitting model to feature number 1056, c7aa42e4a69d17cdb7a5374e8c510136
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.392631 INFO::Fitting model to feature number 1057, f3f9d3323e4d2f15298aaed5439af56a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.423736 INFO::Fitting model to feature number 1058, X765c815c5453e9ed5dc420c8119c5d01
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.455609 INFO::Fitting model to feature number 1059, d36986d02dbdc4507e68cd32d674ac05
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.486132 INFO::Fitting model to feature number 1060, d4c443b446021feb34aad85d999b11cb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.515921 INFO::Fitting model to feature number 1061, X08b4991c00fa39be14bc8f7327e213a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.546234 INFO::Fitting model to feature number 1062, X4d076914382017f211091f07aa8db632
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.576015 INFO::Fitting model to feature number 1063, X944c1fe867c513e7c1e38241d24352c1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.609131 INFO::Fitting model to feature number 1064, b7af501666cfd4d95b15a259a0d33534
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.641704 INFO::Fitting model to feature number 1065, d49f0e0482144b9ac84c3a0dad81f3dc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.679901 INFO::Fitting model to feature number 1066, ea17d1fba365a308f925eb11b72264bf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.722648 INFO::Fitting model to feature number 1067, X53f66ce1db04a14261e921e5f8dd43db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.754512 INFO::Fitting model to feature number 1068, X54ea53f5747082b754171686efe7ae35
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.785191 INFO::Fitting model to feature number 1069, X93628c71a8cacf127ce7fae55d40c1d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.814874 INFO::Fitting model to feature number 1070, a0d658715dbad152018ebec95f5b2a13
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.845893 INFO::Fitting model to feature number 1071, faeb2051952acacd1ee1deff131c3bb0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.878092 INFO::Fitting model to feature number 1072, X58d7505e7eb30d2e9a47244a921e3fae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.908861 INFO::Fitting model to feature number 1073, X85b0ab39a82b695768422d01e88a25d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:43.939411 INFO::Fitting model to feature number 1074, cbb718ece292de4ff64d526aa5604372
## 2026-03-04 18:37:43.970545 INFO::Fitting model to feature number 1075, X5ab675c9e714d085b6b113d8892cd685
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.000386 INFO::Fitting model to feature number 1076, X46b863013a51fc2faa3e8cbf8c77271a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.031094 INFO::Fitting model to feature number 1077, cd2bb0b70d041ae5bb64067b2ab990ec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.062222 INFO::Fitting model to feature number 1078, b1c849febb33944267c7047a39244974
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.091668 INFO::Fitting model to feature number 1079, X5749de9d92f48616e003d42550b8f79d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.121567 INFO::Fitting model to feature number 1080, X51ad6d6b27b7f1551bd818f9f253fb61
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.151408 INFO::Fitting model to feature number 1081, e78437188624ce623ae9bcd812e72d2b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.188984 INFO::Fitting model to feature number 1082, X0baa3c3e2e6042eede193ac9f25d53a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.226155 INFO::Fitting model to feature number 1083, X9c3e3ed907988deb778c5f0cdec6de0a
## 2026-03-04 18:37:44.262474 INFO::Fitting model to feature number 1084, X0f5e024743e654ca7872273094fbed5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.298698 INFO::Fitting model to feature number 1085, X784e1de605fa06ef3726907163c5027e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.331992 INFO::Fitting model to feature number 1086, X9e5b25a18c0bf19a5b288654dc4ea008
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.362654 INFO::Fitting model to feature number 1087, X1004ca7e65c7a7d359941799d33a67b3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.393018 INFO::Fitting model to feature number 1088, X8f57a9d96db41075fa1ba4c4328e26d7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.425196 INFO::Fitting model to feature number 1089, X04ede807f45e5ac5b121a53a237f830b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.458108 INFO::Fitting model to feature number 1090, X3ead1fa5a380796e8a1f7c74c5212bc0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.509582 INFO::Fitting model to feature number 1091, d8edc720feb6c43730bcfe12a789b5e8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.543151 INFO::Fitting model to feature number 1092, X5c10c0918fb25e77e97086722698f1ba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.575738 INFO::Fitting model to feature number 1093, X2d84007ca1c2c3f513156f6f520716fc
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.607601 INFO::Fitting model to feature number 1094, X9f60b24425f186b088fe6280c81a8388
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.63832 INFO::Fitting model to feature number 1095, cfc746d1cf45fc692c0e53c4c7f671d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.668429 INFO::Fitting model to feature number 1096, d874902a8d6182ca9798bc486ebb4a99
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.705407 INFO::Fitting model to feature number 1097, f392db5a807da4e8a2e708bd84e04ca0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.74058 INFO::Fitting model to feature number 1098, X6d3c164b889d28be2ec50d08fb62322e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.774423 INFO::Fitting model to feature number 1099, c18ca53cea68a2a7144fe36bce89340b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.806001 INFO::Fitting model to feature number 1100, ce68fc44b24d1212ae41371700862b72
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.837071 INFO::Fitting model to feature number 1101, X455814e0df838e53217a329758d86744
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.86827 INFO::Fitting model to feature number 1102, X68612badfa6c422975b70cf092a1b873
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.898626 INFO::Fitting model to feature number 1103, X7d5e0a18434e6d972b5c2dc28b8bf51d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.929326 INFO::Fitting model to feature number 1104, X86a4904762e51b23c5eac255935f0534
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.959787 INFO::Fitting model to feature number 1105, X9060dd97d853d275ec247f2db98c265f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:44.991136 INFO::Fitting model to feature number 1106, X49b65a3ca234740512a8dc5b92ccf5a0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.02259 INFO::Fitting model to feature number 1107, X4f9ed6ab77896830974fd4df52525f76
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.054309 INFO::Fitting model to feature number 1108, X674f8b7dba9f080323fcec4292eb16e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.091676 INFO::Fitting model to feature number 1109, X1dbdd4ba3c7ae9b85cea8f0f4cfa8a10
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.122297 INFO::Fitting model to feature number 1110, X9590c032d57e4f5568b0eafa67708b3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.159614 INFO::Fitting model to feature number 1111, b39e9d4cb76e6b56dc9999a92c5283d4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.19177 INFO::Fitting model to feature number 1112, e978d95f90a6f389da9b57cf0868b847
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.223111 INFO::Fitting model to feature number 1113, X198065e0411746d3c31cbc810a1c8b6b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.253541 INFO::Fitting model to feature number 1114, bc25945dde85e8f2f55c642669695015
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.296553 INFO::Fitting model to feature number 1115, ce3868416d546ecf20e2c204da9a1828
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.331814 INFO::Fitting model to feature number 1116, X8910dfaae22f25dfe7f0679ca7c94618
## 2026-03-04 18:37:45.364522 INFO::Fitting model to feature number 1117, e39c0f2f2bf6163e4f2da67826b09132
## 2026-03-04 18:37:45.395543 INFO::Fitting model to feature number 1118, X7ad2b0d2a4533eb2ce608b7852890f98
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.427689 INFO::Fitting model to feature number 1119, X597728f24a67cc2770bedadf83591fdd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.459605 INFO::Fitting model to feature number 1120, X8cfb2b6f68d3ac591ad8d67f325aecc5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.489969 INFO::Fitting model to feature number 1121, a48784576958e537b25f9fa5d4e3f98e
## 2026-03-04 18:37:45.520454 INFO::Fitting model to feature number 1122, X267edef74cf84d7548b0d3f88cb5ff89
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.551833 INFO::Fitting model to feature number 1123, cc9703e82e512181f882f6bf229cb138
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.583656 INFO::Fitting model to feature number 1124, d91a7f73b26a7d13a51d9f3e4fe3cc06
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.615847 INFO::Fitting model to feature number 1125, X6331ea1bb0c3c736beaf34fd6284fd31
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.64769 INFO::Fitting model to feature number 1126, X87eb2526a517ff1a76336ae23315c926
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.679318 INFO::Fitting model to feature number 1127, X2e27f1f31ebbd4b4c0ff20a7feb121a8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.710659 INFO::Fitting model to feature number 1128, bcfe7158c79ee0946dc47793ffbf0549
## 2026-03-04 18:37:45.74137 INFO::Fitting model to feature number 1129, X27a395b590a278e9085d7157656a37ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.772575 INFO::Fitting model to feature number 1130, f6fd2a76558c3266f11fe4c856f929e0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.804347 INFO::Fitting model to feature number 1131, afde775f8b52d66edb513773b710a8e0
## 2026-03-04 18:37:45.835168 INFO::Fitting model to feature number 1132, X9939ab0bed6fd6fa3c489cfa467cf9c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.866542 INFO::Fitting model to feature number 1133, X11f98106eb45262aa2b864ed085f3369
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.897301 INFO::Fitting model to feature number 1134, X19609c541e304d9414eeed6a276f038a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.928123 INFO::Fitting model to feature number 1135, X33151a1928fbfcff924af5d627bc553b
## 2026-03-04 18:37:45.959206 INFO::Fitting model to feature number 1136, a2224b2e402597ee34ee3fdcd9fb7298
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:45.99 INFO::Fitting model to feature number 1137, X29087b19e0622aa58bc91d8e521126e4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.021583 INFO::Fitting model to feature number 1138, fcde5c46c10eec800e6b38a39ea9e0b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.052943 INFO::Fitting model to feature number 1139, X219c33e8cda1089efcbe919a3694dbf6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.098863 INFO::Fitting model to feature number 1140, X29d5f5579a251818c2a9c49b911782ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.131432 INFO::Fitting model to feature number 1141, X6e72215f2c39c515db5fa57a00ba473b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.16221 INFO::Fitting model to feature number 1142, a6583feb07439f47775cf6216ebe3e9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.192575 INFO::Fitting model to feature number 1143, X2d5028d34b33f496c8ea9d664e3d9f24
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.223328 INFO::Fitting model to feature number 1144, X9c0ac0fed41ad08dcd91f61b0440f4bf
## 2026-03-04 18:37:46.253907 INFO::Fitting model to feature number 1145, ff5ab5918e8b9c763a2f1c3cf33d06f1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.283974 INFO::Fitting model to feature number 1146, X67e28f844816795eb9d42efcff09c29b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.314304 INFO::Fitting model to feature number 1147, X79fdee6143db1981342333ae4a977269
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.345966 INFO::Fitting model to feature number 1148, f5aa4ea6105b38599df1b4c9f2c526ab
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.377279 INFO::Fitting model to feature number 1149, d315b9630a2ca429ed87f6ac6856f505
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.409875 INFO::Fitting model to feature number 1150, X0186f5ab59b4d03b51ca5bb5d69fec54
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.442257 INFO::Fitting model to feature number 1151, X334079baf11f432d414d30bd1e0f09e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.47362 INFO::Fitting model to feature number 1152, X960304d4b9cb07256cceb3a66af3234f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.504479 INFO::Fitting model to feature number 1153, X9385f532468f13f7aeef2e694329d709
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.534693 INFO::Fitting model to feature number 1154, f5a9b9639b6ba9307b5166734c9a69fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.565415 INFO::Fitting model to feature number 1155, X26c053114523080cbc3b791b07e8b937
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.596042 INFO::Fitting model to feature number 1156, X54d7f7e9bf7bd4e00d059e3a9f80b832
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.627396 INFO::Fitting model to feature number 1157, X91237a8dd758dacb643569fe824f48b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.658015 INFO::Fitting model to feature number 1158, ca7f6db958f83d9a7daac017d413bcc8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.68896 INFO::Fitting model to feature number 1159, X586c90b4d2a884eeb3bbdfb2c129e908
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.720168 INFO::Fitting model to feature number 1160, X71ed20fd789eccdde477d5b95c58ce98
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.751073 INFO::Fitting model to feature number 1161, X76f64e5e9cbda0476bbb42ad34e4cdf9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.781922 INFO::Fitting model to feature number 1162, c46904bcacd273e0a9a37384d6e45a26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.818238 INFO::Fitting model to feature number 1163, cfcd2110898f3d3051894f249b744767
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.865416 INFO::Fitting model to feature number 1164, X10096fb55ea863485a884e2bf9542163
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.900633 INFO::Fitting model to feature number 1165, c7e8689997df61fe02151acdef28ef63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.932944 INFO::Fitting model to feature number 1166, dac90bcd762785b630843d62a7da6c9a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.963813 INFO::Fitting model to feature number 1167, eed08fcc06d959756c2dd08f4ae201c9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:46.994552 INFO::Fitting model to feature number 1168, f5d15fd366b9bb809bbf1a53a126b694
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.024548 INFO::Fitting model to feature number 1169, X045f2c5846fb21d9f2dfbc9ec3bf64e1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.056032 INFO::Fitting model to feature number 1170, X1cb08a81bec75cf1cb6a5347637f32e7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.087116 INFO::Fitting model to feature number 1171, X393a57d6c6b02ca59aa23b7a88a4b773
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.117911 INFO::Fitting model to feature number 1172, X478a6c05d04db6a58eaae7bb27eb04ee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.149379 INFO::Fitting model to feature number 1173, X522fc4ccc38089e8f0d8d4a57d4e82f8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.180035 INFO::Fitting model to feature number 1174, b909bd96ae3a423bbef71d7add0c86df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.211711 INFO::Fitting model to feature number 1175, c937f4adda1f009333dd511ae285fff7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.242917 INFO::Fitting model to feature number 1176, d11d08c6a579ddd3cc910338ebbe6547
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.273758 INFO::Fitting model to feature number 1177, ec054373da44610fd5b9a258e5f793a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.304729 INFO::Fitting model to feature number 1178, f0af9b7703a95bfa52c7ebcfc5a94b82
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.335802 INFO::Fitting model to feature number 1179, f1f7e45805f74287c1564249d80cafa2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.366672 INFO::Fitting model to feature number 1180, f70c8600851d3d162e26ad9cf59675fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.398332 INFO::Fitting model to feature number 1181, f81c4bee45cfd34ed1fe726d0c844fb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.430076 INFO::Fitting model to feature number 1182, X3a323cf260c90c3306639d009f11b7fc
## 2026-03-04 18:37:47.469724 INFO::Fitting model to feature number 1183, b7513d5b69299f329da60ff370d4acf9
## 2026-03-04 18:37:47.500284 INFO::Fitting model to feature number 1184, X086c111d7310df390ca2ce0ac8144156
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.530176 INFO::Fitting model to feature number 1185, X462983a97faff888d34ac78f5f92825f
## 2026-03-04 18:37:47.560068 INFO::Fitting model to feature number 1186, X3a50580ca187917dbd5d8dd58f5c6a48
## 2026-03-04 18:37:47.5898 INFO::Fitting model to feature number 1187, X62ebbf17f633335b2fc737ebb581df16
## 2026-03-04 18:37:47.619597 INFO::Fitting model to feature number 1188, X727da14bd94f15b7dee15dbe146c0da9
## 2026-03-04 18:37:47.665759 INFO::Fitting model to feature number 1189, X3e737e299b0235bc5df1bc95b9d551cc
## 2026-03-04 18:37:47.698116 INFO::Fitting model to feature number 1190, X6680fc094b1ae654169e3f6596065387
## 2026-03-04 18:37:47.729223 INFO::Fitting model to feature number 1191, X54339f8b1100dc127f86478b017555bb
## 2026-03-04 18:37:47.759862 INFO::Fitting model to feature number 1192, X609000883c374a04e356bda18b419104
## 2026-03-04 18:37:47.790065 INFO::Fitting model to feature number 1193, X95b88fdb11b48228c744237ea2277420
## 2026-03-04 18:37:47.819922 INFO::Fitting model to feature number 1194, X9e2a90893c63429f7974a93724f58473
## 2026-03-04 18:37:47.850426 INFO::Fitting model to feature number 1195, X3646994e928de61adcc5460ca0df9cd0
## 2026-03-04 18:37:47.880392 INFO::Fitting model to feature number 1196, X4dbb7bbc94b14db1ffad5c6e3de8341c
## 2026-03-04 18:37:47.911192 INFO::Fitting model to feature number 1197, f521618f99a90bb8c936cb315d5b33cd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:47.942758 INFO::Fitting model to feature number 1198, a473147383e13716f9365ff1926aa6a6
## 2026-03-04 18:37:47.974694 INFO::Fitting model to feature number 1199, f2c1249d1644237604788e8198b48d31
## 2026-03-04 18:37:48.005087 INFO::Fitting model to feature number 1200, fe3bfbd4d173653203a9d6d256d92240
## 2026-03-04 18:37:48.035582 INFO::Fitting model to feature number 1201, X077b2ad66183e1ab78bd0a27b5170970
## 2026-03-04 18:37:48.070344 INFO::Fitting model to feature number 1202, X39a893e34f9c30e78354a90670c464ff
## 2026-03-04 18:37:48.099975 INFO::Fitting model to feature number 1203, X89c032289bb7d50431817cd8ad8fb3e3
## 2026-03-04 18:37:48.130377 INFO::Fitting model to feature number 1204, X65970e707cca45512e17863477cbf9fd
## 2026-03-04 18:37:48.160855 INFO::Fitting model to feature number 1205, X9e7f26a7544795db8dd3f0dcad58d1d9
## 2026-03-04 18:37:48.191372 INFO::Fitting model to feature number 1206, X738e89f33bb2d9a90cc65b6e38395094
## 2026-03-04 18:37:48.229266 INFO::Fitting model to feature number 1207, e8d8f54456fddc89b8d3d072dfd0dae2
## 2026-03-04 18:37:48.259703 INFO::Fitting model to feature number 1208, X8e068f286e4459a2f92d83bf8fee4c75
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:48.29062 INFO::Fitting model to feature number 1209, ff3ef656c9c5c1e8d9c808185a095229
## 2026-03-04 18:37:48.322912 INFO::Fitting model to feature number 1210, X44816a4453af39d938da26fe98668bf6
## 2026-03-04 18:37:48.353557 INFO::Fitting model to feature number 1211, dced475a8c58af383bb71156ede21a28
## 2026-03-04 18:37:48.38415 INFO::Fitting model to feature number 1212, X1fb31804045eed5be6772950d92bbea1
## 2026-03-04 18:37:48.425716 INFO::Fitting model to feature number 1213, X2d4eff8bb4b08b2ad545641114a1f95d
## 2026-03-04 18:37:48.461689 INFO::Fitting model to feature number 1214, X9b1e1af81c490f7167ba1f0f3a5c27f0
## 2026-03-04 18:37:48.494204 INFO::Fitting model to feature number 1215, X64c9855e3fb0a940f7160dc982b493aa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:48.525023 INFO::Fitting model to feature number 1216, X8375b729fbec2ab3e94c4e34a6b89601
## 2026-03-04 18:37:48.556393 INFO::Fitting model to feature number 1217, adfe3978e202f9a097a1f8cce2fd35b7
## 2026-03-04 18:37:48.587446 INFO::Fitting model to feature number 1218, c52e54c93b3cc50a746340ebdd80aa55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:48.618244 INFO::Fitting model to feature number 1219, X5b1b6becbd52352e6a9618ddadf0d8d9
## 2026-03-04 18:37:48.648514 INFO::Fitting model to feature number 1220, a8fda41c9782acfa04a7d7df17761e59
## 2026-03-04 18:37:48.67971 INFO::Fitting model to feature number 1221, X4e4884420b31400d35db2a1417fa4ec7
## 2026-03-04 18:37:48.717316 INFO::Fitting model to feature number 1222, f236ac69ec7755b0247018f9b997b504
## 2026-03-04 18:37:48.748083 INFO::Fitting model to feature number 1223, X5bcdc43e09669d72979a94199b86a102
## 2026-03-04 18:37:48.778949 INFO::Fitting model to feature number 1224, dd0a2ab3c468bf12f8a8de815ab5dd63
## 2026-03-04 18:37:48.809561 INFO::Fitting model to feature number 1225, X9de219990e1666e0ba94f08161b49e4c
## 2026-03-04 18:37:48.840133 INFO::Fitting model to feature number 1226, ac6cc9a3a10683c803042c539a9caea1
## 2026-03-04 18:37:48.870587 INFO::Fitting model to feature number 1227, b683bb4bfec4e83e90550a48ad78ef9a
## 2026-03-04 18:37:48.909086 INFO::Fitting model to feature number 1228, e720001f2f9616203ce0cd35245a6218
## 2026-03-04 18:37:48.938735 INFO::Fitting model to feature number 1229, f50f30b5af2b1813e816f5aeedb1f22f
## 2026-03-04 18:37:48.968848 INFO::Fitting model to feature number 1230, X596e2f7be3a79400444a18242148fe31
## 2026-03-04 18:37:48.999477 INFO::Fitting model to feature number 1231, X8ccc17b289af3286eacf8da751dd6ffc
## 2026-03-04 18:37:49.029477 INFO::Fitting model to feature number 1232, fd7a5206567863e6bbc88ec1431ce83a
## 2026-03-04 18:37:49.060147 INFO::Fitting model to feature number 1233, X58448ec65533e295cd04dcfce9c81137
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:49.090482 INFO::Fitting model to feature number 1234, d6960b8def5f8b8cd6e93a90054cd2cb
## 2026-03-04 18:37:49.120309 INFO::Fitting model to feature number 1235, X083129e51191b2bf89f87bfafec40b48
## 2026-03-04 18:37:49.149803 INFO::Fitting model to feature number 1236, d923157d249ba936676cb5b0e460ba21
## 2026-03-04 18:37:49.17931 INFO::Fitting model to feature number 1237, X8b00ab07a4a5072a42e9ab7ab540313a
## 2026-03-04 18:37:49.224729 INFO::Fitting model to feature number 1238, f2480634f9b255a885f0577d0c84663d
## 2026-03-04 18:37:49.256722 INFO::Fitting model to feature number 1239, X2ce08b3a482c56941e6838df1a1d0db3
## 2026-03-04 18:37:49.288484 INFO::Fitting model to feature number 1240, X0c6a4054470a42e000bb0b69e5aa1a4e
## 2026-03-04 18:37:49.327049 INFO::Fitting model to feature number 1241, f1d223031e540e85dc5b60c84bbf80c6
## 2026-03-04 18:37:49.357947 INFO::Fitting model to feature number 1242, d52068433df4203a90b8e9da9e6e4dd2
## 2026-03-04 18:37:49.390721 INFO::Fitting model to feature number 1243, X1244fb68a2ab9cacab3a01719181ef35
## 2026-03-04 18:37:49.420502 INFO::Fitting model to feature number 1244, X25410db3bc0d1fa7ca52ce54385cb0af
## 2026-03-04 18:37:49.451538 INFO::Fitting model to feature number 1245, d11910302f2ea2a2dbca76d22dc4bcae
## 2026-03-04 18:37:49.481391 INFO::Fitting model to feature number 1246, X838980cc7edfa38dbcd8a24eeb97f9f0
## 2026-03-04 18:37:49.512918 INFO::Fitting model to feature number 1247, X19855f95f67ddb3a43733ce20d0e6d62
## 2026-03-04 18:37:49.551032 INFO::Fitting model to feature number 1248, X1dee71f32a863cdbaac0e10e0c89d5ec
## 2026-03-04 18:37:49.587483 INFO::Fitting model to feature number 1249, X618fd6e0f9e5467a143b628465cd3c8b
## 2026-03-04 18:37:49.619002 INFO::Fitting model to feature number 1250, X45a02d8c560f708f60e2eb952d45974d
## 2026-03-04 18:37:49.649525 INFO::Fitting model to feature number 1251, X820869c7c18af78f1bfe8f498e88e952
## 2026-03-04 18:37:49.680066 INFO::Fitting model to feature number 1252, dfaeccc139705834e879974f5b88c539
## 2026-03-04 18:37:49.710735 INFO::Fitting model to feature number 1253, X21cd87ba0ed75295c17dd698ef342af5
## 2026-03-04 18:37:49.74361 INFO::Fitting model to feature number 1254, X3c45cc52efc8e822effcb8d8fbec21d9
## 2026-03-04 18:37:49.774537 INFO::Fitting model to feature number 1255, X56e95bebfea3a564631df46966c413f8
## 2026-03-04 18:37:49.805393 INFO::Fitting model to feature number 1256, X1f12408f6bee875e10fcacbb5869bad0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:49.836022 INFO::Fitting model to feature number 1257, X5586fc56392c9262889f6ed54a2be971
## 2026-03-04 18:37:49.86652 INFO::Fitting model to feature number 1258, X7aaec5dc635a1d42d497edfe99e32c73
## 2026-03-04 18:37:49.897078 INFO::Fitting model to feature number 1259, X9b77571525020f2fb7a3d30e0161ce49
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:49.927494 INFO::Fitting model to feature number 1260, X9e3e8224c26b1c7c08dfa6e01de4ca5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:49.958358 INFO::Fitting model to feature number 1261, c886b40eb27052d8337856e81521ffbe
## 2026-03-04 18:37:49.988762 INFO::Fitting model to feature number 1262, e005c3a476aa46f8a7f918f704ef9cb6
## 2026-03-04 18:37:50.034564 INFO::Fitting model to feature number 1263, X13ca585ef9477f5bca479ffd2848097c
## 2026-03-04 18:37:50.067255 INFO::Fitting model to feature number 1264, X2242fe61584902366620d2a2096534f2
## 2026-03-04 18:37:50.098118 INFO::Fitting model to feature number 1265, X8f6c96e5604e4f77d84378b3d03782d0
## 2026-03-04 18:37:50.128681 INFO::Fitting model to feature number 1266, b71116fa177f4b07f53bf5d7996356a4
## 2026-03-04 18:37:50.158926 INFO::Fitting model to feature number 1267, e069f67ca919b8c627480de636066f4c
## 2026-03-04 18:37:50.18991 INFO::Fitting model to feature number 1268, X450cceb1cd3f129094c17ac29feb8d05
## 2026-03-04 18:37:50.220885 INFO::Fitting model to feature number 1269, X8ab70daaf7ddf3ee4baacee5f8c12693
## 2026-03-04 18:37:50.259115 INFO::Fitting model to feature number 1270, a0b24ff1cf7c2fc50822ba074982bfdd
## 2026-03-04 18:37:50.289731 INFO::Fitting model to feature number 1271, f11d3a89c1dae0dcdf1ff1e444155592
## 2026-03-04 18:37:50.320657 INFO::Fitting model to feature number 1272, e14d6eef3ea1efe7f99f0d1cfec2af99
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.352185 INFO::Fitting model to feature number 1273, X08030887b139d81f555f8f8352faf277
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.383015 INFO::Fitting model to feature number 1274, X1125d77abfbe71bd6698e9d969084307
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.414464 INFO::Fitting model to feature number 1275, X5862b017b7c3a0d098e0212f786ab5f0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.447203 INFO::Fitting model to feature number 1276, X7e3607fbbecab3ea931b4089a87ef58e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.478513 INFO::Fitting model to feature number 1277, c3027c24f295e7851649e1cd6e8b5013
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.509065 INFO::Fitting model to feature number 1278, e21e8275b16b5eb4781c76c5711b664c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.53936 INFO::Fitting model to feature number 1279, f222d3ac2aa36e3b93a61fcb909dcb28
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.57027 INFO::Fitting model to feature number 1280, X360335c2e90ac4b8274ba79192d02952
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.600768 INFO::Fitting model to feature number 1281, d0cf180c4dbb461db1db9e31f98e9017
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.631122 INFO::Fitting model to feature number 1282, b448327fe9304f5a67c5e8d7e63b47af
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.661883 INFO::Fitting model to feature number 1283, X265d49e10b1beb36a822e75f8499e466
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.698033 INFO::Fitting model to feature number 1284, c86f8eb8159950eceaddf95a4ab5b57c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.72972 INFO::Fitting model to feature number 1285, ad7793e14bb253aa0e97e853e8b015f2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.762768 INFO::Fitting model to feature number 1286, bcc404be7315530aea16b5a5f3d1c976
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.809871 INFO::Fitting model to feature number 1287, X5f5ee11493e3bb7a92450ba053af4591
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.84268 INFO::Fitting model to feature number 1288, X12b74af5dfbb936b04b4e8852792a241
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.875711 INFO::Fitting model to feature number 1289, X71e1a15d4596752eaeba8282f55eda63
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.907833 INFO::Fitting model to feature number 1290, a5638847d1da6845322b4c096586adf3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.93988 INFO::Fitting model to feature number 1291, X97e1d178e0060f5631143e0b6997db73
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:50.971437 INFO::Fitting model to feature number 1292, a32f6818da54b1c56a3f44425cbbccfd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.002372 INFO::Fitting model to feature number 1293, dcac11c3e60c3e677bb956a2e668718c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.03423 INFO::Fitting model to feature number 1294, f02787626d812f40886dff5a3b0e0fad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.066412 INFO::Fitting model to feature number 1295, X2554ec50a1739d46c1932d4b6d83bc24
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.098088 INFO::Fitting model to feature number 1296, X538b7a57f1012017b54549f162276e5c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.129612 INFO::Fitting model to feature number 1297, X7fe642a71440319d58bcb8d4b05ebe59
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.161708 INFO::Fitting model to feature number 1298, d6f3ce9ec26d83ee4d0c9615049f8056
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.193426 INFO::Fitting model to feature number 1299, f4cafac796fc9f62ae36a2571742ac56
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.225076 INFO::Fitting model to feature number 1300, X370b434407b901f334f8f70fcde243b8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.256505 INFO::Fitting model to feature number 1301, d74a8b7d8302fd484fce300c7f7cb90b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.28769 INFO::Fitting model to feature number 1302, X6138c8580647a54710e38072430b97b6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.318187 INFO::Fitting model to feature number 1303, X71dbe01e632b84071394a25599cc1b41
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.348401 INFO::Fitting model to feature number 1304, X86b78e8356db22087c189f4fb82f3cf5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.378913 INFO::Fitting model to feature number 1305, b6bf63a9790b48aa1c1dcf66baae1882
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.408945 INFO::Fitting model to feature number 1306, ae999de82bf29bf047a871241ccfbe50
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.448123 INFO::Fitting model to feature number 1307, c460b8c9a6aeac17873f4e562d3c631e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.47929 INFO::Fitting model to feature number 1308, d60fd4f4bf69a870d2f8db25d29f6486
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.510537 INFO::Fitting model to feature number 1309, dbadde65054417dedcb3c441117232e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.541578 INFO::Fitting model to feature number 1310, f44cf3ad8a5561421c5c86750b907568
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.583713 INFO::Fitting model to feature number 1311, X61ecc4bd7227a33c05851dfaef146197
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.619489 INFO::Fitting model to feature number 1312, f4c831a96e731d1e318dc1f4043c2f0b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.651838 INFO::Fitting model to feature number 1313, X067720f804b765e567ec36a0554015c3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.684493 INFO::Fitting model to feature number 1314, X2559f65fb87e56ab705a049961c3b596
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.715053 INFO::Fitting model to feature number 1315, e770005755e733e70992d93622110bac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.745838 INFO::Fitting model to feature number 1316, faa8d1a49c539195a242b95a3e78959a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.777539 INFO::Fitting model to feature number 1317, X0aba666ff3f58e07402855543abb0d77
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.809042 INFO::Fitting model to feature number 1318, X2ab2ed82aa11a2e91348a6323587cad0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.84059 INFO::Fitting model to feature number 1319, X7fce3acaf39b24f7c3eceaf8f4a509f1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.872194 INFO::Fitting model to feature number 1320, fc12d8847c55802fd660c17f96a4d239
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.909393 INFO::Fitting model to feature number 1321, X1237e3d469a445359f6344c4d7134df1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.941605 INFO::Fitting model to feature number 1322, be6eee56cb6557d3fde8a15a30c5f522
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:51.974284 INFO::Fitting model to feature number 1323, d978f3ca8b4c0a612f4a89fac040f433
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.012144 INFO::Fitting model to feature number 1324, X2d190136e77dd7bf9e1b271f6417bbbf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.048457 INFO::Fitting model to feature number 1325, X2eb915427713a21123065d4ac7d77481
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.080193 INFO::Fitting model to feature number 1326, df1d25de03331abaee085336549fccbf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.111606 INFO::Fitting model to feature number 1327, X2b33504d2a090ab8b7f93846f188fc84
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.149722 INFO::Fitting model to feature number 1328, X4514fcefdb82718d97b412184a893dc2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.19383 INFO::Fitting model to feature number 1329, X8605ed355fda390c4dbe1d350506b9bd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.227833 INFO::Fitting model to feature number 1330, e716a797ad8788f8c02701addcb9d49a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.258698 INFO::Fitting model to feature number 1331, f1c055f26db635e0551344b8e1c405ae
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.290498 INFO::Fitting model to feature number 1332, bf02fe466b9048dea2d3515ac8765402
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.321089 INFO::Fitting model to feature number 1333, d02e995ace556a85c3f1b34687268a29
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.351393 INFO::Fitting model to feature number 1334, X4289d6cbcecaffdc69d6ec33e18c1fcd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.382227 INFO::Fitting model to feature number 1335, X6dfdaa590e7618e0eaa870da258d29c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.434297 INFO::Fitting model to feature number 1336, X71faecc40533f09ce827e9ccd43f4d39
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.470382 INFO::Fitting model to feature number 1337, c5d7cdc0cd5fc338196e70876edc6839
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.501806 INFO::Fitting model to feature number 1338, X779d376dde67553423bf61868a8153ec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.535877 INFO::Fitting model to feature number 1339, c05488ee8eec8ca5ec58e73157e44cee
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.568933 INFO::Fitting model to feature number 1340, X1477502529efa4a2c43b1146a94e4a1f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.607093 INFO::Fitting model to feature number 1341, X50ec7bb8c7892c4a2d8d687da6beefef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.63795 INFO::Fitting model to feature number 1342, X8a9cb38167b3fcd87381c8985bc6819a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.675569 INFO::Fitting model to feature number 1343, X9856ac7a656cb6596ba6928448ded890
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.712438 INFO::Fitting model to feature number 1344, fbc2f9405994e3c4ca9c528d362faba1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.744702 INFO::Fitting model to feature number 1345, X2d4fb4e16e2ec5372a78961fcc33e48c
## 2026-03-04 18:37:52.779274 INFO::Fitting model to feature number 1346, X5d2b612c3a59a5919ab0c68e518a0921
## 2026-03-04 18:37:52.812741 INFO::Fitting model to feature number 1347, cac45687d42dfeaa1ad41cada5a86692
## 2026-03-04 18:37:52.846165 INFO::Fitting model to feature number 1348, X44d1c2442ea542df4a75d9eda9ffc130
## 2026-03-04 18:37:52.879931 INFO::Fitting model to feature number 1349, X9a9a3e95155b49f25bd224150719c489
## 2026-03-04 18:37:52.914878 INFO::Fitting model to feature number 1350, f7f661315c061fabf9c0d780d805a534
## 2026-03-04 18:37:52.949803 INFO::Fitting model to feature number 1351, d2dfb8391439482e10d58726a5e4fd8a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:52.991482 INFO::Fitting model to feature number 1352, X20e8adbd5fcb3bd31b9a6367c226e7db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.026875 INFO::Fitting model to feature number 1353, X9ecff22f4870ac57e0897f807ee34989
## 2026-03-04 18:37:53.06361 INFO::Fitting model to feature number 1354, d87e16aa9a2cd32ddf77082b4b182168
## 2026-03-04 18:37:53.097749 INFO::Fitting model to feature number 1355, X7aa1908721c25de6fbb0de8dde870b06
## 2026-03-04 18:37:53.138602 INFO::Fitting model to feature number 1356, X5ab5afd02358060e35e6b47024d2487f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.174325 INFO::Fitting model to feature number 1357, X8b955d9fba9b023835776d9b85cf1152
## 2026-03-04 18:37:53.215915 INFO::Fitting model to feature number 1358, ab718d064dc90341d0d9feab5795cede
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.251273 INFO::Fitting model to feature number 1359, X16cab7aa6e20994e74feb61b82b156a8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.302393 INFO::Fitting model to feature number 1360, d3e5fcec01188b4f3131fd3783b39618
## 2026-03-04 18:37:53.348254 INFO::Fitting model to feature number 1361, X8299fbf1188f45f3816d1b83d91a293c
## 2026-03-04 18:37:53.385073 INFO::Fitting model to feature number 1362, d52967973c245a507023695bc1681673
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.426444 INFO::Fitting model to feature number 1363, ff1991d802d8c09b2afedb3900b6677e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.462434 INFO::Fitting model to feature number 1364, X8831ebf2293100dea886fc93618e58ca
## 2026-03-04 18:37:53.498688 INFO::Fitting model to feature number 1365, X4bea2aea04a20f8b9ff601799f76f4a7
## 2026-03-04 18:37:53.534684 INFO::Fitting model to feature number 1366, b2a039758e7ec54d4dc52cb1385b3086
## 2026-03-04 18:37:53.576936 INFO::Fitting model to feature number 1367, X5b025d73c766db00d0bfc9384a36c699
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.612759 INFO::Fitting model to feature number 1368, d78f0ae1704d75a9a51b5d3203b46a18
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.65114 INFO::Fitting model to feature number 1369, X12711049f5f0b0fe67416dba7d01910b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.687433 INFO::Fitting model to feature number 1370, X2d04f5a36847a364b32939c0c7cc27de
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.724113 INFO::Fitting model to feature number 1371, X5da0eb7fa58b1b776bfe7a7c5f62bcdd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.761315 INFO::Fitting model to feature number 1372, X04ce44e19521cecbaa8d4f042e0f0154
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.797354 INFO::Fitting model to feature number 1373, X806807af7f9c8e247c7a5246ab1073ba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.83892 INFO::Fitting model to feature number 1374, f2df8d8cc211effb39df7e76235b0acb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.875066 INFO::Fitting model to feature number 1375, ddb238e065dcd145239865183f17b8a1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.911067 INFO::Fitting model to feature number 1376, X23b44691db0da7bcc193c4abd5e2994f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.946584 INFO::Fitting model to feature number 1377, X78bd9cd77859e67bf1cee19881c72185
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:53.98359 INFO::Fitting model to feature number 1378, fa53e063ccb39a71e3528a3246ab05e0
## 2026-03-04 18:37:54.017864 INFO::Fitting model to feature number 1379, X4e13eb363f4ee6b7744b0f7af637898e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.059891 INFO::Fitting model to feature number 1380, b0022a49e8ef1e2c089251ff260e88a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.094653 INFO::Fitting model to feature number 1381, X0761fbdedb674ca859a2e6464cd54682
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.129545 INFO::Fitting model to feature number 1382, f6d259ec62674092083cb3537a21ad15
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.164316 INFO::Fitting model to feature number 1383, dd233ad31d8d9dac73d102d3e32b503e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.199177 INFO::Fitting model to feature number 1384, X375132e5c4ba0f3c06a8b7a28c0a11d5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.258523 INFO::Fitting model to feature number 1385, X481f2e633138bf123a336f9c6809e5a5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.295168 INFO::Fitting model to feature number 1386, ae38dc32794bf8dbb622878b34e90fd8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.330064 INFO::Fitting model to feature number 1387, b845c8e9fe940562b0a61fd239b9f6a2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.363829 INFO::Fitting model to feature number 1388, X76d335ac41d3f074d0c27f5cfced7bd5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.396886 INFO::Fitting model to feature number 1389, X4d2f31f9e9a66db5837edbc35d6212ac
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.429722 INFO::Fitting model to feature number 1390, X6b644754f9e17055aa0fb066a1f39660
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.465527 INFO::Fitting model to feature number 1391, eab94e28ac20551638369c32c1890678
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.498522 INFO::Fitting model to feature number 1392, e4c6d550da4a82c067f9fe2d23ad7d87
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.532842 INFO::Fitting model to feature number 1393, X6d643d243d7154690d17a90f34214a15
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.566476 INFO::Fitting model to feature number 1394, X858ef2ff37b4306d93fdc9491470e18b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.60005 INFO::Fitting model to feature number 1395, a4bf35ce0ea6e33b1ca5fce93f2c72b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.633447 INFO::Fitting model to feature number 1396, a61c2f470b263001d328b36c16ee6a8f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.666481 INFO::Fitting model to feature number 1397, X947d4cca69b99f31480944eb79193ce7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.700035 INFO::Fitting model to feature number 1398, X271cb325b611cad7db3939e8e11c5541
## 2026-03-04 18:37:54.732121 INFO::Fitting model to feature number 1399, X823d7476979a0b8fc9479ae9113b7bf7
## 2026-03-04 18:37:54.766186 INFO::Fitting model to feature number 1400, X5ddfcafe5718c9c9f989fe6d946e741c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.806819 INFO::Fitting model to feature number 1401, X6cdde78777913ddf95c98995c3b4f62a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.841263 INFO::Fitting model to feature number 1402, X1bd40086437793e1ba33f5cbd83d6fe5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.880824 INFO::Fitting model to feature number 1403, X9abe6885f490b0b13978b8e993cf1e95
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.920821 INFO::Fitting model to feature number 1404, X316ab514dc5c338e49f89b0732fc5fca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.961074 INFO::Fitting model to feature number 1405, X648bee2fa5f6c2bee7662df9d276a906
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:54.997376 INFO::Fitting model to feature number 1406, X4b2e82cadba9d973d121de999f3199f9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.030869 INFO::Fitting model to feature number 1407, X21cd78cabd2ea1730146030a8ee5bb3a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.0659 INFO::Fitting model to feature number 1408, X5a2d5c1a794e467db1cf2588418be645
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.1216 INFO::Fitting model to feature number 1409, X169969cb9914e40725adc23728995292
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.162222 INFO::Fitting model to feature number 1410, c6a0cea6e8edca58bc207bfd93291981
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.198505 INFO::Fitting model to feature number 1411, fc1f59f65cb9dbb055d1b5e6c36ee5fb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.235053 INFO::Fitting model to feature number 1412, X34c29a3124e99390fe16cb7e3cf1542b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.267705 INFO::Fitting model to feature number 1413, X43690cf299862046d6c4f3c14cc26fbf
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.301019 INFO::Fitting model to feature number 1414, X6e151e6f762a3ddd31df62acda1fa509
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.339194 INFO::Fitting model to feature number 1415, X776f5ce8233c756ba7ae51611d4021c7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.372279 INFO::Fitting model to feature number 1416, X8123f26005a06511ef36459cd14c4be7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.405992 INFO::Fitting model to feature number 1417, X9498050fda592707a74ae42cb43cf376
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.439644 INFO::Fitting model to feature number 1418, b7824ad5efd327ec9aa70463ffe41313
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.474283 INFO::Fitting model to feature number 1419, be122441d4bfbf1b799a85d8765fb831
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:55.510496 INFO::Fitting model to feature number 1420, f522ba24c2cefcb9eaecc3bd442fc926
## 2026-03-04 18:37:55.553068 INFO::Fitting model to feature number 1421, X7d5650596629e727155c043cefe80839
## 2026-03-04 18:37:55.591387 INFO::Fitting model to feature number 1422, a57c1c3550d34aabb91ec1aaeccc3da0
## 2026-03-04 18:37:55.624553 INFO::Fitting model to feature number 1423, a6c825f96532401aa96ebbb3ca8cc4d5
## 2026-03-04 18:37:55.655537 INFO::Fitting model to feature number 1424, e6100b80a38222d7f4b53640a2813c92
## 2026-03-04 18:37:55.686935 INFO::Fitting model to feature number 1425, X5e298cb51b192a2586c265bd13ecddec
## 2026-03-04 18:37:55.718146 INFO::Fitting model to feature number 1426, d6bbba5be34ced0b2ddcf6d419a290fd
## 2026-03-04 18:37:55.750736 INFO::Fitting model to feature number 1427, X1e5afa7c69a7d5bc728c0b1c178825a4
## 2026-03-04 18:37:55.782424 INFO::Fitting model to feature number 1428, a72f4146638b9975a2f70a429c7888aa
## 2026-03-04 18:37:55.813907 INFO::Fitting model to feature number 1429, a6dd7f11774533a463aa01d6c837dc00
## 2026-03-04 18:37:55.845556 INFO::Fitting model to feature number 1430, X394106ad4bcc01e28e654f8b904a43f0
## 2026-03-04 18:37:55.876896 INFO::Fitting model to feature number 1431, X3323196417bad3bec4e7a302d5dfdd61
## 2026-03-04 18:37:55.907216 INFO::Fitting model to feature number 1432, X568794df14ca80a20df326e942094f60
## 2026-03-04 18:37:55.938605 INFO::Fitting model to feature number 1433, dc70f2dde3a4e16635013a39d5bc4738
## 2026-03-04 18:37:55.984523 INFO::Fitting model to feature number 1434, X598d0ff1000e5b59358afb0c03dd8f91
## 2026-03-04 18:37:56.016947 INFO::Fitting model to feature number 1435, X48d757c9fc57146957b1f61c94e1cc19
## 2026-03-04 18:37:56.04865 INFO::Fitting model to feature number 1436, c7c955897cbfe6837c99e7cd44d0887e
## 2026-03-04 18:37:56.079576 INFO::Fitting model to feature number 1437, e369521c3e10197580daa61efd80c436
## 2026-03-04 18:37:56.110585 INFO::Fitting model to feature number 1438, f7e47a7f64694d3180b261250205ce4a
## 2026-03-04 18:37:56.141075 INFO::Fitting model to feature number 1439, X3481a13605bd74eb2f7cc589b20486b5
## 2026-03-04 18:37:56.171801 INFO::Fitting model to feature number 1440, X7e75bd65b31ecb2641e7317dac525229
## 2026-03-04 18:37:56.202358 INFO::Fitting model to feature number 1441, b70927f642e4add511f7bd662c667a58
## 2026-03-04 18:37:56.233273 INFO::Fitting model to feature number 1442, X224b3eccb629d9e4a636f3af73fb0268
## 2026-03-04 18:37:56.265128 INFO::Fitting model to feature number 1443, eac0c8b1f85bfc7dba195c75ea8f3561
## 2026-03-04 18:37:56.296898 INFO::Fitting model to feature number 1444, X56c30e736c54be6e83da0343c057e762
## 2026-03-04 18:37:56.327987 INFO::Fitting model to feature number 1445, X9d1c840bc95d7a1d599e26ec32d6be01
## 2026-03-04 18:37:56.359092 INFO::Fitting model to feature number 1446, X49f668308b55c43b37f983fc8bc7ce4b
## 2026-03-04 18:37:56.38993 INFO::Fitting model to feature number 1447, X8fa3bb6a74294c2f87ffa923221aaef1
## 2026-03-04 18:37:56.420729 INFO::Fitting model to feature number 1448, X7c54d19b47ce29e4e80f24ec6e8abcfc
## 2026-03-04 18:37:56.457174 INFO::Fitting model to feature number 1449, X4202fe0cef4d4c2bf106ba2027a23bb8
## 2026-03-04 18:37:56.488838 INFO::Fitting model to feature number 1450, X22bf5bfb7975b9c37934c7334e313408
## 2026-03-04 18:37:56.522047 INFO::Fitting model to feature number 1451, X388c2c1d2035cce97331699e6dcf29b4
## 2026-03-04 18:37:56.55478 INFO::Fitting model to feature number 1452, X439a7765de6fa0674816e837c3810090
## 2026-03-04 18:37:56.586418 INFO::Fitting model to feature number 1453, bd536adea9f0edde6e5eafb1da9be283
## 2026-03-04 18:37:56.618097 INFO::Fitting model to feature number 1454, e3690c2df119bc765abb284d23301e5c
## 2026-03-04 18:37:56.653455 INFO::Fitting model to feature number 1455, X86e184a98df6da6c11e688e0e794d017
## 2026-03-04 18:37:56.685934 INFO::Fitting model to feature number 1456, X20c0706bcdf2adf7708e01585bf70288
## 2026-03-04 18:37:56.724045 INFO::Fitting model to feature number 1457, X2414414f9af0314343395f48bca12ec7
## 2026-03-04 18:37:56.769944 INFO::Fitting model to feature number 1458, X0616a2d8a3a5c92df8a0acd7c90be22a
## 2026-03-04 18:37:56.80763 INFO::Fitting model to feature number 1459, X03d7aab491a1a55630aa6c485fafe72c
## 2026-03-04 18:37:56.8426 INFO::Fitting model to feature number 1460, b0031e259b262ad8ac0ff64393f013fa
## 2026-03-04 18:37:56.872981 INFO::Fitting model to feature number 1461, X46ef44d3a2f2d113352a55c8d8ffad04
## 2026-03-04 18:37:56.902809 INFO::Fitting model to feature number 1462, b713d99eb533b299629c1394f75b6926
## 2026-03-04 18:37:56.93438 INFO::Fitting model to feature number 1463, X14dbca841e6bb9568be475a3360d5a58
## 2026-03-04 18:37:56.964931 INFO::Fitting model to feature number 1464, X4e4c5f7c2c776156be82a7e1d24cdffb
## 2026-03-04 18:37:56.996333 INFO::Fitting model to feature number 1465, X34c3988ed3a5e3630b8488505b7f9f98
## 2026-03-04 18:37:57.028111 INFO::Fitting model to feature number 1466, X72ee97f042fd5f6508c8f55419c8862e
## 2026-03-04 18:37:57.060575 INFO::Fitting model to feature number 1467, X7b3fe7dea742fd319815e05d8c5f9fe1
## 2026-03-04 18:37:57.091912 INFO::Fitting model to feature number 1468, X5006a0d2b1a915d327889a08c9e8885f
## 2026-03-04 18:37:57.123637 INFO::Fitting model to feature number 1469, a4bd3b3c5909749d8f08456f744defa0
## 2026-03-04 18:37:57.15422 INFO::Fitting model to feature number 1470, bf2eac04920ee7f7f1d299a292613df2
## 2026-03-04 18:37:57.184306 INFO::Fitting model to feature number 1471, e79c3000d572833b5371f40ba4a529ea
## 2026-03-04 18:37:57.214456 INFO::Fitting model to feature number 1472, X539fab4811fbd687bbeca6e06a832e8e
## 2026-03-04 18:37:57.244839 INFO::Fitting model to feature number 1473, ca1c86768a0a1d127988348bb1a84a4e
## 2026-03-04 18:37:57.2749 INFO::Fitting model to feature number 1474, X5f9c85a0e24285e0d751487138070fcb
## 2026-03-04 18:37:57.305066 INFO::Fitting model to feature number 1475, X0432cc28a194dcd588a5a751fac80832
## 2026-03-04 18:37:57.338178 INFO::Fitting model to feature number 1476, a75f41af38c6efe6ff870067af0eecd5
## 2026-03-04 18:37:57.369692 INFO::Fitting model to feature number 1477, X5712fb6a5e9c90125f5029a8cc18f538
## 2026-03-04 18:37:57.400303 INFO::Fitting model to feature number 1478, e35b85bcf133e5c1406ade391f27f35e
## 2026-03-04 18:37:57.437626 INFO::Fitting model to feature number 1479, X377d5ee705688586a17936fb3804994a
## 2026-03-04 18:37:57.475123 INFO::Fitting model to feature number 1480, X5b73b77b0b395c24b331446c95c75533
## 2026-03-04 18:37:57.508337 INFO::Fitting model to feature number 1481, de4b3ada2cc38daf9fcdc0741de5dca7
## 2026-03-04 18:37:57.549565 INFO::Fitting model to feature number 1482, X9d9c864dccdebdc7a243f28d6510f2ef
## 2026-03-04 18:37:57.58698 INFO::Fitting model to feature number 1483, X47cf08ac4ef1b9bbfb862292206a9f7c
## 2026-03-04 18:37:57.620195 INFO::Fitting model to feature number 1484, X531af51c5105e6170ef66e92168256c1
## 2026-03-04 18:37:57.651413 INFO::Fitting model to feature number 1485, acd2936ee50c99be98115e15414cb216
## 2026-03-04 18:37:57.682677 INFO::Fitting model to feature number 1486, c0484a4bafa44a0dd28c9d0880a3d591
## 2026-03-04 18:37:57.713062 INFO::Fitting model to feature number 1487, X211046170b10ef39b5c78eae5a7249d3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:57.743769 INFO::Fitting model to feature number 1488, X75c94a2e7f736b8dd4894c788c15ead7
## 2026-03-04 18:37:57.774899 INFO::Fitting model to feature number 1489, b2c82de8806b29b0f4f22aff7ee128c6
## 2026-03-04 18:37:57.805757 INFO::Fitting model to feature number 1490, X44cbcbce5883cdcee679699310064855
## 2026-03-04 18:37:57.836517 INFO::Fitting model to feature number 1491, X28d4648bb079932a475e8b4e2cb7b203
## 2026-03-04 18:37:57.867511 INFO::Fitting model to feature number 1492, X73ed71d087d17b73df11c5ef4d509a48
## 2026-03-04 18:37:57.899412 INFO::Fitting model to feature number 1493, aee63dacdb65a4e383d739b9a823c6d5
## 2026-03-04 18:37:57.930745 INFO::Fitting model to feature number 1494, X02d8636ae9a2d97764ca155d429a3573
## 2026-03-04 18:37:57.962113 INFO::Fitting model to feature number 1495, X7c2c90b4dd80a6c5e232c200e12bb135
## 2026-03-04 18:37:57.992657 INFO::Fitting model to feature number 1496, X0de79de542d618ec0a3640507d175b58
## 2026-03-04 18:37:58.022933 INFO::Fitting model to feature number 1497, b6d0a5d1db20bb06207235f7ca1c80e6
## 2026-03-04 18:37:58.053998 INFO::Fitting model to feature number 1498, X61bedb706ed0ebf66136c2de5e894ae0
## 2026-03-04 18:37:58.084779 INFO::Fitting model to feature number 1499, X85e51deb7d526a41feb684d724a249f7
## 2026-03-04 18:37:58.11483 INFO::Fitting model to feature number 1500, X50e3f160c939685dd9fe2516743de428
## 2026-03-04 18:37:58.144945 INFO::Fitting model to feature number 1501, X9be98a384b63040c28b7ba80f5bcfaa6
## 2026-03-04 18:37:58.174909 INFO::Fitting model to feature number 1502, bb4981f12957c9153404475002597fcf
## 2026-03-04 18:37:58.205184 INFO::Fitting model to feature number 1503, X55b6147730ff0b48500535882c1c90c8
## 2026-03-04 18:37:58.235475 INFO::Fitting model to feature number 1504, X70cec6470fc234fda2298dd42df999c8
## 2026-03-04 18:37:58.266108 INFO::Fitting model to feature number 1505, ab7683701ff9cd914c709a432024c95a
## 2026-03-04 18:37:58.298145 INFO::Fitting model to feature number 1506, ebae2a99fbe85d6f7192a78fccdb8fdb
## 2026-03-04 18:37:58.344368 INFO::Fitting model to feature number 1507, X82949f67dacd7431521ab84ca247a5d3
## 2026-03-04 18:37:58.376443 INFO::Fitting model to feature number 1508, a2a7e83f01568ca2ffc7aaeb75be4515
## 2026-03-04 18:37:58.407773 INFO::Fitting model to feature number 1509, X7cb5a8b160490076fb6549561c3423f1
## 2026-03-04 18:37:58.438734 INFO::Fitting model to feature number 1510, X003efbd0f8008eb1a10db7323ae5b9cd
## 2026-03-04 18:37:58.469949 INFO::Fitting model to feature number 1511, X658685f5b6770693b0e2a9758ffee903
## 2026-03-04 18:37:58.501455 INFO::Fitting model to feature number 1512, X124e851f1317df1c44f69576ded77eea
## 2026-03-04 18:37:58.534511 INFO::Fitting model to feature number 1513, X247c6962f97e1327653c49a55baaf0ef
## 2026-03-04 18:37:58.565839 INFO::Fitting model to feature number 1514, X95acc8353370e4c6473aabd257f89145
## 2026-03-04 18:37:58.596414 INFO::Fitting model to feature number 1515, a6be673d83361f0db6b04f4f6368aac8
## 2026-03-04 18:37:58.626651 INFO::Fitting model to feature number 1516, e78522b1e525cf73de0a832161386cfd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:58.667064 INFO::Fitting model to feature number 1517, X40ead826bf0cabf57669891968fb8028
## 2026-03-04 18:37:58.697844 INFO::Fitting model to feature number 1518, X4eafbc59f43c8a0357c8c6c01196f799
## 2026-03-04 18:37:58.728225 INFO::Fitting model to feature number 1519, b69359c339f22ecb5db0d3f40004badc
## 2026-03-04 18:37:58.75861 INFO::Fitting model to feature number 1520, X462a59f090f13da057d0e2bcf7113ffc
## 2026-03-04 18:37:58.790652 INFO::Fitting model to feature number 1521, df8cea54149f075a0cf136ebd86633f1
## 2026-03-04 18:37:58.821014 INFO::Fitting model to feature number 1522, fb647843e3c65da0b4fcddbc5afa2900
## 2026-03-04 18:37:58.851444 INFO::Fitting model to feature number 1523, ff1b4e1105cab28d1b9274d9aad00686
## 2026-03-04 18:37:58.889803 INFO::Fitting model to feature number 1524, X350fce626781822ad4f661f8ab50950e
## 2026-03-04 18:37:58.920128 INFO::Fitting model to feature number 1525, X516edf87db3e79ac6cb028e4866128c4
## 2026-03-04 18:37:58.950537 INFO::Fitting model to feature number 1526, c6fba20dabb5158b66a8c275f1972c9c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:58.991866 INFO::Fitting model to feature number 1527, cea3551436db23c969c3846aaca74006
## 2026-03-04 18:37:59.022759 INFO::Fitting model to feature number 1528, X40ee6ad12f5a84a745007123217b2f14
## 2026-03-04 18:37:59.054469 INFO::Fitting model to feature number 1529, X41bfe90c919674be43e40cc1ed15f397
## 2026-03-04 18:37:59.085787 INFO::Fitting model to feature number 1530, X5888c34571979386d88d40e153020962
## 2026-03-04 18:37:59.137165 INFO::Fitting model to feature number 1531, X9d1c414e4e43f9035786f41fe51a2c37
## 2026-03-04 18:37:59.171515 INFO::Fitting model to feature number 1532, X2ed61e783e248ef47514084aba18efe5
## 2026-03-04 18:37:59.203787 INFO::Fitting model to feature number 1533, X85eb6d9d8178170857155c18c4747e83
## 2026-03-04 18:37:59.234467 INFO::Fitting model to feature number 1534, X192171db9d9d7eef504462e639aff14e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:59.265288 INFO::Fitting model to feature number 1535, X382cb4846b04b8a16502cbccaa8fd78b
## 2026-03-04 18:37:59.295507 INFO::Fitting model to feature number 1536, X3bc397e5ed09b0d3152b2715837be559
## 2026-03-04 18:37:59.326822 INFO::Fitting model to feature number 1537, X45d5c0bb1e1bb6300d395413db2070c8
## 2026-03-04 18:37:59.357907 INFO::Fitting model to feature number 1538, X712381dc322c52191a6601dfa1ea3d5b
## 2026-03-04 18:37:59.387459 INFO::Fitting model to feature number 1539, X8a5fcf814785cd96fc7be2582c438c7c
## 2026-03-04 18:37:59.417932 INFO::Fitting model to feature number 1540, ae2ef6affcdf43d6f0f516c8212da640
## 2026-03-04 18:37:59.451965 INFO::Fitting model to feature number 1541, X17122257dd2d67ac0e4167bd03651e51
## 2026-03-04 18:37:59.48371 INFO::Fitting model to feature number 1542, X30d593c68137f5c78b3750e68af4b9f7
## 2026-03-04 18:37:59.513911 INFO::Fitting model to feature number 1543, X4540e89c3c1314208895f374fe7c9a6a
## 2026-03-04 18:37:59.545672 INFO::Fitting model to feature number 1544, X6277a075fc7865b1e27d287c43a7e57c
## 2026-03-04 18:37:59.576911 INFO::Fitting model to feature number 1545, X7a7377a0454382bbad229ff7c5d80546
## 2026-03-04 18:37:59.607572 INFO::Fitting model to feature number 1546, X7c3929432548e1d0644f5b531b396ad0
## 2026-03-04 18:37:59.63889 INFO::Fitting model to feature number 1547, b40189cb15dd5644b0adea1ff4afb315
## 2026-03-04 18:37:59.670385 INFO::Fitting model to feature number 1548, X95d765a3b3dfd08fdb458a3c6a7054b9
## 2026-03-04 18:37:59.700811 INFO::Fitting model to feature number 1549, e8019eaaf310c14bc4b7be16e7d3bfa1
## 2026-03-04 18:37:59.731766 INFO::Fitting model to feature number 1550, X7fbb704ff0cbacc0253188c453d69e36
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:37:59.764301 INFO::Fitting model to feature number 1551, X9be5e759a3fdd983eedea94bfca182f9
## 2026-03-04 18:37:59.795379 INFO::Fitting model to feature number 1552, X13bc51f05823bc48301a5d1c4bde9f5a
## 2026-03-04 18:37:59.833352 INFO::Fitting model to feature number 1553, X3e82bcc53ad83cc493c60475fbcc727d
## 2026-03-04 18:37:59.864294 INFO::Fitting model to feature number 1554, X632aa2401824dfbff3dccab0c136de74
## 2026-03-04 18:37:59.89474 INFO::Fitting model to feature number 1555, X6fe547e6fe340a40c701e84a3d8eb1cb
## 2026-03-04 18:37:59.939042 INFO::Fitting model to feature number 1556, X8387bba230a7f9e754d83c9171c28b31
## 2026-03-04 18:37:59.970772 INFO::Fitting model to feature number 1557, a0a9b9c3660fa6ae48f562a494208215
## 2026-03-04 18:38:00.001741 INFO::Fitting model to feature number 1558, X56afd101747172bfbeae0ac912dfc469
## 2026-03-04 18:38:00.033969 INFO::Fitting model to feature number 1559, X84bfbbf1dcab3c3f73b0089d3bb32c7c
## 2026-03-04 18:38:00.066819 INFO::Fitting model to feature number 1560, X42d802e1f056d33a6fc7e8db5a5bbff7
## 2026-03-04 18:38:00.098479 INFO::Fitting model to feature number 1561, X6e330d70c19eb8bf817a65230895f026
## 2026-03-04 18:38:00.129896 INFO::Fitting model to feature number 1562, a1b0ba01f038cc733b4bfe022a6525dc
## 2026-03-04 18:38:00.161249 INFO::Fitting model to feature number 1563, bc8fd2265d620b8bb9eea8d62bf623c6
## 2026-03-04 18:38:00.192493 INFO::Fitting model to feature number 1564, X5091bd9b1943a71c946ae11b4e07b7e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.223958 INFO::Fitting model to feature number 1565, X5d3a339376909059f6ad187fdc28d6df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.255715 INFO::Fitting model to feature number 1566, X5deed8fedf1ec0ed7efc98743a7dbb25
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.286529 INFO::Fitting model to feature number 1567, d77a370884869d02025c707a1059ebf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.317487 INFO::Fitting model to feature number 1568, X96535061b17238b675b7c8397d18b1db
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.351076 INFO::Fitting model to feature number 1569, c94467ba9b103e3de9d85eaca888a027
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.38811 INFO::Fitting model to feature number 1570, X0b7c20e95c7d503306b44bc853c0f948
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.423013 INFO::Fitting model to feature number 1571, X589a69391f6c119b48b9884552436536
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.45511 INFO::Fitting model to feature number 1572, X922ebffab737973fc23feaefe8df04e8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.485775 INFO::Fitting model to feature number 1573, X183006dbbfa96d0df78d911109617c04
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.516158 INFO::Fitting model to feature number 1574, X47add6d359e281d9307ef7b818b1132c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.548607 INFO::Fitting model to feature number 1575, accdc77f3b67a008ff30692229aa28ef
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.580052 INFO::Fitting model to feature number 1576, X76215b5328d849f5bb5aa09fee01861f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.611792 INFO::Fitting model to feature number 1577, fe585ab842d898d8e35d20f9005df2d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.643635 INFO::Fitting model to feature number 1578, X52c94f53dc5b8610cb0a96ec22b853eb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.676381 INFO::Fitting model to feature number 1579, X7c5b049887df626ca87a19fbbe1609d1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.722425 INFO::Fitting model to feature number 1580, e8bb3d4f052a5add7027c46cbaead30b
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.754972 INFO::Fitting model to feature number 1581, X5708f680f66bc7b8047b515f2b4c4bc3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.78773 INFO::Fitting model to feature number 1582, X1948dc0f10daf853102b0731a6186af4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.82016 INFO::Fitting model to feature number 1583, X43286dda280419b7f392dd6ce424759a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.852078 INFO::Fitting model to feature number 1584, ef67f9a0ad8f89623991c3f3d2ec7d16
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.883876 INFO::Fitting model to feature number 1585, bd17fc6e6fee81cc3e9ec0905eb366b3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.915672 INFO::Fitting model to feature number 1586, X5c0401ad6b1ee6bb92b542bd10210f51
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.94756 INFO::Fitting model to feature number 1587, X4807347adae205610568f73ccc96c6a7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:00.979242 INFO::Fitting model to feature number 1588, f4bc6194719dec9e61ae79c9b52fcbf3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.011703 INFO::Fitting model to feature number 1589, X1b76fcb2cecd29cdfff9a0af888d4169
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.050312 INFO::Fitting model to feature number 1590, X6158e79837bbf56a578b6d8b44ef7adb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.088018 INFO::Fitting model to feature number 1591, d71613d583910ae8267587f2b140ceb7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.118879 INFO::Fitting model to feature number 1592, X8ceade32cdc05781325d8ad50f31f2b0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.150307 INFO::Fitting model to feature number 1593, X58c9c50389bbfe8ad20f1d089a45172d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.181297 INFO::Fitting model to feature number 1594, X9f4c8ff7db4a28f9d8a9d36c2c773eba
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.211934 INFO::Fitting model to feature number 1595, eb740349dd1e6ad53ae89a46ec50ecb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.242967 INFO::Fitting model to feature number 1596, X915e1c30b64ef9a338d6c10600e4a617
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.274364 INFO::Fitting model to feature number 1597, f027bfc4feea257cddbf035993ba545a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.304487 INFO::Fitting model to feature number 1598, X6079966a4d2ff4c7cdc74e1741c30fd6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.335442 INFO::Fitting model to feature number 1599, c6ce5406ec57d4e6377b8c641e8975da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.367408 INFO::Fitting model to feature number 1600, X56c5746a26ab97c041720b362ad8661e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.39806 INFO::Fitting model to feature number 1601, X6baf37ed88a786881d547a08161907e4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.430457 INFO::Fitting model to feature number 1602, a34ec679cd800d55496935d18635358d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.46319 INFO::Fitting model to feature number 1603, cc01d85478bd68c5a0abdb6b92a6bb57
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.960881 INFO::Fitting model to feature number 1604, d10cef946e0c66aff6e1b5b7976792ad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:01.992876 INFO::Fitting model to feature number 1605, X1596a42ede5440896d482a43e80cf4be
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.023875 INFO::Fitting model to feature number 1606, X6174b3382fc068d14dd5ff3e1b5f7657
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.055042 INFO::Fitting model to feature number 1607, X80ab0c108d89e0bea0795b22fdd78c45
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.092607 INFO::Fitting model to feature number 1608, X84ad86e97445e7787a94091968391d2d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.129711 INFO::Fitting model to feature number 1609, bc148193d702e654dd5beafc48ed3b3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.16054 INFO::Fitting model to feature number 1610, X5819fc035f12f49f114cdc42ebd7f872
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.19661 INFO::Fitting model to feature number 1611, c1ad43ca73c671841937ae7a7f12bd00
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.234464 INFO::Fitting model to feature number 1612, X3857b4ad66f5288bdfaee996a22c2dd5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.271447 INFO::Fitting model to feature number 1613, X5ec337477fdc4a209b4d08ba19627c12
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.306817 INFO::Fitting model to feature number 1614, X88735c01a26d318d22ab194b70d82670
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.33948 INFO::Fitting model to feature number 1615, X993b533f1f2a19c778a99ee2d73d6401
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.37621 INFO::Fitting model to feature number 1616, c50ece6688d9455f2089d10c593dc5e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.409751 INFO::Fitting model to feature number 1617, X1f07789f2a30d8406f4e1e6c4e59d512
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.451282 INFO::Fitting model to feature number 1618, X12957d1b9a2b84263f7c2577a0be5ffe
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.487772 INFO::Fitting model to feature number 1619, X8ee9e03c0b6df43238a2bac06e8157c6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.5204 INFO::Fitting model to feature number 1620, X7c1a951cbb503325aa9c776b49db8e62
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.555097 INFO::Fitting model to feature number 1621, X37af96a311f2b6d6d9e6f1aeac8adcbb
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.588432 INFO::Fitting model to feature number 1622, aa0d9d7fa21622e919bc657a1c27acea
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.621372 INFO::Fitting model to feature number 1623, a8d44a386b74d5efeb2ab88dd7b14d81
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.654877 INFO::Fitting model to feature number 1624, X6e2c98b75514f3cec57b18d85e847d56
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.694111 INFO::Fitting model to feature number 1625, b6e3cadd45df3d2e2e465e4b2e49fbf6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.72672 INFO::Fitting model to feature number 1626, X8e714ede8d3de600e6e1e85b93ffff3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.759705 INFO::Fitting model to feature number 1627, e36960ffa2d5b1880b964b90d5e01ad1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.79297 INFO::Fitting model to feature number 1628, X1012dfedbd7005d166522a23a70e17e2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.833491 INFO::Fitting model to feature number 1629, ee8bdf48063e7dee314d7106122de356
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.868507 INFO::Fitting model to feature number 1630, X73741e0f0fc53215dfc000d6bcd03bf2
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.905459 INFO::Fitting model to feature number 1631, b4f013a316dd9db71487186d4e5f8c48
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.942161 INFO::Fitting model to feature number 1632, d80f4c95bba89e950a462cde932b4a19
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:02.978323 INFO::Fitting model to feature number 1633, ec281d5415aa63e552133483f7b7f263
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.009043 INFO::Fitting model to feature number 1634, X3889da5feeb2f015f4e43df76d1ba1c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.040404 INFO::Fitting model to feature number 1635, X0fdfc03c69957857ea0aa99afb974be6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.077699 INFO::Fitting model to feature number 1636, bb71589a895e6aac9ac57ef82e01389f
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.108799 INFO::Fitting model to feature number 1637, X6f8c98e3a54279e65fb5b3d6c4596b06
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.148589 INFO::Fitting model to feature number 1638, ef7d9398b53a4c43e226717ea5256aa8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.179261 INFO::Fitting model to feature number 1639, X0744af08f13a84f475cacdfe0a7eb563
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.210978 INFO::Fitting model to feature number 1640, b933d7cfe0a4a8671012717ad8f3d467
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.24314 INFO::Fitting model to feature number 1641, X8e898fac58bba1bbf6b1d86cda3d468a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.275781 INFO::Fitting model to feature number 1642, f90ae11968e525ce0b184a7d096010b5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.306822 INFO::Fitting model to feature number 1643, X6149b30955bdf96ec6892f0d385ddeb4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.33773 INFO::Fitting model to feature number 1644, X878ab77c15d89ce76f14512b57b66554
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.370032 INFO::Fitting model to feature number 1645, X97fb8735d375018eef0bc292ced80c3c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.400269 INFO::Fitting model to feature number 1646, X122ae66707c8158bbd1a25e11191df1e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.432273 INFO::Fitting model to feature number 1647, X08d50335ad948012223e0bcb0e02a98c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.471427 INFO::Fitting model to feature number 1648, X65450912952dfde7e9482ff90cb96c14
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.506855 INFO::Fitting model to feature number 1649, X9a11808482e84a5875187726140985df
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.544862 INFO::Fitting model to feature number 1650, bcab38bdf56dea2067767829a85b6081
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.577249 INFO::Fitting model to feature number 1651, X70aa06947222e61a0342014baad5b4fd
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.609358 INFO::Fitting model to feature number 1652, X75df890f90e4d0da4c0bad5bfe8344c8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.645918 INFO::Fitting model to feature number 1653, X36a5741f495aca0b66ad97f97dc7a194
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.691925 INFO::Fitting model to feature number 1654, X10415149c4af2e9597ec437422950ed1
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.722931 INFO::Fitting model to feature number 1655, X57f59e0ae6b1c8ce77629f2b8cb9b6ec
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.756172 INFO::Fitting model to feature number 1656, X839baf767eec4eed3854e3ac4ff072c5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.788318 INFO::Fitting model to feature number 1657, X2b799a142091659df01cda80e96b81ca
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.820981 INFO::Fitting model to feature number 1658, X947eb89bd4aa5e20d8bf953ea6c0d081
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.852075 INFO::Fitting model to feature number 1659, cc61a14c8cb9f6513a9103c1a21d0b55
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.889641 INFO::Fitting model to feature number 1660, c65985845943a827744fa6d3de55535d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.920319 INFO::Fitting model to feature number 1661, X88117a1cd850c7fafeeef4fa07f31d41
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.95735 INFO::Fitting model to feature number 1662, bedc4e430c931cfaab616b853069f814
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:03.995133 INFO::Fitting model to feature number 1663, X54559e1264fb8ce1cc649277c8326fad
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.026977 INFO::Fitting model to feature number 1664, X10c2b175dcb8d7348ee54da806d75dde
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.05991 INFO::Fitting model to feature number 1665, X5774e92983913d64432d9cee33e6fd3e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.091986 INFO::Fitting model to feature number 1666, X8f6915518fcb472c94ac31408b61a0da
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.135739 INFO::Fitting model to feature number 1667, d74a6273b1480e0d17fe444efd8c855d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.167428 INFO::Fitting model to feature number 1668, ebb85b3d0073e8567c7187e4d074f323
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.198422 INFO::Fitting model to feature number 1669, ff0b4d81a135df7e6b677c1e0e371d26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.229328 INFO::Fitting model to feature number 1670, c7f5deefec0090691af7e3564a82429d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.260594 INFO::Fitting model to feature number 1671, d624453d3dcd59833f711d70d8be5e2a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.290926 INFO::Fitting model to feature number 1672, X6a7065ab25b547f46bbb229e5ed5a8d9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.321015 INFO::Fitting model to feature number 1673, X7b7686f9f714d532db91900814a6ffe9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.352889 INFO::Fitting model to feature number 1674, X858cadaa6bed41000fb485d316ba52a4
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.390415 INFO::Fitting model to feature number 1675, eb56c11fff2b5d82b6918720792b3902
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.420239 INFO::Fitting model to feature number 1676, a693465c86c54e2943f05ad39652a938
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.453954 INFO::Fitting model to feature number 1677, f197afa977ca78606df5de118bbe0bb3
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.487163 INFO::Fitting model to feature number 1678, X9e1f8fa2e00d442b56a9262b77fdd34c
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.52893 INFO::Fitting model to feature number 1679, X466af77afd608fc85dfc61d3a08010fa
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.559751 INFO::Fitting model to feature number 1680, X615e8c85484fde9ca7fc10a9ca8d2aa5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.592434 INFO::Fitting model to feature number 1681, b34dcecfc1eb9380c919f8566479a79e
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.625212 INFO::Fitting model to feature number 1682, X459037d1dd08629797b49e758ec86da7
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.656629 INFO::Fitting model to feature number 1683, X30d07efaf39113af71b777723b986086
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.693812 INFO::Fitting model to feature number 1684, cf22fbcb43140202c9d69acac1d2d7d6
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.726772 INFO::Fitting model to feature number 1685, X708f08717de48853cf9afe00b1198a27
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.759841 INFO::Fitting model to feature number 1686, a1a62a39c30efb6136daa1c51cba279a
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.79058 INFO::Fitting model to feature number 1687, a48690290fc79a4f0d211eb587077fa8
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.828759 INFO::Fitting model to feature number 1688, f67744afabc94b4f4c0a3b49ac8d862d
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.866302 INFO::Fitting model to feature number 1689, X05269a738f261937bd7740c6501454e5
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.904526 INFO::Fitting model to feature number 1690, X2eb6747a7aebca6666a096a37071d7e9
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.936806 INFO::Fitting model to feature number 1691, X3db4777463dd9b345036bcdb868ef547
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:04.969559 INFO::Fitting model to feature number 1692, X53b62dd9d227d3bc7297f0dbd48a4989
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:05.000103 INFO::Fitting model to feature number 1693, X628547bdac555ee7526d0216152bcf26
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:05.032302 INFO::Fitting model to feature number 1694, X6f1fbf83f07a443433b7f1e88a3dd392
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:05.064895 INFO::Fitting model to feature number 1695, X9bb4fc212f3b982613792c538ea52694
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:05.094944 INFO::Fitting model to feature number 1696, a0d57e53d8bc50f10b7857f69a8c5fa0
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 2026-03-04 18:38:05.20706 INFO::Counting total values for each feature
## 2026-03-04 18:38:05.273379 INFO::Writing filtered data to file ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/features/filtered_data.tsv
## 2026-03-04 18:38:05.354021 INFO::Writing filtered, normalized data to file ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/features/filtered_data_norm.tsv
## 2026-03-04 18:38:05.46105 INFO::Writing filtered, normalized, transformed data to file ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/features/filtered_data_norm_transformed.tsv
## 2026-03-04 18:38:05.704345 INFO::Writing residuals to file ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/fits/residuals.rds
## 2026-03-04 18:38:05.729396 INFO::Writing fitted values to file ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/fits/fitted.rds
## 2026-03-04 18:38:05.748605 INFO::Writing extracted random effects to file ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/fits/ranef.rds
## 2026-03-04 18:38:05.76059 INFO::Writing all results to file (ordered by increasing q-values): ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/all_results.tsv
## 2026-03-04 18:38:05.805105 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.050000 ) to file (ordered by increasing q-values): ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/significant_results.tsv
## 2026-03-04 18:38:05.812014 INFO::Writing heatmap of significant results to file: ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations/heatmap.pdf
## [1] "There are no associations to plot!"
## 2026-03-04 18:38:05.822526 INFO::Writing association plots (one for each significant association) to output folder: ./output/maaslin2/pooled-output/rare/20260304_18_37_08-NMR-OTU-sex-234-ref-female-with-relations
## [1] "There are no associations to plot!"
```

Save the fit data object as an rds file.


``` r
# 20260224_21_21_02
# saveRDS(maaslin.fit_data.sex,
#         file = file.path("output/rdafiles",
#                          paste(
#                            paste(format(Sys.time(),format="%Y%m%d"),
#                                  format(Sys.time(),format = "%H_%M_%S"),
#                                  sep = "_"),"maaslin.fit_data.sex",
#                            nmr.sex.data.for_test$output.filename,".rds",sep = "-")))
# "20260215_19_30_37" for both signif and signif.decreased tsv
maaslin.processed_output.sex <-
  process_maaslin2_output(agglom.rank = agglom.rank,
                          maaslin.fit_data = maaslin.fit_data.sex,
                          ps.q.agg = nmr.sex.data.for_test$ps.q.agg,
                          sample.groups = nmr.sex.data.for_test$sample.groups)
```

Downstream analysis and plotting of differentially abundant features. 


``` r
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
```

```
## Joining with `by = join_by(Sample, class)`
```

<img src="007-diffabund-tests_files/figure-html/unnamed-chunk-59-1.png" alt="" width="768" />

``` r
sessionInfo()
```

```
## R version 4.4.3 (2025-02-28 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 22631)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: Asia/Tokyo
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] doRNG_1.8.6.3         rngtools_1.5.2        foreach_1.5.2        
##  [4] Polychrome_1.5.4      phyloseq_1.50.0       ANCOMBC_2.8.1        
##  [7] ALDEx2_1.38.0         latticeExtra_0.6-31   lattice_0.22-6       
## [10] zCompositions_1.5.0-5 survival_3.8-3        truncnorm_1.0-9      
## [13] MASS_7.3-65           Maaslin2_1.20.0       lubridate_1.9.4      
## [16] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
## [19] purrr_1.0.4           readr_2.1.5           tidyr_1.3.1          
## [22] tibble_3.2.1          ggplot2_4.0.0         tidyverse_2.0.0      
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.4.3               cellranger_1.1.0           
##   [3] rpart_4.1.24                lifecycle_1.0.5            
##   [5] Rdpack_2.6.4                doParallel_1.0.17          
##   [7] backports_1.5.0             magrittr_2.0.3             
##   [9] Hmisc_5.2-3                 sass_0.4.10                
##  [11] rmarkdown_2.30              jquerylib_0.1.4            
##  [13] yaml_2.3.12                 otel_0.2.0                 
##  [15] gld_2.6.8                   pbapply_1.7-4              
##  [17] DBI_1.2.3                   minqa_1.2.8                
##  [19] RColorBrewer_1.1-3          ade4_1.7-23                
##  [21] multcomp_1.4-29             abind_1.4-8                
##  [23] zlibbioc_1.52.0             Rtsne_0.17                 
##  [25] expm_1.0-0                  quadprog_1.5-8             
##  [27] GenomicRanges_1.58.0        BiocGenerics_0.52.0        
##  [29] hash_2.2.6.4                nnet_7.3-20                
##  [31] TH.data_1.1-5               sandwich_3.1-1             
##  [33] GenomeInfoDbData_1.2.13     IRanges_2.40.1             
##  [35] S4Vectors_0.44.0            pheatmap_1.0.13            
##  [37] vegan_2.6-4                 microbiome_1.28.0          
##  [39] commonmark_2.0.0            permute_0.9-8              
##  [41] codetools_0.2-20            getopt_1.20.4              
##  [43] DelayedArray_0.32.0         xml2_1.3.8                 
##  [45] ggtext_0.1.2                energy_1.7-12              
##  [47] tidyselect_1.2.1            UCSC.utils_1.2.0           
##  [49] farver_2.1.2                lme4_1.1-37                
##  [51] gmp_0.7-5                   matrixStats_1.5.0          
##  [53] stats4_4.4.3                base64enc_0.1-3            
##  [55] jsonlite_2.0.0              multtest_2.62.0            
##  [57] e1071_1.7-16                Formula_1.2-5              
##  [59] iterators_1.0.14            tools_4.4.3                
##  [61] DescTools_0.99.60           Rcpp_1.0.14                
##  [63] glue_1.8.0                  gridExtra_2.3              
##  [65] SparseArray_1.6.2           mgcv_1.9-1                 
##  [67] xfun_0.56                   MatrixGenerics_1.18.1      
##  [69] GenomeInfoDb_1.42.3         withr_3.0.2                
##  [71] numDeriv_2016.8-1.1         fastmap_1.2.0              
##  [73] rhdf5filters_1.18.1         boot_1.3-31                
##  [75] litedown_0.9                digest_0.6.37              
##  [77] timechange_0.3.0            R6_2.6.1                   
##  [79] colorspace_2.1-2            gtools_3.9.5               
##  [81] markdown_2.0                jpeg_0.1-11                
##  [83] utf8_1.2.4                  generics_0.1.4             
##  [85] data.table_1.17.8           robustbase_0.99-6          
##  [87] class_7.3-23                CVXR_1.0-15                
##  [89] httr_1.4.7                  htmlwidgets_1.6.4          
##  [91] S4Arrays_1.6.0              scatterplot3d_0.3-44       
##  [93] pkgconfig_2.0.3             gtable_0.3.6               
##  [95] Exact_3.3                   Rmpfr_1.1-1                
##  [97] S7_0.2.0                    XVector_0.46.0             
##  [99] pcaPP_2.0-5                 htmltools_0.5.8.1          
## [101] bookdown_0.46               biomformat_1.34.0          
## [103] zigg_0.0.2                  logging_0.10-108           
## [105] scales_1.4.0                Biobase_2.66.0             
## [107] lmom_3.2                    png_0.1-8                  
## [109] optparse_1.7.5              reformulas_0.4.4           
## [111] knitr_1.51                  rstudioapi_0.18.0          
## [113] tzdb_0.5.0                  reshape2_1.4.4             
## [115] checkmate_2.3.3             nlme_3.1-167               
## [117] nloptr_2.2.1                rhdf5_2.50.2               
## [119] proxy_0.4-27                cachem_1.1.0               
## [121] zoo_1.8-14                  rootSolve_1.8.2.4          
## [123] parallel_4.4.3              foreign_0.8-90             
## [125] pillar_1.11.1               grid_4.4.3                 
## [127] vctrs_0.6.5                 cluster_2.1.8.1            
## [129] htmlTable_2.4.3             evaluate_1.0.5             
## [131] mvtnorm_1.3-3               cli_3.6.4                  
## [133] compiler_4.4.3              rlang_1.1.5                
## [135] crayon_1.5.3                labeling_0.4.3             
## [137] interp_1.1-6                plyr_1.8.9                 
## [139] fs_1.6.6                    stringi_1.8.4              
## [141] viridisLite_0.4.3           deldir_2.0-4               
## [143] BiocParallel_1.40.2         lmerTest_3.2-0             
## [145] Biostrings_2.74.1           gsl_2.1-8                  
## [147] Matrix_1.7-4                hms_1.1.4                  
## [149] bit64_4.6.0-1               Rhdf5lib_1.28.0            
## [151] SummarizedExperiment_1.36.0 haven_2.5.5                
## [153] rbibutils_2.3               gridtext_0.1.5             
## [155] Rfast_2.1.5.1               igraph_2.1.4               
## [157] RcppParallel_5.1.11-1       bslib_0.10.0               
## [159] biglm_0.9-3                 DEoptimR_1.1-4             
## [161] directlabels_2025.6.24      bit_4.6.0                  
## [163] readxl_1.4.5                ape_5.8-1
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  9915627 529.6   17299862  924.0  17299862  924.0
## Vcells 17692180 135.0  286184191 2183.5 357730090 2729.3
```

