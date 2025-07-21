# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("Maaslin2","ALDEx2","ANCOMBC","phyloseq"))
# BiocManager::install("ALDEx2")
# BiocManager::install("ANCOMBC")
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","Polychrome"))
# Differential abundance tests ####
library(tidyverse)
library(Maaslin2)
library(ALDEx2)
library(ANCOMBC)
library(phyloseq)
library(Polychrome)
ps.q.df.preprocessed.date_time<-"20240809_13_18_49" 
# "20240426_22_00_04" rarefied table for all hosts, genus level
# "20240524_13_58_11" rarefied table file for NMR, OTU level
# 20240809_13_18_49 rarefied table file for NMR, genus level
inside.host<-TRUE
if(inside.host=="TRUE"){
  # choose what to compare
  comparison<-"age"
  # comparison<-"sex"
  # comparison<-"strain"
  ref.level<-"agegroup0_10" # choose the reference level
  custom.levels<-"NMR"
  # custom.levels<-c("B6mouse",
  #                  "MSMmouse",    
  #                  "FVBNmouse")
}else{
  comparison<-"host"
  ref.level<-"NMR"
  custom.levels<-c("NMR",
                   "DMR",
                   "B6mouse",
                   "MSMmouse",
                   "FVBNmouse",
                   "spalax",
                   "pvo",
                   "hare",
                   "rabbit")
}
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
authorname<-"pooled"
ps.q.agg.date_time<-"20240613_21_47_12"
# 20240613_21_47_12 phyloseq OTU table
# 20240613_21_42_48 phyloseq Genus table

output.filename<-paste(paste(custom.levels,collapse = '-'),
                       agglom.rank,comparison,
                       truncationlvl,"ref",ref.level,
                       sep="-")


# Import data ####
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

if(setequal(custom.levels,"NMR")){
  custom.md<-readRDS("./output/rdafiles/custom.md.ages.rds")
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class=="NMR",Abundance!=0)
}else if(setequal(custom.levels,c("B6mouse","MSMmouse","FVBNmouse"))){
  # select mice and add age groups: B6, old, or young
  # B6 are separate
  # mice born before 2020 are old
  # after 2023 are young
  custom.md<-readRDS("./output/rdafiles/custom.md.rds")
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,Abundance!=0)
  custom.md<-custom.md%>%
    filter(class%in%custom.levels)
  custom.md$agegroup<-ifelse(custom.md$class=="B6mouse","B6",
                                        ifelse(grepl("2020",custom.md$birthday),"old","young"))
}else{
  custom.md<-readRDS("./output/rdafiles/custom.md.rds")
}

# Preparing the dataset ####
# filter the dataset
ps.q.df.wide<-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c("Sample",agglom.rank,"Abundance","class")))%>%
  filter(Abundance!=0)%>%
  pivot_wider(names_from = agglom.rank, # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()%>%
  column_to_rownames("Sample")%>%
  select(-class)
# colnames are OTUs and rownames are sample IDs


## MaAsLin 2 ####
### Create reference levels  ####
if (comparison=="age"){
  maaslin.reference<-paste("agegroup",ref.level,sep = ",")
  maaslin.comparison<-"agegroup"
}else if(comparison=="sex"){
  maaslin.reference<-paste("sex","F",sep = ",")
  maaslin.comparison<-"sex"
}else if(comparison=="strain"){
  maaslin.reference<-paste("class",ref.level,sep = ",")
  maaslin.comparison<-"class"
}

if(setequal(custom.levels,"NMR")){
  relations<-read.table("./data/metadata/pooled-metadata/nmr-relations.tsv",
                        header = T,
                        sep = "\t")
  custom.md<-custom.md%>%
    left_join(relations,by="Sample")
  rownames(custom.md)<-custom.md$Sample
}

### Run Maaslin2 ####
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
             output = file.path("./output/maaslin2",paste0(authorname,"-output"),
                                rare.status,paste(
                                  paste(format(Sys.time(),format="%Y%m%d"),
                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                  output.filename,sep = "-")), 
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
             output = file.path("./output/maaslin2",paste0(authorname,"-output"),
                                rare.status,paste(
                                  paste(format(Sys.time(),format="%Y%m%d"),
                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                  output.filename,sep = "-")), 
             fixed_effects = c("class"),
             reference = paste0("class,",ref.level),
             max_significance = 0.05)
  
}
### Save the workspace ####
save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2",output.filename,"workspace.RData",sep="-")))


## Downstream processing of Maaslin2 output ####
source("./code/r-scripts/make_features_maaslin.R")
# Import ps.q.agg
ps.q.agg<-read.table(file.path("./output/rtables",authorname,paste(
  paste(ps.q.agg.date_time,
        "phyloseq-qiime",authorname,agglom.rank,read.end.type,truncationlvl,
        "table.tsv",sep="-"))
),header = T)

### Extract features with qvalue<0.05 ####
if(min(maaslin.fit_data$results$qval)<0.05){
  maaslin.signif.features<-maaslin.fit_data$results%>%
    filter(qval<0.05) # should be qval
}else{
  maaslin.signif.features<-maaslin.fit_data$results%>%
    arrange(qval)%>%
    head(n = 100) # if no significant results found
}

### Make features pretty ####
# We have to make this exchange because maaslin output treats space and hyphen
# as the same thing
if(agglom.rank=="OTU"){
  maaslin.signif.features$feature<-gsub("^X","",maaslin.signif.features$feature)
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
  maaslin.signif.features<-subset(maaslin.signif.features, select=-get(agglom.rank))
}

### Extract features that are downregulated in all hosts simultaneously ####
if(setequal(custom.levels,"NMR")){
  if(comparison=="age"){
    sample.groups<-c("agegroup0_10",
                     "agegroup10_16")
  }else if(comparison=="sex"){
    sample.groups<-c("F",
                     "M")
  }
}
# Downreglated features have coef<0.
# We also add association strength column
maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by(feature)%>%
  filter(n()==length(sample.groups)-1)%>% 
  arrange(pval,feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

# maaslin.signif.increased<-maaslin.signif.features%>%
#   as_tibble()%>%
#   filter(coef>0)%>%
#   group_by(feature)%>%
#   filter(n()==length(sample.groups)-1)%>% 
#   arrange(feature)%>%
#   mutate(n=n())%>%
#   mutate(assoc.str=-log(qval)*sign(coef))#%>%
# # select(feature,assoc.str,name)
# Check if all features (OTU/taxa/etc) are also present in the ps.q.agg dataframe
if(agglom.rank=="OTU"){
  table(maaslin.signif.decreased$feature%in%ps.q.agg$OTU)
}else{
  table(maaslin.signif.decreased$feature%in%pull(ps.q.agg[,agglom.rank]))
}
write.table(maaslin.signif.features,
            file=file.path("./output/rtables",authorname,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "maaslin2",output.filename,
              "signif.tsv",sep="-")),
            row.names = F,sep = "\t")
write.table(maaslin.signif.decreased,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "maaslin.signif.decreased",
                                 output.filename,
                                 "signif.tsv",sep="-")),
            row.names = F,sep = "\t")




## ALDEx2 ####
### Create covariates: take the samples in the ps.q.df.wide, and find  ####
# the positions of these samples in the custom.md.
# Ref: match returns a vector of the positions of (first) matches of its first argument in its second.
# We're searching the elements of the first vector in the second vector
if(setequal(custom.levels,"NMR")& comparison=="age"){
  covariates<-custom.md$agegroup[match(rownames(ps.q.df.wide),rownames(custom.md))]
}else if (inside.host==FALSE){
  covariates<-custom.md$class[match(rownames(ps.q.df.wide),rownames(custom.md))]
}
# model matrix is covariates and the reference group is the first covariate
mm <- model.matrix(~ covariates-1)
# reorder model.matrix to put ref.level as first column

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
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == aldex.reference), which(colnames(mm) !=aldex.reference))]]
  
}else{
  # Same reordering except this is the case when we have animal hosts. We 
  # directly put ref.level into the string
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == paste0("covariates",ref.level)), 
                          which(colnames(mm) != paste0("covariates",ref.level)))]]
}


### Aldex glm for a complex case ####
set.seed(1)
ps.q.aldex.clr <- aldex.clr(ps.q.df.wide, mm, mc.samples=1000, denom="all", verbose=F)
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
save.image(paste0("./rdafiles/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "aldex2",output.filename,
                        "workspace-test.RData",sep="-")))
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
save.image(paste0("./rdafiles/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "aldex2",output.filename,
                        "workspace-effect.RData",sep = "-")))


### Extract significant features ####
aldex.signif.features<-list()
for (i in 1:length(ps.q.glm.effect)){
  # take all features that have good CI (not overlapping zero) and high effect size
  # identify features with significant effect size and good CI
  sig<-which((ps.q.glm.effect[[i]]$effect.low>0 & ps.q.glm.effect[[i]]$effect.high>0)|
               (ps.q.glm.effect[[i]]$effect.low<0 & ps.q.glm.effect[[i]]$effect.high<0))
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
# check for duplicates
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

write.table(aldex.neg.effect,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "aldex.neg.effect",
                                 output.filename,
                                 "signif.tsv",sep="-")),
            row.names = F,sep = "\t")

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "aldex2",output.filename,"workspace.RData",sep="-")))
write.table(aldex.signif.features,
            file=file.path("./output/rtables",authorname,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "aldex2",output.filename,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")



## ANCOMBC ####
### Load the ps.q.agg object from 001-phyloseq-qiime2.R ####
ps.q.agg<-readRDS(file=file.path("./output/rdafiles",paste(
  ps.q.agg.date_time,
  "phyloseq-qiime",authorname,agglom.rank,read.end.type,truncationlvl,
  "table.rds",sep="-")))


### Extract a taxonomic table ####
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
### Create a phyloseq object for ANCOMBC ####
ps.q.OTU<-t(ps.q.df.wide)
ps.q.OTU<-otu_table(ps.q.OTU,taxa_are_rows = T)
ps.q.TAX<-tax_table(taxmat)
ps.q.phyloseq.new<-phyloseq(otu_table(ps.q.OTU),
                            tax_table(ps.q.TAX),
                            sample_data(custom.md))

### Relevel the comparison vector. The first level will be the reference ####
# for custom leveling
ancombc.levels<-c(ref.level,sample.groups[sample.groups!=ref.level])
if(comparison=="host"){
  sample_data(ps.q.phyloseq.new)$class<-factor(sample_data(ps.q.phyloseq.new)$class,
                                               levels = ancombc.levels)
  ancombc.comparison<-"class"                                           
}else if (comparison=="age"){
  sample_data(ps.q.phyloseq.new)$agegroup<-factor(sample_data(ps.q.phyloseq.new)$agegroup,
                                                  levels = ancombc.levels)
  ancombc.comparison<-"agegroup"
}else if(comparison=="sex"){
  sample_data(ps.q.phyloseq.new)$sex<-factor(sample_data(ps.q.phyloseq.new)$sex,
                                             levels = ancombc.levels)
  ancombc.comparison<-"sex"
}else if(comparison=="strain"){
  sample_data(ps.q.phyloseq.new)$class<-factor(sample_data(ps.q.phyloseq.new)$class,
                                               levels = ancombc.levels)
  ancombc.comparison<-"class"
}

### Perform differential abundance test ####
ancombc.out<-ancombc(
  phyloseq = ps.q.phyloseq.new,
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
ancombc.res<-ancombc.out$res

### Find differentially abundant taxa by multiplying fold change with TRUE/FALSE ####
df_lfc = data.frame(ancombc.res$lfc[, -1] * ancombc.res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancombc.res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
ancombc.signif.features<-subset(df_lfc,rowSums(df_lfc[,-c(1,2)]!=0)==ncol(df_lfc[,-c(1,2)]))

ancombc.signif.decreased<-subset(ancombc.signif.features,
                                 rowSums(ancombc.signif.features[,-c(1,2)]<0)==ncol(ancombc.signif.features[,-c(1,2)]))


save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "ancombc",output.filename,"workspace.RData",sep="-")))
write.table(ancombc.signif.features,
            file=file.path("./rtables",authorname,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "ancombc",output.filename,
              "signif.tsv",sep="-")), 
            row.names = F,sep = "\t")

write.table(ancombc.signif.decreased,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "ancombc.signif.decreased",
                                 output.filename,
                                 "signif.tsv",sep="-")),
            row.names = F,sep = "\t")


# Downstream processing of the test output ####
authorname<-"pooled"
inside.host<-TRUE
if(inside.host=="TRUE"){
  # choose what to compare
  comparison<-"age"
  # comparison<-"sex"
  # comparison<-"strain"
  custom.levels<-"NMR"
  if(comparison=="age"){
    sample.groups<-c("agegroup0_10",
                     "agegroup10_16")
  }else if(comparison=="sex"){
    sample.groups<-c("F",
                     "M")
  }
  # NOT workspace dates but signif table dates
  maaslin2.signif_all.date_time<-c("agegroup0_10"="20240621_17_53_06")
  # maaslin2.signif_all.date_time<-c("F"="20240621_17_54_53")
  ref.level<-"agegroup0_10" # choose the reference level
  # ref.level<-"F"
  if(comparison=="age"){
    sample.groups<-c("agegroup0_10",
                     "agegroup10_16")
  }else if(comparison=="sex"){
    sample.groups<-c("F",
                     "M")
  }
  
}else{
  comparison<-"host"
  ref.level<-"NMR"
  custom.levels<-c("NMR",
                   "B6mouse",
                   "MSMmouse",
                   "FVBNmouse",
                   "DMR",
                   "hare",
                   "rabbit",
                   "spalax",
                   "pvo")
  # Tables of differentially abundant taxa
  maaslin2.signif_all.date_time<-"20240427_16_22_09"
  maaslin2.signif_decreased.date_time<-"20240429_17_15_32"
  aldex.signif_all.date_time<-"20240427_17_39_41"
  aldex.neg.effect.date_time<-"20240429_17_15_56"
  ancombc.signif_all.date_time<-"20240427_17_12_27"
  ancombc.signif.decreased.date_time<-"20240429_17_16_00"
}

phyloseq.workspace.date_time<-"20240524_13_54_21"
# 20240524_13_54_21 for all hosts, OTU level
# 20240426_21_44_30 for all hosts, genus level
truncationlvl<-"234"
agglom.rank<-"OTU"
authorname<-"pooled"
rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"
rtables.directory<-file.path("./output/rtables",authorname)
image.formats<-c("png","tiff")
# This is in all files
output.filename<-paste(paste(custom.levels,collapse = '-'),
                       agglom.rank,comparison,
                       truncationlvl,"ref",ref.level,
                       sep="-")

### Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R) ####
load(file.path("./output/rdafiles",paste(
  phyloseq.workspace.date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)

### Load tables ####
if(inside.host==TRUE){
  maaslin.signif.features<-
    read.table(file.path(rtables.directory,
                         paste(maaslin2.signif_all.date_time,"maaslin2",
                               output.filename,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  maaslin.signif.decreased<-
    read.table(file.path(rtables.directory,
                         paste(maaslin2.signif_decreased.date_time,"maaslin.signif.decreased",
                               output.filename, 
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  
  # aldex.signif.features<-
  #   read.table(file.path(rtables.directory,
  #                        paste(aldex.signif_all.date_time,"aldex2",
  #                              signif.tsv.filename,sep="-")),
  #              header = T,sep = "\t")
  # ancombc.signif.features<-
  #   read.table(file.path(rtables.directory,
  #                        paste(ancombc.signif_all.date_time,"ancombc",
  #                              signif.tsv.filename,sep="-")),
  #              header = T,sep = "\t")
  # ancombc.signif.features<-lvl.ancombc.signif.features%>%
  #   dplyr::select(-X.Intercept.)%>%
  #   pivot_longer(!taxon_id,names_to = "class",values_to = "coef")
}else{
  maaslin.signif.features<-
    read.table(file.path(rtables.directory,
                         paste(maaslin2.signif_all.date_time,"maaslin2",
                               output.filename,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  maaslin.signif.decreased<-
    read.table(file.path(rtables.directory,
                         paste(maaslin2.signif_decreased.date_time,"maaslin.signif.decreased",
                               output.filename, 
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  
  aldex.signif.features<-
    read.table(file.path(rtables.directory,
                         paste(aldex.signif_all.date_time,"aldex2",
                               output.filename, 
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  aldex.neg.effect<-
    read.table(file.path(rtables.directory,
                         paste(aldex.neg.effect.date_time,"aldex.neg.effect",
                               output.filename,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  ancombc.signif.features<-
    read.table(file.path(rtables.directory,
                         paste(ancombc.signif_all.date_time,"ancombc",
                               output.filename,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  ancombc.signif.decreased<-
    read.table(file.path(rtables.directory,
                         paste(ancombc.signif.decreased.date_time,"ancombc.signif.decreased",
                               output.filename,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
}


if(inside.host!=TRUE){
  ### Find common significant features between three tools ####
  print(Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$Taxon,
                        ancombc.signif.features$taxon_id)))
  
  ### Find common significantly decreased features between three tools ####
  print(Reduce(intersect,list(maaslin.signif.decreased$feature,aldex.neg.effect$Taxon,
                        ancombc.signif.decreased$taxon_id)))
  
  # Common significant and decreased features between tools
  common.signif<-Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$Taxon,
                                       ancombc.signif.features$taxon_id))
  ### Only Maaslin2 and ANCOM-BC: they're decreased in other hosts, not in ref ####
  common.decreased<-Reduce(intersect,list(maaslin.signif.decreased$feature,
                                          ancombc.signif.decreased$taxon_id))
  
  ps.q.agg%>%
    filter(class==ref.level,Genus%in%common.decreased)%>%
    distinct(Genus,.keep_all = T)%>%
    ungroup()%>%
    arrange(-MeanRelativeAbundance)%>%
    dplyr::select(Genus,MeanRelativeAbundance)%>%
    head(n=10)%>%
    print()
  
  write.table(common.decreased,
              file=file.path(rtables.directory,
                             paste(paste(format(Sys.time(),format="%Y%m%d"),
                                         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                   "significant-features",
                                   output.filename,
                                   "signif.tsv",sep="-")),
              row.names = F,sep = "\t",col.names = F)
}


# Plot differentially abundant species ####
if(setequal(custom.levels,"NMR")){
  # select nmr and add age groups
  custom.md<-custom.md%>%
    filter(class=="NMR")%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))%>%
    mutate(agegroup=cut(age, breaks =c(0,10,16),
                        right = FALSE))
}else if(setequal(custom.levels,c("B6mouse","MSMmouse","FVBNmouse"))){
  # select mice and add age groups
  custom.md<-custom.md%>%
    filter(Abundance!=0)
  custom.md$agegroup<-ifelse(custom.md$class=="B6mouse","B6",
                                        ifelse(grepl("2020",custom.md$birthday),"old","young"))
}

if(comparison=="host"){
  pretty.level.names<-
    c("NMR" = "*Heterocephalus glaber*", # better labels for facets
      "B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse",
      "DMR" = "*Fukomys Damarensis*",
      "hare" = "*Lepus europaeus*",
      "rabbit" = "*Oryctolagus cuniculus*",
      "spalax" = "*Nannospalax leucodon*",
      "pvo" = "*Pteromys volans orii*")
  ggplot.levels<-names(pretty.level.names)
  gg.labs.name<-"Age group"
  gg.title.groups<-"age groups"
}else if (comparison=="age"){
  pretty.level.names<-names(table(custom.md$agegroup))
  custom.md<-custom.md%>%
    ungroup()%>%
    mutate(old_agegroup=agegroup)%>%
    mutate(agegroup = paste0("agegroup", agegroup))%>%
    mutate(agegroup = gsub("\\(|\\)|\\[|\\]","",agegroup))%>%
    mutate(agegroup = gsub("\\,","_",agegroup))
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
    c("F" = "Females",
      "M" = "Males")
  ggplot.levels<-names(pretty.level.names)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%ggplot.levels)]
  gg.labs.name<-"Host sex"
  gg.title.groups<-"groups"
}else if(comparison=="strain"){
  pretty.level.names<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  ggplot.levels<-intersect(names(pretty.level.names),custom.md$class)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%ggplot.levels)]
  gg.labs.name<-"Strain"
  gg.title.groups<-"strains"
}

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#BF3EFF","#5CACEE","#00CD66",
                                          "#FF8C00","#00EE00","#EEC900", "#00FFFF",
                                          "#FF6EB4",
                                          "#FFA07A"))
names(custom.fill)<-ggplot.levels
swatch(custom.fill)

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
  
  feature.plot<-ps.q.agg%>%
    filter(get(agglom.rank)%in%maaslin.signif.features$feature,class=="NMR")%>%
    left_join(custom.md)%>%
    group_by_at(c(agglom.rank,"agegroup"))%>%
    ggplot(aes(x=factor(agegroup,levels=names(pretty.level.names)),
               y=RelativeAbundance,
               fill=factor(agegroup)))+
    geom_boxplot(show.legend = FALSE)+
    facet_wrap(~OTU,scales = "free_y",
               ncol = 2,
               labeller = as_labeller(pretty.asv.names))+
    theme_bw()+
    labs(x="",
         y="Relative abundance (%)")+
    scale_x_discrete(labels=pretty.level.names,
                     limits=ggplot.levels)+ # rename boxplot labels (x axis)
    scale_fill_manual(values = custom.fill)+
    theme(axis.title = element_text(size = 20),
          axis.text.y = ggtext::element_markdown(size=18),
          axis.text.x = element_text(size=20),
          strip.text.x= ggtext::element_markdown(size=20),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "none")+
    ggtitle(paste0("Relative abundance of differentially abundant ASVs \nin different naked mole-rat groups"))
  for(image.format in image.formats){
    ggsave(paste0("./images/taxaboxplots/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        ref.level,"specific-bacteria","NMR",comparison,
                        sep = "-"),".",image.format),
           plot=feature.plot,
           width = 4000,height = 2000,
           units = "px",dpi=300,device = image.format)
}
}else{
  feature.plot<-ps.q.agg%>%
    filter(get(agglom.rank)%in%common.decreased)%>%
    group_by_at(c("class",agglom.rank))%>%
    ggplot(aes(x=factor(class,level=rev(ggplot.levels)),
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
                     limits=rev(ggplot.levels))+ # rename boxplot labels (x axis)
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
    ggtitle(paste0("Relative abundance of naked mole-rat-specific taxa"))
  for(image.format in image.formats){
    ggsave(paste0("./images/taxaboxplots/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "NMR-specific-bacteria",
                        sep = "-"),".",image.format),
           plot=feature.plot,
           width = 4000,height = 12000,
           units = "px",dpi=300,device = image.format)
  }

}
