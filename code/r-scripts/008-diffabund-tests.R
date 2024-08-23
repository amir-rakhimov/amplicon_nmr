# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("Maaslin2","ALDEx2","ANCOMBC","phyloseq"))
# BiocManager::install("ALDEx2")
# BiocManager::install("ANCOMBC")
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse"))
# Differential abundance tests ####
library(tidyverse)
library(Maaslin2)
library(ALDEx2)
library(ANCOMBC)
library(phyloseq)
# Import custom.md, ps.q.df.wide, custom.levels ####
# input_data_date_time is Rdata workspace from diffabund-input.R
# with a rarefied table in wide format and metadata
input_data_date_time<-"20240809_14_40_39"
# 20240809_14_52_10 for all hosts, genus level
# 20240809_14_40_39 for NMR, OTU level
# 20240809_15_38_22 for NMR, genus level
ps.q.agg.date_time<-"20240613_21_47_12"
# 20240613_21_47_12 phyloseq OTU table
# 20240613_21_42_48 phyloseq Genus table
authorname<-"pooled"
inside.host<-TRUE
if(inside.host=="TRUE"){
  # choose what to compare
  comparison<-"age"
  # comparison<-"sex"
  # comparison<-"strain"
  ref.level<-"agegroup0_10" # choose the reference level
  custom.levels<-"NMR"
  # custom.levels<-c("agegroup0_10",
  #                                 "agegroup10_16")
  # custom.levels<-c("F",
  #                  "M")
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
  
}

truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
output.filename<-paste(paste(custom.levels,collapse = '-'),
                           agglom.rank,comparison,
                           truncationlvl,"ref",ref.level,
                           sep="-")



load(file.path("./output/rdafiles",paste(
  input_data_date_time,
  "diffabund-input-data",rare.status,filter.status,agglom.rank,
  truncationlvl,
  paste(custom.levels,collapse = '-'),"workspace.RData",sep="-")))

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
