# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("Maaslin2"))
# install.packages(c("tidyverse"))
library(tidyverse)
library(Maaslin2)
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

load(file.path("./output/rdafiles",paste(
  input_data_date_time,
  "diffabund-input-data",rare.status,filter.status,agglom.rank,
  truncationlvl,
  paste(custom.levels,collapse = '-'),"workspace.RData",sep="-")))

# Running MaAsLin 2 ####
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
                                  paste(custom.levels,collapse = '-'),
                                  agglom.rank,comparison,
                                  "ref",ref.level,sep = "-")), 
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
                                  paste(custom.levels,collapse = '-'),
                                  agglom.rank,comparison,
                                  "ref",ref.level,sep = "-")), 
             fixed_effects = c("class"),
             reference = paste0("class,",ref.level),
             max_significance = 0.05)
  
}
# Save the workspace ####
save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2",paste(custom.levels,collapse = '-'),
  agglom.rank,comparison,
  truncationlvl,"ref",
  ref.level,"workspace.RData",sep="-")))


# Downstream processing ####
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

# Extract features that are downregulated in all hosts simultaneously ####
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
              "maaslin2",paste(custom.levels,collapse = '-'),
              agglom.rank,comparison,
              truncationlvl,"ref",
              ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")
