# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("Maaslin2"))
# install.packages(c("tidyverse"))
library(tidyverse)
library(Maaslin2)
# Import custom.md, ps.q.df.wide, custom.levels
# input_data_date_time is Rdata workspace from diffabund-input.R
# with a rarefied table in wide format and metadata
input_data_date_time<-"20240809_15_38_22"
# 20240809_14_52_10 for all hosts, genus level
# 20240809_15_38_22 for NMR, genus level
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
                   "pvo"
                   # ,
                   # "NMRwt"
                   )
  
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

# 3.2 Running MaAsLin 2 ####
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

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2",paste(custom.levels,collapse = '-'),
  agglom.rank,comparison,
  truncationlvl,"ref",
  ref.level,"workspace.RData",sep="-")))


q()