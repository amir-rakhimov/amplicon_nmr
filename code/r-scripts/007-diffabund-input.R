# install.packages(c("tidyverse"))
library(tidyverse)
# We use this script to prepare data for differential microbial abundance test.
# The final workspace image has custom.md, ps.q.df.wide,
# custom.levels, rare.status,filter.status,agglom.rank,
# truncationlvl,read.end.type.
# ps.q.df.preprocessed.date_time is a rarefied table from 
# 004-phyloseq-rarefaction-filtering.R
ps.q.df.preprocessed.date_time<-"20240809_13_18_49" 
# "20240426_22_00_04" rarefied table for all hosts, genus level
# "20240524_13_58_11" rarefied table file for NMR, OTU level
# 20240809_13_18_49 rarefied table file for NMR, genus level

agglom.rank<-"Genus"
authorname<-"pooled"
truncationlvl<-"234"
read.end.type<-"single"
# load(file.path("./output/rdafiles",paste(
#   date_time,
#   authorname,read.end.type,"qiime2",
#   truncationlvl,agglom.rank,
#   "phyloseq-workspace.RData",sep = "-")))

rare.status<-"rare"
filter.status<-"nonfiltered"


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
custom.levels<-"NMR"
custom.levels<-c("B6mouse",
                 "MSMmouse",    
                 "FVBNmouse")
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
    filter(class=="NMR",Abundance!=0)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))%>%
    mutate(class=as.factor(class),
           sex=as.factor(sex))
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    left_join(custom.md)
}else if(setequal(custom.levels,c("B6mouse","MSMmouse","FVBNmouse"))){
  # select mice and add age groups: B6, old, or young
  # B6 are separate
  # mice born before 2020 are old
  # after 2023 are young
  custom.md<-readRDS("./output/rdafiles/custom.md.rds")
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,Abundance!=0)
  ps.q.df.preprocessed$agegroup<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                        ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
  # also to metadata
  custom.md<-custom.md%>%
    filter(class%in%custom.levels)
  custom.md$agegroup<-ifelse(custom.md$class=="B6mouse","B6",
                             ifelse(grepl("2020",custom.md$birthday),"old","young"))
}else{
  custom.md<-readRDS("./output/rdafiles/custom.md.rds")
}


# Preparing the dataset ####
# filter the dataset
ps.q.df <-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c("Sample",agglom.rank,"Abundance","class")))%>%
  filter(Abundance!=0)
ps.q.df.wide<-ps.q.df%>%
  pivot_wider(names_from = agglom.rank, # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()%>%
  column_to_rownames("Sample")%>%
  select(-class)
# colnames are OTUs and rownames are sample IDs

objects.to.keep<-c("rare.status","filter.status","agglom.rank",
                   "truncationlvl","read.end.type",
                   "custom.levels","ps.q.df.wide","custom.md")

objects.to.keep<-which(ls()%in%objects.to.keep)
rm(list = ls()[-objects.to.keep])
save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "diffabund-input-data",rare.status,filter.status,agglom.rank,
  truncationlvl,
  paste(custom.levels,collapse = '-'),
  "workspace.RData",sep="-")))

