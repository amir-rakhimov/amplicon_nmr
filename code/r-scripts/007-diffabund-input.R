library(tidyverse)
# The final workspace image has custom.md, ps.q.df.wide,
# custom.levels, ref.level, comparison
# ps.q.df.preprocessed.date_time is a rarefied table
# ps.q.df.preprocessed.date_time<-"20240426_22_00_04" # for all hosts, genus level
ps.q.df.preprocessed.date_time<-"20240524_13_58_11" # for NMR, OTU level
# "20240426_21_44_30" workspace Rdata file for all hosts, genus level
# "20240524_13_58_11" workspace Rdata file for NMR # TODO: NOT FOUND

agglom.rank<-"OTU"
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
ref.level<-"NMR"


### If comparing inside host
custom.md<-readRDS("output/rdafiles/custom.md.rds")
# choose what to compare
comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"
# choose the host of interest
host<-"NMR"
# host<-"mice"
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}


# Import data ####
if(agglom.rank=="OTU"){
  ps.q.df.preprocessed<-read.table(
    file.path("./output/rtables",authorname,paste0(
      paste(
        ps.q.df.preprocessed.date_time,
        paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
        paste(names(host.labels),collapse = '-'),sep = "-"),
      ".tsv")),
    header = T)
}else{
  ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)}

if(agglom.rank=="OTU"){
  if(host=="NMR"){
    # select nmr and add age groups
    ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
      filter(class=="NMR",Abundance!=0)%>%
      group_by(Sample)%>%
      mutate(birthday=as.Date(birthday))%>%
      mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))
    # minimum age and maximum age
    min_boundary <- floor(min(ps.q.df.preprocessed$age)/5) * 5
    max_boundary <- ceiling(max(ps.q.df.preprocessed$age)/5) * 5
    
    # add age group to the dataset of abundances
    # each group is 5 years
    # ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    #   mutate(agegroup=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
    #                       include.lowest = TRUE))
    ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
      mutate(agegroup=cut(age, breaks =c(0,10,16),
                          right = FALSE))
    # we create these new levels because maaslin is itsy bitsy
    unique_levels <- ps.q.df.preprocessed %>%
      ungroup()%>%
      distinct(agegroup)%>%
      arrange(agegroup) %>%
      mutate(new_agegroup = paste0("agegroup", agegroup))%>%
      mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
      mutate(new_agegroup = gsub("\\,","_",new_agegroup))
    ps.q.df.preprocessed <- ps.q.df.preprocessed %>%
      left_join(unique_levels, by = "agegroup")
    colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="agegroup")]<-"old_agegroup"
    colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="new_agegroup")]<-"agegroup"
    # add age group to metadata
    custom.md$Sample<-rownames(custom.md)
    custom.md<-custom.md%>% 
      filter(class=="NMR")%>%
      group_by(Sample)%>%
      mutate(birthday=as.Date(birthday))%>%
      mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))%>%
      left_join(unique(ps.q.df.preprocessed[,c("Sample","agegroup")]),by="Sample")%>%
      as.data.frame()
    rownames(custom.md)<-custom.md$Sample
  }else if(host=="mice"){
    # select mice and add age groups: B6, old, or young
    # B6 are separate
    # mice born before 2020 are old
    # after 2023 are young
    custom.levels<-c("B6mouse",
                     "MSMmouse",
                     "FVBNmouse")
    ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
      filter(class%in%custom.levels,Abundance!=0)
    ps.q.df.preprocessed$agegroup<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                          ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
    # also to metadata
    custom.md<-custom.md%>%
      filter(class%in%custom.levels)
    custom.md$agegroup<-ifelse(custom.md$class=="B6mouse","B6",
                               ifelse(grepl("2020",custom.md$birthday),"old","young"))
  }
  # Creating custom levels ####
  if (comparison=="age"){
    # names for levels are age groups
    pretty.facet.labels<-names(table(ps.q.df.preprocessed$agegroup))
    names(pretty.facet.labels)<-names(table(ps.q.df.preprocessed$agegroup))
    custom.levels<-names(pretty.facet.labels)
    
  }else if (comparison=="sex"){
    pretty.facet.labels<-
      c("F" = "Females",
        "M" = "Males")
    custom.levels<-names(pretty.facet.labels)
    
  }else if(comparison=="strain"){
    pretty.facet.labels<-
      c("B6mouse" = "B6 mouse",
        "MSMmouse" = "MSM/Ms mouse",
        "FVBNmouse" = "FVB/N mouse"
      )
    custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
    pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]
  }
  
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
if(agglom.rank=="OTU"){
  objects.to.keep<-c(objects.to.keep,"host","comparison")
  objects.to.keep<-which(ls()%in%objects.to.keep)
  rm(list = ls()[-objects.to.keep])
  save.image(file.path("./output/rdafiles",paste(
    paste(format(Sys.time(),format="%Y%m%d"),
          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
    "diffabund-input-data",host,rare.status,filter.status,agglom.rank,
    comparison,truncationlvl,
    paste(custom.levels,collapse = '-'),
    "workspace.RData",sep="-")))
}else{
  objects.to.keep<-which(ls()%in%objects.to.keep)
  rm(list = ls()[-objects.to.keep])
  save.image(file.path("./output/rdafiles",paste(
    paste(format(Sys.time(),format="%Y%m%d"),
          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
    "diffabund-input-data",rare.status,filter.status,agglom.rank,
    truncationlvl,
    paste(custom.levels,collapse = '-'),
    "workspace.RData",sep="-")))
}


