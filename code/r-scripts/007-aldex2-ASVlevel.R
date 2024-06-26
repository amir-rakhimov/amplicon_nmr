library(tidyverse)
library(phyloseq)
library(ALDEx2)
library(vegan)
authorname<-"pooled"
# choose what to compare
comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"
# choose the host of interest
# host<-"NMR"
host<-"mice"
ref.level<-"FVBNmouse" # choose the reference level
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}
truncationlvl<-"234"
agglom.rank<-"OTU"
read.end.type<-"single"
load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))
rare.status<-"rare"
filter.status<-"nonfiltered"
# Import data ####
ps.q.df.preprocessed<-read.table(paste0("./rtables/",authorname,"/ps.q.df.",
                                        rare.status,".",filter.status,"-",agglom.rank,"-",
                                        paste(names(host.labels),collapse = '-'),".tsv"),
                                 header = T,sep = "\t")
if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class=="NMR",Abundance!=0)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=year(as.period(interval(birthday,now()))))
  # minimum age and maximum age
  min_boundary <- floor(min(ps.q.df.preprocessed$age)/5) * 5
  max_boundary <- ceiling(max(ps.q.df.preprocessed$age)/5) * 5
  # add age group to the dataset of abundances
  # each group is 5 years
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(age_group=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                          include.lowest = TRUE))
  # add age group to metadata
  custom.md$Sample<-rownames(custom.md)
  custom.md<-custom.md%>% 
    filter(class=="NMR")%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=year(as.period(interval(birthday,now()))))
  min_boundary <- floor(min(custom.md$age)/5) * 5
  max_boundary <- ceiling(max(custom.md$age)/5) * 5
  
  custom.md<-custom.md%>%
    mutate(agegroup=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                        include.lowest = TRUE))%>%
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
  ps.q.df.preprocessed$age_group<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
}
# creating custom levels ####
if (comparison=="age"){
  # names for levels are age groups
  pretty.facet.labels<-names(table(ps.q.df.preprocessed$age_group))
  names(pretty.facet.labels)<-names(table(ps.q.df.preprocessed$age_group))
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
# Preparing the dataset ####
# filter the dataset
ps.q.df <-ps.q.df.preprocessed%>%
  dplyr::select(Sample,OTU,Abundance,class,age_group,sex)%>%
  filter(Abundance!=0)
# convert the data frames into wide format
ps.q.df.aldex.input.wide<-ps.q.df%>%
  dplyr::select(Sample,Abundance,OTU)%>%
  pivot_wider(names_from = "OTU", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.aldex.input.wide)<-ps.q.df.aldex.input.wide$Sample
ps.q.df.aldex.input.wide<-ps.q.df.aldex.input.wide[,-1] 
ps.q.df.aldex.input.wide<-t(ps.q.df.aldex.input.wide)
# create covariates: age groups
if(host=="NMR" & comparison=="age"){
  covariates<-custom.md$agegroup[match(colnames(ps.q.df.aldex.input.wide),rownames(custom.md))]
}else{
  covariates<-custom.md$class[match(colnames(ps.q.df.aldex.input.wide),rownames(custom.md))]
}
# model matrix is covariates and the reference group is the first covariate
mm <- model.matrix(~ covariates-1)
# reorder model.matrix to put ref.level as first column
if (comparison=="age"){
  aldex.reference<-paste0("covariates",ref.level)
}else if(comparison=="sex"){
  aldex.reference<-"covaritatesF"
}else if(comparison=="strain"){
  aldex.reference<-paste0("covariates",ref.level)
}
mm<-mm[, colnames(mm)[c(which(colnames(mm) == aldex.reference), which(colnames(mm) !=aldex.reference))]]
# aldex glm for a complex case ####
set.seed(1)
ps.q.aldex.clr <- aldex.clr(ps.q.df.aldex.input.wide, mm, mc.samples=1000, denom="all", verbose=F)
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
save.image(paste0("./rdafiles/",
                  paste("aldex2",rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,ref.level,
                        "workspace-test.RData",sep="-")))
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
save.image(paste0("./rdafiles/",
                  paste("aldex2",rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,ref.level,
                        "workspace-effect.RData",sep = "-")))
# Extract significant features ####
aldex.signif.features<-list()
for (i in 1:length(ps.q.glm.effect)){
  # take all features that have good CI and high effect size
  # identify features with significant effect size and good CI
  sig<-which((ps.q.glm.effect[[i]]$effect.low>0 & ps.q.glm.effect[[i]]$effect.high>0)|
               (ps.q.glm.effect[[i]]$effect.low<0 & ps.q.glm.effect[[i]]$effect.high<0))
  sig<-intersect(which(abs(ps.q.glm.effect[[i]]$effect)>1), sig)
  
  if(length(sig)!=0){
    signif.df<-ps.q.glm.effect[[i]][sig,]
    signif.df$OTU<-rownames(signif.df)
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
aldex.signif.features%>%
  dplyr::group_by(OTU)%>%
  summarise(n=n())%>%
  arrange(-n)
save.image(paste0("./rdafiles/",
                  paste("aldex2",rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,ref.level,
                        "workspace.RData",sep = "-")))
write.table(aldex.signif.features,
            file=paste0("./rtables/",authorname,"/",
                        paste("aldex2",rare.status,filter.status,host,agglom.rank,
                              comparison,truncationlvl,ref.level,
                              "signif.tsv",sep="-")), 
            row.names = F,sep = "\t")
q()
