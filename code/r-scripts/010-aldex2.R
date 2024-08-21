library(tidyverse)
library(ALDEx2)
# Import custom.md, ps.q.df.wide, custom.levels
# input_data_date_time is Rdata workspace from diffabund-input.R
# with a rarefied table in wide format and metadata
input_data_date_time<-"20240809_15_38_22"
# 20240809_14_40_39 for NMR, OTU level
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


# create covariates: age groups
if(host=="NMR" & comparison=="age"){
  covariates<-custom.md$agegroup[match(colnames(ps.q.df.wide),rownames(custom.md))]
}else{
  covariates<-custom.md$class[match(colnames(ps.q.df.wide),rownames(custom.md))]
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
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == aldex.reference), which(colnames(mm) !=aldex.reference))]]
  
}else{
  mm<-mm[, colnames(mm)[c(which(colnames(mm) == paste0("covariates",ref.level)), 
                          which(colnames(mm) != paste0("covariates",ref.level)))]]
}


# aldex glm for a complex case ####
set.seed(1)
ps.q.aldex.clr <- aldex.clr(ps.q.df.wide, mm, mc.samples=1000, denom="all", verbose=F)
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
save.image(paste0("./rdafiles/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "aldex2",paste(custom.levels,collapse = '-'),
                        agglom.rank,comparison,
                        "ref",ref.level,
                        "workspace-test.RData",sep="-")))
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
save.image(paste0("./rdafiles/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "aldex2",paste(custom.levels,collapse = '-'),
                        agglom.rank,comparison,
                        "ref",ref.level,
                        "workspace-effect.RData",sep = "-")))



# Extract significant features ####
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

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "aldex2",paste(custom.levels,collapse = '-'),
  agglom.rank,
  comparison,truncationlvl,"ref",
  ref.level,"workspace.RData",sep="-")))
write.table(aldex.signif.features,
            file=file.path("./output/rtables",authorname,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "aldex2",paste(custom.levels,collapse = '-'),
              agglom.rank,
              comparison,truncationlvl,"ref",
              ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")
