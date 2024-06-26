library(tidyverse)
library(phyloseq)
library(ALDEx2)
library(vegan)
ref.level<-"NMR"
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
# Import data ####
date_time<-"20240426_21_44_30"
load(file.path("./output/rdafiles",paste(
  date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

custom.levels<-c("NMR",
                 "B6mouse",
                 "MSMmouse",
                 "FVBNmouse",
                 "DMR",
                 "hare",
                 "rabbit",
                 "spalax",
                 "pvo"#,
                 # "NMRwt"
                 )
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)

ps.q.df.aldex.input<- read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      "20240426_22_00_04",
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)
# Prepare the dataset ####
# convert the data frame into wide format
ps.q.df.aldex.input.wide<-ps.q.df.aldex.input%>%
  filter(class%in%custom.levels,Abundance!=0)%>%
  dplyr::select(Sample,Abundance,all_of(agglom.rank))%>%
  pivot_wider(names_from = all_of(agglom.rank),
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.aldex.input.wide)<-ps.q.df.aldex.input.wide$Sample
ps.q.df.aldex.input.wide<-ps.q.df.aldex.input.wide[,-1]  
ps.q.df.aldex.input.wide<-t(ps.q.df.aldex.input.wide)

covariates<-custom.md$class[match(colnames(ps.q.df.aldex.input.wide),rownames(custom.md))]
mm <- model.matrix(~ covariates-1)
# reorder model.matrix to put ref.level as first column
mm<-mm[, colnames(mm)[c(which(colnames(mm) == paste0("covariates",ref.level)), 
                        which(colnames(mm) != paste0("covariates",ref.level)))]]

# aldex glm for a complex case ####
set.seed(1)
ps.q.aldex.clr <- aldex.clr(ps.q.df.aldex.input.wide, mm, mc.samples=1000, denom="all", verbose=F)
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)
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
    signif.df$Taxon<-rownames(signif.df)
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
  dplyr::group_by(Taxon)%>%
  summarise(n=n())%>%
  arrange(-n)

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "aldex2",rare.status,filter.status,agglom.rank,
  paste(custom.levels,collapse = '-'),
  truncationlvl,ref.level,"workspace.RData",sep="-")))

write.table(aldex.signif.features,
            file=file.path("./output/rtables",authorname,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "aldex2",rare.status,
              filter.status,agglom.rank,
              paste(custom.levels,collapse = '-'),truncationlvl,
              ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

q()
