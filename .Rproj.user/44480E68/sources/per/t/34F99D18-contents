library(tidyverse)
library(phyloseq)
library(ALDEx2)
library(vegan)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"
# load("./rdafiles/yashliu-dada2-284-203-Genus-138-1-phyloseq-workspace.RData")

load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))


custom.levels<-c("NMR","SPFmouse","spalax","FukomysDamarensis","rabbit")

# Import data ####
rare.status<-"rarefied"
rare.status<-"nonrarefied"
if (rare.status=="rarefied"){
  ps.q.df.aldex.input<-read.table(paste0("./rtables/alldir/ps.q.df.rare.filtered-",
                                           paste(custom.levels,collapse = '-'),".tsv"),
                                    header = T)
} else if(rare.status=="nonrarefied"){
  ps.q.df.aldex.input<- read.table(paste0("./rtables/alldir/ps.q.df.norare.filtered-",
                                            paste(custom.levels,collapse = '-'),".tsv"),
                                     header = T)
}

# convert the data frames into wide format
ps.q.df.aldex.input.wide<-ps.q.df.aldex.input%>%
  dplyr::select(Sample,Abundance,Taxon)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()


# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.aldex.input.wide)<-ps.q.df.aldex.input.wide$Sample
ps.q.df.aldex.input.wide<-ps.q.df.aldex.input.wide[,-1]  
ps.q.df.aldex.input.wide<-t(ps.q.df.aldex.input.wide)

covariates<-custom.md$class[match(colnames(ps.q.df.aldex.input.wide),rownames(custom.md))]
mm <- model.matrix(~ covariates-1)
# reorder model.matrix to put NMR as first column
# custom.levels[c(which(custom.levels == "NMR"), which(custom.levels != "NMR"))]
mm<-mm[, colnames(mm)[c(which(colnames(mm) == "covariatesNMR"), which(colnames(mm) != "covariatesNMR"))]]


# aldex glm for a complex case ####
set.seed(1)
ps.q.aldex.clr <- aldex.clr(ps.q.df.aldex.input.wide, mm, mc.samples=1000, denom="all", verbose=F)
ps.q.glm.test <- aldex.glm(ps.q.aldex.clr,mm)
ps.q.glm.effect <- aldex.glm.effect(ps.q.aldex.clr, CI= T)


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

save.image(paste0("./rdafiles/aldex2-",rare.status,"-","filtered-",
                  paste(custom.levels,collapse = '-'),
                  "-workspace.RData"))
write.table(aldex.signif.features,
            file=paste0("./rtables/alldir/aldex2-",rare.status,"-","filtered-",
                        paste(custom.levels,collapse = '-'),"-",
                        "nmr-signif",".tsv"), # no rare.status
            row.names = F,sep = "\t")
