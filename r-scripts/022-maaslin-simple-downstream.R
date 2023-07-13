library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"
source("r-scripts/make_features_pretty.R")

custom.levels<-c("NMR","SPFmouse")
rare.status<-"rarefied"
rare.status<-"nonrarefied"
load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
load(paste0("./rdafiles/maaslin-",rare.status,"-","filtered-",
                  paste(custom.levels,collapse = '-'),
                  "-workspace.RData"))


maaslin.fit.df<-maaslin.fit_data$results%>%
  filter(qval<0.05)

maaslin.signif.decreased<-maaslin.fit.df%>%
  as_tibble()%>%
  filter(coef<0)%>%
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

maaslin.signif.increased<-maaslin.fit.df%>%
  as_tibble()%>%
  filter(coef>0)%>%
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

maaslin.signif.decreased<-make_features_pretty(maaslin.signif.decreased,"feature")
maaslin.signif.increased<-make_features_pretty(maaslin.signif.increased,"feature")

table(maaslin.signif.decreased$feature%in%ps.q.agg.rel$Taxon)
table(maaslin.signif.increased$feature%in%ps.q.agg.rel$Taxon)
write.table(maaslin.signif.decreased,
            file=paste0("./rtables/alldir/maaslin-",
                        paste(custom.levels,collapse = '-'),"-",
                        "nmr-signif-",rare.status,".tsv"),
            row.names = F,sep = "\t")
