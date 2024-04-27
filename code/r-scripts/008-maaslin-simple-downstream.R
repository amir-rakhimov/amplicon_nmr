# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("phyloseq","Maaslin2"))
# install.packages(c("tidyverse","Polychrome"))
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)

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


maaslin.signif.features<-maaslin.fit_data$results%>%
  filter(qval<0.05)

maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

maaslin.signif.increased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef>0)%>%
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

maaslin.signif.decreased<-make_features_pretty(maaslin.signif.decreased,"feature")
maaslin.signif.increased<-make_features_pretty(maaslin.signif.increased,"feature")

table(maaslin.signif.decreased$feature%in%ps.q.agg$Taxon)
table(maaslin.signif.increased$feature%in%ps.q.agg$Taxon)

write.table(maaslin.signif.decreased,
            file=file.path("./rtables/pooled",
                           paste("maaslin",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 "nmr-signif.tsv",sep="-")),
            row.names = F,sep = "\t")