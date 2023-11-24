library(tidyverse)
library(phyloseq)
library(vegan)
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
custom.levels<-c("NMR",
                 "B6mouse",
                 "MSMmouse",
                 "FVBNmouse",
                 "DMR",
                 "hare",
                 "rabbit",
                 "spalax",
                 "pvo")
# Import data ####
rare.status<-"rare"
filter.status<-"nonfiltered"

ancom.dir<-"./ANCOM-II"
source(file.path(ancom.dir,"rel_ancom.R"))

ps.q.df.ancom_ii.input<-read.table(paste0("./rtables/pooled/ps.q.df.",
                                         rare.status,".",filter.status,"-",agglom.rank,"-",
                                         paste(custom.levels,collapse = '-'),".tsv"),
                                  header = T)
# ANCOM-II ####
# Input data and metadata ####
# convert the data frames into wide format
custom.md.foo<-custom.md
custom.md<-custom.md.foo

ps.q.df.ancom_ii.input.wide<-ps.q.df.ancom_ii.input%>%
  dplyr::select(Sample,Abundance,Taxon)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()%>%
  rename(Sample.ID=Sample)



custom.md$Sample.ID<-rownames(custom.md)
ps.q.df.ancom_ii.input.wide<-
  ps.q.df.ancom_ii.input.wide[match(custom.md$Sample.ID,ps.q.df.ancom_ii.input.wide$Sample.ID),]

rownames(custom.md)<-NULL
rownames(ps.q.df.ancom_ii.input.wide)<-NULL
custom.md<-custom.md%>%
  select(Sample.ID,class)

ps.q.df.ancom_ii.input.wide$Sample.ID<-custom.md$Sample.ID<-seq_along(custom.md$Sample.ID)



# Example 1. using logarithm of geometric mean as normalizing quantity.
ancom.out<-rel_ancom(OTUdat = ps.q.df.ancom_ii.input.wide,Vardat = custom.md,main.var = "class")
ancom.out$detected_microbes #a list of all microbes that are detected as being differentially abundant across
# the population of interest
ancom.out$structural_zeros
# Example 2. using logarithm of "Muribaculaceae (Muribaculaceae)" as normalizing quantity.
rel_ancom(OTUdat=ps.q.df.ancom_ii.input.wide, Vardat=custom.md,main.var="class",
          ref_name="Muribaculaceae..Muribaculaceae." )
