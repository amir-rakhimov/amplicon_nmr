library(dplyr)
ancom.dir<-"./ANCOM-II"
otu<-read.csv(file.path(ancom.dir,"otu.csv"),dec = ",")
meta<-read.csv(file.path(ancom.dir,"meta.csv"))
source(file.path(ancom.dir,"rel_ancom.R"))
# Example 1. using logarithm of geometric mean as normalizing quantity.
ancom.out<-rel_ancom(OTUdat = otu,Vardat = meta,main.var = "country")
ancom.out$detected_microbes #a list of all microbes that are detected as being differentially abundant across
# the population of interest
ancom.out$structural_zeros
# Example 2. using logarithm of "p__Actinobacteria.g__Bifidobacterium" as normalizing quantity.
rel_ancom(OTUdat=otu, Vardat=meta,main.var="country",
          ref_name="p__Actinobacteria.g__Bifidobacterium" )
