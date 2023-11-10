library(PERFect)
library(ggplot2)
library(knitr)
library(kableExtra)
library(tidyverse)
# Workflow

# The filtering methods for this package are wrapped into two main functions,
# *PERFect_sim()* (performing simultaneous filtering) and *PERFect_perm()* 
# (performing permutation filtering). First, we load the OTU matrix with 240 samples and 46 taxa.
split.ec.levels<-function(ec.df){
  ec.levels<-data.frame(do.call('rbind', 
                                strsplit(as.character(ec.df$predicted.function),'.',
                                         fixed=TRUE)))
  colnames(ec.levels)<-paste0("EC.level.",1:ncol(ec.levels))
  ec.levels$EC.level.1<-gsub("EC\\:","",ec.levels$EC.level.1)
  ec.df.with.lvls<-cbind(ec.df,ec.levels)
  return(ec.df.with.lvls)
}
# Load metadata as a tibble ####
metadata <- 
  read.table(paste0("./data/pooled-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  # filter(class=="NMR"|class=="SPFmouse")
  filter(class=="NMR")
# filter(class=="SPFmouse")

# Load KO abundance as a data.frame ####
# data(ko_abundance)
enz_abundance <- 
  read.table("./data/pooled-data/picrust/supercomp/NMR_SPFmouse/pred_metagenome_unstrat.tsv",
                            header = T,sep="\t")
colnames(enz_abundance)[which(colnames(enz_abundance)=="function.")]<-"predicted.function"
colnames(enz_abundance)<-gsub("X2D10", "2D10",colnames(enz_abundance))
colnames(enz_abundance)<-gsub("X2D14", "2D14",colnames(enz_abundance))
# enz_abundance<-enz_abundance%>%
#   column_to_rownames("predicted.function")

enzymes.long<-enz_abundance%>%
  pivot_longer(!predicted.function,names_to = "Sample",values_to = "Abundance")

enzymes.with.EC<-split.ec.levels(enzymes.long)

filtered.enzymes.long<-enzymes.with.EC%>%
  inner_join(metadata[,c("sample.id","class")],by=c("Sample"="sample.id"))%>%
  mutate(Abundance=as.numeric(Abundance))%>%
  # filter(Sample=="DCG6",EC.level.1==1,EC.level.2==1,EC.level.3==1)%>%
  # summarise(sum=sum(Abundance))
  # group_by(predicted.function)%>%
  group_by(Sample,EC.level.1,EC.level.2,EC.level.3)%>%
  summarise(Abundance=sum(Abundance))%>%
  filter(EC.level.1==1 | EC.level.1==3)%>%
  unite(predicted.function,c(EC.level.1,EC.level.2,EC.level.3),sep = ".")

filtered.enzymes.long$predicted.function<-
  paste0("EC:",filtered.enzymes.long$predicted.function)
filtered.enz_abundance<-filtered.enzymes.long%>%
  pivot_wider(names_from = Sample,values_from = Abundance,values_fill = 0)%>%
  column_to_rownames("predicted.function")%>%
  t()


# By default, the function *PERFect_sim()* takes the data table, $X$, as a 
# matrix or data frame, orders it by the taxa abundance, uses 10%, 25% and 50% 
# quantiles for matching the log of DFL to a Skew-Normal distribution and then 
# calculates the p-value for each taxon at the significance level of
# $\alpha$ = 0.1. The function *PERFect_sim()* only needs a taxa table as the 
# input, and other parameters are set to default.


enz.res_sim <- PERFect_sim(X = filtered.enz_abundance,alpha = 0.4)

# The "filtX" object from the result stores the filtered OTU matrix, which 
# consists of 240 samples and 10 taxa in this example. We can identify the
# remaining taxa by looking into the column of the OTU matrix. 

dim(enz.res_sim$filtX)      
colnames(enz.res_sim$filtX) # signal taxa


# The p-values for all taxa can be extracted as
head(enz.res_sim$pvals)
enz.res_sim[["pvals"]][which(enz.res_sim[["pvals"]]<0.4)]


# and they can be plot using the function *pvals_Plots()*.
p <- pvals_Plots(PERFect = enz.res_sim, X = filtered.enz_abundance, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.4)
p$plot + ggtitle("Simultanenous Filtering")


# Alternatively, we can use permutation filtering *PERFect_perm()* which is 
# more robust than simultaneous filtering. By default, this function generates
# k = 10000 permutations for each taxa, thus it can be computationally expensive 
# for a large OTU matrix. We offer user a fast algorithm which employs an 
# unbalanced binary search that optimally finds the cutoff taxon without
# building the permutation distribution for all taxa. The codes for these 
# options are shown below.


set.seed(1)
enz.res_perm <- PERFect_perm(X = filtered.enz_abundance, 
                             Order = "pvals", 
                             pvals_sim = enz.res_sim, 
                             algorithm = "full",
                             alpha = 0.4)
enz.res_perm[["pvals"]][which(enz.res_perm[["pvals"]]<0.4)]

set.seed(1)
enz.res_perm2 <- PERFect_perm(X = filtered.enz_abundance, 
                              Order = "pvals", 
                              pvals_sim = enz.res_sim, 
                              algorithm = "fast",
                              alpha = 0.4)
enz.res_perm2[["pvals"]][which(enz.res_perm2[["pvals"]]<0.4)]

p1 <- pvals_Plots(enz.res_perm, filtered.enz_abundance)
p1 <- p1$plot + ggtitle("Full Algorithm")
p2 <- pvals_Plots(enz.res_perm2, filtered.enz_abundance)
p2 <- p2$plot + ggtitle("Fast Algorithm")
ggpubr::ggarrange(p1,p2,ncol = 2)

