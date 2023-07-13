library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"

load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

custom.levels<-c("NMR","SPFmouse")
custom.levels<-c("NMR","SPFmouse","spalax","FukomysDamarensis","rabbit")

# Import data ####
rare.status<-"rarefied"
rare.status<-"nonrarefied"
if (rare.status=="rarefied"){
  ps.q.df.maaslin.input<-read.table(paste0("./rtables/alldir/ps.q.df.rare.filtered-",
                                           paste(custom.levels,collapse = '-'),".tsv"),
                                    header = T)
} else if(rare.status=="nonrarefied"){
  ps.q.df.maaslin.input<- read.table(paste0("./rtables/alldir/ps.q.df.norare.filtered-",
                                            paste(custom.levels,collapse = '-'),".tsv"),
                                     header = T)
}


# Maaslin2 ####
# Input data and metadata ####
# convert the data frames into wide format
ps.q.df.maaslin.input.wide<-ps.q.df.maaslin.input%>%
  select(Sample,Abundance,Taxon)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()


# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.maaslin.input.wide)<-ps.q.df.maaslin.input.wide$Sample
ps.q.df.maaslin.input.wide<-ps.q.df.maaslin.input.wide[,-1]  


# 3.2 Running MaAsLin 2 ####
set.seed(1)
maaslin.fit_data = 
  Maaslin2(input_data = ps.q.df.maaslin.input.wide, 
           input_metadata = custom.md, 
           min_prevalence = 0,
           normalization = "TSS",
           transform = "AST",
           analysis_method = "LM",
           random_effects = NULL,
           standardize = FALSE,
           output = paste0("./output/maaslin2/alldir-output/",
                           rare.status,"/",paste(custom.levels,collapse = '-')), 
           fixed_effects = c("class"),
           reference = c("class,NMR"),
           max_significance = 0.05)

save.image(paste0("./rdafiles/maaslin-",rare.status,"-","filtered-",
                  paste(custom.levels,collapse = '-'),
                  "-workspace.RData"))
