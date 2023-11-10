library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"

load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))

custom.levels<-c("NMR",
                 "B6mouse",
                 # "MSMmouse",
                 # "FVBNmouse",
                 "DMR",
                 "hare",
                 "rabbit",
                 "spalax",
                 "pvo"#,
                 # "NMRwt"
                 )
ref.level<-"NMR"
# Import data ####
rare.status<-"rare"
filter.status<-"nonfiltered"
ps.q.df.maaslin.input<-read.table(paste0("./rtables/",authorname,"/ps.q.df.",
                                         rare.status,".",filter.status,"-",agglom.rank,"-",
                                         paste(custom.levels,collapse = '-'),".tsv"),
                                  header = T)
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)
ps.q.total<-ps.q.total%>%
  filter(Sample%in%rownames(custom.md))
ps.q.1pc<-ps.q.1pc%>%
  filter(class%in%custom.levels)
# Maaslin2 ####
# Input data and metadata ####
# convert the data frames into wide format
ps.q.df.maaslin.input.wide<-ps.q.df.maaslin.input%>%
  filter(class%in%custom.levels,Abundance!=0)%>%
  dplyr::select(Sample,Abundance,Taxon)%>%
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
           output = paste0("./output/maaslin2/",authorname,"-output/",
                           rare.status,"/",paste(paste(custom.levels,collapse = '-'),
                           truncationlvl,"ref",ref.level,sep = "-")), 
           fixed_effects = c("class"),
           reference = paste0("class,",ref.level),
           max_significance = 0.05)

save.image(paste0("./rdafiles/",
                  paste("maaslin",rare.status,filter.status,agglom.rank,
                        paste(custom.levels,collapse = '-'),
                        truncationlvl,ref.level,"workspace.RData",sep="-")))
q()