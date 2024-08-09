# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("phyloseq","Maaslin2"))
# install.packages(c("tidyverse","Polychrome"))
library(phyloseq)
library(tidyverse)
library(Maaslin2)
# library(vegan)
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
phyloseq.workspace.date_time<-"20240426_21_44_30"

load(file.path("./output/rdafiles",paste(
  phyloseq.workspace.date_time,
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
                 "pvo"
                 # ,
                 # "NMRwt"
                 )
ref.level<-"NMR"
# Import data ####
rare.status<-"rare"
filter.status<-"nonfiltered"
ps.q.df.maaslin.input<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      "20240426_22_00_04",
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)

ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)
# Maaslin2 ####
# Input data and metadata ####
# convert the data frames into wide format
ps.q.df.maaslin.input.wide<-ps.q.df.maaslin.input%>%
  filter(class%in%custom.levels,Abundance!=0)%>%
  dplyr::select(Sample,Abundance,all_of(agglom.rank))%>%
  pivot_wider(names_from = all_of(agglom.rank),
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
           analysis_method = "LM",
           normalization = "TSS",
           transform = "LOG",
           random_effects = NULL,
           standardize = FALSE,
           output = file.path("./output/maaslin2",paste0(authorname,"-output"),
                           rare.status,paste(
                             paste(format(Sys.time(),format="%Y%m%d"),
                                   format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                             paste(custom.levels,collapse = '-'),
                           truncationlvl,"ref",ref.level,sep = "-")), 
           fixed_effects = c("class"),
           reference = paste0("class,",ref.level),
           max_significance = 0.05)
save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2",rare.status,filter.status,agglom.rank,
                        paste(custom.levels,collapse = '-'),
                        truncationlvl,ref.level,"workspace.RData",sep="-")))
q()