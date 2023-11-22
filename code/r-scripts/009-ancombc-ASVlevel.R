library(eeptools)
library(tidyverse)
library(phyloseq)
library(vegan)
library(ANCOMBC)
truncationlvl<-"234"
agglom.rank<-"OTU"
read.end.type<-"single"
authorname<-"pooled"
load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))


rare.status<-"rare"
filter.status<-"nonfiltered"

host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
comparison<-"sex"
comparison<-"strain"

host.labels<-
  c("NMR" = "*Heterocephalus glaber*")

host.labels<-
  c("B6mouse" = "B6 mouse",
    "MSMmouse" = "MSM/Ms mouse",
    "FVBNmouse" = "FVB/N mouse")

ref.level<-"FVBNmouse"
ps.q.df.preprocessed<-read.table(paste0("./rtables/",authorname,"/ps.q.df.",
                                        rare.status,".",filter.status,"-",agglom.rank,"-",
                                        paste(names(host.labels),collapse = '-'),".tsv"),
                                 header = T,sep = "\t")

if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class=="NMR",Abundance!=0)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=age_calc(birthday,units = "years"))%>%
    mutate(age=round(age))
  
  min_boundary <- floor(min(ps.q.df.preprocessed$age)/5) * 5
  max_boundary <- ceiling(max(ps.q.df.preprocessed$age)/5) * 5
  # create a column with age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(age_group=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                         include.lowest = TRUE))
  
  # Metadata
  custom.md$Sample<-rownames(custom.md)
  custom.md<-custom.md%>% 
    filter(class=="NMR")%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=age_calc(birthday,units = "years"))%>%
    mutate(age=round(age))
  min_boundary <- floor(min(custom.md$age)/5) * 5
  max_boundary <- ceiling(max(custom.md$age)/5) * 5
  
  custom.md<-custom.md%>%
    mutate(age_group=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                        include.lowest = TRUE))%>%
    as.data.frame()
  rownames(custom.md)<-custom.md$Sample
}else if(host=="mice"){
  # select mice and add age groups
  custom.levels<-c("B6mouse",
                   "MSMmouse",
                   "FVBNmouse")
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,Abundance!=0)
  ps.q.df.preprocessed$age_group<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
  # also to metadata
  custom.md<-custom.md%>%
    filter(class%in%custom.levels)
  custom.md$age_group<-ifelse(custom.md$class=="B6mouse","B6",
                             ifelse(grepl("2020",custom.md$birthday),"old","young"))
}


if (comparison=="age"){
  pretty.facet.labels<-names(table(ps.q.df.preprocessed$age_group))
  names(pretty.facet.labels)<-names(table(ps.q.df.preprocessed$age_group))
  custom.levels<-names(pretty.facet.labels)
  
}else if (comparison=="sex"){
  pretty.facet.labels<-
    c("F" = "Females",
      "M" = "Males")
  custom.levels<-names(pretty.facet.labels)
  
}else if(comparison=="strain"){
  pretty.facet.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
  pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]
}
ps.q.df <-ps.q.df.preprocessed%>%
  dplyr::select(Sample,OTU,Abundance,class,age_group,sex)%>%
  filter(Abundance!=0)


# ANCOMBC ####
# Input data and metadata ####
# convert the data frames into wide format
ps.q.df.ancombc.input.wide<-ps.q.df%>%
  dplyr::select(Sample,Abundance,OTU)%>%
  pivot_wider(names_from = "OTU", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.ancombc.input.wide)<-ps.q.df.ancombc.input.wide$Sample
ps.q.df.ancombc.input.wide<-ps.q.df.ancombc.input.wide[,-1] 

taxmat<-ps.q.agg%>%
  dplyr::select(Kingdom,Phylum,Class,Order,Family,Genus,OTU)%>%
  distinct()%>%
  column_to_rownames(var = "OTU")%>%
  as.matrix()
ps.q.OTU<-t(ps.q.df.ancombc.input.wide)
ps.q.OTU<-otu_table(ps.q.OTU,taxa_are_rows = T)
ps.q.TAX<-tax_table(taxmat)
ps.q.phyloseq.new<-phyloseq(otu_table(ps.q.OTU),
                            tax_table(ps.q.TAX),
                            sample_data(custom.md))

# relevel our comparison vector. The first level will be the reference
# for custom leveling
ancombc.levels<-c(ref.level,custom.levels[custom.levels!=ref.level])
if (comparison=="age"){
  sample_data(ps.q.phyloseq.new)$age_group<-factor(sample_data(ps.q.phyloseq.new)$age_group,
                                                   levels = ancombc.levels)
  ancombc.comparison<-"age_group"
}else if(comparison=="sex"){
  sample_data(ps.q.phyloseq.new)$sex<-factor(sample_data(ps.q.phyloseq.new)$sex,
                                                   levels = ancombc.levels)
  ancombc.comparison<-"sex"
}else if(comparison=="strain"){
  sample_data(ps.q.phyloseq.new)$class<-factor(sample_data(ps.q.phyloseq.new)$class,
                                                   levels = ancombc.levels)
  ancombc.comparison<-"class"
}

# perform differential abundance test
ancombc.out<-ancombc(
  phyloseq = ps.q.phyloseq.new,
  # tax_level = ,
  formula = ancombc.comparison,
  p_adj_method = "fdr", 
  prv_cut = 0, # by default prevalence filter of 10% is applied
  lib_cut = 0, 
  group = ancombc.comparison,
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
ancombc.res<-ancombc.out$res

# find differentially abundant taxa by multiplying fold change with TRUE/FALSE
# for diff abund
df_lfc = data.frame(ancombc.res$lfc[, -1] * ancombc.res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancombc.res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
ancombc.signif.features<-subset(df_lfc,rowSums(df_lfc[,-c(1,2)]!=0)==ncol(df_lfc[,-c(1,2)]))

save.image(file.path("./rdafiles",
                  paste("ancombc",rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,ref.level,
                        "workspace.RData",sep="-")))
write.table(ancombc.signif.features,
            file=paste0("./rtables/",authorname,"/",
                        paste("ancombc",rare.status,filter.status,host,agglom.rank,
                              comparison,truncationlvl,ref.level,
                              "signif.tsv",sep="-")), 
            row.names = F,sep = "\t")
