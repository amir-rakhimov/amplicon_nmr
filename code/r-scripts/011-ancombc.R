# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("ANCOMBC")
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse"))
library(tidyverse)
library(ANCOMBC)
library(phyloseq)
# Import custom.md, ps.q.df.wide, custom.levels
# input_data_date_time is Rdata workspace from diffabund-input.R
# with a rarefied table in wide format and metadata
input_data_date_time<-"20240809_15_38_22"
# 20240809_14_40_39 for NMR, OTU level
# 20240809_14_52_10 for all hosts, genus level
# 20240809_15_38_22 for NMR, genus level
ps.q.agg.date_time<-"20240620_12_40_41"
# 20240620_12_40_41 for all hosts, genus level
# 20240620_12_38_18 for all hosts, OTU level
authorname<-"pooled"
inside.host<-TRUE
if(inside.host=="TRUE"){
  # choose what to compare
  comparison<-"age"
  # comparison<-"sex"
  # comparison<-"strain"
  ref.level<-"agegroup0_10" # choose the reference level
  custom.levels<-"NMR"
  if(comparison=="age"){
    sample.groups<-c("agegroup0_10",
                     "agegroup10_16")
  }else if(comparison=="sex"){
    sample.groups<-c("F",
                     "M")
  }
}else{
  comparison<-"host"
  ref.level<-"NMR"
  custom.levels<-c("NMR",
                   "B6mouse",
                   "MSMmouse",
                   "FVBNmouse",
                   "DMR",
                   "hare",
                   "rabbit",
                   "spalax",
                   "pvo")
  
}

truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
# Load the workspace from 007-diffabund-input.R
load(file.path("./output/rdafiles",paste(
  input_data_date_time,
  "diffabund-input-data",rare.status,filter.status,agglom.rank,
  truncationlvl,
  paste(custom.levels,collapse = '-'),"workspace.RData",sep="-")))

# Load the ps.q.agg object from 001-phyloseq-qiime2.R
ps.q.agg<-readRDS(file=file.path("./output/rdafiles",paste(
         ps.q.agg.date_time,
          "phyloseq-qiime",authorname,agglom.rank,read.end.type,truncationlvl,
          "table.rds",sep="-")))


# ANCOMBC ####
### Extract a taxonomic table ####
if(agglom.rank=="OTU"){
  taxmat<-ps.q.agg%>%
    dplyr::select(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>%
    distinct()%>%
    column_to_rownames(var = "OTU")%>%
    as.matrix()
}else{
  all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
  agglom.rank.index<-match(agglom.rank,all.ranks)
  custom.ranks<-all.ranks[1:agglom.rank.index]
  
  taxmat<-ps.q.agg%>%
    ungroup()%>%
    dplyr::select(all_of(custom.ranks))%>%
    distinct()%>%
    column_to_rownames(var = all_of(agglom.rank))%>%
    as.matrix()
}
# Create a phyloseq object for ANCOMBC
ps.q.OTU<-t(ps.q.df.wide)
ps.q.OTU<-otu_table(ps.q.OTU,taxa_are_rows = T)
ps.q.TAX<-tax_table(taxmat)
ps.q.phyloseq.new<-phyloseq(otu_table(ps.q.OTU),
                            tax_table(ps.q.TAX),
                            sample_data(custom.md))

# relevel our comparison vector. The first level will be the reference
# for custom leveling
ancombc.levels<-c(ref.level,sample.groups[sample.groups!=ref.level])
if(comparison=="host"){
  sample_data(ps.q.phyloseq.new)$class<-factor(sample_data(ps.q.phyloseq.new)$class,
                                               levels = ancombc.levels)
  ancombc.comparison<-"class"                                           
}else if (comparison=="age"){
  sample_data(ps.q.phyloseq.new)$agegroup<-factor(sample_data(ps.q.phyloseq.new)$agegroup,
                                                   levels = ancombc.levels)
  ancombc.comparison<-"agegroup"
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

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "ancombc",paste(custom.levels,collapse = '-'),
  agglom.rank,
  comparison,truncationlvl,"ref",
  ref.level,"workspace.RData",sep="-")))
write.table(ancombc.signif.features,
            file=file.path("./rtables",authorname,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "ancombc",paste(custom.levels,collapse = '-'),
              agglom.rank,
              comparison,truncationlvl,"ref",ref.level,
              "signif.tsv",sep="-")), 
            row.names = F,sep = "\t")


