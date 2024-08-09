library(tidyverse)
library(phyloseq)
library(Maaslin2)
# library(vegan)
# Import custom.md, ps.q.df.wide, custom.levels, ref.level, comparison
# input_data_date_time is Rdata workspace from diffabund-input.R
input_data_date_time<-"20240613_21_16_08"
# input_data_date_time<-"20240613_21_16_08" for nmr ages, OTU level
# 20240615_17_33_25 for nmr sexes, OTU level
authorname<-"pooled"
# choose what to compare
comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"
# choose the host of interest
host<-"NMR"
# host<-"mice"
ref.level<-"agegroup0_10" # choose the reference level
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}

truncationlvl<-"234"
agglom.rank<-"OTU"
read.end.type<-"single"

rare.status<-"rare"
filter.status<-"nonfiltered"
custom.levels<-c("agegroup0_10",
                                "agegroup10_16")
# custom.levels<-c("F",
#                  "M")

load(file.path("./output/rdafiles",paste(
  input_data_date_time,
  "diffabund-input-data",host,rare.status,filter.status,agglom.rank,
  comparison,truncationlvl,
  paste(sort(custom.levels),collapse = '-'),"workspace.RData",sep="-")))

# 3.2 Running MaAsLin 2 ####
if (comparison=="age"){
  maaslin.reference<-paste("agegroup",ref.level,sep = ",")
  maaslin.comparison<-"agegroup"
}else if(comparison=="sex"){
  maaslin.reference<-paste("sex","F",sep = ",")
  maaslin.comparison<-"sex"
}else if(comparison=="strain"){
  maaslin.reference<-paste("class",ref.level,sep = ",")
  maaslin.comparison<-"class"
}

relations<-read.table("./data/metadata/pooled-metadata/nmr-relations.tsv",
                      header = T,
                      sep = "\t")
custom.md<-custom.md%>%
  left_join(relations,by="Sample")
rownames(custom.md)<-custom.md$Sample


set.seed(1)
maaslin.fit_data = 
  Maaslin2(input_data = ps.q.df.wide, 
           input_metadata = custom.md, 
           min_prevalence = 0,
           normalization = "TSS",
           transform = "LOG",
           analysis_method = "LM",
           random_effects = c("relation"), 
           standardize = FALSE,
           output = file.path("./output/maaslin2",paste0(authorname,"-output"),
                              rare.status,paste(
                                paste(format(Sys.time(),format="%Y%m%d"),
                                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                host,filter.status,agglom.rank,
                                comparison,truncationlvl,
                                paste(custom.levels,collapse = '-'),
                                "ref",ref.level,sep = "-")), 
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.05)


save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2",host,rare.status,filter.status,agglom.rank,
  comparison,truncationlvl,
  paste(custom.levels,collapse = '-'),"ref",
  ref.level,"workspace.RData",sep="-")))
q()
