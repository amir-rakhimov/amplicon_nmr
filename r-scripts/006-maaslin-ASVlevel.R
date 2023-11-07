library(eeptools)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
# choose what to compare
# comparison<-"age"
# comparison<-"sex"
comparison<-"strain"
# choose the host of interest
# host<-"NMR"
host<-"mice"
ref.level<-"FVBNmouse" # choose the reference level
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
load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
rare.status<-"rare"
filter.status<-"nonfiltered"

# Import data ####
ps.q.df.preprocessed<-read.table(paste0("./rtables/alldir/ps.q.df.",
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
  
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(agegroup=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                        include.lowest = TRUE))
  # we create these new levels because maaslin is itsy bitsy
  unique_levels <- ps.q.df.preprocessed %>%
    ungroup()%>%
    distinct(agegroup)%>%
    arrange(agegroup) %>%
    mutate(new_agegroup = paste0("agegroup", row_number()))
  ps.q.df.preprocessed <- ps.q.df.preprocessed %>%
    left_join(unique_levels, by = "agegroup")
  colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="agegroup")]<-"old_agegroup"
  colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="new_agegroup")]<-"agegroup"
  
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
    mutate(agegroup=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                        include.lowest = TRUE))%>%
    as.data.frame()
  # Extract unique levels from the original_vector
  unique_levels <- custom.md %>%
    distinct(agegroup) %>%
    arrange(agegroup) %>%
    mutate(new_agegroup = paste0("agegroup", row_number()))
  custom.md <- custom.md %>%
    left_join(unique_levels, by = "agegroup")
  colnames(custom.md)[which(colnames(custom.md)=="agegroup")]<-"old_agegroup"
  colnames(custom.md)[which(colnames(custom.md)=="new_agegroup")]<-"agegroup"
  rownames(custom.md)<-custom.md$Sample
}else if(host=="mice"){
  # select mice and add age groups
  custom.levels<-c("B6mouse",
                   "MSMmouse",
                   "FVBNmouse")
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,Abundance!=0)
  ps.q.df.preprocessed$agegroup<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
  # also to metadata
  custom.md<-custom.md%>%
    filter(class%in%custom.levels)
  custom.md$agegroup<-ifelse(custom.md$class=="B6mouse","B6",
                             ifelse(grepl("2020",custom.md$birthday),"old","young"))
}


if (comparison=="age"){
  pretty.facet.labels<-names(table(ps.q.df.preprocessed$agegroup))
  names(pretty.facet.labels)<-names(table(ps.q.df.preprocessed$agegroup))
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
  dplyr::select(Sample,OTU,Abundance,class,agegroup,sex)%>%
  filter(Abundance!=0)
ps.q.df.maaslin.input.wide<-ps.q.df%>%
  pivot_wider(names_from = "OTU", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.maaslin.input.wide)<-ps.q.df.maaslin.input.wide$Sample
ps.q.df.maaslin.input.wide<-ps.q.df.maaslin.input.wide[,-c(1:4)]  


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
                          rare.status,"/",paste(host,filter.status,agglom.rank,
                                                comparison,truncationlvl,
                                                ref.level,sep="-")), 
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.05)


save.image(paste0("./rdafiles/",
                  paste("maaslin",rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,ref.level,
                  "workspace.RData",sep="-")))
q()
