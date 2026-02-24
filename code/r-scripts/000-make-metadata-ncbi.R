library(tidyverse)
authorname<-"pooled" # name of the folder with QIIME2 output
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
metadatadir.16s<-file.path("./data/metadata",
                       paste(authorname,"metadata",sep = "-")) # directory with metadata
metadatadir.wms<-file.path("../metagenome/data/metadata") # directory with metadata


metadata.filename.16s<-file.path(metadatadir.16s,
                             paste("filenames",read.end.type,
                                   authorname,"raw-supercomp.tsv", sep = "-"))
metadata.filename.wms<-file.path(metadatadir.wms,
                             paste("samples.tsv", sep = "-"))
custom.md.16s<-read.table(metadata.filename.16s, header = T)
custom.md.wms<-read.table(metadata.filename.wms, header = F)%>%
  as_tibble
custom.md.wms<-custom.md.wms%>%
  select(V2)%>%
  rename("sample_file_name"="V2")%>%
  mutate(sample_name=gsub("_L[0-9]*_[0-9]*\\.fq\\.gz","",sample_file_name),
         host_subject_id=gsub("_wms","",sample_name))%>%
  distinct(sample_name,.keep_all = TRUE)%>%
  left_join(custom.md.16s,by=join_by("host_subject_id"=="sample.id"))%>%
  select(-sample_file_name,-absolute.filepath)


# Biosample metadata ####
biosample.md<-custom.md.16s%>%
  as_tibble()%>%
  select(-absolute.filepath)%>%
  filter(class%in%c("NMR","B6mouse","FVBNmouse","MSMmouse"))%>%
  add_row(sample.id = "Y29",
          class = "NMR",
          animal = "naked mole rat",
          sex = "F",
          birthday = "2017-12-25")%>%
  arrange(desc(class),sample.id)%>%
  mutate(sample_name=paste(sample.id,"16S",sep="_"))%>%
  rename("host_subject_id"="sample.id")%>%
  relocate(sample_name)%>%
  rbind(custom.md.wms)%>%
  mutate(class_general=case_when(class=="NMR" ~ "NMR",
                                 class=="B6mouse"|
                                   class=="FVBNmouse"|
                                   class=="MSMmouse"~ "mouse"))%>%
  mutate(organism=case_when(class_general=="NMR" ~ "gut metagenome",
                            class_general=="mouse"~ "mouse gut metagenome"))%>%
  mutate(collection_date=case_when(
    host_subject_id == "2D10" & grepl("16S",sample_name) ~ "2022-07-07/2022-11-17",
    host_subject_id == "2D14" & grepl("16S",sample_name) ~ "2022-07-07/2022-11-15",
    host_subject_id == "CG33" & grepl("16S",sample_name) ~ "2022-06-01/2022-10-13",
    host_subject_id == "D27" & grepl("16S",sample_name) ~ "2022-07-29/2022-11-17",
    host_subject_id == "DCG6" & grepl("16S",sample_name) ~ "2022-08-12", # "2022-06-30" in 220809.xlsx !!!
    host_subject_id == "GA17" & grepl("16S",sample_name) ~ "2022-06-21/2022-10-05",
    host_subject_id == "GA5" & grepl("16S",sample_name) ~ "2022-06-29/2022-11-15", 
    host_subject_id == "GH36" & grepl("16S",sample_name) ~ "2022-07-06/2022-11-15",
    host_subject_id == "GH6" & grepl("16S",sample_name)  ~ "2022-07-06/2022-11-15",
    host_subject_id == "L122" & grepl("16S",sample_name) ~ "2022-07-07/2022-11-15",
    host_subject_id == "L128" & grepl("16S",sample_name) ~ "2022-07-28/2022-11-15",
    host_subject_id == "L133" & grepl("16S",sample_name) ~ "2022-07-28/2022-09-29",
    host_subject_id == "M40" & grepl("16S",sample_name) ~ "2022-07-29/2022-11-17",
    host_subject_id == "O15" & grepl("16S",sample_name) ~ "2022-07-27",
    host_subject_id == "R10" & grepl("16S",sample_name) ~ "2022-06-23/2022-11-17",
    host_subject_id == "R35" & grepl("16S",sample_name) ~ "2022-07-29",
    host_subject_id == "Y29" & grepl("16S",sample_name) ~ "2022-07-07/2022-11-15",
    host_subject_id == "Y51b" & grepl("16S",sample_name) ~ "2022-07-27",
    host_subject_id == "Y66b" & grepl("16S",sample_name) ~ "2022-07-27/2022-12-08",
    
    class == "B6mouse" ~ "2022-08-01/2022-09-30",
    
    host_subject_id %in% c("MSM339","MSM340","MSM342") ~ "2023-06-02/2023-06-07",
    host_subject_id %in% c("MSM341", "MSM343", "MSM344", "MSM345") ~ "2023-06-02/2023-06-22",
    host_subject_id %in% c("MSM346", "FVBN549","FVBN550","FVBN551") ~ "2023-06-22/2023-06-23",
    
    host_subject_id == "2D10" & grepl("wms",sample_name) ~ "2023-05-22",
    host_subject_id == "2D14" & grepl("wms",sample_name) ~ "2023-05-22/2023-07-25",
    host_subject_id == "G14" ~ "2023-07-12/2023-07-26",
    host_subject_id == "G18" ~ "2023-05-31/2023-07-25",
    
    host_subject_id == "H3"  ~ "2023-07-26",
    host_subject_id == "H4" ~ "2023-07-25/2023-07-26",
    host_subject_id == "H15" ~ "2023-07-12/2023-07-25",
    host_subject_id == "H21" ~ "2023-05-22/2023-05-31",
    host_subject_id == "O15" & grepl("wms",sample_name) ~ "2023-07-25/2023-07-26",
    host_subject_id == "Y51b" & grepl("wms",sample_name) ~ "2023-07-25",
    host_subject_id == "Y66b" & grepl("wms",sample_name) ~ "2023-06-02/2023-07-25"))%>%
  mutate(env_broad_scale="laboratory environment [ENVO:01001405]",
         env_local_scale="animal-associated environment",
         env_medium="fecal material",
         geo_loc_name=case_when(class== "NMR"| 
                                  class=="B6mouse" ~ "Japan:Kumamoto, Kumamoto",
                                class=="FVBNmouse"|
                                  class=="MSMmouse"~ "Japan:Chiba, Chiba"),
         host = case_when(class_general =="NMR" ~ "Heterocephalus glaber",
                        class_general =="mouse"~ "Mus musculus"),
         host_common_name = case_when(class_general =="NMR" ~ "Naked mole-rat",
                                      class_general =="mouse" ~ "Mouse"),
         lat_lon=case_when(class== "NMR"|class=="B6mouse" ~ "32.793913 N 130.712695 E",
                           class=="FVBNmouse"|
                             class=="MSMmouse"~ "35.58128625 N 140.16152814 E"))%>%
  rename(host_birthday = birthday)%>%
  mutate(host_birthday = ifelse(host_birthday == "-", NA,host_birthday))%>%
  mutate(colony_id = ifelse(class_general =="NMR",
                            gsub("[0-9]*[a-z]*$","",host_subject_id),
                            "not applicable"))%>%
  
  select(-class,-animal,-class_general)%>%
  relocate(sex,.after=last_col())%>%
  rename(host_sex=sex)%>%
  mutate(host_sex=case_when(host_sex == "M" ~ "male",
                            host_sex == "F" ~ "female"))%>%
  mutate(first_collection_date = gsub("\\/[0-9-]*","",collection_date))%>%
  mutate(host_age=year(as.period(interval(host_birthday,first_collection_date))),
         host_age=as.character(host_age))%>%
  mutate(caste = case_when(
    host_subject_id == "2D10" ~ "pair",
    host_subject_id == "2D14"~ "pair",
    host_subject_id == "O15" ~ "pair",
    host_subject_id == "Y29" ~ "pair",
    host_subject_id == "Y51b" ~ "worker",
    host_subject_id == "Y66b" ~ "worker",
    host_subject_id == "H4" ~ "pair",
    host_subject_id == "H21" ~ "pair",
    host_subject_id == "H3" ~ "pair",
    host_subject_id == "G18" ~ "pair",
    host_subject_id == "G14" ~ "pair",
    host_subject_id == "H15" ~ "pair",
    host_subject_id == "CG33" ~ "pair",
    host_subject_id == "D27" ~ "worker",
    host_subject_id == "DCG6" ~ "worker",
    host_subject_id == "GA17" ~ "single",
    host_subject_id == "GA5" ~ "pair",
    host_subject_id == "GH36" ~ "pair",
    host_subject_id == "GH6" ~ "pair",
    host_subject_id == "L122" ~ "pair",
    host_subject_id == "L128" ~ "pair",
    host_subject_id == "L133" ~ "pair",
    host_subject_id == "M40" ~ "pair",
    host_subject_id == "R10" ~ "worker",
    host_subject_id == "R35" ~ "worker"))%>%
  mutate(sequencing_type = case_when(grepl("16S", sample_name) ~ 
                                       paste(host_common_name, "16S rRNA gene sequencing"),
                                     grepl("wms", sample_name) ~
                                       paste(host_common_name, "whole metagenome sequencing")))%>%
  replace(is.na(.),"missing")
  
write.table(biosample.md,file="data/metadata/pooled-metadata/biosample_metadata_for_ncbi.tsv",
            sep = "\t",col.names = T,row.names = F)


# SRA metadata ####
sra.md<-biosample.md%>%
  select(sample_name,sequencing_type)%>%
  mutate(library_ID=sample_name,
         title = sequencing_type,
         library_strategy = case_when(grepl("16S",sample_name) ~ "AMPLICON",
                                      grepl("wms",sample_name) ~ "WGS"),
         library_source = "METAGENOMIC",
         library_selection =  case_when(grepl("16S",sample_name) ~ "PCR",
                                        grepl("wms",sample_name) ~ "RANDOM"), 
         library_layout = "PAIRED",
         platform = "ILLUMINA",
         instrument_model = case_when(grepl("16S",sample_name) ~ "Illumina MiSeq",
                                      grepl("wms",sample_name) ~ "Illumina Hiseq 2500"),
         
         design_description = case_when(grepl("16S",sample_name) ~ "DNA extraction from the feces and 16S rRNA amplicon sequencing were performed by Noster Corporation (Kyoto, Japan) with the QIAamp DNA Microbiome Kit. For the paired-end sequencing with the MiSeq system (Illumina), the MiSeq Reagent Kit v3 and PhiX Control Kit v3 (Illumina) were used",
                                        grepl("wms",sample_name) ~ "DNA was extracted with the QIAmpFast DNA Stool Mini Kit (Qiagen), DNA concentration was measured using NanoDrop, and sequencing was performed with Illumina Hiseq 2500 by Novogene (Japan)"),
         filetype = "fastq",
         filename = case_when(
           library_strategy == "AMPLICON" ~ paste(sample_name,"R1.fastq.gz", sep = "_"),
           library_strategy == "WGS" ~ paste(sample_name,"R1.fq.gz", sep = "_")),
         filename2 = case_when(
           library_strategy == "AMPLICON" ~ paste(sample_name,"R2.fastq.gz", sep = "_"),
           library_strategy == "WGS" ~ paste(sample_name,"R2.fq.gz", sep = "_")))%>%
  select(-sequencing_type)



write.table(sra.md,file="data/metadata/pooled-metadata/sra_metadata_for_ncbi.tsv",
            sep = "\t",col.names = T,row.names = F)
