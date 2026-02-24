#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---


#' ```{r, setup 001-phyloseq-qiime2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("code/r-scripts/001-phyloseq-qiime2.R", knit = FALSE)
#' #file.rename("code/r-scripts/001-phyloseq-qiime2.Rmd", 
#' # "markdown/001-phyloseq-qiime2.Rmd")
#' #rmarkdown::render('./markdown/001-phyloseq-qiime2.Rmd',
#' # 'html_document',
#' # knit_root_dir="/home/rakhimov/projects/amplicon_nmr/")
#' ```
#' 
#+ echo=FALSE
# Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank ####
#' # Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank
#'  
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' Once you produce the feature table, taxonomic classification, 
#' and the phylogenetic tree in QIIME2, it's time to perform 
#' downstream processing in R. First, we need to import the 
#' QZA files using `qiime2R` package.
#' We will convert the QZA files directly into phyloseq objects.
#'
#' The final output is the dataframe ps.q.agg and metadata custom.md.
#' ps.q.agg and custom.md are saved as tsv files and rds files.
#'
#' 1) ps.q.agg is an ASV table with 7-13 columns   
#' (7 if agglomerating at Phylum, 13 if agglomerating at ASV level):  
#' * `Sample`: samples that were sequenced.   
#' * `Abundance`: Absolute abundance of taxa.   
#' * `class`: short names of animal hosts. The variable is factor 
#' with 9 levels at most (B6mouse, DMR, FVBNmouse, hare, 
#' MSMmouse, NMR, pvo, rabbit, spalax).   
#' The next seven columns may not all be in the table. 
#' If you agglomerate by Genus, you don't see the Species column. 
#' And if you agglomerate by Family, you don't see Genus and 
#' Species. But these are taxonomic ranks for ASVs that we got 
#' from QIIME2.  
#' * `Kingdom`  
#' * `Phylum`  
#' * `Class`  
#' * `Order`  
#' * `Family`  
#' * `Genus`  
#' * `Species`  
#' * `OTU`: ASV IDs from QIIME2. phyloseq uses OTU, so we keep it 
#' as it is. Not included if you are not aggomerating by ASVs.  
#' * `RelativeAbundance`: Relative abundance of taxa in each sample.
#' We calculate it by summing the Abundance of a taxon in each sample
#' and dividing that sum by the sum of reads in that sample.  
#' * `MeanRelativeAbundance`: Average relative abundance of a taxon 
#' in each host. We calculate it by summing the absolute 
#' abundance of a taxon from all samples in a host and 
#' dividing by the sum of reads in that host.  
#'
#' 2) custom.md is a dataframe with metadata. It has 5 variables:  
#' * `class`: same as in ps.q.agg  
#' * `animal`: full names of animal hosts.  The variable is 
#' factor with 9 levels at most ("Fukomys Damarensis", 
#' "FVB/N mouse", "Lepus europaeus", "MSM/Ms mouse", 
#' "naked mole rat", "Nannospalax leucodon", "Oryctolagus cuniculus", 
#' "Pteromys volans orii", "SPF mouse, B6").
#' "Fukomys Damarensis" is DMR, "FVB/N mouse" is FVBNmouse, "Lepus europaeus" is hare, 
#' "MSM/Ms mouse" is MSMmouse, "naked mole rat" is NMR, 
#' "Nannospalax leucodon" is spalax, "Oryctolagus cuniculus" is 
#' rabbit, "Pteromys volans orii" is pvo, "SPF mouse, B6" 
#' is B6mouse.  
#' * `sex`: sex of tested samples. Not all samples have it. It is a factor with
#' four levels (F, M, NR, -)  
#' * `birthday`: date of birth of samples. Not all samples have it. 
#' It is a Date format variable.  
#' * `Sample`: same as in ps.q.agg  
#'
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(microViz)

#+ echo=FALSE
## 2. Import data from QIIME2. #### 
#'
#' ## Import data from QIIME2.  
truncationlvl<-"234" # truncation level that we chose in QIIME2
authorname<-"pooled" # name of the folder with QIIME2 output
# qza_file_date_time<-"20240425_02_57_13"
qza_file_date_time<-"20260209_16_33_25"
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime"),
                    paste(qza_file_date_time,read.end.type,truncationlvl,sep="-")) # directory with QZA files

metadatadir<-file.path("./data/metadata",
                       paste(authorname,"metadata",sep = "-")) # directory with metadata

#' Specify the name of your metadata file.
metadata.filename<-file.path(metadatadir,
                          paste("filenames",read.end.type,
                                authorname,"raw-supercomp.tsv", 
                                sep = "-"))
biosample.md<-read.table("./data/metadata/pooled-metadata/biosample_metadata_for_ncbi.tsv",
                         sep = "\t", header= T)

#+ echo=FALSE
## 3. Import qza files and convert them into a phyloseq object. ####
#'
#' ## Import qza files and convert them into a phyloseq object.
ps.q<-qza_to_phyloseq(
  features = file.path(qiimedir, paste0(paste(qza_file_date_time,
                                              authorname,
                                              read.end.type,
                                              "trimmed-dada2-table",
                                              truncationlvl,
                                              "filtered",
                                              sep="-"),".qza")), # feature table
  taxonomy = file.path(qiimedir,paste0(paste(qza_file_date_time,
                                             authorname,
                                             read.end.type,
                                             "trimmed-dada2",
                                             truncationlvl,
                                             "filtered-taxonomy",
                                             sep="-"),".qza")), # taxonomy
  tree = file.path(qiimedir,paste0(paste(qza_file_date_time,
                                         authorname,
                                         read.end.type,
                                         "trimmed-dada2",
                                         truncationlvl,
                                         "filtered-rooted-tree",
                                         sep="-"),
                                   ".qza")) # rooted tree
)

#' Change the name d__Kingdom to Kingdom.
ps.q.taxtab<-as.data.frame(tax_table(ps.q))
ps.q.taxtab$Kingdom<-
  gsub("d__","",ps.q.taxtab$Kingdom)
tax_table(ps.q)<-as.matrix(ps.q.taxtab)
rm(ps.q.taxtab)

#+ echo=FALSE
## 4. Add custom metadata. ####
#'
#' ## Add custom metadata. 
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
#' Convert the Sample column into row names because phyloseq 
#' needs samples as rownames.
#' 
#' Remove absolute.filepath column.
custom.md<-custom.md%>%
  dplyr::select(-absolute.filepath)%>%
  mutate(class= as.factor(class),
         sex = as.factor(sex),
         birthday = as.Date(birthday),
         animal = as.factor(animal))
# birthday=as.Date(ifelse(class=="B6mouse",sampling_date-weeks(12),birthday)))
rownames(custom.md)<-custom.md$Sample

#' Assign the custom metadata as your phyloseq object's metadata.
sample_data(ps.q)<-custom.md

#+ echo=FALSE
### 4.1 For NMR, we create a separate metadata object with age groups. ####
#'
#' ### For NMR, we create a separate metadata object with age groups. ####
custom.md.ages<-custom.md%>%
  dplyr::select(Sample,class)%>%
  filter(class=="NMR")%>%
  left_join(biosample.md,by = join_by("Sample" =="host_subject_id"))%>%
  rename("sex"= host_sex,
         "age" = host_age)%>%
  dplyr::select(-organism,-env_broad_scale,
         -env_local_scale, -env_medium,
         -geo_loc_name, -lat_lon,
         -host)%>%
  mutate(age = as.numeric(age),
         agegroup=cut(age, breaks =c(0,10,16),
                      right = FALSE))
#' We create these new levels for differential microbial abundance.
unique_levels <-custom.md.ages %>%
  ungroup()%>%
  distinct(agegroup)%>%
  arrange(agegroup) %>%
  mutate(new_agegroup = paste0("agegroup", agegroup))%>%
  mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
  mutate(new_agegroup = gsub("\\,","_",new_agegroup))
custom.md.ages <- custom.md.ages %>%
  left_join(unique_levels, by = "agegroup")
#' We preserve the old group names for visualisation.
colnames(custom.md.ages)[which(colnames(custom.md.ages)=="agegroup")]<-"old_agegroup"
colnames(custom.md.ages)[which(colnames(custom.md.ages)=="new_agegroup")]<-"agegroup"

custom.md.ages<-custom.md.ages%>%
  as.data.frame()
rownames(custom.md.ages)<-custom.md.ages$sample_name

# saveRDS(custom.md.ages,file="./output/rdafiles/custom.md.ages.rds")
# write.table(custom.md.ages,file="./output/rtables/pooled/custom.md.ages.tsv",
#             row.names = F,sep = "\t")

#' You can exclude some samples based on class. Specify the excluded classes
#' in a vector, then use the `%in%` operator. It will remove entries 
#' of the `class` column (animal hosts) from the `custom.md` object (metadata).
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
#                                               'ntccontrol','rabbitcontrol',
#                                               'harecontrol'),]
custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
#' You can exclude samples based on their library size (total number of reads).
custom.md<-custom.md[!rownames(custom.md) %in%
                       intersect(names(which(colSums(ps.q@otu_table)<20000)),
                                 rownames(custom.md)),]
# saveRDS(custom.md,file="./output/rdafiles/custom.md.rds")
# write.table(custom.md,file="./output/rtables/pooled/custom.md.tsv",
#             row.names = F,sep = "\t")

#+ echo=FALSE
### 4.2 Construct the phyloseq object directly from QIIME2 output. ####
#'
#' ### Construct the phyloseq object directly from QIIME2 output. ####
#' We combine the phyloseq object with new metadata (if we excluded samples).
ps.foo <- phyloseq(otu_table(ps.q),
                   sample_data(custom.md),
                   tax_table(ps.q),
                   phy_tree(ps.q))
ps.q<-ps.foo
rm(ps.foo)

#' Number of features in the unfiltered dataset:
length(rownames(ps.q@tax_table@.Data))

#' Total frequency in the unfiltered dataset:
sum(colSums(ps.q@otu_table@.Data))

#' Summary statistics (min, median, max, quartiles) of the unfiltered dataset:
ps.q@otu_table@.Data%>%
  colSums()%>%
  summary()

#' Select only Bacteria. Remove chloroplast and mitochondria
ps.q<-ps.q %>%
  subset_taxa(Kingdom%in%"Bacteria")%>%
  subset_taxa(!Order %in% "Chloroplast")%>%
  subset_taxa(!Family %in% "Mitochondria")

#+ echo=FALSE
### 4.3 Fix empty taxa with higher rank taxon. ####
#'
#' ### Fix empty taxa with higher rank taxon. ####
#' Because we want to remove NA values and make ambiguous "uncultured" or 
#' "unclassified" taxa more understandable.
ps.q<-tax_fix(ps.q,unknowns = c("NA","uncultured","Unassigned",
                                "uncultured_bacterium","uncultured_rumen",
                                "gut_metagenome","human_gut","mouse_gut",
                                "wallaby_gut","uncultured_soil", 
                                "uncultured_organism","uncultured_prokaryote"))
#+ echo=FALSE
## 5. Convert the phyloseq object into a dataframe. ####
#'
#' ## Convert the phyloseq object into a dataframe.
ps.q.agg<-ps.q %>%
  psmelt() 
ps.q.agg.phylum<-ps.q %>%
  tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by agglom.rank
  psmelt()  # transform the phyloseq object into an R dataframe
ps.q.agg.family<-ps.q %>%
  tax_glom("Family",NArm = FALSE) %>% # agglomerate by agglom.rank
  psmelt()  # transform the phyloseq object into an R dataframe
ps.q.agg.genus<-ps.q %>%
  tax_glom("Genus",NArm = FALSE) %>% # agglomerate by agglom.rank
  psmelt()  # transform the phyloseq object into an R dataframe
ps.list <- list("OTU" = ps.q.agg, 
               "Phylum" = ps.q.agg.phylum,
               "Family" = ps.q.agg.family, 
               "Genus" = ps.q.agg.genus)
for (ps.df.index in names(ps.list)){
  ps.df <- ps.list[[ps.df.index]]
  print(paste("Parsing data from the", ps.df.index, "table"))
  print(paste("Number of rows in the dataframe:", nrow(ps.df)))
  
  # Remove entries with zero Abundance.
  print("Removing rows with zero Abundance")
  ps.df <- ps.df %>%
    filter(Abundance!=0)%>%
    dplyr::select(-sample_Sample)%>% # remove the duplicate column
    dplyr::select(-sex,-birthday,-animal)
  print(paste("Number of rows in the filtered dataset:",nrow(ps.df)))
  ### 5.1 Number of samples in the filtered dataset. ####
  print(paste("Number of samples in the filtered dataset:"))
  ps.df%>%
    distinct(Sample)%>%
    tally()%>%
    print()
  ### 5.2 Number of features in the filtered dataset. ####
  print(paste("Number of features in the filtered dataset:"))
  ps.df%>%
    distinct(OTU)%>%
    tally()%>%
    print()
  ### 5.3 Total frequency in the filtered dataset. ####
  print(paste("Total frequency in the filtered dataset:"))
  ps.df%>%
    summarise(TotalAbundance=sum(Abundance))%>%
    print()
  ### 5.4 Summary statistics (min, median, max, quartiles) of the filtered dataset. ####
  print(paste("Summary statistics (min, median, max, quartiles) of the filtered dataset:"))
  ps.df%>%
    dplyr::select(Sample,Abundance)%>%
    group_by(Sample)%>%
    summarise(FrequencyPerSample=sum(Abundance))%>%
    dplyr::select(FrequencyPerSample)%>%
    summary()%>%
    print()
  
  ## 6. Add relative abundance column: Abundance divided by total abundance in a sample. ####
  ps.df<-ps.df%>%
    group_by(class,Sample)%>%
    mutate(TotalSample=sum(Abundance))%>%
    group_by_at(c("class","Sample",ps.df.index))%>%
    mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
    ungroup()
  
  # Sanity check: is total relative abundance of each sample 100%?
  print(paste("Sanity check: is total relative abundance of each sample 100%?"))
  ps.df %>%
    group_by(Sample) %>% # Group by sample id
    summarise(sumRelativeAbundance = sum(RelativeAbundance)) %>% # Sum all abundances
    mutate(diff_from_100 = sumRelativeAbundance-100, # compare each value to 100
           is_different = as.logical(round(diff_from_100,digits = 10)))%>% 
    arrange(desc(diff_from_100))%>% # show the most deviating samples
    print()
  
  ps.df %>%
    group_by(Sample)%>%
    mutate(sumRelativeAbundance = sum(RelativeAbundance)) %>%
    ungroup()%>%
    distinct(Sample,.keep_all = T)%>%
    ggplot(aes(x=Sample,y=sumRelativeAbundance))+
    geom_bar(stat="identity")+
    coord_flip()
  print(last_plot())
  ### 6.1 Add mean relative abundance data. ####
  # We will group the dataset by three columns: class (animal host), 
  # two taxonomic ranks (e.g Genus, Family), and maybe OTU (actually ASV)
  # if we agglomerate at ASV level.
  
  # Group the dataframe by classes (animal hosts).
  # First, we calculate the library size per sample.
  # Then, inside each class, we take a agglom.rank, sum its abundances from all samples,
  # then take a mean. This will be our MeanRelativeAbundance.
  ps.df<-ps.df%>%
    group_by(class)%>% # group by class (animal host),
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",ps.df.index))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
    ungroup()
  
  
  if(ps.df.index != "OTU"){
    ps.df<-ps.df%>%
      dplyr::select(-OTU)
  }
  ps.df<-ps.df%>%
    ungroup()%>%
    dplyr::select(-TotalClass,-TotalSample,-TotalAgglomRank)
  # Save the tables in TSV format and as an RDS object
  # write.table(ps.df,
  #             file=file.path("./output/rtables",authorname,paste(
  #               paste(format(Sys.time(),format="%Y%m%d"),
  #                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #               "phyloseq-qiime",authorname,ps.df.index,read.end.type,truncationlvl,
  #               "table.tsv",sep="-")),
  #             row.names = F,sep = "\t")
  # saveRDS(ps.df,
  #         file=file.path("./output/rdafiles",paste(
  #           paste(format(Sys.time(),format="%Y%m%d"),
  #                 format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #           "phyloseq-qiime",authorname,ps.df.index,read.end.type,truncationlvl,
  #           "table.rds",sep="-")))
  
}

sessionInfo()
rm(list = ls(all=TRUE))
gc()