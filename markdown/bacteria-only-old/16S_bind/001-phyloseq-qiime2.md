---
output: 
  bookdown::html_document2:
     toc: true
---




# Processing QIIME2 output into phyloseq format, agglomeration by taxonomic rank
 




## Introduction
Once you produce the feature table, taxonomic classification, 
and the phylogenetic tree in QIIME2, it's time to perform 
downstream processing in R. First, we need to import the 
QZA files using `qiime2R` package.
We will convert the QZA files directly into phyloseq objects.

The final output is the dataframe ps.q.agg and metadata custom.md.
ps.q.agg and custom.md are saved as tsv files and rds files.

1) ps.q.agg is an ASV table with 7-13 columns   
(7 if agglomerating at Phylum, 13 if agglomerating at ASV level):  
* `Sample`: samples that were sequenced.   
* `Abundance`: Absolute abundance of taxa.   
* `class`: short names of animal hosts. The variable is factor 
with 9 levels at most (B6mouse, DMR, FVBNmouse, hare, 
MSMmouse, NMR, pvo, rabbit, spalax).   
The next seven columns may not all be in the table. 
If you agglomerate by Genus, you don't see the Species column. 
And if you agglomerate by Family, you don't see Genus and 
Species. But these are taxonomic ranks for ASVs that we got 
from QIIME2.  
* `Kingdom`  
* `Phylum`  
* `Class`  
* `Order`  
* `Family`  
* `Genus`  
* `Species`  
* `OTU`: ASV IDs from QIIME2. phyloseq uses OTU, so we keep it 
as it is. Not included if you are not aggomerating by ASVs.  
* `RelativeAbundance`: Relative abundance of taxa in each sample.
We calculate it by summing the Abundance of a taxon in each sample
and dividing that sum by the sum of reads in that sample.  
* `MeanRelativeAbundance`: Average relative abundance of a taxon 
in each host. We calculate it by summing the absolute 
abundance of a taxon from all samples in a host and 
dividing by the sum of reads in that host.  

2) custom.md is a dataframe with metadata. It has 5 variables:  
* `class`: same as in ps.q.agg  
* `animal`: full names of animal hosts.  The variable is 
factor with 9 levels at most ("Fukomys Damarensis", 
"FVB/N mouse", "Lepus europaeus", "MSM/Ms mouse", 
"naked mole rat", "Nannospalax leucodon", "Oryctolagus cuniculus", 
"Pteromys volans orii", "SPF mouse, B6").
"Fukomys Damarensis" is DMR, "FVB/N mouse" is FVBNmouse, "Lepus europaeus" is hare, 
"MSM/Ms mouse" is MSMmouse, "naked mole rat" is NMR, 
"Nannospalax leucodon" is spalax, "Oryctolagus cuniculus" is 
rabbit, "Pteromys volans orii" is pvo, "SPF mouse, B6" 
is B6mouse.  
* `sex`: sex of tested samples. Not all samples have it. It is a factor with
four levels (F, M, NR, -)  
* `birthday`: date of birth of samples. Not all samples have it. 
It is a Date format variable.  
* `Sample`: same as in ps.q.agg  





## Load necessary libraries.


``` r
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
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   4.0.0     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.0.4     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ purrr::%||%()   masks base::%||%()
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(microViz)
```

```
## microViz version 0.12.1 - Copyright (C) 2021-2024 David Barnett
## ! Website: https://david-barnett.github.io/microViz
## ✔ Useful?  For citation details, run: `citation("microViz")`
## ✖ Silence? `suppressPackageStartupMessages(library(microViz))`
```



## Import data from QIIME2.  


``` r
truncationlvl<-"234" # truncation level that we chose in QIIME2
authorname<-"pooled" # name of the folder with QIIME2 output
# qza_file_date_time<-"20240425_02_57_13"
qza_file_date_time<-"20260209_16_33_25"
read.end.type<-"single" # single reads or paired reads: decided in QIIME2
qiimedir<-file.path("./output/qiime",paste0(authorname,"-qiime"),
                    paste(qza_file_date_time,read.end.type,truncationlvl,sep="-")) # directory with QZA files

metadatadir<-file.path("./data/metadata",
                       paste(authorname,"metadata",sep = "-")) # directory with metadata
```

Specify the name of your metadata file.


``` r
metadata.filename<-file.path(metadatadir,
                          paste("filenames",read.end.type,
                                authorname,"raw-supercomp.tsv", 
                                sep = "-"))
biosample.md<-read.table("./data/metadata/pooled-metadata/biosample_metadata_for_ncbi.tsv",
                         sep = "\t", header= T)
```



## Import qza files and convert them into a phyloseq object.


``` r
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
```

Change the name d__Kingdom to Kingdom.


``` r
ps.q.taxtab<-as.data.frame(tax_table(ps.q))
ps.q.taxtab$Kingdom<-
  gsub("d__","",ps.q.taxtab$Kingdom)
tax_table(ps.q)<-as.matrix(ps.q.taxtab)
rm(ps.q.taxtab)
```



## Add custom metadata. 


``` r
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
```

Convert the Sample column into row names because phyloseq 
needs samples as rownames.

Remove absolute.filepath column.


``` r
custom.md<-custom.md%>%
  dplyr::select(-absolute.filepath)%>%
  mutate(class= as.factor(class),
         sex = as.factor(sex),
         birthday = as.Date(birthday),
         animal = as.factor(animal))
# birthday=as.Date(ifelse(class=="B6mouse",sampling_date-weeks(12),birthday)))
rownames(custom.md)<-custom.md$Sample
```

Assign the custom metadata as your phyloseq object's metadata.


``` r
sample_data(ps.q)<-custom.md
```



### For NMR, we create a separate metadata object with age groups. ####


``` r
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
```

We create these new levels for differential microbial abundance.


``` r
unique_levels <-custom.md.ages %>%
  ungroup()%>%
  distinct(agegroup)%>%
  arrange(agegroup) %>%
  mutate(new_agegroup = paste0("agegroup", agegroup))%>%
  mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
  mutate(new_agegroup = gsub("\\,","_",new_agegroup))
custom.md.ages <- custom.md.ages %>%
  left_join(unique_levels, by = "agegroup")
```

We preserve the old group names for visualisation.


``` r
colnames(custom.md.ages)[which(colnames(custom.md.ages)=="agegroup")]<-"old_agegroup"
colnames(custom.md.ages)[which(colnames(custom.md.ages)=="new_agegroup")]<-"agegroup"

custom.md.ages<-custom.md.ages%>%
  as.data.frame()
rownames(custom.md.ages)<-custom.md.ages$sample_name

# saveRDS(custom.md.ages,file="./output/rdafiles/custom.md.ages.rds")
# write.table(custom.md.ages,file="./output/rtables/pooled/custom.md.ages.tsv",
#             row.names = F,sep = "\t")
```

You can exclude some samples based on class. Specify the excluded classes
in a vector, then use the `%in%` operator. It will remove entries 
of the `class` column (animal hosts) from the `custom.md` object (metadata).


``` r
# custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','pvo','tx',
#                                               'ntccontrol','rabbitcontrol',
#                                               'harecontrol'),]
custom.md<-custom.md[!custom.md$class  %in% c('pal','ppg','tx'),]
```

You can exclude samples based on their library size (total number of reads).


``` r
custom.md<-custom.md[!rownames(custom.md) %in%
                       intersect(names(which(colSums(ps.q@otu_table)<20000)),
                                 rownames(custom.md)),]
# saveRDS(custom.md,file="./output/rdafiles/custom.md.rds")
# write.table(custom.md,file="./output/rtables/pooled/custom.md.tsv",
#             row.names = F,sep = "\t")
```



### Construct the phyloseq object directly from QIIME2 output. ####
We combine the phyloseq object with new metadata (if we excluded samples).


``` r
ps.foo <- phyloseq(otu_table(ps.q),
                   sample_data(custom.md),
                   tax_table(ps.q),
                   phy_tree(ps.q))
ps.q<-ps.foo
rm(ps.foo)
```

Number of features in the unfiltered dataset:


``` r
length(rownames(ps.q@tax_table@.Data))
```

```
## [1] 13162
```

Total frequency in the unfiltered dataset:


``` r
sum(colSums(ps.q@otu_table@.Data))
```

```
## [1] 13027170
```

Summary statistics (min, median, max, quartiles) of the unfiltered dataset:


``` r
ps.q@otu_table@.Data%>%
  colSums()%>%
  summary()
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   21924   58453   88966  131588  154586  416954
```

Select only Bacteria. Remove chloroplast and mitochondria


``` r
ps.q<-ps.q %>%
  subset_taxa(Kingdom%in%"Bacteria")%>%
  subset_taxa(!Order %in% "Chloroplast")%>%
  subset_taxa(!Family %in% "Mitochondria")
```



### Fix empty taxa with higher rank taxon. ####
Because we want to remove NA values and make ambiguous "uncultured" or 
"unclassified" taxa more understandable.


``` r
ps.q<-tax_fix(ps.q,unknowns = c("NA","uncultured","Unassigned",
                                "uncultured_bacterium","uncultured_rumen",
                                "gut_metagenome","human_gut","mouse_gut",
                                "wallaby_gut","uncultured_soil", 
                                "uncultured_organism","uncultured_prokaryote"))
```



## Convert the phyloseq object into a dataframe.


``` r
ps.q.agg<-ps.q %>%
  psmelt() 
```

```
## Warning in psmelt(.): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

``` r
ps.q.agg.phylum<-ps.q %>%
  tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by phylum
  psmelt()  # transform the phyloseq object into an R dataframe
```

```
## Warning in psmelt(.): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

``` r
ps.q.agg.family<-ps.q %>%
  tax_glom("Family",NArm = FALSE) %>% # agglomerate by family
  psmelt()  # transform the phyloseq object into an R dataframe
```

```
## Warning in psmelt(.): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

``` r
ps.q.agg.genus<-ps.q %>%
  tax_glom("Genus",NArm = FALSE) %>% # agglomerate by genus
  psmelt()  # transform the phyloseq object into an R dataframe
```

```
## Warning in psmelt(.): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

``` r
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
```

```
## [1] "Parsing data from the OTU table"
## [1] "Number of rows in the dataframe: 1233639"
## [1] "Removing rows with zero Abundance"
## [1] "Number of rows in the filtered dataset: 41826"
## [1] "Number of samples in the filtered dataset:"
##    n
## 1 99
## [1] "Number of features in the filtered dataset:"
##       n
## 1 12461
## [1] "Total frequency in the filtered dataset:"
##   TotalAbundance
## 1       12911847
## [1] "Summary statistics (min, median, max, quartiles) of the filtered dataset:"
##  FrequencyPerSample
##  Min.   : 21302    
##  1st Qu.: 58453    
##  Median : 88959    
##  Mean   :130423    
##  3rd Qu.:154575    
##  Max.   :416947    
## [1] "Sanity check: is total relative abundance of each sample 100%?"
## # A tibble: 99 × 4
##    Sample sumRelativeAbundance diff_from_100 is_different
##    <chr>                 <dbl>         <dbl> <lgl>       
##  1 2D10                    100             0 FALSE       
##  2 2D14                    100             0 FALSE       
##  3 CG33                    100             0 FALSE       
##  4 D27                     100             0 FALSE       
##  5 DCG6                    100             0 FALSE       
##  6 F198                    100             0 FALSE       
##  7 F200                    100             0 FALSE       
##  8 F201                    100             0 FALSE       
##  9 F202                    100             0 FALSE       
## 10 F212                    100             0 FALSE       
## # ℹ 89 more rows
```

<img src="001-phyloseq-qiime2_files/figure-html/unnamed-chunk-30-1.png" alt="" width="672" />

```
## [1] "Parsing data from the Phylum table"
## [1] "Number of rows in the dataframe: 2475"
## [1] "Removing rows with zero Abundance"
## [1] "Number of rows in the filtered dataset: 1132"
## [1] "Number of samples in the filtered dataset:"
##    n
## 1 99
## [1] "Number of features in the filtered dataset:"
##    n
## 1 25
## [1] "Total frequency in the filtered dataset:"
##   TotalAbundance
## 1       12911847
## [1] "Summary statistics (min, median, max, quartiles) of the filtered dataset:"
##  FrequencyPerSample
##  Min.   : 21302    
##  1st Qu.: 58453    
##  Median : 88959    
##  Mean   :130423    
##  3rd Qu.:154575    
##  Max.   :416947    
## [1] "Sanity check: is total relative abundance of each sample 100%?"
## # A tibble: 99 × 4
##    Sample sumRelativeAbundance diff_from_100 is_different
##    <chr>                 <dbl>         <dbl> <lgl>       
##  1 F229                    100      1.42e-14 FALSE       
##  2 F514                    100      1.42e-14 FALSE       
##  3 F566                    100      1.42e-14 FALSE       
##  4 H15                     100      1.42e-14 FALSE       
##  5 MF_143                  100      1.42e-14 FALSE       
##  6 MF_144                  100      1.42e-14 FALSE       
##  7 PVO_17                  100      1.42e-14 FALSE       
##  8 R10                     100      1.42e-14 FALSE       
##  9 S15                     100      1.42e-14 FALSE       
## 10 2D10                    100      0        FALSE       
## # ℹ 89 more rows
```

<img src="001-phyloseq-qiime2_files/figure-html/unnamed-chunk-30-2.png" alt="" width="672" />

```
## [1] "Parsing data from the Family table"
## [1] "Number of rows in the dataframe: 21978"
## [1] "Removing rows with zero Abundance"
## [1] "Number of rows in the filtered dataset: 5010"
## [1] "Number of samples in the filtered dataset:"
##    n
## 1 99
## [1] "Number of features in the filtered dataset:"
##     n
## 1 222
## [1] "Total frequency in the filtered dataset:"
##   TotalAbundance
## 1       12911847
## [1] "Summary statistics (min, median, max, quartiles) of the filtered dataset:"
##  FrequencyPerSample
##  Min.   : 21302    
##  1st Qu.: 58453    
##  Median : 88959    
##  Mean   :130423    
##  3rd Qu.:154575    
##  Max.   :416947    
## [1] "Sanity check: is total relative abundance of each sample 100%?"
## # A tibble: 99 × 4
##    Sample sumRelativeAbundance diff_from_100 is_different
##    <chr>                 <dbl>         <dbl> <lgl>       
##  1 S11                     100      1.42e-14 FALSE       
##  2 2D10                    100      0        FALSE       
##  3 2D14                    100      0        FALSE       
##  4 CG33                    100      0        FALSE       
##  5 D27                     100      0        FALSE       
##  6 DCG6                    100      0        FALSE       
##  7 F198                    100      0        FALSE       
##  8 F200                    100      0        FALSE       
##  9 F201                    100      0        FALSE       
## 10 F202                    100      0        FALSE       
## # ℹ 89 more rows
```

<img src="001-phyloseq-qiime2_files/figure-html/unnamed-chunk-30-3.png" alt="" width="672" />

```
## [1] "Parsing data from the Genus table"
## [1] "Number of rows in the dataframe: 51876"
## [1] "Removing rows with zero Abundance"
## [1] "Number of rows in the filtered dataset: 9103"
## [1] "Number of samples in the filtered dataset:"
##    n
## 1 99
## [1] "Number of features in the filtered dataset:"
##     n
## 1 524
## [1] "Total frequency in the filtered dataset:"
##   TotalAbundance
## 1       12911847
## [1] "Summary statistics (min, median, max, quartiles) of the filtered dataset:"
##  FrequencyPerSample
##  Min.   : 21302    
##  1st Qu.: 58453    
##  Median : 88959    
##  Mean   :130423    
##  3rd Qu.:154575    
##  Max.   :416947    
## [1] "Sanity check: is total relative abundance of each sample 100%?"
## # A tibble: 99 × 4
##    Sample sumRelativeAbundance diff_from_100 is_different
##    <chr>                 <dbl>         <dbl> <lgl>       
##  1 2D10                    100             0 FALSE       
##  2 2D14                    100             0 FALSE       
##  3 CG33                    100             0 FALSE       
##  4 D27                     100             0 FALSE       
##  5 DCG6                    100             0 FALSE       
##  6 F198                    100             0 FALSE       
##  7 F200                    100             0 FALSE       
##  8 F201                    100             0 FALSE       
##  9 F202                    100             0 FALSE       
## 10 F212                    100             0 FALSE       
## # ℹ 89 more rows
```

<img src="001-phyloseq-qiime2_files/figure-html/unnamed-chunk-30-4.png" alt="" width="672" />

``` r
sessionInfo()
```

```
## R version 4.4.3 (2025-02-28 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 22631)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: Asia/Tokyo
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] microViz_0.12.1 lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
##  [5] dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1    
##  [9] tibble_3.2.1    ggplot2_4.0.0   tidyverse_2.0.0 phyloseq_1.50.0
## [13] qiime2R_0.99.6 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-9            gridExtra_2.3           permute_0.9-8          
##  [4] rlang_1.1.5             magrittr_2.0.3          ade4_1.7-23            
##  [7] otel_0.2.0              compiler_4.4.3          mgcv_1.9-1             
## [10] vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3        
## [13] crayon_1.5.3            fastmap_1.2.0           backports_1.5.0        
## [16] XVector_0.46.0          labeling_0.4.3          utf8_1.2.4             
## [19] rmarkdown_2.30          tzdb_0.5.0              UCSC.utils_1.2.0       
## [22] xfun_0.56               zlibbioc_1.52.0         cachem_1.1.0           
## [25] GenomeInfoDb_1.42.3     jsonlite_2.0.0          biomformat_1.34.0      
## [28] rhdf5filters_1.18.1     Rhdf5lib_1.28.0         parallel_4.4.3         
## [31] cluster_2.1.8.1         R6_2.6.1                bslib_0.10.0           
## [34] stringi_1.8.4           RColorBrewer_1.1-3      zCompositions_1.5.0-5  
## [37] rpart_4.1.24            jquerylib_0.1.4         Rcpp_1.0.14            
## [40] bookdown_0.46           iterators_1.0.14        knitr_1.51             
## [43] base64enc_0.1-3         IRanges_2.40.1          Matrix_1.7-4           
## [46] splines_4.4.3           nnet_7.3-20             igraph_2.1.4           
## [49] timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.18.0      
## [52] yaml_2.3.12             vegan_2.6-4             codetools_0.2-20       
## [55] lattice_0.22-6          plyr_1.8.9              Biobase_2.66.0         
## [58] withr_3.0.2             S7_0.2.0                evaluate_1.0.5         
## [61] foreign_0.8-90          survival_3.8-3          Biostrings_2.74.1      
## [64] pillar_1.11.1           checkmate_2.3.3         DT_0.34.0              
## [67] foreach_1.5.2           stats4_4.4.3            generics_0.1.4         
## [70] RCurl_1.98-1.17         truncnorm_1.0-9         S4Vectors_0.44.0       
## [73] hms_1.1.4               scales_1.4.0            glue_1.8.0             
## [76] Hmisc_5.2-3             tools_4.4.3             data.table_1.17.8      
## [79] rhdf5_2.50.2            grid_4.4.3              ape_5.8-1              
## [82] colorspace_2.1-2        nlme_3.1-167            GenomeInfoDbData_1.2.13
## [85] htmlTable_2.4.3         Formula_1.2-5           cli_3.6.4              
## [88] gtable_0.3.6            sass_0.4.10             digest_0.6.37          
## [91] BiocGenerics_0.52.0     htmlwidgets_1.6.4       farver_2.1.2           
## [94] htmltools_0.5.8.1       multtest_2.62.0         lifecycle_1.0.5        
## [97] httr_1.4.7              MASS_7.3-65
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##            used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells  5760375 307.7   10085427 538.7 10085427 538.7
## Vcells 10996046  83.9   67674822 516.4 84593527 645.4
```

