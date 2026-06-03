---
output: 
  bookdown::html_document2:
     toc: true
---





# Compare our data with results from the original papers.





## Introduction
In this script, we will compare the results from original papers on DMR,
flying squirrels (PVO), and hares and rabbits with our 
re-analysis. We will also create a bar plot of dominant families in our data
vs the original wild NMR data.




## Load necessary libraries and scripts.


``` r
# install.packages(c("tidyverse","ggtext"))
# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
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
library(phyloseq)
library(ggtext)
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
```



## Specifying parameters and directory/file names. 
Name of the folder with QIIME2 output:


``` r
authorname<-"pooled" 
```

Name of the folder with QIIME2 output:


``` r
rdafiles.directory<-"./output/rdafiles"
```

Truncation level that we chose in QIIME2:


``` r
truncationlvl<-"234" 
```

Single reads or paired reads: decided in QIIME2.


``` r
read.end.type<-"single"
```



## Import datasets.
Import datasets as rds files.


``` r
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(
    "20260211_17_01_07",
        "phyloseq-qiime",authorname,"OTU",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20260211_17_01_10",
        "phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
```

Import metadata:


``` r
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
```



## Re-analysis of the DMR dataset.
Data was taken from github of the original study.


``` r
dmr.asvs <- read_tsv("./data/dmr-original-data/asv_table_FreezedriedVsFrozen.tsv", 
                     col_types =cols(asv = col_character(),
                                     sample = col_character(),
                                     count = col_double(),
                                     relab = col_double() ))
dmr.metadata <- read_csv("./data/dmr-original-data/FDvsFrozenMetadata.csv", 
                         col_types = cols(.default = col_character(),
                                          SampleDate = col_date(format = ""),
                                          SampleOrder = col_number())) %>% 
  suppressWarnings()
```

Rename SampleNumber


``` r
dmr.metadata  <- dmr.metadata %>% 
  dplyr::rename(OldSampleNumber = SampleNumber, 
                SampleNumber = NewSampleNumber )
```

Make a character of the sample number


``` r
dmr.metadata <- dmr.metadata %>% 
  mutate(Asample = paste("A", sample, sep= "_"))
dmr.taxonomy <- read_tsv("./data/dmr-original-data/taxonomy_FreezedriedVsFrozen.tsv", 
                         col_types =  cols(asv = col_character(),
                                           kingdom = col_character(), 
                                           phylum = col_character(),
                                           class = col_character(),
                                           order = col_character(),
                                           family = col_character(),
                                           genus = col_character(),
                                           species = col_character()))
```

Rename unclassified taxa:


``` r
dmr.taxonomy  <- dmr.taxonomy %>%
  mutate(
    phylum = ifelse(is.na(phylum), sprintf("%s unclassified", kingdom), 
                    phylum),
    class = ifelse(is.na(class), sprintf("%s unclassified", 
                                         str_remove(phylum, "unclassified")), 
                   class),
    order = ifelse(is.na(order), sprintf("%s unclassified", 
                                         str_remove(class, "unclassified")), 
                   order),
    family = ifelse(is.na(family), sprintf("%s unclassified", 
                                           str_remove(order, "unclassified")), 
                    family),
    genus = ifelse(is.na(genus), sprintf("%s unclassified", 
                                         str_remove(family, "unclassified")), 
                   genus),
    species = ifelse(is.na(species), sprintf("%s unclassified", 
                                             str_remove(genus, "unclassified")), 
                     species))
```



### Filter to keep only frozen samples.
Also, join the asv data with taxonomy table.


``` r
dmr.asv.tax.df<-dmr.asvs%>%
  left_join(dmr.taxonomy)%>%
  rename("Sample"=sample)%>%
  rename("Abundance"=count)%>%
  rename("Kingdom"=kingdom)%>%
  rename("Phylum"=phylum)%>%
  rename("Class"=class)%>%
  rename("Order"=order)%>%
  rename("Family"=family)%>%
  rename("Genus"=genus)%>%
  rename("Species"=species)%>%
  mutate(class="DMR")
```

```
## Joining with `by = join_by(asv)`
```

``` r
dmr.asv.tax.df.frozen<-dmr.asv.tax.df%>%
  filter(Sample%in%subset(dmr.metadata,Treatment=="Frozen")$sample)
dmr.taxonomy.frozen<-dmr.taxonomy%>%
  filter(asv%in%dmr.asv.tax.df.frozen$asv)

## ====== Results after Table 1: ====== ####
```



### Compare numbers of phyla, families, genera, and asv between datasets.
Phyla in original frozen (19), original total (19), and my data (21).


``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Phylum")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR      19
```

``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Phylum")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR      19
```

``` r
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg.phylum, class=="DMR"),
                         tax.rank = "Phylum")
```

```
## # A tibble: 1 × 2
##   class     n
##   <fct> <int>
## 1 DMR      21
```

Families in original frozen (97), original total (117), and my data (156).


``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Family")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR      97
```

``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Family")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR     117
```

``` r
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg.family, class=="DMR"),
                         tax.rank = "Family")
```

```
## # A tibble: 1 × 2
##   class     n
##   <fct> <int>
## 1 DMR     156
```

Genera in original frozen (165), original total (210), and my data (257).


``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "Genus")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR     165
```

``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "Genus")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR     210
```

``` r
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg.genus, class=="DMR"),
                         tax.rank = "Genus")
```

```
## # A tibble: 1 × 2
##   class     n
##   <fct> <int>
## 1 DMR     257
```

ASVs in original frozen (1368), original total (1768), and my data (1826).


``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df.frozen,tax.rank = "asv")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR    1368
```

``` r
get_n_uniq_taxa_per_host(tax.df = dmr.asv.tax.df,tax.rank = "asv")
```

```
## # A tibble: 1 × 2
##   class     n
##   <chr> <int>
## 1 DMR    1768
```

``` r
get_n_uniq_taxa_per_host(tax.df = subset(ps.q.agg, class=="DMR"),
                         tax.rank = "OTU")
```

```
## # A tibble: 1 × 2
##   class     n
##   <fct> <int>
## 1 DMR    1826
```



### Classified and unclassified genera.
41 unclassified genera in frozen samples:


``` r
dmr.taxonomy.frozen[grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    41
```

51 unclassified genera in all samples:


``` r
dmr.taxonomy[grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    51
```

124 classified genera in frozen samples:


``` r
dmr.taxonomy.frozen[!grepl("uncultured|unclassified",dmr.taxonomy.frozen$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally 
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1   124
```

159 classified genera in all samples:


``` r
dmr.taxonomy[!grepl("uncultured|unclassified",dmr.taxonomy$genus),]%>%
  distinct(genus,.keep_all = T)%>%
  select(genus)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1   159
```



### How many total reads in frozen data.
1,853,295 reads in frozen data:


``` r
sum(dmr.asv.tax.df.frozen$Abundance)
```

```
## [1] 1853295
```

3,626,584 reads in total:


``` r
sum(dmr.asv.tax.df$Abundance)
```

```
## [1] 3626584
```



### Library size in frozen data.
Mean 92665, SD 44202


``` r
dmr.asv.tax.df.frozen%>%
  group_by(Sample)%>%
  summarise(LibrarySize=sum(Abundance))%>%
  summarise(MeanLibrarySize=mean(LibrarySize),
            SDLibrarySize=sd(LibrarySize))
```

```
## # A tibble: 1 × 2
##   MeanLibrarySize SDLibrarySize
##             <dbl>         <dbl>
## 1          92665.        44202.
```



### Dominant phyla in the frozen dataset vs my data (Results).
Original data:


``` r
dmr.asv.tax.df.frozen %>% 
  group_by(Sample, Phylum) %>% 
  summarise(sample_relab = sum(relab)*100) %>%
  group_by(Phylum) %>% 
  summarise(MeanRelativeAbundance =mean(sample_relab), 
            MinRelativeAbundance= min(sample_relab), 
            MaxRelativeAbundance= max(sample_relab))%>%
  arrange(desc(MeanRelativeAbundance)) %>% 
  mutate(MeanRelativeAbundance = round(MeanRelativeAbundance,digits = 5),
         MinRelativeAbundance = round(MinRelativeAbundance,digits = 5),
         MaxRelativeAbundance = round(MaxRelativeAbundance,digits = 5))
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 19 × 4
##    Phylum        MeanRelativeAbundance MinRelativeAbundance MaxRelativeAbundance
##    <chr>                         <dbl>                <dbl>                <dbl>
##  1 Bacteroidetes              68.9                 31.4                 86.7    
##  2 Firmicutes                 20.3                  5.20                47.8    
##  3 Proteobacter…               4.07                 0.319               60.9    
##  4 Cyanobacteria               3.37                 0.121                6.93   
##  5 Actinobacter…               2.46                 0.788                8.82   
##  6 Spirochaetes                0.406                0.0339               1.90   
##  7 Patescibacte…               0.272                0.00997              0.755  
##  8 Synergistetes               0.137                0.0177               0.260  
##  9 Lentisphaerae               0.102                0.00207              0.559  
## 10 Bacteria unc…               0.0621               0.0031               0.192  
## 11 Tenericutes                 0.0124               0.00374              0.0201 
## 12 Deferribacte…               0.00957              0.00171              0.0300 
## 13 Verrucomicro…               0.00816              0.00109              0.0310 
## 14 Euryarchaeota               0.0062               0.0062               0.0062 
## 15 Chloroflexi                 0.00378              0.00378              0.00378
## 16 Gemmatimonad…               0.00303              0.00303              0.00303
## 17 Chlamydiae                  0.00292              0.00292              0.00292
## 18 Fusobacteria                0.0029               0.0029               0.0029 
## 19 Planctomycet…               0.00179              0.00151              0.00207
```

My data:


``` r
ps.q.agg.phylum%>%
  filter(class=="DMR")%>% 
  group_by(Phylum) %>% 
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) %>% 
  mutate(MeanRelativeAbundance = round(MeanRelativeAbundance,digits = 5),
         MinRelativeAbundance = round(MinRelativeAbundance,digits = 5),
         MaxRelativeAbundance = round(MaxRelativeAbundance,digits = 5))%>%
  print(n = 21)
```

```
## # A tibble: 21 × 4
##    Phylum        MeanRelativeAbundance MinRelativeAbundance MaxRelativeAbundance
##    <chr>                         <dbl>                <dbl>                <dbl>
##  1 Bacteroidota               69.7                 41.8                 90.4    
##  2 Firmicutes                 22.5                  4.01                40.4    
##  3 Cyanobacteria               2.43                 0.0855               5.19   
##  4 Proteobacter…               2.22                 0.0299              40.0    
##  5 Actinobacter…               1.68                 0.734                5.32   
##  6 Rs-K70_termi…               0.374                0.0950               1.38   
##  7 Spirochaetota               0.347                0.0315               2.13   
##  8 Desulfobacte…               0.252                0.0386               1.24   
##  9 Patescibacte…               0.157                0.00726              0.478  
## 10 Bacteria Kin…               0.130                0.0230               0.394  
## 11 Synergistota                0.0947               0.0153               0.166  
## 12 Verrucomicro…               0.0649               0.00178              0.387  
## 13 Deferribacte…               0.00717              0.00138              0.0194 
## 14 Gemmatimonad…               0.00595              0.00595              0.00595
## 15 Chloroflexi                 0.00494              0.00093              0.0109 
## 16 Acidobacteri…               0.00285              0.00093              0.00793
## 17 Elusimicrobi…               0.00226              0.00226              0.00226
## 18 Fusobacterio…               0.00178              0.00178              0.00178
## 19 Methylomirab…               0.00149              0.00149              0.00149
## 20 Myxococcota                 0.00138              0.00062              0.00181
## 21 Planctomycet…               0.00126              0.00062              0.00198
```



## Re-analysis of the PVO dataset.
We don't have the original data, so I'm showing my results. The comparison
is provided in the manuscript.




### Dominant phyla (Results).


``` r
ps.q.agg.phylum%>%
  filter(class=="pvo")%>%
  ungroup()%>%
  group_by(Phylum)%>%
  summarise(MeanRelativeAbundance = round(mean(RelativeAbundance), digits = 5), 
            MinRelativeAbundance= round(min(RelativeAbundance), digits = 5), 
            MaxRelativeAbundance= round(max(RelativeAbundance), digits = 5),
            sd_relab = round(sd(RelativeAbundance), digits =5))%>%
  arrange(desc(MeanRelativeAbundance))
```

```
## # A tibble: 13 × 5
##    Phylum        MeanRelativeAbundance MinRelativeAbundance MaxRelativeAbundance
##    <chr>                         <dbl>                <dbl>                <dbl>
##  1 Firmicutes                 62.7                 27.8                 75.4    
##  2 Bacteroidota               20.6                 15.0                 31.4    
##  3 Proteobacter…               5.69                 0.357               44.9    
##  4 Bacteria Kin…               5.50                 1.69                 9.84   
##  5 Actinobacter…               3.44                 0.848                7.95   
##  6 Verrucomicro…               0.965                0.00436              5.59   
##  7 Cyanobacteria               0.870                0.341                1.73   
##  8 Elusimicrobi…               0.217                0.0155               0.749  
##  9 Desulfobacte…               0.207                0.0286               0.515  
## 10 Campylobacte…               0.0335               0.0152               0.0568 
## 11 Gemmatimonad…               0.0128               0.0128               0.0128 
## 12 Planctomycet…               0.0128               0.0128               0.0128 
## 13 Patescibacte…               0.00827              0.00782              0.00873
## # ℹ 1 more variable: sd_relab <dbl>
```



### Dominant families (Results).


``` r
get_dominant_taxa_in_host(ps.q.agg.family,
                          tax.rank = "Family",
                          host = "pvo",
                          nonbacterial.table = F)
```

```
## # A tibble: 96 × 3
## # Groups:   Family [96]
##    Phylum           Family                                MeanRelativeAbundance
##    <chr>            <chr>                                                 <dbl>
##  1 Firmicutes       Lachnospiraceae                                       36.6 
##  2 Bacteroidota     Muribaculaceae                                        15.3 
##  3 Firmicutes       Ruminococcaceae                                        9.22
##  4 Firmicutes       Clostridia_UCG-014                                     7.57
##  5 Firmicutes       Oscillospiraceae                                       7.12
##  6 Bacteria Kingdom Bacteria Kingdom                                       6.05
##  7 Actinobacteriota Eggerthellaceae                                        2.12
##  8 Bacteroidota     Rikenellaceae                                          1.60
##  9 Firmicutes       Clostridia_vadinBB60_group                             1.60
## 10 Firmicutes       [Eubacterium]_coprostanoligenes_group                  1.01
## # ℹ 86 more rows
```



## Re-analysis of the Hares and rabbits dataset.




### Remove pregnant and lactating samples.


``` r
selected.hare_rabbit.samples<-c("Hare 1","Hare 2", "Hare 3","Hare 4","Hare 5",
                                "Hare 6","Hare 7","Hare 8",
                                "Rabbit 2","Rabbit 3", "Rabbit 6","Rabbit 7",
                                "Rabbit 10", "Rabbit 11", "Rabbit 12")
```

Import abundances from the original paper (supplementary data). Only 
Illumina samples.


``` r
hare_rabbit.phylum<-read.table("./data/hares-rabbits-original-data/phyla.tsv",sep = "\t",
                               header = T)%>%
  as_tibble()
```

Change to long format and remove the "Unassigned" taxon.


``` r
hare_rabbit.phylum<-hare_rabbit.phylum%>%
  filter(Sample%in%selected.hare_rabbit.samples)%>%
  pivot_longer(!Sample,names_to = "Phylum",values_to = "Abundance")%>%
  filter(Phylum!="Unassigned")%>%
  filter(Abundance!=0)
```

Import family abundances.


``` r
hare_rabbit.family<-read.table("./data/hares-rabbits-original-data/families.tsv",
                               sep = "\t",
                               header = T)%>%
  as_tibble()%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Family",values_to = "Abundance")%>%
  filter(Family!="unassigned")%>%
  filter(Abundance!=0)
```

Import genus abundances. 


``` r
hare_rabbit.genus<-read.table("./data/hares-rabbits-original-data/genera.tsv",sep = "\t",
                              header = T)%>%
  as_tibble()%>%
  filter(Sample_ID%in%selected.hare_rabbit.samples)%>%
  rename(Sample=Sample_ID)%>%
  pivot_longer(!Sample,names_to = "Genus",values_to = "Abundance")%>%
  filter(Abundance!=0)
```



### Compare numbers of phyla, families, genera, and asv between datasets.
14 phyla in the original data


``` r
hare_rabbit.phylum%>%
  distinct(Phylum)%>%
  tally()
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    14
```

17 phyla in the my data


``` r
ps.q.agg.phylum%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Phylum)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    17
```

83 families in the original data


``` r
hare_rabbit.family%>%
  distinct(Family)%>%
  tally()
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    83
```

120 families in my data


``` r
ps.q.agg.family%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Family)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1   120
```

70 genera in the original data


``` r
hare_rabbit.genus%>%
  distinct(Genus)%>%
  tally()
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    70
```

241 genera in my data


``` r
ps.q.agg.genus%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(Genus)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1   241
```

4642 ASV in my data


``` r
ps.q.agg%>%
  filter(class%in%c("hare","rabbit"))%>%
  ungroup%>%
  distinct(OTU)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1  4642
```



### Compare dominant taxa between datasets.
Phyla in the original dataset.


``` r
hare_rabbit.phylum%>%
  mutate(class = ifelse(grepl("Hare",Sample), "hare", "rabbit"))%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(class,Phylum)%>%
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) 
```

```
## `summarise()` has grouped output by 'class'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 22 × 5
## # Groups:   class [2]
##    class  Phylum MeanRelativeAbundance MinRelativeAbundance MaxRelativeAbundance
##    <chr>  <chr>                  <dbl>                <dbl>                <dbl>
##  1 rabbit Firmi…                 82.7               64.9                   95.5 
##  2 hare   Firmi…                 62.3               49.0                   69.1 
##  3 hare   Bacte…                 31.0               22.4                   44.5 
##  4 rabbit Bacte…                  6.03               1.15                  18.6 
##  5 rabbit Verru…                  4.85               0.0343                21.3 
##  6 rabbit Cyano…                  2.33               0.0806                 8.11
##  7 rabbit Actin…                  1.83               0.528                  4.18
##  8 hare   Lenti…                  1.61               0.309                  3.45
##  9 hare   Prote…                  1.35               0.886                  2.29
## 10 rabbit Tener…                  1.32               0.259                  3.01
## # ℹ 12 more rows
```

Phyla in my dataset.


``` r
get_dominant_taxa_in_host(ps.q.agg.phylum,tax.rank = "Phylum",
                          host = c("hare","rabbit"),
                          nonbacterial.table = F)
```

```
## Adding missing grouping variables: `class`
```

```
## # A tibble: 27 × 3
## # Groups:   class, Phylum [27]
##    class  Phylum            MeanRelativeAbundance
##    <fct>  <chr>                             <dbl>
##  1 rabbit Firmicutes                       84.0  
##  2 hare   Firmicutes                       60.6  
##  3 hare   Bacteroidota                     32.9  
##  4 rabbit Bacteroidota                      7.15 
##  5 rabbit Verrucomicrobiota                 5.05 
##  6 rabbit Actinobacteriota                  2.03 
##  7 hare   Verrucomicrobiota                 1.36 
##  8 hare   Actinobacteriota                  1.18 
##  9 hare   Patescibacteria                   1.01 
## 10 hare   Synergistota                      0.967
## # ℹ 17 more rows
```

Families in the original dataset.


``` r
hare_rabbit.family%>%
  mutate(class = ifelse(grepl("Hare",Sample), "hare", "rabbit"))%>%
  group_by(Sample)%>%
  mutate(RelativeAbundance=Abundance/sum(Abundance)*100)%>%
  group_by(class,Family)%>%
  summarise(MeanRelativeAbundance =mean(RelativeAbundance ), 
            MinRelativeAbundance= min(RelativeAbundance ), 
            MaxRelativeAbundance= max(RelativeAbundance ))%>%
  arrange(desc(MeanRelativeAbundance)) %>%
  mutate(Family=sub(".*\\.\\.","",Family))%>%
  head
```

```
## `summarise()` has grouped output by 'class'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 6 × 5
## # Groups:   class [2]
##   class  Family  MeanRelativeAbundance MinRelativeAbundance MaxRelativeAbundance
##   <chr>  <chr>                   <dbl>                <dbl>                <dbl>
## 1 rabbit Rumino…                  42.8             5.64                     68.4
## 2 rabbit Clostr…                  37.2             0.0114                   87.5
## 3 hare   Rumino…                  33.8            28.3                      42.3
## 4 hare   Bacter…                  18.9            12.3                      26.8
## 5 hare   Lachno…                  17.2             8.61                     27.4
## 6 rabbit Peptos…                  16.0             0.000978                 32.0
```

Families in my dataset.


``` r
get_dominant_taxa_in_host(ps.q.agg.family,
                          tax.rank = "Family",
                          host = c("hare","rabbit"),
                          nonbacterial.table = F)
```

```
## Adding missing grouping variables: `class`
```

```
## # A tibble: 188 × 4
## # Groups:   class, Family [188]
##    class  Phylum       Family                MeanRelativeAbundance
##    <fct>  <chr>        <chr>                                 <dbl>
##  1 rabbit Firmicutes   Clostridiaceae                        22.9 
##  2 hare   Firmicutes   Lachnospiraceae                       17.3 
##  3 hare   Firmicutes   Ruminococcaceae                       17.3 
##  4 rabbit Firmicutes   Ruminococcaceae                       16.5 
##  5 rabbit Firmicutes   Oscillospiraceae                       9.93
##  6 hare   Bacteroidota Bacteroidaceae                         9.72
##  7 hare   Bacteroidota Bacteroidales Order                    9.52
##  8 hare   Firmicutes   Christensenellaceae                    8.54
##  9 hare   Bacteroidota Rikenellaceae                          5.99
## 10 rabbit Firmicutes   Peptostreptococcaceae                  5.44
## # ℹ 178 more rows
```



## Bar plot of dominant families reported in the wild NMR paper.


``` r
wild.nmr.families<-data.frame(
  Family=c("Lachnospiraceae",
           "Prevotellaceae",
           "Paraprevotellaceae",
           "Bacteroidales Order",
           "Clostridiales Order",
           "Oscillospiraceae",
           "Veillonellaceae",
           "Clostridiaceae",
           "Muribaculaceae",
           "Porphyromonadaceae",
           "Spirochaetaceae"),
  MeanRelativeAbundance=c(17.6,
                          11,
                          8.8,
                          6.2,
                          6.1,
                          5.7,
                          4.7,
                          4.1,
                          4,
                          3,
                          10.9)
)

ps.q.agg.dominant.families.nmr<-ps.q.agg.family%>%
  filter(class=="NMR")%>%
  distinct(class,Phylum,Family, MeanRelativeAbundance)%>%
  arrange(class,desc(MeanRelativeAbundance))
  
# Add a column if the family was reported as dominant or not
ps.q.agg.dominant.families.nmr.with_domin<-ps.q.agg.dominant.families.nmr%>%
  mutate(row.index=as.numeric(rownames(ps.q.agg.dominant.families.nmr)),
         dominant=ifelse(row.index<=10,"dominant","not_dominant"))

debebe.comparison.plot<-ps.q.agg.dominant.families.nmr.with_domin%>%
  select(Family,MeanRelativeAbundance,dominant)%>%
  filter(Family%in%wild.nmr.families$Family)%>%
  full_join(ps.q.agg.dominant.families.nmr.with_domin[1:10,
                                                      c("Family",
                                                        "MeanRelativeAbundance",
                                                        "dominant")])%>%
  full_join(wild.nmr.families[,c("Family","MeanRelativeAbundance")],
            by="Family")%>%
  rename(NMR=MeanRelativeAbundance.x,NMRwt=MeanRelativeAbundance.y)%>%
  pivot_longer(cols = c("NMR","NMRwt"),
               names_to = "class",
               values_to = "MeanRelativeAbundance")%>%
  arrange(-MeanRelativeAbundance)%>%
  mutate(Family = ifelse(!grepl(" ",Family),
                         paste0("<i>",Family,"</i>"),
                         Family),
         Family = gsub ("Bacteroidales Order", "<i>Bacteroidales Order</i>", Family),
         Family = gsub ("Clostridiales Order", "<i>Clostridiales Order</i>", Family)
         )%>%
  mutate(Family=factor(Family,levels=unique(Family)))%>%
  group_by(Family)%>%
  mutate(Family_id=cur_group_id(),
         Family_id=factor(Family_id,levels=unique(sort(Family_id))))%>%
  ungroup%>%
  mutate(animal=ifelse(class=="NMR","Lab-bred naked mole-rats","Wild naked mole-rats"))%>%
  replace_na(replace =list( MeanRelativeAbundance=0))%>%
  mutate(Family=gsub("_", " ", Family),
         Family = stringr::str_wrap(Family,width=30),
         Family = stringr::str_replace_all(Family, "\n", "<br>"),
         Family=paste(Family_id,Family,sep = ": "),
         Family=factor(Family, levels=unique(Family)))%>%
  arrange(animal,desc(MeanRelativeAbundance))%>%
  ggplot(aes(x=Family_id,
             y=MeanRelativeAbundance,
             fill=Family))+
  geom_bar(stat = "identity",
           position = position_dodge2(),
           width=0.8 # distance between bars
           )+
  theme_bw()+
  facet_grid(~animal,
             scales="free_x",  # each species will have its own bars inside
             # facet (instead of all bars)
             space = "free_x")+
  labs(x="",
       y="Average relative abundance (%)")+
  coord_cartesian(expand = c("bottom"=FALSE))+
  scale_fill_viridis_d(option = "C")+
  theme(axis.text.x = element_text(angle=0,size=11,hjust=0.5,colour = "black"),# shift 
        # the x-axis labels to the right
        # legend.position = "inside",
        axis.title = element_text(size = 15), # size of axis names
        axis.text.y = element_text(size=11,color="black"),
        strip.text.x = ggtext::element_markdown(size = 14),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        legend.position.inside = c(0.9,0.8),
        legend.text = ggtext::element_markdown(size=13),
        legend.title = element_text(size = 16), # size of legend title
        plot.caption = element_text(size=15), # size of plot caption
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
        )
```

```
## Joining with `by = join_by(Family, MeanRelativeAbundance, dominant)`
```

``` r
print(debebe.comparison.plot)
```

<img src="003-compare-data_files/figure-html/unnamed-chunk-65-1.png" alt="" width="1056" />

``` r
# for(image.format in c("png","tiff")){
#   ggsave(file.path("./images/barplots",
#                  paste0(paste(format(Sys.time(),format="%Y%m%d"),
#                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                        "-barplot-wild-vs-lab-nmr-from-paper.",image.format)),
#        plot=debebe.comparison.plot,
#        width=11, height=6,units="in",
#        # width = 5700,height =2800, units = "px",
#        dpi=300,device = image.format)
# }
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
##  [1] ggtext_0.1.2    phyloseq_1.50.0 lubridate_1.9.4 forcats_1.0.0  
##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5    
##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_4.0.0   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] ade4_1.7-23             tidyselect_1.2.1        viridisLite_0.4.3      
##  [4] farver_2.1.2            Biostrings_2.74.1       S7_0.2.0               
##  [7] fastmap_1.2.0           digest_0.6.37           timechange_0.3.0       
## [10] lifecycle_1.0.5         cluster_2.1.8.1         survival_3.8-3         
## [13] magrittr_2.0.3          compiler_4.4.3          rlang_1.1.5            
## [16] sass_0.4.10             tools_4.4.3             utf8_1.2.4             
## [19] igraph_2.1.4            yaml_2.3.12             data.table_1.17.8      
## [22] knitr_1.51              labeling_0.4.3          bit_4.6.0              
## [25] xml2_1.3.8              plyr_1.8.9              RColorBrewer_1.1-3     
## [28] withr_3.0.2             BiocGenerics_0.52.0     grid_4.4.3             
## [31] stats4_4.4.3            multtest_2.62.0         biomformat_1.34.0      
## [34] colorspace_2.1-2        Rhdf5lib_1.28.0         scales_1.4.0           
## [37] iterators_1.0.14        MASS_7.3-65             cli_3.6.4              
## [40] rmarkdown_2.30          vegan_2.6-4             crayon_1.5.3           
## [43] generics_0.1.4          otel_0.2.0              rstudioapi_0.18.0      
## [46] httr_1.4.7              reshape2_1.4.4          tzdb_0.5.0             
## [49] commonmark_2.0.0        ape_5.8-1               cachem_1.1.0           
## [52] rhdf5_2.50.2            zlibbioc_1.52.0         splines_4.4.3          
## [55] parallel_4.4.3          XVector_0.46.0          vctrs_0.6.5            
## [58] Matrix_1.7-4            jsonlite_2.0.0          litedown_0.9           
## [61] bookdown_0.46           IRanges_2.40.1          hms_1.1.4              
## [64] S4Vectors_0.44.0        bit64_4.6.0-1           foreach_1.5.2          
## [67] jquerylib_0.1.4         glue_1.8.0              codetools_0.2-20       
## [70] stringi_1.8.4           gtable_0.3.6            GenomeInfoDb_1.42.3    
## [73] UCSC.utils_1.2.0        pillar_1.11.1           htmltools_0.5.8.1      
## [76] rhdf5filters_1.18.1     GenomeInfoDbData_1.2.13 R6_2.6.1               
## [79] vroom_1.7.0             evaluate_1.0.5          lattice_0.22-6         
## [82] Biobase_2.66.0          markdown_2.0            gridtext_0.1.5         
## [85] bslib_0.10.0            Rcpp_1.0.14             nlme_3.1-167           
## [88] permute_0.9-8           mgcv_1.9-1              xfun_0.56              
## [91] pkgconfig_2.0.3
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 5621895 300.3    8825951 471.4  8701554 464.8
## Vcells 9733391  74.3   16146932 123.2 13389110 102.2
```

