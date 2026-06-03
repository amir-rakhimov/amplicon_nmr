---
output: 
  bookdown::html_document2:
     toc: true
---





# Analysing QIIME2 data with phyloseq.





## Introduction
In this script, we will explore the imported dataset from QIIME2 using qiime2R
and phyloseq.

We will also rarefy the data for future analyses.

We will use the data from 001-phyloseq-qiime2.R script (ps.q.agg
agglomerated tables at phylum, family, genus, and OTU level).




## Load necessary libraries and scripts.


``` r
# install.packages(c("tidyverse","vegan","ggtext","Polychrome","ggrepel"))
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
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.6-4
```

``` r
library(Polychrome)
library(ggtext)
library(ggrepel)
```

Load necessary scripts.


``` r
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
source("./code/r-scripts/create_summary_stats_table.R")
source("./code/r-scripts/add_relab_to_tax_df.R")
source("./code/r-scripts/add_agegroup_to_tax_df.R")
source("./code/r-scripts/get_unclassified_summary_stats.R")
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/ggplot_species.R")
source("./code/r-scripts/add_zero_rows.R")
```



## Specifying parameters and directory/file names. 
Name of the folder with QIIME2 output:


``` r
authorname<-"pooled"
```

Directories with input files:


``` r
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)
```

Truncation level that we chose in QIIME2:


``` r
truncationlvl<-"234" 
```

Single reads or paired reads (decided in QIIME2):


``` r
read.end.type<-"single"
```

Import datasets as rds files.


``` r
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_38_18","phyloseq-qiime",authorname,"OTU",read.end.type,# old
  paste("20260211_17_01_07","phyloseq-qiime",authorname,"OTU",read.end.type, # new 
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240620_12_40_41","phyloseq-qiime",authorname,"Genus",read.end.type,
  paste("20260211_17_01_10","phyloseq-qiime",authorname,"Genus",read.end.type,
        truncationlvl,"table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_10_54_12-phyloseq-qiime",authorname,"Family",read.end.type,
  paste("20260211_17_01_09-phyloseq-qiime",authorname,"Family",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  # paste("20240917_21_29_36-phyloseq-qiime",authorname,"Phylum",read.end.type,
  paste("20260211_17_01_08-phyloseq-qiime",authorname,"Phylum",read.end.type,
        truncationlvl,"table.rds",sep = "-")))
```

Import metadata:


``` r
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
```

Specify paths and image formats:


``` r
barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/"
image.formats<-c("png","tiff")
```



## Setup plots.


``` r
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "DMR" = "*Fukomys damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*")
```

Use only the taxa that are present in the workspace
(custom.md is metadata from the rdafile).


``` r
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
```

Setup general ggplot theme


``` r
mytheme<-theme(axis.text.y = element_text(size=10), # size of y axis ticks
               axis.title = element_text(size = 10), # size of axis names
               legend.text = element_text(size = 10), # size of legend text
               legend.title = element_text(size = 15), # size of legend title
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank()
)
```



## Calculating summary statistics.





### Check the total number of unique ASV/phyla/families/genera per class.
We will use it for the summary table in the next section.


``` r
n.asv.per.host<-get_n_uniq_taxa_per_host(ps.q.agg,"OTU")
n.phylum.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.phylum,"Phylum")
n.family.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.family,"Family")
n.genus.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.genus,"Genus")
```



### Create a summary table.
The columns are:  
- Total reads  
- Library size (mean Abundance ± SD)  
- Number of ASV per host  
- Number of phyla per host  
- Number of families per host  
- Number of genera per host  


``` r
summary.stats.table<-create_summary_stats_table(ps.q.agg,
                                                n.asv.per.host,
                                                n.phylum.per.host,
                                                n.family.per.host,
                                                n.genus.per.host)
```

```
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
```

``` r
print(summary.stats.table)
```

```
## # A tibble: 9 × 9
## # Groups:   class [9]
##   class     TotalSamplesPerHost TotalReadsPerHost MeanLibrarySize SDLibrarySize
##   <fct>                   <int>             <dbl>           <dbl>         <dbl>
## 1 B6mouse                     4            350539           87635          9606
## 2 DMR                        20           3226720          161336         68602
## 3 FVBNmouse                   3            236930           78977          9116
## 4 hare                        8           2923880          365485         48701
## 5 MSMmouse                    8            613902           76738          8144
## 6 NMR                        24           1844935           76872         26287
## 7 pvo                        10            758468           75847         55176
## 8 rabbit                      7           2417378          345340         70861
## 9 spalax                     15            539095           35940         11383
## # ℹ 4 more variables: ASVPerHost <int>, PhylaPerHost <int>,
## #   FamiliesPerHost <int>, GeneraPerHost <int>
```

``` r
# write.table(summary.stats.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



## Add relative abundance and average relative abundance columns.


``` r
ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")
ps.q.agg.relab<-add_relab_to_tax_df(ps.q.agg,"OTU")
head(ps.q.agg.phylum.relab)
```

```
## # A tibble: 6 × 11
##   Sample Abundance class  Kingdom Phylum RelativeAbundance MeanRelativeAbundance
##   <chr>      <dbl> <fct>  <chr>   <chr>              <dbl>                 <dbl>
## 1 MF_153    391461 rabbit Bacter… Firmi…              94.4                  84.0
## 2 MF_154    363129 rabbit Bacter… Firmi…              89.6                  84.0
## 3 MF_144    305772 rabbit Bacter… Firmi…              87.3                  84.0
## 4 MF_151    287896 hare   Bacter… Firmi…              69.0                  60.6
## 5 MF_140    275649 rabbit Bacter… Firmi…              84.7                  84.0
## 6 MF_136    267945 hare   Bacter… Firmi…              64.5                  60.6
## # ℹ 4 more variables: TotalSample <dbl>, TotalClass <dbl>,
## #   TotalAgglomRank <dbl>, sdRelativeAbundance <dbl>
```

``` r
head(ps.q.agg.family.relab)
```

```
## # A tibble: 6 × 14
##   Sample Abundance class  Kingdom  Phylum   Class Order Family RelativeAbundance
##   <chr>      <dbl> <fct>  <chr>    <chr>    <chr> <chr> <chr>              <dbl>
## 1 MF_153    351733 rabbit Bacteria Firmicu… Clos… Clos… Clost…              84.8
## 2 MF_144    192256 rabbit Bacteria Firmicu… Clos… Clos… Clost…              54.9
## 3 MF_143    131605 rabbit Bacteria Firmicu… Clos… Pept… Pepto…              35.8
## 4 F568      122474 DMR    Bacteria Bactero… Bact… Bact… Murib…              38.1
## 5 MF_140    113241 rabbit Bacteria Firmicu… Clos… Osci… Rumin…              34.8
## 6 MF_136    112693 hare   Bacteria Firmicu… Clos… Lach… Lachn…              27.1
## # ℹ 5 more variables: MeanRelativeAbundance <dbl>, TotalSample <dbl>,
## #   TotalClass <dbl>, TotalAgglomRank <dbl>, sdRelativeAbundance <dbl>
```

``` r
head(ps.q.agg.genus.relab)
```

```
## # A tibble: 6 × 15
##   Sample Abundance class  Kingdom  Phylum       Class       Order   Family Genus
##   <chr>      <dbl> <fct>  <chr>    <chr>        <chr>       <chr>   <chr>  <chr>
## 1 MF_153    351692 rabbit Bacteria Firmicutes   Clostridia  Clostr… Clost… Clos…
## 2 MF_144    163176 rabbit Bacteria Firmicutes   Clostridia  Clostr… Clost… Clos…
## 3 MF_143    131605 rabbit Bacteria Firmicutes   Clostridia  Peptos… Pepto… Pept…
## 4 F568      120752 DMR    Bacteria Bacteroidota Bacteroidia Bacter… Murib… Muri…
## 5 MF_136    103379 hare   Bacteria Firmicutes   Clostridia  Lachno… Lachn… Lach…
## 6 F541       93198 DMR    Bacteria Bacteroidota Bacteroidia Bacter… Murib… Muri…
## # ℹ 6 more variables: RelativeAbundance <dbl>, MeanRelativeAbundance <dbl>,
## #   TotalSample <dbl>, TotalClass <dbl>, TotalAgglomRank <dbl>,
## #   sdRelativeAbundance <dbl>
```

``` r
head(ps.q.agg.relab)
```

```
## # A tibble: 6 × 17
##   OTU     Sample Abundance class Kingdom Phylum Class Order Family Genus Species
##   <chr>   <chr>      <dbl> <fct> <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>  
## 1 29e99f… MF_153    317376 rabb… Bacter… Firmi… Clos… Clos… Clost… Clos… Clostr…
## 2 c731b1… MF_144    144914 rabb… Bacter… Firmi… Clos… Clos… Clost… Clos… Clostr…
## 3 dc17e8… MF_143    109556 rabb… Bacter… Firmi… Clos… Pept… Pepto… Pept… Peptos…
## 4 a45df8… MF_143     87032 rabb… Bacter… Verru… Verr… Verr… Akker… Akke… Akkerm…
## 5 2972be… MF_140     67455 rabb… Bacter… Firmi… Clos… Osci… Rumin… Rumi… uncult…
## 6 4514fc… F519       65842 DMR   Bacter… Prote… Gamm… Ente… Enter… Esch… Escher…
## # ℹ 6 more variables: RelativeAbundance <dbl>, MeanRelativeAbundance <dbl>,
## #   TotalSample <dbl>, TotalClass <dbl>, TotalAgglomRank <dbl>,
## #   sdRelativeAbundance <dbl>
```



## Add agegroup variable to NMR data (must run for plotting).


``` r
custom.md.ages<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")

ps.q.agg.genus.relab.nmr<-ps.q.agg.genus.relab%>%
  filter(class=="NMR")

ps.q.agg.relab.nmr<-ps.q.agg.relab%>%
  filter(class=="NMR")
```

Add the age groups.


``` r
ps.q.agg.genus.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.genus.relab.nmr,"Genus",
                                             custom.md.ages)
ps.q.agg.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.relab.nmr,"OTU",
                                       custom.md.ages)
head(ps.q.agg.genus.relab.nmr)
```

```
## # A tibble: 6 × 18
##   Sample Abundance class Kingdom  Phylum       Class       Order    Family Genus
##   <chr>      <dbl> <fct> <chr>    <chr>        <chr>       <chr>    <chr>  <chr>
## 1 CG33       22665 NMR   Bacteria Firmicutes   Clostridia  Eubacte… Eubac… Euba…
## 2 2D10       17907 NMR   Bacteria Firmicutes   Bacilli     Erysipe… Erysi… Allo…
## 3 D27        16698 NMR   Bacteria Firmicutes   Clostridia  Eubacte… Eubac… Euba…
## 4 2D14       16108 NMR   Bacteria Firmicutes   Bacilli     Erysipe… Erysi… Allo…
## 5 L133       15999 NMR   Bacteria Bacteroidota Bacteroidia Bactero… Murib… Muri…
## 6 Y66b       15163 NMR   Bacteria Bacteroidota Bacteroidia Bactero… Murib… Muri…
## # ℹ 9 more variables: RelativeAbundance <dbl>, MeanRelativeAbundance <dbl>,
## #   TotalSample <dbl>, TotalClass <dbl>, TotalAgglomRank <dbl>,
## #   sdRelativeAbundance <dbl>, agegroup <chr>, old_agegroup <fct>,
## #   MeanRelativeAbundanceAgegroup <dbl>
```

``` r
head(ps.q.agg.relab.nmr)
```

```
## # A tibble: 6 × 20
##   OTU     Sample Abundance class Kingdom Phylum Class Order Family Genus Species
##   <chr>   <chr>      <dbl> <fct> <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>  
## 1 775f1d… D27        12318 NMR   Bacter… Firmi… Clos… Euba… Eubac… Euba… Eubact…
## 2 fbcb14… 2D10       12293 NMR   Bacter… Firmi… Baci… Erys… Erysi… Allo… Alloba…
## 3 5a8180… DCG6       11288 NMR   Bacter… Bacte… Bact… Bact… p-251… p-25… p-251-…
## 4 2a9b76… L122       10857 NMR   Bacter… Fibro… Fibr… Fibr… Fibro… Fibr… Fibrob…
## 5 fbcb14… 2D14       10674 NMR   Bacter… Firmi… Baci… Erys… Erysi… Allo… Alloba…
## 6 5a8180… GA17       10346 NMR   Bacter… Bacte… Bact… Bact… p-251… p-25… p-251-…
## # ℹ 9 more variables: RelativeAbundance <dbl>, MeanRelativeAbundance <dbl>,
## #   TotalSample <dbl>, TotalClass <dbl>, TotalAgglomRank <dbl>,
## #   sdRelativeAbundance <dbl>, agegroup <chr>, old_agegroup <fct>,
## #   MeanRelativeAbundanceAgegroup <dbl>
```



## Calculate summary stats of unclassified genera in each animal (unrarefied). ####
Here, we are interested in unclassified genera, but you can also try Families,
Orders, Classes, etc.


``` r
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
agglom.rank<-"Genus"
unclassified.genus.summary.stats.table<-
  get_unclassified_summary_stats(ps.q.agg.genus.relab,"Genus")
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
```

``` r
unclassified.genus.summary.stats.table
```

```
## # A tibble: 9 × 10
##   class     MeanTotalUnclassifie…¹ SDTotalUnclassifiedP…² minTotalUnclassified…³
##   <fct>                      <dbl>                  <dbl>                  <dbl>
## 1 pvo                           30                     11                     12
## 2 NMR                           29                      6                     17
## 3 hare                          28                      8                     18
## 4 DMR                           26                     11                      6
## 5 rabbit                        23                     15                      6
## 6 MSMmouse                      13                      9                      7
## 7 B6mouse                       11                      6                      4
## 8 spalax                        10                      7                      2
## 9 FVBNmouse                      9                      7                      4
## # ℹ abbreviated names: ¹​MeanTotalUnclassifiedPercent,
## #   ²​SDTotalUnclassifiedPercent, ³​minTotalUnclassifiedPercent
## # ℹ 6 more variables: maxTotalUnclassifiedPercent <dbl>,
## #   MedianTotalUnclassifiedPercent <dbl>, NumUnclassifiedTaxa <int>,
## #   NumClassifiedTaxa <int>, MinTotalTaxa <int>, MaxTotalTaxa <int>
```

``` r
# write.table(unclassified.genus.summary.stats.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



## Rarefy the table and check the percentage of unclassified taxa. ####
Convert the data frame into wide format: rows are samples and columns
are taxa


``` r
get_rarefied_table<-function(tax.df,tax.rank,host.classes){
  tax.df.wide<-tax.df%>%
    filter(class %in% host.classes,Abundance!=0)%>%
    dplyr::select(Sample,Abundance,class,all_of(tax.rank))%>%
    filter(Abundance!=0)%>%
    pivot_wider(names_from = all_of(tax.rank),
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("Sample")%>% # Set sample names as row names
    dplyr::select(-class)
  # Find the smallest sample size
  min.n_seqs.all<-tax.df%>%
    filter(class %in% host.classes,Abundance!=0)%>%
    dplyr::select(Sample,all_of(tax.rank),Abundance)%>%
    group_by(Sample)%>%
    summarize(n_seqs=sum(Abundance))%>%
    summarize(min=min(n_seqs))%>%
    pull(min)
  print(paste("Smallest sample size:", min.n_seqs.all))
  
  ### Rarefied asv table with vegan ####
  set.seed(1)
  tax.df.rare<-rrarefy(tax.df.wide,sample=min.n_seqs.all)
  tax.df.rare<-tax.df.rare%>%
    as_tibble(rownames="Sample")%>%
    pivot_longer(-Sample)%>%
    as.data.frame()%>%
    left_join(unique(tax.df[,c("Sample","class")]),
              by="Sample")
  if(tax.rank=="OTU"){
    tax.df.rare<-tax.df.rare%>%
      rename(OTU=name,Abundance=value)%>%
      filter(Abundance!=0)  
  }else{
    # rename the 'name' column corresponding to the tax.rank
    tax.df.rare[,paste(tax.rank)]<-tax.df.rare$name
    tax.df.rare<-tax.df.rare%>%
      dplyr::select(-name)%>%
      rename(Abundance=value)%>%
      filter(Abundance!=0)
  }
  # write.table(tax.df.rare,
  #             file = file.path(rtables.directory,paste0(
  #               paste(
  #                 paste(format(Sys.time(),format="%Y%m%d"),
  #                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #                 "ps.q.df.rare-nonfiltered",tax.rank,
  #                 paste(host.classes,collapse = '-'),sep = "-"),
  #               ".tsv")),
  #             row.names = F,
  #             sep = "\t")
  # saveRDS(tax.df.rare,
  #         file = file.path(rdafiles.directory,paste0(
  #           paste(
  #             paste(format(Sys.time(),format="%Y%m%d"),
  #                   format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #             "ps.q.df.rare-nonfiltered",tax.rank,
  #             paste(host.classes,collapse = '-'),sep = "-"),
  #           ".rds")))
  return(tax.df.rare)
}

ps.q.agg.rare<-get_rarefied_table(ps.q.agg,"OTU",custom.levels)
```

```
## [1] "Smallest sample size: 21302"
```

```
## Warning in rrarefy(tax.df.wide, sample = min.n_seqs.all): function should be
## used for observed counts, but smallest count is 2
```

``` r
ps.q.agg.phylum.rare<-get_rarefied_table(ps.q.agg.phylum,"Phylum",custom.levels)
```

```
## [1] "Smallest sample size: 21302"
```

```
## Warning in rrarefy(tax.df.wide, sample = min.n_seqs.all): function should be
## used for observed counts, but smallest count is 2
```

``` r
ps.q.agg.family.rare<-get_rarefied_table(ps.q.agg.family,"Family",custom.levels)
```

```
## [1] "Smallest sample size: 21302"
```

```
## Warning in rrarefy(tax.df.wide, sample = min.n_seqs.all): function should be
## used for observed counts, but smallest count is 2
```

``` r
ps.q.agg.genus.rare<-get_rarefied_table(ps.q.agg.genus,"Genus",custom.levels)
```

```
## [1] "Smallest sample size: 21302"
```

```
## Warning in rrarefy(tax.df.wide, sample = min.n_seqs.all): function should be
## used for observed counts, but smallest count is 2
```

``` r
ps.q.agg.genus.nmr.rare<-get_rarefied_table(ps.q.agg.genus.relab.nmr,"Genus","NMR")
```

```
## [1] "Smallest sample size: 25421"
```

```
## Warning in rrarefy(tax.df.wide, sample = min.n_seqs.all): function should be
## used for observed counts, but smallest count is 2
```

``` r
ps.q.agg.nmr.rare<-get_rarefied_table(ps.q.agg.relab.nmr,"OTU","NMR")
```

```
## [1] "Smallest sample size: 25421"
```

```
## Warning in rrarefy(tax.df.wide, sample = min.n_seqs.all): function should be
## used for observed counts, but smallest count is 2
```

``` r
head(ps.q.agg.relab.nmr)
```

```
## # A tibble: 6 × 20
##   OTU     Sample Abundance class Kingdom Phylum Class Order Family Genus Species
##   <chr>   <chr>      <dbl> <fct> <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>  
## 1 775f1d… D27        12318 NMR   Bacter… Firmi… Clos… Euba… Eubac… Euba… Eubact…
## 2 fbcb14… 2D10       12293 NMR   Bacter… Firmi… Baci… Erys… Erysi… Allo… Alloba…
## 3 5a8180… DCG6       11288 NMR   Bacter… Bacte… Bact… Bact… p-251… p-25… p-251-…
## 4 2a9b76… L122       10857 NMR   Bacter… Fibro… Fibr… Fibr… Fibro… Fibr… Fibrob…
## 5 fbcb14… 2D14       10674 NMR   Bacter… Firmi… Baci… Erys… Erysi… Allo… Alloba…
## 6 5a8180… GA17       10346 NMR   Bacter… Bacte… Bact… Bact… p-251… p-25… p-251-…
## # ℹ 9 more variables: RelativeAbundance <dbl>, MeanRelativeAbundance <dbl>,
## #   TotalSample <dbl>, TotalClass <dbl>, TotalAgglomRank <dbl>,
## #   sdRelativeAbundance <dbl>, agegroup <chr>, old_agegroup <fct>,
## #   MeanRelativeAbundanceAgegroup <dbl>
```

``` r
head(ps.q.agg.genus.nmr.rare)
```

```
##   Sample Abundance class                      Genus
## 1   CG33      5236   NMR      Eubacteriaceae Family
## 2   CG33      2828   NMR                Allobaculum
## 3   CG33      2468   NMR             Muribaculaceae
## 4   CG33      1349   NMR                Fibrobacter
## 5   CG33       357   NMR     Lachnospiraceae Family
## 6   CG33      1906   NMR Erysipelotrichaceae Family
```

``` r
head(ps.q.agg.nmr.rare)
```

```
##   Sample                              OTU Abundance class
## 1    D27 775f1d61e616978239ee36c337b0403c      3284   NMR
## 2    D27 5a818007a6452d5db9223ef1388e451c      1580   NMR
## 3    D27 2a9b76ddc8fc1e5f0e10f1fd40de1c79       784   NMR
## 4    D27 c513b4ef037af47125244161fb1eab50      1356   NMR
## 5    D27 72aa784b6a3f05f45910fe33902caa81       353   NMR
## 6    D27 37d6420f6843a99bbb892a3b37178db1      1000   NMR
```



### Add relative abundances and taxonomic information to the rarefied dataframe. 
All hosts (genus)


``` r
ps.q.agg.genus.rare.relab<-add_relab_to_tax_df(ps.q.agg.genus.rare,"Genus")
head(ps.q.agg.genus.rare.relab)
```

```
## # A tibble: 6 × 10
##   Sample Abundance class  Genus         TotalSample RelativeAbundance TotalClass
##   <chr>      <dbl> <fct>  <chr>               <dbl>             <dbl>      <dbl>
## 1 MF_153     18071 rabbit Clostridium_…       21302           84.8        149114
## 2 MF_153        42 rabbit Lachnospirac…       21302            0.197      149114
## 3 MF_153       314 rabbit Ruminococcus        21302            1.47       149114
## 4 MF_153        11 rabbit Akkermansia         21302            0.0516     149114
## 5 MF_153        10 rabbit Escherichia-…       21302            0.0469     149114
## 6 MF_153        16 rabbit Bacteroidale…       21302            0.0751     149114
## # ℹ 3 more variables: TotalAgglomRank <dbl>, MeanRelativeAbundance <dbl>,
## #   sdRelativeAbundance <dbl>
```

Add other taxonomic ranks to the dataframe


``` r
ps.q.agg.genus.rare.relab<-ps.q.agg.genus.rare.relab%>%
  left_join(unique(ps.q.agg.genus[,c("Kingdom","Phylum","Class","Order","Family","Genus")]))
```

```
## Joining with `by = join_by(Genus)`
```

``` r
head(ps.q.agg.genus.rare.relab)
```

```
## # A tibble: 6 × 15
##   Sample Abundance class  Genus         TotalSample RelativeAbundance TotalClass
##   <chr>      <dbl> <fct>  <chr>               <dbl>             <dbl>      <dbl>
## 1 MF_153     18071 rabbit Clostridium_…       21302           84.8        149114
## 2 MF_153        42 rabbit Lachnospirac…       21302            0.197      149114
## 3 MF_153       314 rabbit Ruminococcus        21302            1.47       149114
## 4 MF_153        11 rabbit Akkermansia         21302            0.0516     149114
## 5 MF_153        10 rabbit Escherichia-…       21302            0.0469     149114
## 6 MF_153        16 rabbit Bacteroidale…       21302            0.0751     149114
## # ℹ 8 more variables: TotalAgglomRank <dbl>, MeanRelativeAbundance <dbl>,
## #   sdRelativeAbundance <dbl>, Kingdom <chr>, Phylum <chr>, Class <chr>,
## #   Order <chr>, Family <chr>
```

Add relative abundances to NMR dataframe (ASV level)


``` r
ps.q.agg.nmr.rare.relab<-add_relab_to_tax_df(ps.q.agg.nmr.rare,"OTU")
ps.q.agg.nmr.rare.relab<-ps.q.agg.nmr.rare.relab%>%
  left_join(unique(ps.q.agg[,c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")]))
```

```
## Joining with `by = join_by(OTU)`
```

``` r
### 8.2 Plot a rarefaction curve. ####
# ps.q.mat<-as(t(otu_table(ps.q)),"matrix") # from phyloseq
# ps.q.genus.mat<-ps.q.agg.genus%>%
#   filter(class %in% custom.levels,Abundance!=0)%>%
#   dplyr::select(Sample,Abundance,class,all_of(agglom.rank))%>%
#   filter(Abundance!=0)%>%
#   pivot_wider(names_from = all_of(agglom.rank),
#               values_from = "Abundance",
#               values_fill = 0)%>%
#   as.data.frame()%>%
#   column_to_rownames("Sample")%>% # Set sample names as row names
#   dplyr::select(-class)%>%
#   as.matrix() # convert to matrix
# set.seed(1)
# rare.df<-rarecurve(ps.q.genus.mat,step = 100,sample=min(rowSums(ps.q.genus.mat)),tidy = TRUE)
# rare.df%>%
#   # filter(Sample<=100000)%>%
#   group_by(Site)%>% # site is sample name
#   mutate(label=if_else(Sample==max(Sample),as.character(Site),NA_character_))%>%
#   # filter(Site%in%rownames(custom.md[which(custom.md$class=="NMR"),]))%>%
#   filter(Site%in%unique(ps.q.agg.genus$Sample[ps.q.agg.genus$class%in%custom.levels]))%>%
#   ggplot(.,aes(x=Sample,y=Species,col=Site))+
#   geom_line()+
#   # coord_cartesian(xlim=c(0,100000))+
#   geom_vline(xintercept = min(rowSums(ps.q.genus.mat)))+
#   annotate("text",
#            x=min(rowSums(ps.q.genus.mat))+2000,
#            y=10,
#            label=min(rowSums(ps.q.genus.mat)))+
#   geom_label_repel(aes(label = label),
#                    nudge_x = 1,
#                    na.rm = TRUE) +
#   theme_bw()+
#   labs(x="Sample size",
#        y="ASV")+
#   theme(legend.position = "none")
# ggsave(paste0("./images/lineplots/",
#               paste(paste(format(Sys.time(),format="%Y%m%d"),
#                           format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                     "rarecurve",
#                     truncationlvl,agglom.rank,
#                     sep = "-"),".png"),
#        plot=last_plot(),
#        width = 4500,height = 3000,
#        units = "px",dpi=300,device = "png")
```



### Create a summary stats table for the rarefied dataframe.


``` r
n.asv.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.rare,"OTU")
n.phylum.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.phylum.rare,"Phylum")
n.family.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.family.rare,"Family")
n.genus.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.genus.rare,"Genus")

summary.stats.table.rare<-create_summary_stats_table(ps.q.agg.rare,
                                                n.asv.per.host.rare,
                                                n.phylum.per.host.rare,
                                                n.family.per.host.rare,
                                                n.genus.per.host.rare)
```

```
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
```

``` r
summary.stats.table.rare
```

```
## # A tibble: 9 × 9
## # Groups:   class [9]
##   class     TotalSamplesPerHost TotalReadsPerHost MeanLibrarySize SDLibrarySize
##   <fct>                   <int>             <dbl>           <dbl>         <dbl>
## 1 B6mouse                     4             85208           21302             0
## 2 DMR                        20            426040           21302             0
## 3 FVBNmouse                   3             63906           21302             0
## 4 hare                        8            170416           21302             0
## 5 MSMmouse                    8            170416           21302             0
## 6 NMR                        24            511248           21302             0
## 7 pvo                        10            213020           21302             0
## 8 rabbit                      7            149114           21302             0
## 9 spalax                     15            319530           21302             0
## # ℹ 4 more variables: ASVPerHost <int>, PhylaPerHost <int>,
## #   FamiliesPerHost <int>, GeneraPerHost <int>
```



## Calculate summary stats of unclassified genera in rarefied data.


``` r
unclassified.genus.summary.stats.table.rare<-get_unclassified_summary_stats(ps.q.agg.genus.rare.relab,
                                                                            "Genus")
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
## Joining with `by = join_by(class)`
```

``` r
unclassified.genus.summary.stats.table.rare
```

```
## # A tibble: 9 × 10
##   class     MeanTotalUnclassifie…¹ SDTotalUnclassifiedP…² minTotalUnclassified…³
##   <fct>                      <dbl>                  <dbl>                  <dbl>
## 1 pvo                           30                     11                     12
## 2 NMR                           29                      6                     17
## 3 hare                          28                      8                     18
## 4 DMR                           26                     11                      5
## 5 rabbit                        23                     15                      6
## 6 MSMmouse                      13                      9                      7
## 7 B6mouse                       11                      6                      4
## 8 spalax                        10                      7                      2
## 9 FVBNmouse                      9                      7                      4
## # ℹ abbreviated names: ¹​MeanTotalUnclassifiedPercent,
## #   ²​SDTotalUnclassifiedPercent, ³​minTotalUnclassifiedPercent
## # ℹ 6 more variables: maxTotalUnclassifiedPercent <dbl>,
## #   MedianTotalUnclassifiedPercent <dbl>, NumUnclassifiedTaxa <int>,
## #   NumClassifiedTaxa <int>, MinTotalTaxa <int>, MaxTotalTaxa <int>
```

``` r
# write.table(unclassified.genus.summary.stats.table.rare,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "unclassified-genus-summary-table-rarefied.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



## Check the most abundant phyla, families, genera in NMR and other hosts.
Phyla:


``` r
ps.q.agg.dominant.phyla.all_hosts<-ps.q.agg.phylum.relab%>%
  group_by(class,Phylum)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class,Phylum, MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class,Phylum)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.phyla.all_hosts)
```

```
## # A tibble: 6 × 7
##   class   Phylum  MeanRelativeAbundance sdRelativeAbundance     min    max     n
##   <fct>   <chr>                   <dbl>               <dbl>   <dbl>  <dbl> <int>
## 1 B6mouse Firmic…                51.6                 6.51  43.3    59.0       4
## 2 B6mouse Bacter…                45.5                 7.53  36.2    54.3       4
## 3 B6mouse Actino…                 1.02                0.652  0.341   1.89      4
## 4 B6mouse Proteo…                 0.722               0.291  0.427   1.06      4
## 5 B6mouse Desulf…                 0.535               0.812  0.0649  1.70      4
## 6 B6mouse Patesc…                 0.266               0.130  0.0955  0.389     4
```

``` r
ps.q.agg.dominant.phyla.nmr<-ps.q.agg.dominant.phyla.all_hosts%>%
  filter(class=="NMR")
head(ps.q.agg.dominant.phyla.nmr)
```

```
## # A tibble: 6 × 7
##   class Phylum     MeanRelativeAbundance sdRelativeAbundance     min   max     n
##   <fct> <chr>                      <dbl>               <dbl>   <dbl> <dbl> <int>
## 1 NMR   Bacteroid…                 47.0               10.4   23.8    61.0     24
## 2 NMR   Firmicutes                 38.7               10.5   20.6    60       24
## 3 NMR   Fibrobact…                  5.71               5.01   0.0185 18.6     23
## 4 NMR   Spirochae…                  2.63               2.65   0.633  11.6     23
## 5 NMR   Synergist…                  1.16               0.648  0.168   3.31    24
## 6 NMR   Actinobac…                  1.10               0.569  0.0661  2.64    24
```

Families


``` r
ps.q.agg.dominant.families.all_hosts<-ps.q.agg.family.relab%>%
  group_by(class,Family)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class, Phylum, Family, MeanRelativeAbundance, sdRelativeAbundance, 
           min,max,n)%>%
  group_by(class,Phylum,Family)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.families.all_hosts)
```

```
## # A tibble: 6 × 8
##   class   Phylum   Family MeanRelativeAbundance sdRelativeAbundance    min   max
##   <fct>   <chr>    <chr>                  <dbl>               <dbl>  <dbl> <dbl>
## 1 B6mouse Bactero… Murib…                 35.0                 6.09 29.6   43.4 
## 2 B6mouse Firmicu… Lachn…                 18.3                12.7   5.52  35.4 
## 3 B6mouse Firmicu… Lacto…                 14.6                13.2   1.50  30.4 
## 4 B6mouse Firmicu… Erysi…                  8.13                7.79  0.138 17.1 
## 5 B6mouse Firmicu… Clost…                  3.83                7.88  0.219 14.3 
## 6 B6mouse Bactero… Prevo…                  3.75                3.21  1.85   8.69
## # ℹ 1 more variable: n <int>
```

``` r
ps.q.agg.dominant.families.nmr<-ps.q.agg.dominant.families.all_hosts%>%
  filter(class=="NMR")
head(ps.q.agg.dominant.families.nmr) 
```

```
## # A tibble: 6 × 8
##   class Phylum    Family MeanRelativeAbundance sdRelativeAbundance     min   max
##   <fct> <chr>     <chr>                  <dbl>               <dbl>   <dbl> <dbl>
## 1 NMR   Bacteroi… Prevo…                 19.5                 7.61 0.0242   31.7
## 2 NMR   Firmicut… Erysi…                 11.5                 6.87 0.834    24.4
## 3 NMR   Bacteroi… Murib…                 11.1                 3.14 5.97     17.9
## 4 NMR   Firmicut… Eubac…                  8.09                5.68 0.519    20.2
## 5 NMR   Firmicut… Lachn…                  6.49                5.20 1.67     23.5
## 6 NMR   Bacteroi… p-251…                  6.32                3.73 0.00441  12.6
## # ℹ 1 more variable: n <int>
```

Genera:


``` r
ps.q.agg.dominant.genera.all_hosts<-ps.q.agg.genus.relab%>%
  group_by(class, Genus)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class,Phylum,Family,Genus, MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class,Phylum,Family,Genus)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.genera.all_hosts)
```

```
## # A tibble: 6 × 9
##   class   Phylum   Family Genus MeanRelativeAbundance sdRelativeAbundance    min
##   <fct>   <chr>    <chr>  <chr>                 <dbl>               <dbl>  <dbl>
## 1 B6mouse Bactero… Murib… Muri…                 34.0                 6.32 28.0  
## 2 B6mouse Firmicu… Lachn… Lach…                  9.26                4.91  3.00 
## 3 B6mouse Firmicu… Lacto… Ligi…                  6.74                7.51  0.215
## 4 B6mouse Firmicu… Lachn… Lach…                  5.57                6.34  0.731
## 5 B6mouse Firmicu… Lacto… Lact…                  5.20                6.45  1.09 
## 6 B6mouse Bactero… Bacte… Bact…                  3.65                3.60  1.22 
## # ℹ 2 more variables: max <dbl>, n <int>
```

``` r
ps.q.agg.dominant.genera.nmr<-ps.q.agg.dominant.genera.all_hosts%>%
  filter(class=="NMR")
head(ps.q.agg.dominant.genera.nmr)
```

```
## # A tibble: 6 × 9
##   class Phylum    Family Genus MeanRelativeAbundance sdRelativeAbundance     min
##   <fct> <chr>     <chr>  <chr>                 <dbl>               <dbl>   <dbl>
## 1 NMR   Bacteroi… Murib… Muri…                 11.1                 3.14 5.97   
## 2 NMR   Firmicut… Eubac… Euba…                  7.95                5.73 0.438  
## 3 NMR   Bacteroi… p-251… p-25…                  6.32                3.73 0.00441
## 4 NMR   Fibrobac… Fibro… Fibr…                  5.71                5.01 0.0185 
## 5 NMR   Firmicut… Erysi… Erys…                  5.53                4.57 0.783  
## 6 NMR   Bacteroi… Prevo… Prev…                  5.48                2.66 0.00551
## # ℹ 2 more variables: max <dbl>, n <int>
```

``` r
# write.table(ps.q.agg.dominant.phyla.nmr,
#             file=file.path(rtables.directory,
#                            # paste("20240523_12_02_33", 
#                             paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-phyla-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
# 
# write.table(ps.q.agg.dominant.families.nmr,
#             file=file.path(rtables.directory,
#                            # paste("20240523_12_07_02",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-families-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
# 
# write.table(ps.q.agg.dominant.genera.nmr,
#             file=file.path(rtables.directory,
#                            # paste("20240523_12_31_19",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-genera-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
# # For all hosts
# write.table(ps.q.agg.dominant.phyla.all_hosts,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-phyla-table-all_hosts.tsv",sep="-")),
#             row.names = F,sep = "\t")
# 
# write.table(ps.q.agg.dominant.families.all_hosts,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-families-table-all_hosts.tsv",sep="-")),
#             row.names = F,sep = "\t")
# 
# write.table(ps.q.agg.dominant.genera.all_hosts,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "dominant-genera-table-all_hosts.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



## Check how much Bacteroidaceae are in NMR.


``` r
bacteroidaceae.nmr<-ps.q.agg.dominant.families.nmr%>%
  filter(Family=="Bacteroidaceae")
bacteroidaceae.nmr
```

```
## # A tibble: 1 × 8
##   class Phylum     Family MeanRelativeAbundance sdRelativeAbundance    min   max
##   <fct> <chr>      <chr>                  <dbl>               <dbl>  <dbl> <dbl>
## 1 NMR   Bacteroid… Bacte…                  1.11                1.69 0.0923  8.68
## # ℹ 1 more variable: n <int>
```

``` r
# write.table(bacteroidaceae.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "bacteroidaceae-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



### Check the most dominant Bacteroidota families in NMR.


``` r
bacteroidota.nmr<-ps.q.agg.dominant.families.nmr%>%
  filter(Phylum=="Bacteroidota")
head(bacteroidota.nmr)
```

```
## # A tibble: 6 × 8
##   class Phylum    Family MeanRelativeAbundance sdRelativeAbundance     min   max
##   <fct> <chr>     <chr>                  <dbl>               <dbl>   <dbl> <dbl>
## 1 NMR   Bacteroi… Prevo…                 19.5                 7.61 0.0242  31.7 
## 2 NMR   Bacteroi… Murib…                 11.1                 3.14 5.97    17.9 
## 3 NMR   Bacteroi… p-251…                  6.32                3.73 0.00441 12.6 
## 4 NMR   Bacteroi… Riken…                  3.59                2.21 1.31    11.1 
## 5 NMR   Bacteroi… Palud…                  2.17                1.51 0.616    6.82
## 6 NMR   Bacteroi… Bacte…                  1.96                1.25 0.517    4.89
## # ℹ 1 more variable: n <int>
```

``` r
# write.table(bacteroidota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "bacteroidota-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



## Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. 


``` r
spirochaetaceae.nmr<-ps.q.agg.dominant.families.nmr%>%
  filter(Family=="Spirochaetaceae",class=="NMR")
spirochaetaceae.nmr
```

```
## # A tibble: 1 × 8
##   class Phylum      Family MeanRelativeAbundance sdRelativeAbundance   min   max
##   <fct> <chr>       <chr>                  <dbl>               <dbl> <dbl> <dbl>
## 1 NMR   Spirochaet… Spiro…                  2.63                2.65 0.633  11.6
## # ℹ 1 more variable: n <int>
```

``` r
spirochaetota.nmr<-ps.q.agg.dominant.phyla.nmr%>%
  filter(Phylum=="Spirochaetota")
spirochaetota.nmr
```

```
## # A tibble: 1 × 7
##   class Phylum       MeanRelativeAbundance sdRelativeAbundance   min   max     n
##   <fct> <chr>                        <dbl>               <dbl> <dbl> <dbl> <int>
## 1 NMR   Spirochaeto…                  2.63                2.65 0.633  11.6    23
```

``` r
# write.table(spirochaetaceae.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "spirochaetaceae-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")

treponema.nmr<-ps.q.agg.dominant.genera.nmr%>%
  filter(Genus=="Treponema")
treponema.nmr
```

```
## # A tibble: 1 × 9
##   class Phylum      Family Genus MeanRelativeAbundance sdRelativeAbundance   min
##   <fct> <chr>       <chr>  <chr>                 <dbl>               <dbl> <dbl>
## 1 NMR   Spirochaet… Spiro… Trep…                  2.20                2.71 0.159
## # ℹ 2 more variables: max <dbl>, n <int>
```

``` r
# write.table(treponema.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "treponema-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



### Check the number of ASVs in Treponema from NMR.


``` r
ps.q.agg%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1   109
```



## Check Mogibacteriaceae (renamed to Anaerovoracaceae) in all hosts.


``` r
mogibacteriaceae_anaerovoracaceae.all<-ps.q.agg.dominant.families.all_hosts%>%
  filter(Family=="Anaerovoracaceae")
head(mogibacteriaceae_anaerovoracaceae.all)
```

```
## # A tibble: 6 × 8
##   class     Phylum Family MeanRelativeAbundance sdRelativeAbundance    min   max
##   <fct>     <chr>  <chr>                  <dbl>               <dbl>  <dbl> <dbl>
## 1 B6mouse   Firmi… Anaer…                 0.136              0.0392 0.0872 0.183
## 2 DMR       Firmi… Anaer…                 0.194              0.190  0.0279 0.829
## 3 FVBNmouse Firmi… Anaer…                 0.371              0.103  0.286  0.486
## 4 hare      Firmi… Anaer…                 0.499              0.146  0.339  0.721
## 5 MSMmouse  Firmi… Anaer…                 0.361              0.117  0.212  0.528
## 6 NMR       Firmi… Anaer…                 0.740              0.510  0.0744 1.71 
## # ℹ 1 more variable: n <int>
```

``` r
# write.table(mogibacteriaceae_anaerovoracaceae.all,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



## Analyse sulfur-metabolising bacteria in NMR. 


``` r
desulfobacterota.nmr<-ps.q.agg.dominant.genera.nmr%>%
  filter(Phylum=="Desulfobacterota",class=="NMR")
head(desulfobacterota.nmr)
```

```
## # A tibble: 5 × 9
##   class Phylum    Family Genus MeanRelativeAbundance sdRelativeAbundance     min
##   <fct> <chr>     <chr>  <chr>                 <dbl>               <dbl>   <dbl>
## 1 NMR   Desulfob… Desul… Desu…               0.385               0.307   0.0107 
## 2 NMR   Desulfob… Desul… Desu…               0.152               0.179   0.00661
## 3 NMR   Desulfob… Desul… Bilo…               0.0398              0.0905  0.00305
## 4 NMR   Desulfob… Desul… Desu…               0.00835             0.0176  0.00367
## 5 NMR   Desulfob… Desul… Mail…               0.00488             0.00805 0.00392
## # ℹ 2 more variables: max <dbl>, n <int>
```

Desulfobacterota in other hosts


``` r
desulfobacterota.all<-ps.q.agg.dominant.genera.all_hosts%>%
  filter(Phylum=="Desulfobacterota")
head(desulfobacterota.all)
```

```
## # A tibble: 6 × 9
##   class   Phylum  Family Genus MeanRelativeAbundance sdRelativeAbundance     min
##   <fct>   <chr>   <chr>  <chr>                 <dbl>               <dbl>   <dbl>
## 1 B6mouse Desulf… Desul… Desu…                0.424              1.03    0.0218 
## 2 B6mouse Desulf… Desul… Desu…                0.0628             0.0580  0.0746 
## 3 B6mouse Desulf… Desul… Bilo…                0.0485             0.00396 0.0606 
## 4 DMR     Desulf… Desul… Desu…                0.168              0.243   0.0137 
## 5 DMR     Desulf… Desul… Desu…                0.0493             0.0452  0.00528
## 6 DMR     Desulf… Desul… Bilo…                0.0294             0.0273  0.00356
## # ℹ 2 more variables: max <dbl>, n <int>
```



### Total Desulfobacterota in NMR.


``` r
ps.q.agg.dominant.phyla.nmr%>%
  filter(Phylum=="Desulfobacterota")
```

```
## # A tibble: 1 × 7
##   class Phylum      MeanRelativeAbundance sdRelativeAbundance    min   max     n
##   <fct> <chr>                       <dbl>               <dbl>  <dbl> <dbl> <int>
## 1 NMR   Desulfobac…                 0.590               0.315 0.0827  1.38    24
```

``` r
# write.table(desulfobacterota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "desulfobacterota-table-nmr.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(desulfobacterota.all,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "desulfobacterota-table-all.tsv",sep="-")),
#             row.names = F,sep = "\t")
```



### Plot Desulfobacterota.


``` r
# It's a flipped plot
desulfobacterota.plot<-ps.q.agg.genus.relab%>%
  mutate(class=factor(class,levels=custom.levels))%>%
  filter(Phylum=="Desulfobacterota")%>%
  group_by(class,Genus)%>%
  ggplot(aes(x=factor(class,levels=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE,
               outliers = FALSE)+
  geom_jitter(aes(color=class),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=2,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
  facet_wrap(~Genus,scales = "free_x",
             ncol = 2)+
  theme_bw()+
  coord_flip()+
  labs(x="",
       y="Relative abundance (%)")+
  # VVV scale_x_discrete works on x axis but the final labels are on y because the coord_flip() flipped the plot.
  # So, we use rev() to address the flipping VVV
  scale_x_discrete(labels=rev(pretty.level.names), # new labels (named vector) on the axis
                   limits=rev(custom.levels) # limits adjust which levels (and in what order) are displayed
  )+ # rename boxplot labels 
  scale_fill_viridis_d(option="C")+
  stat_summary(fun=median, geom="point", shape=23, size=2, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  # ggtitle(paste0("Relative abundance of Desulfobacterota phylum members"))+
  mytheme+
  theme(axis.text.y = ggtext::element_markdown(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(size=10),
        strip.text.x = element_text(size=10),
        plot.title = element_text(size = 8), # size of plot title
        plot.caption = element_text(size=8), # size of plot caption
        legend.position = "none")
```

``` r
print(desulfobacterota.plot + 
  ggtitle(paste0("Relative abundance of Desulfobacterota phylum members"))+
  theme(plot.title = element_text(size = 14))
)
```

<img src="002-summary-stats-qiime2_files/figure-html/unnamed-chunk-61-1.png" alt="" width="768" />

``` r
# for (image.format in c("png","tiff")){
#   ggsave(paste0("./images/taxaboxplots/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "Desulfobacterota phylum members-all-hosts",
#                       sep = "-"),".",image.format),
#          plot=desulfobacterota.plot,
#          width=8, height=10,units="in",
#          # width = 4000,height = 6000,
#          # width = 1200,
#          # units = "px",
#          dpi=300,device = image.format)
# }
```



## Plot Treponema and other Spirochaetota.


``` r
spirochaetota.plot<-ps.q.agg.genus.relab%>%
  mutate(class=factor(class,levels=custom.levels))%>%
  filter(Phylum=="Spirochaetota")%>%
  group_by(class,Genus)%>%
  ggplot(aes(x=factor(class,levels=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE,
               outliers = FALSE)+
  geom_jitter(aes(color=class),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=2,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
  facet_wrap(~Genus,scales = "free_x",
             ncol = 2)+
  theme_bw()+
  coord_flip()+
  labs(x="",
       y="Relative abundance (%)")+
  # VVV scale_x_discrete works on x axis but the final labels are on y because the coord_flip() flipped the plot.
  # So, we use rev() to address the flipping VVV
  scale_x_discrete(labels=rev(pretty.level.names), # new labels (named vector) on the axis
                   limits=rev(custom.levels) # limits adjust which levels (and in what order) are displayed
  )+ # rename boxplot labels 
  scale_fill_viridis_d(option="C")+
  stat_summary(fun=median, geom="point", shape=23, size=2, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  # ggtitle(paste0("Relative abundance of Spirochaetota phylum members"))+
  mytheme+
  theme(axis.text.y = ggtext::element_markdown(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(size=10),
        strip.text.x = element_text(size=10),
        plot.title = element_text(size = 8), # size of plot title
        plot.caption = element_text(size=8), # size of plot caption
        legend.position = "none")
```

``` r
print(spirochaetota.plot + 
        ggtitle(paste0("Relative abundance of Spirochaetota phylum members"))+
        theme(plot.title = element_text(size = 14))
        )
```

<img src="002-summary-stats-qiime2_files/figure-html/unnamed-chunk-64-1.png" alt="" width="768" />

``` r
# for (image.format in c("png","tiff")){
#   ggsave(paste0("./images/taxaboxplots/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "Spirochaetota phylum members-all-hosts",
#                       sep = "-"),".",image.format),
#          plot=spirochaetota.plot,
#          width=8, height=10,units="in",
#          # width = 4000,height = 6000,
#          # width = 1200,
#          # units = "px",
#          dpi=300,device = image.format)
# }
```



## Analysis of naked mole-rat data ASVs.
Give ASVs shorter names: Genus, "ASV", OTU, OTU number. For example, Allobaculum_ASV_22.


``` r
nmr.asv.names<-ps.q.agg.relab.nmr%>%
  dplyr::select(OTU,Genus)%>%
  group_by(Genus)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(Genus,OTU)%>%
  mutate(row.index=row_number())%>%
  mutate(ASV_name=paste(Genus,row.index,sep="_ASV_"))
head(nmr.asv.names)
```

```
## # A tibble: 6 × 4
## # Groups:   Genus [3]
##   OTU                              Genus                      row.index ASV_name
##   <chr>                            <chr>                          <int> <chr>   
## 1 1630c97bc99db59912c92545f184a2c8 ASF356                             1 ASF356_…
## 2 a5d3cef1b84bffa024ae0259ebbe7819 ASF356                             2 ASF356_…
## 3 fc224cd4a1ffe36b74bdec848bbba352 ASF356                             3 ASF356_…
## 4 3c23ed3c668cb5db5fe54e0920da4c77 Absconditabacteriales_(SR…         1 Abscond…
## 5 446609f70ba4143bd817f44c210a0a4c Acetanaerobacterium                1 Acetana…
## 6 ba4f8c14f1a6f625c90f05848475b43a Acetanaerobacterium                2 Acetana…
```

The new name becomes the OTU column. The old name becomes OTU_old_name column.


``` r
ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")])%>%
  rename("OTU_old_name"="OTU",
         "OTU"="ASV_name")%>%
  relocate(OTU,.before = Sample)%>%
  relocate(OTU_old_name,.after = Genus)
```

```
## Joining with `by = join_by(OTU, Genus)`
```


### How many ASVs are shared between two age groups?
First, find ASVs in young samples


``` r
otu.young<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  dplyr::select(OTU,agegroup,MeanRelativeAbundanceAgegroup, sdRelativeAbundance)
```

Next, find ASVs in old samples


``` r
otu.old<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  dplyr::select(OTU,agegroup,MeanRelativeAbundanceAgegroup)
```

1745 ASVs in young individuals:


``` r
nrow(otu.young) 
```

```
## [1] 1745
```

771 ASVs in old individuals:


``` r
nrow(otu.old) 
```

```
## [1] 771
```

668 shared ASVs:


``` r
shared.otu<-intersect(otu.young$OTU,otu.old$OTU)
length(shared.otu)
```

```
## [1] 668
```

ASVs unique to old samples: 771 - 668 = 103


``` r
length(otu.old$OTU)-length(shared.otu)
```

```
## [1] 103
```



### Check 103 ASV only in old individuals.


``` r
ps.q.agg.relab.nmr%>%
  group_by(OTU,agegroup)%>%
  mutate(n_samples=n_distinct(Sample))%>% # Find the number of old samples the ASV was found in
  filter(OTU %in% otu.old$OTU,
         !OTU %in% shared.otu)%>%
  filter(n_samples >= 3)%>% 
  arrange(desc(MeanRelativeAbundanceAgegroup))%>%
  dplyr::select(OTU,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)
```

```
## Adding missing grouping variables: `agegroup`
```

```
## # A tibble: 0 × 4
## # Groups:   OTU, agegroup [0]
## # ℹ 4 variables: agegroup <chr>, OTU <chr>, MeanRelativeAbundance <dbl>,
## #   MeanRelativeAbundanceAgegroup <dbl>
```

No ASVs with at least 3 samples! The 103 ASVs are individual-specific.
Which genera do the 103 old-specific ASVs belong to?


``` r
ps.q.agg.relab.nmr%>%
  filter(OTU %in% otu.old$OTU, # ASV in old but not shared vector
         !OTU %in% shared.otu)%>%
  dplyr::select(OTU,Family,Genus, MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
  group_by(Family, Genus)%>%
  summarise(n=n_distinct(OTU))%>%
  arrange(desc(n))%>%
  head
```

```
## `summarise()` has grouped output by 'Family'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 6 × 3
## # Groups:   Family [6]
##   Family              Genus                      n
##   <chr>               <chr>                  <int>
## 1 Lachnospiraceae     Lachnospiraceae Family    14
## 2 Bacteria Kingdom    Bacteria Kingdom          11
## 3 Prevotellaceae      Prevotellaceae Family      7
## 4 Spirochaetaceae     Treponema                  5
## 5 Bacteroidales Order Bacteroidales Order        4
## 6 Muribaculaceae      Muribaculaceae             4
```



### How much % do shared ASVs take on average?


``` r
ps.q.agg.relab.nmr%>%
  # no separation by age
  filter(OTU%in%shared.otu)%>%
  group_by(Sample)%>%
  summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
  summarise(MeanRelAbSharedASVTotal = mean(SumRelAbSharedASV),
            sdRelaAbSharedASVTotal = sd(SumRelAbSharedASV))
```

```
## # A tibble: 1 × 2
##   MeanRelAbSharedASVTotal sdRelaAbSharedASVTotal
##                     <dbl>                  <dbl>
## 1                    94.2                   10.9
```

On average 94.2% ± 10.9%





### How much % do shared ASVs take in each age group?


``` r
ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  # separation by age
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbSharedASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbSharedASVTotalAge=mean(SumRelAbSharedASV),
            sdRelaAbSharedASVTotal = sd(SumRelAbSharedASV))
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 2 × 3
##   agegroup      MeanRelAbSharedASVTotalAge sdRelaAbSharedASVTotal
##   <chr>                              <dbl>                  <dbl>
## 1 agegroup0_10                        92.6                  12.2 
## 2 agegroup10_16                       98.8                   1.77
```

Higher variation in young individuals





### Are the shared ASVs enriched in certain genera? 
First, get the table of shared ASVs and their genera.


``` r
shared.otu.genera<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  # keep unique rows
  distinct(Genus,OTU)%>%
  group_by(Genus)
head(shared.otu.genera)
```

```
## # A tibble: 6 × 2
## # Groups:   Genus [5]
##   Genus                      OTU                              
##   <chr>                      <chr>                            
## 1 Eubacteriaceae Family      Eubacteriaceae Family_ASV_16     
## 2 Allobaculum                Allobaculum_ASV_22               
## 3 p-251-o5                   p-251-o5_ASV_2                   
## 4 Fibrobacter                Fibrobacter_ASV_2                
## 5 Erysipelotrichaceae Family Erysipelotrichaceae Family_ASV_11
## 6 Fibrobacter                Fibrobacter_ASV_5
```

Calculate the ASVs in each genus with cumsum() and pull the most numerous genera.


``` r
shared.otu.genera.cumsum<-shared.otu.genera%>%
  group_by(Genus)%>%
  # count rows (ASVs) for each genus; add as a column for sorting
  summarise(num_asvs =n_distinct(OTU))%>%
  # genera with the highest number of ASVs will be on top
  arrange(-num_asvs,Genus)%>%
  ungroup()%>%
  # cumulative sum shows how many ASVs the top genera take
  mutate(cum_sum=cumsum(num_asvs))

head(shared.otu.genera.cumsum)
```

```
## # A tibble: 6 × 3
##   Genus                   num_asvs cum_sum
##   <chr>                      <int>   <int>
## 1 Lachnospiraceae Family        55      55
## 2 Muribaculaceae                54     109
## 3 Treponema                     42     151
## 4 Bacteria Kingdom              30     181
## 5 Oscillospiraceae Family       22     203
## 6 Parabacteroides               22     225
```

Six genera account for 1/3 shared ASVs (225 out of 668 ASVS)




### How many ASVs of the top 5 genera in shared.otu.genera are actually found there?
The five most common genera according to the shared.otu.genera dataframe are 
Lachnospiraceae Family (55 ASVs), Muribaculaceae (54 ASVS), 
Treponema (42 ASVs), Bacteria Kingdom (30 ASVs), and
Oscillospiraceae Family (22 ASVs).  


``` r
ps.q.agg.relab.nmr%>%
  filter(Genus%in%c("Lachnospiraceae Family","Treponema",
                    "Muribaculaceae","Bacteria Kingdom",
                    "Oscillospiraceae Family"))%>%
  distinct(Genus,OTU)%>%
  mutate(is_shared=ifelse(OTU %in%shared.otu.genera$OTU, TRUE, FALSE))%>%
  group_by(Genus,is_shared)%>%
  summarise(n_asvs=n())%>%
  ungroup()%>%
  arrange(Genus,desc(is_shared))
```

```
## `summarise()` has grouped output by 'Genus'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 10 × 3
##    Genus                   is_shared n_asvs
##    <chr>                   <lgl>      <int>
##  1 Bacteria Kingdom        TRUE          30
##  2 Bacteria Kingdom        FALSE        125
##  3 Lachnospiraceae Family  TRUE          55
##  4 Lachnospiraceae Family  FALSE        108
##  5 Muribaculaceae          TRUE          54
##  6 Muribaculaceae          FALSE         71
##  7 Oscillospiraceae Family TRUE          22
##  8 Oscillospiraceae Family FALSE         25
##  9 Treponema               TRUE          42
## 10 Treponema               FALSE         67
```

Lachnospiraceae Family: 55 ASVs are shared, 108 are not.
Muribaculaceae: 54 are shared, 71 are not.
Treponema: 42 ASVs are shared, 67 are not.
Bacteria Kingdom: 30 ASVs are shared, 125 are not.
Oscillospiraceae Family: 22 ASVs are shared, 25 are not.




### What is the average relative abundance of each genus in the top 5 of shared.otu.genera?
TODO: This might be the correction of our paper


``` r
ps.q.agg.relab.nmr%>%
  # the subset command relies on the fact that we sorted the shared.otu.genera
  # by the number of ASVs. So, the unique(shared.otu.genera$Genus)[1:5] has 
  # genera with the highest number of ASVs
  filter(OTU%in%subset(shared.otu.genera, 
                       Genus %in% unique(shared.otu.genera.cumsum$Genus)[1:5])$OTU)%>%
  group_by(Genus,Sample)%>%
  summarise(TotalGenus=sum(RelativeAbundance))%>%
  group_by(Genus)%>%
  summarise(MeanRelAbGenus = mean(TotalGenus),
            sdRelAbGenus = sd(TotalGenus))%>%
  arrange(-MeanRelAbGenus)
```

```
## `summarise()` has grouped output by 'Genus'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 5 × 3
##   Genus                   MeanRelAbGenus sdRelAbGenus
##   <chr>                            <dbl>        <dbl>
## 1 Muribaculaceae                  10.6          3.17 
## 2 Lachnospiraceae Family           1.93         1.91 
## 3 Treponema                        1.91         2.26 
## 4 Oscillospiraceae Family          0.426        0.237
## 5 Bacteria Kingdom                 0.282        0.201
```



## Plot ASVs in NMR. 
Setup sample levels for NMR for barplots.


``` r
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  dplyr::select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))
```

Bar plot of abundances: high variability


``` r
ps.q.agg.relab.nmr%>%
  left_join(custom.md.ages)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  filter(OTU%in%subset(shared.otu.genera,
                       Genus%in%unique(shared.otu.genera.cumsum$Genus)[1:5])$OTU)%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=Genus))+
  geom_bar(stat="identity")+
  labs(x = "Sample")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
```

```
## Joining with `by = join_by(Sample, class, agegroup, old_agegroup)`
```

<img src="002-summary-stats-qiime2_files/figure-html/unnamed-chunk-91-1.png" alt="" width="1056" />



### 10 most abundant shared ASVs on average account for 30-40% of samples (barplot).
Find the most abundant ASVs on average.


``` r
top10.asv.average<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,MeanRelativeAbundance)%>%
  head(n=10)
top10.asv.average%>%
  arrange(Genus)
```

```
## # A tibble: 10 × 3
## # Groups:   OTU [10]
##    Genus                      OTU                          MeanRelativeAbundance
##    <chr>                      <chr>                                        <dbl>
##  1 Allobaculum                Allobaculum_ASV_22                            3.28
##  2 Erysipelotrichaceae Family Erysipelotrichaceae Family_…                  4.15
##  3 Eubacteriaceae Family      Eubacteriaceae Family_ASV_16                  1.74
##  4 Eubacteriaceae Family      Eubacteriaceae Family_ASV_4                   1.67
##  5 Fibrobacter                Fibrobacter_ASV_2                             4.16
##  6 Prevotella                 Prevotella_ASV_7                              2.35
##  7 Prevotellaceae Family      Prevotellaceae Family_ASV_31                  2.48
##  8 Prevotellaceae_UCG-003     Prevotellaceae_UCG-003_ASV_…                  2.68
##  9 Prevotellaceae_UCG-003     Prevotellaceae_UCG-003_ASV_2                  1.53
## 10 p-251-o5                   p-251-o5_ASV_2                                6.31
```

Bar plot shows the 10 most abundant shared ASVs on average account for 30-40% of samples


``` r
ps.q.agg.relab.nmr%>%
  filter(OTU%in%top10.asv.average$OTU)%>%
  mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=new_OTU))+
  geom_bar(stat="identity")+
  facet_grid(~agegroup,
             scales="free",
             space = "free")+
  labs(x = "Sample")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
```

<img src="002-summary-stats-qiime2_files/figure-html/unnamed-chunk-94-1.png" alt="" width="1056" />

``` r
### 17.2 Most abundant shared ASVs in each age group. ####
```


### Most abundant shared ASVs in each age group.


``` r
top.shared.asvs.by_age<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU,agegroup)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,agegroup,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
  ungroup()
head(top.shared.asvs.by_age)
```

```
## # A tibble: 6 × 5
##   Genus              OTU   agegroup MeanRelativeAbundance MeanRelativeAbundanc…¹
##   <chr>              <chr> <chr>                    <dbl>                  <dbl>
## 1 p-251-o5           p-25… agegrou…                  6.31                   6.66
## 2 p-251-o5           p-25… agegrou…                  6.31                   4.02
## 3 Fibrobacter        Fibr… agegrou…                  4.16                   3.90
## 4 Fibrobacter        Fibr… agegrou…                  4.16                   5.84
## 5 Erysipelotrichace… Erys… agegrou…                  4.15                   3.73
## 6 Erysipelotrichace… Erys… agegrou…                  4.15                   6.82
## # ℹ abbreviated name: ¹​MeanRelativeAbundanceAgegroup
```

``` r
top10.shared.asv.young<-top.shared.asvs.by_age%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  head(n=10)

top10.shared.asv.old<-top.shared.asvs.by_age%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  head(n=10)
```

Are top 10 most abundant ASVs same in two age groups?


``` r
setequal(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)
```

```
## [1] FALSE
```

No. How many are common?


``` r
intersect(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)
```

```
## [1] "p-251-o5_ASV_2"                    "Fibrobacter_ASV_2"                
## [3] "Erysipelotrichaceae Family_ASV_11" "Allobaculum_ASV_22"               
## [5] "Prevotellaceae_UCG-003_ASV_13"     "Prevotellaceae Family_ASV_31"     
## [7] "Prevotella_ASV_7"
```

Seven ASVs

Union of the top 10 ASVs in each of the two age groups.


``` r
top10.shared.asv.union<-sort(union(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU))
```

Make the names shorter for the barplot


``` r
top10.shared.asv.union<-top10.shared.asv.union%>%
  as_tibble()%>%
  rename("OTU"="value")%>%
  mutate(new_OTU = OTU,
         new_OTU = gsub("Allobaculum_", "Allob. ",new_OTU),
         new_OTU = gsub("Erysipelotrichaceae Family_", "Erysip. F. ",new_OTU),
         new_OTU = gsub("Eubacteriaceae Family_", "Eubac. F. ",new_OTU),
         new_OTU = gsub("Fibrobacter_", "Fibrob. ",new_OTU),
         new_OTU = gsub("Muribaculaceae_", "Murib. ",new_OTU),
         new_OTU = gsub("Paludibacteraceae Family_", "Palud. F. ",new_OTU),
         new_OTU = gsub("o5_ASV", "o5 ASV",new_OTU),
         new_OTU = gsub("Prevotella_", "Prev. ",new_OTU),
         new_OTU = gsub("Prevotellaceae Family_", "Prevot. F. ",new_OTU),
         new_OTU = gsub("Prevotellaceae_UCG", "Prevot. UCG",new_OTU),
         new_OTU = gsub("UCG-001_", "UCG-001 ",new_OTU),
         new_OTU = gsub("UCG-003_", "UCG-003 ",new_OTU)
         )
```

Prepare a custom fill with Polychrome package


``` r
set.seed(1)
otu.fill<-createPalette(nrow(top10.shared.asv.union),
                        seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(otu.fill)<-top10.shared.asv.union$new_OTU
```



### Barplot of the most abundant ASVs.


``` r
top10.shared.asv.plot<-ps.q.agg.relab.nmr%>%
  left_join(custom.md.ages)%>%
  filter(OTU%in%top10.shared.asv.union$OTU) %>%
  # mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  left_join(top10.shared.asv.union)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=new_OTU))+
  geom_bar(stat="identity")+
  scale_fill_manual(labels=names(otu.fill),
                    values=otu.fill)+
  theme_bw()+
  coord_cartesian(expand = FALSE)+
  labs(x="Sample",
       y="Relative abundance (%)",
       # title="Top 10 most abundant ASVs across age",
       fill="ASV"
       )+
  mytheme +
  theme(
    legend.title = element_text(size=10),
    legend.text = element_text(size=9),
    legend.key.size = unit(0.3, 'cm'), #change legend key size
    legend.key.spacing.y = unit(0.1, "lines"), # distant between key text
    legend.box.spacing = unit(0.1,"lines"), # space between the plot and the legend box
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle=45,size=10,hjust=1),# rotate 
    strip.text.x = ggtext::element_markdown(size=10),
    panel.spacing = unit(0.8, "cm"), # increase distance between facets
    plot.title = element_text(size = 8), # size of plot title
    plot.caption = element_text(size=8), # size of plot caption
    legend.position = "bottom")
```

```
## Joining with `by = join_by(Sample, class, agegroup, old_agegroup)`
## Joining with `by = join_by(OTU)`
```

``` r
print(top10.shared.asv.plot+
  ggtitle("Top 10 most abundant ASVs across age")+
  theme(plot.title = element_text(size = 14)))
```

<img src="002-summary-stats-qiime2_files/figure-html/unnamed-chunk-103-1.png" alt="" width="768" />

``` r
# for(image.format in image.formats){
#   ggsave(file.path("./images/barplots",
#                  paste0(paste(format(Sys.time(),format="%Y%m%d"),
#                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                        "-top10-asv", ".",image.format)),
#        plot=top10.shared.asv.plot,
#        width=8, height=6,units="in",
#        # width = 5000,height = 3500,
#        # units = "px",
#        dpi=300,device = image.format)
# }
```



### M40 sample is very different. 


``` r
m40.asvs<-ps.q.agg.relab.nmr%>% 
  filter(Sample=="M40")%>%
  dplyr::select(OTU,Genus,RelativeAbundance,MeanRelativeAbundance,
         MeanRelativeAbundanceAgegroup)%>%
  arrange(-RelativeAbundance)%>%
  head(n=10)

set.seed(1)
m40.otu.fill<-createPalette(length(m40.asvs$OTU),
                            seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                            range=c(30, 80))
names(m40.otu.fill)<-sort(m40.asvs$OTU)
m40.asv.plot<-ps.q.agg.relab.nmr%>%
  left_join(custom.md.ages)%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample),
         NewSample=paste0(Sample," (",age,")"),
         NewSample=factor(NewSample,levels=sample.levels$NewSample))%>%
  filter(OTU%in%m40.asvs$OTU)%>%
  # mutate(new_OTU=paste0(OTU," (",Genus,")"))%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=OTU))+
  geom_bar(stat="identity")+
  scale_fill_manual(labels=names(m40.otu.fill),
                    values=m40.otu.fill)+
  theme_bw()+
  coord_cartesian(expand = FALSE)+
  labs(x="Sample",
       y="Relative abundance (%)",
       fill="ASV",
       title="Top 10 most abundant ASVs in M40 sample")+
  guides(fill = guide_legend(ncol = 3))+
  theme(
    legend.title = element_text(size=10),
    legend.text = element_text(size=9),
    legend.key.size = unit(0.3, 'cm'), #change legend key size
    legend.key.spacing.y = unit(0.1, "lines"), # distant between key text
    legend.box.spacing = unit(0.1,"lines"), # space between the plot and the legend box
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle=45,size=10,hjust=1),# rotate 
    strip.text.x = ggtext::element_markdown(size=10),
    panel.spacing = unit(0.8, "cm"), # increase distance between facets
    plot.title = element_text(size = 14), # size of plot title
    plot.caption = element_text(size=8), # size of plot caption
    legend.position = "bottom")
```

```
## Joining with `by = join_by(Sample, class, agegroup, old_agegroup)`
```

``` r
print(m40.asv.plot)
```

<img src="002-summary-stats-qiime2_files/figure-html/unnamed-chunk-106-1.png" alt="" width="960" />

``` r
# ggsave(file.path("./images/barplots",
#                  paste(paste(format(Sys.time(),format="%Y%m%d"),
#                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                        "top10-asv-m40.png")),
#        plot=last_plot(),
#        width=8, height=6,units="in",
#        units = "px",dpi=300,device = "png")
```



## Import the rarefied dataframe 


``` r
# Between NMR (ASV level)
# ps.q.df.preprocessed.date_time<-"20240524_13_58_11"
ps.q.df.preprocessed.date_time<-"20260211_17_14_21"
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables/",authorname,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      paste0("ps.q.df.","rare"),"nonfiltered","OTU",
      paste("NMR",collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)
```

Or between species (Genus level)


``` r
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables/",authorname,paste0(
    paste(
      # "20240426_22_00_04",
      "20260211_17_14_18",
      paste0("ps.q.df.","rare"),"nonfiltered","Genus",
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)
```



### Let's find which major ASVs are specific to one age group ####


``` r
young.ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

old.ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)

length(unique(young.ps.q.agg.relab.nmr$OTU))
```

```
## [1] 1745
```

``` r
length(unique(old.ps.q.agg.relab.nmr$OTU))
```

```
## [1] 771
```



### How many genera are shared between two age groups ####


``` r
genera.young<-ps.q.agg.genus.relab.nmr%>%
  # add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  # left_join(custom.md.ages)%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)
genera.old<-ps.q.agg.genus.relab.nmr%>%
  # add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  # left_join(custom.md.ages)%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(Genus,.keep_all = T)%>%
  arrange(-MeanRelativeAbundanceAgegroup)

length(genera.young$Genus)
```

```
## [1] 229
```

``` r
length(genera.old$Genus)
```

```
## [1] 157
```

``` r
shared.genera<-intersect(genera.young$Genus,genera.old$Genus)
length(shared.genera)
```

```
## [1] 155
```



### Show 5 genera found in old but not young NMR ####


``` r
unique.to_old.genera<-setdiff(genera.old$Genus,genera.young$Genus)
ps.q.agg.genus.relab.nmr%>%
  # add_agegroup_to_tax_df(.,"Genus",custom.md.ages)%>%
  filter(Genus %in% unique.to_old.genera)%>%
  dplyr::select(Sample,Family,Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup )
```

```
## # A tibble: 3 × 5
##   Sample Family          Genus      MeanRelativeAbundance MeanRelativeAbundanc…¹
##   <chr>  <chr>           <chr>                      <dbl>                  <dbl>
## 1 G18    Lachnospiraceae Roseburia               0.00141                 0.0106 
## 2 G14    Veillonellaceae Veillonel…              0.000434                0.00326
## 3 H15    Veillonellaceae Veillonel…              0.000434                0.00326
## # ℹ abbreviated name: ¹​MeanRelativeAbundanceAgegroup
```

``` r
# How much % do common genera take in each age group on average?
ps.q.agg.relab.nmr%>%
  filter(Genus%in%shared.genera)%>%
  group_by(Sample,agegroup)%>%
  summarise(SumRelAbCommonASV=sum(RelativeAbundance))%>%
  arrange(agegroup)%>%
  group_by(agegroup)%>%
  summarise(MeanRelAbCommonASVTotalAge=mean(SumRelAbCommonASV))
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 2 × 2
##   agegroup      MeanRelAbCommonASVTotalAge
##   <chr>                              <dbl>
## 1 agegroup0_10                        99.6
## 2 agegroup10_16                      100.0
```

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
##  [1] ggrepel_0.9.6    ggtext_0.1.2     Polychrome_1.5.4 vegan_2.6-4     
##  [5] lattice_0.22-6   permute_0.9-8    lubridate_1.9.4  forcats_1.0.0   
##  [9] stringr_1.5.1    dplyr_1.1.4      purrr_1.0.4      readr_2.1.5     
## [13] tidyr_1.3.1      tibble_3.2.1     ggplot2_4.0.0    tidyverse_2.0.0 
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.6         xfun_0.56            bslib_0.10.0        
##  [4] tzdb_0.5.0           vctrs_0.6.5          tools_4.4.3         
##  [7] generics_0.1.4       parallel_4.4.3       cluster_2.1.8.1     
## [10] pkgconfig_2.0.3      Matrix_1.7-4         RColorBrewer_1.1-3  
## [13] S7_0.2.0             scatterplot3d_0.3-44 lifecycle_1.0.5     
## [16] compiler_4.4.3       farver_2.1.2         litedown_0.9        
## [19] htmltools_0.5.8.1    sass_0.4.10          yaml_2.3.12         
## [22] pillar_1.11.1        jquerylib_0.1.4      MASS_7.3-65         
## [25] cachem_1.1.0         nlme_3.1-167         commonmark_2.0.0    
## [28] tidyselect_1.2.1     digest_0.6.37        stringi_1.8.4       
## [31] bookdown_0.46        labeling_0.4.3       splines_4.4.3       
## [34] fastmap_1.2.0        grid_4.4.3           colorspace_2.1-2    
## [37] cli_3.6.4            magrittr_2.0.3       utf8_1.2.4          
## [40] withr_3.0.2          scales_1.4.0         timechange_0.3.0    
## [43] rmarkdown_2.30       otel_0.2.0           hms_1.1.4           
## [46] evaluate_1.0.5       knitr_1.51           viridisLite_0.4.3   
## [49] markdown_2.0         mgcv_1.9-1           rlang_1.1.5         
## [52] gridtext_0.1.5       Rcpp_1.0.14          glue_1.8.0          
## [55] xml2_1.3.8           rstudioapi_0.18.0    jsonlite_2.0.0      
## [58] R6_2.6.1
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 2968951 158.6    4379735 234.0  4379735 234.0
## Vcells 5240330  40.0   21892305 167.1 21324459 162.7
```

