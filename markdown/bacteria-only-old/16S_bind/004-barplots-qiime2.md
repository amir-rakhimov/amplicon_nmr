---
output: 
  bookdown::html_document2:
     toc: true
---




# Barplots





## Introduction
After importing the QZA files into R and processing them with phyloseq,
it's time to explore the taxonomic composition of our data.
We will use the `Polychrome` package to create a custom palette for the 
barplots.
We will build a combined barplot for all hosts.





## Load necessary libraries.


``` r
# install.packages(c("tidyverse","ggtext","Polychrome"))
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
library(Polychrome)
library(ggtext)
```



## Specifying parameters and directory/file names. 
Name of the folder with QIIME2 output:


``` r
authorname<-"pooled" 
```

The taxonomic rank that was used for agglomeration:


``` r
agglom.rank<-"Genus"
```

Truncation level that we chose in QIIME2:


``` r
truncationlvl<-"234" 
```

Single reads or paired reads: decided in QIIME2.


``` r
read.end.type<-"single"
```

Specify paths and image formats:


``` r
barplot.directory<-"./images/barplots/" 
rdafiles.directory<-"./output/rdafiles"
```

Path for custom metadata:


``` r
metadata.directory<-"./output/rdafiles" 
image.formats<-c("png","tiff")
```

Import abundance table as an rds file (NOT rarefied): 


``` r
ps.q.agg.date_time<-"20260211_17_01_10"
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))

# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
```

Import metadata:


``` r
custom.md<-readRDS(file.path(metadata.directory,"custom.md.rds"))
custom.levels<-c("NMR","B6mouse","MSMmouse","FVBNmouse",
                 "DMR","hare","rabbit","spalax","pvo")
```

Set "pretty" labels for barplot facets that correspond to animal hosts. 
Here, the left side of the vector (values) is taken from the metadata, 
while the right side (names) are the pretty labels that will be shown on 
the final barplot.


``` r
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "DMR" = "*Fukomys damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N<br>mouse",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*")
```



## Set up the barplot.
First, find where the agglom.rank column is located (what number the
`Genus` columns is located at, for example). In case of ASV table, the 
column is actually called `OTU`.  
We also need to know the location of the preceding taxonomic rank 
(e.g. `Family`).
Usually, it's to the left of the agglom.rank column.


``` r
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}
custom_order <- c("Kingdom", "Phylum", "Class", "Order", "Family",
                  "Genus","Species")
```

If we agglomerate by higher level (Order, Class, etc), need to adjust 
the rank by
removing the levels below. In case of Order-level data, remove the 
`Family`, and `Genus`


``` r
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}
```






### Tidy up the taxonomic names. ####
First, merge the taxonomic columns into the "taxon" column to find 
unclassified taxa.
If a taxon has "Kingdom", "Phylum", "Class", "Order", etc. in one of those 
columns, this taxon is unclassified. These will be left as they are 
("taxon" value becomes the agglom.rank), but classified taxa 
should have the preceding rank in the brackets (agglom.rank and 
preceding.rank). 
For example, "Prevotella (Prevotellaceae)" is classified, 
while "Lachnospiraceae Family" is not. The barplot needs to reflect this. At 
this point, the "taxon" column is changed. Moreover, if the average 
relative abundance of a taxon is less than 1%, the "taxon" value becomes 
"Remainder (mean relative abundance < 1%)" to avoid too many names in the 
legend.

Unclassified taxa need to be separate from classified ones: 
we put unclassified ones at the top of the 
legend. Then, we sort both groups alphabetically.  
Another thing we need to do is wrapping long names using `str_wrap`. We 
convert the "taxon" column into factor and set the levels as follows: 
the first level is the "Remainder (mean relative abundance < 1%)" string. 
Next, unclassified taxa, then classified taxa. 
We use the str_detect to find unclassified taxa because these have 
Kingdom or Phylum in the name.


``` r
new.df<-ps.q.agg%>%
  unite("taxon",Kingdom:all_of(agglom.rank),sep = ";",remove = FALSE)%>%
  mutate(is_unclassified = grepl
         ("Kingdom|Phylum|Class|Order|Family|Genus|Species",taxon))%>% # add column to show that a taxon was unclassified
  mutate(taxon=ifelse(is_unclassified,
                       get(agglom.rank),
                       paste0(get(agglom.rank), " (",get(preceding.rank),
                              ")")),
         taxon=ifelse(MeanRelativeAbundance<1,"Remainder (mean relative abundance < 1%)",taxon))%>% # if a taxon is classified, we change it to show the preceding rank, e.g Family (Genus)
  group_by_at(c("is_unclassified",preceding.rank))%>%
  arrange(desc(is_unclassified),
          get(preceding.rank),
          taxon,
          .by_group = TRUE)%>% # sort within unclassified and classified taxa (by group)
  ungroup()%>%
  mutate(taxon = gsub("_", " ", taxon),
         taxon = str_wrap(taxon,width=32))%>%
  mutate(taxon=factor(taxon,levels=unique(taxon)),
         taxon=taxon%>%
           fct_relevel(c("Remainder (mean relative\nabundance < 1%)",
                         unique(unlist(lapply(custom_order,
                                              function(pat){levels(.)[str_detect(levels(.),
                                                                                 fixed(pat))]
                                                            })))
                         ),
                       after=0))%>% # change taxon column to factor. Then, we set "Remainder" as the first level, unclassified taxa next, then classified taxa. we use the str_detect to find unclassified taxa because these have Kingdom or Phylum in the name
  ungroup()
```

Create a color palette


``` r
set.seed(1)
plot.cols<-createPalette(length(levels(unique(new.df$taxon))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(new.df$taxon))
```



## Plot the barplot.
We will concatenate the sample name on the x-axis with the 
number of reads in that sample. It will look like "Sample 1 (n = 25000)". 
Actually, this kind of string will be stored in a new column 
called NewSample.

Then, we will convert the "class" column (host names) into factors that we 
defined in custom.levels. This will order our bars according to the vector 
of levels. It must be factor because this allows us ordering 
custom.levels as we want, otherwise it would be alphabetic. And in 
our vector, the first level is naked mole-rats. So, the first facet 
will also be naked mole-rats. Facets are panels that correspond 
to animal host.


``` r
barplot.all<-new.df%>%
  group_by(Sample)%>%
  mutate(TotalAbundance=sum(Abundance))%>% # add total counts per sample, 
  # so we can have info about sample size
  ungroup()%>%
  mutate(NewSample=paste0(Sample," (n = ", TotalAbundance, ")"))%>% # add a 
  # column where sample names are together with sample sizes
  mutate(class=fct_relevel(class,
                           custom.levels,
                           after = 0))%>%
  # filter(class%in%c("NMR","B6mouse"))%>%
  ggplot(aes(x=NewSample,y=RelativeAbundance,fill=taxon))+
  geom_bar(stat="identity")+
  facet_grid(~class, # separate animal hosts
             scales="free",  # each species will have its own bars inside
             # facet (instead of all bars)
             space = "free",
             labeller = labeller(class=pretty.level.names))+
  theme(legend.position = "bottom")+
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  # custom palette
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon",
       caption="Mean Relative Abundance was calculated for each host separately")+
  theme_bw()+
  guides(fill = guide_legend(ncol = 6))+
  # ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different rodents"))+
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.5), 'cm'),
        axis.text.y = element_text(size=5), # size of y axis ticks
        axis.title = element_text(size = 7), # size of axis names
        strip.text.x = ggtext::element_markdown(size = 5), # the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        axis.text.x = element_text(angle=45,size=6,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        plot.caption = element_text(size=5), # size of plot caption
        plot.title = element_text(size = 5), # size of plot title
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title= element_text(size = 7),
        legend.box.spacing = unit(0.1,"lines"), # space between the plot and the legend box
        legend.key.spacing.y = unit(2, "pt"), # distance between the legend's key text
        legend.key.size = unit(0.3, 'cm'), # change legend key size
        legend.text = element_text(margin = margin(l=1), size = 5), # size of the legend text
        legend.position = "bottom" # place the legend under the plot
  )
```

Information about padding and legend size was obtained from:
https://stackoverflow.com/a/78067650   
Show the plot:


``` r
print(barplot.all+
        ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different rodents"))+
        theme(plot.title = element_text(size=10))
      )
```

<img src="004-barplots-qiime2_files/figure-html/unnamed-chunk-23-1.png" alt="" width="1056" />

``` r
# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "barplot",paste(custom.levels,collapse = '-'),
#                       truncationlvl,agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=barplot.all,
#          width=11, height=8,units="in",
#          
#          # width = 13500,height = 5200,
#          # units = "px",
#          dpi=800,device = image.format)
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
##  [1] ggtext_0.1.2     Polychrome_1.5.4 lubridate_1.9.4  forcats_1.0.0   
##  [5] stringr_1.5.1    dplyr_1.1.4      purrr_1.0.4      readr_2.1.5     
##  [9] tidyr_1.3.1      tibble_3.2.1     ggplot2_4.0.0    tidyverse_2.0.0 
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.10          generics_0.1.4       xml2_1.3.8          
##  [4] stringi_1.8.4        hms_1.1.4            digest_0.6.37       
##  [7] magrittr_2.0.3       evaluate_1.0.5       grid_4.4.3          
## [10] timechange_0.3.0     RColorBrewer_1.1-3   bookdown_0.46       
## [13] fastmap_1.2.0        jsonlite_2.0.0       scales_1.4.0        
## [16] jquerylib_0.1.4      cli_3.6.4            rlang_1.1.5         
## [19] scatterplot3d_0.3-44 litedown_0.9         commonmark_2.0.0    
## [22] withr_3.0.2          cachem_1.1.0         yaml_2.3.12         
## [25] otel_0.2.0           tools_4.4.3          tzdb_0.5.0          
## [28] colorspace_2.1-2     vctrs_0.6.5          R6_2.6.1            
## [31] lifecycle_1.0.5      pkgconfig_2.0.3      pillar_1.11.1       
## [34] bslib_0.10.0         gtable_0.3.6         glue_1.8.0          
## [37] Rcpp_1.0.14          xfun_0.56            tidyselect_1.2.1    
## [40] rstudioapi_0.18.0    knitr_1.51           farver_2.1.2        
## [43] htmltools_0.5.8.1    rmarkdown_2.30       labeling_0.4.3      
## [46] compiler_4.4.3       S7_0.2.0             markdown_2.0        
## [49] gridtext_0.1.5
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##           used (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 1800928 96.2    3468485 185.3  2488569 133.0
## Vcells 3459548 26.4    8388608  64.0  8324211  63.6
```

