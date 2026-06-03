---
output: 
  bookdown::html_document2:
     toc: true
---





# Alpha diversity analysis of QIIME2 output





## Introduction
This script analyses alpha diversity among hosts. We will compare 
hosts using three metrics: the number of observed genera, 
Shannon diversity index, and inverse Simpson index.  





## Load necessary libraries and scripts.


``` r
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan","Polychrome"))
# library(phyloseq)
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



## Specifying parameters and directory/file names. 
Name of the folder with QIIME2 output:


``` r
authorname<-"pooled"
```

Specify paths and image formats:


``` r
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)
image.formats<-c("png","tiff")
```

Truncation level that we chose in QIIME2:


``` r
truncationlvl<-"234"
```

The taxonomic rank that was used for agglomeration:


``` r
agglom.rank<-"Genus"
```

Single reads or paired reads (decided in QIIME2):


``` r
read.end.type<-"single"
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


Import the rarefied abundance table (rds  file):


``` r
ps.q.df.preprocessed.date_time<-"20260211_17_14_19"
# 20260211_17_14_19 rarefied table for all hosts, genus level
# 20260211_17_14_20 rarefied table file for NMR, genus level
# 20260211_17_14_21 rarefied table file for NMR, ASV level
```

Import metadata:


``` r
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
```

Set "pretty" labels


``` r
pretty.level.names<-c("NMR" = "*H. glaber*", # better labels for facets
                      "DMR" = "*F. damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*N. leucodon*",
                      "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
```

Facet labels


``` r
metric.labs <- c('sobs'= case_when(agglom.rank == "OTU" ~ 
                                     "Richness \n(Observed species)",
                                agglom.rank == "Genus" ~ 
                                  "Richness \n(Observed genera)") ,
              'shannon' = "Shannon",
              'invsimpson' = "Inverse \nSimpson")
```

Metrics to plot


``` r
plot.metrics<-c("sobs","shannon","invsimpson")
div.indices<-c("sobs","shannon", "invsimpson")

mytheme<-theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_text(size=15),
  axis.title = element_text(size = 15),
  strip.text.x = element_text(size=15),
  plot.title = element_text(size = 15),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 18),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())
```

Remove rows with 0 Abundance, if there's any.


``` r
ps.q.agg <- ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
```

Import rarefied dataframe and select Sample, Abundance, 
class, and agglom.rank columns.


``` r
ps.q.df.preprocessed<-readRDS(
  file.path(rdafiles.directory,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".rds")))%>%
  filter(class%in%custom.levels, Abundance!=0)%>%
  dplyr::select(Sample,Abundance,class,all_of(agglom.rank),)
```



## Alpha diversity analysis.




### Compute alpha diversity metrics.


``` r
all.div<-ps.q.df.preprocessed%>%
  group_by(Sample)%>%
  reframe(sobs=specnumber(Abundance), # richness (num of species)
          shannon=diversity(Abundance,index = "shannon"),
          invsimpson=diversity(Abundance, index="invsimpson"), # inverse simpson
          tot=sum(Abundance),
          class=class)%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()
```



### Alpha diversity tests. ####
The Kruskal-Wallis rank sum test analyzes differences in each 
alpha diversity metric across animal hosts.


``` r
kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
colnames(kt.results)<-div.indices
rownames(kt.results)<-c("statistic", "pvalue")
```

All unique pairwise combinations:


``` r
combinations<-combn(custom.levels,2) 
w.results<-data.frame(matrix(nrow = ncol(combinations),
                             ncol=length(div.indices))) # ncol(combinations) pairwise comparisons
colnames(w.results)<-div.indices
w.results<-array(dim = c(length(table(combinations[1,])), # nrows
                         length(table(combinations[2,])), # ncols
                         length(div.indices)), # num of 2D arrays (stacking) 
                 dimnames = list(NULL, NULL, div.indices))



for (div.metric in div.indices) {
  metric.ds<-all.div%>%
    filter(metric==div.metric)%>%
    distinct()
  # perform kruskal-wallis test
  kt<-kruskal.test(value~class,data=metric.ds)
  
  kt.results["statistic",div.metric]<- kt$statistic
  kt.results["pvalue",div.metric]<- kt$p.value
  # low pvalue indicates statistical difference between some groups
  # But which groups are different from others?
  # Perform pairwise wilcoxon test
  
  # NB: You must not run pairwise wilcoxon test unless kruskal wallis results 
  # are significant
  if(kt$p.value<0.05){
    w.test<-pairwise.wilcox.test(metric.ds$value,
                                 metric.ds$class,
                                 p.adjust.method = "BH"
                                 ,exact=FALSE
                                 )
    w.results[,,div.metric]<-w.test$p.value
    
  }else(
    w.results[,,div.metric]<-matrix(data = "n.s.",
                                    nrow = nrow(w.test$p.value),
                                    ncol = ncol(w.test$p.value))
    
  )
  dimnames(w.results)[[1]]<-dimnames(w.test$p.value)[[1]] # change rownames of w.results
  dimnames(w.results)[[2]]<-dimnames(w.test$p.value)[[2]] # change colnames of w.results
}

kt.results
```

```
##                   sobs      shannon   invsimpson
## statistic 7.280958e+01 7.646242e+01 7.374686e+01
## pvalue    1.352624e-12 2.511894e-13 8.787165e-13
```

``` r
w.results
```

```
## , , sobs
## 
##               B6mouse          DMR FVBNmouse         hare     MSMmouse
## DMR       0.691993266           NA        NA           NA           NA
## FVBNmouse 0.691993266 2.057970e-01        NA           NA           NA
## hare      0.236012189 1.740406e-03 0.8273413           NA           NA
## MSMmouse  0.981025330 5.162308e-01 0.1830412 0.0040781266           NA
## NMR       0.016498264 3.560095e-06 0.1320516 0.0312766917 0.0004875759
## pvo       0.638688814 3.262471e-01 0.2360122 0.0204299011 0.6625897294
## rabbit    1.000000000 7.172796e-01 0.4371681 0.0741462942 0.9810253298
## spalax    0.008726373 7.392181e-06 0.0204299 0.0005478033 0.0005478033
##                    NMR          pvo       rabbit
## DMR                 NA           NA           NA
## FVBNmouse           NA           NA           NA
## hare                NA           NA           NA
## MSMmouse            NA           NA           NA
## NMR                 NA           NA           NA
## pvo       3.996673e-04           NA           NA
## rabbit    2.590009e-03 0.6625897294           NA
## spalax    3.929571e-06 0.0004875759 0.0009793901
## 
## , , shannon
## 
##               B6mouse          DMR  FVBNmouse        hare     MSMmouse
## DMR       0.425803710           NA         NA          NA           NA
## FVBNmouse 0.521666471 6.050488e-02         NA          NA           NA
## hare      0.019068304 4.168804e-04 0.63327805          NA           NA
## MSMmouse  0.091437337 5.452710e-03 0.82250895 0.007101269           NA
## NMR       0.005452710 6.815077e-07 0.02432257 0.000441647 0.0002152878
## pvo       0.814132811 1.170594e-01 0.79599888 0.096156624 0.7245959403
## rabbit    0.822508946 3.850466e-01 1.00000000 0.795998883 0.9811102146
## spalax    0.007626574 7.499778e-06 0.01938005 0.000441647 0.0004416470
##                    NMR          pvo      rabbit
## DMR                 NA           NA          NA
## FVBNmouse           NA           NA          NA
## hare                NA           NA          NA
## MSMmouse            NA           NA          NA
## NMR                 NA           NA          NA
## pvo       1.608806e-04           NA          NA
## rabbit    1.014765e-01 0.7959988831          NA
## spalax    3.957312e-06 0.0002152878 0.006747595
## 
## , , invsimpson
## 
##               B6mouse          DMR  FVBNmouse         hare     MSMmouse
## DMR       1.000000000           NA         NA           NA           NA
## FVBNmouse 0.798870448 3.807495e-01         NA           NA           NA
## hare      0.019068304 6.682285e-04 0.09446310           NA           NA
## MSMmouse  0.063951580 1.044555e-02 0.88195688 0.0271745461           NA
## NMR       0.005948411 5.944136e-07 0.01479553 0.0004907189 0.0002469123
## pvo       0.415935888 2.811976e-01 0.98750821 0.1718127525 0.8819568759
## rabbit    0.818487813 7.988704e-01 1.00000000 0.8508953578 0.9639463193
## spalax    0.008799893 7.499778e-06 0.01938005 0.0004907189 0.0004907189
##                    NMR          pvo      rabbit
## DMR                 NA           NA          NA
## FVBNmouse           NA           NA          NA
## hare                NA           NA          NA
## MSMmouse            NA           NA          NA
## NMR                 NA           NA          NA
## pvo       2.469123e-04           NA          NA
## rabbit    2.673735e-01 0.9639463193          NA
## spalax    3.957312e-06 0.0002469123 0.008799893
```

``` r
stopifnot(all(kt.results[2,]<0.05))

# convert multidimensional array into data frame
# w.results.df<-w.results%>%
#   as.data.frame.table()%>%
#   drop_na()
# colnames(w.results.df)<-c("class1","class2","div.metric","p.value")

max.values<-all.div%>%
  group_by(metric)%>%
  summarise(max_val=max(value))

all.div$class<-factor(all.div$class,levels=custom.levels)
```



## Plot alpha diversity metrics. ####


``` r
div.plot<-all.div%>%
  filter(metric%in%plot.metrics)%>%
  mutate(class = factor(class, levels = custom.levels))%>%
  group_by(class)%>%
  mutate(class_id = cur_group_id(), # add numeric id
         class_letter = LETTERS[class_id], # then, turn it into a letter
         class_name = paste(class_letter, pretty.level.names[class_id],sep = ": "))%>%
  ggplot(aes(x = class_letter, 
             y = value,
             fill = class_name))+
  geom_boxplot(show.legend = T,
               outliers = F)+ # outliers won't be visible anyway
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             ncol=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+
  # scale_color_manual(breaks = unname(pretty.level.names),
  #                    labels=unname(pretty.level.names))+
  # scale_x_discrete(labels=pretty.level.names,
  #                  limits=custom.levels)+ # rename boxplot labels (x axis)
  scale_fill_viridis_d(option="C")+
  stat_summary(fun=median, geom="point", shape=23, size=3, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  labs(fill = "Host")+
  # ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
  mytheme + # general theme
  theme(axis.text.x = ggtext::element_markdown(size=12),
        legend.text = ggtext::element_markdown(size=15))

div.plot.jitter<-div.plot+
  geom_jitter(aes(colour=class_letter),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=1.8,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
  guides(color = "none")
```

``` r
print(div.plot.jitter + 
        ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14))
)
```

<img src="006-alpha-diversity_files/figure-html/unnamed-chunk-27-1.png" alt="" width="1056" />

``` r
div.plot.dots<-div.plot+
  geom_dotplot(aes(colour=class),
               fill="white",
               color="black",
               binaxis='y',
               stackdir='center', 
               dotsize=0.4) # add dots
```

``` r
print(div.plot.dots + 
        ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14))
)
```

```
## Bin width defaults to 1/30 of the range of the data. Pick better value with
## `binwidth`.
```

<img src="006-alpha-diversity_files/figure-html/unnamed-chunk-28-1.png" alt="" width="1056" />

``` r
# for(image.format in image.formats){
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha-jitter",
#                       paste(plot.metrics,collapse = "-"),
#                       paste(custom.levels,collapse = '-'),
#                       agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=div.plot.jitter,
#          width=11, height=5,units="in",
#          # width = 4500,height = 2000,units = "px",
#          dpi=300,device = image.format)
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha-dots",
#                       paste(plot.metrics,collapse = "-"),
#                       paste(custom.levels,collapse = '-'),
#                       agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=div.plot.dots,
#          width=11, height=5,units="in",
#          # width = 4500,height = 2000,units = "px",
#          dpi=300,device = image.format)
# }
```



### Add significance stars. 


``` r
coord.combs<-combn(seq_along(custom.levels),2)

# First, find significant results (stars)
stars.list<-matrix(NA,nrow=length(div.indices),ncol = ncol(coord.combs))
for (k in 1:dim(w.results)[3]) { # loop over each sub array
  for (i in 1:dim(w.results)[1]) { # then each row
    for (j in 1:dim(w.results)[2]) { # then each column
      if (!is.na(w.results[i,j,k])) { # extract non-NA values
        ith.row<-rownames(w.results[,,k])[i] 
        jth.col<-colnames(w.results[,,k])[j]
        w.val<- w.results[i,j,k]
        
        # e.g. ith.row="NMR"
        # jth.col="hare"
        # find in custom.levels the positions that correspond to c("NMR","hare")
        # the result: level.position = c(1,4)
        levels.position<-which(custom.levels%in%c(ith.row,jth.col))
        # find these positions in the coord.combs to get the column 
        # that stores levels.position
        # result: level.col=3 
        level.col<-which(apply(coord.combs,2,function(x) 
          return(all(x==levels.position))))
        # assign the w.val to stars.list:kth row is diversity metric
        # level.col column is the combination of levels
        # we must assign either "*" or "n.s." depeding on w.result[i,j,k]
        stars.list[k,level.col]<-ifelse(w.val!="n.s.",ifelse(w.val<0.05,"*","n.s."),w.val)
      }
    }
  }
}

stars.list<-as.vector(t(stars.list))

stars<-tibble(
  metric=factor(rep(plot.metrics,each=ncol(coord.combs)), # each plot metric repeated by the number of combinations
                levels = plot.metrics), # names of our metrics
  label=stars.list
)
star.indices<-which(stars$label=="*") # only significant results

freqs<-as.data.frame(table(stars))%>%filter(label=="*")
yvalues<-c()
for (i in seq_along(freqs)){
  vec<-seq(from=1, by=0.08,length.out=freqs[i,"Freq"])*
    as.numeric(max.values[max.values$metric==freqs[i,"metric"],"max_val"])
  yvalues<-c(yvalues,vec)
}

# the horizontal lines
# merge two vectors: two rows
# c() turns them into a vector

# y coordinates
# multiply vectors to get a matrix: nrows is the number of comparisons
# ncols is the number of plot.metrics
# seq is the sequence of values that will multiply maxvalues

# start and end values can be found from pairwise combinations
# x values have dimensions: num of metrics * num of pairwise comparisons
xvalues<-rep(coord.combs[1,],
             length=length(div.indices)*ncol(coord.combs)) # first row is x start values
xendvalues<-rep(coord.combs[2,],
                length=length(div.indices)*ncol(coord.combs)) # second row is x end values

# select only significant results
xvalues<-xvalues[star.indices]
xendvalues<-xendvalues[star.indices]

horizontal.lines<-tibble(
  metric=factor(rep(plot.metrics,each=ncol(coord.combs)),
                levels = plot.metrics)[star.indices], # names of our metrics
  x = xvalues,       #1 2 1 1 2 1 1 2 1 1 2 1
  xend = xendvalues, #2 3 3 2 3 3 2 3 3 2 3 3
  y =yvalues,
  yend = yvalues
)

# labels depend on our tests (statistical significance)
xstars<-(xvalues+xendvalues)/2
stars<-stars%>% filter(label=="*")%>%
  mutate(x= xstars,
         y =c(yvalues)*1.01)

newplot<-div.plot+
  geom_segment(data=horizontal.lines, # add horizontal.lines of significance
               aes(x=x, xend=xend, y=y, yend=yend),
               inherit.aes = FALSE)+# no conflict with different fills
  geom_text(data = stars,aes(x=x, y=y, label=label),
            inherit.aes = FALSE,size=10) # add stars 
```

``` r
print(newplot)
```

<img src="006-alpha-diversity_files/figure-html/unnamed-chunk-31-1.png" alt="" width="1056" />

``` r
# save plot with significance bars
# Use the table of wilcoxon tests, check if there are any pairwise comparisons
# that were significant
# if yes, save the plot
# if(table(w.results<0.05)[2]>0){
#   for(image.format in image.formats){
#     ggsave(paste0("./images/diversity/alpha/",
#                   paste(paste(format(Sys.time(),format="%Y%m%d"),
#                               format(Sys.time(),format = "%H_%M_%S"),
#                               sep = "_"),
#                         "alpha",
#                         paste(plot.metrics,collapse = "-"),
#                         paste(custom.levels,collapse = '-'),
#                         agglom.rank,truncationlvl,"signif",
#                         sep = "-"),".",image.format),
#            plot=newplot,
#            width=11, height=8,units="in",
#            # width = 4500,height = 2000,units = "px",
#            dpi=300,device = image.format)
#   }
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
##  [1] vegan_2.6-4     lattice_0.22-6  permute_0.9-8   lubridate_1.9.4
##  [5] forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4    
##  [9] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_4.0.0  
## [13] tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.10        generics_0.1.4     xml2_1.3.8         stringi_1.8.4     
##  [5] hms_1.1.4          digest_0.6.37      magrittr_2.0.3     evaluate_1.0.5    
##  [9] grid_4.4.3         timechange_0.3.0   RColorBrewer_1.1-3 bookdown_0.46     
## [13] fastmap_1.2.0      Matrix_1.7-4       jsonlite_2.0.0     ggtext_0.1.2      
## [17] mgcv_1.9-1         viridisLite_0.4.3  scales_1.4.0       jquerylib_0.1.4   
## [21] cli_3.6.4          rlang_1.1.5        litedown_0.9       commonmark_2.0.0  
## [25] splines_4.4.3      withr_3.0.2        cachem_1.1.0       yaml_2.3.12       
## [29] otel_0.2.0         tools_4.4.3        parallel_4.4.3     tzdb_0.5.0        
## [33] vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.5    MASS_7.3-65       
## [37] cluster_2.1.8.1    pkgconfig_2.0.3    pillar_1.11.1      bslib_0.10.0      
## [41] gtable_0.3.6       Rcpp_1.0.14        glue_1.8.0         xfun_0.56         
## [45] tidyselect_1.2.1   rstudioapi_0.18.0  knitr_1.51         farver_2.1.2      
## [49] nlme_3.1-167       htmltools_0.5.8.1  labeling_0.4.3     rmarkdown_2.30    
## [53] compiler_4.4.3     S7_0.2.0           markdown_2.0       gridtext_0.1.5
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 2860997 152.8    4347284 232.2  4347284 232.2
## Vcells 5015809  38.3   10146329  77.5  8002247  61.1
```

