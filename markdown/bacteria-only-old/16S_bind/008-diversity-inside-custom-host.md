---
output: 
  bookdown::html_document2:
     toc: true
---




# Alpha and beta diversity in NMR samples 




## Introduction
This script performs alpha diversity analysis and PCA on NMR samples




## Load necessary libraries.


``` r
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-4
```

``` r
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   4.0.0     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.0.4
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ purrr::%||%()   masks base::%||%()
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(ggrepel)
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
rare.status<-"rare"
filter.status<-"nonfiltered"
```

The taxonomic rank that was used for agglomeration:


``` r
agglom.rank<-"OTU"
```

Import rarefied data (rds).


``` r
ps.q.df.preprocessed.date_time<-"20260211_17_14_21" # ASV NMR
```

Import abundance table as an rds file (NOT rarefied): 


``` r
ps.q.agg.date_time<-"20260211_17_01_07" # ASV
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))
```

Import metadata.


``` r
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
```

Specify paths and image formats.


``` r
image.formats<-c("png","tiff")
```

Labels for plots:


``` r
metric.labs=c('sobs'= case_when(agglom.rank == "OTU" ~ "Richness (Observed species)",
                                agglom.rank == "Genus" ~ "Richness (Observed genera)") ,
              'shannon' = "Shannon",
              'invsimpson' = "Inverse Simpson")
plot.metrics<-c("sobs","shannon",
                "invsimpson") # metrics to plot
div.indices<-c("sobs","shannon",
               "invsimpson")
```

If you compare NMRs or mice (inside host):


``` r
host<-"NMR"
# host.class<-c("NMR"="naked mole-rat",
#               "mice"="mouse")
# comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"
gg.title.taxon <- ifelse(agglom.rank=="OTU","(ASV level)",
                       paste0("(",agglom.rank," level)"))

host.labels <- c("NMR" = "*Heterocephalus glaber*")
```

Load the rarefied dataframe:


``` r
ps.q.df.preprocessed<-readRDS(
  file.path(rdafiles.directory,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      paste0("ps.q.df.",rare.status),filter.status,agglom.rank,
      paste(names(host.labels),collapse = '-'),sep = "-"),
    ".rds")))
```

Load metadata with age information:


``` r
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
```

Theme for alpha diversity plots:


``` r
alpha.plot.theme<-theme(#plot.margin=unit(c(1,1,1,2), 'cm'),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = ggtext::element_markdown(size=18),
  axis.text.y = element_text(size=18),
  axis.title = element_text(size = 20),
  strip.text.x = element_text(size=20),
  plot.title = element_text(size = 27),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 25),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())
```

Theme for PCA scatter plots:


``` r
scatter.plot.theme<-theme(
  axis.text.x = element_text(angle=0,size=20,hjust=1),
  axis.text.y = element_text(size=20),
  axis.title = element_text(size = 20),
  plot.title = element_text(size = 27),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 25),
  legend.position = "right",
  plot.caption = ggtext::element_markdown(hjust = 0, size=20),
  plot.caption.position = "plot",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())
```



## Alpha diversity analysis. 




### Compute alpha diversity metrics.


``` r
all.div<-ps.q.df.preprocessed%>%
  group_by(Sample)%>%
  reframe(sobs=specnumber(Abundance), # richness (num of species)
          shannon=diversity(Abundance,index = "shannon"),
          # simpson=diversity(Abundance, index="simpson"),
          invsimpson=diversity(Abundance, index="invsimpson"), # inverse simpson
          tot=sum(Abundance))%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()
```



### Alpha diversity tests. 
Function that runs Kruskal-Wallis test and pairwise Wilcoxon tests:


``` r
test_alpha_diversity<-function(div.indices, 
                               all.div, 
                               comparison,
                               custom.levels,
                               custom.md){
  kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
  colnames(kt.results)<-div.indices
  rownames(kt.results)<-c("statistic", "pvalue")
  combinations<-combn(custom.levels,2) # all unique pairwise combinations
  w.results<-data.frame(matrix(nrow = ncol(combinations),
                               ncol=length(div.indices))) # ncol(combinations) pairwise comparisons
  colnames(w.results)<-div.indices
  w.results<-array(dim = c(length(table(combinations[1,])), # nrows
                           length(table(combinations[2,])), # ncols
                           length(div.indices)), # num of 2D arrays (stacking) 
                   dimnames = list(NULL, NULL, div.indices))
  for (div.metric in div.indices) {
    metric.ds<-all.div%>%
      left_join(custom.md)%>%
      filter(metric==div.metric)%>%
      distinct()
    # perform kruskal-wallis test
    if(comparison=="age"){
      kt<-kruskal.test(value~agegroup,data=metric.ds)
      
    }else if (comparison=="sex"){
      kt<-kruskal.test(value~sex,data=metric.ds)
      
    }else if(comparison=="strain"){
      kt<-kruskal.test(value~class,data=metric.ds)
      
    }
    
    kt.results["statistic",div.metric]<- kt$statistic
    kt.results["pvalue",div.metric]<- kt$p.value
    # low pvalue indicates statistical difference between some groups
    # But which groups are different from others?
    # Perform pairwise wilcoxon test
    
    # NB: You must not run pairwise wilcoxon test unless kruskal wallis results 
    # are significant
    if(kt$p.value<0.05){
      if(comparison=="age"){
        w.test<-pairwise.wilcox.test(metric.ds$value,
                                     metric.ds$agegroup,
                                     p.adjust.method = "BH",
                                     exact=FALSE)
      }else if (comparison=="sex"){
        w.test<-pairwise.wilcox.test(metric.ds$value,
                                     metric.ds$sex,
                                     p.adjust.method = "BH",
                                     exact=FALSE)
        
      }else if(comparison=="strain"){
        w.test<-pairwise.wilcox.test(metric.ds$value,
                                     metric.ds$class,
                                     p.adjust.method = "BH",
                                     exact=FALSE)
      }
      
      w.results[,,div.metric]<-w.test$p.value
      
    }else(
      w.results[,,div.metric]<-matrix(data = "n.s.",
                                      nrow = length(table(combinations[1,])), # nrows,
                                      ncol = length(table(combinations[2,])))
                                      # nrow = nrow(w.test$p.value),
                                      # ncol = ncol(w.test$p.value))
      
    )
    # dimnames(w.results)[[1]]<-dimnames(w.test$p.value)[[1]] # change rownames of w.results
    # dimnames(w.results)[[2]]<-dimnames(w.test$p.value)[[2]] # change colnames of w.results
    dimnames(w.results)[[1]]<-custom.levels[1] # change rownames of w.results
    dimnames(w.results)[[2]]<-custom.levels[2] # change colnames of w.results
  }
  return(list(kt.results = kt.results,
         w.results = w.results))
}
```



### Test on ages:  


``` r
custom.levels <- names(table(custom.md$agegroup))
alpha.test.age<-test_alpha_diversity(div.indices = div.indices, 
                     all.div = all.div, 
                     comparison = "age",
                     custom.levels = custom.levels,
                     custom.md = custom.md)
```

```
## Joining with `by = join_by(Sample)`
## Joining with `by = join_by(Sample)`
## Joining with `by = join_by(Sample)`
```

``` r
alpha.test.age$kt.results
```

```
##                   sobs      shannon  invsimpson
## statistic 1.179804e+01 10.671111111 6.760000000
## pvalue    5.929321e-04  0.001088217 0.009322376
```

``` r
alpha.test.age$w.results
```

```
## , , sobs
## 
##              agegroup10_16
## agegroup0_10  0.0006702222
## 
## , , shannon
## 
##              agegroup10_16
## agegroup0_10   0.001223547
## 
## , , invsimpson
## 
##              agegroup10_16
## agegroup0_10    0.01026813
```

``` r
# stopifnot(all(alpha.test.age$kt.results[2,]<0.05))
all(alpha.test.age$kt.results[2,]<0.05)
```

```
## [1] TRUE
```



### Test on sexes:


``` r
custom.levels <- names(table(custom.md$sex))
alpha.test.sex<-test_alpha_diversity(div.indices = div.indices, 
                     all.div = all.div, 
                     comparison = "sex",
                     custom.levels = custom.levels,
                     custom.md = custom.md)
```

```
## Joining with `by = join_by(Sample)`
## Joining with `by = join_by(Sample)`
## Joining with `by = join_by(Sample)`
```

``` r
alpha.test.sex$kt.results
```

```
##                sobs   shannon invsimpson
## statistic 0.0569384 2.3120000 0.04355556
## pvalue    0.8114020 0.1283788 0.83468269
```

``` r
alpha.test.sex$w.results
```

```
## , , sobs
## 
##        male  
## female "n.s."
## 
## , , shannon
## 
##        male  
## female "n.s."
## 
## , , invsimpson
## 
##        male  
## female "n.s."
```

``` r
# stopifnot(all(alpha.test.sex$kt.results[2,]<0.05))
all(alpha.test.sex$kt.results[2,]<0.05)
```

```
## [1] FALSE
```

For plots.


``` r
max.values<-all.div%>%
  group_by(metric)%>%
  summarise(max_val=max(value))
```



### Prepare data for plotting.


``` r
custom.md<-custom.md%>%
  mutate(agegroup = factor(agegroup,levels = c("agegroup0_10", 
                                               "agegroup10_16")))%>%
  mutate(sex = case_when(grepl("Female",custom.md$sex,ignore.case = T) ~ "female",
                         grepl("Male",custom.md$sex,ignore.case = T) ~ "male"),
         sex = factor(sex,levels = c("female", "male")))
```



### Plot diversity on sexes.


``` r
pretty.level.names.sex <-
  c("female" = "Females",
    "male" = "Males")
custom.levels.sex <- names(pretty.level.names.sex)
pretty.level.names.sex <-pretty.level.names.sex[which(names(pretty.level.names.sex)%in%custom.levels.sex)]
gg.labs.name.sex <-"Host sex"
gg.title.groups.sex <-"groups"

div.plot.sex<-all.div%>%
  filter(metric%in%plot.metrics)%>%
  left_join(custom.md)%>%
  ggplot(aes(x=factor(sex, levels=custom.levels.sex),
             y=value,
             fill=factor(sex)))+
  geom_boxplot(show.legend = FALSE,outliers = F)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             ncol=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+ 
  labs(fill = gg.labs.name.sex)+
  scale_fill_viridis_d(name = NULL,
                       option = "C", 
                       labels = unname(pretty.level.names.sex),
                       end = 0.9,
                       direction = (-1),
                       alpha = 0.5)+
  scale_color_manual(breaks = unname(pretty.level.names.sex),
                     labels=unname(pretty.level.names.sex))+
  scale_x_discrete(labels=pretty.level.names.sex,
                   limits=custom.levels.sex)+ # rename boxplot labels (x axis)
  # scale_fill_manual(values = custom.fill)+
  # ggtitle(paste("Alpha diversity of the gut microbiota of different",as.character(host.class[host]),gg.title.groups,
  #               gg.title.taxon))+
  alpha.plot.theme + # theme for alpha diversity plots
  theme(legend.position = "none",
        strip.text.x = element_text(size=15))
```

```
## Joining with `by = join_by(Sample)`
```

``` r
div.plot.sex.jitter<-div.plot.sex+
  geom_jitter(width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=3,
              show.legend = FALSE)
```

``` r
print(div.plot.sex.jitter)
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-35-1.png" alt="" width="960" />

``` r
div.plot.sex.dots<-div.plot.sex+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) # add dots
```

``` r
print(div.plot.sex.dots)
```

```
## Bin width defaults to 1/30 of the range of the data. Pick better value with
## `binwidth`.
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-36-1.png" alt="" width="960" />

``` r
# for(image.format in image.formats){
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha-jitter",
#                       paste(plot.metrics,collapse = "-"),
#                       host,"sex",agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=div.plot.sex.jitter,
#          width=10, height=6,units="in",
#          # width = 6000,height = 3000, units = "px",
#          dpi=300,device = image.format)
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha-dots",
#                       paste(plot.metrics,collapse = "-"),
#                       host,"sex",agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=div.plot.sex.dots,
#          width=10, height=6,units="in",
#          # width = 6000,height = 3000, units = "px",
#          dpi=300,device = image.format)
# }
```



### Plot diversity on ages.


``` r
pretty.level.names.age<-names(table(custom.md$old_agegroup))
names(pretty.level.names.age)<-names(table(custom.md$agegroup))
custom.levels<-names(pretty.level.names.age)
gg.labs.name.age<-"Age group"
gg.title.groups.age<-"age groups"

nmr.sample.levels<-custom.md%>%
  filter(class=="NMR")%>%
  # mutate(age=round(time_length(sampling_date-birthday,unit="years")))%>%
  dplyr::select(Sample,age,class)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age," y)"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))

nmr.age.div.plot<-all.div%>%
  filter(metric%in%plot.metrics)%>%
  left_join(nmr.sample.levels)%>%
  ggplot(aes(x=age,y=value,fill=age))+
  geom_point(size=4,stat = "identity",colour = "white",shape=21)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             nrow=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  scale_x_continuous(breaks = round(seq(min(nmr.sample.levels$age), 
                                        max(nmr.sample.levels$age), by = 1),1)) +
  scale_fill_viridis_c(option="D")+
  labs(fill="Age (years)",
       x="Age (years)")+
  theme_bw()+
  # ggtitle(paste0("Alpha diversity of the naked mole-rat gut microbiota across age\n",gg.title.taxon))+
  alpha.plot.theme + # general theme
  theme(axis.title.x = element_text(size=20),
        legend.title = element_text(vjust = 2))
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(nmr.age.div.plot)
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-39-1.png" alt="" width="1056" />

``` r
# for(image.format in image.formats){
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha-per-age",
#                       paste(plot.metrics,collapse = "-"),
#                       "NMR","age",agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=nmr.age.div.plot,
#          width=11, height=8,units="in",
#          # width = 4000,height = 3000,units = "px",
#          dpi=300,device = image.format)
# }
```



## PCA
Convert df into wide format.


``` r
ps.q.df.wide.pca<-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
```

Colnames are OTUs and rownames are sample IDs.


``` r
rownames(ps.q.df.wide.pca)<-ps.q.df.wide.pca$Sample
ps.q.df.wide.pca<-ps.q.df.wide.pca[,-1] # remove sample names
ps.q.df.wide.pca<-as.matrix(ps.q.df.wide.pca) # convert to matrix
# The first crucial step is to transform the compositional data:
# Use a log-ratio transformation, typically the centered log-ratio (clr)
# or isometric log-ratio (ilr) transformation.
# The clr transformation is often preferred for interpretability, 
# while ilr is used for statistical properties.
# We will use robust CLR (RCLR) with natural log
ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,method="rclr",logbase = exp)

#### >if you want to exclude specific samples
if(host=="NMR"){
  ps.q.df.wide.pca<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("M40")),]
}else if(host=="mice"){
  ps.q.df.wide.pca<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("mf_1","MSM343")),]
}
# remove ASVs that are zero because their samples are removed
ps.q.df.wide.pca<-ps.q.df.wide.pca[,which(colSums(ps.q.df.wide.pca)!=0)]
ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,
                                method="rclr",
                                logbase = exp)
####<

# calculate principal components
pca.q<-prcomp(ps.q.df.wide.pca.tfm)
str(pca.q)
```

```
## List of 5
##  $ sdev    : num [1:23] 8.06 6.92 6.52 6.26 5.7 ...
##  $ rotation: num [1:1568, 1:23] -0.1123 -0.1153 -0.0275 -0.0522 -0.0422 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:1568] "775f1d61e616978239ee36c337b0403c" "5a818007a6452d5db9223ef1388e451c" "2a9b76ddc8fc1e5f0e10f1fd40de1c79" "c513b4ef037af47125244161fb1eab50" ...
##   .. ..$ : chr [1:23] "PC1" "PC2" "PC3" "PC4" ...
##  $ center  : Named num [1:1568] 1.9 4.27 3.95 4 2.22 ...
##   ..- attr(*, "names")= chr [1:1568] "775f1d61e616978239ee36c337b0403c" "5a818007a6452d5db9223ef1388e451c" "2a9b76ddc8fc1e5f0e10f1fd40de1c79" "c513b4ef037af47125244161fb1eab50" ...
##  $ scale   : logi FALSE
##  $ x       : num [1:23, 1:23] -9.12 17.6 -0.38 4.2 12.72 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:23] "D27" "2D10" "DCG6" "L122" ...
##   .. ..$ : chr [1:23] "PC1" "PC2" "PC3" "PC4" ...
##  - attr(*, "class")= chr "prcomp"
```

``` r
dim(pca.q$x)
```

```
## [1] 23 23
```

``` r
# reverse the signs
pca.q$rotation<- -1*pca.q$rotation

# display principal components (loadings)
head(pca.q$rotation)
```

```
##                                         PC1          PC2         PC3
## 775f1d61e616978239ee36c337b0403c 0.11229656 -0.017578273 -0.09731180
## 5a818007a6452d5db9223ef1388e451c 0.11532601  0.114262592  0.09972222
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.02749886  0.082882525 -0.06280918
## c513b4ef037af47125244161fb1eab50 0.05222805  0.015134460  0.05147067
## 72aa784b6a3f05f45910fe33902caa81 0.04216571  0.005725382  0.05682973
## 37d6420f6843a99bbb892a3b37178db1 0.04183904  0.161195002 -0.09502338
##                                           PC4          PC5          PC6
## 775f1d61e616978239ee36c337b0403c  0.122695874  0.167475081 -0.032721054
## 5a818007a6452d5db9223ef1388e451c -0.043714745  0.036566765 -0.039876348
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79  0.045960763 -0.042518501 -0.002630799
## c513b4ef037af47125244161fb1eab50  0.006356976  0.064962189  0.004214299
## 72aa784b6a3f05f45910fe33902caa81  0.073374998 -0.171138646  0.016795174
## 37d6420f6843a99bbb892a3b37178db1  0.058140902  0.006065293  0.073430044
##                                          PC7          PC8          PC9
## 775f1d61e616978239ee36c337b0403c -0.11717519 -0.045617847 -0.083021792
## 5a818007a6452d5db9223ef1388e451c  0.04046028  0.017928239  0.030937646
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 -0.11652036 -0.012588960  0.003286941
## c513b4ef037af47125244161fb1eab50  0.01093424 -0.004909519 -0.004681348
## 72aa784b6a3f05f45910fe33902caa81 -0.16446031 -0.144817935 -0.011862911
## 37d6420f6843a99bbb892a3b37178db1 -0.04704658  0.037128538  0.011259495
##                                           PC10         PC11        PC12
## 775f1d61e616978239ee36c337b0403c  0.0383658583  0.020069367  0.04353561
## 5a818007a6452d5db9223ef1388e451c  0.0376631275  0.014802475  0.05937540
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 -0.0037557361 -0.083190555 -0.04558988
## c513b4ef037af47125244161fb1eab50 -0.0292490722  0.026864649 -0.05064640
## 72aa784b6a3f05f45910fe33902caa81  0.0435671311 -0.070474637 -0.07968657
## 37d6420f6843a99bbb892a3b37178db1  0.0004739799 -0.002513979  0.11465579
##                                          PC13        PC14        PC15
## 775f1d61e616978239ee36c337b0403c -0.006839302  0.02007494  0.01162128
## 5a818007a6452d5db9223ef1388e451c  0.008339569 -0.01905280 -0.01149095
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 -0.043099213 -0.08045531  0.06478452
## c513b4ef037af47125244161fb1eab50  0.018698856 -0.05148711 -0.03875526
## 72aa784b6a3f05f45910fe33902caa81  0.013854348 -0.02403359  0.03607152
## 37d6420f6843a99bbb892a3b37178db1 -0.070535234  0.04356465 -0.08276010
##                                         PC16         PC17         PC18
## 775f1d61e616978239ee36c337b0403c  0.01367695 -0.087460717  0.011184639
## 5a818007a6452d5db9223ef1388e451c -0.04734302 -0.011121974 -0.005378269
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 -0.02120797  0.004878748 -0.066324578
## c513b4ef037af47125244161fb1eab50  0.07366780 -0.092287378  0.039717193
## 72aa784b6a3f05f45910fe33902caa81  0.07931191 -0.009267001 -0.018493649
## 37d6420f6843a99bbb892a3b37178db1  0.06217742  0.019878644 -0.037783028
##                                          PC19        PC20        PC21
## 775f1d61e616978239ee36c337b0403c -0.020680912 -0.04138453  0.03381683
## 5a818007a6452d5db9223ef1388e451c  0.007288279 -0.04247859 -0.03036506
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79  0.053202434 -0.02443008  0.03851377
## c513b4ef037af47125244161fb1eab50  0.015137012  0.07618898  0.06793726
## 72aa784b6a3f05f45910fe33902caa81  0.021320854 -0.08810553 -0.02724575
## 37d6420f6843a99bbb892a3b37178db1 -0.001549329 -0.03450310  0.00673757
##                                         PC22        PC23
## 775f1d61e616978239ee36c337b0403c  0.00189574  0.24098421
## 5a818007a6452d5db9223ef1388e451c  0.05407447  0.09451974
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 -0.03761066  0.30037245
## c513b4ef037af47125244161fb1eab50  0.04605305 -0.14573999
## 72aa784b6a3f05f45910fe33902caa81  0.03780556 -0.05277378
## 37d6420f6843a99bbb892a3b37178db1 -0.02267975 -0.66066385
```

``` r
# reverse the signs of the scores
pca.q$x<- -1*pca.q$x

# display the first six scores
head(pca.q$x)
```

```
##              PC1       PC2       PC3       PC4       PC5        PC6        PC7
## D27    9.1165614  3.830851 -1.394304  3.628042  3.085807  4.4753124 -3.0408370
## 2D10 -17.6006016 -9.521624  1.832852  8.022424 -4.917435 13.2895285 -1.9252541
## DCG6   0.3795246  9.586394  3.326733  3.487412 -4.198388  2.9296244  7.3156171
## L122  -4.1960798 17.600548 -8.894850  5.101943 -2.397220  3.1729576  3.1199405
## 2D14 -12.7223992 -4.773917  9.502256  6.132527  3.493594  0.6996189  5.9602129
## GA17   3.0186730 -4.922015 10.932941 -7.632667  1.556425 -1.0139729 -0.6676987
##            PC8        PC9      PC10       PC11      PC12       PC13      PC14
## D27  -4.194607 -4.6789000 -3.003502 -2.5608388  4.448412 -11.058679  7.349089
## 2D10  1.614644 -1.0894454 -1.735064  0.2942293 -3.633317  -7.066023 -1.639834
## DCG6 -8.444869 13.9456527  0.688520 -2.4436653 -6.977537   0.752624  5.659803
## L122  4.813584 -0.1427333 -1.109490  0.8493992 11.212605   3.266524 -5.482350
## 2D14  4.053665 -1.6934375  4.055454  4.6870298  5.569492   2.373997  2.064765
## GA17 -3.310625  8.4489643 -9.109883 -2.9686274  7.736234  -1.179429 -4.432115
##           PC15      PC16      PC17       PC18      PC19       PC20        PC21
## D27  -7.348624  6.315821 -5.075570 -1.4876272 -4.490656 -3.4013034 -0.08831742
## 2D10  6.348807 -3.142479 -1.486999  2.0039660 -1.699199 -1.3699633  1.44599395
## DCG6 -3.276201 -3.984205  2.816253 -1.4714331  0.885116 -2.3610802  1.17287116
## L122  1.217440 -2.151518  2.953707 -1.0055685 -5.890760  0.8976717  2.62012761
## 2D14 -8.674243  6.089183  5.548628 -1.7303252  6.675852  0.6605632 -0.81230231
## GA17  7.323373  7.252182  2.808790 -0.6281696 -0.312759  0.2458381  2.33191937
##            PC22          PC23
## D27   0.3787286 -2.051278e-14
## 2D10  0.5179635 -2.037400e-14
## DCG6  0.3412885 -3.161501e-14
## L122 -0.8600773  6.521151e-15
## 2D14 -1.2614174  4.939083e-15
## GA17  1.0857881  6.805526e-14
```

``` r
#calculate total variance explained by each principal component
pca.q$sdev^2 / sum(pca.q$sdev^2)
```

```
##  [1] 1.066934e-01 7.851020e-02 6.980386e-02 6.435556e-02 5.334621e-02
##  [6] 5.014471e-02 4.781295e-02 4.531818e-02 4.297095e-02 4.176699e-02
## [11] 4.099048e-02 4.050578e-02 4.006243e-02 3.693788e-02 3.566756e-02
## [16] 3.368975e-02 3.358310e-02 3.265505e-02 3.099356e-02 2.964324e-02
## [21] 2.547015e-02 1.907800e-02 3.953098e-31
```

``` r
#calculate total variance explained by each principal component
var_explained = pca.q$sdev^2 / sum(pca.q$sdev^2)

#create scree plot
qplot(seq_along(1:nrow(ps.q.df.wide.pca.tfm)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
```

```
## Warning: `qplot()` was deprecated in ggplot2 3.4.0.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-42-1.png" alt="" width="672" />



### PCA Plot.


``` r
PC1<-pca.q$x[,1]
PC2<-pca.q$x[,2]
perc.var<-round(100*summary(pca.q)$importance[2,1:2],2)

# PC1.df<-as.data.frame(PC1)%>%rownames_to_column("Sample")
PC1.df<-as.data.frame(PC1)%>%mutate(Sample=rownames(ps.q.df.wide.pca))
# PC2.df<-as.data.frame(PC2)%>%rownames_to_column("Sample")
PC2.df<-as.data.frame(PC2)%>%mutate(Sample=rownames(ps.q.df.wide.pca))

pca.plot<-PC1.df%>%
  left_join(PC2.df)%>%
  # left_join(custom.md)
  left_join(custom.md)
```

```
## Joining with `by = join_by(Sample)`
## Joining with `by = join_by(Sample)`
```

``` r
pca.plot.sex<-pca.plot%>%
  ggplot(aes(x = PC1, y = PC2,
             color = sex,
             fill = sex,
             shape = sex))+
  scale_shape_manual(values = 21:22,
                     labels = unname(pretty.level.names.sex))+
  scale_color_viridis_d(labels = unname(pretty.level.names.sex),
                        option = "C",
                        end = 0.9,
                        direction = (-1))+
  scale_fill_viridis_d(name = NULL,
                       option = "C", 
                       labels = unname(pretty.level.names.sex),
                       end = 0.9,
                       direction = (-1))+
  guides(fill = "none",
         color = "none")+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  labs(shape=gg.labs.name.sex)+
  geom_point(size = 5,
             color = "black" 
  ) +
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"))+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  # ggtitle(paste("PCA between different",as.character(host.class[host]),gg.title.groups, "\n",
  #               gg.title.taxon))+
  scatter.plot.theme+  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(legend.title = element_text(vjust = 2))

pca.plot.age<-pca.plot%>%
    ggplot(aes(x = PC1, 
               y = PC2, 
               color = age,
               shape= class,
               fill = age))+
    scale_shape_manual(values = 21)+
    scale_fill_viridis_c(option = "D")+
    scale_color_viridis_c(option = "D")+
    guides(color = "none",
           shape = "none")+
    labs(fill="Age (years)")+
  geom_point(size = 5,
             color = "black" 
  ) +
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"))+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  # ggtitle(paste("PCA between different",as.character(host.class[host]),gg.title.groups, "\n",
  #               gg.title.taxon))+
  scatter.plot.theme+  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(legend.title = element_text(vjust = 2))
```

``` r
print(pca.plot.sex)
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-45-1.png" alt="" width="1056" />

``` r
print(pca.plot.age)
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-46-1.png" alt="" width="1056" />

Label the samples


``` r
pca.sex.labeled<-pca.plot.sex +
  ggrepel::geom_text_repel(aes(label=Sample),
                           show.legend = FALSE,size=7) # add labels to samples

pca.age.labeled<-pca.plot.age+
  # ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
  ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
                           max.overlaps = Inf,
                           show.legend = FALSE,
                           size=7) +# add labels to samples
  theme(legend.title = element_text(vjust = 2))
```

``` r
print(pca.sex.labeled)
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-48-1.png" alt="" width="1056" />

``` r
print(pca.age.labeled)
```

<img src="008-diversity-inside-custom-host_files/figure-html/unnamed-chunk-49-1.png" alt="" width="1056" />

``` r
# for(image.format in image.formats){
#   ggsave(paste0("./images/diversity/pca/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "pca",
#                       host,
#                       comparison,agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=pca.plot,
#          width=11, height=8,units="in",
#          # width = 4000,height = 2500,units = "px",
#          dpi=300,device = image.format)
#   ggsave(paste0("./images/diversity/pca/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "pca-labeled",
#                       host,
#                       comparison,agglom.rank,truncationlvl,
#                       sep = "-"),".",image.format),
#          plot=pca.labeled,
#          width=11, height=8,units="in",
#          # width = 4000,height = 2500,units = "px",units = "px",
#          dpi=300,device = image.format)
# }
```



### Find ASVs that contribute to PCs.


``` r
aload<-abs(pca.q$rotation)
head(sweep(aload,2,colSums(aload),"/"))
```

```
##                                          PC1          PC2         PC3
## 775f1d61e616978239ee36c337b0403c 0.004524987 0.0006723835 0.003554045
## 5a818007a6452d5db9223ef1388e451c 0.004647059 0.0043706390 0.003642079
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.001108066 0.0031703254 0.002293932
## c513b4ef037af47125244161fb1eab50 0.002104528 0.0005789057 0.001879824
## 72aa784b6a3f05f45910fe33902caa81 0.001699066 0.0002190006 0.002075549
## 37d6420f6843a99bbb892a3b37178db1 0.001685903 0.0061658426 0.003470467
##                                           PC4          PC5          PC6
## 775f1d61e616978239ee36c337b0403c 0.0045175530 0.0059800566 1.148978e-03
## 5a818007a6452d5db9223ef1388e451c 0.0016095381 0.0013056947 1.400231e-03
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.0016922344 0.0015182142 9.237876e-05
## c513b4ef037af47125244161fb1eab50 0.0002340582 0.0023196141 1.479823e-04
## 72aa784b6a3f05f45910fe33902caa81 0.0027016022 0.0061108721 5.897514e-04
## 37d6420f6843a99bbb892a3b37178db1 0.0021406963 0.0002165743 2.578447e-03
##                                           PC7          PC8          PC9
## 775f1d61e616978239ee36c337b0403c 0.0041751382 0.0015762670 0.0030054454
## 5a818007a6452d5db9223ef1388e451c 0.0014416642 0.0006194876 0.0011199639
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.0041518057 0.0004349956 0.0001189895
## c513b4ef037af47125244161fb1eab50 0.0003896043 0.0001696422 0.0001694680
## 72aa784b6a3f05f45910fe33902caa81 0.0058599822 0.0050040004 0.0004294455
## 37d6420f6843a99bbb892a3b37178db1 0.0016763442 0.0012829296 0.0004076014
##                                          PC10         PC11        PC12
## 775f1d61e616978239ee36c337b0403c 0.0012987309 7.048972e-04 0.001505151
## 5a818007a6452d5db9223ef1388e451c 0.0012749426 5.199079e-04 0.002052778
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.0001271362 2.921905e-03 0.001576173
## c513b4ef037af47125244161fb1eab50 0.0009901166 9.435681e-04 0.001750991
## 72aa784b6a3f05f45910fe33902caa81 0.0014748003 2.475283e-03 0.002754993
## 37d6420f6843a99bbb892a3b37178db1 0.0000160448 8.829858e-05 0.003963979
##                                          PC13         PC14         PC15
## 775f1d61e616978239ee36c337b0403c 0.0002414194 0.0006850522 0.0003894717
## 5a818007a6452d5db9223ef1388e451c 0.0002943770 0.0006501718 0.0003851040
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.0015213520 0.0027455169 0.0021711677
## c513b4ef037af47125244161fb1eab50 0.0006600478 0.0017569844 0.0012988314
## 72aa784b6a3f05f45910fe33902caa81 0.0004890423 0.0008201400 0.0012088894
## 37d6420f6843a99bbb892a3b37178db1 0.0024898115 0.0014866326 0.0027735954
##                                          PC16         PC17         PC18
## 775f1d61e616978239ee36c337b0403c 0.0004676719 0.0030642453 0.0003901128
## 5a818007a6452d5db9223ef1388e451c 0.0016188550 0.0003896659 0.0001875905
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.0007251887 0.0001709302 0.0023133572
## c513b4ef037af47125244161fb1eab50 0.0025190090 0.0032333507 0.0013853093
## 72aa784b6a3f05f45910fe33902caa81 0.0027120045 0.0003246756 0.0006450462
## 37d6420f6843a99bbb892a3b37178db1 0.0021261049 0.0006964617 0.0013178469
##                                          PC19         PC20         PC21
## 775f1d61e616978239ee36c337b0403c 7.288666e-04 0.0015435505 0.0012587331
## 5a818007a6452d5db9223ef1388e451c 2.568640e-04 0.0015843566 0.0011302512
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 1.875037e-03 0.0009111874 0.0014335632
## c513b4ef037af47125244161fb1eab50 5.334805e-04 0.0028416785 0.0025287667
## 72aa784b6a3f05f45910fe33902caa81 7.514203e-04 0.0032861394 0.0010141437
## 37d6420f6843a99bbb892a3b37178db1 5.460369e-05 0.0012868884 0.0002507864
##                                          PC22        PC23
## 775f1d61e616978239ee36c337b0403c 0.0000902954 0.024757589
## 5a818007a6452d5db9223ef1388e451c 0.0025756037 0.009710516
## 2a9b76ddc8fc1e5f0e10f1fd40de1c79 0.0017914212 0.030858859
## c513b4ef037af47125244161fb1eab50 0.0021935380 0.014972644
## 72aa784b6a3f05f45910fe33902caa81 0.0018007045 0.005421731
## 37d6420f6843a99bbb892a3b37178db1 0.0010802518 0.067873508
```

``` r
colSums(sweep(aload, 2, colSums(aload), "/"))
```

```
##  PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 
##    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
## PC17 PC18 PC19 PC20 PC21 PC22 PC23 
##    1    1    1    1    1    1    1
```

``` r
colSums(aload)
```

```
##       PC1       PC2       PC3       PC4       PC5       PC6       PC7       PC8 
## 24.816989 26.143224 27.380576 27.159808 28.005601 28.478399 28.064985 28.940432 
##       PC9      PC10      PC11      PC12      PC13      PC14      PC15      PC16 
## 27.623789 29.541038 28.471340 28.924416 28.329548 29.304250 29.838561 29.244756 
##      PC17      PC18      PC19      PC20      PC21      PC22      PC23 
## 28.542335 28.670271 28.374071 26.811259 26.865766 20.994872  9.733751
```

``` r
aload[1,1]/colSums(aload)[1]
```

```
##         PC1 
## 0.004524987
```

``` r
pc.df<-as.data.frame(sweep(aload,2,colSums(aload),"/")[,1:2])
lapply(pc.df,max)
```

```
## $PC1
## [1] 0.00967505
## 
## $PC2
## [1] 0.006165843
```

``` r
max.ind<-lapply(pc.df,which.max)
pc.df[max.ind$PC1,]
```

```
##                                         PC1         PC2
## fbcb14cd27216d317b57152441189a35 0.00967505 0.002296274
```

``` r
pc.df[max.ind$PC2,]
```

```
##                                          PC1         PC2
## 37d6420f6843a99bbb892a3b37178db1 0.001685903 0.006165843
```

``` r
pc1.loadings<-pc.df%>%
  rownames_to_column(var="OTU")%>%
  arrange(-PC1)%>%
  head(n=3)%>%
  left_join(ps.q.df.preprocessed)%>%
  distinct()%>%
  right_join(ps.q.agg[,c("Sample","OTU","Species")])%>%
  filter(class=="NMR")%>%
  arrange(-PC1,-Abundance)
```

```
## Joining with `by = join_by(OTU)`
## Joining with `by = join_by(OTU, Sample)`
```

``` r
pc2.loadings<-pc.df%>%
  rownames_to_column(var="OTU")%>%
  arrange(-PC2)%>%
  head(n=3)%>%
  left_join(ps.q.df.preprocessed)%>%
  distinct()%>%
  right_join(ps.q.agg[,c("Sample","OTU","Species")])%>%
  filter(class=="NMR")%>%
  arrange(-PC2,-Abundance)
```

```
## Joining with `by = join_by(OTU)`
## Joining with `by = join_by(OTU, Sample)`
```

``` r
rbind(pc1.loadings,pc2.loadings)%>%
  dplyr::select(Sample,OTU,Species,Abundance)%>%
  group_by(Species,OTU)%>%
  mutate(ASV= paste(Species,cur_group_id()))%>%
  ungroup%>%
  mutate(OTU=ASV,
         OTU=gsub(" ","_",OTU))%>%
  head
```

```
## # A tibble: 6 × 5
##   Sample OTU                 Species           Abundance ASV                
##   <chr>  <chr>               <chr>                 <dbl> <chr>              
## 1 2D10   Allobaculum_Genus_1 Allobaculum Genus      3558 Allobaculum Genus 1
## 2 2D14   Allobaculum_Genus_1 Allobaculum Genus      2991 Allobaculum Genus 1
## 3 GA5    Allobaculum_Genus_1 Allobaculum Genus      2158 Allobaculum Genus 1
## 4 O15    Allobaculum_Genus_1 Allobaculum Genus      2109 Allobaculum Genus 1
## 5 G14    Allobaculum_Genus_1 Allobaculum Genus      1817 Allobaculum Genus 1
## 6 CG33   Allobaculum_Genus_1 Allobaculum Genus      1621 Allobaculum Genus 1
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
##  [1] ggrepel_0.9.6   lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
##  [5] dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1    
##  [9] tibble_3.2.1    ggplot2_4.0.0   tidyverse_2.0.0 vegan_2.6-4    
## [13] lattice_0.22-6  permute_0.9-8  
## 
## loaded via a namespace (and not attached):
##  [1] utf8_1.2.4         sass_0.4.10        generics_0.1.4     xml2_1.3.8        
##  [5] stringi_1.8.4      hms_1.1.4          digest_0.6.37      magrittr_2.0.3    
##  [9] timechange_0.3.0   evaluate_1.0.5     grid_4.4.3         RColorBrewer_1.1-3
## [13] bookdown_0.46      fastmap_1.2.0      jsonlite_2.0.0     Matrix_1.7-4      
## [17] ggtext_0.1.2       mgcv_1.9-1         viridisLite_0.4.3  scales_1.4.0      
## [21] jquerylib_0.1.4    cli_3.6.4          rlang_1.1.5        litedown_0.9      
## [25] commonmark_2.0.0   splines_4.4.3      withr_3.0.2        cachem_1.1.0      
## [29] yaml_2.3.12        otel_0.2.0         tools_4.4.3        parallel_4.4.3    
## [33] tzdb_0.5.0         vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.5   
## [37] MASS_7.3-65        cluster_2.1.8.1    pkgconfig_2.0.3    bslib_0.10.0      
## [41] pillar_1.11.1      gtable_0.3.6       Rcpp_1.0.14        glue_1.8.0        
## [45] xfun_0.56          tidyselect_1.2.1   rstudioapi_0.18.0  knitr_1.51        
## [49] farver_2.1.2       htmltools_0.5.8.1  nlme_3.1-167       labeling_0.4.3    
## [53] rmarkdown_2.30     compiler_4.4.3     S7_0.2.0           markdown_2.0      
## [57] gridtext_0.1.5
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used (Mb)
## Ncells 2908681 155.4    4941454 264.0  4941454  264
## Vcells 5086098  38.9   10146329  77.5  8388608   64
```

