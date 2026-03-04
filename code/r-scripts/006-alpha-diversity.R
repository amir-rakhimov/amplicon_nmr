#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 006-alpha-diversity.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# Alpha diversity analysis of QIIME2 output ####
#' # Alpha diversity analysis of QIIME2 output
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script analyses alpha diversity among hosts. We will compare 
#' hosts using three metrics: the number of observed genera, 
#' Shannon diversity index, and inverse Simpson index.  
#' 
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries and scripts.
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan","Polychrome"))
# library(phyloseq)
library(tidyverse)
library(vegan)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Name of the folder with QIIME2 output:
authorname<-"pooled"
#' Specify paths and image formats:
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-file.path("./output/rtables",authorname)
image.formats<-c("png","tiff")
#' Truncation level that we chose in QIIME2:
truncationlvl<-"234"
#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"Genus"
#' Single reads or paired reads (decided in QIIME2):
read.end.type<-"single"
#' Import abundance table as an rds file (NOT rarefied): 
ps.q.agg.date_time<-"20260211_17_01_10"
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#'
#' Import the rarefied abundance table (rds  file):
ps.q.df.preprocessed.date_time<-"20260211_17_14_19"
# 20260211_17_14_19 rarefied table for all hosts, genus level
# 20260211_17_14_20 rarefied table file for NMR, genus level
# 20260211_17_14_21 rarefied table file for NMR, ASV level

#' Import metadata:
custom.md<-readRDS(file.path(rdafiles.directory,"custom.md.rds"))
#' Set "pretty" labels
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

#' Facet labels
metric.labs <- c('sobs'= case_when(agglom.rank == "OTU" ~ 
                                     "Richness \n(Observed species)",
                                agglom.rank == "Genus" ~ 
                                  "Richness \n(Observed genera)") ,
              'shannon' = "Shannon",
              'invsimpson' = "Inverse \nSimpson")
#' Metrics to plot
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
  
#' Remove rows with 0 Abundance, if there's any.
ps.q.agg <- ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)

#' Import rarefied dataframe and select Sample, Abundance, 
#' class, and agglom.rank columns.
ps.q.df.preprocessed<-readRDS(
  file.path(rdafiles.directory,paste0(
    paste(
      ps.q.df.preprocessed.date_time,
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".rds")))%>%
  filter(class%in%custom.levels, Abundance!=0)%>%
  dplyr::select(Sample,Abundance,class,all_of(agglom.rank),)

#+ echo=FALSE
## 3. Alpha diversity analysis. ####
#' 
#' ## Alpha diversity analysis.
#+ echo=FALSE
### 3.1 Compute alpha diversity metrics. ####
#' 
#' ### Compute alpha diversity metrics.
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

#+ echo=FALSE
### 3.2 Alpha diversity tests. ####
#' 
#' ### Alpha diversity tests. ####
#' The Kruskal-Wallis rank sum test analyzes differences in each 
#' alpha diversity metric across animal hosts.
kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
colnames(kt.results)<-div.indices
rownames(kt.results)<-c("statistic", "pvalue")
#' All unique pairwise combinations:
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
w.results
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

#+ echo=FALSE
## 4. Plot alpha diversity metrics. ####
#'
#' ## Plot alpha diversity metrics. ####
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
#+ fig.height=5, fig.width=11
print(div.plot.jitter + 
        ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14))
)

div.plot.dots<-div.plot+
  geom_dotplot(aes(colour=class),
               fill="white",
               color="black",
               binaxis='y',
               stackdir='center', 
               dotsize=0.4) # add dots
#+ fig.height=5, fig.width=11
print(div.plot.dots + 
        ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14))
)
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

#+ echo=FALSE
### 4.1 Add significance stars. ####
#' 
#' ### Add significance stars. 
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
#+ fig.width=11, fig.height=8
print(newplot)

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
rm(list = ls(all=TRUE))
gc()