#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
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
library(phyloseq)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
message("The analysis focus is: ", active.analysis)

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
source(file.path(util.functions.r,"get_average_rarefied_diversity.R"))
#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"Genus"
ps.q<-readRDS(file= ps.q.filtered.fname)
custom.md <- readRDS(custom.md.path)
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

#' Import rarefied dataframe 
ps.q.agg.asv.rare <-readRDS(file.path(community.composition.rdafiles,paste0(
  paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
)
#+ echo=FALSE
## 3. Alpha diversity analysis. ####
#' 
#' ## Alpha diversity analysis.
#+ echo=FALSE
### 3.1 Compute alpha diversity metrics. ####
#' 
#' ### Compute alpha diversity metrics.
all.div.fname <-file.path(diversity.rdafiles,paste0(
  paste("ps.q.rare.diversity.average",agglom.rank,paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
if(!file.exists(all.div.fname)){
  all.div<-get_average_rarefied_diversity(ps.q.agg.asv.rare,agglom.rank)%>%
    left_join(custom.md, by= "Sample")
  saveRDS(all.div,
          file = all.div.fname)
}else{
  all.div<- readRDS(all.div.fname)
}

#+ echo=FALSE
### 3.2 Alpha diversity tests. ####
#' 
#' ### Alpha diversity tests. ####
#' Function that runs Kruskal-Wallis test and pairwise Wilcoxon tests:
test_alpha_diversity<-function(div.indices, 
                               div.df, 
                               comparison,
                               comparison_colname,
                               custom.levels,
                               metadata){
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
    metric.ds <- div.df%>%
      left_join(metadata)%>%
      filter(metric==div.metric)%>%
      distinct()
    # perform kruskal-wallis test
    if(comparison=="age"){
      kt<-kruskal.test(value ~ agegroup,data = metric.ds)
      
    }else if (comparison == "sex"){
      kt<-kruskal.test(value ~ sex,data = metric.ds)
      
    }else if(comparison == "strain" | comparison == "host" ){
      kt<-kruskal.test(value ~ class, data = metric.ds)
      
    }
    
    kt.results["statistic",div.metric]<- kt$statistic
    kt.results["pvalue",div.metric]<- kt$p.value
    # low pvalue indicates statistical difference between some groups
    # But which groups are different from others?
    # Perform pairwise wilcoxon test
    
    # NB: You must not run pairwise wilcoxon test unless kruskal wallis results 
    # are significant
    if(kt$p.value<0.05){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds[[comparison_colname]],
                                   p.adjust.method = "BH",
                                   exact=FALSE)
      w.results[,,div.metric]<-w.test$p.value
      dimnames(w.results)[[1]]<-dimnames(w.test$p.value)[[1]] # change rownames of w.results
      dimnames(w.results)[[2]]<-dimnames(w.test$p.value)[[2]] # change colnames of w.results
    }else{
      w.results[,,div.metric]<-matrix(data = "n.s.",
                                      nrow = length(table(combinations[1,])), # nrows,
                                      ncol = length(table(combinations[2,])))
      # nrow = nrow(w.test$p.value),
      # ncol = ncol(w.test$p.value))
      dimnames(w.results)[[1]]<-custom.levels[1] # change rownames of w.results
      dimnames(w.results)[[2]]<-custom.levels[2] # change colnames of w.results
    }
    

  }
  return(list(kt.results = kt.results,
              w.results = w.results))
}

alpha.test.host <- test_alpha_diversity(div.indices = div.indices, 
                                     div.df = all.div, 
                                     comparison = "host",
                                     comparison_colname = "class",
                                     custom.levels = custom.levels,
                                     metadata = custom.md)
alpha.test.host$kt.results
alpha.test.host$w.results
# stopifnot(all(alpha.test.age$kt.results[2,]<0.05))
#' Are all p-values less than 0.05?
all(alpha.test.host$kt.results[2,]<0.05)

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
  scale_fill_viridis_d(option="C")

saveRDS(div.plot,
        file = file.path(manuscript.figures,"div.plot.rds"))

div.plot <- div.plot +
  stat_summary(fun=median, geom="point", shape=23, size=3, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  labs(fill = "Host")+
  # ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
  project_theme + # general theme
  alpha.plot.theme +
  theme(axis.text.x = ggtext::element_markdown(size=12),
    axis.text.y = element_text(size=15),
    axis.title = element_text(size = 15),
    strip.text.x = element_text(size=15),
    plot.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18),
    legend.position = "right")

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
for(image.format in image.formats){
  ggsave(paste0(paste("alpha-jitter",
                      paste(plot.metrics,collapse = "-"),
                      "all_hosts",
                      agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.jitter,
         path = file.path(diversity.figures,"alpha"),
         width=11, height=5,units="in",
         dpi=300,device = image.format)
  ggsave(paste0(paste("alpha-dots",
                      paste(plot.metrics,collapse = "-"),
                      "all_hosts",
                      agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.dots,
         path = file.path(diversity.figures,"alpha"),
         width=11, height=5,units="in",
         dpi=300,device = image.format)
}


#+ echo=FALSE
### 4.1 Add significance stars. ####
#' 
#' ### Add significance stars. 
coord.combs<-combn(seq_along(custom.levels),2)

# First, find significant results (stars)
stars.list<-matrix(NA,nrow=length(div.indices),ncol = ncol(coord.combs))
for (k in 1:dim(alpha.test.host$w.results)[3]) { # loop over each sub array
  for (i in 1:dim(alpha.test.host$w.results)[1]) { # then each row
    for (j in 1:dim(alpha.test.host$w.results)[2]) { # then each column
      if (!is.na(alpha.test.host$w.results[i,j,k])) { # extract non-NA values
        ith.row<-rownames(alpha.test.host$w.results[,,k])[i] 
        jth.col<-colnames(alpha.test.host$w.results[,,k])[j]
        w.val<- alpha.test.host$w.results[i,j,k]
        
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



#+ echo=FALSE
# Alpha and beta diversity in NMR samples ####
#' # Alpha and beta diversity in NMR samples 
#+ echo=FALSE
# 5. Analyse naked mole-rat data #####
agglom.rank <- "OTU"
#' The taxonomic rank that was used for agglomeration:
ps.q.nmr <- readRDS(file= ps.q.filtered.fname)%>%
  subset_samples(class == "NMR")
#' Load metadata with age information:
custom.md<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

#' Facet labels
metric.labs <- c('sobs'= case_when(agglom.rank == "OTU" ~ 
                                     "Richness \n(Observed species)",
                                   agglom.rank == "Genus" ~ 
                                     "Richness \n(Observed genera)") ,
                 'shannon' = "Shannon",
                 'invsimpson' = "Inverse \nSimpson")
#' Metrics to plot
gg.title.taxon <- ifelse(agglom.rank=="OTU","(ASV level)",
                         paste0("(",agglom.rank," level)"))
host.labels <- c("NMR" = "*Heterocephalus glaber*")

#' Import rarefied dataframe 
ps.q.agg.asv.nmr.rare <-readRDS(file.path(community.composition.rdafiles,paste0(
  paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
)
#+ echo=FALSE
## 6. Alpha diversity analysis. ####
#' 
#' ## Alpha diversity analysis.
#+ echo=FALSE
### 6.1 Compute alpha diversity metrics. ####
#' 
#' ### Compute alpha diversity metrics.
all.div.nmr.fname <-file.path(diversity.rdafiles,paste0(
  paste("ps.q.rare.diversity.average",agglom.rank,paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
if(!file.exists(all.div.nmr.fname)){
  all.div.nmr <- get_average_rarefied_diversity(ps.q.agg.asv.nmr.rare,
                                                agglom.rank)%>%
    left_join(custom.md, by= "Sample")
  saveRDS(all.div.nmr,
          file = all.div.nmr.fname)
}else{
  all.div.nmr <- readRDS(all.div.nmr.fname)
}

#+ echo=FALSE
### 6.2 Alpha diversity tests. ####
#'
#' ### Alpha diversity tests. 

#+ echo=FALSE
#### 6.2.1 Test on ages: ####
#' 
#' #### Test on ages:  
custom.levels <- names(table(custom.md$agegroup))
alpha.test.age<-test_alpha_diversity(div.indices = div.indices, 
                                     div.df = all.div.nmr, 
                                     comparison = "age",
                                     comparison_colname = "agegroup",
                                     custom.levels = custom.levels,
                                     metadata = custom.md)
alpha.test.age$kt.results
alpha.test.age$w.results
# stopifnot(all(alpha.test.age$kt.results[2,]<0.05))
all(alpha.test.age$kt.results[2,]<0.05)
#+ echo=FALSE
#### 6.2.2 Test on sexes: ####
#' 
#' #### Test on sexes:
custom.levels <- names(table(custom.md$sex))
alpha.test.sex<-test_alpha_diversity(div.indices = div.indices, 
                                     div.df = all.div.nmr, 
                                     comparison = "sex",
                                     comparison_colname = "sex",
                                     custom.levels = custom.levels,
                                     metadata = custom.md)
alpha.test.sex$kt.results
alpha.test.sex$w.results
# stopifnot(all(alpha.test.sex$kt.results[2,]<0.05))
all(alpha.test.sex$kt.results[2,]<0.05)
#' For plots.
max.values<-all.div.nmr%>%
  group_by(metric)%>%
  summarise(max_val=max(value))
#+ echo=FALSE
### 6.3 Prepare data for plotting. ####
#'
#' ### Prepare data for plotting.
custom.md <- custom.md%>%
  mutate(agegroup = factor(agegroup,levels = c("agegroup0_10", 
                                               "agegroup10_16")))%>%
  mutate(sex = case_when(grepl("Female",custom.md$sex,ignore.case = T) ~ "female",
                         grepl("Male",custom.md$sex,ignore.case = T) ~ "male"),
         sex = factor(sex,levels = c("female", "male")))

#+ echo=FALSE
### 6.4 Plot diversity on ages. ####
#' 
#' ### Plot diversity on ages.
pretty.level.names.age<-names(table(custom.md$agegroup))
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

nmr.age.div.plot<-all.div.nmr%>%
  filter(metric%in%plot.metrics)%>%
  left_join(nmr.sample.levels)%>%
  ggplot(aes(x=age,y=value,fill=age))+
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
  project_theme 

saveRDS(nmr.age.div.plot,
        file = file.path(manuscript.figures,"nmr.age.div.plot.rds"))

nmr.age.div.plot <- nmr.age.div.plot +
  geom_point(size=4,stat = "identity",colour = "white",shape=21)+
  alpha.plot.theme + # general theme
  theme(axis.title.x = element_text(size=20),
        legend.title = element_text(vjust = 2))
#+ fig.width=11, fig.height=8
print(nmr.age.div.plot)
for(image.format in image.formats){
  ggsave(paste0(paste("alpha-per-age",
                      paste(plot.metrics,collapse = "-"),
                      "NMR","age",agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot = nmr.age.div.plot,
         path = file.path(diversity.figures,"alpha"),
         width=11, height=8,units="in",
         dpi=300,device = image.format)
}


#+ echo=FALSE
### 6.5 Plot diversity on sexes. ####
#' 
#' ### Plot diversity on sexes.
pretty.level.names.sex <-
  c("female" = "Females",
    "male" = "Males")
custom.levels.sex <- names(pretty.level.names.sex)
pretty.level.names.sex <-pretty.level.names.sex[which(names(pretty.level.names.sex)%in%custom.levels.sex)]
gg.labs.name.sex <-"Host sex"
gg.title.groups.sex <-"groups"

div.plot.sex <- all.div.nmr %>%
  filter(metric %in% plot.metrics)%>%
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
  project_theme +
  alpha.plot.theme + # theme for alpha diversity plots
  theme(legend.position = "none",
        strip.text.x = element_text(size=15))

div.plot.sex.jitter<-div.plot.sex+
  geom_jitter(width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=3,
              show.legend = FALSE)
#+ fig.width=10, fig.height=6
print(div.plot.sex.jitter)

div.plot.sex.dots<-div.plot.sex+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6) # add dots
#+ fig.width=10, fig.height=6
print(div.plot.sex.dots)

for(image.format in image.formats){
  ggsave(paste0(paste("alpha-jitter",
                      paste(plot.metrics,collapse = "-"),
                      "NMR","sex",agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.sex.jitter,
         path = file.path(diversity.figures,"alpha"),
         width=10, height=6,units="in",
         dpi=300,device = image.format)
  ggsave(paste0(paste("alpha-dots",
                      paste(plot.metrics,collapse = "-"),
                      "NMR","sex",agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.sex.dots,
         path = file.path(diversity.figures,"alpha"),
         width=10, height=6,units="in",
         # width = 6000,height = 3000, units = "px",
         dpi=300,device = image.format)
}


sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()