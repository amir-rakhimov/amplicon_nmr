#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---

#' ```{r, setup 008-diversity-inside-custom-host.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# Alpha and beta diversity in NMR samples ####
#' # Alpha and beta diversity in NMR samples 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script performs alpha diversity analysis and PCA on NMR samples

#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
library(vegan)
library(tidyverse)
library(ggrepel)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"OTU"
rare.status<-"rare"
filter.status<-"nonfiltered"
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
# comparison<- cmdargs[2] #"age" or "sex" 
# host <- cmdargs[3] # NMR
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
# comparison<- params$comparison #"age" or "sex"
# host <-"NMR"# params$host # NMR
print(paste("The analysis focus is:", active.analysis, comparison, host))
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

ps.q.agg<-readRDS(file=file.path(
  community.composition.rdafiles,
  paste("phyloseq-qiime",authorname,agglom.rank,read.end.type,
        truncationlvl,"table.rds",sep = "-")))%>%
  filter(class == "NMR",Abundance!=0)
custom.md <- readRDS(custom.md.path)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
gg.title.taxon <- ifelse(agglom.rank=="OTU","(ASV level)",
                       paste0("(",agglom.rank," level)"))
host.labels <- c("NMR" = "*Heterocephalus glaber*")
#' Load the rarefied dataframe:
ps.q.df.preprocessed<-readRDS(
  file.path(community.composition.rdafiles,
            paste0(
              paste("ps.q.df.rare-nonfiltered",agglom.rank,
                    paste(custom.levels,collapse = '-'),sep = "-"),
              ".rds")))%>%
  filter(class == "NMR",Abundance!=0)
#' Load metadata with age information:
custom.md<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
#' Labels for plots:
metric.labs=c('sobs'= case_when(agglom.rank == "OTU" ~ "Richness (Observed species)",
                                agglom.rank == "Genus" ~ "Richness (Observed genera)") ,
              'shannon' = "Shannon",
              'invsimpson' = "Inverse Simpson")
plot.metrics<-c("sobs","shannon",
                "invsimpson") # metrics to plot
div.indices<-c("sobs","shannon",
               "invsimpson")

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
          # simpson=diversity(Abundance, index="simpson"),
          invsimpson=diversity(Abundance, index="invsimpson"), # inverse simpson
          tot=sum(Abundance))%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()

#+ echo=FALSE
### 3.2 Alpha diversity tests. ####
#'
#' ### Alpha diversity tests. 
#' Function that runs Kruskal-Wallis test and pairwise Wilcoxon tests:
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

#+ echo=FALSE
### 3.3 Test on ages: ####
#' 
#' ### Test on ages:  
custom.levels <- names(table(custom.md$agegroup))
alpha.test.age<-test_alpha_diversity(div.indices = div.indices, 
                     all.div = all.div, 
                     comparison = "age",
                     custom.levels = custom.levels,
                     custom.md = custom.md)
alpha.test.age$kt.results
alpha.test.age$w.results
# stopifnot(all(alpha.test.age$kt.results[2,]<0.05))
all(alpha.test.age$kt.results[2,]<0.05)
#+ echo=FALSE
### 3.4 Test on sexes: ####
#' 
#' ### Test on sexes:
custom.levels <- names(table(custom.md$sex))
alpha.test.sex<-test_alpha_diversity(div.indices = div.indices, 
                     all.div = all.div, 
                     comparison = "sex",
                     custom.levels = custom.levels,
                     custom.md = custom.md)
alpha.test.sex$kt.results
alpha.test.sex$w.results
# stopifnot(all(alpha.test.sex$kt.results[2,]<0.05))
all(alpha.test.sex$kt.results[2,]<0.05)
#' For plots.
max.values<-all.div%>%
  group_by(metric)%>%
  summarise(max_val=max(value))

#+ echo=FALSE
### 3.5 Prepare data for plotting. ####
#'
#' ### Prepare data for plotting.
custom.md<-custom.md%>%
  mutate(agegroup = factor(agegroup,levels = c("agegroup0_10", 
                                               "agegroup10_16")))%>%
  mutate(sex = case_when(grepl("Female",custom.md$sex,ignore.case = T) ~ "female",
                         grepl("Male",custom.md$sex,ignore.case = T) ~ "male"),
         sex = factor(sex,levels = c("female", "male")))


#+ echo=FALSE
### 3.6 Plot diversity on sexes. ####
#' 
#' ### Plot diversity on sexes.
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

#+ echo=FALSE
### 3.7 Plot diversity on ages. ####
#' 
#' ### Plot diversity on ages.
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
  project_theme +
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
## 4. PCA ####
#' 
#' ## PCA
#' Convert df into wide format.
ps.q.df.wide.pca<-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

#' Colnames are OTUs and rownames are sample IDs.
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
ps.q.df.wide.pca<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("M40")),]
# if(host=="NMR"){
#   ps.q.df.wide.pca<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("M40")),]
# }else if(host=="mice"){
#   ps.q.df.wide.pca<-ps.q.df.wide.pca[-which(rownames(ps.q.df.wide.pca)%in%c("mf_1","MSM343")),]
# }
# remove ASVs that are zero because their samples are removed
ps.q.df.wide.pca<-ps.q.df.wide.pca[,which(colSums(ps.q.df.wide.pca)!=0)]
ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,
                                method="rclr",
                                logbase = exp)
####<

#' Calculate principal components:
pca.q<-prcomp(ps.q.df.wide.pca.tfm)
str(pca.q)
dim(pca.q$x)

#' Reverse the signs:
pca.q$rotation<- -1*pca.q$rotation

#' Display principal components (loadings):
head(pca.q$rotation)

#' Reverse the signs of the scores:
pca.q$x<- -1*pca.q$x

#' Display the first six scores:
head(pca.q$x)

#' Calculate total variance explained by each principal component:
pca.q$sdev^2 / sum(pca.q$sdev^2)

#' Calculate total variance explained by each principal component:
var_explained = pca.q$sdev^2 / sum(pca.q$sdev^2)

#' Create scree plot:
qplot(seq_along(1:nrow(ps.q.df.wide.pca.tfm)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

#+ echo=FALSE
### 4.1 PCA Plot. ####
#' 
#' ### PCA Plot.
PC1<-pca.q$x[,1]
PC2<-pca.q$x[,2]
perc.var<-round(100*summary(pca.q)$importance[2,1:2],2)

# PC1.df<-as.data.frame(PC1)%>%rownames_to_column("Sample")
PC1.df<-as.data.frame(PC1)%>%mutate(Sample=rownames(ps.q.df.wide.pca))
# PC2.df<-as.data.frame(PC2)%>%rownames_to_column("Sample")
PC2.df<-as.data.frame(PC2)%>%mutate(Sample=rownames(ps.q.df.wide.pca))

pca.plot<-PC1.df%>%
  left_join(PC2.df)%>%
  left_join(biosample.md%>%filter(sequencing_type=="Naked mole-rat 16S rRNA gene sequencing"),
            by = join_by("Sample" == "host_subject_id"))%>%
  mutate(first_collection_date = as.Date(first_collection_date),
         collection_year = year(first_collection_date))
  # left_join(custom.md)

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
  project_theme +
  pca.plot.theme +  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(legend.text = ggtext::element_markdown(size = 20),
        legend.title = element_text(vjust = 2))

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
           shape = "none"
           )+
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
  project_theme +
  pca.plot.theme +  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(legend.title = element_text(vjust = 2))

pca.plot.batch<-pca.plot%>%
  ggplot(aes(x = PC1, 
             y = PC2, 
             color = collection_year,
             fill = collection_year))+
  scale_fill_viridis_c(option = "D")+
  scale_color_viridis_c(option = "D")+
  guides(color = "none"
  )+
  labs(fill="Collection year")+
  geom_point(size = 5  ) +
  labs(x=paste0("PC1 (", perc.var[1], "%)"),
       y=paste0("PC2 (", perc.var[2], "%)"))+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")


#+ fig.width=11, fig.height=8
print(pca.plot.sex)
#+ fig.width=11, fig.height=8
print(pca.plot.age)

#' Label the samples
pca.sex.labeled<-pca.plot.sex +
  ggrepel::geom_text_repel(aes(label=Sample),
                           show.legend = FALSE,size=7) # add labels to samples

pca.age.labeled<-pca.plot.age+
  ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
                           max.overlaps = Inf,
                           show.legend = FALSE,
                           size=7) +# add labels to samples
  theme(legend.title = element_text(vjust = 2))

#+ fig.width=11, fig.height=8
print(pca.sex.labeled)
#+ fig.width=11, fig.height=8
print(pca.age.labeled)

pca.plots<-list("age" = c(pca.plot.age, pca.age.labeled),
                "sex" = c(pca.plot.sex, pca.sex.labeled)
                )
for (pca_comparison in names(pca.plots)){
  for(image.format in image.formats){
    ggsave(paste0(paste("pca",
                        "NMR",
                        pca_comparison,agglom.rank,truncationlvl,
                        sep = "-"),".",image.format),
           plot = pca.plots[[pca_comparison]][1],
           path = file.path(diversity.figures,"pca"),
           width=11, height=8,units="in",
           dpi=300,device = image.format)
    ggsave(paste0(paste("pca-labeled",
                        "NMR",
                        pca_comparison,agglom.rank,truncationlvl,
                        sep = "-"),".",image.format),
           plot = pca.plots[[pca_comparison]][2],
           path = file.path(diversity.figures,"pca"),
           width=11, height=8,units="in",
           dpi=300,device = image.format)
  }
}

#+ echo=FALSE
### 4.2 Find ASVs that contribute to PCs. ####
#' 
#' ### Find ASVs that contribute to PCs.
aload<-abs(pca.q$rotation)
head(sweep(aload,2,colSums(aload),"/"))
colSums(sweep(aload, 2, colSums(aload), "/"))
colSums(aload)
aload[1,1]/colSums(aload)[1]

pc.df<-as.data.frame(sweep(aload,2,colSums(aload),"/")[,1:2])
lapply(pc.df,max)
max.ind<-lapply(pc.df,which.max)
pc.df[max.ind$PC1,]
pc.df[max.ind$PC2,]

pc1.loadings<-pc.df%>%
  rownames_to_column(var="OTU")%>%
  arrange(-PC1)%>%
  head(n=3)%>%
  left_join(ps.q.df.preprocessed)%>%
  distinct()%>%
  right_join(ps.q.agg[,c("Sample","OTU","Species")])%>%
  filter(class=="NMR")%>%
  arrange(-PC1,-Abundance)

pc2.loadings<-pc.df%>%
  rownames_to_column(var="OTU")%>%
  arrange(-PC2)%>%
  head(n=3)%>%
  left_join(ps.q.df.preprocessed)%>%
  distinct()%>%
  right_join(ps.q.agg[,c("Sample","OTU","Species")])%>%
  filter(class=="NMR")%>%
  arrange(-PC2,-Abundance)

rbind(pc1.loadings,pc2.loadings)%>%
  dplyr::select(Sample,OTU,Species,Abundance)%>%
  group_by(Species,OTU)%>%
  mutate(ASV= paste(Species,cur_group_id()))%>%
  ungroup%>%
  mutate(OTU=ASV,
         OTU=gsub(" ","_",OTU))%>%
  distinct(OTU,.keep_all = )%>%
  head


sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), "active.analysis"))
gc()