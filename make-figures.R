#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---
#' 
#' ```{r, setup make-figures.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("make-figures.R", knit = FALSE)
#' #file.rename("make-figures.Rmd", "markdown/make-figures.Rmd")
#' #rmarkdown::render('./markdown/make-figures.Rmd', 'html_document',
#' # knit_root_dir="/home/rakhimov/projects/amplicon_nmr/")
#' ```
#' 
#+ echo=FALSE
# Generating figures ####
#' # Generating figures
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script generates all figures from 16S rRNA gene sequencing data:
#' 
#' 1) Alpha diversity on hosts (genus level)
#' 
#' 2) PCA plot on hosts (genus level)
#' 
#' 3) PCA plot on NMR ages (ASV level)
#' 
#' 4) PCA plot on NMR sexes (ASV level)
#' 
#' 5) Alpha diversity on NMR ages (ASV level)
#' 
#' 6) Barplot of the most abundant ASVs in NMR
#'
#' 7) Comparison of our NMR data with NMR data from Debebe et al.
#' 
#' We will use the data from 001-phyloseq-qiime2.R script (ps.q.agg
#' agglomerated tables at phylum, family, genus, and OTU level).

#+ echo=FALSE
## 1. Load necessary libraries and scripts. ####
#'
#' ## Load necessary libraries and scripts.
# install.packages(c("tidyverse","vegan","ggtext","Polychrome","ggrepel"))
library(tidyverse)
library(vegan)
library(Polychrome)
library(ggtext)
library(ggrepel)
library(cowplot)

#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Use only the taxa that are present in the workspace
#' (custom.md is metadata from the rdafile).
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

plot.metrics<-c("sobs","shannon","invsimpson")
div.indices<-c("sobs","shannon", "invsimpson")

#+ echo=FALSE
## 3. Plots. ####
#'
#' ## Plots.
#+ echo=FALSE
### Alpha diversity on hosts ####
#' ### Alpha diversity on hosts ####
div.plot.jitter <- readRDS(file = file.path(manuscript.figures,"div.plot.rds"))+
  stat_summary(fun=median, geom="point", shape=23, size=2, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  labs(fill = "Host",
       y="",
       x="")+
  # ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
  project_theme + 
  manuscript.alpha.plot.theme +
  theme(
    # axis.title = element_text(size = 10),
    # plot.title = element_text(size = 10),
    # legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    legend.box.spacing = unit(0.01,units="cm"),
    legend.position = "right",
    legend.text = ggtext::element_markdown(size=10))+  
  geom_jitter(aes(colour=class_letter),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=0.7,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
  guides(color = "none")

#+ echo=FALSE
### PCA on hosts ####
#' ### PCA on hosts ####
pca.hosts.plot <- readRDS(file.path(manuscript.figures,"pca.hosts.plot.rds")) + 
  geom_point(aes(shape = class_name),
             size = 1.5,
             color = "black")+ 
  theme_bw() +
  project_theme +
  manuscript.pca.plot.theme +
  theme(axis.text.x = element_text(angle=0,size=13),
        axis.text.y = element_text(size=13),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.key.spacing = unit(0.01,"cm"),
        legend.box.spacing = unit(0.01,units="cm"),
        plot.caption = ggtext::element_markdown(hjust = 0,
                                                size=12),
        legend.position = "bottom")+
  labs(shape ="")

#+ echo=FALSE
### PCA on NMR ages ####
#' ### PCA on NMR ages ####
pca.nmr.plot.age.labeled <- readRDS(file.path(manuscript.figures,"pca.nmr.plot.age.rds")) + 
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  geom_point(size = 2,
           color = "black" ) +
  project_theme +
  ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
                           max.overlaps = Inf,
                           show.legend = FALSE,
                           size=1.8) +# add labels to samples
  # manuscript.pca.plot.theme +  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(
    legend.title = element_text(vjust = 5, size = 10),
    legend.key.spacing = unit(0.01,"cm"),
    legend.box.spacing = unit(0.01,units="cm")
  )


#+ echo=FALSE
### PCA on NMR sexes ####
#' ### PCA on NMR sexes ####
pca.nmr.plot.sex <- readRDS(file.path(manuscript.figures,"pca.nmr.plot.sex.rds")) + 
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  geom_point(size = 2,
               color = "black" 
  ) +
  #               gg.title.taxon))+
  project_theme +
  manuscript.pca.plot.theme + 
  
  theme(legend.title = element_text(vjust = 2),
        legend.key.spacing = unit(0.01,"cm"),
        legend.box.spacing = unit(0.01,units="cm")
  )

#+ echo=FALSE
### Alpha diversity on NMR ages ####
#' ### Alpha diversity on NMR ages ####
agglom.rank <- "OTU"
metric.labs <- c('sobs'= case_when(agglom.rank == "OTU" ~ 
                                     "Richness (Observed species)",
                                   agglom.rank == "Genus" ~ 
                                     "Richness (Observed genera)") ,
                 'shannon' = "Shannon",
                 'invsimpson' = "Inverse Simpson")

nmr.age.div.plot <- readRDS(file = file.path(manuscript.figures,"nmr.age.div.plot.rds")) +
  facet_wrap(~factor(metric, # reorder facets
                     levels=names(metric.labs)),
             nrow=length(names(metric.labs)),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  geom_point(size=2,stat = "identity",colour = "white",shape=21)+
  labs(fill="Age (years)",
       x="Age (years)")+
  manuscript.alpha.plot.theme + # general theme
  theme(axis.title.x = element_text(size=11),
        axis.text.y=element_text(size=8),
        legend.box.spacing = unit(0.01,units="cm"),
        legend.title = element_text(vjust = 2))

#+ echo=FALSE
### NMR ASV barplot ####
#' ### NMR ASV barplot ####
top10.shared.asv.plot.filtered <- readRDS(file = file.path(manuscript.figures,"top10.shared.asv.plot.raw.rds"))+
  manuscript.asv.barplot.theme +
  theme(legend.key.spacing.x = unit(0.5, "cm"))

#+ echo=FALSE
### Comparison of my NMR data with Debebe et al. (Figure 6) ####
#' ### Comparison of my NMR data with Debebe et al. (Figure 6) ####
debebe.comparison.plot <- readRDS(file = file.path(manuscript.figures,
                                                   "debebe.comparison.plot.rds"))
debebe.comparison.plot <- debebe.comparison.plot + 
  theme(axis.text.x = element_text(angle=0,size=7,
                                   hjust=0.5,colour = "black"),# shift 
        axis.title.x=element_blank(),
        # the x-axis labels to the right
        # legend.position = "inside",
        axis.title.y = element_text(size = 10), # size of axis names
        axis.text.y = element_text(size=10,color="black"),
        strip.text.x = ggtext::element_markdown(size = 11),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        # legend.position.inside = c(0.9,0.8),
        legend.text = ggtext::element_markdown(size=9,margin = margin(l=0.02,
                                                                      unit="cm")),
        legend.title = element_text(size = 12), # size of legend title
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.spacing.y = unit(0.01, "cm"), # distant between key text
        legend.key.spacing.x = unit(0.05, "cm"), # distant between key text
        legend.box.spacing = unit(0.004,"cm"), # space between the plot and the legend box
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
  )


#+ echo=FALSE
## Save figures ####
#' ## Save figures ####
figure2<-plot_grid(div.plot.jitter,
                   pca.hosts.plot,
                   ncol=1,
                   labels = "auto",
                   label_size = 15,
                   rel_heights = c(1,1.4))

f3.top_row <- plot_grid(pca.nmr.plot.age.labeled, 
                        pca.nmr.plot.sex, 
                        labels = c('a', 'd'), label_size = 15)
figure3 <-plot_grid(f3.top_row, 
                    nmr.age.div.plot, 
                    top10.shared.asv.plot.filtered,
                    labels = c('', 'b','c'), label_size = 15, ncol = 1,
                    rel_heights = c(0.8,1.1,1))

for(image.format in image.formats){
  ggsave(filename = paste("figure2",image.format,sep = "."),
         path = manuscript.figures,
         plot = figure2,
         device = image.format,
         dpi=300,
         width = 170,
         height = 190,
         units = "mm",
         scale = 1)
}

for(image.format in image.formats){
  ggsave(filename = paste("figure3",image.format,sep = "."),
         path = manuscript.figures,
         plot = figure3,
         device = image.format,
         dpi=300,
         width = 170,
         height = 220,
         units = "mm",
         scale = 1)
}

for(image.format in image.formats){
  ggsave(filename = paste("figure6",image.format,sep = "."),
         path = manuscript.figures,
         plot = debebe.comparison.plot,
         device = image.format,
         dpi=300,
         width = 170, 
         height = 90,
         units = "mm")
}

sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()
