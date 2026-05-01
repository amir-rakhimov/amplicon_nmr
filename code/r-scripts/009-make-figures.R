#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 002-summary-stats-qiime2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("code/r-scripts/002-summary-stats-qiime2.R", knit = FALSE)
#' #file.rename("code/r-scripts/002-summary-stats-qiime2.Rmd", "markdown/002-summary-stats-qiime2.Rmd")
#' #rmarkdown::render('./markdown/002-summary-stats-qiime2.Rmd', 'html_document',
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
#' This script generates all figures
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
#' Load necessary scripts.
source("./code/r-scripts/get_n_uniq_taxa_per_host.R")
source("./code/r-scripts/create_summary_stats_table.R")
source("./code/r-scripts/add_relab_to_tax_df.R")
source("./code/r-scripts/add_agegroup_to_tax_df.R")
source("./code/r-scripts/get_unclassified_summary_stats.R")
source("./code/r-scripts/get_dominant_taxa_in_host.R")
source("./code/r-scripts/ggplot_species.R")

#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Use only the taxa that are present in the workspace
#' (custom.md is metadata from the rdafile).
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

#' Import the rarefied abundance table (rds  file):
# 20260211_17_14_19 rarefied table for all hosts, genus level
# 20260211_17_14_20 rarefied table file for NMR, genus level
# 20260211_17_14_21 rarefied table file for NMR, ASV level
#' Import rarefied dataframe and select Sample, Abundance, 
#' class, and agglom.rank columns.
ps.q.df.preprocessed.genus<-readRDS(
  file.path(rdafiles.directory,paste0(
    paste(
      "20260211_17_14_19",
      "ps.q.df.rare-nonfiltered","Genus",
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".rds")))%>%
  filter(class%in%custom.levels, Abundance!=0)%>%
  dplyr::select(Sample,Abundance,class,Genus)

ps.q.df.preprocessed.asv<-readRDS(
  file.path(rdafiles.directory,paste0(
    paste(
      "20260211_17_14_21",
      "ps.q.df.rare-nonfiltered","OTU",
      "NMR",sep = "-"),
    ".rds")))%>%
  filter(Abundance!=0)%>%
  dplyr::select(Sample,Abundance,class,OTU)





## Figure 2 ####
#+ echo=FALSE
## 3. Alpha diversity analysis. ####
#' 
#' ## Alpha diversity analysis.
#+ echo=FALSE
### 3.1 Compute alpha diversity metrics. ####
#' 
#' ### Compute alpha diversity metrics.
all.div<-ps.q.df.preprocessed.genus%>%
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
## 4. Plot alpha diversity metrics. ####
#'
#' ## Plot alpha diversity metrics. ####
#' Facet labels
agglom.rank<-"Genus"
metric.labs <- c('sobs'= case_when(agglom.rank == "OTU" ~ 
                                     "Richness \n(Observed species)",
                                   agglom.rank == "Genus" ~ 
                                     "Richness \n(Observed genera)") ,
                 'shannon' = "Shannon",
                 'invsimpson' = "Inverse \nSimpson")

div.plot<-all.div%>%
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
                     levels=names(metric.labs)),
             ncol=length(names(metric.labs)),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+
  # scale_color_manual(breaks = unname(pretty.level.names),
  #                    labels=unname(pretty.level.names))+
  # scale_x_discrete(labels=pretty.level.names,
  #                  limits=custom.levels)+ # rename boxplot labels (x axis)
  scale_fill_viridis_d(option="C")+
  stat_summary(fun=median, geom="point", shape=23, size=2, color="black",fill="white",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single"))+
  labs(fill = "Host",
       y="",
       x="")+
  # ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=10),
    axis.title = element_text(size = 10),
    strip.text.x = element_text(size=10),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    legend.box.spacing = unit(0.01,units="cm"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) + # general theme
  theme(axis.text.x = ggtext::element_markdown(size=10),
        legend.text = ggtext::element_markdown(size=10))

div.plot.jitter<-div.plot+
  geom_jitter(aes(colour=class_letter),
              width = 0.1, # jitter spread
              shape = 1, # empty dots
              size=0.7,
              show.legend = FALSE)+
  scale_color_viridis_d(direction = -1,option="C")+
  guides(color = "none")

#+ echo=FALSE
## 3. PCA ####
#'
#' ## PCA ####
#' Convert into wide format.
ps.q.df.wide.pca<-ps.q.df.preprocessed.genus%>%
  dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
#' Column names are ASVs and row names are sample IDs.
rownames(ps.q.df.wide.pca)<-ps.q.df.wide.pca$Sample
#' Remove the Sample column (the first column)
ps.q.df.wide.pca<-ps.q.df.wide.pca[,-1] 
ps.q.df.wide.pca<-as.matrix(ps.q.df.wide.pca)

ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,
                                method="rclr",
                                logbase = exp)
#' Calculate principal components
pca.q<-prcomp(ps.q.df.wide.pca.tfm)
#+ echo=FALSE
## 4. PCA Plot. ####
#' 
#' ## PCA Plot. ####
PC1<-pca.q$x[,1]
PC2<-pca.q$x[,2]
#' Check how much variance is explained by the first two principal
#' components:
perc.var<-round(100*summary(pca.q)$importance[2,1:2],2)
perc.var
#' Combine the principal component vectors into a dataframe:
PC1.df <- as.data.frame(PC1)%>%
  mutate(Sample = rownames(ps.q.df.wide.pca))
PC2.df <- as.data.frame(PC2)%>%
  mutate(Sample = rownames(ps.q.df.wide.pca))
#' Plot the PCA results:
pca.plot<-PC1.df%>%
  left_join(PC2.df)%>%
  left_join(custom.md)%>%
  mutate(class = factor(class, levels=custom.levels))%>%
  group_by(class)%>%
  mutate(class_id = cur_group_id(), # add numeric id
         class_name = pretty.level.names[class_id])%>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = class_name,
             fill = class_name)) +
  geom_point(aes(shape = class_name),
             size = 1.5,
             color = "black")+ 
  scale_shape_manual(values = 0:9,
                     breaks = pretty.level.names,
                     labels = unname(pretty.level.names))+
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha = 0.2,
               show.legend = FALSE,
               color = "black")+
  scale_fill_viridis_d(option = "C", 
                       labels = unname(pretty.level.names))+
  theme_bw()+
  # ggtitle(paste0("PCA between different rodents (",agglom.rank, " level)"))+
  labs(x = paste0("PC1 (", perc.var[1], "%)"),
       y = paste0("PC2 (", perc.var[2], "%)"),
       shape = "Host")+
  guides(fill="none")+
  theme(axis.text.x = element_text(angle=0,size=13),
        axis.text.y = element_text(size=13),
        axis.title = element_text(size = 15),
        legend.text = ggtext::element_markdown(size = 10),
        legend.title = element_text(size = 12),
        legend.key.spacing = unit(0.01,"cm"),
        legend.box.spacing = unit(0.01,units="cm"),
        plot.caption = ggtext::element_markdown(hjust = 0,
                                                size=12),
        plot.caption.position = "plot",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+ # general theme
  theme(legend.position = "bottom")+
  labs(shape ="")





## Figure 3 ####

#+ echo=FALSE
## 6. Add agegroup variable to NMR data (must run for plotting). ####
#'
#' ## Add agegroup variable to NMR data (must run for plotting).
custom.md.ages<-readRDS(file.path(rdafiles.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")

ps.q.agg.relab.nmr<-ps.q.agg.asv%>%
  filter(class=="NMR")
source("./code/r-scripts/add_agegroup_to_tax_df.R")
ps.q.agg.relab.nmr<-add_agegroup_to_tax_df(ps.q.agg.relab.nmr,"OTU",
                                           custom.md.ages)

#+ echo=FALSE
## 16. Analysis of naked mole-rat data ASVs. ####
#' 
#' ## Analysis of naked mole-rat data ASVs.
#' Give ASVs shorter names: Genus, "ASV", OTU, OTU number. For example, Allobaculum_ASV_22.
nmr.asv.names<-ps.q.agg.relab.nmr%>%
  dplyr::select(OTU,Genus)%>%
  group_by(Genus)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(Genus,OTU)%>%
  mutate(row.index=row_number())%>%
  mutate(ASV_name=paste(Genus,row.index,sep="_ASV_"))

#' The new name becomes the OTU column. The old name becomes OTU_old_name column.
ps.q.agg.relab.nmr<-ps.q.agg.relab.nmr%>%
  left_join(nmr.asv.names[,c("Genus","OTU","ASV_name")])%>%
  rename("OTU_old_name"="OTU",
         "OTU"="ASV_name")%>%
  relocate(OTU,.before = Sample)%>%
  relocate(OTU_old_name,.after = Genus)
#+ echo=FALSE
## 17. Plot ASVs in NMR. ####
#' 
#' ## Plot ASVs in NMR. 
#' Setup sample levels for NMR for barplots.
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  dplyr::select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))

#+ echo=FALSE
### 16.1 How many ASVs are shared between two age groups? ####
#' ### How many ASVs are shared between two age groups?
#' First, find ASVs in young samples
otu.young<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup0_10",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)
#' Next, find ASVs in old samples
otu.old<-ps.q.agg.relab.nmr%>%
  filter(agegroup=="agegroup10_16",Abundance!=0)%>%
  distinct(OTU,.keep_all = T)
#' 668 shared ASVs:
shared.otu<-intersect(otu.young$OTU,otu.old$OTU)

### 17.2 Most abundant shared ASVs in each age group. ####
#'
#' ### Most abundant shared ASVs in each age group.
top.shared.asvs.by_age<-ps.q.agg.relab.nmr%>%
  filter(OTU%in%shared.otu)%>%
  group_by(OTU,agegroup)%>%
  distinct(OTU,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  dplyr::select(Genus,OTU,agegroup,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
  ungroup()
head(top.shared.asvs.by_age)

top10.shared.asv.young<-top.shared.asvs.by_age%>%
  filter(agegroup=="agegroup0_10")%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  head(n=10)

top10.shared.asv.old<-top.shared.asvs.by_age%>%
  filter(agegroup=="agegroup10_16")%>%
  arrange(-MeanRelativeAbundanceAgegroup)%>%
  head(n=10)

#' Are top 10 most abundant ASVs same in two age groups?
setequal(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)
#' No. How many are common?
intersect(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU)
#' Seven ASVs
#'
#' Union of the top 10 ASVs in each of the two age groups.
top10.shared.asv.union<-sort(union(top10.shared.asv.young$OTU,top10.shared.asv.old$OTU))
#' Make the names shorter for the barplot
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
#' Prepare a custom fill with Polychrome package
set.seed(1)
otu.fill<-createPalette(nrow(top10.shared.asv.union),
                        seedcolors = c("#FF0000", "#00FF00", "#0000FF"),
                        range=c(30, 80))
names(otu.fill)<-top10.shared.asv.union$new_OTU

#+ echo=FALSE
### 17.3 Barplot of the most abundant ASVs. ####
#' 
#' ### Barplot of the most abundant ASVs.
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
  theme(axis.text.y = element_text(size=8), # size of y axis ticks
    axis.text.x = element_text(angle=45,size=8,hjust=1),# rotate 
    axis.title.x = element_text(size = 10), # size of axis names
    axis.title.y = element_text(size = 9, hjust= 0.2),
    strip.text.x = ggtext::element_markdown(size=10),
    legend.text = element_text(size = 7.2), # size of legend text
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size=10),
    legend.key.size = unit(0.3, 'cm'), #change legend key size
    legend.key.spacing.y = unit(0.01, "cm"), # distant between key text
    legend.key.spacing.x = unit(0.05, "cm"), # distant between key text
    legend.box.spacing = unit(0.01,"cm"), # space between the plot and the legend box
    panel.spacing = unit(0.8, "cm"), # increase distance between facets
    plot.title = element_text(size = 8), # size of plot title
    plot.caption = element_text(size=8), # size of plot caption
    legend.position = "bottom")



#+ echo=FALSE
## 4. PCA for NMR ####
#' 
#' ## PCA for NMR
#' Convert df into wide format.
agglom.rank<-"OTU"
metric.labs <- c('sobs'= case_when(agglom.rank == "OTU" ~ 
                                     "Richness (Observed species)",
                                   agglom.rank == "Genus" ~ 
                                     "Richness \n(Observed genera)") ,
                 'shannon' = "Shannon",
                 'invsimpson' =case_when(agglom.rank == "OTU" ~ 
                                           "Inverse Simpson",
                                         agglom.rank == "Genus" ~ 
                                           "Inverse \nSimpson"))

ps.q.df.wide.pca.nmr<-ps.q.df.preprocessed.asv%>%
  dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
  pivot_wider(names_from = all_of(agglom.rank), 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

#' Colnames are OTUs and rownames are sample IDs.
rownames(ps.q.df.wide.pca.nmr)<-ps.q.df.wide.pca.nmr$Sample
ps.q.df.wide.pca.nmr<-ps.q.df.wide.pca.nmr[,-1] # remove sample names
ps.q.df.wide.pca.nmr<-as.matrix(ps.q.df.wide.pca.nmr) # convert to matrix
ps.q.df.wide.pca.nmr.tfm<-decostand(ps.q.df.wide.pca.nmr,method="rclr",logbase = exp)

#### >if you want to exclude specific samples
ps.q.df.wide.pca.nmr<-ps.q.df.wide.pca.nmr[-which(rownames(ps.q.df.wide.pca.nmr)%in%c("M40")),]
# remove ASVs that are zero because their samples are removed
ps.q.df.wide.pca.nmr<-ps.q.df.wide.pca.nmr[,which(colSums(ps.q.df.wide.pca.nmr)!=0)]
ps.q.df.wide.pca.nmr.tfm<-decostand(ps.q.df.wide.pca.nmr,
                                method="rclr",
                                logbase = exp)
####<

# calculate principal components
pca.q.nmr<-prcomp(ps.q.df.wide.pca.nmr.tfm)
str(pca.q.nmr)
dim(pca.q.nmr$x)

# reverse the signs
pca.q.nmr$rotation<- -1*pca.q.nmr$rotation

# reverse the signs of the scores
pca.q.nmr$x<- -1*pca.q.nmr$x

#calculate total variance explained by each principal component
var_explained.nmr = pca.q.nmr$sdev^2 / sum(pca.q.nmr$sdev^2)

#+ echo=FALSE
### 4.1 PCA Plot for NMR. ####
#' 
#' ### PCA Plot for NMR.
PC1.nmr<-pca.q.nmr$x[,1]
PC2.nmr<-pca.q.nmr$x[,2]
perc.var<-round(100*summary(pca.q.nmr)$importance[2,1:2],2)

PC1.df.nmr<-as.data.frame(PC1.nmr)%>%mutate(Sample=rownames(ps.q.df.wide.pca.nmr))
PC2.df.nmr<-as.data.frame(PC2.nmr)%>%mutate(Sample=rownames(ps.q.df.wide.pca.nmr))

pca.plot.nmr<-PC1.df.nmr%>%
  left_join(PC2.df.nmr)%>%
  # left_join(custom.md)
  left_join(custom.md.ages)

#' Theme for alpha diversity plots:
alpha.plot.theme<-theme(#plot.margin=unit(c(1,1,1,2), 'cm'),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = ggtext::element_markdown(size=10),
  axis.text.y = element_text(size=10),
  axis.title = element_text(size = 13),
  strip.text.x = element_text(size=10),
  plot.title = element_text(size = 14),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 10),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())
#' Theme for PCA scatter plots:
scatter.plot.theme<-theme(
  axis.text.x = element_text(angle=0,size=10,hjust=1),
  axis.text.y = element_text(size=10),
  axis.title = element_text(size = 10),
  plot.title = element_text(size = 14),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  legend.position = "right",
  plot.caption = ggtext::element_markdown(hjust = 0, size=13),
  plot.caption.position = "plot",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank())

#+ echo=FALSE
### 3.6 Plot PCA on sexes. ####
#' 
#' ### Plot PCA on sexes.
pretty.level.names.sex <-
  c("female" = "Females",
    "male" = "Males")
custom.levels.sex <- names(pretty.level.names.sex)
pretty.level.names.sex <-pretty.level.names.sex[which(names(pretty.level.names.sex)%in%custom.levels.sex)]
gg.labs.name.sex <-"Host sex"
gg.title.groups.sex <-"groups"
pca.plot.nmr.sex<-pca.plot.nmr%>%
  ggplot(aes(x = PC1.nmr, y = PC2.nmr,
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
  geom_point(size = 2,
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
  theme(legend.title = element_text(vjust = 2),
        legend.key.spacing = unit(0.01,"cm"),
        legend.box.spacing = unit(0.01,units="cm")
  )
#+ echo=FALSE
### 3.7 Plot PCA on ages. ####
#' 
#' ### Plot PCA on ages.
pca.plot.nmr.age<-pca.plot.nmr%>%
  ggplot(aes(x = PC1.nmr, 
             y = PC2.nmr, 
             color = age,
             shape= class,
             fill = age))+
  scale_shape_manual(values = 21)+
  scale_fill_viridis_c(option = "D")+
  scale_color_viridis_c(option = "D")+
  guides(color = "none",
         shape = "none")+
  labs(fill="Age (years)")+
  geom_point(size = 2,
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
  theme(
    legend.title = element_text(vjust = 5),
        legend.key.spacing = unit(0.01,"cm"),
        legend.box.spacing = unit(0.01,units="cm")
  )

pca.plot.nmr.age.labeled<-pca.plot.nmr.age+
  # ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
  ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
                           max.overlaps = Inf,
                           show.legend = FALSE,
                           size=1.8) +# add labels to samples
  theme(legend.title = element_text(vjust = 5))



#+ echo=FALSE
## 3. Alpha diversity analysis on NMR. ####
#' 
#' ## Alpha diversity analysis on NMR.
#+ echo=FALSE
### 3.1 Compute alpha diversity metrics. ####
#' 
#' ### Compute alpha diversity metrics.
pretty.level.names.age<-names(table(custom.md.ages$old_agegroup))
names(pretty.level.names.age)<-names(table(custom.md.ages$agegroup))
custom.levels<-names(pretty.level.names.age)
gg.labs.name.age<-"Age group"
gg.title.groups.age<-"age groups"

nmr.sample.levels<-custom.md.ages%>%
  filter(class=="NMR")%>%
  # mutate(age=round(time_length(sampling_date-birthday,unit="years")))%>%
  select(Sample,age,class)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age," y)"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))
all.div.nmr<-ps.q.df.preprocessed.asv%>%
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

nmr.age.div.plot<-all.div.nmr%>%
  filter(metric%in%names(metric.labs))%>%
  left_join(nmr.sample.levels)%>%
  ggplot(aes(x=age,y=value,fill=age))+
  geom_point(size=2,stat = "identity",colour = "white",shape=21)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=names(metric.labs)),
             nrow=length(names(metric.labs)),
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
  theme(axis.title.x = element_text(size=11),
        axis.text.y=element_text(size=8),
        legend.box.spacing = unit(0.01,units="cm"),
        legend.title = element_text(vjust = 2))


# Save figures ####
library(cowplot)
figure2<-plot_grid(div.plot.jitter,pca.plot,ncol=1,
                   labels = "auto",label_size = 20,
                   rel_heights = c(1,1.4))
for(image.format in image.formats){
  ggsave(filename = paste(paste(format(Sys.time(),format="%Y%m%d"),
                                format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                          paste("figure2",image.format,sep = "."),sep = "-"),
         path = "images",
         plot = figure2,
         device = image.format,
         dpi=300,
         width = 170,
         height = 190,
         units = "mm",
         scale = 1)
}

f3.top_row <- plot_grid(pca.plot.nmr.age.labeled, 
                        pca.plot.nmr.sex, 
                        labels = c('a', 'd'), label_size = 20)

figure3<-plot_grid(f3.top_row, 
                   nmr.age.div.plot, 
                   top10.shared.asv.plot,
          labels = c('', 'b','c'), label_size = 20, ncol = 1,
          rel_heights = c(0.8,1.1,1))



# figure3<-plot_grid(top10.shared.asv.plot, ,
#                    nmr.age.div.plot, ,ncol=1,
#                    labels = "auto",label_size = 20)
for(image.format in image.formats){
  ggsave(filename = paste(paste(format(Sys.time(),format="%Y%m%d"),
                                format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                          paste("figure3",image.format,sep = "."),sep = "-"),
         path = "images",
         plot = figure3,
         device = image.format,
         dpi=300,
         width = 170,
         height = 220,
         units = "mm",
         scale = 1)
}


