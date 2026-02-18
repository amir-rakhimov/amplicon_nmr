#' ---
#' title: "Beta diversity analysis"
#' output: 
#'   html_document:
#'      toc: true
#'      toc-location: left
#' ---
#' 
#' ```{r, setup 005-beta-diversity.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("code/r-scripts/005-beta-diversity.R", knit = FALSE)
#' #file.rename("code/r-scripts/005-beta-diversity.Rmd", "markdown/005-beta-diversity.Rmd")
#' #rmarkdown::render('./markdown/005-beta-diversity.Rmd', 'html_document',
#' # knit_root_dir="/home/rakhimov/projects/amplicon_nmr/")
#' ```
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script performs PCA to see how similar our samples are, 
#' then plot the results using PCA biplot. 
#' 
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## 1. Load necessary libraries and scripts.
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan","Polychrome"))
library(tidyverse)
library(vegan) # Need v2.6-4 instead of current (2.7-1)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## 2. Specifying parameters and directory/file names. 
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
#' Truncation level that we chose in QIIME2:
read.end.type<-"single"
#' Import abundance table from 001-phyloseq-qiime2.R as rds file
#' (NOT rarefied):
ps.q.agg.date_time<-"20260211_17_01_10"
ps.q.agg<-readRDS(file=file.path(
  "./output/rdafiles",
  paste(ps.q.agg.date_time,"phyloseq-qiime",authorname,agglom.rank,
        read.end.type, truncationlvl,"table.rds",sep = "-")))
# 20260211_17_01_07 ps.q.agg.date_time ASV level, all hosts
# 20260211_17_01_10 ps.q.agg.date_time genus level, all hosts
# 20260211_17_01_09 ps.q.agg.date_time family level, all hosts
# 20260211_17_01_08 ps.q.agg.date_time phylum level, all hosts
#'
#' Import the rarefied abundance table:
# ps.q.df.pca.input.date_time<-"20240426_22_00_04"
ps.q.df.pca.input.date_time<-"20260211_17_14_18"
#' Import metadata:
custom.md<-readRDS("./output/rdafiles/custom.md.rds")
#' This is for plot filenames:
plot.types<-c("plot"="",
              "plot.with.labels"="-with-labels")
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
#' Set a general theme for ggplot2. Add more parameters depending on the plot.
mytheme <- theme(plot.title = element_text(size = 27),
                 axis.text.x = element_text(angle=0,size=20),
                 axis.text.y = element_text(size=20),
                 axis.title = element_text(size = 20),
                 legend.text = ggtext::element_markdown(size = 15),
                 legend.title = element_text(size = 25),
                 legend.position = "right",
                 plot.caption = ggtext::element_markdown(hjust = 0,
                                                         size=20),
                 plot.caption.position = "plot",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()
)
#' Remove rows with 0 Abundance, if there's any.
ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
#+ echo=FALSE
## 3. PCA ####
#'
#' ## 3. PCA ####
#' Import rarefied dataframe.
ps.q.df.pca.input<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      ps.q.df.pca.input.date_time,
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
           header = T)
#' Convert into wide format.
ps.q.df.wide.pca<-ps.q.df.pca.input%>%
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
#' The first crucial step is to transform the compositional data:
#' Use a log-ratio transformation, typically the centered log-ratio 
#' (clr) or isometric log-ratio (ilr) transformation.
#' The clr transformation is often preferred for interpretability, 
#' while ilr is used for statistical properties.
#' We will use robust CLR (RCLR) with natural log.
ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,
                                method="rclr",
                                logbase = exp)
#' Calculate principal components
pca.q<-prcomp(ps.q.df.wide.pca.tfm)
str(pca.q)
#' * `sdev` is standard deviations of principal components. Used to 
#' construct the scree plot   
#' * `rotation` is loadings  
#' * `center` is values used for centering  
#' `x` is scores (original data in the new coordinate system)  s
dim(pca.q$x)
#'
#+ echo=FALSE
### 3.1 Scree plot. ####
#'
#'### 3.1 Scree plot. 
plot(pca.q)
#' Calculate total variance explained by each principal component. 
pca.q$sdev^2 / sum(pca.q$sdev^2)
#' Calculate total variance explained by each principal component. 
var_explained = pca.q$sdev^2 / sum(pca.q$sdev^2)
#' Create a scree plot:
qplot(seq_along(1:nrow(ps.q.df.wide.pca.tfm)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
#+ echo=FALSE
### 3.2 Loading plot. #### 
#' 
#' ### 3.2 Loading plot. #### 
par(mar = c(10, 4, 2, 2) + 0.2) #
barplot(pca.q$rotation[1:10,1],las=2)
#' Display principal components (loadings):
head(pca.q$rotation[,1:10])
#' Display the first six scores:
head(pca.q$x[,1:10])

#+ echo=FALSE
## 4. PCA Plot. ####
#' 
#' ## 4. PCA Plot. ####
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
             size = 2,
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
  mytheme+ # general theme
  theme(legend.position = "bottom")+
  labs(shape ="")

#' Label each point by sample name:
pca.plot.with.labels<-pca.plot+
  ggrepel::geom_text_repel(aes(label = Sample),
                           show.legend = FALSE,
                           color = "black") 
#+ fig.height=8, fig.width=11
print(pca.plot + 
        ggtitle(paste0("PCA between different rodents (",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14)))
#+ fig.height=8, fig.width=11
print(pca.plot.with.labels + 
        ggtitle(paste0("PCA between different rodents (",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14)))
# Save the plots
# for(plot.type in seq_along(plot.types)){
#   for(image.format in image.formats){
#     ggsave(filename = paste0(paste(paste(format(Sys.time(), format="%Y%m%d"),
#                               format(Sys.time(), format = "%H_%M_%S"),
#                               sep = "_"),
#                         "pca",
#                         paste(custom.levels,collapse = '-'),
#                         agglom.rank,truncationlvl,sep = "-"),
#                   plot.types[plot.type], ".",image.format),
#            path = "./images/diversity/pca",
#            plot=get(paste0("pca.",names(plot.types[plot.type]))),
#            width=11, height=8,units="in",
#            # width = 4500,height = 3000,units = "px",
#            dpi=300,device = image.format)
#   }
# }
#+ echo=FALSE
## 5. Find ASVs that contribute to PCs. ####
#' 
#' ## 5. Find ASVs that contribute to PCs. ####
aload<-abs(pca.q$rotation)
#' Below, 2 means columns; sweep divides each column of loadings and 
#' divides by the column sum.
head(sweep(aload,2,colSums(aload),"/"))[1:5,1:5] 
colSums(sweep(aload, 2, colSums(aload), "/"))
colSums(aload)
aload[1,1]/colSums(aload)[1]

pc.df<-as.data.frame(sweep(aload,2,colSums(aload),"/")[,1:2])
lapply(pc.df,max)
max.ind<-lapply(pc.df,which.max)
pc.df[max.ind$PC1,]
pc.df[max.ind$PC2,]

pc1.loadings<-pc.df%>%
  rownames_to_column(var=agglom.rank)%>%
  arrange(-PC1)%>%
  mutate(PC_num="PC1")%>%
  head()
pc2.loadings<-pc.df%>%
  rownames_to_column(var=agglom.rank)%>%
  arrange(-PC2)%>%
  mutate(PC_num="PC2")%>%
  head()

pc1.loadings$Genus
pc2.loadings$Genus

loadings.df<-pc1.loadings%>%
  rbind(pc2.loadings)%>%
  select(-PC1,-PC2)%>%
  add_count(Genus)%>%
  mutate(PC_num=ifelse(n==2,"PC1_and_PC2",PC_num))%>%
  distinct(Genus,PC_num)
loadings.df  
#' Check if the loadings are found in the data:
source("./code/r-scripts/add_relab_to_tax_df.R")
ps.q.agg.pc_loadings<-add_relab_to_tax_df(ps.q.agg,"Genus")%>%
  distinct(class,Genus,MeanRelativeAbundance,sdRelativeAbundance)%>%
  right_join(loadings.df)
# write.table(ps.q.agg.pc_loadings,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "pca-loadings-mean_relative_abundance.tsv",sep="-")),
#             row.names = F,sep = "\t")
ps.q.agg.pc_loadings