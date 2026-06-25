#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---
#' 
#' ```{r, setup 005-pca.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#+ echo=FALSE
# PCA ####
#' # PCA
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
#' ## Load necessary libraries.
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan","Polychrome"))
library(phyloseq)
library(tidyverse)
library(vegan) # Need v2.6-4 instead of current (2.7-1)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
# cmdargs <- commandArgs(trailingOnly = TRUE)
# active.analysis <- cmdargs[1]
unlockBinding("params", env = .GlobalEnv)
active.analysis <- params$active.analysis
print(paste("The analysis focus is:", active.analysis))

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

dir.create(diversity.tables,recursive = TRUE)
dir.create(diversity.rdafiles,recursive = TRUE)
dir.create(file.path(diversity.figures, "pca"),recursive = TRUE)
dir.create(file.path(diversity.figures, "alpha"),recursive = TRUE)
plot.types<-c("plot"="",
              "plot.with.labels"="-with-labels")

#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"Genus"
ps.q<-readRDS(file= ps.q.filtered.fname)

custom.md <- readRDS(custom.md.path)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

ps.q.agg.asv.rare <-readRDS(file.path(community.composition.rdafiles,paste0(
  paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
)

ps.q.agg.agglom_rank.rare<-map(ps.q.agg.asv.rare, ~tax_glom(.x, taxrank = agglom.rank))

ps.q.agg.agglom_rank.rare.melt <-map(ps.q.agg.agglom_rank.rare, ~psmelt(.x) %>% 
                           filter(Abundance !=0) %>% 
                             dplyr::select(Sample, Abundance, all_of(agglom.rank))%>%
                             pivot_wider(names_from = all_of(agglom.rank),
                                         values_from = Abundance,
                                         values_fill = 0)%>%
                             column_to_rownames("Sample"))


rclr_list<-lapply(ps.q.agg.agglom_rank.rare.melt, function(wide.df){
  # Samples are rows, taxa are columns
  wide.df.rclr <-decostand(wide.df, method = "rclr", logbase = exp(1))
  
  missing.taxa <- setdiff(ps.q@tax_table%>%as.data.frame()%>%distinct(get(agglom.rank))%>%pull, 
                          colnames(wide.df.rclr))
  
  missing.taxa.df <- matrix(nrow= nrow(wide.df.rclr),
                            ncol = length(missing.taxa),
                            data = 0,
                            dimnames = list(rownames(wide.df.rclr),missing.taxa))
  wide.df.rclr <- cbind(wide.df.rclr, missing.taxa.df)
  wide.df.rclr <- wide.df.rclr%>%select(order(colnames(.)))
  
})

rclr_avg <- Reduce ("+", rclr_list)/ length(rclr_list)
prcomp.rclr_avg<- prcomp(rclr_avg, center = TRUE, scale. = FALSE)
pca.q <-prcomp.rclr_avg

#' #### Alternative: single rarefaction 
#' ps.q.df.wide.pca <- ps.q.agg.agglom_rank.rare[[1]]
#' ps.q.df.wide.pca.tfm <- rclr_list[[1]]
#' #### OR raw????
#' ps.q.df.pca.input <- ps.q.agg
#' ####
#' 
#' if(agglom.rank=="OTU"){
#'   ps.q.agg <- ps.q%>%
#'     psmelt()
#' }else{
#'   ps.q.agg<-ps.q%>%
#'     tax_glom(taxrank = agglom.rank)%>%
#'     psmelt()%>%
#'     select(-OTU)
#' }
#' #' Remove rows with 0 Abundance, if there's any.
#' ps.q.agg<-ps.q.agg%>%
#'     filter(class%in%custom.levels,Abundance!=0)
#+ echo=FALSE
## 3. PCA ####
#'
#' ## PCA ####
#' Import rarefied dataframe.
# ps.q.df.pca.input <- readRDS(
#   file.path(community.composition.rdafiles,
#             paste0(
#     paste("ps.q.df.rare-nonfiltered",agglom.rank,
#       paste(custom.levels,collapse = '-'),sep = "-"),
#     ".rds")))
#' #' Convert into wide format.
#' ps.q.df.wide.pca<-ps.q.df.pca.input%>%
#'   dplyr::select(all_of(c(agglom.rank,"Sample","Abundance")))%>%
#'   pivot_wider(names_from = all_of(agglom.rank), 
#'               values_from = "Abundance",
#'               values_fill = 0)%>%
#'   as.data.frame()
#' #' Column names are ASVs and row names are sample IDs.
#' ps.q.df.wide.pca<-ps.q.df.wide.pca%>%
#'   column_to_rownames("Sample")%>%
#'   as.matrix()
#' #' The first crucial step is to transform the compositional data:
#' #' Use a log-ratio transformation, typically the centered log-ratio 
#' #' (clr) or isometric log-ratio (ilr) transformation.
#' #' The clr transformation is often preferred for interpretability, 
#' #' while ilr is used for statistical properties.
#' #' We will use robust CLR (RCLR) with natural log.
#' ps.q.df.wide.pca.tfm<-decostand(ps.q.df.wide.pca,
#'                                 method="rclr",
#'                                 logbase = exp(1))
#' #' Calculate principal components
#' pca.q<-prcomp(ps.q.df.wide.pca.tfm)
#' str(pca.q)
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
#'### Scree plot. 
plot(pca.q)
#' Calculate total variance explained by each principal component. 
pca.q$sdev^2 / sum(pca.q$sdev^2)
#' Calculate total variance explained by each principal component. 
var_explained = pca.q$sdev^2 / sum(pca.q$sdev^2)
#' Create a scree plot:
# qplot(seq_along(1:nrow(ps.q.df.wide.pca.tfm)), var_explained) + 
qplot(seq_along(1:nrow(rclr_avg)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
#+ echo=FALSE
### 3.2 Loading plot. #### 
#' 
#' ### Loading plot. #### 
par(mar = c(10, 4, 2, 2) + 0.2) #
barplot(pca.q$rotation[1:10,1],las=2)
#' Display principal components (loadings):
head(pca.q$rotation[,1:10])
#' Display the first six scores:
head(pca.q$x[,1:10])

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
  # mutate(Sample = rownames(ps.q.df.wide.pca))
  mutate(Sample = rownames(rclr_avg))
PC2.df <- as.data.frame(PC2)%>%
  # mutate(Sample = rownames(ps.q.df.wide.pca))
  mutate(Sample = rownames(rclr_avg))
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
  project_theme + # general theme
  pca.plot.theme+
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
for(plot.type in seq_along(plot.types)){
  for(image.format in image.formats){
    ggsave(filename = paste0(paste(
                        "pca",
                        paste(custom.levels,collapse = '-'),
                        agglom.rank,truncationlvl,sep = "-"),
                  plot.types[plot.type], ".",image.format),
           path = file.path(diversity.figures,"pca"),
           plot=get(paste0("pca.",names(plot.types[plot.type]))),
           width=11, height=8,units="in",
           dpi=300,device = image.format)
  }
}
#+ echo=FALSE
## 5. Find ASVs that contribute to PCs. ####
#' 
#' ## Find ASVs that contribute to PCs. ####
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
  dplyr::select(-PC1,-PC2)%>%
  add_count(Genus)%>%
  mutate(PC_num=ifelse(n==2,"PC1_and_PC2",PC_num))%>%
  distinct(Genus,PC_num)
loadings.df  
#' Check mean relative abundance of loadings.
mean_sd_relab.all_hosts<- read.table(file.path(community.composition.tables,
                                               paste0("mean-max-min-sd-relab-all_hosts-",agglom.rank,".tsv")),
                                     sep = "\t", header = T)

ps.q.agg.pc_loadings<-loadings.df %>%
  left_join(mean_sd_relab.all_hosts, by ="Genus")
ps.q.agg.pc_loadings%>%
  knitr::kable(format = "simple")

ps.q.agg.pc_loadings.fname <-file.path(diversity.tables,
                                       paste("pca-loadings-mean_relative_abundance.tsv",sep="-"))
if(!file.exists(ps.q.agg.pc_loadings.fname)){
  write.table(ps.q.agg.pc_loadings,
              file = ps.q.agg.pc_loadings.fname,
              row.names = F,sep = "\t")
}

sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()