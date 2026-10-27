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
message("The analysis focus is: ", active.analysis)
#' Suppress summarise info
options(dplyr.summarise.inform = FALSE)

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

dir.create(diversity.tables,recursive = TRUE)
dir.create(diversity.rdafiles,recursive = TRUE)
dir.create(file.path(diversity.figures, "pca"),recursive = TRUE)
dir.create(file.path(diversity.figures, "alpha"),recursive = TRUE)
dir.create(file.path(manuscript.figures),recursive = TRUE)
plot.types<-c("plot"="",
              "plot.with.labels"="-with-labels")

#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"Genus"
ps.q<-readRDS(file= ps.q.filtered.fname)

custom.md <- readRDS(custom.md.path)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

rclr_avg.fname <-file.path(diversity.rdafiles,paste0(
  paste("ps.q.rare.rclr.average",agglom.rank,paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))

ps.q.agg.asv.rare <-readRDS(file.path(community.composition.rdafiles,paste0(
  paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
)

#+ echo=FALSE
## 3. PCA on a list with rarefied subsamples. ####
#'
#' ## PCA on a list with rarefied subsamples. ####
get_average_rclr <- function(ps.rare.list,
                             ps,
                             agglom.rank,
                             out.fname){
  if(file.exists(out.fname)){
    message("Average RCLR dataset exists. Skipping")
    rclr_avg <-readRDS(out.fname)
  }else{
    
    if(agglom.rank == "OTU"){
      ps.q.agg.agglom_rank.rare<-ps.list
    }else{
      ps.q.agg.agglom_rank.rare<-map(ps.list, ~tax_glom(.x, taxrank = agglom.rank))
      
    }
    # Melt each subsample and convert into wide format
    # Column names are ASVs and row names are sample IDs.
    ps.q.agg.agglom_rank.rare.melt <-map(ps.q.agg.agglom_rank.rare, ~psmelt(.x) %>% 
                                           filter(Abundance !=0) %>% 
                                           dplyr::select(Sample, Abundance, all_of(c(agglom.rank)))%>%
                                           pivot_wider(names_from = all_of(agglom.rank),
                                                       values_from = Abundance,
                                                       values_fill = 0)%>%
                                           column_to_rownames("Sample"))
    
    
    rclr_list <- lapply(ps.q.agg.agglom_rank.rare.melt, function(wide.df){
      # The crucial step is to transform the compositional data:
      # Use a log-ratio transformation, typically the centered log-ratio (clr)
      # or isometric log-ratio (ilr) transformation.
      # The clr transformation is often preferred for interpretability, 
      # while ilr is used for statistical properties.
      # We will use robust CLR (RCLR) with natural log
      
      # Transform the raw abundances with rclr
      # Add taxa that are missing in each subsample
      # Bind columns
      # Sort the columns by name
      # Samples are rows, taxa are columns
      wide.df.rclr <-decostand(wide.df, method = "rclr", logbase = exp(1))
      
      if(agglom.rank == "OTU"){
        missing.taxa <- setdiff(unique(rownames(ps@tax_table)), 
                                colnames(wide.df.rclr))
      }else{
        missing.taxa <- setdiff(ps@tax_table%>%as.data.frame()%>%distinct(get(agglom.rank))%>%pull, 
                                colnames(wide.df.rclr))
      }
      
      
      missing.taxa.df <- matrix(nrow= nrow(wide.df.rclr),
                                ncol = length(missing.taxa),
                                data = 0,
                                dimnames = list(rownames(wide.df.rclr),missing.taxa))
      wide.df.rclr <- cbind(wide.df.rclr, missing.taxa.df)
      wide.df.rclr <- wide.df.rclr%>%
        dplyr::select(order(colnames(.)))
      
    })
    # Sum the RCLR and divide by length of the list
    rclr_avg <- Reduce ("+", rclr_list)/ length(rclr_list)
    saveRDS(rclr_avg,
            file = out.fname)
  }
  return(rclr_avg)
}

rclr_avg.hosts <- get_average_rclr(ps.q.agg.asv.rare, 
                                   ps.q,
                             agglom.rank,
                             rclr_avg.fname)
#' #### Alternative: single rarefaction 
#' ps.q.df.wide.pca <- ps.q.agg.agglom_rank.rare[[1]]
#' ps.q.df.wide.pca.tfm <- rclr_list[[1]]
#' #### OR raw????
#' ps.q.df.pca.input <- ps.q.agg
#' ####
#' 

make_pca_df <- function(wide.df){
  pca.result <- prcomp(wide.df, center = TRUE, scale. = FALSE)
  message("Calculating principal components:")
  print(str(pca.result))
  #' * `sdev` is standard deviations of principal components. Used to 
  #' construct the scree plot   
  #' * `rotation` is loadings  
  #' * `center` is values used for centering  
  #' `x` is scores (original data in the new coordinate system)  s
  message("dim() of the prcomp()$x:")
  print(dim(pca.result$x))
  
  message("Scree plot:")
  plot(pca.result)
  
  message("Total variance explained by each principal component:")
  var_explained = pca.result$sdev^2 / sum(pca.result$sdev^2)
  print(var_explained)
  
  message("Scree plot:")
  qplot(seq_along(1:nrow(wide.df)), var_explained) + 
    geom_line() + 
    xlab("Principal Component") + 
    ylab("Variance Explained") +
    ggtitle("Scree Plot") +
    ylim(0, 1)%>%
    print()
  
  message("Loading plot:")
  par(mar = c(10, 4, 2, 2) + 0.2) #
  print(barplot(pca.result$rotation[1:10,1],las=2))
  
  #' Reverse the signs:
  pca.result$rotation<- -1*pca.result$rotation
  
  message("Principal components (loadings):")
  head(pca.result$rotation[,1:10])%>%
    knitr::kable(format = "simple")%>%
    print()
  
  #' Reverse the signs of the scores:
  pca.result$x<- -1*pca.result$x
  
  message("The first six scores:")
  head(pca.result$x[,1:10])%>%
    knitr::kable(format = "simple")%>%
    print()
  
  #+ echo=FALSE
  ## 4. PCA Plot. ####
  #' 
  #' ## PCA Plot. ####
  PC1<-pca.result$x[,1]
  PC2<-pca.result$x[,2]
  
  message("Check how much variance is explained by the first two principal components:")
  perc.var<-round(100*summary(pca.result)$importance[2,1:2],2)
  print(perc.var)
  
  #' Combine the principal component vectors into a dataframe:
  PC1.df <- as.data.frame(PC1)%>%
    # mutate(Sample = rownames(ps.q.df.wide.pca))
    mutate(Sample = rownames(wide.df))
  PC2.df <- as.data.frame(PC2)%>%
    # mutate(Sample = rownames(ps.q.df.wide.pca))
    mutate(Sample = rownames(wide.df))
  PC.1_2 <-PC1.df%>%
    left_join(PC2.df, by ="Sample")%>%
    relocate(PC1, .after = Sample)
  return(list(pca.result = pca.result,
              PC.1_2 = PC.1_2,
              perc.var = perc.var ))
}

pca.hosts <- make_pca_df(rclr_avg.hosts)

#' Plot the PCA results:
pca.hosts.plot <- pca.hosts$PC.1_2%>%
  left_join(custom.md)%>%
  mutate(class = factor(class, levels=custom.levels))%>%
  group_by(class)%>%
  mutate(class_id = cur_group_id(), # add numeric id
         class_name = pretty.level.names[class_id])%>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = class_name,
             fill = class_name)) +
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
  labs(x = paste0("PC1 (", pca.hosts$perc.var[1], "%)"),
       y = paste0("PC2 (", pca.hosts$perc.var[2], "%)"),
       shape = "Host")+
  guides(fill="none")

#' Label each point by sample name:
pca.hosts.plot.with.labels <- pca.hosts.plot +
  ggrepel::geom_text_repel(aes(label = Sample),
                           show.legend = FALSE,
                           color = "black") 

#' Save plot as RDS
saveRDS(pca.hosts.plot, 
        file = file.path(manuscript.figures, "pca.hosts.plot.rds"))

pca.hosts.plot <- pca.hosts.plot+
  geom_point(aes(shape = class_name),
             size = 2,
             color = "black")+ 
  project_theme + # general theme
  pca.plot.theme +
  theme(legend.position = "bottom")+
  labs(shape ="")

#+ fig.height=8, fig.width=11
print(pca.hosts.plot + 
        ggtitle(paste0("PCA between different rodents (",agglom.rank, " level)"))+
        theme(plot.title = element_text(size = 14)))
#+ fig.height=8, fig.width=11
print(pca.hosts.plot.with.labels + 
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
           plot=get(paste0("pca.hosts.",names(plot.types[plot.type]))),
           width=11, height=8,units="in",
           dpi=300,device = image.format)
  }
}
#+ echo=FALSE
## 5. Find ASVs that contribute to PCs. ####
#' 
#' ## Find ASVs that contribute to PCs. ####
aload.hosts <- abs(pca.hosts$pca.result$rotation)
#' Below, 2 means columns; sweep divides each column of loadings and 
#' divides by the column sum.
head(sweep(aload.hosts,2,colSums(aload.hosts),"/"))[1:5,1:5] 
colSums(sweep(aload.hosts, 2, colSums(aload.hosts), "/"))
colSums(aload.hosts)
aload.hosts[1,1]/colSums(aload.hosts)[1]

pc.hosts.df<-as.data.frame(sweep(aload.hosts,2,colSums(aload.hosts),"/")[,1:2])
lapply(pc.hosts.df,max)
max.ind<-lapply(pc.hosts.df,which.max)
pc.hosts.df[max.ind$PC1,]
pc.hosts.df[max.ind$PC2,]

#' Extracting the first few loadings in each component
pc1.hosts.loadings<-pc.hosts.df%>%
  rownames_to_column(var=agglom.rank)%>%
  arrange(-PC1)%>%
  mutate(PC_num="PC1")%>%
  head()
pc2.hosts.loadings<-pc.hosts.df%>%
  rownames_to_column(var=agglom.rank)%>%
  arrange(-PC2)%>%
  mutate(PC_num="PC2")%>%
  head()

pc1.hosts.loadings[[agglom.rank]]
pc2.hosts.loadings[[agglom.rank]]

loadings.hosts.df <- pc1.hosts.loadings%>%
  rbind(pc2.hosts.loadings)%>%
  dplyr::select(-PC1,-PC2)%>%
  add_count(Genus)%>%
  mutate(PC_num=ifelse(n==2,"PC1_and_PC2",PC_num))%>%
  distinct(Genus,PC_num)
loadings.hosts.df  
#' Check mean relative abundance of loadings.
mean_sd_relab.all_hosts<- read.table(file.path(community.composition.tables,
                                               paste0("mean-max-min-sd-relab-all_hosts-",agglom.rank,".tsv")),
                                     sep = "\t", header = T)

loadings.hosts.df.rel_ab <- loadings.hosts.df %>%
  left_join(mean_sd_relab.all_hosts, by ="Genus")
loadings.hosts.df.rel_ab %>%
  knitr::kable(format = "simple")

loadings.hosts.df.rel_ab.fname <-file.path(diversity.tables,
                                       paste("pca-loadings-hosts-mean_relative_abundance.tsv",sep="-"))
if(!file.exists(loadings.hosts.df.rel_ab.fname)){
  write.table(loadings.hosts.df.rel_ab,
              file = loadings.hosts.df.rel_ab.fname,
              row.names = F,sep = "\t")
}


#+ echo=FALSE
## 6. Bonus: PCA on each RCLR subsample (no averaging). ####
#' 
#' ## Bonus: PCA on each RCLR subsample (no averaging).
pc.list.fname <- file.path(diversity.rdafiles,paste0(
  paste("pca_plot_list","Genus",paste0(rare.num_samples,"_iter"),
        "all_hosts",sep = "-"),".rds"))
rclr_list_100.fname <- file.path(diversity.rdafiles,paste0(
  paste("rclr_wide_table_list","Genus",paste0(rare.num_samples,"_iter"),
        "all_hosts",sep = "-"),".rds"))

if(file.exists(pc.list.fname) & file.exists(rclr_list_100.fname)){
  message("PCA plots exist. Skipping")
  pc.list <- readRDS(pc.list.fname)
  rclr_list_100 <- readRDS(rclr_list_100.fname)
}else{
  genus.rarefied.fname.all <- file.path(community.composition.rdafiles,paste0(
    paste("ps.q.rare","Genus",paste0(rare.num_samples,"_iter"),
          paste(custom.levels,collapse = '-'),sep = "-"),".rds"))
  if(file.exists(genus.rarefied.fname.all)){
    message("Rarefied data exists. Skipping")
    ps.rare_100.genus<- readRDS(genus.rarefied.fname.all)
  }else{
    ps.rare_100.genus<-map(ps.q.agg.asv.rare, ~tax_glom(.x, taxrank = "Genus"))
  }
  #' Convert phyloseq objects into wide table.
  ps.rare_100.genus.melt <-map(ps.rare_100.genus, ~psmelt(.x) %>% 
                                 filter(Abundance !=0) %>% 
                                 dplyr::select(Sample, Abundance, all_of(c(agglom.rank)))%>%
                                 pivot_wider(names_from = all_of(agglom.rank),
                                             values_from = Abundance,
                                             values_fill = 0)%>%
                                 column_to_rownames("Sample"))
  #' Perform RCLR transformation and make all tables have the same columns (same
  #' names in the same order) even if they have 0 abundance
  rclr_list_100 <- lapply(ps.rare_100.genus.melt, function(wide.df){
    wide.df.rclr <-decostand(wide.df, method = "rclr", logbase = exp(1))
    missing.taxa <- setdiff(ps.q@tax_table%>%as.data.frame()%>%distinct(get(agglom.rank))%>%pull, 
                            colnames(wide.df.rclr))
    #' Fill missing columns
    missing.taxa.df <- matrix(nrow= nrow(wide.df.rclr),
                              ncol = length(missing.taxa),
                              data = 0,
                              dimnames = list(rownames(wide.df.rclr),missing.taxa))
    wide.df.rclr <- cbind(wide.df.rclr, missing.taxa.df)
    wide.df.rclr <- wide.df.rclr%>%
      dplyr::select(order(colnames(.)))
  })
  #' Perform PCA on all tables
  pca.rclr_list_100 <- map(rclr_list_100, ~prcomp(.x, center = TRUE, scale. = FALSE))
  #' Convert prcomp objects into a dataframe
  pc.list <- lapply (pca.rclr_list_100, function (pca.res){
    PC1<-pca.res$x[,1]
    PC2<-pca.res$x[,2]
    perc.var<-round(100*summary(pca.res)$importance[2,1:2],2)
    print(perc.var)
    
    #' Combine the principal component vectors into a dataframe:
    PC1.df <- as.data.frame(PC1)%>%
      # mutate(Sample = rownames(ps.q.df.wide.pca))
      mutate(Sample = rownames(pca.res$x))
    PC2.df <- as.data.frame(PC2)%>%
      # mutate(Sample = rownames(ps.q.df.wide.pca))
      mutate(Sample = rownames(pca.res$x))
    PC.1_2 <-PC1.df%>%
      left_join(PC2.df, by ="Sample")%>%
      relocate(PC1, .after = Sample)
    pca.plot <- PC.1_2%>%
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
                   color = "black")
    return(list(pca.plot = pca.plot,
                perc.var = perc.var))
  })
  
}

#' Finally, create PCA plots.
pc.list.final <- map(pc.list, function (pc.list.item){
  pc.list.item$pca.plot+
    scale_fill_viridis_d(option = "C", 
                         labels = unname(pretty.level.names))+
    theme_bw()+
    # ggtitle(paste0("PCA between different rodents (",agglom.rank, " level)"))+
    labs(x = paste0("PC1 (", pc.list.item$perc.var[1], "%)"),
         y = paste0("PC2 (", pc.list.item$perc.var[2], "%)"),
         shape = "Host")+
    guides(fill="none")
  
})

dir.create(file.path(diversity.figures, "pca","100_iter" ))
for (i in seq_along(pc.list.final)){
  ggsave(plot =  pc.list.final[[i]],
         filename = paste0(i,".png"),  
         path = file.path(diversity.figures, "pca","100_iter" ),
         device = "png",
         width = 5,height = 6)
  
}

#+ echo=FALSE
## 7. PCA on naked mole-rat samples. ####
#' 
#' ## PCA on naked mole-rat samples.
#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"OTU"
ps.q.nmr <- readRDS(file= ps.q.filtered.fname)%>%
  subset_samples(class == "NMR")
#' Load metadata with age information:
custom.md<-readRDS(file.path(new.metadata.dir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat 16S rRNA gene sequencing")
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

ps.q.agg.asv.nmr.rare <- readRDS(file.path(community.composition.rdafiles,paste0(
  paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
        paste("NMR",collapse = '-'),sep = "-"),
  ".rds")))

rclr_avg.nmr.fname <- file.path(diversity.rdafiles,paste0(
  paste("ps.q.rare.rclr.average",agglom.rank,paste0(rare.num_samples,"_iter"),
        paste(custom.levels,collapse = '-'),sep = "-"),".rds"))

# Get an average RCLR dataset
rclr_avg.nmr <- get_average_rclr(ps.q.agg.asv.nmr.rare,
                                 ps.q,
                                 agglom.rank,
                                 rclr_avg.nmr.fname)
# Calculate principal components
pca.nmr <- make_pca_df(rclr_avg.nmr)

pca.nmr.plot <- pca.nmr$PC.1_2 %>%
  left_join(custom.md%>%filter(sequencing_type=="Naked mole-rat 16S rRNA gene sequencing"),
            by = "Sample")%>%
  mutate(first_collection_date = as.Date(first_collection_date),
         collection_year = year(first_collection_date))
# left_join(custom.md)


#+ echo=FALSE
### 7.1 PCA plot on sexes. ####
#' ### PCA plot on sexes. ####
pretty.level.names.sex <- c("female" = "Females",
                            "male" = "Males")
pretty.level.names.age <- c("agegroup0_10" = "agegroup0_10",
                            "agegroup10_16" = "agegroup10_16")

pca.nmr.plot.sex <- pca.nmr.plot%>%
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
  labs(shape= "Host sex")+
  labs(x=paste0("PC1 (", pca.nmr$perc.var[1], "%)"),
       y=paste0("PC2 (", pca.nmr$perc.var[2], "%)"))+
  theme_bw() 

saveRDS(pca.nmr.plot.sex,
        file.path(manuscript.figures,"pca.nmr.plot.sex.rds"))

pca.nmr.plot.sex <- pca.nmr.plot.sex +
  geom_point(size = 5,
             color = "black" 
  ) +
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  # ggtitle(paste("PCA between different",as.character(host.class[host]),gg.title.groups, "\n",
  #               gg.title.taxon))+
  project_theme +
  pca.plot.theme +  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(legend.text = ggtext::element_markdown(size = 20),
        legend.title = element_text(vjust = 2))

#+ echo=FALSE
### 7.2 PCA plot on ages. ####
#' ### PCA plot on ages. ####
pca.nmr.plot.age <- pca.nmr.plot%>%
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
  labs(x=paste0("PC1 (", pca.nmr$perc.var[1], "%)"),
       y=paste0("PC2 (", pca.nmr$perc.var[2], "%)"))+
  theme_bw()

saveRDS(pca.nmr.plot.age,
        file.path(manuscript.figures,"pca.nmr.plot.age.rds"))

pca.nmr.plot.age <- pca.nmr.plot.age +
  geom_point(size = 5,
             color = "black" 
  ) +
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")+
  # ggtitle(paste("PCA between different",as.character(host.class[host]),gg.title.groups, "\n",
  #               gg.title.taxon))+
  project_theme +
  pca.plot.theme +  # theme for scatter plots (PCA, PCoA, nMDS)
  theme(legend.title = element_text(vjust = 2))

pca.nmr.plot.batch <- pca.nmr.plot%>%
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
  labs(x=paste0("PC1 (", pca.nmr$perc.var[1], "%)"),
       y=paste0("PC2 (", pca.nmr$perc.var[2], "%)"))+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dashed", color = "grey")+
  geom_vline(xintercept=0,linetype="dashed", color = "grey")


#+ fig.width=11, fig.height=8
print(pca.nmr.plot.sex)
#+ fig.width=11, fig.height=8
print(pca.nmr.plot.age)

#' Label the samples
pca.nmr.plot.sex.labeled<-pca.nmr.plot.sex +
  ggrepel::geom_text_repel(aes(label=Sample),
                           show.legend = FALSE,size=7) # add labels to samples

pca.nmr.plot.age.labeled<-pca.nmr.plot.age+
  ggrepel::geom_text_repel(aes(label=paste0(Sample," (",age,")")),
                           max.overlaps = Inf,
                           show.legend = FALSE,
                           size=7) +# add labels to samples
  theme(legend.title = element_text(vjust = 2))

saveRDS(pca.nmr.plot.age.labeled,
        file.path(manuscript.figures,"pca.nmr.plot.age.labeled.rds"))

#+ fig.width=11, fig.height=8
print(pca.nmr.plot.sex.labeled)
#+ fig.width=11, fig.height=8
print(pca.nmr.plot.age.labeled)

pca.nmr.plots<-list("age" = c(pca.nmr.plot.age, pca.nmr.plot.age.labeled),
                "sex" = c(pca.nmr.plot.sex, pca.nmr.plot.sex.labeled)
)
for (pca_comparison in names(pca.nmr.plots)){
  for(image.format in image.formats){
    ggsave(paste0(paste("pca",
                        "NMR",
                        pca_comparison,agglom.rank,truncationlvl,
                        sep = "-"),".",image.format),
           plot = pca.nmr.plots[[pca_comparison]][1],
           path = file.path(diversity.figures,"pca"),
           width=11, height=8,units="in",
           dpi=300,device = image.format)
    ggsave(paste0(paste("pca-labeled",
                        "NMR",
                        pca_comparison,agglom.rank,truncationlvl,
                        sep = "-"),".",image.format),
           plot = pca.nmr.plots[[pca_comparison]][2],
           path = file.path(diversity.figures,"pca"),
           width=11, height=8,units="in",
           dpi=300,device = image.format)
  }
}

#+ echo=FALSE
### 7.3 Find ASVs that contribute to PCs. ####
#' 
#' ### Find ASVs that contribute to PCs.
mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-asv.tsv")

mean_sd_relab.nmr_ages.asv <- read.table(mean_sd_relab.nmr_ages.asv.fname,
                                         header = TRUE,
                                         sep = "\t")
aload.nmr<-abs(pca.nmr$pca.result$rotation)
head(sweep(aload.nmr,2,colSums(aload.nmr),"/"))
colSums(sweep(aload.nmr, 2, colSums(aload.nmr), "/"))
colSums(aload.nmr)
aload.nmr[1,1]/colSums(aload.nmr)[1]

pc.nmr.df<-as.data.frame(sweep(aload.nmr,2,colSums(aload.nmr),"/")[,1:2])
lapply(pc.nmr.df,max)
max.ind<-lapply(pc.nmr.df,which.max)
pc.nmr.df[max.ind$PC1,]
pc.nmr.df[max.ind$PC2,]

pc1.nmr.loadings <- pc.nmr.df%>%
  rownames_to_column(var="OTU")%>%
  arrange(-PC1)%>%
  head(n=3)%>%
  left_join(mean_sd_relab.nmr_ages.asv, by = "OTU")%>%
  distinct()%>%
  arrange(-PC1,-MeanRelativeAbundance)

pc2.nmr.loadings <- pc.nmr.df%>%
  rownames_to_column(var="OTU")%>%
  arrange(-PC2)%>%
  head(n=3)%>%
  left_join(mean_sd_relab.nmr_ages.asv, by = "OTU")%>%
  distinct()%>%
  arrange(-PC2,-MeanRelativeAbundance)

rbind(pc1.nmr.loadings, pc2.nmr.loadings)%>%
  dplyr::select(OTU,Species,MeanRelativeAbundance, sdRelativeAbundance,
                minRelativeAbundance, 
                maxRelativeAbundance, 
                n_samples)%>%
  group_by(Species,OTU)%>%
  mutate(ASV= paste(Species,cur_group_id()))%>%
  ungroup%>%
  mutate(OTU=ASV,
         OTU=gsub(" ","_",OTU))%>%
  distinct(OTU,.keep_all = )%>%
  knitr::kable(format = "simple")%>%
  print()


sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("active.analysis","markdown.dir")))
gc()