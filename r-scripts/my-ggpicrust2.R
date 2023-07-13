library(ggpicrust2)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)

# foo<-KEGGREST::keggGet("ko04940")[[1]]

# Load necessary data: abundance data and metadata
abundance_file <- "./data/alldir-data/picrust/supercomp/NMR_FD/pred_metagenome_unstrat.tsv"

metadata <- read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  filter(class=="NMR"|class=="FukomysDamarensis")

# Run ggpicrust2 with input file path
# results_file_input <- ggpicrust2(file = abundance_file,
#                                  metadata = metadata,
#                                  group = "class", # For example dataset, group = "Environment"
#                                  pathway = "KO",
#                                  daa_method = "LinDA",
#                                  ko_to_kegg = TRUE,
#                                  order = "pathway_class",
#                                  p_values_bar = TRUE,
#                                  x_lab = "pathway_name")
# 
# # Run ggpicrust2 with imported data.frame
# abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
# 
# # Run ggpicrust2 with input data
# results_data_input <- ggpicrust2(data = abundance_data,
#                                  metadata = metadata,
#                                  group = "class", # For example dataset, group = "Environment"
#                                  pathway = "KO",
#                                  daa_method = "LinDA",
#                                  ko_to_kegg = TRUE,
#                                  order = "pathway_class",
#                                  p_values_bar = TRUE,
#                                  x_lab = "pathway_name")
# 
# # Access the plot and results dataframe for the first DA method
# example_plot <- results_file_input[[1]]$plot
# example_results <- results_file_input[[1]]$results

# KO2KEGG ####
# Load metadata as a tibble
# data(metadata)
metadata <- read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  filter(class=="NMR"|class=="FukomysDamarensis")

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("./data/alldir-data/picrust/supercomp/NMR_FD/pred_metagenome_unstrat.tsv") 

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance, 
                              metadata = metadata, 
                              group = "class", 
                              daa_method = "ALDEx2", 
                              select = NULL, 
                              reference = NULL) 

# Filter results for ALDEx2_Wilcoxon's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
daa_sub_method_results_df<-daa_sub_method_results_df%>%
  filter(p_adjust<0.05)%>%
  arrange(p_adjust)%>%
  slice(1:25)


pathways.test<-daa_sub_method_results_df
batches<-split(1:nrow(pathways.test), ceiling(seq_along(1:nrow(pathways.test))/50))
pathways.list<-list()

for(batch_index in 1:length(batches)){
  batch.item<-pathways.test[batches[[batch_index]],]
  for(i in batches[[batch_index]]){
    tryCatch({
      message('try part starts')
      pathway.kegg<-KEGGREST::keggGet(pathways.test$feature[i])
      
    }, error=function(cond){
      message("KEGG ",pathways.test$feature[i]," doesn't exist")
      message("Here's the original warning message:")
      message(cond)
      return(NA)
    })
    if(!is.na(pathway.kegg)){
      pathway.kegg<-pathway.kegg[[1]]
      batch.item$pathway.name[i]<-pathway.kegg$NAME
      # pathways.test$pathway.description[i]<-pathway.kegg$DESCRIPTION[i]
      batch.item$pathway.class[i]<-pathway.kegg$CLASS
      batch.item$pathway.map[i]<-pathway.kegg$PATHWAY_MAP[[1]]
    }else{
      batch.item$pathway.name[i]<-NA
      # pathways.test$pathway.description[i]<-pathway.kegg$DESCRIPTION[i]
      batch.item$pathway.class[i]<-NA
      batch.item$pathway.map[i]<-NA
    }
  }
  pathways.list[[batch_index]]<-batch.item
}

# (KEGGREST::keggGet("md:stm_M00841"))
# class(KEGGREST::keggGet("k0000000"))
# 
# kegg.foo<-"md:stm_M00841"
# tryCatch({
#   message('try part starts')
#   KEGGREST::keggGet(kegg.foo)
# }, error=function(cond){
#   message("KEGG doesn't exist")
#   message("Here's the original warning message:")
#   message(cond)
#   return(NA)
# }, warning=function(cond){
#   message(paste("URL caused a warning:", kegg.foo))
#   message("Here's the original warning message:")
#   message(cond)
#   # Choose a return value in case of warning
#   return(NULL)
# },finally={
#   # NOTE:
#   # Here goes everything that should be executed at the end,
#   # regardless of success or error.
#   # If you want more than one expression to be executed, then you 
#   # need to wrap them in curly brackets ({...}); otherwise you could
#   # just have written 'finally=<expression>' 
#   message(paste("Processed URL:", kegg.foo,"finish"))
#   message("Some other message at the end")
# })

b1<-daa_sub_method_results_df[batches[[1]],]
b2<-daa_sub_method_results_df[batches[[2]],]
b3<-daa_sub_method_results_df[batches[[3]],]
b4<-daa_sub_method_results_df[batches[[4]],]


b1.a<-pathway_annotation(pathway = "KO", 
                         daa_results_df = b1, 
                         ko_to_kegg = TRUE)
b2.a<-pathway_annotation(pathway = "KO", 
                         daa_results_df = b2, 
                         ko_to_kegg = TRUE)
b3.a<-pathway_annotation(pathway = "KO", 
                         daa_results_df = b3, 
                         ko_to_kegg = TRUE)
b4.a<-pathway_annotation(pathway = "KO", 
                         daa_results_df = b4, 
                         ko_to_kegg = TRUE)

b1234.a<-rbind(b1.a,b2.a,b3.a,b4.a)
b.foo<-b1234.a[grep("Metabolism",b1234.a$pathway_class),]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- 
  pathway_annotation(pathway = "KO", 
                     daa_results_df = daa_sub_method_results_df, 
                     ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = kegg_abundance, 
                      # daa_results_df = daa_annotated_sub_method_results_df, 
                      daa_results_df = b.foo,
                      Group = metadata$class, 
                      p_values_threshold = 0.05, 
                      order = "pathway_class", 
                      select = NULL, 
                      ko_to_kegg = TRUE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "pathway_name")

p



#  -----------------------------------------------------------------#
# EC, MetaCyc, and KO without conversions: turn ko_to_kegg to FALSE. ####

# Load metadata as a tibble
metadata <- read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  filter(class=="NMR"|class=="FukomysDamarensis")

# Load KO abundance as a data.frame
# data(ko_abundance)
ko_abundance <- read.table("./data/alldir-data/picrust/supercomp/NMR_FD/pred_metagenome_unstrat.tsv",
                           header = T,sep="\t")
colnames(ko_abundance)[which(colnames(ko_abundance)=="function.")]<-"predicted.function"
colnames(ko_abundance)<-gsub("X2D10", "2D10",colnames(ko_abundance))
colnames(ko_abundance)<-gsub("X2D14", "2D14",colnames(ko_abundance))

# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = ko_abundance %>%column_to_rownames("predicted.function"), 
                              metadata = metadata, 
                              group = "class", 
                              daa_method = "ALDEx2", 
                              select = NULL, 
                              reference = NULL)

# Filter results for ALDEx2_Kruskal-Wallace test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
daa_sub_method_results_df<-daa_sub_method_results_df%>%
  arrange(p_adjust)%>%
  slice(1:20)

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", 
                                                          daa_results_df = daa_sub_method_results_df, 
                                                          ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = ko_abundance %>% column_to_rownames("predicted.function"), 
                      daa_results_df = daa_annotated_sub_method_results_df, 
                      Group = metadata$class, 
                      p_values_threshold = 0.05, 
                      order = "group",
                      select = daa_annotated_sub_method_results_df %>% arrange(p_adjust) %>% slice(1:20) %>% dplyr::select(feature) %>% pull(), 
                      ko_to_kegg = FALSE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "description")

p

# --------------------------------------- #
# Workflow for MetaCyc Pathway and EC ####

# Load MetaCyc pathway abundance and metadata
metadata <- read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  filter(class=="NMR"|class=="FukomysDamarensis")

metacyc_abundance<-read.table("./data/alldir-data/picrust/supercomp/NMR_FD/path_abun_unstrat.tsv",
                              header = T,sep="\t")
colnames(metacyc_abundance)<-gsub("X2D10", "2D10",colnames(metacyc_abundance))
colnames(metacyc_abundance)<-gsub("X2D14", "2D14",colnames(metacyc_abundance))

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), 
                                      metadata = metadata, 
                                      group = "class", 
                                      daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = metacyc_daa_results_df, 
                                                       ko_to_kegg = FALSE)
metacyc_daa_annotated_results_df<-metacyc_daa_annotated_results_df%>%
  filter(p_adjust<0.05)%>%
  slice(1:25)
# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), 
                 daa_results_df = metacyc_daa_annotated_results_df, 
                 Group = metadata$class, 
                 ko_to_kegg = FALSE, 
                 p_values_threshold = 0.05, 
                 order = "group", 
                 select = NULL, 
                 p_value_bar = TRUE, 
                 colors = NULL, 
                 x_lab = "description")
library(ggh4x)
# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = metacyc_abundance %>% 
                  filter(pathway %in% feature_with_p_0.05$feature) %>%
                  column_to_rownames("pathway"), 
                metadata = metadata, 
                group = "class")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% 
              column_to_rownames("pathway"), 
            metadata = metadata, 
            group = "class")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "Maaslin2", "LinDA")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), 
              metadata = metadata, 
              group = "class",
              daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, 
                                          method_names = c("ALDEx2_Welch's t test", 
                                                           "ALDEx2_Wilcoxon rank test", 
                                                           "Maaslin2", 
                                                           "LinDA"))




#--------------------------------#
# Load metadata as a tibble
# data(metadata)
metadata <- read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  filter(class=="NMR")

# Load KEGG pathway abundance
# data(kegg_abundance)
# kegg_abundance <- ko2kegg_abundance("./data/alldir-data/picrust/supercomp/NMR/path_abun_unstrat.tsv") 

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
# daa_results_df <- pathway_daa(abundance = kegg_abundance, 
#                               metadata = metadata, 
#                               group = "class", 
#                               daa_method = "ALDEx2", 
#                               select = NULL, 
#                               reference = NULL) 

# Filter results for ALDEx2_Wilcoxon's t test method
# Please check the unique(daa_results_df$method) and choose one
# daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
# daa_sub_method_results_df<-daa_sub_method_results_df%>%
#   arrange(p_adjust)%>%
#   slice(1:20)

# Annotate pathway results using KO to KEGG conversion
annotated_results_df <- 
  pathway_annotation(file="./data/alldir-data/picrust/supercomp/NMR/path_abun_unstrat.tsv",
                     pathway = "MetaCyc", 
                     daa_results_df = kegg_abundance, 
                     ko_to_kegg = FALSE)
foo<-KEGGREST::keggGet("ko00010")[[1]]









# Just stuff ####
rep.seqs.q<-read_qza(paste0(qiimedir,"pooled","-rep-seqs-trimmed-dada2-",
                            truncationlvl,".qza"))
tree.q<-read_qza(paste0(qiimedir,"pooled","-rooted-tree-trimmed-dada2-",
                        truncationlvl,".qza"))

colnames(ko_abundance)[which(colnames(ko_abundance)=="function.")]<-"predicted.function"
foo<-ko_abundance[rowSums(ko_abundance[,-1] == 0) <= 3, ]

foo%>%
  # column_to_rownames("function")%>%
  pivot_longer(!predicted.function,names_to = "Sample",values_to = "abundance")%>%
  group_by(predicted.function)%>%
  mutate(abundance=as.numeric(abundance))%>%
  mutate(mean.abundance=mean(abundance))%>%
  mutate(var.abundance=var(abundance))%>%
  select(predicted.function,mean.abundance,var.abundance)%>%
  distinct()%>%
  # filter(mean.abundance>1000,var.abundance>100)%>%
  mutate(mean.int=mean.abundance %/% 1000 +1)%>%
  mutate(var.int=var.abundance %/% 1000 +1)%>%
  arrange(desc(mean.int),var.int)%>%View()

foo%>%
  pivot_longer(!predicted.function,names_to = "Sample",values_to = "abundance")%>%
  # filter(predicted.function=="K03088")%>%
  select(-predicted.function)%>%
  ggplot(aes(x=Sample,y=abundance))+
  geom_boxplot()




metacyc.foo<-metacyc_abundance[rowSums(metacyc_abundance[,-1] == 0) <= 3, ]

metacyc.foo%>%
  # column_to_rownames("function")%>%
  pivot_longer(!pathway,names_to = "Sample",values_to = "abundance")%>%
  group_by(pathway)%>%
  mutate(abundance=as.numeric(abundance))%>%
  mutate(mean.abundance=mean(abundance))%>%
  mutate(var.abundance=var(abundance))%>%
  select(pathway,mean.abundance,var.abundance)%>%
  distinct()%>%
  # filter(mean.abundance>1000,var.abundance>100)%>%
  mutate(mean.int=mean.abundance %/% 1000 +1)%>%
  mutate(var.int=var.abundance %/% 1000 +1)%>%
  arrange(desc(mean.int),var.int)%>%View()

metacyc.foo%>%
  pivot_longer(!pathway,names_to = "Sample",values_to = "abundance")%>%
  # filter(predicted.function=="K03088")%>%
  select(-pathway)%>%
  ggplot(aes(x=Sample,y=abundance))+
  geom_boxplot()






library(biomaRt)
mart <- useMart("ensembl")
View(listDatasets(mart))
mart <- useDataset("hgfemale_gene_ensembl", mart)
enzymes <-  getBM(attributes = c("go_id", "description"), 
                  
                  values = "GO:0005975", mart = mart)
