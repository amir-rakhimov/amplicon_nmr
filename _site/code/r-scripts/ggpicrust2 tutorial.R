library(ggpicrust2)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)


# Load necessary data: abundance data and metadata
abundance_file <- "./data/pooled-data/picrust/pred_metagenome_unstrat.tsv"

metadata <- read_delim(
  "./data/pooled-data/picrust/metadata.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Run ggpicrust2 with input file path
results_file_input <- ggpicrust2(file = abundance_file,
                                 metadata = metadata,
                                 group = "Environment", # For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Run ggpicrust2 with imported data.frame
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Run ggpicrust2 with input data
results_data_input <- ggpicrust2(data = abundance_data,
                                 metadata = metadata,
                                 group = "Environment", # For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Access the plot and results dataframe for the first DA method
example_plot <- results_file_input[[1]]$plot
example_results <- results_file_input[[1]]$results

# Use the example data in ggpicrust2 package
data(ko_abundance)
data(metadata)
results_file_input <- ggpicrust2(data = ko_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Analyze the EC or MetaCyc pathway
data(metacyc_abundance)
results_file_input <- ggpicrust2(data = metacyc_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "MetaCyc",
                                 daa_method = "LinDA",
                                 ko_to_kegg = FALSE,
                                 order = "group",
                                 p_values_bar = TRUE,
                                 x_lab = "description")
results_file_input[[1]]$plot
results_file_input[[1]]$results


# Step by step ####
# Load metadata as a tibble
# data(metadata)
metadata <- read_delim(
  "./data/pooled-data/picrust/metadata.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("./data/pooled-data/picrust/pred_metagenome_unstrat.tsv") 

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance, 
                              metadata = metadata, 
                              group = "Environment", 
                              daa_method = "ALDEx2", 
                              select = NULL, 
                              reference = NULL) 

# Filter results for ALDEx2_Wilcoxon's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- 
  pathway_annotation(pathway = "KO", 
                     daa_results_df = daa_sub_method_results_df, 
                     ko_to_kegg = TRUE)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = kegg_abundance, 
                      daa_results_df = daa_annotated_sub_method_results_df, 
                      Group = metadata$Environment, 
                      p_values_threshold = 0.05, 
                      order = "pathway_class", 
                      select = NULL, 
                      ko_to_kegg = TRUE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "pathway_name")

# If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE. ####

# Load metadata as a tibble
# data(metadata)
metadata <- read_delim(
  "./data/pooled-data/picrust/metadata.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Load KO abundance as a data.frame
# data(ko_abundance)
ko_abundance <- read.delim("./data/pooled-data/picrust/pred_metagenome_unstrat.tsv")

# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = ko_abundance %>%column_to_rownames("X.NAME"), 
                              metadata = metadata, 
                              group = "Environment", 
                              daa_method = "ALDEx2", 
                              select = NULL, 
                              reference = NULL)

# Filter results for ALDEx2_Kruskal-Wallace test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", 
                                                          daa_results_df = daa_sub_method_results_df, 
                                                          ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = ko_abundance %>% column_to_rownames("X.NAME"), 
                      daa_results_df = daa_annotated_sub_method_results_df, 
                      Group = metadata$Environment, 
                      p_values_threshold = 0.05, 
                      order = "group",
                      select = daa_annotated_sub_method_results_df %>% arrange(p_adjust) %>% slice(1:20) %>% dplyr::select(feature) %>% pull(), 
                      ko_to_kegg = FALSE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "description")

# Workflow for MetaCyc Pathway and EC ####

# Load MetaCyc pathway abundance and metadata
data("metacyc_abundance")
data("metadata")

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), 
                                      metadata = metadata, 
                                      group = "Environment", 
                                      daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = metacyc_daa_results_df, 
                                                       ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), 
                 daa_results_df = metacyc_daa_annotated_results_df, 
                 Group = metadata$Environment, 
                 ko_to_kegg = FALSE, 
                 p_values_threshold = 0.05, 
                 order = "group", 
                 select = NULL, 
                 p_value_bar = TRUE, 
                 colors = NULL, 
                 x_lab = "description")

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))
