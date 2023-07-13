library(ggpicrust2)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)


# EC without conversions: turn ko_to_kegg to FALSE. ####

# Load metadata as a tibble
metadata <- read.table(paste0("./data/alldir-data/","filenames-single-pooled-raw-supercomp.tsv"),
                       header = T)
metadata<-metadata%>%
  filter(class=="NMR"|class=="SPFmouse")

# Load KO abundance as a data.frame
# data(ko_abundance)
enz_abundance <- read.table("./data/alldir-data/picrust/supercomp/NMR_SPFmouse/pred_metagenome_unstrat.tsv",
                           header = T,sep="\t")
colnames(enz_abundance)[which(colnames(enz_abundance)=="function.")]<-"predicted.function"
colnames(enz_abundance)<-gsub("X2D10", "2D10",colnames(enz_abundance))
colnames(enz_abundance)<-gsub("X2D14", "2D14",colnames(enz_abundance))
# foo<-enz_abundance[rowSums(enz_abundance[,-1] == 0) <= 3, ]

# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = enz_abundance %>%column_to_rownames("predicted.function"), 
                              metadata = metadata, 
                              group = "class", 
                              daa_method = "ALDEx2", 
                              select = NULL, 
                              reference = NULL)

# Filter results for ALDEx2_Kruskal-Wallace test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
daa_sub_method_results_df<-daa_sub_method_results_df%>%
  filter(p_adjust<0.05)%>%
  arrange(p_adjust)%>%
  slice(1:20)

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "EC", 
                                                          daa_results_df = daa_sub_method_results_df, 
                                                          ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = enz_abundance %>% column_to_rownames("predicted.function"), 
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



# df[c('First', 'Last')] <- str_split_fixed(df$player, '_', 2)

# foo <- str_split_fixed(enz_abundance$predicted.function, '\\.', 4)

enzymes.long<-enz_abundance%>%
  pivot_longer(!predicted.function,names_to = "Sample",values_to = "Abundance")


split.ec.levels<-function(ec.df){
  ec.levels<-data.frame(do.call('rbind', 
                     strsplit(as.character(ec.df$predicted.function),'.',
                              fixed=TRUE)))
  colnames(ec.levels)<-paste0("EC.level.",1:ncol(ec.levels))
  ec.levels$EC.level.1<-gsub("EC\\:","",ec.levels$EC.level.1)
  ec.df.with.lvls<-cbind(ec.df,ec.levels)
  return(ec.df.with.lvls)
}
# ec.levels<-split.ec.levels(enzymes.long)
# ec.levels<-split.ec.levels(enz_abundance)
# ec.levels<-
#   data.frame(do.call('rbind', 
#                      strsplit(as.character(enz_abundance$predicted.function),'.',
#                               fixed=TRUE)))

# foo<-within(enz_abundance, FOO<-data.frame(do.call('rbind', strsplit(as.character(enz_abundance$predicted.function), '.', fixed=TRUE))))
enzymes.with.EC<-split.ec.levels(enzymes.long)

enzymes.with.EC%>%
  left_join(metadata[,c("sample.id","class")],by=c("Sample"="sample.id"))%>%
  mutate(Abundance=as.numeric(Abundance))%>%
  # filter(Sample=="DCG6",EC.level.1==1,EC.level.2==1,EC.level.3==1)%>%
  # summarise(sum=sum(Abundance))
  # group_by(predicted.function)%>%
  group_by(Sample,EC.level.1,EC.level.2,EC.level.3)%>%
  summarise(Abundance=sum(Abundance))%>%
  filter(EC.level.1==1 | EC.level.1==3)%>%  
  View()

enz.mean.var<-enzymes.with.EC%>%
  left_join(metadata[,c("sample.id","class")],by=c("Sample"="sample.id"))%>%
  # filter(EC.level.1==3,EC.level.2==6,EC.level.3==4,EC.level.4==12)%>%
  group_by(class)%>%
  # summarise(sum=sum(Abundance))
  
  
  filter(EC.level.1==1 | EC.level.1==3)%>%  
  group_by(class,predicted.function)%>%
  mutate(mean.abundance=mean(Abundance))%>%
  mutate(var.abundance=var(Abundance))%>%
  select(predicted.function,mean.abundance,var.abundance)%>%
  distinct()%>%
  arrange(-mean.abundance)%>%
  left_join(daa_sub_method_results_df,by=c("predicted.function"="feature"))

top10.signif<-enz.mean.var%>%
  filter(!is.na(p_adjust))%>%
  arrange(p_adjust,-mean.abundance)%>%
  group_by(class)%>%
  slice(1:10)%>%
  ungroup()

colnames(top10.signif)[which(colnames(top10.signif)=="predicted.function")]<-"feature"

daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "EC", 
                                                          daa_results_df = top10.signif, 
                                                          ko_to_kegg = FALSE)
p <- pathway_errorbar(abundance = enz_abundance %>% column_to_rownames("predicted.function"), 
                      daa_results_df = daa_annotated_sub_method_results_df, 
                      Group = metadata$class, 
                      p_values_threshold = 0.05, 
                      order = "group",
                      select = daa_annotated_sub_method_results_df %>% 
                        arrange(p_adjust) %>% slice(1:20) %>% 
                        dplyr::select(feature) %>%
                        pull(), 
                      ko_to_kegg = FALSE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "description")

p