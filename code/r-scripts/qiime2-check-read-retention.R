library(tidyverse)
active.analysis<-"all_domains"
source("config/R/config.R")
custom.md<-readRDS(file = "data/metadata/shared/custom.md.rds")
qiime2.output.dir <-"results/1-qiime2/20260521_191423-shared-paired_end/01-qiime2_output"
# list.files(stats.dir)
# read.table(file.path(stats.dir,list.files(stats.dir)[1],"stats.tsv"), 
#            sep ="\t",header = T)
# authors <- list.files(qiime2.output.dir)
# authors <- c("yasuda", "bensch", "liu", "shanmuganandam", "sibai")

read_stats_files <- function(filepath){
  stats.table <- read.table(filepath, sep = "\t", header = T)
  stats.table$filename <- dirname (filepath)
  return (stats.table)
}

files <-list.files(qiime2.output.dir, full.names = T, recursive = T,pattern = "stats\\.tsv$")
df <- do.call(rbind, lapply(files, read_stats_files))
df$filename <- gsub(paste0(qiime2.output.dir,"/"), "", df$filename)
df$filename <- gsub("exported-stats\\/trimmed-dada2-stats-", "", df$filename)
all.stats <- df %>% 
  as_tibble()%>%
  separate_wider_delim(filename, delim = "/", 
                     names = c("author", "params"))%>%
  rename("Sample"= "sample.id")%>%
  left_join(custom.md[,c("Sample","class")])%>%
  mutate(class = as.character(class),
         class = ifelse(Sample == "Y29", "NMR", class))%>%
  # mutate(retained = non.chimeric/input *100)%>%
  mutate(class = factor(class, levels = names(pretty.level.names)))%>%
  arrange(params,class,percentage.of.input.non.chimeric,Sample)%>%
  mutate(params = factor(params, levels = unique (params)),
         Sample = factor(Sample, levels = unique(Sample)))#%>%
# group_by(params, class)%>%
# mutate(mean_retention = mean(percentage.of.input.non.chimeric))%>%
# ungroup()


non_chimeric.retention.plot <- all.stats%>%
  filter(params != "284_203")%>%
  separate_wider_delim(params,
                       delim = "_",names = c("Fwd","Rev"))%>%
  mutate(Fwd = paste("Forward:", Fwd),
         Rev = paste("Reverse:", Rev))%>%
  # filter(params %in% c("265_220", "270_210","259_199", 
  #                      "225_225", "230_225"))%>%
  ggplot(aes(x = Sample, 
             y = percentage.of.input.non.chimeric,
             # color = params,
             shape = class))+
  geom_point(size = 2)+
  scale_shape_manual(values = c(7:10, 12:20),
                     labels = pretty.level.names)+
  scale_color_viridis_d()+
  facet_wrap( ~ Fwd + Rev, ncol = 4)+
  geom_hline(yintercept = c(0,50,75),
             linetype = "dashed",
             colour = rep(viridis::magma(3, end = 0.8), 
                          length(unique(all.stats$params))-1))+
  theme_bw()+
  # facet_wrap(~ params, ncol = 3)+
  labs(y = "Percentage of input non-chimeric",
       shape = "Host",
       # color = "Truncation parameters\n(Forward _ Reverse)",
       x = ""
    )+
  theme(axis.text.x = element_text(angle = 45, hjust =1, size = 4),
        legend.text = ggtext::element_markdown(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(filename = "non_chimeric-retention-plot.png",
       plot = non_chimeric.retention.plot,
       path = "results/1-qiime2/20260521_191423-shared-paired_end",
       device = "png",
       width = 30, height = 13, units = "in",
       dpi = 300)


all.stats%>%
  ggplot(aes (x = params,
              y = percentage.of.input.non.chimeric,
              # color = Sample,
              shape = class))+
  # geom_point()+
  geom_boxplot()+
  scale_shape_manual(values = c(12:20))+
  facet_wrap(~ class)+
  guides(color="none")+
  theme(axis.text.x = element_text(angle = 45, hjust =1))


pipeline.loss.boxplot <- all.stats%>%
  filter(params != "284_203")%>%
  separate_wider_delim(params,
                       delim = "_",names = c("Fwd","Rev"),cols_remove = F)%>%
  mutate(Fwd = paste("Forward:", Fwd),
         Rev = paste("Reverse:", Rev))%>%
  pivot_longer(cols = c("percentage.of.input.passed.filter",
                        "percentage.of.input.merged",
                        "percentage.of.input.non.chimeric"),
               names_to = "pipeline_step",
               values_to = "percentage")%>%
  mutate(pipeline_step = factor(pipeline_step, 
                                levels = c("percentage.of.input.passed.filter",
                                           "percentage.of.input.merged",
                                           "percentage.of.input.non.chimeric")),
         pipeline_step_for_plot = 
           case_when(pipeline_step == "percentage.of.input.passed.filter" ~"Passed filter",
                     pipeline_step == "percentage.of.input.merged" ~ "Merged",
                     pipeline_step == "percentage.of.input.non.chimeric" ~ "Non-chimeric"),
         pipeline_step_for_plot= factor(pipeline_step_for_plot, 
                                        levels = c("Passed filter",
                                                   "Merged",
                                                   "Non-chimeric")))%>%
  ggplot(aes(x = params, 
             y = percentage,
             # color = params,
             fill = pipeline_step_for_plot 
             # shape = class
             ))+
  geom_boxplot(position = position_dodge(),
               outliers = FALSE)+
  geom_point(
             position = position_jitterdodge(jitter.width = 0.12,),
             size = 1,
             shape = 1)+
  scale_fill_viridis_d(begin = 0.5,
                       option = "C")+
  facet_wrap(~ class,
             ncol =3,
             labeller = as_labeller(pretty.level.names))+
  theme_bw()+
  labs(x = "Truncation parameters",
       y = "Percentage of input",
       fill = "Pipeline step"
  )+
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        legend.position = "bottom",
        strip.text.x = ggtext::element_markdown(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(filename = "pipeline-loss-boxplot.png",
       plot = pipeline.loss.boxplot,
       path = "results/1-qiime2/20260521_191423-shared-paired_end",
       device = "png",
       width = 21, height = 12, units = "in",
       dpi = 300)

####

all.stats%>%
  filter(params %in% c("225_225", "230_225"))%>%
  select(Sample, percentage.of.input.non.chimeric, params)%>%
  pivot_wider(names_from = "params",
              values_from = percentage.of.input.non.chimeric)%>%
  mutate(perc_diff = `225_225` - `230_225`)%>%
  arrange(desc(perc_diff))

##### 

dfsingle<-read.table("results/1-qiime2/20260518_214705-shared-paired_end/01-qiime2_output/single-trimmed-dada2-stats-234.tsv",sep="\t",header = T)
df259_199<-read.table("results/1-qiime2/20260518_214705-shared-paired_end/01-qiime2_output/trimmed-dada2-stats-259_199.tsv",sep="\t", header=T)
df273_214<-read.table("results/1-qiime2/20260518_214705-shared-paired_end/01-qiime2_output/trimmed-dada2-stats-273_214.tsv",sep="\t", header=T)

dfsingle<-dfsingle%>%
  as_tibble()%>%
  rename("Sample"= "sample.id")%>%
  left_join(custom.md[,c("Sample","class")])%>%
  filter(class!="spalax")
df259_199<-df259_199%>%
  as_tibble()%>%
  rename("Sample"= "sample.id")%>%
  left_join(custom.md[,c("Sample","class")])%>%
  filter(class!="spalax")
df273_214<-df273_214%>%
  as_tibble()%>%
  rename("Sample"= "sample.id")%>%
  left_join(custom.md[,c("Sample","class")])%>%
  filter(class!="spalax")


dfsingle%>%
  mutate(retained = non.chimeric/input *100)%>%
  arrange(class,retained,Sample)%>%
  mutate(class = factor(class, levels = unique(class)),
         Sample = factor(Sample, levels = Sample))%>%
  ggplot(aes(x=Sample, 
             y=retained,
             color = class))+
  geom_point()

df273_214%>%
  mutate(retained = non.chimeric/input *100)%>%
  arrange(class,retained,Sample)%>%
  mutate(class = factor(class, levels = unique(class)),
         Sample = factor(Sample, levels = Sample))%>%
  ggplot(aes(x=Sample, 
             y=retained,
             color = class))+
  geom_point()

df259_199%>%
  mutate(author = "pooled",
         params = "259_199")%>%
  relocate(class,.after = last_col())%>%
  # group_by(params, class)%>%
  # mutate(mean_retention = mean(percentage.of.input.non.chimeric))%>%
  # ungroup()%>%
  # rbind(all.stats)%>%
  arrange(params,class,percentage.of.input.non.chimeric,Sample)%>%
  mutate(params = factor(params, levels = unique (params)),
         class = factor(class, levels = unique(class)),
         Sample = factor(Sample, levels = unique(Sample)))%>%
  ggplot(aes(x=Sample, 
           y=percentage.of.input.non.chimeric,
           color = params,
           shape = class))+
  geom_point(size = 2)+
  scale_shape_manual(values = c(12:19))+
  scale_color_viridis_d()+
  theme_bw()+
  facet_wrap(~params)
# facet_wrap(author~params,
#            scales =  "free_x")




df1<- df259_199%>%
  filter(class %in%c("NMR","B6mouse","FVBNmouse","MSMmouse"))%>%
  mutate(params="259_199")%>%
  select(-class)%>%
  arrange(Sample)

df2<- stats.list[["259_199"]]%>%  
  as_tibble()%>%
  rename("Sample"= "sample.id")%>%
  filter(Sample !="Y29")%>%
  arrange(Sample)

table(df1$percentage.of.input.passed.filter == df2$percentage.of.input.passed.filter)
table(round(df1$percentage.of.input.merged - df2$percentage.of.input.merged,digits = 3))
table(round(df1$percentage.of.input.non.chimeric - df2$percentage.of.input.non.chimeric,digits = 3))
