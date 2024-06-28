# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install(c("phyloseq","Maaslin2"))
# install.packages(c("tidyverse","Polychrome","pheatmap"))
library(tidyverse)
library(phyloseq)
# library(Maaslin2)
# library(vegan)
library(pheatmap)
ps.q.agg.date_time<-"20240613_21_47_12"
agglom.rank<-"Genus"
# 20240613_21_47_12 phyloseq OTU table
# 20240613_21_42_48 phyloseq Genus table

# Maaslin2 workspace
# 20240621_17_32_57 for NMR ages
# 20240621_17_44_54 for NMR sexes
# Workspace dates for multiple comparisons
# maaslin.dates<-c("agegroup0_5"="20240613_15_16_22",
#                  "agegroup5_10"="20240613_15_18_13",
#                  "agegroup10_15"="20240613_15_09_12")
maaslin.dates<-c("agegroup0_10"="20240621_17_32_57")
maaslin.dates<-c("F"="20240621_17_44_54")

ref.level<-"agegroup0_10"
ref.level<-"NMR"

truncationlvl<-"234"
authorname<-"pooled"
source("./code/r-scripts/make_features_maaslin.R")

rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"

# Import ps.q.agg
ps.q.agg<-read.table(file.path("./output/rtables",authorname,paste(
  paste(ps.q.agg.date_time,
  "phyloseq-qiime",authorname,agglom.rank,read.end.type,truncationlvl,
  "table.tsv",sep="-"))
  ),header = T)
# custom.levels<-c("NMR",
#                  "B6mouse",
#                  "MSMmouse",
#                  "FVBNmouse",
#                  "DMR",
#                  "hare",
#                  "rabbit",
#                  "spalax",
#                  "pvo"#,
#                  # "NMRwt"
#                  )

# custom.levels<-c("agegroup0_5",
#                  "agegroup5_10",
#                  "agegroup10_15")
custom.levels<-c("agegroup0_10",
                 "agegroup10_16")
custom.levels<-c("F",
                 "M")
host<-"NMR"
# host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"

if(length(custom.levels)==2){
  forloop.length<-length(custom.levels)-1
}else{
  forloop.length<-length(custom.levels)
}


for(i in seq(1, forloop.length)){
  ref.level<-custom.levels[i]
  if(agglom.rank=="OTU"){
    load(file.path("./output/rdafiles",paste(
      maaslin.dates[ref.level],
      "maaslin2",host,rare.status,filter.status,agglom.rank,
      comparison,truncationlvl,
      paste(sort(custom.levels),collapse = '-'),
      "ref",ref.level,"workspace.RData",sep="-")))
  }else{
    load(file.path("./output/rdafiles",paste(
      maaslin.dates,
      "maaslin2",rare.status,filter.status,agglom.rank,
      paste(custom.levels,collapse = '-'),
      truncationlvl,"ref",ref.level,"workspace.RData",sep="-")))
  }
  # extract features with qvalue<0.05
  if(min(maaslin.fit_data$results$qval)<0.05){
    maaslin.signif.features<-maaslin.fit_data$results%>%
      filter(qval<0.05) # should be qval
  }else{
    maaslin.signif.features<-maaslin.fit_data$results%>%
      arrange(qval)%>%
      head(n = 100) # if no significant results found
  }
  
  # make features pretty
  if(agglom.rank=="OTU"){
    maaslin.signif.features$feature<-gsub("^X","",maaslin.signif.features$feature)
  }else{
    foo<-ps.q.agg
    foo<-foo%>%mutate("maaslin"=get(agglom.rank))
    foo<-make_features_maaslin(foo,"maaslin")
    foo<-unique(foo[,c("maaslin",agglom.rank)])
    maaslin.signif.features<-maaslin.signif.features%>%
      left_join(foo[,c("maaslin",agglom.rank)],by=c("feature"="maaslin"))%>%
      distinct()
    rm(foo)
    maaslin.signif.features$feature<-maaslin.signif.features[,agglom.rank]
    maaslin.signif.features<-subset(maaslin.signif.features, select=-get(agglom.rank))
  }
  # we had to make this exchange because maaslin output treats space and hyphen
  # as the same thing
  
  # extract features that are downregulated in all hosts simultaneously
  maaslin.signif.decreased<-maaslin.signif.features%>%
    as_tibble()%>%
    filter(coef<0)%>%
    group_by(feature)%>%
    filter(n()==length(custom.levels)-1)%>% 
    arrange(pval,feature)%>%
    mutate(n=n())%>%
    mutate(assoc.str=-log(qval)*sign(coef))#%>%
  # select(feature,assoc.str,name)
  
  # maaslin.signif.increased<-maaslin.signif.features%>%
  #   as_tibble()%>%
  #   filter(coef>0)%>%
  #   group_by(feature)%>%
  #   filter(n()==length(custom.levels)-1)%>% 
  #   arrange(feature)%>%
  #   mutate(n=n())%>%
  #   mutate(assoc.str=-log(qval)*sign(coef))#%>%
  # # select(feature,assoc.str,name)
  if(agglom.rank=="OTU"){
    table(maaslin.signif.decreased$feature%in%ps.q.agg$OTU)
  }else{
    table(maaslin.signif.decreased$feature%in%pull(ps.q.agg[,agglom.rank]))
  }
  
  if(agglom.rank=="OTU"){
    write.table(maaslin.signif.features,
                file=file.path("./output/rtables",authorname,paste(
                  paste(format(Sys.time(),format="%Y%m%d"),
                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                  "maaslin2",host,rare.status,
                  filter.status,agglom.rank,
                  comparison,truncationlvl,
                  paste(sort(custom.levels),collapse = '-'),
                  "ref",ref.level,"signif.tsv",sep="-")),
                row.names = F,sep = "\t")
  }else{
    write.table(maaslin.signif.features,
                file=file.path("./output/rtables",authorname,paste(
                  paste(format(Sys.time(),format="%Y%m%d"),
                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                  "maaslin2",rare.status,
                  filter.status,agglom.rank,
                  paste(custom.levels,collapse = '-'),truncationlvl,
                  "ref",ref.level,"signif.tsv",sep="-")),
                row.names = F,sep = "\t")
  }
}


# Heatmap ####
heatmap.df<-maaslin.signif.features%>%
  mutate(assoc.str=-log(qval)*sign(coef))%>%
  select(feature,assoc.str,value)%>%
  distinct(.,feature,.keep_all = TRUE)%>%
  pivot_wider(names_from = "feature",
              values_from = "assoc.str",
              values_fill = 0)%>%
  as.data.frame()
heatmap.df$value<-gsub("class","",heatmap.df$value)
heatmap.df<-heatmap.df%>%
  arrange(factor(value, levels=custom.levels))

rownames(heatmap.df)<-heatmap.df$value
heatmap.df<-heatmap.df%>%select(-value)

max_value <- ceiling(max(heatmap.df))
min_value <- ceiling(min(heatmap.df))
range_value <- max(c(abs(max_value),abs(min_value)))
breaks <- seq(-1*range_value, range_value, by = 1)
pheatmap(t(heatmap.df)[1:50,],
         color=colorRampPalette(c("blue", "grey90", "red"))(range_value * 2) ,
         # cellwidth = 5,
         #   cellheight = 5,
         # changed to 3
         # fontsize = 6,
         kmeans_k = NA,
         border = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         legend = TRUE,
         breaks = breaks, 
         display_numbers=TRUE)

#####

heatmap(as.matrix(t(heatmap.df)))

ps.q.df.maaslin.input%>%
  filter(get(agglom.rank)=="Prevotellaceae Family",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

ps.q.agg%>%
  filter(get(agglom.rank)=="Eubacteriaceae Family",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Uncultured (Eubacteriaceae)")


