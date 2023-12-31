library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(pheatmap)
ref.level<-"NMR"

truncationlvl<-"234"
agglom.rank<-"Genus"
agglom.rank<-"OTU"
source("r-scripts/make_features_maaslin.R")

rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"

custom.levels<-c("NMR",
                 "B6mouse",
                 # "MSMmouse",
                 # "FVBNmouse",
                 "DMR",
                 "hare",
                 "rabbit",
                 "spalax",
                 "pvo"#,
                 # "NMRwt"
                 )


host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
comparison<-"sex"
comparison<-"strain"

# load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
#                                 truncationlvl,agglom.rank,
#                                 "phyloseq-workspace.RData",sep = "-")))

if(agglom.rank=="OTU"){
  load(file.path("./rdafiles",
              paste("maaslin",rare.status,filter.status,host,agglom.rank,
                    comparison,truncationlvl,ref.level,
                    "workspace.RData",sep="-")))
}else{
  load(file.path("./rdafiles",
              paste("maaslin",rare.status,filter.status,agglom.rank,
                    paste(custom.levels,collapse = '-'),
                    truncationlvl,ref.level,"workspace.RData",sep="-")))
}

# extract features with qvalue<0.05
maaslin.signif.features<-maaslin.fit_data$results%>%
  filter(qval<0.05)
# make features pretty
if(agglom.rank=="OTU"){
  maaslin.signif.features$feature<-gsub("^X","",maaslin.signif.features$feature)
}else{
  foo<-ps.q.agg
  foo$maaslin<-foo$Taxon
  foo<-make_features_maaslin(foo,"maaslin")
  foo<-unique(foo[,c("maaslin","Taxon")])
  maaslin.signif.features<-maaslin.signif.features%>%
    left_join(foo[,c("maaslin","Taxon")],by=c("feature"="maaslin"))%>%
    distinct()
  rm(foo)
  maaslin.signif.features$feature<-maaslin.signif.features$Taxon
  maaslin.signif.features<-subset(maaslin.signif.features, select=-Taxon)
}
# we had to make this exchange because maaslin output treats space and hyphen
# as the same thing

# extract features that are downregulated in all hosts simultaneously
maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by(feature)%>%
  filter(n()==length(custom.levels)-1)%>% 
  arrange(feature)%>%
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
  table(maaslin.signif.decreased$feature%in%ps.q.agg$Taxon)
}

if(agglom.rank=="OTU"){
  write.table(maaslin.signif.features,
              file=file.path("./rtables",authorname,
                          paste("maaslin",rare.status,
                                filter.status,host,agglom.rank,
                                comparison,truncationlvl,
                                ref.level,"signif.tsv",sep="-")),
              row.names = F,sep = "\t")
}else{
  write.table(maaslin.signif.features,
              file=file.path("./rtables",authorname,
                          paste("maaslin",rare.status,
                                filter.status,agglom.rank,
                                paste(custom.levels,collapse = '-'),truncationlvl,
                                ref.level,"signif.tsv",sep="-")),
              row.names = F,sep = "\t")
}



# Heatmap ####
heatmap.df<-maaslin.signif.features%>%
  mutate(assoc.str=-log(qval)*sign(coef))%>%
  select(feature,assoc.str,name)%>%
  pivot_wider(names_from = "feature",
              values_from = "assoc.str",
              values_fill = 0)%>%
  as.data.frame()
heatmap.df$name<-gsub("class","",heatmap.df$name)
heatmap.df<-heatmap.df%>%
  arrange(factor(name, levels=custom.levels))

rownames(heatmap.df)<-heatmap.df$name
heatmap.df<-heatmap.df%>%select(-name)

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
  filter(Taxon=="Prevotellaceae UCG-001 (Prevotellaceae)",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

ps.q.agg%>%
  filter(Taxon=="Uncultured (Eubacteriaceae)",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Uncultured (Eubacteriaceae)")

