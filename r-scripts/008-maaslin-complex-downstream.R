library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(pheatmap)

truncationlvl<-"234"
agglom.rank<-"Genus"
agglom.rank<-"OTU"
source("r-scripts/make_features_pretty.R")

custom.levels<-c("NMR",
                 "B6mouse",
                 "MSMmouse",
                 "FVBNmouse",
                 "DMR",
                 "hare",
                 "rabbit",
                 "spalax",
                 "pvo")
rare.status<-"rare"
filter.status<-"nonfiltered"
host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

read.end.type<-"single"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

if(asvlevel==TRUE){
  load(paste0("./rdafiles/maaslin-",rare.status,"-",filter.status,"-",
              paste(custom.levels,collapse = '-'),
              "-workspace-ASVlevel.RData"))
}else{
  load(paste0("./rdafiles/maaslin-",rare.status,"-",filter.status,"-",
            paste(custom.levels,collapse = '-'),
            "-workspace.RData"))
}
paste0("./rdafiles/",
       paste("maaslin",rare.status,filter.status,
             paste(custom.levels,collapse = '-'),agglom.rank,
             truncationlvl,"workspace.RData",sep="-"))
paste0("./rdafiles/",
       paste("maaslin",rare.status,filter.status,host,agglom.rank,
             comparison,truncationlvl,
             "workspace.RData",sep="-"))


maaslin.fit.df<-maaslin.fit_data$results%>%
  filter(qval<0.05)
if(asvlevel==TRUE){
  maaslin.fit.df$feature<-gsub("^X","",maaslin.fit.df$feature)
}else{
  maaslin.fit.df<-make_features_pretty(maaslin.fit.df, "feature")
}


maaslin.signif.decreased<-maaslin.fit.df%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by(feature)%>%
  filter(n()==length(custom.levels)-1)%>% 
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

# maaslin.signif.increased<-maaslin.fit.df%>%
#   as_tibble()%>%
#   filter(coef>0)%>%
#   group_by(feature)%>%
#   filter(n()==4)%>% 
#   arrange(feature)%>%
#   mutate(n=n())%>%
#   mutate(assoc.str=-log(qval)*sign(coef))#%>%
# # select(feature,assoc.str,name)
if(asvlevel==TRUE){
  table(maaslin.signif.decreased$feature%in%ps.q.agg.abs$OTU)
}else{
  table(maaslin.signif.decreased$feature%in%ps.q.agg.rel$Taxon)
}


if(asvlevel==TRUE){
  write.table(maaslin.signif.decreased,
              file=paste0("./rtables/alldir/maaslin-",rare.status,"-",
                          filter.status,"-",
                          paste(custom.levels,collapse = '-'),"-",
                          "nmr-signif-ASVlevel.tsv"),
              row.names = F,sep = "\t")
}else{
  write.table(maaslin.signif.decreased,
              file=paste0("./rtables/alldir/maaslin-",rare.status,"-",
                          filter.status,"-",
                          paste(custom.levels,collapse = '-'),"-",
                          "nmr-signif.tsv"),
              row.names = F,sep = "\t")
}

# Heatmap ####
heatmap.df<-maaslin.fit.df%>%
  mutate(assoc.str=-log(qval)*sign(coef))%>%
  select(feature,assoc.str,name)%>%
  pivot_wider(names_from = "feature",
              values_from = "assoc.str",
              values_fill = 0)%>%
  as.data.frame()
heatmap.df$name<-gsub("class","",heatmap.df$name)
heatmap.df<-heatmap.df%>%
  arrange(factor(name, levels=c("rabbit","B6mouse","spalax","DMR")))

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

ps.q.df.maaslin.input%>%
  filter(Taxon=="Prevotellaceae_UCG-001 (Prevotellaceae)",
         class%in%c("NMR","B6mouse","spalax","DMR","rabbit"))%>%
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

