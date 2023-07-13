library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(pheatmap)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"
source("r-scripts/make_features_pretty.R")

custom.levels<-c("NMR","SPFmouse","spalax","FukomysDamarensis","rabbit")
rare.status<-"rarefied"
rare.status<-"nonrarefied"
load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
load(paste0("./rdafiles/maaslin-",rare.status,"-","filtered-",
            paste(custom.levels,collapse = '-'),
            "-workspace.RData"))

maaslin.fit.df<-maaslin.fit_data$results%>%
  filter(qval<0.05)
maaslin.fit.df<-make_features_pretty(maaslin.fit.df, "feature")

maaslin.signif.decreased<-maaslin.fit.df%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by(feature)%>%
  filter(n()==4)%>% 
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

table(maaslin.signif.decreased$feature%in%ps.q.agg.rel$Taxon)
write.table(maaslin.signif.decreased,
            file=paste0("./rtables/alldir/maaslin-",
                        paste(custom.levels,collapse = '-'),"-",
                        "nmr-signif-",rare.status,".tsv"),
            row.names = F,sep = "\t")

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
  arrange(factor(name, levels=c("rabbit","SPFmouse","spalax","FukomysDamarensis")))

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
         class%in%c("NMR","SPFmouse","spalax","FukomysDamarensis","rabbit"))%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

ps.q.agg.abs%>%
  filter(Taxon=="Uncultured (Eubacteriaceae)",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Uncultured (Eubacteriaceae)")





aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)

intersect(maaslin.fit.df$feature,aldex.signif.features$Taxon)
match(aldex.signif.features$Taxon,maaslin.fit.df$feature)

ggplots<-list()
for (feature_index in seq_along(aldex.neg.effect$Taxon)){
  feature.plot<-ps.q.agg.abs%>%
    filter(Taxon==aldex.neg.effect$Taxon[[feature_index]],
           class%in%custom.levels)%>%
    ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
    geom_boxplot()+
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
    ggtitle(aldex.neg.effect$Taxon[[feature_index]])
  ggplots[[aldex.neg.effect$Taxon[[feature_index]]]]<-feature.plot
  print(feature.plot)
}

# magic with dplyr
maaslin.decreased.foo<-maaslin.signif.decreased%>%
  arrange(coef)%>%
  dplyr::select(feature)%>%
  distinct()%>%
  head()
subcols<-c("OTU",
           "Kingdom",          
           "Phylum",         
           "Class",           
           "Order",          
           "Family",          
           "Genus",          
           "Taxon")
foo.ps.q.agg.abs<-ps.q.agg.abs

# for (feature_index in seq_along(aldex.neg.effect$Taxon)){
for (feature_index in seq_along(maaslin.decreased.foo$feature)){
  existing.classes<-ps.q.agg.abs%>%
    filter(Taxon==maaslin.decreased.foo$feature[[feature_index]],
           class%in%custom.levels)%>%
    pull(class) %>%
    unique()
  missing.classes<-setdiff(custom.levels,existing.classes)
  if (length(missing.classes)!=0){
    existing.dummy<-ps.q.agg.abs[ps.q.agg.abs$Taxon==maaslin.decreased.foo$feature[feature_index],][1,]
    for (class_index in seq_along(missing.classes)){
      
      dummy.row<-ps.q.agg.abs[ps.q.agg.abs$class==missing.classes[[class_index]],][1,]
      dummy.row[subcols]<-existing.dummy[subcols]
      dummy.row$Abundance<-0
      
      foo.ps.q.agg.abs<-rbind(foo.ps.q.agg.abs,dummy.row)
    }
  }
  
}

foo.ps.q.agg.abs%>%
  filter(Taxon %in% maaslin.decreased.foo$feature,
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  facet_wrap(Taxon~.,scales = "free",ncol=2)
  # ggtitle(maaslin.decreased.foo$feature[[feature_index]])




##########
foo.ps.q.df.maaslin.input<-ps.q.df.maaslin.input
subcols<-c("Family",          
           "Genus",          
           "Taxon")
# for (feature_index in seq_along(aldex.neg.effect$Taxon)){
for (feature_index in seq_along(maaslin.decreased.foo$feature)){
  existing.classes<-ps.q.df.maaslin.input%>%
    filter(Taxon==maaslin.decreased.foo$feature[[feature_index]],
           class%in%custom.levels)%>%
    pull(class) %>%
    unique()
  missing.classes<-setdiff(custom.levels,existing.classes)
  if (length(missing.classes)!=0){
    existing.dummy<-ps.q.df.maaslin.input[ps.q.df.maaslin.input$Taxon==maaslin.decreased.foo$feature[feature_index],][1,]
    for (class_index in seq_along(missing.classes)){
      
      dummy.row<-ps.q.df.maaslin.input[ps.q.df.maaslin.input$class==missing.classes[[class_index]],][1,]
      dummy.row[subcols]<-existing.dummy[subcols]
      dummy.row$Abundance<-0
      dummy.row$observed_samples<-0
      dummy.row$PercentageSamples<-0
      dummy.row$MeanAbundance<-0
      foo.ps.q.df.maaslin.input<-rbind(foo.ps.q.df.maaslin.input,dummy.row)
    }
  }
  
}


foo.ps.q.df.maaslin.input%>%
  filter(Taxon %in% maaslin.decreased.foo$feature,
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  facet_wrap(Taxon~.,scales = "free",ncol=2)
# ggtitle(maaslin.decreased.foo$feature[[feature_index]])


