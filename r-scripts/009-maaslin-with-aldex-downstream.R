library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)

truncationlvl<-"234"
agglom.rank<-"Genus"
source("r-scripts/make_features_pretty.R")

custom.levels<-c("NMR","B6mouse","DMR","hare","rabbit",
                 "spalax","pvo")
rare.status<-"rare"
filter.status<-"nonfiltered"

read.end.type<-"single"
asvlevel=TRUE
load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
if(asvlevel==TRUE){
  load(paste0("./rdafiles/maaslin-",rare.status,"-",filter.status,"-",
              paste(custom.levels,collapse = '-'),
              "-workspace-ASVlevel.RData"))
  maaslin.fit.df<-maaslin.fit_data$results%>%
    filter(qval<0.05)
  maaslin.fit.df$feature<-gsub("^X","",maaslin.fit.df$feature)
  maaslin.signif.decreased<-
    read.table(paste0("./rtables/alldir/maaslin-",
                      paste(custom.levels,collapse = '-'),"-",
                      "nmr-signif-",rare.status,"-ASVlevel.tsv"),
               header = T,sep = "\t")
}else{
  load(paste0("./rdafiles/maaslin-",rare.status,"-",filter.status,"-",
              paste(custom.levels,collapse = '-'),
              "-workspace.RData"))
  maaslin.fit.df<-maaslin.fit_data$results%>%
    filter(qval<0.05)
  maaslin.fit.df<-make_features_pretty(maaslin.fit.df, "feature")
  maaslin.signif.decreased<-
    read.table(paste0("./rtables/alldir/maaslin-",
                      paste(custom.levels,collapse = '-'),"-",
                      "nmr-signif-",rare.status,".tsv"),
               header = T,sep = "\t")
}


if(asvlevel==TRUE){
  aldex.signif.features<-
    read.table(paste0("./rtables/alldir/aldex2-",
                      rare.status,"-",filter.status,"-",
                      paste(custom.levels,collapse = '-'),"-",
                      "nmr-signif","-ASVlevel.tsv"), # no rare.status
               header = T,sep = "\t")
}else{
  aldex.signif.features<-
    read.table(paste0("./rtables/alldir/aldex2-",
                      rare.status,"-",filter.status,"-",
                      paste(custom.levels,collapse = '-'),"-",
                      "nmr-signif",".tsv"), # no rare.status
               header = T,sep = "\t")
}



aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)

intersect(maaslin.fit.df$feature,aldex.signif.features$Taxon)
intersect(maaslin.fit.df$feature,aldex.signif.features$OTU)


match(aldex.signif.features$Taxon,maaslin.fit.df$feature)
match(aldex.signif.features$OTU,maaslin.fit.df$feature)

ggplots<-list()
for (feature_index in seq_along(aldex.neg.effect$Taxon)){
  # take each taxon one by one from aldex.neg.effect
  # make a boxplot of abundances
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

