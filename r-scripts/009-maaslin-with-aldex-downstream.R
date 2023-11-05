library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)

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
read.end.type<-"single"

host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
comparison<-"sex"
comparison<-"strain"
# load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
#             "-phyloseq-workspace.RData"))

if(agglom.rank=="OTU"){
  load(file.path("./rdafiles",
              paste("maaslin",rare.status,filter.status,host,agglom.rank,
                    comparison,truncationlvl,
                    "workspace.RData",sep="-")))
}else{
  load(file.path("./rdafiles",
              paste("maaslin",rare.status,filter.status,agglom.rank,
                    paste(custom.levels,collapse = '-'),
                    truncationlvl,"workspace.RData",sep="-")))
}
maaslin.fit.df<-maaslin.fit_data$results%>%
  filter(qval<0.05)
if(agglom.rank=="OTU"){
  maaslin.fit.df$feature<-gsub("^X","",maaslin.fit.df$feature)
}else{ #TODO: fix substitutions dash and space
  maaslin.fit.df<-make_features_pretty(maaslin.fit.df, "feature")
}
maaslin.signif.decreased<-
  read.table(file.path("./rtables/alldir",
                    paste("maaslin",rare.status,
                          filter.status,agglom.rank,
                          paste(custom.levels,collapse = '-'),truncationlvl,
                          "nmr-signif.tsv",sep="-")),
             header = T,sep = "\t")




if(agglom.rank=="OTU"){
  maaslin.signif.decreased<-
    read.table(file.path("./rtables/alldir",
                          paste("maaslin",rare.status,
                                filter.status,host,agglom.rank,
                                comparison,truncationlvl,custom.levels[1],#pre-loaded
                                "signif.tsv",sep="-")),
               header = T,sep = "\t")
  aldex.signif.features<-
    read.table(file.path("./rtables/alldir/",
                      paste("aldex2",rare.status,
                            filter.status,host,agglom.rank,
                            comparison,truncationlvl,
                            "signif.tsv",sep="-")),
               header = T,sep = "\t")
}else{
  maaslin.signif.decreased<-
    read.table(file.path("./rtables/alldir/",
                          paste("maaslin",rare.status,
                                filter.status,agglom.rank,
                                paste(custom.levels,collapse = '-'),truncationlvl,
                                "nmr-signif.tsv",sep="-")),
               header = T,sep = "\t")
  aldex.signif.features<-
    read.table(file.path("./rtables/alldir/",
                      paste("aldex2",rare.status,
                            filter.status,agglom.rank,
                            paste(custom.levels,collapse = '-'),truncationlvl,
                            "nmr-signif.tsv",sep="-")),
               header = T,sep = "\t")
}

aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)

if(agglom.rank=="OTU"){
  intersect(maaslin.fit.df$feature,aldex.signif.features$OTU)
}else{ 
  intersect(maaslin.fit.df$feature,aldex.signif.features$Taxon)
}

if(agglom.rank=="OTU"){
  match(aldex.signif.features$OTU,maaslin.fit.df$feature)
}else{ 
  match(aldex.signif.features$Taxon,maaslin.fit.df$feature)
}

ggplots<-list()
if(agglom.rank=="OTU"){
  for (feature_index in seq_along(aldex.neg.effect$Taxon)){
    # take each taxon one by one from aldex.neg.effect
    # make a boxplot of abundances
    feature.plot<-ps.q.agg%>%
      filter(Taxon==aldex.neg.effect$Taxon[[feature_index]],
             class%in%custom.levels)%>%
      ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
      geom_boxplot()+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      ggtitle(aldex.neg.effect$Taxon[[feature_index]])
    ggplots[[aldex.neg.effect$Taxon[[feature_index]]]]<-feature.plot
    print(feature.plot)
  }
}else{ 
  for (feature_index in seq_along(aldex.neg.effect$OTU)){
    # take each otu one by one from aldex.neg.effect
    # make a boxplot of abundances
    feature.plot<-ps.q.agg%>%
      filter(OTU==aldex.neg.effect$OTU[[feature_index]],
             class%in%names(host.labels))%>%
      ggplot(aes(x=factor(class,level=host.labels),y=Abundance,fill=class))+
      geom_boxplot()+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      ggtitle(paste(aldex.neg.effect$OTU[[feature_index]],
                    as.character(ps.q.agg[ps.q.agg$OTU==aldex.neg.effect$OTU[[feature_index]],
                                          "Taxon"][1,1])))
    ggplots[[aldex.neg.effect$OTU[[feature_index]]]]<-feature.plot
    print(feature.plot)
  }
}
ggplots
