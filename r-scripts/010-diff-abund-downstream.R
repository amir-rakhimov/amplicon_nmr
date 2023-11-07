library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)

truncationlvl<-"234"
agglom.rank<-"Genus"
agglom.rank<-"OTU"
source("r-scripts/make_features_maaslin.R")

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
load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))
ref.level="NMR"
if(agglom.rank=="OTU"){
  maaslin.signif.features<-
    read.table(file.path("./rtables/alldir",
                          paste("maaslin",rare.status,
                                filter.status,host,agglom.rank,
                                comparison,truncationlvl,ref.level,
                                "signif.tsv",sep="-")),
               header = T,sep = "\t")
  aldex.signif.features<-
    read.table(file.path("./rtables/alldir",
                      paste("aldex2",rare.status,
                            filter.status,host,agglom.rank,
                            comparison,truncationlvl,ref.level,
                            "signif.tsv",sep="-")),
               header = T,sep = "\t")
  ancombc.signif.features<-read.table(file.path("./rtables/alldir",
                      paste("ancombc",rare.status,
                            filter.status,host,agglom.rank,
                            comparison,truncationlvl,ref.level,
                            "signif.tsv",sep="-")),
               header = T,sep = "\t")
}else{
  maaslin.signif.features<-
    read.table(file.path("./rtables/alldir",
                          paste("maaslin",rare.status,
                                filter.status,agglom.rank,
                                paste(custom.levels,collapse = '-'),
                                truncationlvl,ref.level,
                                "signif.tsv",sep="-")),
               header = T,sep = "\t")
  aldex.signif.features<-
    read.table(file.path("./rtables/alldir",
                      paste("aldex2",rare.status,
                            filter.status,agglom.rank,
                            paste(custom.levels,collapse = '-'),
                            truncationlvl,ref.level,
                            "signif.tsv",sep="-")),
               header = T,sep = "\t")
  ancombc.signif.features<-
    read.table(file.path("./rtables/alldir",
                         paste("ancombc",rare.status,
                               filter.status,agglom.rank,
                               paste(custom.levels,collapse = '-'),
                               truncationlvl,ref.level,
                               "signif.tsv",sep="-")),
               header = T,sep = "\t")
  
}

maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  group_by(feature)%>%
  filter(n()==length(custom.levels)-1)%>% 
  arrange(feature)%>%
  mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))

aldex.neg.effect<-aldex.signif.features%>%
  filter(effect<0)

ancombc.signif.decreased<-subset(ancombc.signif.features,
                       rowSums(ancombc.signif.features[,-c(1,2)]<0)==ncol(ancombc.signif.features[,-c(1,2)]))

# find common significant features between three tools
if(agglom.rank=="OTU"){
  Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$OTU,
                        ancombc.signif.features$taxon_id))
}else{ 
  Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$Taxon,
                        ancombc.signif.features$taxon_id))
}

# find common significantly decreased features between three tools
if(agglom.rank=="OTU"){
  Reduce(intersect,list(maaslin.signif.decreased$feature,aldex.neg.effect$OTU,
                        ancombc.signif.decreased$taxon_id))
}else{ 
  Reduce(intersect,list(maaslin.signif.decreased$feature,aldex.neg.effect$Taxon,
                        ancombc.signif.decreased$taxon_id))
}

if(agglom.rank=="OTU"){
  match(aldex.signif.features$OTU,maaslin.signif.features$feature)
}else{ 
  match(aldex.signif.features$Taxon,maaslin.signif.features$feature)
}

# common significant and decreased features between two tools
common.signif<-Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$Taxon,
                                     ancombc.signif.features$taxon_id))
# they're decreased in other hosts, not in ref
common.decreased<-Reduce(intersect,list(maaslin.signif.decreased$feature,
                                     ancombc.signif.decreased$taxon_id))


write.table(common.decreased,
            file=file.path("./rtables/alldir",
                           paste("significant-features",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t",col.names = F)

