library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)
ref.level="NMR"

truncationlvl<-"234"
agglom.rank<-"Genus"
source("r-scripts/make_features_maaslin.R")

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
rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"

authorname<-"pooled"
load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)
ps.q.total<-ps.q.total%>%
  filter(Sample%in%rownames(custom.md))
ps.q.1pc<-ps.q.1pc%>%
  filter(class%in%custom.levels)

maaslin.signif.features<-
  read.table(file.path("./rtables",authorname,
                        paste("maaslin",rare.status,
                              filter.status,agglom.rank,
                              paste(custom.levels,collapse = '-'),
                              truncationlvl,ref.level,
                              "signif.tsv",sep="-")),
             header = T,sep = "\t")
aldex.signif.features<-
  read.table(file.path("./rtables",authorname,
                    paste("aldex2",rare.status,
                          filter.status,agglom.rank,
                          paste(custom.levels,collapse = '-'),
                          truncationlvl,ref.level,
                          "signif.tsv",sep="-")),
             header = T,sep = "\t")
ancombc.signif.features<-
  read.table(file.path("./rtables",authorname,
                       paste("ancombc",rare.status,
                             filter.status,agglom.rank,
                             paste(custom.levels,collapse = '-'),
                             truncationlvl,ref.level,
                             "signif.tsv",sep="-")),
             header = T,sep = "\t")

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

write.table(maaslin.signif.decreased,
            file=file.path("./rtables",authorname,
                           paste("maaslin.signif.decreased",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(aldex.neg.effect,
            file=file.path("./rtables",authorname,
                           paste("aldex.neg.effect",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ancombc.signif.decreased,
            file=file.path("./rtables",authorname,
                           paste("ancombc.signif.decreased",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

# find common significant features between three tools
Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$Taxon,
                        ancombc.signif.features$taxon_id))


# find common significantly decreased features between three tools
Reduce(intersect,list(maaslin.signif.decreased$feature,aldex.neg.effect$Taxon,
                        ancombc.signif.decreased$taxon_id))

match(aldex.signif.features$Taxon,maaslin.signif.features$feature)

# common significant and decreased features between two tools
common.signif<-Reduce(intersect,list(maaslin.signif.features$feature,aldex.signif.features$Taxon,
                                     ancombc.signif.features$taxon_id))
# they're decreased in other hosts, not in ref
common.decreased<-Reduce(intersect,list(maaslin.signif.decreased$feature,
                                     ancombc.signif.decreased$taxon_id))

write.table(common.decreased,
            file=file.path("./rtables",authorname,
                           paste("significant-features",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t",col.names = F)
