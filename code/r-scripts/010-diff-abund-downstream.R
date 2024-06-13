library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)
# Import data ####
ref.level="NMR"
truncationlvl<-"234"
agglom.rank<-"Genus"
source("./code/r-scripts/make_features_maaslin.R")
phyloseq.workspace.date_time<-"20240426_21_44_30"
maaslin.date_time<-"20240427_16_22_09"
aldex2.date_time<-"20240427_17_39_41"
ancombc.date_time<-"20240427_17_12_27"
custom.levels<-c("NMR",
                 "B6mouse",
                 "MSMmouse",
                 "FVBNmouse",
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
# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
load(file.path("./output/rdafiles",paste(
  phyloseq.workspace.date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)

maaslin.signif.features<-
  read.table(file.path("./output/rtables",authorname,
                        paste(maaslin.date_time,"maaslin",rare.status,
                              filter.status,agglom.rank,
                              paste(custom.levels,collapse = '-'),
                              truncationlvl,ref.level,
                              "signif.tsv",sep="-")),
             header = T,sep = "\t")
aldex.signif.features<-
  read.table(file.path("./output/rtables",authorname,
                    paste(aldex2.date_time,"aldex2",rare.status,
                          filter.status,agglom.rank,
                          paste(custom.levels,collapse = '-'),
                          truncationlvl,ref.level,
                          "signif.tsv",sep="-")),
             header = T,sep = "\t")
ancombc.signif.features<-
  read.table(file.path("./output/rtables",authorname,
                       paste(ancombc.date_time,"ancombc",rare.status,
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
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "maaslin.signif.decreased",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(aldex.neg.effect,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "aldex.neg.effect",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ancombc.signif.decreased,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "ancombc.signif.decreased",rare.status,
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
# Only Maaslin2 and ANCOM-BC: they're decreased in other hosts, not in ref
common.decreased<-Reduce(intersect,list(maaslin.signif.decreased$feature,
                                     ancombc.signif.decreased$taxon_id))

write.table(common.decreased,
            file=file.path("./output/rtables",authorname,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "significant-features",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            row.names = F,sep = "\t",col.names = F)


maaslin.signif.features<-c("Desulfovibrio","Campylobacter")
for (i in seq(maaslin.signif.features)){
  taxon.plot<-ps.q.agg%>%
    filter(class=="NMR")%>%
    left_join(custom.md,by="Sample")%>%
    filter(Genus==maaslin.signif.features[i])%>%
    group_by(Sample,old_agegroup,Genus)%>%
    summarise(RelativeAbundance=sum(RelativeAbundance))%>%
    ggplot(aes(x=Sample,y=RelativeAbundance))+
    geom_bar(stat="identity")+
    ylab("Relative Abundance (%)")+
    
    facet_grid(~old_agegroup,scales = "free_x")+
    ggtitle(paste(maaslin.signif.features[i],"relative abundance"))+
    theme_bw()+
    theme(axis.line = element_blank(),
          strip.text.x = ggtext::element_markdown(size = 20),# the name of
          # each facet will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          axis.text.x = element_text(size=20),# rotate
          # the x-axis labels by 45 degrees and shift to the right
          axis.text.y = element_text(size=20), # size of y axis ticks
          axis.title = element_text(size = 20), # size of axis names
          plot.title = ggtext::element_markdown(size = 25), # the plot
          # title will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          plot.caption = element_text(size=23),# size of plot caption
          legend.text = element_text(size = 20),# size of legend text
          legend.title = element_text(size = 25), # size of legend title
          legend.position = "right")
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "barplot","NMR",maaslin.signif.features[i],
                      agglom.rank,sep = "-"),".png"),
         plot=taxon.plot,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = "png")
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "barplot","NMR",maaslin.signif.features[i],
                      agglom.rank,sep = "-"),".tiff"),
         plot=taxon.plot,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = "tiff")
}

