library(tidyverse)
# library(phyloseq)
# library(Maaslin2)
# library(vegan)
# library(ALDEx2)
# Import data ####
ref.level="NMR"
truncationlvl<-"234"
agglom.rank<-"Genus"
source("./code/r-scripts/make_features_maaslin.R")
phyloseq.workspace.date_time<-"20240426_21_44_30"
# Tables of differentially abundant taxa
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
                        paste(maaslin.date_time,"maaslin2",rare.status,
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

# write.table(maaslin.signif.decreased,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "maaslin.signif.decreased",rare.status,
#                                  filter.status,agglom.rank,
#                                  paste(custom.levels,collapse = '-'),truncationlvl,
#                                  ref.level,"signif.tsv",sep="-")),
#             row.names = F,sep = "\t")
# 
# write.table(aldex.neg.effect,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "aldex.neg.effect",rare.status,
#                                  filter.status,agglom.rank,
#                                  paste(custom.levels,collapse = '-'),truncationlvl,
#                                  ref.level,"signif.tsv",sep="-")),
#             row.names = F,sep = "\t")
# 
# write.table(ancombc.signif.decreased,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "ancombc.signif.decreased",rare.status,
#                                  filter.status,agglom.rank,
#                                  paste(custom.levels,collapse = '-'),truncationlvl,
#                                  ref.level,"signif.tsv",sep="-")),
#             row.names = F,sep = "\t")

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

# write.table(common.decreased,
#             file=file.path("./output/rtables",authorname,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "significant-features",rare.status,
#                                  filter.status,agglom.rank,
#                                  paste(custom.levels,collapse = '-'),truncationlvl,
#                                  ref.level,"signif.tsv",sep="-")),
#             row.names = F,sep = "\t",col.names = F)



####
library(Polychrome)
library(ggtext)
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "DMR" = "*Fukomys Damarensis*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*"#,
                      # "NMRwt"="Wild *Heterocephalus glaber*"
)
set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900",
                                          "#FFA07A"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

ps.q.agg%>%
  filter(get(agglom.rank)%in%common.decreased)%>%
  group_by_at(c("class",agglom.rank))%>%
  ggplot(aes(x=factor(class,level=rev(custom.levels)),
             y=RelativeAbundance,
             fill=factor(class)))+
  geom_boxplot(show.legend = FALSE)+
  # facet_wrap(~Genus,scales = "free_x",
  #            ncol = 2)+
  theme_bw()+
  coord_flip()+
  labs(x="",
       y="Relative abundance (%)")+
  scale_color_manual(breaks = rev(unname(pretty.level.names)),
                     labels=rev(unname(pretty.level.names)))+
  scale_x_discrete(labels=rev(pretty.level.names),
                   limits=rev(custom.levels))+ # rename boxplot labels (x axis)
  scale_fill_manual(values = rev(custom.fill))+
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "none")+
  ggtitle(paste0("Relative abundance of naked mole-rat-specific taxa"))
ggsave(paste0("./images/taxaboxplots/",
              paste(paste(format(Sys.time(),format="%Y%m%d"),
                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                    "NMR-specific-bacteria",
                    sep = "-"),".","png"),
       plot=last_plot(),
       width = 4000,height = 12000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/taxaboxplots/",
              paste(paste(format(Sys.time(),format="%Y%m%d"),
                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                    "NMR-specific-bacteria",
                    sep = "-"),".","tiff"),
       plot=last_plot(),
       width = 4000,height = 12000,
       units = "px",dpi=300,device = "tiff")
