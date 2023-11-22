library(phyloseq)
library(vegan)
library(tidyverse)
library(patchwork)
library(Polychrome)
authorname<-"pooled"
agglom.rank<-"Genus"
truncationlvl<-"234"
read.end.type<-"single"

load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))

pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "B6mouse" = "B6 mouse",
                       "MSMmouse" = "MSM/Ms mouse",
                       "FVBNmouse" = "FVB/N <br>mouse",
                       "DMR" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus <br>cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pvo" = "*Pteromys volans orii*",
                       "NMRwt"="Wild *Heterocephalus glaber*"
)
# merge metadata
custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)

# "Clean" column: Strip families from "Unclassified" ####
# Order the data frame by the higher taxonomic rank 
# (which is the "Clean" column)
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}
ps.q.agg$Clean<-
  gsub("Unclassified \\(|Uncultured \\(", "", ps.q.agg[[agglom.rank.col-1]])
ps.q.agg$Clean<-
  gsub("\\)", "", ps.q.agg$Clean)

# Convert taxa with abundance<1% into Remainder ####
ps.q.agg$Clean<-
  ifelse(ps.q.agg$MeanRelativeAbundance<1,
         "Remainder (Mean abundance < 1%)",
         ps.q.agg$Clean)

# Separate Unclassified taxa from the rest and sort by agglom.rank-1 (higher rank) ####
# unclas is agg.rel.pooled dataset with Unclassified taxa
ps.q.agg.unclas<-
  ps.q.agg[grep("Unclassified|Uncultured", 
                          ps.q.agg[[agglom.rank.col]]),]
# but it doesn't have Remainder taxa 
ps.q.agg.unclas<-
  ps.q.agg.unclas[!grepl("Remainder", ps.q.agg.unclas$Clean),]

# clas is agg.rel.pooled dataset without unclassified taxa
ps.q.agg.clas<-
  ps.q.agg[!grepl("Unclassified|Uncultured", 
                            ps.q.agg[[agglom.rank.col]]),]
# But no Remainders
ps.q.agg.clas<-
  ps.q.agg.clas[!grepl("Remainder", ps.q.agg.clas$Clean),]


# Remainders ####
# only remainder taxa
ps.q.agg.rem<-
  ps.q.agg[grep("Remainder", ps.q.agg$Clean),]
# Taxon.bp is for the barplot
ps.q.agg.rem$Taxon<-ps.q.agg.rem$Clean


## Order by the Clean column ####
ps.q.agg.unclas<-
  ps.q.agg.unclas[order(ps.q.agg.unclas$Clean),]

ps.q.agg.clas<-
  ps.q.agg.clas[order(ps.q.agg.clas$Clean),]


# Merge them back
ps.q.agg<-rbind(ps.q.agg.rem,ps.q.agg.unclas,
                           ps.q.agg.clas)

# Create a common legend ####
ps.q.legend.unclas<-ps.q.agg.unclas
ps.q.legend.unclas<-ps.q.legend.unclas[!duplicated(ps.q.legend.unclas$Taxon),]
ps.q.legend.unclas<-ps.q.legend.unclas[order(ps.q.legend.unclas$Clean),]

ps.q.legend.clas<-ps.q.agg.clas
ps.q.legend.clas<-ps.q.legend.clas[!duplicated(ps.q.legend.clas$Taxon),]
ps.q.legend.clas<-ps.q.legend.clas[order(ps.q.legend.clas$Clean),]

ps.q.legend.rem<-ps.q.agg.clas[1,]
ps.q.legend.rem[which(sapply(ps.q.legend.rem,is.character))]<-"Remainder (Mean abundance < 1%)"
ps.q.legend.rem[which(sapply(ps.q.legend.rem,is.numeric))]<-0

ps.q.legend<-rbind(ps.q.legend.rem,ps.q.legend.unclas,
                   ps.q.legend.clas)
ps.q.legend<-ps.q.legend%>%
  select(Taxon,Taxon.bp,Clean)

# Plot the barplots ####
set.seed(1)
plot.cols<-createPalette(nrow(ps.q.legend),
                         seedcolors = rainbow(7))
# plot.cols<-createPalette(nrow(ps.q.legend),
#                          seedcolors = c("#B22222", "#0000FF","#006400", 
#                                         "#FF8C00", "#68228B"))
col.vec<-setNames(plot.cols,ps.q.legend$Taxon)


mainplot<-ps.q.agg%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ", TotalAbundance, ")"))%>% # add a column where sample names are together with sample sizes
  # filter(class=="NMR",Abundance!=0)%>%
  ggplot(aes(x=NewSample, y=RelativeAbundance,  
             fill=factor(Taxon.bp, levels=ps.q.legend$Taxon)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  
  # geom_text(mapping=aes(label=paste(Sample,Abundance.y),y=-1),
  #           size=3,
  #           angle=45,
  #           hjust=1)+
  
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
                )+
  # guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec)+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon",
       caption="Mean Relative Abundance was calculated for each host separately")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different rodents"))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        
        axis.text.x = element_text(angle=45,size=20,hjust=1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25),
        plot.caption = element_text(size=23),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "bottom")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),paste(custom.levels,collapse = '-'),
                    "barplot",truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=mainplot,
       width = 13500,height = 5200,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),paste(custom.levels,collapse = '-'),
                    "barplot",truncationlvl,
                    agglom.rank,sep = "-"),".tiff"),
       plot=mainplot,
       width = 13500,height = 5200,
       units = "px",dpi=300,device = "tiff")


# library(gridExtra)
# # TODO: arrange a grid
# ps.q.agg.rel<-ps.q.agg.rel%>%
#   left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
#   # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
#   mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
#   mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")")) # add a column where sample names are together with sample sizes
# 
# classes<-custom.levels
# classes<-c("NMR","B6mouse","hare","rabbit")
# plots<-list()
# for (i in 1:length(classes)) {
#   # Subset data for each class
#   df_sub <- subset(ps.q.agg.rel, class == classes[i])
#   
#   # Create plot for each class
#   p <- ggplot(df_sub, aes(x=NewSample, y=Abundance,  
#                           fill=factor(Taxon.bp, levels=ps.q.legend$Taxon))) +
#     geom_bar(stat = "identity") +
#     ggtitle(paste0("Class ", classes[i]))+
#     scale_fill_manual(values = col.vec)+
#     theme(
#       axis.text.x = element_text(angle=45,hjust=1),
#       axis.line = element_blank()) +
#  
#     coord_cartesian(expand=FALSE) +
#     xlab("") +
#     ylab("Mean Relative Abundance (%)")+
#     
#     theme_bw()+
#     theme(
#       axis.text.x = element_text(angle=45,hjust=1),
#       axis.line = element_blank()) +
#     # ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats and SPF mice"))+
#     theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
#           axis.line = element_blank(), 
#           strip.text.x = ggtext::element_markdown(size = 20),
#           
#           axis.text.x = element_text(angle=45,size=20,hjust=1),
#           axis.text.y = element_text(size=20),
#           axis.title = element_text(size = 20),
#           plot.title = element_text(size = 25),
#           plot.caption = element_text(size=23))+
#     guides(fill=guide_legend(ncol=1))+
#     theme(legend.position = "right")
#   # Add plot to list of plots
#   plots[[i]] <- p
# }
# plots[[length(plots)]] <- plots[[length(plots)]] + 
#   guides(fill=guide_legend(ncol=1))+
#   theme(legend.position = "right",
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 25))
# 
# # Arrange plots in a grid
# grid.arrange(grobs = plots, ncol = 4)
# 
# plots[[2]]

# Separate barplots for each host
for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.agg%>%
    left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
    filter(class==custom.levels[i],Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[i])
  lvl.name<-gsub("<br>"," ", lvl.name)
  # the total legend is big, we need to narrow down to our host
  host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))]
  
  lvl.plot<-lvl.df%>%ggplot(aes(x=NewSample, y=RelativeAbundance,  
                      fill=factor(Taxon.bp, levels=host.legend)))+
    # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
    geom_bar(stat = "identity")+ # barplot
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    # facet_grid(~class, # separate species
    #            scales="free",  # each species will have its own bars inside facet (instead of all bars)
    #            space = "free", # bars will have same widths
    #            labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
    # )+
    guides(fill=guide_legend(ncol=1))+ # legend as one column
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values = col.vec,breaks = names(col.vec))+
    xlab("") +
    ylab("Relative Abundance (%)")+
    labs(fill="Taxon")+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle=45,hjust=1),
      axis.line = element_blank()) +
    ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", lvl.name))+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          axis.line = element_blank(), 
          strip.text.x = ggtext::element_markdown(size = 20),
          
          axis.text.x = element_text(angle=45,size=20,hjust=1),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size = 20),
          plot.title = ggtext::element_markdown(size = 25),
          plot.caption = element_text(size=23),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  ggsave(paste0("./images/barplots/",
                paste(Sys.Date(),custom.levels[i],"barplot",truncationlvl,
                      agglom.rank,sep = "-"),".png"),
         plot=lvl.plot,
         width = 8000,height = 6000,
         units = "px",dpi=300,device = "png")
}

# Barplot for NMR and mice ####
lvl.df<-ps.q.agg%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
  filter(class%in%c("NMR","B6mouse"),Abundance!=0)
lvl.name<-unname(pretty.facet.labels[names(pretty.facet.labels)%in%c("NMR","B6mouse")])
lvl.name<-gsub("<br>"," ", lvl.name)
# the total legend is big, we need to narrow down to our host
host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))]

lvl.plot<-lvl.df%>%ggplot(aes(x=NewSample, y=RelativeAbundance,  
                              fill=factor(Taxon.bp, levels=host.legend)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
  )+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,breaks = names(col.vec))+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", paste(lvl.name,collapse = " and ")))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        plot.title = ggtext::element_markdown(size = 25),
        plot.caption = element_text(size=23),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")

ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"NMR-B6mouse","barplot",truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=lvl.plot,
       width = 8000,height = 6000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"NMR-B6mouse","barplot",truncationlvl,
                    agglom.rank,sep = "-"),".tiff"),
       plot=lvl.plot,
       width = 8000,height = 6000,
       units = "px",dpi=300,device = "tiff")



# Barplot for NMR and NMR wt ####
lvl.df<-ps.q.agg%>%
  left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
  # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
  mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
  mutate(NewSample=paste0(Sample," (n = ",TotalAbundance, ")"))%>%
  filter(class%in%c("NMR","NMRwt"),Abundance!=0)
lvl.name<-unname(pretty.facet.labels[names(pretty.facet.labels)%in%c("NMR","NMRwt")])
lvl.name<-gsub("<br>"," ", lvl.name)
# the total legend is big, we need to narrow down to our host
host.legend<-ps.q.legend$Taxon[ps.q.legend$Taxon%in%names(table(lvl.df$Taxon.bp))]

lvl.plot<-lvl.df%>%ggplot(aes(x=NewSample, y=RelativeAbundance,  
                              fill=factor(Taxon.bp, levels=host.legend)))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  facet_grid(~class, # separate species
             scales="free",  # each species will have its own bars inside facet (instead of all bars)
             space = "free", # bars will have same widths
             labeller = labeller(class=pretty.facet.labels) # labeller will change facet labels to custom
  )+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,breaks = names(col.vec))+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=45,hjust=1),
    axis.line = element_blank()) +
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from", paste(lvl.name,collapse = " and ")))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        plot.title = ggtext::element_markdown(size = 25),
        plot.caption = element_text(size=23),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")

ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"NMR-NMRwt","barplot",truncationlvl,
                    agglom.rank,sep = "-"),".png"),
       plot=lvl.plot,
       width = 9000,height = 6000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"NMR-NMRwt","barplot",truncationlvl,
                    agglom.rank,sep = "-"),".tiff"),
       plot=lvl.plot,
       width = 9000,height = 6000,
       units = "px",dpi=300,device = "tiff")

