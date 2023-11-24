library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)

authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"OTU"

# Import data ####
read.end.type<-"single"
rare.status<-"rare"
# rare.status<-"nonrare"
# filter.status<-"filtered"
filter.status<-"nonfiltered"

host<-"NMR"
host<-"mice"
host.class<-c("NMR"="naked mole-rat",
              "mice"="mouse")

comparison<-"age"
comparison<-"sex"
comparison<-"strain"

load(paste0("./rdafiles/",paste(authorname,read.end.type,"qiime2",
                                truncationlvl,agglom.rank,
                                "phyloseq-workspace.RData",sep = "-")))

gg.title.taxon<-ifelse(agglom.rank=="OTU","(ASV level)",
                       paste0("(",agglom.rank," level)"))

if(host=="NMR"){
  host.labels<-
    c("NMR" = "*Heterocephalus glaber*")
}else if(host=="mice"){
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}

ps.q.df.preprocessed<-ps.q.agg

if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(Abundance!=0)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=age_calc(birthday,units = "years"))%>%
    mutate(age=round(age))
  
  min_boundary <- floor(min(ps.q.df.preprocessed$age)/5) * 5
  max_boundary <- ceiling(max(ps.q.df.preprocessed$age)/5) * 5
  
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(age_group=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                         include.lowest = TRUE))
}else if(host=="mice"){
  # select mice and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(Abundance!=0)
  ps.q.df.preprocessed$age_group<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
}




if (comparison=="age"){
  pretty.axis.labels<-names(table(ps.q.df.preprocessed$age_group))
  names(pretty.axis.labels)<-names(table(ps.q.df.preprocessed$age_group))
  custom.levels<-names(pretty.axis.labels)
  gg.labs.name<-"Age group"
  gg.title.groups<-"age groups"
  
}else if (comparison=="sex"){
  pretty.axis.labels<-
    c("F" = "Females",
      "M" = "Males")
  custom.levels<-names(pretty.axis.labels)
  pretty.axis.labels<-pretty.axis.labels[which(names(pretty.axis.labels)%in%custom.levels)]
  gg.labs.name<-"Host sex"
  gg.title.groups<-"groups"
}else if(comparison=="strain"){
  pretty.axis.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.axis.labels),custom.md$class)
  pretty.axis.labels<-pretty.axis.labels[which(names(pretty.axis.labels)%in%custom.levels)]
  gg.labs.name<-"Strain"
  gg.title.groups<-"strains"
}

scale.color.labels<-unname(pretty.axis.labels)
scale.color.breaks<-unname(pretty.axis.labels)

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

# load the list of significant features
load(file=paste0("./rdafiles/",
                 paste("significant-features-all-groups",
                       rare.status,filter.status,host,agglom.rank,
                       comparison,truncationlvl,
                       sep="-"),".RData"))
ggplots<-list()
for(ref.level in custom.levels){
  ref.level.plots<-list()
  for (feature_index in seq_along(signif.all.groups[[ref.level]])){
    # take each taxon one by one from the list of significant features
    # make a boxplot of abundances
    feature.plot<-ps.q.agg%>%
      filter(OTU==signif.all.groups[[ref.level]][[feature_index]],
             class%in%custom.levels)%>%
      ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
      geom_boxplot(show.legend = FALSE)+
      # scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      ggtitle(paste0(signif.all.groups[[ref.level]][[feature_index]],"\n(",
                    as.character(ps.q.agg[ps.q.agg$OTU==signif.all.groups[[ref.level]][[feature_index]],
                                          "Taxon"][1,1]),")\n
                    Mean relative abundance: ", round(ps.q.agg$MeanRelativeAbundance[ps.q.agg$OTU==signif.all.groups[[ref.level]][[feature_index]]][1],2), "%"))+
      theme_bw()+
      xlab("")+
      scale_color_manual(breaks = scale.color.breaks,
                         labels=scale.color.labels)+
      scale_x_discrete(labels=pretty.axis.labels,
                       limits=custom.levels)+ # rename boxplot labels (x axis)
      scale_fill_manual(values = custom.fill)+
      theme(plot.margin=unit(c(1,1,1,2), 'cm'),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = ggtext::element_markdown(hjust=1,size=10),
            axis.text.y = element_text(size=10),
            axis.title = element_text(size = 10),
            strip.text.x = element_text(size=10),
            plot.title = element_text(size = 15),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 15),
            legend.position = "right")
    ref.level.plots[[signif.all.groups[[ref.level]][[feature_index]]]]<-feature.plot
    ggsave(paste0("./images/taxaboxplots/",
                  paste(Sys.Date(),"taxaboxplot",ref.level,signif.all.groups[[ref.level]][[feature_index]],
                        rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,sep="-"),".png"),
           plot=feature.plot,
           width = 2000,height = 1500,
           units = "px",dpi=300,device = "png")
    
    ggsave(paste0("./images/taxaboxplots/",
                  paste(Sys.Date(),"taxaboxplot",ref.level,signif.all.groups[[ref.level]][[feature_index]],
                        rare.status,filter.status,host,agglom.rank,
                        comparison,truncationlvl,sep="-"),".tiff"),
           plot=feature.plot,
           width = 2000,height = 1500,
           units = "px",dpi=300,device = "tiff")
  }
  ggplots[[ref.level]]<-ref.level.plots
}
ggplots

for(i in seq_along(ggplots)){
  tax.plot<-ggplots[i]
  ggsave(paste0("./images/taxaboxplots/",
                paste(Sys.Date(),"taxaboxplot",ref.level,names(ggplots[i]),
                      rare.status,filter.status,host,agglom.rank,
                      comparison,truncationlvl,sep="-"),".png"),
         plot=tax.plot[[1]],
         width = 2000,height = 1500,
         units = "px",dpi=300,device = "png")
  
  ggsave(paste0("./images/taxaboxplots/",
                paste(Sys.Date(),"taxaboxplot",ref.level,names(ggplots[i]),
                      rare.status,filter.status,host,agglom.rank,
                      comparison,truncationlvl,sep="-"),".tiff"),
         plot=tax.plot[[1]],
         width = 2000,height = 1500,
         units = "px",dpi=300,device = "tiff")
}

signif.all.groups.df<-stack(signif.all.groups)
colnames(signif.all.groups.df)<-c("OTU","class")

signif.all.groups.df<-signif.all.groups.df%>%
  left_join(ps.q.agg,by="OTU")%>%
  dplyr::select(OTU,Taxon)%>%
  distinct()%>%
  left_join(signif.all.groups.df,by="OTU")%>%
  group_by(class,Taxon)%>%
  summarize(n=n())

signif.legend<-signif.all.groups.df%>%
  left_join(ps.q.agg)%>%
  dplyr::select(class,Taxon,n,Family)%>%
  distinct()%>%
  arrange(Family)

set.seed(1)
plot.cols<-createPalette(length(unique(signif.legend$Taxon)),
                         seedcolors = rainbow(7))
col.vec<-setNames(plot.cols,unique(signif.legend$Taxon))

signif.proportions.plot<-signif.all.groups.df%>%
  ggplot(aes(x=class, y=n,  
           fill=Taxon))+
  geom_bar(stat = "identity")+ # barplot
  theme(
    axis.text.x = element_text(hjust=1),
    axis.line = element_blank()) +
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,breaks = names(col.vec))+
  xlab("") +
  ylab("Number of ASVs belonging to a Genus")+
  labs(fill="Genus (Mean Relative Abundance >= 0.5%)",
       caption = "Differentially abundant ASVs were identified by MaAsLin2 and ANCOM-BC")+
  theme_bw()+
  theme(
    axis.text.x = element_text(),
    axis.line = element_blank()) +
  ggtitle("Differentially abundant ASVs belong to several genera")+
  # scale_y_continuous(breaks=seq(0,20,1))+# my y scale
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 15),
        
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size = 10),
        plot.title = ggtext::element_markdown(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 13),
        legend.position = "right",
        plot.caption = ggtext::element_markdown(hjust = 0, size=13),
        plot.caption.position = "plot")


ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"diff-proportions-no-filter",
                    rare.status,filter.status,host,agglom.rank,
                    comparison,truncationlvl,sep="-"),".png"),
       plot=signif.proportions.plot,
       width = 3000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"diff-proportions-no-filter",
                    rare.status,filter.status,host,agglom.rank,
                    comparison,truncationlvl,sep="-"),".tiff"),
       plot=signif.proportions.plot,
       width = 3000,height = 3000,
       units = "px",dpi=300,device = "tiff")


foo<-as.data.frame(signif.all.groups$MSMmouse)
colnames(foo)<-"OTU"
foo%>%
  inner_join(ps.q.agg)%>%
  filter(class=="MSMmouse")%>%
  dplyr::select(Taxon,OTU,Abundance)%>%
  # group_by(OTU)%>%
  mutate(tot.rel.ab=Abundance/sum(Abundance))%>%
  arrange(-tot.rel.ab)%>%
  pull(tot.rel.ab)%>%
  sum

ps.q.agg%>%
  filter(class=="MSMmouse")%>%
  dplyr::select(Taxon,OTU,Abundance)%>%
  mutate(tot.rel.ab=Abundance/sum(Abundance))%>%
  arrange(-tot.rel.ab)%>%
  filter(OTU%in%foo$OTU)%>%
  pull(tot.rel.ab)%>%
  sum
