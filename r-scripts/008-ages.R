library(eeptools)



trulvl<-"265-198"
trunclvls<-c("313-229",
             "284-229",
             "284-203",
             "273-203",
             "273-198",
             "265-203",
             "265-198")

load("./rdafiles/biagi-phyloseq-workspace.RData")
truncationlvl<-trulvl

load(paste0("./rdafiles/yasuda-",truncationlvl,"-phyloseq-workspace.RData"))

source("./r-scripts/007-prepare-barplot.R")

set.seed(1)
plot.cols<-createPalette(nrow(all.legends),
                         seedcolors = rainbow(7))

col.vec<-setNames(plot.cols,all.legends$Taxon)
pretty.facet.labels<-c("NMR"= "*Heterocephalus glaber*", # better labels for facets
                       "control"= "SPF mouse,<br>B6",
                       "NMRwt"="*Heterocephalus glaber* wild-type")

ps.total<-ps.total%>%
  bind_rows(ps.b.total)


nmr.age.groups<-nmr.agg.rel%>%
  group_by(Sample)%>%
  mutate(birthday=as.Date(birthday))%>%
  mutate(age=age_calc(birthday,units = "years"))%>%
  mutate(age=round(age))%>%
  mutate(age_group = cut(age, breaks = c(0,5,10)))

age.barplot<-ggplot(nmr.age.groups,aes(x=Sample, y=Abundance,
                       fill=factor(Taxon,
                                   levels=nmr.legend$Taxon)))+
  geom_bar(stat = "identity")+
  facet_grid(~age_group, # separate nmr from mice
             scales="free",  # free scales to separate bars according to species
             space = "free"#, # bars will have same widths
             # labeller = labeller(class=pretty.facet.labels)
  )+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec)+
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different naked mole-rat age groups"))+
  theme(
    axis.line = element_blank(), 
    strip.text.x = ggtext::element_markdown(size = 20),
    
    axis.text.x = element_text(angle=0,size=20),
    axis.text.y = element_text(size=20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25),
    plot.caption = element_text(size=23),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.position = "right",
    panel.spacing = unit(2, "lines")) # increase space between facets

ggsave(paste0("./images/barplots/",Sys.Date(),
              "-age-barplot-",truncationlvl,"-",
              agglom.rank,".png"),plot = age.barplot,
       width = 8000,height = 4500,
       dpi = 300,units = "px",device = "png")
