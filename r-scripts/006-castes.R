library(ggVennDiagram)

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

nmr.pair<-ps.agg.rel %>%
  filter(class=="NMR")%>%
  filter(caste=="pair")
nmr.single<-nmr.agg.rel%>%
  filter(caste=="single")
nmr.worker<-nmr.agg.rel%>%
  filter(caste=="worker")



caste.plot<-ggplot(nmr.agg.rel,aes(x=Sample, y=Abundance,
                       fill=factor(Taxon,
                                   levels=nmr.legend$Taxon)))+
  geom_bar(stat = "identity")+
  facet_grid(~caste, # separate nmr from mice
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
  theme_bw()+ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from different naked mole-rat castes"))+
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
        panel.spacing = unit(2, "lines"))

ggsave(paste0("./images/barplots/",Sys.Date(),
              "-caste-barplot-",truncationlvl,"-",
              agglom.rank,".png"),plot = caste.plot,
       width = 8000,height = 4500,
       dpi = 300,units = "px",device = "png")


pair.genus.freqs<-as.data.frame(table(nmr.pair$Genus))
single.genus.freqs<-as.data.frame(table(nmr.single$Genus))
worker.genus.freqs<-as.data.frame(table(nmr.worker$Genus))

overlap<-calculate.overlap(x<-list("pair"=nmr.pair$Genus,
                                   "single"=nmr.single$Genus,
                                   "worker"=nmr.worker$Genus))
m40<-nmr.pair%>%
  filter(Sample=="M40")

foo<-list(
  pair=nmr.pair$Genus,
  single=nmr.single$Genus,
  worker=nmr.worker$Genus
)

pretty.facet.labels<-c("NMRwt"="*Heterocephalus glaber* \nwild-type",
                       "NMR"= "*Heterocephalus glaber*", # better labels for facets
                       "control"= "SPF mouse,<br>B6")

ggVennDiagram(foo, label_alpha = 0,
              category.names = names(foo)) +
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")
