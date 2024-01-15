custom.md$Sample<-rownames(custom.md)
custom.md<-custom.md%>% 
  filter(class=="NMR")%>%
  group_by(Sample)%>%
  mutate(birthday=as.Date(birthday))%>%
  mutate(age=age_calc(birthday,units = "years"))%>%
  mutate(age=round(age))
min_boundary <- floor(min(custom.md$age)/5) * 5
max_boundary <- ceiling(max(custom.md$age)/5) * 5

custom.md<-custom.md%>%
  mutate(agegroup=cut(age, breaks = seq(min_boundary, max_boundary, by = 5), 
                      include.lowest = TRUE))%>%
  as.data.frame()
rownames(custom.md)<-custom.md$Sample

foo<-custom.md%>%filter(sample.id%in%c("2D10", "2D14", "O15","Y51b",'Y66b','H4','H21','H3','G18','G14','H15'))
ageplot<-foo%>%
  group_by(agegroup)%>%
  summarise(freq=n())%>%
  ggplot(aes(x=agegroup, y=freq))+
  # Taxon.bp is from ps.q.agg.rel, while Taxon is from ps.q.legend
  geom_bar(stat = "identity")+ # barplot
  # coord_cartesian(expand=FALSE) +
  xlab("Age group") +
  ylab("Number of samples")+
  theme_bw()+
  theme(strip.text.x = ggtext::element_markdown(size = 20),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(size=20),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "bottom") 
barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved
ggsave(paste0(barplot.directory,
              paste(Sys.Date(),"age-barplot",sep = "-"),".png"),
       plot=ageplot,
       width = 2000,height = 1500,
       units = "px",dpi=300,device = "png")
