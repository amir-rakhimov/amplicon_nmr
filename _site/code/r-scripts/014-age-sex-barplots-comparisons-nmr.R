library(eeptools)
library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)
library(pals)

truncationlvl<-"234"
agglom.rank<-"Genus"
agglom.rank<-"OTU"
read.end.type<-"single"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))


set.seed(1)
plot.cols<-createPalette(nrow(all.legends),
                         seedcolors = rainbow(7))

col.vec<-setNames(plot.cols,all.legends$Taxon)
pretty.facet.labels<-c("NMR"= "*Heterocephalus glaber*", # better labels for facets
                       "control"= "SPF mouse,<br>B6",
                       "NMRwt"="*Heterocephalus glaber* wild-type")


ps.q.nmr.with.ages<-ps.q.agg.rel%>%
  filter(class=="NMR",Abundance!=0)%>%
  group_by(Sample)%>%
  mutate(birthday=as.Date(birthday))%>%
  mutate(age=age_calc(birthday,units = "years"))%>%
  mutate(age=round(age))%>%
  mutate(age_group = cut(age, breaks = c(0,5,10,15)))

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

# Comparing sexes ####
# Find unique taxa
m.taxa<-ps.q.agg.abs%>%
  filter(class=="NMR",sex=="M",Abundance!=0)%>%
  select(Taxon)%>%
  distinct()%>%
  as_vector()
f.taxa<-ps.q.agg.abs%>%
  filter(class=="NMR",sex=="F",Abundance!=0)%>%
  select(Taxon)%>%
  distinct()%>%
  as_vector()
setdiff(m.taxa,f.taxa)
setdiff(f.taxa,m.taxa)

m.asv<-ps.q.agg.abs%>%
  filter(class=="NMR",sex=="M",Abundance!=0)%>%
  select(OTU)%>%
  distinct()%>%
  as_vector()
f.asv<-ps.q.agg.abs%>%
  filter(class=="NMR",sex=="F",Abundance!=0)%>%
  select(OTU)%>%
  distinct()%>%
  as_vector()
setdiff(m.asv,f.asv)
setdiff(f.asv,m.asv)

# Comparing age groups ####
# Find unique taxa
age_groups<-names(table(ps.q.nmr.with.ages$age_group))
taxa.list<-list()
asv.list<-list()

for (i in seq_along(age_groups)){
  taxa<-ps.q.nmr.with.ages%>%ungroup%>%
    filter(age_group==age_groups[i])%>%
    select(Taxon)%>%
    distinct()%>%
    as.vector()
  asv<-ps.q.nmr.with.ages%>%ungroup%>%
    filter(age_group==age_groups[i])%>%
    select(OTU)%>%
    distinct()%>%
    as.vector()
  taxa.list[[i]]<-taxa
  asv.list[[i]]<-asv
}

combinations<-combn(c(1,2,3),2)
diff.list.taxa<-list()
diff.list.asv<-list()
for (i in 1:ncol(combinations)) {
  diff.list.taxa[[i]]<-setdiff(taxa.list[combinations[1,i]][[1]]$Taxon,taxa.list[combinations[2,i]][[1]]$Taxon)
  diff.list.asv[[i]]<-setdiff(asv.list[combinations[1,i]][[1]]$OTU,asv.list[combinations[2,i]][[1]]$OTU)
}

