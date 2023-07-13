
truncationlvl<-"234"
agglom.rank<-"Genus"
# load("./rdafiles/yashliu-dada2-284-203-Genus-138-1-phyloseq-workspace.RData")

load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "SPFmouse" = "SPF mouse, B6",
                       "FukomysDamarensis" = "*Fukomys Damarensis*",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*",
                       "spalax" = "*Nannospalax leucodon*",
                       "pal" = "*Petaurista alborufus lena*",
                       "pvo" = "*Pteromys volans orii*",
                       "ppg" = "*Petaurista philippensis grandis*",
                       "tx" = "*Trogopterus xanthipes*"
                       
                       # ,
                       # "rabbitcontrol"="*Oryctolagus cuniculus*",
                       # "harecontrol" = "*Lepus europaeus*",
                       # "ntccontrol" = "Non-treatment<br>control"
)
custom.levels<-names(pretty.facet.labels)

detach("package:VennDiagram", unload = TRUE)
library(ggVennDiagram)
grid.newpage()
# create Venn diagram with four sets
draw.quad.venn(area1=72, area2=86, area3=50, 
               area4 =52, n12=44, n23=38, n13=27, 
               n14= 32,n24=32, n34=20, n123=18, 
               n124=17, n234=13, n134=11, n1234=6, 
               category=c("Cricket","Football","Badminton","Table Tennis"),
               col="Green",fill=c("Red","Pink","Blue","Orange"),lty="dashed")
draw.quintuple.venn(area1=72, 
                    area2=86, 
                    area3=50, 
                    area4 =52,
                    area5 = 30,
                    )



pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
              "SPFmouse" = "SPF mouse, B6",
              "FukomysDamarensis" = "*Fukomys Damarensis*",
              "spalax" = "*Nannospalax leucodon*",
              "rabbit" = "*Oryctolagus cuniculus*")

pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "SPFmouse" = "SPF mouse, B6",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*")
pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                       "pal" = "*Petaurista alborufus lena*",
                       "pvo" = "*Pteromys volans orii*",
                       "ppg" = "*Petaurista philippensis grandis*",
                       "tx" = "*Trogopterus xanthipes*")

custom.levels<-names(pretty.facet.labels)

ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)%>%
  filter(Abundance!=0)

# Filtering by prevalence ####
# for each host, calculate the number of samples an ASV was observed in 
ps.q.sample_counts <- ps.q.agg.abs%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class,Family,Genus) %>% 
  summarize(observed_samples = n_distinct(Sample))

# calculate the total number of samples for each host
ps.q.host_counts <- ps.q.agg.abs%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class)%>%
  summarize(total_samples = n_distinct(Sample)) %>%
  ungroup()

ps.q.prevalences <- ps.q.host_counts %>%
  left_join(ps.q.sample_counts, by = "class") %>%
  mutate(PercentageSamples = observed_samples/total_samples*100)%>% # percentage
  # of samples an ASV was observed in
  left_join(ps.q.agg.rel.pooled,by=c("class","Genus","Family"))%>% # add mean
  # relative abundances
  rename(MeanAbundance=Abundance)%>%
  left_join(unique(ps.q.agg.abs[,c("Genus","Family","Taxon","class")]),
            by=c("class","Genus","Family")) # add the Taxon column

## filter by percentage of samples an ASV was observed in ####
ps.q.df.norare.filtered <-ps.q.df%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"))%>%
  filter(PercentageSamples>=10)%>%
  filter(MeanAbundance>=1)

taxa.list<-list()
for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.agg.rel%>%
    left_join(.,ps.q.total,by="Sample",suffix=c("",".y"))%>% # merge our data with the dataset of total abundances, so we can have info about sample size
    # mutate(class=factor(class,levels=c("NMRwt","NMR","control")))%>% # change the order of our class column, so the NMR will be first
    mutate(class=factor(class,levels=custom.levels))%>% # change the order of our class column, so the NMR will be first
    mutate(NewSample=paste0(Sample," (n = ",Abundance.y, ")"))%>%
    filter(class==custom.levels[i],Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[i])
  lvl.name<-gsub("<br>"," ", lvl.name)
  taxa.list[[custom.levels[i]]]<-unique(lvl.df$Taxon)
}

for(i in seq_along(custom.levels)){
  lvl.df<-ps.q.df.norare.filtered%>%
    filter(class==custom.levels[i],Abundance!=0)
  lvl.name<-unname(pretty.facet.labels[i])
  lvl.name<-gsub("<br>"," ", lvl.name)
  taxa.list[[custom.levels[i]]]<-unique(lvl.df$Taxon)
}

vennD<-ggVennDiagram(taxa.list, 
                     label_alpha = 0,
                     category.names = pretty.facet.labels, 
                     label_geom = "label", 
                     label_size = 4) +
  scale_fill_distiller(palette = "Reds", direction = 1)+ 
  scale_x_continuous(expand = expansion(mult = .2))

ggsave(paste0("./images/venn/",Sys.Date(),
              "-venn-",paste(names(taxa.list),collapse = "_"),"-",
              agglom.rank,"-",truncationlvl,".png"),plot = vennD,
       width = 2000,height = 2500,
       dpi = 300,units = "px",device = "png") 
