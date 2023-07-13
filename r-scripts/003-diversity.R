library(vegan)
library(tidyverse)
library(phyloseq)

qzaphyloseq.script.yasuda<-"./r-scripts/001-qza-to-phyloseq-yasuda.R"
qzaphyloseq.script.biagi<-"./r-scripts/001-qza-to-phyloseq-biagi.R"
source(qzaphyloseq.script.yasuda)
source(qzaphyloseq.script.biagi)
set.seed(1)



ps.df <-ps.agg.abs%>%
  select(Sample,OTU,Abundance,class,taxa.full)

ps.df <-ps.agg.abs%>%
  select(Sample,OTU,Abundance,class)
ps.df<-ps.df[!grepl("Unclassified|Uncultured",ps.df$taxa.full),]

# Alpha diversity ####
## Compute alpha diversity metrics ####
all.div<-ps.df%>%
  group_by(Sample)%>%
  summarize(sobs=specnumber(Abundance), # richness (num of species)
            shannon=diversity(Abundance,index = "shannon"),
            simpson=diversity(Abundance, index="simpson"),
            invsimpson=1/simpson, # inverse simpson
            tot=sum(Abundance),
            class=class)%>%
  pivot_longer(cols=c(sobs,shannon,simpson,invsimpson),
               names_to="metric")

## Alpha diversity tests ####
div.indices<-c("sobs","shannon", "simpson", "invsimpson") 
kt.results<-data.frame(matrix(nrow = 2,ncol=4))
w.results<-data.frame(matrix(nrow = 1,ncol=4))

colnames(kt.results)<-div.indices
colnames(w.results)<-div.indices

for (div.metric in div.indices) {
  metric.ds<-all.div%>%
    filter(metric==div.metric)%>%
    distinct()
  kt<-kruskal.test(value~class,data=metric.ds)
  w.test<-pairwise.wilcox.test(metric.ds$value,
                               metric.ds$class,
                               p.adjust.method = "BH")
  kt.results[div.metric]<-t(as.data.frame(lapply(kt[c("statistic","p.value")],as.vector)))
  w.results[div.metric]<-as.vector(w.test$p.value)
}
kt.results
w.results
stopifnot(all(kt.results[2,]<0.05))

## Plot alpha diversity metrics ####
pretty.axis.labels<-c("NMR"= "*Heterocephalus glaber*", # better labels for facets
                       "control"= "SPF mouse, B6")

metric.labs=c('sobs'="Richness \n(Observed species)",
              'shannon' = "Shannon",
              'simpson' = "Simpson",
              'invsimpson' = "Inverse \nSimpson")

div.plot<-ggplot(all.div[all.div$metric %in%
                           c("shannon","sobs",
                             "simpson","invsimpson"),],
                 aes(x=reorder(class,-value),y=value,fill=class))+
  geom_boxplot(show.legend = FALSE)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=c('sobs','shannon','simpson','invsimpson')),
             ncol=4,
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+
  scale_color_manual(breaks = c("Heterocephalus glaber",
                                "SPF mouse, B6"),
                     labels=c("Heterocephalus glaber",
                              "SPF mouse, B6"))+
  scale_x_discrete(labels=pretty.axis.labels)+ # rename boxplot labels (x axis)
  scale_fill_manual(values = c("indianred","seagreen"))+
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = ggtext::element_markdown(angle=45,hjust=1,size=18),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")+
  ggtitle(paste0("Alpha diversity of the gut microbiota of naked mole-rats and SPF mice \n(",agglom.rank, " level)"))


ggsave(paste0("./images/diversity/",Sys.Date(),
              "-alpha-shannon-sobs-",agglom.rank,"-",truncationlvl,
              ".png"),plot=div.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")

ggsave(paste0("./images/diversity/",Sys.Date(),
              "-alpha-shannon-sobs-",agglom.rank,"-",truncationlvl,
              ".tiff"),plot=div.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "tiff")



# rarefy without replacement to the min sample size
set.seed(1)
ps.agg.abs.rarefied = rarefy_even_depth(ps.agg.abs, 
                                        rngseed=1, 
                                        sample.size=min(sample_sums(ps.agg.abs)), 
                                        replace=F)
ps
ps.agg.abs.rarefied
plot_bar(ps.agg.abs.rarefied, fill="Genus")+
  facet_grid(~animal, scales="free",space="free")+  
  guides(fill=guide_legend(ncol=1))

# PCoA plot using the unweighted UniFrac as distance
set.seed(1)
wunifrac_dist = phyloseq::distance(ps.agg.abs.rarefied, 
                                   method="unifrac", 
                                   weighted=F)
jaccard.dist<- phyloseq::distance(ps.agg.abs.rarefied,
                               method = "jaccard",
                               binary=TRUE)

ordination = ordinate(ps.agg.abs.rarefied, 
                      method="PCoA", 
                      distance=jaccard.dist)
plot_ordination(ps.agg.abs.rarefied, 
                ordination, 
                color="animal") + 
  theme(aspect.ratio=1)

# Test whether the seasons differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:
adonis2(jaccard.dist ~ sample_data(ps.agg.abs.rarefied)$animal)
