library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)
library(pals)

# "274-203"
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

# Using Polychrome package
set.seed(1)
# custom.colors<- createPalette(length(custom.levels),
#               seedcolors = rainbow(6))
custom.colors<- createPalette(length(custom.levels),
                              seedcolors = c("#B22222", "#0000FF","#006400", 
                                             "#FF8C00", "#68228B"))
custom.colors<-unname(custom.colors)
swatch(custom.colors)
# custom colors for scale_fill_manual (maybe not needed)
# custom.fill<-c("dodgerblue", "seagreen","indianred")

# custom labels for scale_color_manual
custom.color.labels<-unname(pretty.facet.labels)


set.seed(1)
permut.num<-1000 # number of permutations for PERMANOVA

ps.sampledata<-ps.q.agg.abs%>%
  select(c("Sample","class","origin","relation",
           "sex","caste","birthday"))%>%
  distinct() # metadata



ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)
# convert the data frame into wide format
all.wide<-ps.q.df%>%
  select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs and rownames are sample IDs
rownames(all.wide)<-all.wide$Sample
all.wide<-all.wide[,-c(1,2)]  
all.wide<-as.matrix(all.wide)


# Bray ####
# Meaningful only for integers (counts)
set.seed(1)
# computes the dissimilarity matrix of a dataset MULTIPLE times using vegdist
# while randomly subsampling the dataset each time.
# All of the subsampled iterations are then averaged (mean) to provide a 
# distance matrix that represents the average of multiple subsampling iterations.
b.dist<-avgdist(all.wide,
                dmethod="bray",
                sample=1000,
                iterations = 1000)
set.seed(1)



b.dist.df<-b.dist%>% # convert distances into tibble
  as.matrix()%>%
  as_tibble()
b.dist.df<-cbind(labels(b.dist),b.dist.df) # create a column of sample ids
colnames(b.dist.df)[1]<-"Sample"

meta.dist<-b.dist.df%>%inner_join(ps.sampledata,.,by="Sample") # distances with metadata

all.dist<-meta.dist%>%
  select(all_of(.$Sample))%>% # extract distances (corresponding to Sample vector)
  as.dist()

## PCoA ####
## Bray
b.pcoa<-cmdscale(b.dist,
                 k=2,
                 eig = TRUE,
                 add = TRUE) # PCoA
# k is the num of principal coordinates
# eig allows to see % of variation explained
# add rescales eigenvalues to make them all positive
b.positions<-b.pcoa$points # pcoa values to plot
colnames(b.positions)<-c("pcoa1", "pcoa2")

b.percent_explained<-round(100* b.pcoa$eig / 
                             sum(b.pcoa$eig),1)
# previous won't have 0 after decimal

# format() command will show exact number of digits you need 
b.percent_exp<-format(round(100* b.pcoa$eig / 
                              sum(b.pcoa$eig),1),1)


# Envfit on Bray PCoA ggplot2 ####
fit<-envfit(b.pcoa,all.wide,perm=999)
pvals<-fit$vectors$pvals[fit$vectors$pvals<0.05]
fit_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)
fit_coord_cont$veclen<-sqrt(fit_coord_cont$Dim1^2+fit_coord_cont$Dim2^2)
fit_coord_cont<-fit_coord_cont[rownames(fit_coord_cont)%in%names(pvals),]
fit_coord_cont<-fit_coord_cont%>%arrange(-veclen)

envfit.bray.pcoa.plot<-ggplot(data = b.pcoa.positions, aes(x = pcoa1, y = pcoa2,
                                                           color = class)) + 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  guides(fill="none")+
  theme_bw()+
  geom_point() + 
  scale_colour_manual(breaks = custom.levels,
                      labels=custom.color.labels,
                      values = custom.colors) + 
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+ 
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = fit_coord_cont, aes(x = Dim1, y = Dim2), colour = "grey30", 
            fontface = "bold", label = row.names(fit_coord_cont)) +
  ggtitle(paste0("PCoA on Bray-Curtis distances between different rodents
                 (",agglom.rank, " level)"))+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 25), 
        legend.text = ggtext::element_markdown(size = 20),
        legend.position = "right",
        plot.caption = ggtext::element_markdown(hjust = 0, size=20),
        plot.caption.position = "plot")+
  labs(x=paste0("PCo 1 (", b.percent_exp[1],"%)"),
       y=paste0("PCo 2 (", b.percent_exp[2],"%)"),
       color="Host")
ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-bray-pcoa-all-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=envfit.bray.pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")

## Show only top 10 most significant taxa of envfit bray ####
envfit.bray.pcoa.top10<-ggplot(data = b.pcoa.positions, 
                               aes(x = pcoa1, y = pcoa2,color = class)) + 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  # geom_point(size = 3, alpha = 0.5)+
  geom_point()+
  theme_bw() + 
  guides(fill="none")+
  scale_colour_manual(breaks = custom.levels,
                      labels=custom.color.labels,
                      values = custom.colors) + 
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  # top 10 most significant taxa
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont[1:10,], size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = fit_coord_cont[1:10,], aes(x = Dim1, y = Dim2), colour = "grey30", 
            fontface = "bold", label = row.names(fit_coord_cont[1:10,])) +
  ggtitle(paste0("PCoA on Bray-Curtis distances between different rodents 
  (",agglom.rank, " level). Top 10 most significant taxa (p<0.05)"))+
  
  xlab("PCo 1")+
  ylab("PCo 2")+
  theme(axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 25), 
        legend.text = ggtext::element_markdown(size = 20),
        legend.position = "right",
        plot.caption = ggtext::element_markdown(hjust = 0, size=20),
        plot.caption.position = "plot")+
  labs(colour = "Host")

ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-bray-pcoa-all-top10-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=envfit.bray.pcoa.top10,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")

## Boxplots for significant bacteria that drive differences ####
top10drivers<-as.vector(rownames(fit_coord_cont[which(sort(pvals,decreasing = FALSE)<0.05)[1:10],]))

top10.df<-ps.q.df[ps.q.df$Taxon%in%top10drivers,]
# TODO: fix colors
top10.envfit.uni.unw.plot<-ggplot(top10.df,aes(x=factor(top10.df$class,
                                                        custom.levels),
                                               y=Abundance
                                               # ,fill=class
))+
  geom_boxplot()+
  facet_wrap(Taxon~.,scales = "free",ncol=2)+
  theme_bw()+
  # scale_color_manual(breaks = unname(pretty.facet.labels),
  #                    labels=unname(pretty.facet.labels))+
  
  scale_x_discrete(labels=pretty.facet.labels,
                   limits=custom.levels)+
  # scale_fill_manual(values = custom.colors)+
  coord_flip()+
  ggtitle("Relative abundances of top 10 most significant taxa")+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank())

ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-bray-boxplots-all-top10-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=top10.envfit.uni.unw.plot,
       width = 5500,height = 5000,
       units = "px",dpi=300,device = "png")

# % explained by each axis is the value in eig divided 
# by the sum of eig and multiplied by 100
# we create a vector of % (each axis is divided by sum(eig))
# so, we actually have % for all axes
# first two axes are percent_explained[1:2]
b.pcoa.positions<-b.positions %>% 
  as_tibble(rownames="Sample") %>%
  inner_join(.,ps.sampledata, by="Sample")  



# Unweighted UniFrac ####
# Use relative abundances
ps.q.df <-ps.q.agg.rel%>%
  select(Sample,OTU,Abundance,class,Taxon)
all.wide<-ps.q.df%>%
  select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

# colnames are OTUs and rownames are sample IDs
rownames(all.wide)<-all.wide$Sample
all.wide<-all.wide[,-c(1,2)]  
all.wide<-as.matrix(all.wide)

uni.unw<-UniFrac(ps.q, weighted = FALSE)
uni.unw.df<-uni.unw%>% # convert distances into tibble
  as.matrix()%>%
  as_tibble()

uni.unw.df<-cbind(labels(uni.unw),uni.unw.df) # create a column of sample ids
colnames(uni.unw.df)[1]<-"Sample"

## Unweighted unifrac PCoA ####
uni.unw.pcoa<-cmdscale(uni.unw,
                       k=2,
                       eig = TRUE,
                       add = TRUE)

uni.unw.positions<-uni.unw.pcoa$points # pcoa values to plot
colnames(uni.unw.positions)<-c("pcoa1", "pcoa2")

uni.unw.percent_explained<-round(100* uni.unw.pcoa$eig / 
                                   sum(uni.unw.pcoa$eig),1)
# previous won't have 0 after decimal

# format() command will show exact number of digits you need 
uni.unw.percent_exp<-format(round(100* uni.unw.pcoa$eig / 
                                    sum(uni.unw.pcoa$eig),1),1)

# % explained by each axis is the value in eig divided 
# by the sum of eig and multiplied by 100
# we create a vector of % (each axis is divided by sum(eig))
# so, we actually have % for all axes
# first two axes are percent_explained[1:2]
uni.unw.pcoa.positions<-uni.unw.positions %>% 
  as_tibble(rownames="Sample") %>%
  inner_join(.,ps.sampledata, by="Sample")  

## Envfit PCoA ####
fit<-envfit(uni.unw.pcoa,all.wide,perm=999)
pvals<-fit$vectors$pvals[fit$vectors$pvals<0.05]
fit_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)
fit_coord_cont$veclen<-sqrt(fit_coord_cont$Dim1^2+fit_coord_cont$Dim2^2)
fit_coord_cont<-fit_coord_cont[rownames(fit_coord_cont)%in%names(pvals),]
fit_coord_cont<-fit_coord_cont%>%arrange(-veclen)


### Envfit
envfit.uni.unw.pcoa.plot<-ggplot(data = uni.unw.pcoa.positions, 
                                 aes(x = pcoa1, y = pcoa2,color = class)) + 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  # geom_point(size = 3, alpha = 0.5)+
  geom_point()+
  theme_bw() + 
  guides(fill="none")+
  scale_colour_manual(breaks = custom.levels,
                      labels=custom.color.labels,
                      values = custom.colors) + 
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = fit_coord_cont, aes(x = Dim1, y = Dim2), colour = "grey30", 
            fontface = "bold", label = row.names(fit_coord_cont)) +
  ggtitle(paste0("PCoA on unweighted UniFrac distances between different rodents
                 (",agglom.rank, " level)"))+
  
  xlab("PCo 1")+
  ylab("PCo 2")+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 25), 
        legend.text = ggtext::element_markdown(size = 20),
        legend.position = "right",
        plot.caption = ggtext::element_markdown(hjust = 0, size=20),
        plot.caption.position = "plot")+
  labs(colour = "Host")

ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-uni-unw-pcoa-all-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=envfit.uni.unw.pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-uni-unw-pcoa-all-",agglom.rank,
              "-",truncationlvl,
              ".tiff"),plot=envfit.uni.unw.pcoa.plot,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "tiff")

## Show only top 10 most significant taxa of envfit unweighted unifrac ####
envfit.uni.unw.pcoa.top10<-ggplot(data = uni.unw.pcoa.positions, 
                                  aes(x = pcoa1, y = pcoa2,color = class)) + 
  stat_ellipse(geom = "polygon",
               level = 0.8,
               alpha=0.2,
               show.legend = FALSE)+
  # geom_point(size = 3, alpha = 0.5)+
  geom_point()+
  theme_bw() + 
  guides(fill="none")+
  scale_colour_manual(breaks = custom.levels,
                      labels=custom.color.labels,
                      values = custom.colors) + 
  scale_fill_manual(name=NULL, breaks = custom.levels,
                    labels=custom.levels,
                    values = custom.colors)+
  # top 10 most significant taxa
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               data = fit_coord_cont[1:10,], size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = fit_coord_cont[1:10,], aes(x = Dim1, y = Dim2), colour = "grey30", 
            fontface = "bold", label = row.names(fit_coord_cont[1:10,])) +
  ggtitle(paste0("PCoA on unweighted UniFrac distances between different rodents 
  (",agglom.rank, " level). Top 10 most significant taxa (p<0.05)"))+
  
  xlab("PCo 1")+
  ylab("PCo 2")+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 25), 
        legend.text = ggtext::element_markdown(size = 20),
        legend.position = "right",
        plot.caption = ggtext::element_markdown(hjust = 0, size=20),
        plot.caption.position = "plot")+
  labs(colour = "Host")
ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-uni-unw-pcoa-all-top10-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=envfit.uni.unw.pcoa.top10,
       width = 4000,height = 3000,
       units = "px",dpi=300,device = "png")

## Boxplots for significant bacteria that drive differences ####
top10drivers<-as.vector(rownames(fit_coord_cont[which(pvals<0.05)[1:10],]))

top10.df<-ps.q.df[ps.q.df$Taxon%in%top10drivers,]
# TODO: fix colors
top10.envfit.uni.unw.plot<-ggplot(top10.df,aes(x=factor(top10.df$class,
                                                        custom.levels),
                                               y=Abundance
                                               # ,fill=class
))+
  geom_boxplot()+
  facet_wrap(Taxon~.,scales = "free",ncol=2)+
  theme_bw()+
  # scale_color_manual(breaks = unname(pretty.facet.labels),
  #                    labels=unname(pretty.facet.labels))+
  
  scale_x_discrete(labels=pretty.facet.labels,
                   limits=custom.levels)+
  # scale_fill_manual(values = custom.colors)+
  coord_flip()+
  ggtitle("Relative abundances of top 10 most significant taxa (p<0.05)")+
  theme(axis.title = element_text(size = 20), 
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank())

ggsave(paste0("./images/diversity/",Sys.Date(),"-envfit-uni-unw-boxplots-all-top10-",agglom.rank,
              "-",truncationlvl,
              ".png"),plot=top10.envfit.uni.unw.plot,
       width = 5500,height = 5000,
       units = "px",dpi=300,device = "png")
