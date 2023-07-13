library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)

# trulvl<-"313-229"
# trulvl<-"284-229"
# trulvl<-"284-203"
# trulvl<-"273-203"
# trulvl<-"273-198"
# trulvl<-"265-203"
# trulvl<-"265-198"

truncationlvl<-"284-203"
agglom.rank<-"Genus"
load("./rdafiles/yashliu-qiime2-284-203-Genus-138-1-phyloseq-workspace.RData")


pretty.axis.labels<-c("NMR" = "*Heterocephalus glaber*", # better labels for axes
                       "SPFmouse" = "SPF mouse, B6",
                       "hare" = "*Lepus europaeus*",
                       "rabbit" = "*Oryctolagus cuniculus*"
                       # ,
                       # "rabbitcontrol"="*Oryctolagus cuniculus*",
                       # "harecontrol" = "*Lepus europaeus*",
                       # "ntccontrol" = "Non-treatment<br>control"
                       )
custom.levels<-names(pretty.axis.labels)
metric.labs=c('sobs'="Richness \n(Observed species)",
              'shannon' = "Shannon",
              'simpson' = "Simpson",
              'invsimpson' = "Inverse \nSimpson")
plot.metrics<-c("sobs","shannon","simpson","invsimpson") # metrics to plot
div.indices<-c("sobs","shannon", "simpson", "invsimpson") 
scale.color.labels<-unname(pretty.axis.labels)
scale.color.breaks<-unname(pretty.axis.labels)

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                              seedcolors = c("#FF3030", "#1E90FF","#A2CD5A", 
                                             "#FF8C00", "#BF3EFF"))
custom.fill<-unname(custom.fill)

ps.q.df <-ps.q.agg.abs%>%
  select(Sample,OTU,Abundance,class,Taxon)

# ps.q.df <-ps.q.agg.rel%>%
#   select(Sample,OTU,Abundance,class,Taxon)

# ps.q.df<-ps.q.df[!grepl("Unclassified|Uncultured",ps.q.df$taxa.full),] # optional?


# Alpha diversity ####
## Compute alpha diversity metrics ####
all.div<-ps.q.df%>%
  group_by(Sample)%>%
  summarize(sobs=specnumber(Abundance), # richness (num of species)
            shannon=diversity(Abundance,index = "shannon"),
            simpson=diversity(Abundance, index="simpson"),
            invsimpson=1/simpson, # inverse simpson
            tot=sum(Abundance),
            class=class)%>%
  pivot_longer(cols=c(sobs,shannon,simpson,invsimpson),
               names_to="metric")%>%
  distinct()

## Alpha diversity tests ####
kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
colnames(kt.results)<-div.indices
rownames(kt.results)<-c("statistic", "pvalue")

combinations<-combn(custom.levels,2) # all unique pairwise combinations
w.results<-data.frame(matrix(nrow = ncol(combinations),ncol=4)) # ncol(combinations) pairwise comparisons
colnames(w.results)<-div.indices
w.results<-array(dim = c(length(table(combinations[1,])), # nrows
                         length(table(combinations[2,])), # ncols
                         length(div.indices)), # num of 2D arrays (stacking) 
                 dimnames = list(NULL, NULL, div.indices))



# TODO: Fix wilcoxon!
for (div.metric in div.indices) {
  metric.ds<-all.div%>%
    filter(metric==div.metric)%>%
    distinct()
  # perform kruskal-wallis test
  kt<-kruskal.test(value~class,data=metric.ds)
  
  kt.results["statistic",div.metric]<- kt$statistic
  kt.results["pvalue",div.metric]<- kt$p.value
  # low pvalue indicates statistical difference between some groups
  # But which groups are different from others?
  # Perform pairwise wilcoxon test
  
  # NB: You must not run pairwise wilcoxon test unless kruskal wallis results 
  # are significant
  if(kt$p.value<0.05){
    w.test<-pairwise.wilcox.test(metric.ds$value,
                                 metric.ds$class,
                                 p.adjust.method = "BH"
                                 ,exact=FALSE
                                 )
    w.results[,,div.metric]<-w.test$p.value
    
  }else(
    w.results[,,div.metric]<-matrix(data = "n.s.",
                                    nrow = nrow(w.test$p.value),
                                    ncol = ncol(w.test$p.value))
    
  )
  dimnames(w.results)[[1]]<-dimnames(w.test$p.value)[[1]] # change rownames of w.results
  dimnames(w.results)[[2]]<-dimnames(w.test$p.value)[[2]] # change colnames of w.results
}

kt.results
w.results
stopifnot(all(kt.results[2,]<0.05))

## Plot alpha diversity metrics ####
div.plot<-ggplot(all.div[all.div$metric %in%
                           plot.metrics,],
                 # aes(x=reorder(class,-value),y=value,fill=class))+
                 aes(x=factor(all.div$class,
                              level=custom.levels),y=value,fill=class))+
  geom_boxplot(show.legend = FALSE)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             ncol=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+
  scale_color_manual(breaks = scale.color.breaks,
                     labels=scale.color.labels)+
  scale_x_discrete(labels=pretty.axis.labels,
                   limits=custom.levels)+ # rename boxplot labels (x axis)
  scale_fill_manual(values = custom.fill)+
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
  ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))

maxvalues<-all.div%>%
  group_by(metric)%>%
  summarise(max=max(value))%>%
  arrange(.dots=plot.metrics)%>%
  pull(max)*1.03
# merge two vectors: two rows
# c() turns them into a vector
# y coordinates
# multiply vectors to get a matrix: nrows is the number of comparisons
# ncols is the number of plot.metrics
# seq is the sequence of values that will multiply maxvalues
yvalues<-seq(from=1, by=0.05,length.out=ncol(combinations))%*%t(maxvalues)
yvalues<-c(rbind(yvalues))

# start and end values can be found from pairwise combinations
coord.combs<-combn(seq_along(custom.levels),2)
# x values have dimensions: num of metrics * num of pairwise comparisons
xvalues<-rep(coord.combs[1,],
             length=length(div.indices)*ncol(coord.combs)) # first row is x start values
xendvalues<-rep(coord.combs[2,],
                length=length(div.indices)*ncol(coord.combs)) # second row is x end values


lines<-tibble(
  metric=factor(rep(plot.metrics,each=ncol(coord.combs)),
                levels = plot.metrics), # names of our metrics
  x = xvalues,       #1 2 1 1 2 1 1 2 1 1 2 1
  xend = xendvalues, #2 3 3 2 3 3 2 3 3 2 3 3
  y =yvalues,
  yend = yvalues
)

# labels depend on our tests (statistical significance)
xstars<-(xvalues+xendvalues)/2
stars.list<-matrix(NA,nrow=length(div.indices),ncol = ncol(coord.combs))
for (k in 1:dim(w.results)[3]) { # loop over each sub array
  for (i in 1:dim(w.results)[1]) { # then each row
    for (j in 1:dim(w.results)[2]) { # then each column
      if (!is.na(w.results[i,j,k])) { # extract non-NA values
        ith.row<-rownames(w.results[,,k])[i] 
        jth.col<-colnames(w.results[,,k])[j]
        w.val<- w.results[i,j,k]
        
        # e.g. ith.row="NMR"
        # jth.col="hare"
        # find in custom.levels the positions that correspond to c("NMR","hare")
        # the result: level.position = c(1,4)
        levels.position<-which(custom.levels%in%c(ith.row,jth.col))
        # find these positions in the coord.combs to get the column 
        # that stores levels.position
        # result: level.col=3 
        level.col<-which(apply(coord.combs,2,function(x) 
          return(all(x==levels.position))))
        # assign the w.val to stars.list:kth row is diversity metric
        # level.col column is the combination of levels
        # we must assign either "*" or "n.s." depeding on w.result[i,j,k]
        stars.list[k,level.col]<-ifelse(w.val!="n.s.",ifelse(w.val<0.05,"*","n.s."),w.val)
      }
    }
  }
}

stars.list<-as.vector(t(stars.list))

stars<-tibble(
  metric=factor(rep(plot.metrics,each=ncol(coord.combs)),
                levels = plot.metrics), # names of our metrics
  x = xstars, #
  y =c(yvalues)*1.01,
  
  label=stars.list

)

newplot<-div.plot+
  geom_segment(data=lines, # add lines of significance
               aes(x=x, xend=xend, y=y, yend=yend),
               inherit.aes = FALSE)+# no conflict with different fills
  geom_text(data = stars,aes(x=x, y=y, label=label),
            inherit.aes = FALSE,size=3) # add stars 

# ggsave(paste0("./images/diversity/",Sys.Date(),
#               "-alpha-shannon-sobs-",agglom.rank,"-",truncationlvl,
#               ".png"),plot=div.plot,
#        width = 6000,height = 3000,
#        units = "px",dpi=300,device = "png")
# 
# ggsave(paste0("./images/diversity/",Sys.Date(),
#               "-alpha-shannon-sobs-",agglom.rank,"-",truncationlvl,
#               ".tiff"),plot=div.plot,
#        width = 4000,height = 3000,
#        units = "px",dpi=300,device = "tiff")


ggsave(paste0("./images/diversity/",Sys.Date(),
              "-alpha-shannon-sobs-",agglom.rank,"-total-"
              ,truncationlvl,
              ".png"),plot=newplot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")
