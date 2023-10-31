library(vegan)
library(tidyverse)
library(phyloseq)
library(Polychrome)

truncationlvl<-"234"
agglom.rank<-"Genus"

# Import data ####
read.end.type<-"single"
rare.status<-"rare"
# rare.status<-"nonrare"
# filter.status<-"filtered"
filter.status<-"nonfiltered"

load(paste0("./rdafiles/pooled-",read.end.type,"-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

pretty.axis.labels<-
  c("NMR" = "*Heterocephalus glaber*", # better labels for facets
    "B6mouse" = "B6 mouse",
    "MSMmouse" = "MSM/Ms mouse",
    "FVBNmouse" = "FVB/N mouse",
    "DMR" = "*Fukomys Damarensis*",
    "hare" = "*Lepus europaeus*",
    "rabbit" = "*Oryctolagus cuniculus*",
    "spalax" = "*Nannospalax leucodon*",
    "pvo" = "*Pteromys volans orii*"
  )

custom.levels<-intersect(names(pretty.axis.labels),custom.md$class)

# load the output of 003-phyloseq-rarefaction-filtering.R file
ps.q.df.preprocessed<-read.table(paste0("./rtables/alldir/ps.q.df.",
                                         rare.status,".",filter.status,"-",agglom.rank,"-",
                                         paste(custom.levels,collapse = '-'),".tsv"),
                                  header = T,sep = "\t")

metric.labs=c('sobs'="Richness \n(Observed species)",
              'shannon' = "Shannon",
              # 'simpson' = "Simpson",
              'invsimpson' = "Inverse \nSimpson")
plot.metrics<-c("sobs","shannon", # "simpson",
                "invsimpson") # metrics to plot
div.indices<-c("sobs","shannon",# "simpson",
               "invsimpson")
scale.color.labels<-unname(pretty.axis.labels)
scale.color.breaks<-unname(pretty.axis.labels)

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

ps.q.df <-ps.q.df.preprocessed%>%
  filter(class%in%custom.levels)%>%
  select(Sample,Abundance,class,Taxon,sex)# select(Sample,OTU,Abundance,class,Taxon)

# ps.q.df <-ps.q.agg.rel%>%
#   select(Sample,OTU,Abundance,class,Taxon)

# ps.q.df<-ps.q.df[!grepl("Unclassified|Uncultured",ps.q.df$taxa.full),] # optional?

# Alpha diversity ####
## Compute alpha diversity metrics ####
all.div<-ps.q.df%>%
  group_by(Sample)%>%
  reframe(sobs=specnumber(Abundance), # richness (num of species)
            shannon=diversity(Abundance,index = "shannon"),
            # simpson=diversity(Abundance, index="simpson"),
            invsimpson=diversity(Abundance, index="invsimpson"), # inverse simpson
            tot=sum(Abundance),
            class=class)%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()

## Alpha diversity tests ####
kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
colnames(kt.results)<-div.indices
rownames(kt.results)<-c("statistic", "pvalue")

combinations<-combn(custom.levels,2) # all unique pairwise combinations
w.results<-data.frame(matrix(nrow = ncol(combinations),ncol=length(div.indices))) # ncol(combinations) pairwise comparisons
colnames(w.results)<-div.indices
w.results<-array(dim = c(length(table(combinations[1,])), # nrows
                         length(table(combinations[2,])), # ncols
                         length(div.indices)), # num of 2D arrays (stacking) 
                 dimnames = list(NULL, NULL, div.indices))



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

# convert multidimensional array into data frame
# w.results.df<-w.results%>%
#   as.data.frame.table()%>%
#   drop_na()
# colnames(w.results.df)<-c("class1","class2","div.metric","p.value")


max.values<-all.div%>%
  group_by(metric)%>%
  summarise(max_val=max(value))


# max.values.pairs<-max.values%>%
#   cross_join(max.values)%>% # create pairwise combinations of max values for each metric and class
#   filter(metric.x==metric.y,class.x!=class.y)%>% # remove pairs of same class 
#   rename(class1=class.x,class2=class.y)

# columns_to_concatenate <- c("class1", "class2")

# create class pairs where classes are sorted
# max.values.pairs$class.pair <- 
#   apply(max.values.pairs[columns_to_concatenate], 1, function(row) {
#     sorted_values <- c(pmin(row[1], row[2]), pmax(row[1], row[2]))
#     paste(sorted_values, collapse = ", ")
#   })
# 
# max.values.pairs<-max.values.pairs%>%
#   group_by(class.pair,metric.x)%>%
#   distinct(class.pair,metric.x,.keep_all = TRUE)%>% # removed
#   rename(metric=metric.x)%>% 
#   ungroup()%>%
#   select(-metric.y)%>% # remove unnecessary columns
#   group_by(metric,class1,class2)%>%
#   mutate(max_val=max(max_val.x,max_val.y))%>% # get max quartile for each pair of boxplots
#   
#   select(-max_val.x,-max_val.y)%>% # remove unnecessary columns
#   relocate(class1,class2)%>%
#   ungroup()
# 
# 
# w.results.df$class.pair <- 
#   apply(w.results.df[columns_to_concatenate], 1, function(row) {
#     sorted_values <- c(pmin(row[1], row[2]), pmax(row[1], row[2]))
#     paste(sorted_values, collapse = ", ")
#   })
# 
# 
# w.results.max_val<-w.results.df%>%
#   inner_join(max.values,by=c("div.metric"="metric"))
# 
# w.results.max_val<-w.results.max_val%>%
#   as_tibble()%>%
#   mutate(class1=factor(class1,levels=custom.levels),
#          class2=factor(class2,levels=custom.levels))%>%
#   mutate(pos1=as.numeric(class1),
#          pos2=as.numeric(class2),
#          pos=min(pos1,pos2),
#          x=if_else(pos1==pos, pos1, pos2),
#          xend=if_else(pos1==pos, pos2, pos1),
#          y=max_val,
#          yend=max_val)
# for(metric in plot.metrics){
#   w.results.max_val[w.results.max_val$div.metric==metric,"y"]<-seq(from=1, by=0.08,
#                length.out=nrow(w.results.max_val[w.results.max_val$div.metric==metric,]))*
#     w.results.max_val[w.results.max_val$div.metric==metric,"y"]
#   
# }
# w.results.max_val$yend<-w.results.max_val$y


all.div$class<-factor(all.div$class,levels=custom.levels)

## Plot alpha diversity metrics ####
div.plot<-ggplot(all.div[all.div$metric %in%
                           plot.metrics,],
                 # aes(x=reorder(class,-value),y=value,fill=class))+
                 aes(x=factor(all.div$class,
                              level=custom.levels),y=value,fill=factor(class)))+
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

coord.combs<-combn(seq_along(custom.levels),2)

# First, find significant results (stars)
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
  metric=factor(rep(plot.metrics,each=ncol(coord.combs)), # each plot metric repeated by the number of combinations
                levels = plot.metrics), # names of our metrics
  label=stars.list
)
star.indices<-which(stars$label=="*") # only significant results

freqs<-as.data.frame(table(stars))%>%filter(label=="*")
yvalues<-c()
for (i in seq_along(freqs)){
  vec<-seq(from=1, by=0.08,length.out=freqs[i,"Freq"])*
    as.numeric(max.values[max.values$metric==freqs[i,"metric"],"max_val"])
  yvalues<-c(yvalues,vec)
}


# the horizontal lines
# merge two vectors: two rows
# c() turns them into a vector

# y coordinates
# multiply vectors to get a matrix: nrows is the number of comparisons
# ncols is the number of plot.metrics
# seq is the sequence of values that will multiply maxvalues

# start and end values can be found from pairwise combinations
# x values have dimensions: num of metrics * num of pairwise comparisons
xvalues<-rep(coord.combs[1,],
             length=length(div.indices)*ncol(coord.combs)) # first row is x start values
xendvalues<-rep(coord.combs[2,],
                length=length(div.indices)*ncol(coord.combs)) # second row is x end values

# select only significant results
xvalues<-xvalues[star.indices]
xendvalues<-xendvalues[star.indices]

horizontal.lines<-tibble(
  metric=factor(rep(plot.metrics,each=ncol(coord.combs)),
                levels = plot.metrics)[star.indices], # names of our metrics
  x = xvalues,       #1 2 1 1 2 1 1 2 1 1 2 1
  xend = xendvalues, #2 3 3 2 3 3 2 3 3 2 3 3
  y =yvalues,
  yend = yvalues
)




# labels depend on our tests (statistical significance)
xstars<-(xvalues+xendvalues)/2
stars<-stars%>% filter(label=="*")%>%
  mutate(x= xstars,
         y =c(yvalues)*1.01)



newplot<-div.plot+
  geom_segment(data=horizontal.lines, # add horizontal.lines of significance
               aes(x=x, xend=xend, y=y, yend=yend),
               inherit.aes = FALSE)+# no conflict with different fills
  geom_text(data = stars,aes(x=x, y=y, label=label),
            inherit.aes = FALSE,size=10) # add stars 

ggsave(paste0("./images/diversity/alpha/",Sys.Date(),"-alpha-",
              paste(plot.metrics,collapse = "-"),"-",agglom.rank,"-",truncationlvl,
              ".png"),plot=div.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")

ggsave(paste0("./images/diversity/alpha/",Sys.Date(),"-alpha-",
              paste(plot.metrics,collapse = "-"),"-",agglom.rank,"-",truncationlvl,
              ".tiff"),plot=div.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "tiff")


ggsave(paste0("./images/diversity/alpha/",Sys.Date(),
              "-alpha-shannon-sobs-",agglom.rank,"-signif-"
              ,truncationlvl,
              ".png"),plot=newplot,
       width = 6000,height = 5000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/diversity/alpha/",Sys.Date(),
              "-alpha-shannon-sobs-",agglom.rank,"-signif-"
              ,truncationlvl,
              ".tiff"),plot=newplot,
       width = 6000,height = 5000,
       units = "px",dpi=300,device = "png")
