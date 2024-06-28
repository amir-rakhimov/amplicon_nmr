# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(c("tidyverse","vegan","Polychrome"))
# library(phyloseq)
library(tidyverse)
library(vegan)
library(Polychrome)
authorname<-"pooled"
truncationlvl<-"234"
agglom.rank<-"Genus"

# Import data ####
read.end.type<-"single"
rare.status<-"rare"
filter.status<-"nonfiltered"
date_time<-"20240426_21_44_30"
image.formats<-c("png","tiff")

load(file.path("./output/rdafiles",paste(
  date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

pretty.level.names<-
  c("NMR" = "*Heterocephalus glaber*", # better labels for facets
    "B6mouse" = "B6 mouse",
    "MSMmouse" = "MSM/Ms mouse",
    "FVBNmouse" = "FVB/N mouse",
    "DMR" = "*Fukomys Damarensis*",
    "hare" = "*Lepus europaeus*",
    "rabbit" = "*Oryctolagus cuniculus*",
    "spalax" = "*Nannospalax leucodon*",
    "pvo" = "*Pteromys volans orii*",
    "NMRwt"="Wild *Heterocephalus glaber*"
  )

excluded.samples<-
  c(#"MSMmouse",
    #"FVBNmouse",
    "NMRwt")
custom.levels<-intersect(names(pretty.level.names),custom.md$class)

# facet labels
metric.labs=c('sobs'="Richness \n(Observed species)",
              'shannon' = "Shannon",
              # 'simpson' = "Simpson",
              'invsimpson' = "Inverse \nSimpson")
# metrics to plot
plot.metrics<-c("sobs","shannon", # "simpson",
                "invsimpson")
div.indices<-c("sobs","shannon",# "simpson",
               "invsimpson")

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900",
                                          "#FFA07A"))
names(custom.fill)<-custom.levels
swatch(custom.fill)



# filter your data
if(exists("excluded.samples")){
  custom.levels<-custom.levels[!custom.levels%in%excluded.samples]
  pretty.level.names<-
    pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  pretty.level.names<-pretty.level.names[!names(pretty.level.names)%in%excluded.samples]
  ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,!class%in%excluded.samples,Abundance!=0)
}else{
  pretty.level.names<-
    pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
  ps.q.agg<-ps.q.agg%>%
    filter(class%in%custom.levels,Abundance!=0)
}

# load the output of 003-phyloseq-rarefaction-filtering.R file ####
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      "20240426_22_00_04",
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)


# Exclude samples if you need ####
if(exists("excluded.samples")){
  ps.q.df <-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,!class%in%excluded.samples,Abundance!=0)%>%
    select(Sample,Abundance,class,all_of(agglom.rank),sex)
  custom.fill<-custom.fill[!names(custom.fill)%in%excluded.samples]
}else{
  ps.q.df <-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,Abundance!=0)%>%
    select(Sample,Abundance,class,all_of(agglom.rank),sex)
  
}

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
# The Kruskal-Wallis rank sum test analyzes differences in each alpha diversity
# metric across animal hosts
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
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
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
        legend.position = "none")+
  ggtitle(paste0("Alpha diversity of the gut microbiota of different rodents \n(",agglom.rank, " level)"))

div.plot.with.dots<-div.plot+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4) # add dots
  

for(image.format in image.formats){
  ggsave(paste0("./images/diversity/alpha/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "alpha",
                      paste(plot.metrics,collapse = "-"),
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = image.format)
  ggsave(paste0("./images/diversity/alpha/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "alpha-dots",
                      paste(plot.metrics,collapse = "-"),
                      paste(custom.levels,collapse = '-'),
                      agglom.rank,truncationlvl,
                      sep = "-"),".",image.format),
         plot=div.plot.with.dots,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = image.format)
}

# Add significance stars ####
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


# save plot with significance bars
# Use the table of wilcoxon tests, check if there are any pairwise comparisons
# that were significant
# if yes, save the plot
if(table(w.results<0.05)[2]>0){
  # ggsave(paste0("./images/diversity/alpha/",
  #               paste(Sys.Date(),"alpha",
  #                     paste(plot.metrics,collapse = "-"),
  #                     paste(custom.levels,collapse = '-'),
  #                     agglom.rank,truncationlvl,sep = "-"),
  #               "-signif.png"),
  #        plot=newplot,
  #        width = 6000,height = 5000,
  #        units = "px",dpi=300,device = "png")
  # 
  # ggsave(paste0("./images/diversity/alpha/",
  #               paste(Sys.Date(),"alpha",
  #                     paste(plot.metrics,collapse = "-"),
  #                     paste(custom.levels,collapse = '-'),
  #                     agglom.rank,truncationlvl,sep = "-"),
  #               "-signif.tiff"),
  #        plot=newplot,
  #        width = 6000,height = 5000,
  #        units = "px",dpi=300,device = "tiff")
  for(image.format in image.formats){
    ggsave(paste0("./images/diversity/alpha/",
                  paste(paste(format(Sys.time(),format="%Y%m%d"),
                              format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                        "alpha",
                        paste(plot.metrics,collapse = "-"),
                        paste(custom.levels,collapse = '-'),
                        agglom.rank,truncationlvl,"signif",
                        sep = "-"),".",image.format),
           plot=newplot,
           width = 6000,height = 5000,
           units = "px",dpi=300,device = image.format)
  }
}

