library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
library(ALDEx2)
library(Polychrome)
truncationlvl<-"234"
agglom.rank<-"Genus"
rare.status<-"rare"
filter.status<-"nonfiltered"
read.end.type<-"single"
ref.level<-"NMR"
source("r-scripts/make_features_maaslin.R")
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
# custom.levels are the levels that exist in the current metadata
custom.levels<-intersect(names(pretty.axis.labels),custom.md$class)
# custom color palette based on custom levels from metadata
scale.color.labels<-unname(pretty.axis.labels)
scale.color.breaks<-unname(pretty.axis.labels)

# boxplot fill colors
set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900"))
# colors should be assigned to levels
names(custom.fill)<-custom.levels
swatch(custom.fill)

common.decreased<- 
  read.table(file=file.path("./rtables/alldir",
                           paste("significant-features",rare.status,
                                 filter.status,agglom.rank,
                                 paste(custom.levels,collapse = '-'),truncationlvl,
                                 ref.level,"signif.tsv",sep="-")),
            header = F,sep = "\t")

# filter by mean relative abundance: >=1%
common.decreased.filtered<-ps.q.agg%>%
  group_by(Taxon)%>%
  filter(Taxon %in% common.decreased)%>% # taxa from ps.q.agg that are unique for ref
  dplyr::select(Taxon, MeanRelativeAbundance)%>%
  arrange(-MeanRelativeAbundance)%>%
  distinct(Taxon,.keep_all = T)%>%
  filter(MeanRelativeAbundance>=1)%>%
  pull(Taxon)

common.decreased.filtered<-as.character(common.decreased.filtered$V1)
# boxplots of abundances per feature
ggplots<-list()
for (feature_index in seq_along(common.decreased.filtered)){
    # take each taxon one by one from the list of significant features
    # make a boxplot of abundances
    feature.plot<-ps.q.agg%>%
      filter(Taxon==common.decreased.filtered[[feature_index]],
             class%in%custom.levels)%>%
      ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
      geom_boxplot(show.legend = FALSE)+
      # scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      ggtitle(common.decreased.filtered[[feature_index]])+
      theme_bw()+
      xlab("")+
      scale_color_manual(breaks = scale.color.breaks, # custom color palette
                         labels=scale.color.labels)+ # based on our levels
      scale_x_discrete(labels=pretty.axis.labels,
                       limits=custom.levels)+ # rename boxplot labels (x axis)
      scale_fill_manual(values = custom.fill)+ # boxplot fill colors
      theme(plot.margin=unit(c(1,1,1,2), 'cm'),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = ggtext::element_markdown(angle=45,hjust=1,size=10),
            axis.text.y = element_text(size=10),
            axis.title = element_text(size = 10),
            strip.text.x = element_text(size=10),
            plot.title = element_text(size = 15),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 15),
            legend.position = "right")
    ggplots[[common.decreased.filtered[[feature_index]]]]<-feature.plot
    
    ggsave(paste0("./images/taxaboxplots/",
                  paste(Sys.Date(),"taxaboxplot",ref.level,
                        common.decreased.filtered[[feature_index]],
                        rare.status,filter.status,
                        agglom.rank,paste(custom.levels,collapse = '-'),
                        truncationlvl,sep="-"),".png"),
           plot=feature.plot,
           width = 2000,height = 1500,
           units = "px",dpi=300,device = "png")
    
    ggsave(paste0("./images/taxaboxplots/",
                  paste(Sys.Date(),"taxaboxplot",ref.level,
                        common.decreased.filtered[[feature_index]],
                        rare.status,filter.status,
                        agglom.rank,paste(custom.levels,collapse = '-'),
                        truncationlvl,sep="-"),".tiff"),
           plot=feature.plot,
           width = 2000,height = 1500,
           units = "px",dpi=300,device = "tiff")
}
ggplots


ps.q.agg%>%
  filter(Taxon=="Eisenbergiella (Lachnospiraceae)",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Eisenbergiella (Lachnospiraceae)")

# boxplots of abundances per feature
ggplots<-list()
if(agglom.rank=="OTU"){
  for (feature_index in seq_along(common.decreased.filtered)){
    # take each otu one by one from aldex.neg.effect
    # make a boxplot of abundances
    feature.plot<-ps.q.agg%>%
      filter(OTU==common.decreased.filtered[[feature_index]],
             class%in%names(host.labels))%>%
      ggplot(aes(x=factor(class,level=host.labels),y=Abundance,fill=class))+
      geom_boxplot()+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      ggtitle(paste(common.decreased.filtered[[feature_index]],
                    as.character(ps.q.agg[ps.q.agg$OTU==common.decreased.filtered[[feature_index]],
                                          "Taxon"][1,1])))
    ggplots[[common.decreased.filtered[[feature_index]]]]<-feature.plot
  }
}else{ 
  for (feature_index in seq_along(common.decreased.filtered)){
    # take each taxon one by one from aldex.neg.effect
    # make a boxplot of abundances
    feature.plot<-ps.q.agg%>%
      filter(Taxon==common.decreased.filtered[[feature_index]],
             class%in%custom.levels)%>%
      ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
      geom_boxplot()+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      ggtitle(common.decreased.filtered[[feature_index]])
    ggplots[[common.decreased.filtered[[feature_index]]]]<-feature.plot
  }
}
ggplots


# Visualization for class
df_fig_class = ancombc.signif.features %>% 
  dplyr::select(-X.Intercept.)%>%
  mutate(across(where(is.numeric),\(x) round(x,2)))%>%
  pivot_longer(cols = colnames(ancombc.signif.features)[3:ncol(ancombc.signif.features)], 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_class$value))
up = ceiling(max(df_fig_class$value))
mid = (lo + up)/2
df_decreased = ancombc.signif.features %>% 
  dplyr::select(-X.Intercept.)%>%
  filter(taxon_id%in%common.decreased)%>%
  mutate(across(where(is.numeric),\(x) round(x,2)))%>%
  pivot_longer(cols = colnames(ancombc.signif.features)[3:ncol(ancombc.signif.features)], 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
p_decreased_class = df_decreased %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to NMR") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_decreased_class


