# The function ggplot.species creates a faceted barplot of abundance for each
# taxon per sample (only NMR). It uses a vector with names of taxa, a dataframe with
# abundances, and a taxonomic rank (e.g. Genus). It creates age groups for the 
# barplot and creates levels for the color palette.
# The function keeps only taxa that are found in the taxa.to.plot. vector using
# filter(get(tax.rank)%in%taxa.to.plot). It doesn't matter if we plot only
# one taxon or multiple. Taxonomic rank also doesn't matter.
# The barplot is faceted by taxon and bars are colored by age group. The
# palette is created with Polychrome library. The function returns a ggplot object.

# If we are plotting species from a vector of names, use filter(Species%in%species.to.plot)
# If we are plotting all species from a bigger group like entire genus, use
# filter(get(agglom.rank)==species.to.plot)
ggplot_species<-function(taxa.to.plot,
                         tax.df,
                         tax.rank,
                         sample.order,
                         group.names,
                         grouping.variable,
                         metadata.df,
                         ggplot.fill.name,
                         ggplot.fill.vector){
  ggplot.object<-tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot)
  
  # order the plot by species vector
  taxa.to.plot<-gsub("_"," ",taxa.to.plot)
  taxa.to.plot<-paste0("<i>",taxa.to.plot,"</i>")
  # !!tax.rank:= taxa (bang bang) evaluates the tax.rank before the rest is evaluated
  # to substitute an environment-variable (created with <-) with a data-variable (inside a data frame).
  # https://rlang.r-lib.org/reference/topic-inject.html
  ggplot.object<-ggplot.object%>%
    left_join(metadata.df[c("Sample",grouping.variable)])%>%
    mutate(!!tax.rank:=gsub("_"," ",get(tax.rank)),
           !!tax.rank:=paste0("<i>",get(tax.rank),"</i>"),
           !!tax.rank:=factor(get(tax.rank),levels=taxa.to.plot))%>%
    mutate(Sample=factor(Sample,levels=sample.order$Sample))%>%
    group_by_at(c("class",tax.rank))%>%
    ggplot(aes(x=Sample,
               y=RelativeAbundance,
               fill=factor(get(grouping.variable))))+
    geom_bar(stat="identity")+
    facet_wrap(~get(tax.rank),
               scales = "free",
               ncol = 2)+
    theme_bw()+
    labs(x="",
         y="Relative abundance (%)",
         fill=ggplot.fill.name)+
    scale_color_manual(breaks = unname(group.names),
                       labels=unname(group.names))+
    scale_x_discrete(labels=sample.order$Sample,
                     limits=sample.order$Sample)+
    scale_fill_manual(values = ggplot.fill.vector,
                      labels=group.names)+
    theme(axis.title.y = element_text(size = 25),
          axis.title = element_text(size = 20),
          axis.text.y = ggtext::element_markdown(size=18),
          axis.text.x = element_text(size=20),
          strip.text.x = ggtext::element_markdown(size=20),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  return(ggplot.object)
}
