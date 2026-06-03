project_theme <- theme(axis.text.y = element_text(size=10), # size of y axis ticks
                       axis.title = element_text(size = 10), # size of axis names
                       legend.text = element_text(size = 10), # size of legend text
                       legend.title = element_text(size = 15), # size of legend title
                       panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank())

alpha.plot.theme <- theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = ggtext::element_markdown(size=18),
                          axis.text.y = element_text(size=18),
                          axis.title = element_text(size = 20),
                          strip.text.x = element_text(size=20),
                          plot.title = element_text(size = 27),
                          legend.text = element_text(size = 20),
                          legend.title = element_text(size = 25))

pca.plot.theme <- theme(axis.text.x = element_text(angle=0,size=20),
                        axis.text.y = element_text(size=20),
                        axis.title = element_text(size = 20),
                        plot.title = element_text(size = 27),
                        legend.text = element_text(size = 15),
                        legend.title = element_text(size = 25),
                        legend.position = "right",
                        plot.caption = ggtext::element_markdown(hjust = 0, size=20),
                        plot.caption.position = "plot")

asv.barplot.theme <- theme(
  legend.title = element_text(size=10),
  legend.text = element_text(size=9),
  legend.key.size = unit(0.3, 'cm'), #change legend key size
  legend.key.spacing.y = unit(0.1, "lines"), # distant between key text
  legend.box.spacing = unit(0.1,"lines"), # space between the plot and the legend box
  axis.title.y = element_text(size = 10),
  axis.text.x = element_text(angle=45,size=10,hjust=1),# rotate 
  strip.text.x = ggtext::element_markdown(size=10),
  panel.spacing = unit(0.8, "cm"), # increase distance between facets
  plot.title = element_text(size = 14), # size of plot title
  plot.caption = element_text(size=8), # size of plot caption
  legend.position = "bottom")

taxa.plot.theme<- project_theme +
  theme(axis.text.y = ggtext::element_markdown(size=10),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=12),
      axis.text.x = element_text(size=10),
      strip.text.x = element_text(size=10),
      plot.title = element_text(size = 8), # size of plot title
      plot.caption = element_text(size=8), # size of plot caption
      legend.position = "none")

diffabund.plot.theme <- theme(
  axis.title.y = element_blank(),
  axis.title = element_text(size = 5),
  axis.text.y = ggtext::element_markdown(size=5),
  # axis.text.y = element_text(size=5),
  axis.text.x = element_text(size=5),
  strip.text.x = ggtext::element_markdown(size=5),
  plot.title = element_text(size =5),
  legend.text = element_text(size = 5),
  legend.title = element_text(size = 5),
  legend.position = "none",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank()
)
# 009
#' Theme for alpha diversity plots:
manuscript.alpha.plot.theme<-theme(#plot.margin=unit(c(1,1,1,2), 'cm'),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = ggtext::element_markdown(size=10),
  axis.text.y = element_text(size=10),
  axis.title = element_text(size = 13),
  strip.text.x = element_text(size=10),
  plot.title = element_text(size = 14),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 10))

#' Theme for PCA scatter plots:
manuscript.pca.plot.theme<-theme(
  axis.text.x = element_text(angle=0,size=10),
  axis.text.y = element_text(size=10),
  axis.title = element_text(size = 10),
  plot.title = element_text(size = 14),
  legend.text = ggtext::element_markdown(size = 10),
  legend.title = element_text(size = 10),
  legend.position = "right",
  plot.caption = ggtext::element_markdown(hjust = 0, size=13),
  plot.caption.position = "plot")

manuscript.asv.barplot.theme <- theme(
  axis.text.y = element_text(size=8), # size of y axis ticks
  axis.text.x = element_text(angle=45,size=8,hjust=1),# rotate 
  axis.title.x = element_text(size = 10), # size of axis names
  axis.title.y = element_text(size = 9, hjust= 0.2),
  strip.text.x = ggtext::element_markdown(size=10),
  legend.text = element_text(size = 7.2), # size of legend text
  legend.title = element_text(size=10),
  legend.key.size = unit(0.3, 'cm'), #change legend key size
  legend.key.spacing.y = unit(0.01, "cm"), # distant between key text
  legend.key.spacing.x = unit(0.05, "cm"), # distant between key text
  legend.box.spacing = unit(0.01,"cm"), # space between the plot and the legend box
  panel.spacing = unit(0.8, "cm"), # increase distance between facets
  plot.title = element_text(size = 8), # size of plot title
  plot.caption = element_text(size=8), # size of plot caption
  legend.position = "bottom")
