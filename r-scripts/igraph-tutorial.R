library(tidyverse)
library(igraph)

taxa_table<-ps.q.agg.rel%>%
  select(OTU,Kingdom,Phylum,Class,Order,Family,Genus,Taxon)
metadata<-custom.md
sp_ratio<-ps.q.agg.rel%>%
  select(OTU,Abundance,Sample)

# The taxa_table dataset will contain taxonomic infos about our Otu, 
# metadata refers to informations about our samples, and sp_ratio will 
# contain relative abundance of each Otu per samples

sp_ratio<-sp_ratio%>%
  pivot_wider(names_from = Sample,
              values_from = Abundance)

sp_ratio <- sp_ratio %>% 
  as_tibble() %>% 
  rename_all(tolower) %>% 
  column_to_rownames(var = "otu")

# Next we use prevalence filtering to remove Otus that are not expressed 
# in at least 13 samples

# then we apply Spearman’s correlation on the data set to get 
# a correlation matrix. 

# Afterward, we define a threshold of correlation value to 
# be represented in the network, otherwise we might end up with a network 
# containing hundreds of nodes and edges (which will make it difficult 
# to read and interpret). 

# Next we create the network file and we remove nodes without edges. 

# Finally we have a first look at our plot which is not really 
# informative at this stage.

min.prevalence=13
incidence=sp_ratio
incidence[incidence>0]=1
sp_ratio_filtered <- 
  sp_ratio[which(rowSums(incidence)>=min.prevalence),] ### end of prevalence filtering

sp_correl <- sp_ratio_filtered %>% 
  t() %>% 
  cor(method = "spearman")  ### correlation calculations

sp_correl[abs(sp_correl)<0.65]=0  ### define threshold for correlations

net_work <- graph_from_adjacency_matrix(sp_correl,
                                        mode="lower",
                                        weighted=TRUE, 
                                        diag=FALSE)  ## create network file
net_work <- delete.vertices(net_work,
                            degree(net_work)==0) #remove nodes without edges

plot(net_work, 
     vertex.label = NA, 
     edge.width = 5, 
     vertex.size = 10) ## first plot 


## Adding attributes to the network ####
taxa_tbl <- taxa_table %>% 
  as_tibble() %>% 
  mutate(otu = OTU) #tidy taxonomic infos

net_work_used <- V(net_work)$name %>% 
  as_tibble() %>% 
  mutate(otu = value) %>% 
  select(otu) #extract Taxa represented in the network

v_attr <- sp_ratio %>% ### we create a table of attributes for nodes (vertex)
  rownames_to_column( var = "otu") %>% 
  as_tibble() %>% 
  pivot_longer(-otu, names_to = "sample_id", values_to = "ratio" ) %>% 
  group_by(otu) %>% 
  summarise(rel_abundance = sum(ratio)) %>% 
  inner_join(net_work_used, by = "otu") %>% 
  inner_join(taxa_tbl, by = "otu") %>%  
  mutate(rel_abundance = abs(exp (rel_abundance))) # we join taxonomic infos, 
                      # relative abundance and otu that were only
                    # represented in the network to create the attribute table

network_table <- 
  igraph::as_data_frame(net_work, 'both') ##we convert the network in data frames

network_table$vertices <- network_table$vertices %>%
  as_tibble() %>% 
  inner_join(v_attr, by = c("name"="otu"))  # we add our attribute in the data frames 

net_work1 <- 
  graph_from_data_frame(network_table$edges,
                        directed = F,
                        vertices = network_table$vertices) # we convert back data frames to a network

mom_data <- V(net_work1)$class %>% 
  as_tibble_col(column_name = "species") #formating the class variable as factor (is needed for coloring edges)
mom_data$species <- as_factor(mom_data$species)

color_easy <-  c("pink", "brown", "orange", "dodgerblue", "lightgray")[mom_data$species] #creating a color palette to represent each class levels

V(net_work1)$color <-  color_easy ## we have now the color attributes based on class

E(net_work1)$sign <- E(net_work1)$weight ##we create another attribute from the weight

E(net_work1)$weight <- abs(E(net_work1)$weight) ## we then use absolute value of weight because of the specific layout we will use

### details about the network plot (with class attribute used for coloring vertex)

plot(net_work1, vertex.size= 8, #size of nodes
     edge.width=abs(E(net_work1)$weight)*4, #width of edges (edge is a segment between nodes)
     vertex.label = NA, #remove microbes names
     edge.color=ifelse(E(net_work1)$sign > 0, "blue","red"), # color edges based on weight sign
     layout=layout_with_kk(net_work1)) #specific layout

title_legend1 <- V(net_work1)$class %>% 
  as_factor() %>% 
  levels() #extracting class levels in an object for the legend

legend(x = 0.78, y = 1, title_legend1,
       pch = 21, pt.bg = c("pink", "brown", "orange", "dodgerblue", "lightgray"), 
       pt.cex = 2.5, bty = "n", ncol = 1) #adding legend