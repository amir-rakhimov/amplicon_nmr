library(tidyverse)
library(vegan)
source("read_matrix.R")

dist_matrix <- read_matrix("mice.braycurtis.dist")

dist_tbl <- as_tibble(dist_matrix, rownames="samples")

sample_lookup <- dist_tbl %>% 
  select(samples) %>%
  mutate(delimited = str_replace(samples,
                                 "^(([FM])\\d+)D(\\d+)$",
                                 "\\2-\\1-\\3")) %>%
  separate(col=delimited,
           into=c("sex", "animal", "day"), sep="-",
           convert=TRUE)

days_wanted <- c(0:9, 141:150)

dist_matrix <- dist_tbl %>%
  pivot_longer(cols=-samples, names_to="b", values_to="distances") %>%
  inner_join(., sample_lookup, by="samples") %>%
  inner_join(., sample_lookup, by=c("b" = "samples")) %>%
  filter(day.x %in% days_wanted & day.y %in% days_wanted) %>%
  select(samples, b, distances) %>%
  pivot_wider(names_from="b", values_from="distances") %>%
  select(-samples) %>%
  as.dist()

# PcoA: metric dimensional scaling ####
pcoa<-cmdscale(dist_matrix,
               k=2,
               eig = TRUE,
               add = TRUE) # PCoA
# k is the num of principal coordinates
# eig allows to see % of variation explained
# add rescales eigenvalues to make them all positive
positions<-pcoa$points # pcoa values to plot
colnames(positions)<-c("pcoa1", "pcoa2")

percent_explained<-round(100* pcoa$eig / sum(pcoa$eig),1)
# previous won't have 0 after decimal

# format() command will show exact number of digits you need 
percent_exp<-format(round(100* pcoa$eig / 
                                  sum(pcoa$eig),1),1)

# % explained by each axis is the value in eig divided 
# by the sum of eig and multiplied by 100
# we create a vector of % (each axis is divided by sum(eig))
# so, we actually have % for all axes
# first two axes are percent_explained[1:2]
positions %>% 
  as.tibble(rownames="samples") %>%
  inner_join(.,sample_lookup, by="samples") %>%
  mutate(period=if_else(day<10, "early","late")) %>%
  ggplot(aes(x=pcoa1,y=pcoa2,color=period))+
  geom_point()+
  labs(x=paste0("PCo 1 (", percent_exp[1],"%)"),
       y=paste0("PCo 2 (", percent_exp[2],"%)"))+
  theme_bw()


# Scree plot shows how much variation is explained by axes
tibble(pe=percent_explained,
       axis=1:length(percent_explained)) %>%
  ggplot(aes(x=axis,y=pe)) +
  geom_line()+
  coord_cartesian(xlim = c(1,10))

# Cumulative % explained
tibble(pe=cumsum(percent_explained),
       axis=1:length(percent_explained)) %>%
  ggplot(aes(x=axis,y=pe)) +
  geom_line()+
  coord_cartesian(xlim = c(1,10),ylim=c(0,50))

# NMDS: non-metric multidimensional scaling
# use vegan's metaMDS
set.seed(1)
nmds<-metaMDS(dist_matrix)

scores(nmds) %>%
  as.tibble(rownames="samples") %>%
  inner_join(.,sample_lookup, by="samples") %>%
  mutate(period=if_else(day<10, "early","late")) %>%
  ggplot(aes(x=NMDS1,y=NMDS2,color=period))+
  geom_point()+
  theme_bw()
