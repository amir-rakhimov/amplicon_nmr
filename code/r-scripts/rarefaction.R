library(vegan)
library(tidyverse)
days_wanted<-c(0:9,141:150)
shared<-read_tsv("mice.shared")
#label, Group, numOtus, Otu1, Otu2 etc

shared<-shared%>%
  select(Group,starts_with("Otu"))%>% # Only Group(i.e. Sample) and OTUs 
  mutate(day=str_replace(Group, ".*D",""))%>%
  filter(day%in%days_wanted)%>%
  select(-day)%>%
  pivot_longer(-Group) %>% #Group, name (otus), value
  group_by(Group)%>%
  mutate(total=sum(value))%>%
  filter(total>1800)%>%
  group_by(name)%>%
  mutate(total=sum(value))%>%
  filter(total!=0)%>%
  ungroup()%>%
  select(-total)%>% # for null model, don't run below
  pivot_wider(Group)%>%
  as.data.frame()


# rownames(shared) <- shared$Group
# shared <- shared[, -1]
# shared <- as.matrix(shared)


# null model: repeat each row by the value in the 'value' column
# e.g otu1 will be shown 115 times
# then, randomise the order of otus and recount everything
rand<-shared%>%
  uncount(value)%>% # to repeat each row and get 1.8 mln rows
  mutate(rand_name=sample(name)) %>% # randomise without replacement
  select(-name)%>%
  dplyr::count(Group, rand_name)

# sum of reads for each sample
rand_group_count<-rand%>%
  group_by(Group)%>%
  summarise(n=sum(n))
shared_group_count<-shared%>%
  group_by(Group)%>%
  summarise(n=sum(value))

inner_join(shared_group_count,rand_group_count,by="Group")

# sum of reads for each otu
shared_otu_count<-shared%>%
  group_by(name)%>%
  summarise(n=sum(value))
rand_otu_count<-rand%>%
  group_by(rand_name)%>%
  summarise(n=sum(n))
inner_join(shared_otu_count,rand_otu_count,by=c("name"="rand_name"))


# each sample is a statistical sample of the overall community distro

# we preserve the total number of sequences in each otu and 
# total num of sequences in each sample

rand_group_count%>%
  ggplot(aes(x=n))+
  geom_histogram()

# turn rand into matrix
rand_df<-rand%>%
  pivot_wider(names_from = "rand_name",
              values_from = "n",
              values_fill = 0)%>%
  as.data.frame()
rownames(rand_df)<-rand_df$Group
rand_df<-rand_df[,-1]

rand_matrix<-as.matrix(rand_df)

# no rarefaction distance matrix
norare_dist_matrix<-vegdist(rand_matrix,method="bray")
# rarefaction distance matrix
rare_dist_matrix<-avgdist(rand_matrix,
                          dmethod="bray",
                          sample=1826) # size is the smallest library size

# compare distance calculation approaches
norare_dist_tibble<-norare_dist_matrix%>%
  as.matrix()%>%
  as_tibble(rownames="sample")%>%
  pivot_longer(-sample)%>%
  filter(name<sample) # only unique comparisons (lower triangle)


rare_dist_tibble<-rare_dist_matrix%>%
  as.matrix()%>%
  as_tibble(rownames="sample")%>%
  pivot_longer(-sample)%>%
  filter(name<sample) # only unique comparisons (lower triangle)

# range of rarefied data is small, while for non-rarefied it's huge: 0-1
# rarefied distances are averaging about the same
comparison<-inner_join(norare_dist_tibble,rare_dist_tibble,by=c("sample","name"))%>%
  select(sample,name,norare=value.x,rare=value.y)%>%
  inner_join(.,rand_group_count,by=c("sample"="Group"))%>% # add total counts of left sample
  inner_join(.,rand_group_count,by=c("name"="Group"))%>% # add total counts of right sample (name)
  mutate(n_diff=abs(n.x-n.y))%>%
  select(-n.x,-n.y)

comparison%>%
  ggplot(aes(x=norare,y=rare,color=n_diff))+
  geom_point(size=0.25,alpha=0.25)+
  geom_smooth()

# make a faceted plot where x axis is the number of sequences that's 
# different between two compared samples, and y axis is distances
comparison%>%
  # split the df: for each comparison, there's a row for dist in rare
  # and one row for dist in norare
  pivot_longer(cols=c("norare","rare"),names_to="type",values_to = "dist")%>%
  ggplot(aes(x=n_diff,y=dist))+
  geom_point(size=0.25,alpha=0.25)+
  facet_wrap(~type,nrow = 2)

# rarefaction removes the variability due to sequencing depth


dist<-vegdist(shared,method = "bray")
set.seed(19760620)
nmds<-metaMDS(dist)
nmds<-metaMDS(shared,autotransform = FALSE)

scores(nmds)%>%
  as_tibble(rownames="Group")%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point()

set.seed(19760620)
dist<-avgdist(shared,dmethod = "bray",sample=1800)
set.seed(1)
nmds<-metaMDS(dist)
scores(nmds)%>%
  as_tibble(rownames="Group")%>%
  ggplot(aes(x=NMDS1,y=NMDS2))+
  geom_point()
