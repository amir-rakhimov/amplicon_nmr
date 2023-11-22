library(tidyverse)
library(vegan)
days_wanted<-c(0:9,141:150)
micedf<-read_tsv("mice.shared")
shared<-micedf %>%
  select(Group,starts_with("Otu"))%>%
  mutate(day=str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day)%>%
  pivot_longer(-Group) %>% 
  group_by(Group)%>%
  mutate(total=sum(value)) %>%
  filter(total>1800) %>% # filter by sample size
  group_by(name) %>%
  mutate(total=sum(value)) %>% # count otu frequency
  filter(total!=0) %>% # drop otus with 0 frequency
  ungroup() %>% 
  select(-total)

head(shared)

rand<-shared %>%
  uncount(value)%>%
  mutate(rand_name=sample(name)) %>%
  select(-name) %>%
  count(Group,rand_name) 

richness<-function(x){
  r<-sum(x>0)
  return(r)
}

shannon<-function(x){
  rabund<-x[x>0]/sum(x)
  -sum(rabund*log(rabund))
}

simpson<-function(x){
  n<-sum(x)
  sum(x * (x-1) / (n * (n-1)))
}

shared %>%
  group_by(Group)%>%
  summarize(sobs=richness(value),
            shannon=shannon(value),
            simpson=simpson(value),
            invsimpson=1/simpson,
            n=sum(value))%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson,simpson),
               names_to="metric")%>%
  ggplot(aes(x=n,y=value))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~metric,nrow=4,scales="free_y")+
  theme_bw()


rand %>%
  group_by(Group)%>%
  summarize(sobs=richness(n),
            shannon=shannon(n),
            simpson=simpson(n),
            invsimpson=1/simpson,
            tot=sum(n))%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson,simpson),
               names_to="metric")%>%
  ggplot(aes(x=tot,y=value))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~metric,nrow=4,scales="free_y")+
  theme_bw()


# Using vegan
rand %>%
  group_by(Group)%>%
  summarize(sobs=specnumber(n), # richness
            shannon=diversity(n,index = "shannon"),
            simpson=diversity(n, index="simpson"),
            invsimpson=1/simpson,
            tot=sum(n))%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson,simpson),
               names_to="metric")%>%
  ggplot(aes(x=tot,y=value))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~metric,nrow=4,scales="free_y")+
  theme_bw()
