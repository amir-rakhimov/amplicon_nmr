library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(broom)

metadata <- read_excel("schubert.metadata.xlsx", na="NA") %>%
  select(sample_id, disease_stat) %>%
  drop_na(disease_stat)

otu_counts <- read_tsv("schubert.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")

nseqs_per_sample <- otu_counts %>%
  group_by(sample_id) %>%
  summarize(N = sum(count), .groups="drop") %>%
  count(N) %>%
  pull(N)

stopifnot(length(nseqs_per_sample) == 1)

lod <- 100* 1/nseqs_per_sample

taxonomy <- read_tsv("schubert.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{genus}<br>({pretty_otu})")) %>%
  select(otu, taxon)

otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = 100*count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  mutate(disease_stat = factor(disease_stat,
                               levels=c("Case",
                                        "DiarrhealControl",
                                        "NonDiarrhealControl")))


taxon_pool <- otu_rel_abund %>%
  group_by(disease_stat, taxon) %>%
  summarize(median=median(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(median) < 0.5,
            median = max(median),
            .groups="drop")

otu_disease_rel_abund <- inner_join(otu_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", as.character(taxon))) %>%
  group_by(sample_id, disease_stat, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            median = max(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE))


otu_plot<-otu_disease_rel_abund %>%
  mutate(rel_abund = if_else(rel_abund == 0,
                           2/3 * lod,
                           rel_abund)) %>%
  ggplot(aes(y=taxon, x=rel_abund, color=disease_stat)) +
  geom_vline(xintercept = lod, size=0.2) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6)) +
  coord_trans(x="log10") +
  scale_x_continuous(limits=c(NA, 100),
                     breaks=c(0.1, 1, 10, 100),
                     labels=c(0.1, 1, 10, 100)) +
  scale_color_manual(name=NULL,
                     breaks=c("NonDiarrhealControl",
                              "DiarrhealControl",
                              "Case"),
                     labels=c("Healthy",
                              "Diarrhea,<br>*C. difficile* negative",
                              "Diarrhea,<br>*C. difficile* positive"),
                     values=c("gray", "blue", "red")) +
  labs(y=NULL,
       x="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        # legend.position = c(0.8, 0.6),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )




experiment_significance <- otu_disease_rel_abund %>%
  group_by(taxon) %>%
  nest() %>%
  mutate(experiment_tests = map(.x=data,
                                ~kruskal.test(rel_abund~disease_stat, data=.x) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  mutate(p.experiment = p.adjust(p.value, method="BH")) %>%
  select(taxon, data, p.experiment) %>%
  filter(p.experiment < 0.05)

get_max_quartile<-function(x){
  x%>%group_by(disease_stat)%>%
    summarise(third_q=quantile(rel_abund,prob=0.75), .groups="drop")%>%
    summarise(max_quartile=max(third_q))%>%
    pull(max_quartile)
}

pairwise_test<-experiment_significance %>%
  mutate(max_quartile=map_dbl(.x=data,~get_max_quartile(.x)))%>%
  mutate(pairwise_tests = map(.x=data,
                              ~pairwise.wilcox.test(x=.x$rel_abund,
                                                    g=.x$disease_stat,
                                                    p.adjust.method = "BH") %>%
                                tidy())) %>%
  unnest(pairwise_tests) %>%
  filter(p.value < 0.05) %>%
  select(taxon, group1, group2, p.value,max_quartile) %>%
  mutate(pos=as.numeric(taxon),
         y=if_else(group1=="NonDiarrhealControl", pos+0.2, pos),
         yend=if_else(group2=="DiarrhealControl", pos, pos-0.2),
         x=case_when(group1=="NonDiarrhealControl"&
                       group2=="DiarrhealControl"~max_quartile*1.3, 
                     group1=="NonDiarrhealControl"&
                       group2=="Case"~max_quartile*1.6, 
                     group1=="DiarrhealControl"&
                       group2=="Case"~max_quartile*1.9 ),
         xend=x,
         x_star=1.1*x,
         y_star=case_when(group1=="NonDiarrhealControl"&
                            group2=="DiarrhealControl"~pos+0.05, 
                          group1=="NonDiarrhealControl"&
                            group2=="Case"~pos-0.05, 
                          group1=="DiarrhealControl"&
                            group2=="Case"~pos-0.15 ))
  
otu_plot+
  geom_segment(data=pairwise_test,
               aes(x=x,xend=xend,y=y,yend=yend),inherit.aes = FALSE)+
  geom_text(data=pairwise_test,
            aes(x=x_star,y=y_star),label="*",inherit.aes = FALSE)

ggsave("schubert_otu.tiff", width=7, height=6)
