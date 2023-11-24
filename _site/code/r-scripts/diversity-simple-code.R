ps.foo<-prune_taxa(taxa_sums(ps)>0,ps)

ps.foo = rarefy_even_depth(ps.foo,
                           sample.size=17000, 
                           replace=F)

rich=estimate_richness(ps.foo,
                       measures = c("Shannon", "Observed",
                                    "ACE","Simpson","Fisher",
                                    "Chao1"))
kruskal.test(rich$Shannon,sample_data(ps.foo)$class,
             p.adjust.method='BH')

w.sh<-pairwise.wilcox.test(rich$Shannon,
                           sample_data(ps)$class,
                           p.adjust.method = "BH")

w.sh$p.value

kruskal.test(rich$ACE,sample_data(ps.foo)$class,
             p.adjust.method='BH')

w.sh<-pairwise.wilcox.test(rich$ACE,
                           sample_data(ps)$class,
                           p.adjust.method = "BH")

w.sh$p.value


kruskal.test(rich$Chao1,sample_data(ps.foo)$class,
             p.adjust.method='BH')

w.sh<-pairwise.wilcox.test(rich$Chao1,
                           sample_data(ps)$class,
                           p.adjust.method = "BH")

w.sh$p.value

kruskal.test(rich$Simpson,sample_data(ps.foo)$class,
             p.adjust.method='BH')

w.sh<-pairwise.wilcox.test(rich$Simpson,
                           sample_data(ps)$class,
                           p.adjust.method = "BH")

w.sh$p.value

kruskal.test(rich$Fisher,sample_data(ps.foo)$class,
             p.adjust.method='BH')

w.sh<-pairwise.wilcox.test(rich$Fisher,
                           sample_data(ps)$class,
                           p.adjust.method = "BH")

w.sh$p.value


kruskal.test(rich$Observed,sample_data(ps.foo)$class,
             p.adjust.method='BH')
w.sh<-pairwise.wilcox.test(rich$Observed,
                           sample_data(ps)$class,
                           p.adjust.method = "BH")

w.sh$p.value


ps.foo<-ps.foo%>% # agglomerate by agglom.rank
  subset_taxa(Kingdom=="d__Bacteria")
plot_richness(ps.foo, x="animal", 
              measures=c("Observed", "Shannon","Simpson")) + 
  geom_boxplot()