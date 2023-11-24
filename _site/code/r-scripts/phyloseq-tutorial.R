library(phyloseq)
library(ggplot2)
library(vegan)
ps<-qza_to_phyloseq(features = "./data/table-moving-pics.qza",
                     taxonomy = "./data/taxonomy-moving-pics.qza",
                metadata="./data/sample-metadata.tsv",
                tree = "./data/rooted-tree-moving-pics.qza",
                )
sample_data(ps)
rarecurve(t(matrix(otu_table(ps))), step=50, cex=0.5)

# rarefy without replacement
ps.rarefied = 
  rarefy_even_depth(ps, rngseed=1, 
                    sample.size=0.9*min(sample_sums(ps)), 
                    replace=F)

# Bar plots
plot_bar(ps.rarefied, fill="Phylum")
plot_bar(ps.rarefied, fill="Phylum") + 
  facet_wrap(~body.site, scales="free_x", nrow=1)

ps.phylum = tax_glom(ps.rarefied, 
                     taxrank="Phylum", NArm=FALSE)
ps.phylum

plot_bar(ps.phylum, fill="Phylum")
plot_bar(ps.phylum, fill="Phylum") + 
  facet_wrap(~body.site, scales= "free_x", nrow=1)

plot_richness(ps.rarefied, x="days.since.experiment.start", 
              color="body.site", 
              measures=c("Observed"))
plot_richness(ps.rarefied, x="body.site", 
              measures=c("Observed", "Shannon")) + 
  geom_boxplot()

rich = estimate_richness(ps.rarefied)
rich
pairwise.wilcox.test(rich$Observed, 
                     sample_data(ps.rarefied)$body.site)

# Beta diversity
# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps.rarefied, 
                                   method="unifrac", 
                                   weighted=F)
ordination = ordinate(ps.rarefied, method="PCoA", 
                      distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination, color="body.site") + 
  theme(aspect.ratio=1)

# test whether body sites differ significantly from each other
# Use PERMANOVA
foo<-adonis(wunifrac_dist ~ sample_data(ps.rarefied)$body.site)
