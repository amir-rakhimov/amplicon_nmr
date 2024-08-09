library(tidyverse)
library(phyloseq)
library(vegan)
library(ANCOMBC)
# Import data ####
ref.level<-"NMR"
truncationlvl<-"234"
agglom.rank<-"Genus"
read.end.type<-"single"
authorname<-"pooled"
date_time<-"20240426_21_44_30"

load(file.path("./output/rdafiles",paste(
  date_time,
  authorname,read.end.type,"qiime2",
  truncationlvl,agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

custom.levels<-c("NMR",
                 "B6mouse",
                 "MSMmouse",
                 "FVBNmouse",
                 "DMR",
                 "hare",
                 "rabbit",
                 "spalax",
                 "pvo"#,
                 # "NMRwt"
                 )
rare.status<-"rare"
filter.status<-"nonfiltered"

ps.q.df.ancombc.input<-read.table(
  file.path("./output/rtables",authorname,paste0(
    paste(
      "20240426_22_00_04",
      "ps.q.df.rare-nonfiltered",agglom.rank,
      paste(custom.levels,collapse = '-'),sep = "-"),
    ".tsv")),
  header = T)
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)
custom.md<-custom.md%>%
  filter(class%in%custom.levels)

# ANCOMBC ####
# Input data and metadata ####
# convert the data frames into wide format
ps.q.df.ancombc.input.wide<-ps.q.df.ancombc.input%>%
  filter(class%in%custom.levels,Abundance!=0)%>%
  dplyr::select(Sample,Abundance,all_of(agglom.rank))%>%
  pivot_wider(names_from = all_of(agglom.rank), # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()%>%
  rename(Sample.ID=Sample)
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.ancombc.input.wide)<-ps.q.df.ancombc.input.wide$Sample
ps.q.df.ancombc.input.wide<-ps.q.df.ancombc.input.wide[,-1] 

if(agglom.rank=="OTU"){
  taxmat<-ps.q.agg%>%
    dplyr::select(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>%
    distinct()%>%
    column_to_rownames(var = "OTU")%>%
    as.matrix()
}else{
  all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")
  agglom.rank.index<-match(agglom.rank,all.ranks)
  custom.ranks<-all.ranks[1:agglom.rank.index]
  
  taxmat<-ps.q.agg%>%
    ungroup()%>%
    dplyr::select(all_of(custom.ranks))%>%
    distinct()%>%
    column_to_rownames(var = all_of(agglom.rank))%>%
    as.matrix()
}

ps.q.OTU<-t(ps.q.df.ancombc.input.wide)
ps.q.OTU<-otu_table(ps.q.OTU,taxa_are_rows = T)
ps.q.TAX<-tax_table(taxmat)
ps.q.phyloseq.new<-phyloseq(otu_table(ps.q.OTU),
                            tax_table(ps.q.TAX),
                            sample_data(custom.md))
ancombc.levels<-c(ref.level,custom.levels[custom.levels!=ref.level])
sample_data(ps.q.phyloseq.new)$class<-factor(sample_data(ps.q.phyloseq.new)$class,
                                             levels = ancombc.levels)

ancombc.out<-ancombc(
  phyloseq = ps.q.phyloseq.new,
  # tax_level = ,
  formula = "class",
  p_adj_method = "fdr", 
  prv_cut = 0, # by default prevalence filter of 10% is applied
  lib_cut = 0, 
  group = "class",
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
ancombc.res<-ancombc.out$res
ancombc.res_global<-ancombc.out$res_global

## 4.2 ANCOMBC primary result ####
### LFC ####
tab_lfc = ancombc.res$lfc
# col_name = c("Taxon", "Intercept", "Age", "NE - CE", "SE - CE", 
#              "US - CE", "Overweight - Obese", "Lean - Obese")
# colnames(tab_lfc) = col_name
tab_lfc %>% 
  DT::datatable(caption = "Log Fold Changes from the Primary Result")# %>%
  # formatRound(col_name[-1], digits = 2)

### SE  ####
tab_se = ancombc.res$se
# colnames(tab_se) = col_name
tab_se %>% 
  DT::datatable(caption = "SEs from the Primary Result")# %>%
  # formatRound(col_name[-1], digits = 2)

### Test statistic ####
tab_w = ancombc.res$W
# colnames(tab_w) = col_name
tab_w %>% 
  DT::datatable(caption = "Test Statistics from the Primary Result") #%>%
  # formatRound(col_name[-1], digits = 2)

### P-values  ####
tab_p = ancombc.res$p_val
# colnames(tab_p) = col_name
tab_p %>% 
  DT::datatable(caption = "P-values from the Primary Result")# %>%
  # formatRound(col_name[-1], digits = 2)

### Adjusted p-values ####
tab_q = ancombc.res$q
# colnames(tab_q) = col_name
tab_q %>% 
  DT::datatable(caption = "Adjusted p-values from the Primary Result")# %>%
  # formatRound(col_name[-1], digits = 2)


### Differentially abundant taxa ####
tab_diff = ancombc.res$diff_abn
# colnames(tab_diff) = col_name
tab_diff %>% 
  DT::datatable(caption = "Differentially Abundant Taxa from the Primary Result")

# find differentially abundant taxa by multiplying fold change with TRUE/FALSE
# for diff abund
df_lfc = data.frame(ancombc.res$lfc[, -1] * ancombc.res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancombc.res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
ancombc.signif.features<-subset(df_lfc,rowSums(df_lfc[,-c(1,2)]!=0)==ncol(df_lfc[,-c(1,2)]))

# significantly decreased in all hosts
uniq.decreased<-subset(ancombc.signif.features,
                         rowSums(ancombc.signif.features[,-c(1,2)]<0)==ncol(ancombc.signif.features[,-c(1,2)]))
uniq.increased<-subset(ancombc.signif.features,
                           rowSums(ancombc.signif.features[,-c(1,2)]>0)==ncol(ancombc.signif.features[,-c(1,2)]))

ps.q.agg%>%
  filter(get(agglom.rank)=="Paludicola",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Paludicola (Ruminococcaceae)")

ps.q.agg%>%
  filter(get(agglom.rank)=="Alistipes",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Alistipes (Rikenellaceae)")


# Visualization for class
df_fig_class = ancombc.signif.features %>% 
  select(-`(Intercept)`)%>%
  mutate(across(where(is.numeric),\(x) round(x,2)))%>%
  pivot_longer(cols = colnames(ancombc.signif.features)[3:ncol(ancombc.signif.features)], 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_class$value))
up = ceiling(max(df_fig_class$value))
mid = (lo + up)/2
p_class = df_fig_class %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to NMR") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_class

ps.q.agg%>%
  filter(get(agglom.rank)=="Fibrobacter",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Fibrobacter (Fibrobacteraceae)")


df_decreased = uniq.decreased %>% 
  select(-`(Intercept)`)%>%
  mutate(across(where(is.numeric),\(x) round(x,2)))%>%
  pivot_longer(cols = colnames(uniq.decreased)[3:ncol(uniq.decreased)], 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_class$value))
up = ceiling(max(df_fig_class$value))
mid = (lo + up)/2
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

ps.q.agg%>%
  filter(get(agglom.rank)=="Prevotellaceae_UCG-003",
         class%in%custom.levels)%>%
  ggplot(aes(x=factor(class,level=custom.levels),y=Abundance,fill=class))+
  geom_boxplot()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  ggtitle("Prevotellaceae_UCG-003 (Prevotellaceae)")

save.image(file.path("./output/rdafiles/",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "ancombc",rare.status,filter.status,agglom.rank,
  paste(custom.levels,collapse = '-'),
  truncationlvl, ref.level,"workspace.RData",sep="-")))

write.table(ancombc.signif.features,
            file.path("./output/rtables",authorname,paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "ancombc",rare.status,
  filter.status,agglom.rank,
  paste(custom.levels,collapse = '-'),truncationlvl,
  ref.level,"signif.tsv",sep="-")),
  row.names = F,sep = "\t")

