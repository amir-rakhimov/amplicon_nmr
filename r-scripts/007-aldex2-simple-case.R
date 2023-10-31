library(tidyverse)
library(phyloseq)
library(Polychrome)
library(pals)
library(ALDEx2)
library(vegan)

# "274-203"
truncationlvl<-"234"
agglom.rank<-"Genus"
# load("./rdafiles/yashliu-dada2-284-203-Genus-138-1-phyloseq-workspace.RData")

load(paste0("./rdafiles/pooled-qiime2-",truncationlvl,"-",agglom.rank,
            "-phyloseq-workspace.RData"))

pretty.facet.labels<-
  c("NMR" = "*Heterocephalus glaber*", # better labels for facets
    "B6mouse" = "SPF mouse, B6",
    "DMR" = "*Fukomys Damarensis*",
    "hare" = "*Lepus europaeus*",
    "rabbit" = "*Oryctolagus cuniculus*",
    "spalax" = "*Nannospalax leucodon*",
    "pal" = "*Petaurista alborufus lena*",
    "pvo" = "*Pteromys volans orii*",
    "ppg" = "*Petaurista philippensis grandis*",
    "tx" = "*Trogopterus xanthipes*"
    
    # ,
    # "rabbitcontrol"="*Oryctolagus cuniculus*",
    # "harecontrol" = "*Lepus europaeus*",
    # "ntccontrol" = "Non-treatment<br>control"
  )
custom.levels<-names(pretty.facet.labels)

# Using Polychrome package
set.seed(1)
# custom.colors<- createPalette(length(custom.levels),
#               seedcolors = rainbow(6))
custom.colors<- 
  createPalette(length(custom.levels),
                seedcolors = c("#B22222", "#0000FF","#006400", 
                               "#FF8C00", "#68228B"))
custom.colors<-unname(custom.colors)
swatch(custom.colors)
# custom colors for scale_fill_manual (maybe not needed)
# custom.fill<-c("dodgerblue", "seagreen","indianred")

# custom labels for scale_color_manual
custom.color.labels<-unname(pretty.facet.labels)

# Import data ####
ps.q.df <-ps.q.agg.abs%>%
  dplyr::select(Sample,OTU,Abundance,class,Taxon)%>%
  filter(Abundance!=0)

ps.q.df<-ps.q.df%>%
  filter(class %in% c("NMR","B6mouse"))


# convert the data frame into wide format
ps.q.df.wide<-ps.q.df%>%
  dplyr::select(-OTU)%>%
  pivot_wider(names_from = "Taxon", # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  

# # Rarefaction ####
# # find the smallest sample size
# min.n_seqs.all<-ps.q.agg.abs%>%
#   # filter(class %in% c("NMR","B6mouse","spalax","DMR","rabbit"))%>%
#   filter(class %in% c("NMR","B6mouse"))%>%
#   select(Sample,OTU,Abundance)%>%
#   group_by(Sample)%>%
#   summarize(n_seqs=sum(Abundance))%>%
#   summarize(min=min(n_seqs))%>%
#   pull(min)
# 
# # rarefied asv table with vegan
# set.seed(1)
# ps.q.df.rare<-rrarefy(ps.q.df.wide,sample=min.n_seqs.all)%>%
#   as_tibble(rownames="Sample")%>%
#   pivot_longer(-Sample)%>%
#   as.data.frame()%>%
#   inner_join(unique(ps.q.agg.abs[,c("Sample","class")]),
#              by="Sample")%>%
#   rename(Taxon=name,Abundance=value)
# condsNMR condsrabbit condsspalax condsB6mouse
ps.q.df.wide<-t(ps.q.df.wide)
covariates<-custom.md$class[match(colnames(ps.q.df.wide),rownames(custom.md))]

# aldex2 all
ps.q.aldex.all <- aldex(ps.q.df.wide, # Counts table 
               covariates,  # vector of conditions
               mc.samples=128, # num of monte-carlo (DMC) instances
               test="t", 
               effect=TRUE,
               include.sample.summary=FALSE, 
               denom="all", # indicating if iqlr, zero or all feature are used 
               # as the denominator is required
               verbose=FALSE, 
               paired.test=FALSE)
par(mfrow=c(1,2))
aldex.plot(ps.q.aldex.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(ps.q.aldex.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")


# aldex2 modules ####
## The aldex.clr module ####
# generates random instances of the centred log-ratio transformed values.

# There are three inputs: counts table, a vector of conditions, and the number 
# of Monte-Carlo (DMC) instances; and several parameters: a string indicating 
# if iqlr, zero or all feature are used as the denominator is required, 
# and level of verbosity (TRUE or FALSE). 

# the output is in the S3 object 'x'
ps.q.aldex.clr <- aldex.clr(ps.q.df.wide, covariates, mc.samples=128, denom="all", verbose=F)
# We recommend 128 or more mc.samples for the t-test, 1000 for a rigorous effect 
# size calculation, and at least 16 for ANOVA

# in fact we recommend that the number of samples in the smallest group 
# multiplied by the number of DMC be equal at least 1000 in order to generate a
# reasonably stable estimate of the posterior distribution


## The aldex.ttest module ####
# performs the Welch’s t and Wilcoxon rank test for the situation when there
# are only two conditions. 
# There are only two inputs: the aldex object from  aldex.clr and whether a
# paired test should be conducted or not (TRUE or FALSE).
ps.q.aldex.tt <- aldex.ttest(ps.q.aldex.clr, paired.test=FALSE, verbose=FALSE)

# we.ep - Expected P value of Welch’s t test

# we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test

# wi.ep - Expected P value of Wilcoxon rank test

# wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test

## The aldex.kw and the aldex.glm modules ####
# Alternatively to the t-test, the user can perform the glm and Kruskal Wallace 
# tests for one-way ANOVA of two or more conditions. 
# Here there are only two inputs: the aldex object from aldex.clr, and 
# the vector of conditions
ps.q.aldex.kw <- aldex.kw(ps.q.aldex.clr)
# The aldex.glm module is the preferred method for testing of more than two 
# conditons or for complex study designs.

# kw.ep - Expected P value of Kruskal-Wallace test

# kw.eBH - Expected Benjamini-Hochberg corrected P value of Kruskal-Wallace test

# glm.ep - Expected P value of glm test

# glm.eBH - Expected Benjamini-Hochberg corrected P value of glm test

## The aldex.effect module ####
# Finally, we estimate effect size and the within and between condition values
# in the case of two conditions.
ps.q.aldex.effect <- aldex.effect(ps.q.aldex.clr, CI=T, verbose=FALSE, paired.test=FALSE)
# There is one input: the aldex object from aldex.clr; and several useful 
# parameters. It is possible to include the 95% confidence interval information
# for the effect size estimate with the boolean flag CI=TRUE. This can be helpful
# when deciding whether to include or exclude specific features from consideration


# rab.all - median clr value for all samples in the feature
# 
# rab.win.NS - median clr value for the NS group of samples
# 
# rab.win.S - median clr value for the S group of samples
# 
# rab.X1_BNS.q50 - median expression value of features in sample X1_BNS if include.item.summary=TRUE
# 
# dif.btw - median difference in clr values between S and NS groups
# 
# dif.win - median of the largest difference in clr values within S and NS groups
# 
# effect - median effect size: diff.btw / max(diff.win) for all instances
# 
# overlap - proportion of effect size that overlaps 0 (i.e. no effect)


## The aldex.plot module ####
ps.q.aldex.all <- data.frame(ps.q.aldex.tt,ps.q.aldex.effect)

par(mfrow=c(1,2))
aldex.plot(ps.q.aldex.all, type="MA", test="welch")
aldex.plot(ps.q.aldex.all, type="MW", test="welch")

# CI doesn't cross 0
sig<-which((ps.q.aldex.all$effect.low>0 & ps.q.aldex.all$effect.high>0)|
             (ps.q.aldex.all$effect.low<0 & ps.q.aldex.all$effect.high<0))
sig<-intersect(which(abs(ps.q.aldex.all$effect)>2), sig) 
# which(ps.q.aldex.all$effect.low<0 & ps.q.aldex.all$effect.high<0)

aldex.plot(ps.q.aldex.all, type="MA", test="welch")
points(ps.q.aldex.all$rab.all[sig],
       ps.q.aldex.all$diff.btw[sig],col="blue")

aldex.plot(ps.q.aldex.all, type="MW", test="welch")
points(ps.q.aldex.all$diff.win[sig],
       ps.q.aldex.all$diff.btw[sig],col="blue")


View(ps.q.aldex.effect[sig,])

intersect(rownames(ps.q.aldex.effect[sig,][ps.q.aldex.effect[sig,]$effect<0,]),feature.list$feature)
intersect(rownames(ps.q.aldex.effect[sig,][ps.q.aldex.effect[sig,]$effect>0,]),feature.list$feature)
