library(ALDEx2)
data(selex)
head(selex)
#subset only the last 400 features for efficiency
selex.sub <- selex[1:400,]
conds<-c(rep("NS",7),rep("S",7))
x.all <- aldex(selex.sub, # Counts table 
               conds,  # vector of conditions
               mc.samples=16, # num of monte-carlo (DMC) instances
               test="t", 
               effect=TRUE,
               include.sample.summary=FALSE, 
               denom="all", # indicating if iqlr, zero or all feature are used 
                            # as the denominator is required
               verbose=FALSE, 
               paired.test=FALSE)
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")


# aldex2 modules ####
## The aldex.clr module ####
# generates random instances of the centred log-ratio transformed values.

# There are three inputs: counts table, a vector of conditions, and the number 
# of Monte-Carlo (DMC) instances; and several parameters: a string indicating 
# if iqlr, zero or all feature are used as the denominator is required, 
# and level of verbosity (TRUE or FALSE). 

# the output is in the S3 object 'x'
x <- aldex.clr(selex.sub, conds, mc.samples=16, denom="all", verbose=F)
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
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

## The aldex.kw and the aldex.glm modules ####
# Alternatively to the t-test, the user can perform the glm and Kruskal Wallace 
# tests for one-way ANOVA of two or more conditions. 
# Here there are only two inputs: the aldex object from aldex.clr, and 
# the vector of conditions
x.kw <- aldex.kw(x)
# The aldex.glm module is the preferred method for testing of more than two 
# conditons or for complex study designs.

## The aldex.effect module ####
# Finally, we estimate effect size and the within and between condition values
# in the case of two conditions.
x.effect <- aldex.effect(x, CI=T, verbose=FALSE, paired.test=FALSE)
# There is one input: the aldex object from aldex.clr; and several useful 
# parameters. It is possible to include the 95% confidence interval information
# for the effect size estimate with the boolean flag CI=TRUE. This can be helpful
# when deciding whether to include or exclude specific features from consideration

## The aldex.plot module ####
x.all <- data.frame(x.tt,x.effect)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch")
aldex.plot(x.all, type="MW", test="welch")

sig<-x.all$effect>2
sig<-which(x.all$effect>2 & x.all$effect.low>0)

# CI doesn't cross 0
sig<-which(x.all$effect.low>0 & x.all$effect.high>0)
sig<-intersect(which(x.all$effect>2), sig)
# which(x.all$effect.low<0 & x.all$effect.high<0)

aldex.plot(x.all, type="MA", test="welch")
points(x.all$rab.all[sig],
       x.all$diff.btw[sig],col="blue")

aldex.plot(x.all, type="MW", test="welch")
points(x.all$diff.win[sig],
       x.all$diff.btw[sig],col="blue")

plot(x.all$diff.win,
     x.all$diff.btw,pch=1)
sig<-x.all$effect<0.1
points(x.all$diff.win[sig],
       x.all$diff.btw[sig],col="blue")

# Complex case ####
data(selex)
selex.sub <- selex[1:500, ]
set.seed(1)
covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
                         "B" = c(rep(0, 7), rep(1, 7)),
                         "Z" = sample(c(1,2,3), 14, replace=TRUE))
mm <- model.matrix(~ A + Z + B, covariates) # a, z, b, intercept
mm <- model.matrix(~ B + Z + A, covariates) # a, z, b, intercept
mm <- model.matrix(~ A + B + Z, covariates) # a, z, b, intercept

x <- aldex.clr(selex.sub, mm, mc.samples=4, denom="all", verbose=F)
glm.test <- aldex.glm(x, mm)

glm.effect <- aldex.glm.effect(x,CI=T) # list of A and B. But no Z 
# Z is not supported

aldex.plot(glm.effect[["B"]], test="effect", cutoff=2)
sig <- glm.test[,20]<0.05
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")
sig <- glm.test[,20]<0.2
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")

View(glm.effect[["B"]])
sig<-which(glm.effect[["B"]]$effect.low>0 & glm.effect[["B"]]$effect.high>0)
# intersect(which(glm.effect[["B"]]$effect>2),sig)

aldex.plot(glm.effect[["B"]], test="effect", cutoff=2)
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")


# aldex.glm ####
data(selex)
#subset for efficiency
selex.sub <- selex[1:500,]
set.seed(1)
selex.sub[,9:14]<-selex.sub[,9:14]+round(rnorm(nrow(selex.sub),mean = 500))
covariates <- data.frame("A" = c(rep(1, 4), rep(0, 10)),
                         "B" = c(rep(0, 4), rep(1, 4),rep(0, 6)),
                         "Z" = c(rep(0, 8), rep(1, 6)))
mm <- model.matrix(~ A + B+Z-1, covariates) 

# this is also ok
covariates<-c(rep("A", 4), rep("B", 4), rep("Z", 6))
mm <- model.matrix(~ covariates) 

x <- aldex.clr(selex.sub, mm, mc.samples=1, denom="all")
glm.test <- aldex.glm(x,mm)
glm.effect <- aldex.glm.effect(x, CI= T)

aldex.plot(glm.effect[["B"]], test="effect", cutoff=2)
sig <- glm.test[,15]<0.05
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")
sig <- glm.test[,15]<0.2
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")

# aldex.plot(glm.effect[["A"]], test="effect", cutoff=2)
# sig <- glm.test[,15]<0.05
# points(glm.effect[["A"]]$diff.win[sig],
#        glm.effect[["A"]]$diff.btw[sig], col="blue")
# sig <- glm.test[,15]<0.2
# points(glm.effect[["A"]]$diff.win[sig],
#        glm.effect[["A"]]$diff.btw[sig], col="blue")

aldex.plot(glm.effect[["Z"]], test="effect", cutoff=2)
sig <- glm.test[,15]<0.05
points(glm.effect[["Z"]]$diff.win[sig],
       glm.effect[["Z"]]$diff.btw[sig], col="blue")
sig <- glm.test[,15]<0.2
points(glm.effect[["Z"]]$diff.win[sig],
       glm.effect[["Z"]]$diff.btw[sig], col="blue")


# Effect size matters more
sig<-which((glm.effect[["covariatesB"]]$effect.low>0 & glm.effect[["covariatesB"]]$effect.high>0)|
             (glm.effect[["covariatesB"]]$effect.low<0 & glm.effect[["covariatesB"]]$effect.high<0))
sig<-intersect(which(abs(glm.effect[["covariatesB"]]$effect)>2), sig)
# which(x.all$effect.low<0 & x.all$effect.high<0)

aldex.plot(glm.effect[["covariatesB"]], test="effect", cutoff=2)
points(glm.effect[["covariatesB"]]$diff.win[sig],
       glm.effect[["covariatesB"]]$diff.btw[sig],col="blue")


# features with large effect and CI not overlapping 0
rownames(glm.effect[["covariatesB"]][sig,])

###########

sig<-which(glm.effect[["covariatesZ"]]$effect.low>0 & glm.effect[["covariatesZ"]]$effect.high>0)
sig<-intersect(which(glm.effect[["covariatesZ"]]$effect>2), sig)
aldex.plot(glm.effect[["covariatesZ"]], test="effect", cutoff=2)
points(glm.effect[["covariatesZ"]]$diff.win[sig],
       glm.effect[["covariatesZ"]]$diff.btw[sig],col="blue")


# features with large effect and CI not overlapping 0
rownames(glm.effect[["covariatesZ"]][sig,])
