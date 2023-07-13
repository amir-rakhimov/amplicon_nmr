library(ALDEx2)
data(selex) # samples are columns, while features are rows

#subset only the last 400 features for efficiency
selex.sub <- selex[1:400,]

conds <- c(rep("NS", 7), rep("S", 7))
x.all <- aldex(selex.sub,  # counts table
               conds,  # a vector of conditions
               mc.samples=16,  # the number of Monte-Carlo (DMC) instances
               test="t", # evaluates the data as a two-factor experiment using 
               # both the Welch's t and the Wilcoxon rank test
               effect=TRUE, 
               include.sample.summary=FALSE, 
               denom="all", 
               verbose=FALSE)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
    ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
    ylab="Difference")


# The aldex.clr module ####
# the output is in the S3 object 'x'
x <- aldex.clr(selex.sub, conds, mc.samples=16, 
               denom="all", # a string indicating if iqlr, zero or all 
               #feature are used as the denominator is required,
               verbose=F)

# The aldex.ttest module ####
# Welch and wilcoxon rank test
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
head(x.tt)

# The aldex.kw module ####
# Two or more conditions
# glm and kruskal wallis
x.kw <- aldex.kw(x)
head(x.kw)

# The aldex.effect module ####
# we estimate effect size and the within 
# and between condition values in the case of two conditions
x.effect <- aldex.effect(x, CI=T, verbose=FALSE, paired.test=FALSE)
head(x.effect)

# The aldex.plot module ####
x.all <- data.frame(x.tt,x.effect)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch")
aldex.plot(x.all, type="MW", test="welch")
