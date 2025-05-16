set.seed(4131)
x1 <- c(rnorm(100,68,5),rnorm(100,78,5),rnorm(100,63,4))        ## height
x2 <- c(rnorm(100,215,40),rnorm(100,200,20),rnorm(100,125,15))  ## weight
x3 <- c(rnorm(100,280,50),rnorm(100,180,15),rnorm(100,95,15))   ## exercise
label <- rep(c("Compact Guys","Tall Guys","Women"),each=100)    ## labels


# scale the data
x1 <- scale(x1)
x2 <- scale(x2)
x3 <- scale(x3)

# form the data matrix X and covariance matrix E
x <- cbind(x1,x2,x3)
colnames <- c("Average Guys","Tall Guys","Women")
E <- cov(x)

cor(x)
cov(x)

# eigendecomposition of E by solving det(E - lambda I)=0
evd<-eigen(E)
evd

# Rescale the data X using eigenvectors
ev1 <- evd$vectors[,1]    ## eigenvector 1, corresponding to the largest eigenvalue thus most variation
ev2 <- evd$vectors[,2]    ## eigenvector 2
ev3 <- evd$vectors[,3]    ## eigenvector 3

pc1 <- x %*% ev1
pc2 <- x %*% ev2
pc3 <- x %*% ev3

# Regress the first two PC
df <- data.frame(cbind(pc1,pc2),label)
names(df) <- c('PC1','PC2','Group')
ggplot(df,aes(x=PC1,y=PC2,colour=Group)) + geom_point(alpha=1)



# PCA command
out <- princomp(x)

df <- data.frame(cbind(out$scores[,1],out$scores[,2]),label)
names(df) <- c('PC1','PC2','Group')
ggplot(df,aes(x=PC1,y=PC2,colour=Group)) + geom_point(alpha=1)
cbind(ev1,ev2,ev3)

cbind(ev1,ev2,ev3)
out$loadings

plot(out,type="l")

####
data("USArrests")
head(USArrests)
# calculare PCs
results <-prcomp(USArrests,scale=TRUE)

# reverse the signs
results$rotation<--1*results$rotation

# display PCs
results$rotation
# the first principal component (PC1) has high values for Murder, Assault, 
# and Rape which indicates that this principal component describes the most
# variation in these variables.