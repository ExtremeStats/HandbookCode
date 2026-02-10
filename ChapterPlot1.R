#############################################
# Plot Figure 1, 2 and 3
#############################################

library(graphicalExtremes)
library(rgl)

#############################################
# Examples of linear factor models
#############################################

# Left plot from linear factor model without noise

set.seed(31)
n<-1000
Z1<--1/log(runif(n))
Z2<--1/log(runif(n))
X1<-rep(0,n)
X2<-rep(0,n)
for (i in 1:n){
  X1[i]<-0.02*Z1[i]+0.08*Z2[i]
  X2[i]<-0.08*Z1[i]+0.02*Z2[i]
}
plot(X1,X2,xlab=expression(Y[1]),ylab=expression(Y[2]))

# Right plot from linear factor model with noise

set.seed(43)
n<-1000
Z1<--1/log(runif(n))
Z2<--1/log(runif(n))
eta<--1/log(runif(n))
N1<-abs(rnorm(n))
N2<-abs(rnorm(n))
X1noise<-rep(0,n)
X2noise<-rep(0,n)
for (i in 1:n){
  X1noise[i]<-0.02*Z1[i]+0.08*Z2[i]+0.06*N1[i]*eta[i]
  X2noise[i]<-0.08*Z1[i]+0.02*Z2[i]+0.06*N2[i]*eta[i]
}
plot(X1noise,X2noise,xlab=expression(Y[1]),ylab=expression(Y[2]))



#############################################
# 3D plot of river extremes
#############################################

# Load danube river data and select Stations 27,28,29
data <- danube$data_clustered
n <- nrow(data)
i1 <- 27
i2 <- 28
i3 <- 29
dat <- data[,c(i1,i2,i3)]

# Transform the data to Frechet margin
Frechettrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1)) 
river.ext <- apply(dat,2,Frechettrans)

# Select the extreme observations above a quantile
q <- 0.95
norm_vec <- function(x) sqrt(sum(x^2)) 
norms <- apply(river.ext,1,norm_vec)
river.ext <- river.ext[norms>quantile(norms,q),]
norms <- apply(river.ext,1,norm_vec)
river.ext <- river.ext/norms 

# Plot 3D plot of the extreme observation
open3d(windowRect = 20 + c(0,0,600,600))
# par3d(userMatrix = view)
plot3d(x = river.ext[,2], y = river.ext[,3], z = river.ext[,1],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,cex=1.9)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), 
                                  radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font
# Add axis labels. 'line' specifies how far to set the label from the axis.
mtext3d("Station 28",       edge = "x-+", line = 5)
mtext3d("Station 29",       edge = "y-+", line = 2,pos=c(0,0,1))
mtext3d("Station 27",          edge = "z+-", line = 5)

## Note! Can adjust the angle on the 3D plot manually, then use the following command to save the optimal view.##
# view <- par3d("userMatrix")

## To replot using the selected view, add the following line in the beginning of the plotting command and run everything again. ##
## par3d(userMatrix = view) ##

## Save plots
# rgl.snapshot("RiverAng.png", fmt="png")
# rgl.postscript('RiverAng.pdf', fmt = 'pdf')

##########################################
# K-means clustering on river extremes
##########################################

set.seed(800)
library(skmeans)
k <- 3
fit <- skmeans(river.ext,k,method="pclust",control = list(nruns = 1000, maxchains=100))
fit$cluster

open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(x = river.ext[fit$cluster==1,2], 
       y = river.ext[fit$cluster==1,3], 
       z = river.ext[fit$cluster==1,1],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='red',cex=1.9,pch=0)
plot3d(x = river.ext[fit$cluster==2,2], 
       y = river.ext[fit$cluster==2,3], 
       z = river.ext[fit$cluster==2,1],cex=1.9,pch=1,add=T)
plot3d(x = river.ext[fit$cluster==3,2], 
       y = river.ext[fit$cluster==3,3], 
       z = river.ext[fit$cluster==3,1],col='blue',cex=1.9,pch=2,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
legend3d("topright", c("cluster 1", "cluster 2",'cluster 3'), 
         pch=16,col=c('black','red','blue'))
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font
# Add axis labels. 'line' specifies how far to set the label from the axis.
mtext3d("Station 28",       edge = "x-+", line = 5)
mtext3d("Station 29",       edge = "y-+", line = 2,pos=c(0,0,1))
mtext3d("Station 27",          edge = "z+-", line = 5)

# rgl.snapshot("RiverCluster.png", fmt="png")
# rgl.postscript('RiverCluster.pdf', fmt = 'pdf')


#####################
## Elbow plot ##
#####################
value<-rep(0,10)
for (k in 1:10){
  value[k]<-skmeans(river.ext,k,method="pclust",control = list(nruns = 1000))$value
  print(k)
}
par(mfrow=c(1,1))
plot(value[1:10],ylab='objective function', xlab="k")


