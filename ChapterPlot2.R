###########################
# Plot Figure 5
###########################

# Load data
library(graphicalExtremes)
data <- danube$data_clustered
n <- nrow(data)
dat <- data
d <- ncol(data)

# Transform and truncate
Frechettrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1)) 
river.ext <- apply(dat,2,Frechettrans)
norm_vec <- function(x) sqrt(sum(x^2)) 
norms<-apply(river.ext,1,norm_vec)

q <- 0.9
river.ext<-river.ext[norms>quantile(norms,q),]
norms<-apply(river.ext,1,norm_vec)
river.ext<-river.ext/norms 

###########################
# k-means
###########################

set.seed(800)
library(skmeans)
k <- 6
fit <- skmeans(river.ext,k,method="pclust",control = list(nruns = 500, maxchains=100))
fit$cluster

# Re-order the clusters.  
cluster.id <- rep(NA,length(fit$cluster))
cluster.id[fit$cluster==1] <- 1
cluster.id[fit$cluster==5] <- 2
cluster.id[fit$cluster==3] <- 3
cluster.id[fit$cluster==2] <- 4
cluster.id[fit$cluster==4] <- 5
cluster.id[fit$cluster==6] <- 6

# Plot the observations in heatmap
library('plot.matrix')
test.matrix <- river.ext[order(cluster.id),]
cluster <- sort(cluster.id)
row.names(test.matrix) <- as.character(sort(cluster.id))
colnames(test.matrix) <- as.character(1:31)
pdf('1Kmeans.pdf',12,10)
plot(test.matrix,main='Spherical K-means clustering')
dev.off()

## Elbow plot, if needed ##
# value<-rep(0,30)
# for (k in 2:30){
#   value[k]<-skmeans(river.ext,k,method="pclust",control = list(nruns = 500))$value
#   print(k)
# }
# par(mfrow=c(1,1))
# plot(2:30,value[2:30],ylab='objective function', xlab="k")
# 


###########################
# k-pc
###########################

# Use the KPC algorithm from Fomichov and Ivanov
set.seed(800)
source('PCandSpherical.R')
fitPC <- clusterPC(river.ext,k)

# Extract cluster allocations
n.ext <- nrow(river.ext)
cluster.id.PC <- rep(NA,n.ext)
for (i in 1:n.ext){
  cluster.id.PC[i] <- which.max(fitPC%*%river.ext[i,])
}

# Re-order
cluster.id <- rep(NA,n.ext)
cluster.id[cluster.id.PC==2] <- 1
cluster.id[cluster.id.PC==5] <- 2
cluster.id[cluster.id.PC==3] <- 3
cluster.id[cluster.id.PC==1] <- 4
cluster.id[cluster.id.PC==6] <- 5
cluster.id[cluster.id.PC==4] <- 6

# Plot
# After each cluster, an empty row to distinguish the clusters
library('plot.matrix')
test.matrix <- river.ext[order(cluster.id),]
cluster <- sort(cluster.id)
row.names(test.matrix) <- as.character(sort(cluster.id))
colnames(test.matrix) <- as.character(1:31)
test.matrix1 <- c()
rowname.vec <- c()
for (j in 1:6){
  test.matrix1 <- rbind(test.matrix1,river.ext[cluster.id==j,])
  test.matrix1 <- rbind(test.matrix1,rep(-2,d))
  rowname.vec <- c(rowname.vec,rep(as.character(j),sum(cluster.id==j)),'')
}
row.names(test.matrix1) <- rowname.vec 
colnames(test.matrix1) <- as.character(1:31)
pdf('1Kpc.pdf',12,10)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(test.matrix1,main='Spherical K-means and K-pc clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')
dev.off()

###########################
# spectral clustering
###########################

set.seed(800)

dist_matrix <- as.matrix(dist(river.ext,diag=T,upper=T))
similarity_matrix <- exp(-dist_matrix^2/(2 * 1^2))
k.neighbour <- 10
for (i in 1:n.ext){
  row.i <- similarity_matrix[i,]
  row.i[order(row.i)[1:(n.ext-k.neighbour)]] <- 1
  similarity_matrix[i,] <- row.i
}
W <- matrix(NA,n.ext,n.ext)
for (i in 1:n.ext){
  for (j in 1:n.ext){
    W[i,j] <- sqrt(similarity_matrix[i,j]*similarity_matrix[j,i])
  }
}

D12 <- diag(1/sqrt(rowSums(W)))
L <- diag(1,nrow=n.ext) - D12 %*% W %*% D12
eigenvectors <- eigen(L)$vectors
k <- 6
U <- eigenvectors[, 1:k]
for (i in 1:n.ext){
  U[i,] <- U[i,]/sqrt(sum(U[i,]^2))
}
cluster_assignments <- kmeans(U, centers = k)$cluster


cluster.id <- rep(NA,n.ext)
cluster.id[cluster_assignments ==1] <- 1
cluster.id[cluster_assignments ==5] <- 2
cluster.id[cluster_assignments ==4] <- 3
cluster.id[cluster_assignments ==2] <- 4
cluster.id[cluster_assignments ==3] <- 5
cluster.id[cluster_assignments ==6] <- 6
library('plot.matrix')
test.matrix <- river.ext[order(cluster.id),]
cluster <- sort(cluster.id)
row.names(test.matrix) <- as.character(sort(cluster.id))
colnames(test.matrix) <- as.character(1:31)


test.matrix1 <- c()
rowname.vec <- c()
for (j in 1:6){
  test.matrix1 <- rbind(test.matrix1,river.ext[cluster.id==j,])
  test.matrix1 <- rbind(test.matrix1,rep(-2,d))
  rowname.vec <- c(rowname.vec,rep(as.character(j),sum(cluster.id==j)),'')
}
row.names(test.matrix1) <- rowname.vec 

colnames(test.matrix1) <- as.character(1:31)

pdf('1Spectral.pdf',12,10)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(test.matrix1,main='Spectral clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')
dev.off()
