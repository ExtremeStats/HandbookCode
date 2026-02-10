#read in data
dat <- read.csv("FrenchFama1970.csv", header = T)
d <- 30  #dimension of data
dim(dat)
head(dat)
rtn <- -dat[, -1]  #negate to get negative returns
head(rtn)
pairs(rtn[,1:4])  #exploratory
rtn <- as.matrix(rtn)


#estimate TPDM
normVal <- sqrt(apply(rtn^2, 1, sum)) #find norm of each obs
library(evir)
hill(normVal, end = 2000)  #hill plot of norms
abline(v = 600)
abline(h = 3.35)
#estimate TPDM by taking large values (in terms of norm of entire vector)
q600 <- quantile(normVal, 1-600/13599)  #14.28 threshold for estimating TPDM
keep <- normVal > q600
wMtxAll <- rtn/normVal
wMtx <- wMtxAll[keep,]
n <- dim(wMtxAll)[1]
nExc <- dim(wMtx)[1]
covMtxEst <- matrix(nrow = d, ncol = d)
for(i in seq(1, d)){
  for(j in seq(1, d)){
    covMtxEst[i,j] <- 1/nExc * sum(wMtx[,i] * wMtx[,j])
  }
}


eig1 <- eigen(covMtxEst)
eig1$vectors[,1] <- -eig1$vectors[,1]  
round(eig1$values[1:10], 3)
#[1] 0.743 0.063 0.033 0.026 0.018 0.014 0.011 0.010 0.009 0.008
round(eig1$vectors[,1:6], 3)

#plot to visualize eigenvectors
library(fields)
evenColors <- function (n = 64)
{
  orig <- c("#00008F", "#FFFFFF", "#8F0000")
  if (n == 3)
    return(orig)
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- seq(0, 1, , 3)
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    hold <- splint(x, rgb.tim[, k], xg)
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}

eigColors <- function(mtx, axes.args = NULL, numCols = 9)
{
  lim <- max(abs(mtx))
  zLim <- c(-lim, lim)
  x <- seq(1, dim(mtx)[1])
  y <- seq(1, dim(mtx)[2])
  d <- dim(mtx)[2]
  mtxPlot <- mtx[, d:1]
  z <- mtxPlot
  image.plot(x, y, z, zlim = zLim, col = evenColors(numCols), xlab = "", ylab = "Eigenvector", axes = F)
}

#plot with sectors ordered by name
#pdf("eigenvectorsFinance.pdf", width = 6.5, height = 2.3)
#par(mar = c(3,4,1,1))
eigColors(eig1$vectors[, 1:5], numCols = 9)
box()
axis(2, at = 1:5, labels = 5:1)
axis(1, at = seq(1,d), labels = F)
text(x=seq(1,d), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), labels=colnames(rtn), srt=45, adj=1, xpd=TRUE, cex = .7)
#dev.off()

#plot with sectors clustered by eigenvector 2
order1 <- sort.list(eig1$vectors[,2])
eigColorsOrder <- function(mtx, axes.args = NULL, numCols = 9, orderList = NULL)
{
  lim <- max(abs(mtx))
  zLim <- c(-lim, lim)
  x <- seq(1, dim(mtx)[1])
  y <- seq(1, dim(mtx)[2])
  if(is.null(orderList)){orderList = 1:30}
  d <- dim(mtx)[2]
  print(d)  
  mtxPlot <- mtx[orderList, d:1]
  z <- mtxPlot
  image.plot(x, y, z, zlim = zLim, col = evenColors(numCols), xlab = "", ylab = "Eigenvector", axes = F)
  box()
  axis(2, at = 1:5, labels = 5:1)
  axis(1, at = x, labels = F)
  print(orderList)  
  text(x=x, y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), labels=colnames(rtn)[orderList], srt=45, adj=1, xpd=TRUE, cex = .7)
}
eigColorsOrder(eig1$vectors[, 1:5], numCols = 9, orderList = order1)


#more advanced clustering
lvl1 <- sort.list(eig1$vectors[,2])
break1 <- lvl1[1:9]
break2 <- lvl1[10:26]
break3 <- lvl1[27:30]
lvl2A <- break1[sort.list(eig1$vectors[break1, 3])]
lvl2B <- break2[sort.list(eig1$vectors[break2, 3])]
lvl2C <- break3[sort.list(eig1$vectors[break3, 3])]
length(unique(c(lvl2A, lvl2B, lvl2C)))
lvl2 <- c(lvl2A, lvl2B, lvl2C)
pdf("eigVecFinanceOrder.pdf", width = 6.5, height = 2.3)
par(mar = c(3,4,1,1))
eigColorsOrder(eig1$vectors[, 1:5], numCols = 9, orderList = lvl2)
dev.off()


#time series plots of principal components
pc <- rtn %*% eig1$vectors
plot(pc[,1], type = 'l')
plot(pc[,2], type = 'l')
plot(pc[,3], type = 'l')

#scatterplot of first two principal components
xMax <- max(abs(pc[,1]))
yMax <- max(abs(pc[,2]))
q100 <- quantile(normVal, 1 - 100/13599)
plotBlack <- normVal > q100
colVec= rep("gray", n)
colVec[plotBlack] <- "black"
plot(pc[,1], pc[,2], xlim = c(-xMax, xMax), ylim = c(-yMax, yMax), xlab = expression(z[i1]), ylab = expression(z[i2]), col = colVec)

#adding dates of large pc's and investigating
which(pc[,1] > 75)  #4498
dat[4498,]
#mines anc coal have the smallest negative returns
pc[4498, c(1,2)]
text(91, -16, "1987/10/19", cex = .6)

which(pc[,1] < -50 & pc[,2] < -7)  #12670
dat[9790, 1]  #2008/10/13
pc[9790, c(1,2)]
text(-62, -7.2, "2008/10/13", cex = .6)
sort(rtn[9790,])
dat[9790,]

#compare risk to 5 industry portfolio categories
ff1 <- ff2 <- ff3 <- ff4 <- ff5 <- rep(0, 30)
ff1[c(1,2,3,6,7,26,27)] <- 1
ff1 <- ff1/sum(ff1)
ff1
ff2[c(9,12,13,15,16,20,24)] <- 1
ff2 <- ff2/sum(ff2)
ff2
ff3[c(14,21,22,23)] <- 1
ff3 <- ff3/sum(ff3)
ff3
ff4[8] <- 1
ff5[c(4,5,10,11,17,18,19,25,28,29,30)] <- 1
ff5 <- ff5/sum(ff5)
ff5
1*(cbind(ff1, ff2, ff3, ff4, ff5) > 0)

large100 <- normVal > q100
obs100 <- rtn[large100,]
ang100 <- wMtxAll[large100,]
coefsPC <- ang100 %*% eig1$vectors[,1:5]
projPC <- t(eig1$vectors[, 1:5] %*% t(coefsPC))

FMtx <- cbind(ff1, ff2, ff3, ff4, ff5)
FtF <- t(FMtx) %*% FMtx
hatMtx <- solve(FtF, t(FMtx))  
coefsFF <- hatMtx %*% t(ang100)  
projFF <-  t(FMtx %*% coefsFF) 

plot3 <- function(num){
  rng <- range(c(ang100[num,], projPC[num,], projFF[num,]))
  plot(ang100[num,], ylim = rng)
  points(1:30, projPC[num,], col = 2)
  points(1:30, projFF[num,], col = 3)
}

sqdDiffPC <- apply((ang100-projPC)^2, 1, sum)
mean(sqdDiffPC) #[1] 0.06980505
sqdDiffFF <- apply((ang100-projFF)^2, 1, sum)
mean(sqdDiffFF) #[1] 0.1508683 





