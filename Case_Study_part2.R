rm(list=ls()) # Clear R environment
library(evd)

# Specify path to the Ch14_TSExtremes/ folder
path.to.folder <- ".../Ch14_TSExtremes/"

source(paste0(path.to.folder,"functions/MarkovFitsGaussZ.R"))
source(paste0(path.to.folder,"functions/ForwardSimulationFunctions.R"))
source(paste0(path.to.folder,"functions/GPDfunctions.R"))
source(paste0(path.to.folder,"functions/ClusterStatisticsFunctions.R"))

# Find data fites in Ch14_TSExtremes folder
path.to.data <- paste0(path.to.folder,"Data/")
file.nms <- list.files(path.to.data)
data.file.nms <- file.nms[substr(file.nms,1,3)=="TN_"]
data.file.nms

# Parse lines of .txt data file and store in matrix
data.strt <- 0
txt <- readLines(paste0(path.to.data,data.file.nms))
header <- c("SOUID","DATE","TN","Q_TN")
for(l in 1:length(txt)){
  print(l)
  line.txt <- strsplit(gsub(" ","",txt),",")[[l]]
  if(length(line.txt)==length(header)){
    if(sum(line.txt==header)==length(header)){
      data.strt <- l
      break
    }
  }
}
txt.data <- txt[(data.strt+1):length(txt)]
mat.txt.data <- do.call(rbind, strsplit(gsub(" ","",txt.data),","))

# Select only the DATE and TN columns
ind.cols <- which(header %in% c("DATE","TN"))
df.data <- data.frame(date = mat.txt.data[,ind.cols[1]],
                      X    = as.numeric(mat.txt.data[,ind.cols[2]]))

# Remove missing data
df.data$X[df.data$X == -9999] <- NA
df.data$X[df.data$X == -999] <- NA

# Make temperatures Celsius
df.data$X <- df.data$X/10

# Format the dates in 4 different columns
df.data$year <- as.numeric(substr(df.data$date,1,4))
df.data$month <- as.numeric(substr(df.data$date,5,6))
df.data$day <- as.numeric(substr(df.data$date,7,8))
df.data$dt <- paste0(substr(df.data$date,1,4),"-",
                     substr(df.data$date,5,6),"-",
                     substr(df.data$date,7,8))

############################
### Summer sub-selection ###
############################

# Select months of June, July, August
mnths <- c(6,7,8)
df.data <- df.data[df.data$month %in% mnths,]
df.data <- df.data[!is.na(df.data$X),]

##########################
### Linear de-trending ###
##########################

# Same as in Case_Study_part1.R, see that file for more details

co2 <- read.csv(paste0(path.to.folder,"Data/co2_mm_mlo.csv"))
co2_pre_1958 <- read.csv(paste0(path.to.folder,"Data/co2_pre_1958.csv"))

unique.years <- unique(co2$year)
co2.year.mean <- sapply(unique.years,function(xx) mean(co2$average[co2$year==xx]))
co2.year.mean <- cbind(year=unique.years,co2=co2.year.mean)
co2.year.mean <- rbind(co2_pre_1958,co2.year.mean)

df.data$co2 <- sapply(df.data$year, function(xx) co2.year.mean[co2.year.mean[,1]==xx,2])
df.data$co2.mod <- log(df.data$co2/280)

lin.fit <- lm(X ~ 1+co2.mod,data=df.data)
coefs <- lin.fit$coefficients
summary(lin.fit)

# Detrended series
df.data$Y <- df.data$X-(coefs[1]+df.data$co2.mod*coefs[2])

###############################################
### Asymptotic Independence diagnostic plot ###
###############################################

n <- length(df.data$Y)

# Lagged (1,2,3,4) time series
l1 <- cbind(df.data$Y[1:(n-1)],
            df.data$Y[2:(n)])
l2 <- cbind(df.data$Y[1:(n-2)],
            df.data$Y[3:(n)])
l3 <- cbind(df.data$Y[1:(n-3)],
            df.data$Y[4:(n)])
l4 <- cbind(df.data$Y[1:(n-4)],
            df.data$Y[5:(n)])

# chi-plots
par(mfrow=c(1,4),mar=c(4,4,1.5,0.3),mgp=c(2.3,1,0),pty="m")
chiplot(l1,which=1,main1 = "Lag 1")
abline(h=0)
chiplot(l2,which=1,main1 = "Lag 2")
abline(h=0)
chiplot(l3,which=1,main1 = "Lag 3")
abline(h=0)
chiplot(l4,which=1,main1 = "Lag 4")
abline(h=0)

# m-lag scatter plots
lag <- 12
par(mfrow=c(3,ceiling(lag/3)),mar=c(4,4,1.5,0.3),mgp=c(2.3,1,0),pty="m")
n <- length(df.data$Y)
for(l in 1:lag){
  plot(df.data$Y[1:(n-l)],df.data$Y[(1+l):n],
       main=paste0("Lag l=",l),xlab=expression(Y[t]),ylab=expression(Y[t+l]))
}

#######################
### Standardization ###
#######################

# Create dataframe with dates and detrended observations X
df <- data.frame(dt       = df.data$date,
                 time.brk = df.data$year,
                 X        = df.data$Y)

# Convert data to Laplace margins with semi-parametric model
thresh.conv.lap <- quantile(df$X,0.90) # u above which rank transform becomes parametric GPD
LapTransf <- ToLaplace_df(df=df, thresh = thresh.conv.lap)

# Select same quantile level as in Part 1
q <- 0.94
q.Y.star <- quantile(LapTransf$df$X,q)

par(mfrow=c(1,1),mar=c(6,6,0.2,0.2),mgp=c(4,1.6,0),pty="m")
plot(as.Date(df.data$dt),LapTransf$df$X,xlab = "Date", ylab= "Y_t",
     cex = 0.6, pch = 20,cex.axis = 1, cex.lab = 1,col="grey70")
points(as.Date(df.data$dt)[LapTransf$df$X>q.Y.star],
       LapTransf$df$X[LapTransf$df$X>q.Y.star],cex = 0.6, pch = 20)

# Getting rough idea of cluster sizes
rr <- 2
clust.sizes <- unlist(lapply(clusters(data=df$X,u=q.Y.star,r=rr),length))
hist(clust.sizes)
sum(clust.sizes>7)
length(clust.sizes)

################################################################################
################################### Case Study 2 ###############################
################################################################################

#############################################
### Model selection via cross-validation ####
#############################################

# L
L <- 6

# Form of function h
norming <- "classic"

# Number of forward simulations
n.sim <- 5000

#Filter the data to get blocks of L observations following an exceedance
filt.data <- filterData(df=LapTransf$df, L=L, u=q.Y.star)

# Indices of start of folds for cross validation
start.ind <- seq(1,length(filt.data),by=10)

# Store probability integral transformed observations
U.list <- lapply(0:L, function(kk) matrix(NA, ncol=3,nrow=length(filt.data)))

# Store alphas from model fits
alphas <- list(alphas1 = matrix(NA, ncol=L,nrow=length(start.ind)),
               alphas2 = matrix(NA, ncol=L,nrow=length(start.ind)),
               alphas3 = matrix(NA, ncol=L,nrow=length(start.ind)))

# Store betas from model fits
betas <- list(beta1 = rep(NA, length(start.ind)),
              beta2 = rep(NA, length(start.ind)),
              beta3 = rep(NA, length(start.ind)))

# Store simultaneous prediction intervals coverage
simPred.cover <- matrix(NA,ncol=3,nrow=length(filt.data))

cnt <- 1
set.seed(444)
# Loop over the start indices of the validation sets
for(i in 1:length(start.ind)){
  print(start.ind[i])
  # Isolate validation data
  valid.data <- filt.data[start.ind[i]:min((start.ind[i]+9),length(filt.data))]

  # Isolate training data
  fit.data   <- filt.data[-(start.ind[i]:min((start.ind[i]+9),length(filt.data)))]

  # Fit the models of order 1, 2, and 3
  fitted.mod1 <- fit.Markov(data=fit.data,L=L,pars=c(1,0),
                            norming=norming,orderMarkov=1)
  fitted.mod2 <- fit.Markov(data=fit.data,L=L,pars=c(1,1,0),
                            norming=norming,orderMarkov=2)
  fitted.mod3 <- fit.Markov(data=fit.data,L=L,pars=c(1.5,0.75,0.35,0),
                            norming=norming,orderMarkov=3)

  # Store alpaha an beta values
  alphas[[1]][i,] <- fitted.mod1$alphas
  alphas[[2]][i,] <- fitted.mod2$alphas
  alphas[[3]][i,] <- fitted.mod3$alphas
  betas[[1]][i] <- fitted.mod1$beta
  betas[[2]][i] <- fitted.mod2$beta
  betas[[3]][i] <- fitted.mod3$beta

  # for each block of observation in the validation set
  for(j in 1:length(valid.data)){

    # Perform forward simulation from each model
    frwrd.pred1 <- forwardSim(fitted.mod1,fit.data, q.Y.star, sim.t_0=">",n.sim)
    frwrd.pred1 <- do.call(rbind,frwrd.pred1)
    frwrd.pred2 <- forwardSim(fitted.mod2,fit.data, q.Y.star, sim.t_0=">",n.sim)
    frwrd.pred2 <- do.call(rbind,frwrd.pred2)
    frwrd.pred3 <- forwardSim(fitted.mod3,fit.data, q.Y.star, sim.t_0=">",n.sim)
    frwrd.pred3 <- do.call(rbind,frwrd.pred3)

    # Obtain simultaneous prediction intervals for the forward simulations
    simPred1 <- cluster.simPredict(frwrd.pred1)
    simPred2 <- cluster.simPredict(frwrd.pred2)
    simPred3 <- cluster.simPredict(frwrd.pred3)

    # Store coverage of the simultaneous prediction intervals
    C1 <- ifelse(sum(valid.data[[j]]<simPred1[,1] | valid.data[[j]]>simPred1[,2])>0,1,0)
    C2 <- ifelse(sum(valid.data[[j]]<simPred2[,1] | valid.data[[j]]>simPred2[,2])>0,1,0)
    C3 <- ifelse(sum(valid.data[[j]]<simPred3[,1] | valid.data[[j]]>simPred3[,2])>0,1,0)
    simPred.cover[cnt,] <- c(C1,C2,C3)

    # Probability Integral Transform of the observation from forward simulation
    for(l in 0:L){
      F.t <-ecdf(frwrd.pred1[,l+1]); U1 <- F.t(valid.data[[j]][l+1])
      F.t <-ecdf(frwrd.pred2[,l+1]); U2 <- F.t(valid.data[[j]][l+1])
      F.t <-ecdf(frwrd.pred3[,l+1]); U3 <- F.t(valid.data[[j]][l+1])

      # Store PITs in U.list
      U.list[[l+1]][cnt,] <- c(U1,U2,U3)
    }
    cnt <- cnt + 1
  }
}

# Save mean and sd of cross-validation alpha estimates
mean.alphas1 <- round(sapply(1:L,function(lag) mean(alphas[[1]][,lag])),3)
sd.alphas1 <- round(sapply(1:L,function(lag) sd(alphas[[1]][,lag])),3)

mean.alphas2 <- round(sapply(1:L,function(lag) mean(alphas[[2]][,lag])),3)
sd.alphas2 <- round(sapply(1:L,function(lag) sd(alphas[[2]][,lag])),3)

mean.alphas3 <- round(sapply(1:L,function(lag) mean(alphas[[3]][,lag])),3)
sd.alphas3 <- round(sapply(1:L,function(lag) sd(alphas[[3]][,lag])),3)

# Save mean and sd of cross-validation beta estimates
mean.betas1 <- round(mean(betas[[1]]),3)
sd.betas1 <-  round(sd(betas[[1]]),3)

mean.betas2 <- round(mean(betas[[2]]),3)
sd.betas2 <-  round(sd(betas[[2]]),3)

mean.betas3 <- round(mean(betas[[3]]),3)
sd.betas3 <-  round(sd(betas[[3]]),3)

par(mfrow=c(1,3),mar=c(6,6,0.2,0.2),mgp=c(4,1.6,0),pty="s")
lag <- 4

# PP plots of the forward simulations at 3rd lag ahead
for(mod in 1:3){
  Us <- U.list[[lag]][,mod]; Us <- Us[!is.na(Us)]
  plot(c(1:length(Us))/(length(Us)+1),sort(Us),
       cex = 0.8, pch = 20,col="grey35",
       xlab = "Empirical probability", ylab= "Model probability")
  abline(a=0,b=1)
}

# PIT histograms
for(mod in 1:3){
  Us <- U.list[[lag]][,mod]; Us <- Us[!is.na(Us)]
  hist(Us,main="",
       xlab="Model probability",
       ylab="Frequency")
}

# Coverage study of path X_{t+1}, ..., X_{t+6} (acknowledge sampling variability)
round(1-apply(simPred.cover,2,mean),3)

#################################
### Model fitting on all data ###
#################################

# Filter data as before
filt.data <- filterData(df=LapTransf$df, L=L, u=q.Y.star)

# Model fitting
norming <- "classic"
fitted.mod1 <- fit.Markov(data=filt.data,L=L,pars = c(1,0),norming=norming,orderMarkov=1)
fitted.mod2 <- fit.Markov(data=filt.data,L=L,pars = c(1,1,0),norming=norming,orderMarkov=2)
fitted.mod3 <- fit.Markov(data=filt.data,L=L,pars = c(1.5,0.75,0.35,0),norming=norming,orderMarkov=3)

##########################
### Forward simulation ###
##########################

# Number of forward simulations to run given q.Y.star > u
n.sim <- 5000

# Select a threshold above which to simulate q.Y.star with standard exp.

# Run the forward simulation based on q.Y.star, fitted alphas and beta, residuals
set.seed(444)
clusters1 <- forwardSim(fitted.mod1,filt.data, q.Y.star, sim.t_0=">",n.sim)
clusters2 <- forwardSim(fitted.mod2,filt.data, q.Y.star, sim.t_0=">",n.sim)
clusters3 <- forwardSim(fitted.mod3,filt.data, q.Y.star, sim.t_0=">",n.sim)

# Plot forward simulations on Laplace scale
par(mfrow=c(1,1),mar=c(4.1,5.1,0.5,1.1),mgp=c(2.6,0.8,0))
plot(clusters1[[1]],type="l",ylim=c(min(unlist(clusters1)),max(unlist(clusters1))))
for(i in 2:n.sim){
  lines(clusters1[[i]])
}

# Find marginal quantiles of the forward simulations at all lags
quantsLap1 <- apply(do.call(rbind,clusters1),2,function(xx) quantile(xx,c(0.025,0.25,0.5,0.75,0.975)))
quantsLap2 <- apply(do.call(rbind,clusters2),2,function(xx) quantile(xx,c(0.025,0.25,0.5,0.75,0.975)))
quantsLap3 <- apply(do.call(rbind,clusters3),2,function(xx) quantile(xx,c(0.025,0.25,0.5,0.75,0.975)))

# Plot the marginal quantiles
plot(0:L,quantsLap1[1,],ylim=c(min(quantsLap1)*0.98,max(quantsLap1)*1.02),
     type="l",xlab="Lags ahead",ylab="Y_t")
for(i in 2:nrow(quantsLap1)){
  lines(0:L,quantsLap1[i,],type="l")
}
for(i in 1:nrow(quantsLap2)){
  lines(0:L,quantsLap2[i,],type="l",lty="dashed")
}
for(i in 1:nrow(quantsLap3)){
  lines(0:L,quantsLap3[i,],type="l",lty="dotted")
}


# Store the forward simulation in observations scale
CO2.t <- log(unique(df.data$co2[df.data$year==2023])/280)
frwrd.sim1 <- matrix(NA,nrow=n.sim,ncol=L+1)
frwrd.sim2 <- matrix(NA,nrow=n.sim,ncol=L+1)
frwrd.sim3 <- matrix(NA,nrow=n.sim,ncol=L+1)

# Backtransform the forward simulations from Laplace to original scale
for(i in 1:n.sim){
  frwrd.sim1[i,] <- BackTransform(x.lap = clusters1[[i]],
                                  y.star = LapTransf$thresh,
                                  rate = LapTransf$rate,
                                  sig = LapTransf$sig,
                                  xi = LapTransf$xi,
                                  x.orig = df.data$Y)
  frwrd.sim2[i,] <- BackTransform(x.lap = clusters2[[i]],
                                  y.star = LapTransf$thresh,
                                  rate = LapTransf$rate,
                                  sig = LapTransf$sig,
                                  xi = LapTransf$xi,
                                  x.orig = df.data$Y)
  frwrd.sim3[i,] <- BackTransform(x.lap = clusters3[[i]],
                                  y.star = LapTransf$thresh,
                                  rate = LapTransf$rate,
                                  sig = LapTransf$sig,
                                  xi = LapTransf$xi,
                                  x.orig = df.data$Y)

  frwrd.sim1[i,] <- frwrd.sim1[i,] + (coefs[1]+CO2.t*coefs[2])
  frwrd.sim2[i,] <- frwrd.sim2[i,] + (coefs[1]+CO2.t*coefs[2])
  frwrd.sim3[i,] <- frwrd.sim3[i,] + (coefs[1]+CO2.t*coefs[2])
}

# Plot forward simulations on observation scale
par(mfrow=c(1,1),mar=c(4.1,5.1,0.5,1.1),mgp=c(2.6,0.8,0))

# Compute marginal quantiles of forward simulations on original scale
quants1 <- apply(frwrd.sim1,2,function(xx) quantile(xx,c(0.025,0.25,0.5,0.75,0.975)))
quants2 <- apply(frwrd.sim2,2,function(xx) quantile(xx,c(0.025,0.25,0.5,0.75,0.975)))
quants3 <- apply(frwrd.sim3,2,function(xx) quantile(xx,c(0.025,0.25,0.5,0.75,0.975)))

# Plot marginal quantiles of forward simulations on original scale
plot(0:L,quants1[1,],ylim=c(min(quants1)*0.98,max(quants1)*1.02),
     type="l",xlab="Lags ahead",ylab="X_t")
for(i in 2:nrow(quants1)){
  lines(0:L,quants1[i,],type="l")
}
for(i in 1:nrow(quants2)){
  lines(0:L,quants2[i,],type="l",lty="dashed")
}
for(i in 1:nrow(quants2)){
  lines(0:L,quants3[i,],type="l",lty="dashed")
}


##########################
### Cluster statistics ###
##########################

# Estimated mean cluster size at quantile q
1/ExtrIndex_Intervals(df.data$Y,u=quantile(df.data$Y,q))

# Corresomping threshold on the original scale
u.2023 <- quantile(df.data$Y,q) + (coefs[1]+CO2.t*coefs[2])

# Estimate of mean cluster size from model 1
mean(apply(frwrd.sim1,1,function(xx) sum(xx>u.2023)))
sd(apply(frwrd.sim1,1,function(xx) sum(xx>u.2023)))

# Estimate of mean cluster size from model 2
mean(apply(frwrd.sim2,1,function(xx) sum(xx>u.2023)))
sd(apply(frwrd.sim2,1,function(xx) sum(xx>u.2023)))

# Estimate of mean cluster size from model 3
mean(apply(frwrd.sim3,1,function(xx) sum(xx>u.2023)))
sd(apply(frwrd.sim3,1,function(xx) sum(xx>u.2023)))

# pmf from 3 models
pmf1 <- cluster.size.pmf(frwrd.sim1, u.2023, L+1)
pmf2 <- cluster.size.pmf(frwrd.sim2, u.2023, L+1)
pmf3 <- cluster.size.pmf(frwrd.sim3, u.2023, L+1)

# pmf plots
plot(pmf1,pch=20,cex=3,
     xlab="Cluster size",ylab="PMF")
points(pmf2,pch=20,cex=3,col="grey40")
points(pmf3,pch=20,cex=3,col="grey80")

# cluster max functional
cluster.max(frwrd.sim1,plot.max=TRUE)

# Cluster mean functional
meanPred1 <- cluster.mean(frwrd.sim1,plot.mean=FALSE)
meanPred2 <- cluster.mean(frwrd.sim2,plot.mean=FALSE)
meanPred3 <- cluster.mean(frwrd.sim3,plot.mean=FALSE)

# Simultaneous predictive intervals from forward simulations
simPred1 <- cluster.simPredict(frwrd.sim1,plot.sim=FALSE)
simPred2 <- cluster.simPredict(frwrd.sim2,plot.sim=FALSE)
simPred3 <- cluster.simPredict(frwrd.sim3,plot.sim=FALSE)

# Plots of simultaenous predictive intervals
lwd <- 2.5
plot(0:L,simPred1[,1],ylim=range(c(simPred1,simPred2,simPred3)),type="l",
     xlab="Lags ahead",ylab="Y_t",lwd=lwd)
lines(0:L,simPred1[,2],lwd=lwd)
lines(0:L,simPred2[,1],lty="dashed",lwd=lwd)
lines(0:L,simPred2[,2],lty="dashed",lwd=lwd)
lines(0:L,simPred3[,1],lty="dotted",lwd=lwd)
lines(0:L,simPred3[,2],lty="dotted",lwd=lwd)
lines(0:L,meanPred1,lwd=lwd)
lines(0:L,meanPred2,lty="dashed",lwd=lwd)
lines(0:L,meanPred3,lty="dotted",lwd=lwd)
