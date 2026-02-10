rm(list=ls())

set.seed(123)

#########################################################
### Load functions for bivariate conditional extremes ###
#########################################################
source("BivariateCE - functions.R")

##############################################
### Read in data and remove missing values ###
##############################################
load("PollutionData_Bournemouth_Winter.RData")
air.pollution.cleaned <- na.omit(air.pollution)[,2:3]

####################################
### Transform to Laplace margins ###
####################################
air.pollution.Laplace <- apply(air.pollution.cleaned, 2, laplace.margins)

####################################################################
####################################################################
### Figure 1: Plots on original scale and standard Laplace scale ###
####################################################################
####################################################################
par(mfrow=c(1,2), mar=c(5,5,2,2))
plot(air.pollution.cleaned, col="grey70", cex=0.75, pch=16,
     xlab=expression(NO), ylab=expression(NO[2]), 
     cex.lab=1.5, cex.axis=1.5)
plot(air.pollution.Laplace, col="grey70", cex=0.75, pch=16, asp=1, 
     xlab=expression(Y[1]), ylab=expression(Y[2]), 
     cex.lab=1.5, cex.axis=1.5)

### Add threshold and highlight extreme points ###
u <- qlaplace(0.95)
abline(v=u, lwd=3, lty=2)
air.pollution.Laplace.extreme <- air.pollution.Laplace[air.pollution.Laplace[,1]>u,]
points(air.pollution.Laplace.extreme, cex=1, pch=17)

##################################
### Conditional extremes fit 1 ###
##################################
### With full negative log-likelihood ###
fitted.nll <- optim(rep(0.5,4), fn=CE.nll, data=air.pollution.Laplace.extreme)$par

### With profile negative log-likelihood ###
fitted.PL <- optim(rep(0.5,2), fn=CE.nll.PL, data=air.pollution.Laplace.extreme)$par

##################################
##################################
### Figure 2: Diagnostic plots ###
##################################
##################################
### Calculation of empirical residuals and threshold excesses from the data ###
Y1.ext <- air.pollution.Laplace.extreme[,1]
Y2.ext <- air.pollution.Laplace.extreme[,2]
a.hat <- fitted.PL[1]
b.hat <- fitted.PL[2]

residuals <- (Y2.ext - a.hat*Y1.ext)/Y1.ext^b.hat
excesses <- Y1.ext - u

### Diagnostic plot (i): checking independence of residuals and excesses ###
plot(excesses, residuals, col="grey70", pch=16,
     main="(i)", xlab="Threshold excesses", ylab="Empirical residuals", 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)

### Correlation between residuals and excesses (using Kendall's tau) ###
cor(excesses, residuals, method = "kendall")

### Non-parametric test of independence ###
library(testforDEP)
testforDEP(excesses, residuals, test="HOEFFD")

### Simulation of new observations from the fitted model ###
N <- 10*length(Y1.ext)
y1.sim <- u + rexp(N)
resid.sim <- sample(residuals, N, rep=T)
y2.sim <- a.hat*y1.sim + resid.sim*(y1.sim^b.hat)

### Diagnostic plot (ii): checking similarity between observed and simulated data ###
plot(y1.sim, y2.sim, col="grey70", pch=16,
     main="(ii)", xlab=expression(Y[1]), ylab=expression(Y[2]), 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(air.pollution.Laplace.extreme, cex=1.25, pch=17)

###################################################
###################################################
### Figure 4: Estimation of joint probabilities ###
###################################################
###################################################
### Aiming to estimate P(Y1>v, Y2>r) for some large v ###
### Set the value of v ###
v <- qlaplace(0.99)

### Simulation of new observations with Y1>v ###
M <- 10000
y1.sim <- v + rexp(M)
resid.sim <- sample(residuals, M, rep=T)
y2.sim <- a.hat*y1.sim + resid.sim*(y1.sim^b.hat)

### Plot to demonstrate estimation of P(Y1>v, Y2>r) ###
r <- qlaplace(0.99)
plot(air.pollution.Laplace, col="grey90", cex=0.75, pch=16, asp=1, 
     xlim=range(air.pollution.Laplace, y1.sim, y2.sim),
     ylim=range(air.pollution.Laplace, y1.sim, y2.sim),
     xlab=expression(Y[1]), ylab=expression(Y[2]), 
     cex.lab=1.5, cex.axis=1.5)
points(y1.sim, y2.sim, col="grey70", pch=16)
points(y1.sim[y2.sim>r], y2.sim[y2.sim>r], pch=16)
abline(v=v, lty=3, lwd=3)

### Estimation of P(Y1>v, Y2>r) across different values of r ###
rs <- qlaplace(seq(0.5,0.99999,by=0.00001))
probs <- NULL
for(i in 1:length(rs)){
  probs[i] <- 0.01*mean(y2.sim>rs[i])  ## 0.01 because v was set to the 0.99 quantile of Y1
}

### Compare to empirical estimates from the original data ###
rs.orig <- seq(0, 11, by=0.01)
probs.orig <- NULL
for(i in 1:length(rs.orig)){
  probs.orig[i] <- mean(air.pollution.Laplace[,1]>v & air.pollution.Laplace[,2]>rs.orig[i])
}

### Plot to compare estimates ###
plot(rs, probs, type="l", col="grey50", lwd=4,
     xlab="r", ylab=expression(P(Y[1] > v, Y[2] > r)), 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(rs.orig, probs.orig, pch=16)  

#####################
### Bootstrapping ###
#####################
### First checking for autocorrelation in the margins, to decide on bootstrapping scheme ###
#acf(air.pollution.Laplace[,1]) (Commented out here for plotting purposes)
#acf(air.pollution.Laplace[,2])

### Stationary bootstrap with average block length 10 ###
n.iter <- 1000
prob.boot <- matrix(NA, nrow=n.iter, ncol=length(rs))

for(iter in 1:n.iter){
  ### Generate bootstrap sample from original data and transform to Laplace margins ###
  n <- nrow(air.pollution.Laplace)
  boot.sample <- boot.stat(air.pollution.Laplace, 10)

  ### Estimate model parameters ###
  fitted.PL <- optim(rep(0.5,2), fn=CE.nll.PL, data=boot.sample[boot.sample[,1]>u,])$par
  a.hat <- fitted.PL[1]
  b.hat  <- fitted.PL[2]
  
  ### Estimate required probabilities ###
  Y1.ext <- boot.sample[boot.sample[,1]>u,1]
  Y2.ext <- boot.sample[boot.sample[,1]>u,2]
  residuals <- (Y2.ext - a.hat*Y1.ext)/Y1.ext^b.hat
  y1.sim <- v + rexp(M)
  resid.sim <- sample(residuals, M, rep=T)
  y2.sim <- a.hat*y1.sim + resid.sim*(y1.sim^b.hat)
  for(i in 1:length(rs)){
    prob.boot[iter,i] <- 0.01*mean(y2.sim>rs[i])  
  }
  print(iter)
}

### 95% bootstrapped CIs for probabilities of interest ###
prob.CI.lower <- apply(prob.boot, 2, quantile, 0.025)  
prob.CI.upper <- apply(prob.boot, 2, quantile, 0.975)  

### Add confidence intervals to Figure 4 ###
points(rs, prob.CI.lower, type="l", lwd=4, col="grey50", lty=3)
points(rs, prob.CI.upper, type="l", lwd=4, col="grey50", lty=3)

##########################################
##########################################
### Figure 5: Threshold selection plot ###
##########################################
##########################################
### Fix thresholds to test ###
thresholds <- qlaplace(seq(0.5,0.99,by=0.01))

### Model fit using negative profile log-likelihood for each threshold ###
alphas <- NULL
betas  <- NULL
for(iter in 1:length(thresholds)){
  u <- thresholds[iter]
  air.pollution.Laplace.extreme <- air.pollution.Laplace[air.pollution.Laplace[,1]>u,]
  fitted.PL <- optim(rep(0.5,2), fn=CE.nll.PL, data=air.pollution.Laplace.extreme)$par
  alphas[iter] <- fitted.PL[1]
  betas[iter]  <- fitted.PL[2]
}

### Bootstrapping for parameter estimates across different thresholds ###
alphas.boot <- matrix(NA, nrow=1000, ncol=length(thresholds))
betas.boot  <- matrix(NA, nrow=1000, ncol=length(thresholds))

for(b in 1:1000){
  n <- nrow(air.pollution.Laplace)
  boot.sample <- boot.stat(air.pollution.Laplace, 10)
  
  alph <- NULL
  beta <- NULL
  for(iter in 1:length(thresholds)){
    u <- thresholds[iter]
    boot.sample.extreme <- boot.sample[boot.sample[,1]>u,]
    fitted.PL <- optim(rep(0.5,2), fn=CE.nll.PL, data=boot.sample.extreme)$par
    alph[iter] <- fitted.PL[1]
    beta[iter] <- fitted.PL[2]
  }
  alphas.boot[b,] <- alph
  betas.boot[b,]  <- beta
}

### Threshold stability plots ###
plot(thresholds, alphas, pch=16,
     xlab="u", ylab=expression(hat(alpha)), ylim=range(alphas.boot),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
abline(v=1.514, col="grey80", lwd=3, lty=2)
abline(v=2.303, col="grey50", lwd=3, lty=2)
for(i in 1:length(thresholds)){
  points(rep(thresholds[i],2), quantile(alphas.boot[,i], c(0.025,0.975)), type="l")
}

plot(thresholds, betas, pch=16,
     xlab="u", ylab=expression(hat(beta)), ylim=c(-5,1),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
abline(v=1.514, col="grey80", lwd=3, lty=2)
abline(v=2.303, col="grey50", lwd=3, lty=2)
for(i in 1:length(thresholds)){
  points(rep(thresholds[i],2), quantile(betas.boot[,i], c(0.025,0.975)), type="l")
}

#############################################################
#############################################################
### Figure 6: Comparing results after threshold selection ###
#############################################################
#############################################################
### Calculate new probability estimates and 95% confidence intervals ###
u2 <- qlaplace(0.89)
air.pollution.Laplace.extreme <- air.pollution.Laplace[air.pollution.Laplace[,1]>u2,]

### Parameter estimation with profile negative log-likelihood ###
fitted.PL <- optim(rep(0.5,2), fn=CE.nll.PL, data=air.pollution.Laplace.extreme)$par

### Calculation of empirical residuals ###
Y1.ext <- air.pollution.Laplace.extreme[,1]
Y2.ext <- air.pollution.Laplace.extreme[,2]
a.hat <- fitted.PL[1]
b.hat <- fitted.PL[2]
residuals <- (Y2.ext - a.hat*Y1.ext)/Y1.ext^b.hat
excesses <- Y1.ext - u2

### Estimation of P(Y1>v, Y2>r) across different values of r ###
v <- qlaplace(0.99)
M <- 10000
y1.sim <- v + rexp(M)
resid.sim <- sample(residuals, M, rep=T)
y2.sim <- a.hat*y1.sim + resid.sim*(y1.sim^b.hat)

rs <- qlaplace(seq(0.5,0.99999,by=0.00001))
probs2 <- NULL
for(i in 1:length(rs)){
  probs2[i] <- 0.01*mean(y2.sim>rs[i])  
}

#####################
### Bootstrapping ###
#####################
n.iter <- 1000
prob.boot2 <- matrix(NA, nrow=n.iter, ncol=length(rs))

for(iter in 1:n.iter){
  ### Generate bootstrap sample from original data and transform to Laplace margins ###
  n <- nrow(air.pollution.Laplace)
  boot.sample <- boot.stat(air.pollution.Laplace, 10)
  
  ### Estimate model parameters ###
  fitted.PL <- optim(rep(0.5,2), fn=CE.nll.PL, data=boot.sample[boot.sample[,1]>u2,])$par
  a.hat <- fitted.PL[1]
  b.hat  <- fitted.PL[2]
  
  ### Estimate required probabilities ###
  Y1.ext <- boot.sample[boot.sample[,1]>u2,1]
  Y2.ext <- boot.sample[boot.sample[,1]>u2,2]
  residuals <- (Y2.ext - a.hat*Y1.ext)/Y1.ext^b.hat
  y1.sim <- v + rexp(M)
  resid.sim <- sample(residuals, M, rep=T)
  y2.sim <- a.hat*y1.sim + resid.sim*(y1.sim^b.hat)
  for(i in 1:length(rs)){
    prob.boot2[iter,i] <- 0.01*mean(y2.sim>rs[i])  
  }
  print(iter)
}

### 95% bootstrapped CIs for probabilities of interest ###
prob.CI.lower2 <- apply(prob.boot2, 2, quantile, 0.025)  
prob.CI.upper2 <- apply(prob.boot2, 2, quantile, 0.975)  

### Plot to compare the two sets of estimates ###
par(mfrow=c(1,1))
plot(rs, probs, type="l", col="grey70", lwd=4,
     xlab="r", ylab=expression(P(Y[1] > v, Y[2] > r)), 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(rs, prob.CI.lower, type="l", lwd=4, col="grey70", lty=3)
points(rs, prob.CI.upper, type="l", lwd=4, col="grey70", lty=3)

points(rs, probs2, type="l", lwd=4)
points(rs, prob.CI.lower2, type="l", lwd=4, lty=3)
points(rs, prob.CI.upper2, type="l", lwd=4, lty=3)
