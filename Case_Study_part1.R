# Clear R environment
rm(list=ls())

# Specify path to the Ch14_TSExtremes/ folder
path.to.folder <- ".../Ch14_TSExtremes/"

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

# Select points that exceed the 0.95-quantile of the observed series
q <- quantile(df.data$X,0.95)
dts.exc <- df.data$dt[df.data$X>q]
X.exc <- df.data$X[df.data$X>q]

# Plot TN and exceedances of 0.95-quantile
par(mfrow=c(1,1),mar=c(6,6,0.2,0.2),mgp=c(4,1.6,0),pty="m")
plot(as.Date(df.data$dt),df.data$X,xlab = "XXXXX", ylab= "YYYYYYYYYYYY",
     cex = 0.6, pch = 20,col="grey70",cex.axis = 1, cex.lab = 1)
points(as.Date(dts.exc),X.exc,cex = 0.6, pch = 20,col="grey20")

##########################
### Linear de-trending ###
##########################

# Read the CO_2 .csv files
co2 <- read.csv(paste0(path.to.folder,"Data/co2_mm_mlo.csv"))
co2_pre_1958 <- read.csv(paste0(path.to.folder,"Data/co2_pre_1958.csv"))

# Find yearly mean CO_2
unique.years <- unique(co2$year)
co2.year.mean <- sapply(unique.years,function(xx) mean(co2$average[co2$year==xx]))

# Bind two sources of CO_2 data together
co2.year.mean <- cbind(year=unique.years,co2=co2.year.mean)
co2.year.mean <- rbind(co2_pre_1958,co2.year.mean)

# Add CO_2 data to the df.data dataframe
df.data$co2 <- sapply(df.data$year, function(xx) co2.year.mean[co2.year.mean[,1]==xx,2])

# Add log-scaled CO_2 covariate to the df.data dataframe
df.data$co2.mod <- log(df.data$co2/280)

par(mfrow=c(1,1),mar=c(6,6,0.2,6),mgp=c(4.5,1.6,0),pty="m")
plot(co2.year.mean$year,co2.year.mean$co2,xlab = "Date", ylab= "CO_2 (ppm)",
     type="l",# cex = 0.6, pch = 20,
     cex.axis = 1, cex.lab = 1)
par(new=TRUE)
plot(co2.year.mean$year,log(co2.year.mean$co2/280),xlab = "Date", ylab= "CO_2 (ppm)",
     type="l",axes=FALSE, col="red")
mtext("log(CO_2/280)",side=4,col="red",line=5)
axis(4, ylim=range(log(co2.year.mean$co2/280)), col="red",col.axis="red",las=3)

# ACF plot of the original time series X
acf(df.data$X,lag.max = 200)

# Fit linear model with co2.mod covariate
lin.fit <- lm(X ~ 1+co2.mod,data=df.data)
coefs <- lin.fit$coefficients
summary(lin.fit)

# Create the detrended variable Y (\widetilde{x} in chapter)
df.data$Y <- df.data$X-(coefs[1]+df.data$co2.mod*coefs[2])

par(mfrow=c(1,1),mar=c(6,6,0.2,0.2),mgp=c(4.5,1.6,0),pty="m")
plot(as.Date(df.data$dt),df.data$X,xlab = "Date", ylab= "Y_t (Celcius)",
     cex = 0.6, pch = 20,col="grey40",cex.axis = 1, cex.lab = 1)
lines(as.Date(df.data$dt),coefs[1]+df.data$co2.mod*coefs[2],col="red")

# Find the 0.95-quantile of the detrended series
q <- quantile(df.data$Y,0.95)
dts.exc <- df.data$dt[df.data$Y>q]

# Find exceedances of the 0.95-quantile
Y.exc <- df.data$Y[df.data$Y>q]

par(mfrow=c(1,1),mar=c(6,6,0.2,0.2),mgp=c(4.5,1.6,0),pty="m")
plot(as.Date(df.data$dt),df.data$Y,xlab = "Date", ylab= "\\widetilde{Y}_t",
     cex = 0.6, pch = 20,col="grey70")
points(as.Date(dts.exc),Y.exc,cex = 0.6, pch = 20,col="grey20")

# Produce ACF plot with Chapter template
ACF <- acf(df.data$Y,lag.max=31,xlab="XXXXX",ylab="YYYYYYYYYYYY",main="",plot = F)
get_clim <- function(x, ci=0.95, ci.type="white"){
  #' Gets confidence limit data from acf object `x`
  if (!ci.type %in% c("white", "ma")) stop('`ci.type` must be "white" or "ma"')
  if (class(x) != "acf") stop('pass in object of class "acf"')
  clim0 <- qnorm((1 + ci)/2) / sqrt(x$n.used)
  if (ci.type == "ma") {
    clim <- clim0 * sqrt(cumsum(c(1, 2 * x$acf[-1]^2)))
    return(clim[-length(clim)])
  } else {
    return(clim0)
  }
}
CI <- get_clim(ACF)

# ACF plot of the detrended time series Y
par(mfrow=c(1,1),mar=c(6,6,0.2,0.2),mgp=c(4.5,1.6,0),pty="m")
plot(ACF$acf,type="h",xlab="Lags",ylab="ACF")
abline(h=c(0,CI,-CI),lty=c("solid","dashed","dashed"))

##########################
### Extremal index est ###
##########################

# import the evd package for extremal index + GEV and GP model fits
library(evd)

# Define lower and upper limits for a sequence of thresholds u
tlim <- quantile(df.data$Y, probs = c(0.8,0.999))

# Intervals estimator
thetas <- exiplot(df.data$Y, tlim, r = 0)

# Runs estimator
rs <- c(1,2,3,4,5,10,20) # Various run lengths
thetas.matrix <- matrix(NA,ncol=length(rs),nrow=length(thetas$x))
for(i in 1:length(rs)){
  thetas.matrix[,i] <- exiplot(df.data$Y, tlim, r = rs[i], add = TRUE, lty = 2)$y
}

# Plot estimates of \theta(u,r)
plot(thetas$x,thetas$y,xlab = "Exceedance threshold u", ylab= "Theta",
     ylim=c(0,1),type="l",cex.axis = 1, cex.lab = 1)
abline(v=quantile(df.data$Y, probs = .9975),lty="dotted")
for(i in 1:length(rs)){
  lines(thetas$x,thetas.matrix[,i],lty=2)
}

# Plot estimates of 1/\theta(u,r)
plot(thetas$x,1/thetas$y,ylim=c(0,6),
     xlab = "Exceedance threshold u", ylab= "1/Theta",type="l",
     cex.axis = 1, cex.lab = 1)
abline(v=quantile(df.data$Y, probs = .9975),lty="dotted")
for(i in 1:length(rs)){
  lines(thetas$x,1/thetas.matrix[,i],lty=2)
}

###################
### GEV fitting ###
###################

# Isolate yearly maxima
year.max <- sapply(unique(df.data$year),function(xx) max(df.data$Y[df.data$year==xx]))

# Find dates associated with yearly maxima
dt.max <- sapply(unique(df.data$year),function(xx) df.data$dt[df.data$year==xx][which.max(df.data$Y[df.data$year==xx])])

# Plot yearly maxima
plot(as.Date(dt.max),year.max,xlab = "Date", ylab= "Yearly maxima of wildetilde{y}_t",
     cex = 0.8, pch = 20,col="grey35")

# Fit the gev model using the evd package
fitted.gev <- fgev(year.max)

# Extract the location, scale, and shape parameters
pars <- fitted.gev$estimate
mu <- pars[1]; sig <- pars[2]; xi <- pars[3]

# Find model quantiles
U <- c(1:length(year.max))/(length(year.max)+1)
z <- qgev(U,mu,sig,xi)

# Perform QQ plot
plot(sort(year.max),z,cex = 0.8, pch = 20,col="grey35",
     xlab = "Empirical", ylab= "Model")
abline(a=0,b=1)
Ulow <- sapply(1:length(z), function(i) {
  qbeta(0.025, i, length(z) + 1 - i)
})
Uup <- sapply(1:length(z), function(i) {
  qbeta(0.975, i, length(z) + 1 - i)
})
lines(z, qgev(Ulow,mu,sig,xi), lty = 2)
lines(z, qgev(Uup,mu,sig,xi), lty = 2)

# Obtain profile likelihood estimates of the parameters
nms <- names(fitted.gev$estimate)
prof.pars <- profile(fitted.gev,which=nms,
                     mesh = fitted.gev$std.err[nms]/100)
q <- qchisq(0.95,1)

par(mfrow=c(1,3),mar=c(6,6,1,0.2),mgp=c(4,1.6,0),pty="m")
mat <- prof.pars$loc
plot(mat[,1],-mat[,2],type="l",xaxt = "n",yaxt = "n",
     xlab="mu",ylab="Profile log-likelihood")
axis(1, at = seq(5.2,6.8,by=0.4))
axis(2, at = c(-254,-250,-246))
abline(v=mu,lty="dotted")
abline(h=max(-mat[,2])-q,lty="dashed")

mat <- prof.pars$scale
plot(mat[,1],-mat[,2],type="l",xaxt = "n",yaxt = "n",
     xlab="sigma",ylab="")
ran <- c(floor(min(mat[,1])*10),ceiling(max(mat[,1])*10))/10
axis(1, at = seq(ran[1],ran[2],by=0.4))
axis(2, at = c(-254,-250,-246))
abline(v=sig,lty="dotted")
abline(h=max(-mat[,2])-q,lty="dashed")

mat <- prof.pars$shape
plot(mat[,1],-mat[,2],type="l",xaxt = "n",yaxt = "n",
     xlab="xi",ylab="")
ran <- c(floor(min(mat[,1])*10),ceiling(max(mat[,1])*10))/10
axis(1, at = c(-0.5,-0.3,-0.1,0.1))
axis(2, at = c(-254,-250,-246))
abline(v=xi,lty="dotted")
abline(h=max(-mat[,2])-q,lty="dashed")

#####################
### Return levels ###
#####################

# Define a range of return periods
ret.T <- c(2:10,seq(20,45,by=5),
           seq(50,90,by=10),
           seq(100,1000,by=50))

# Find the associated probability of a return period
ps <- 1-1/ret.T

# Obtain the return level associated with a probability p
z.theta <- qgev(ps,mu,sig,xi)

# Estimate of the upper bound of the marginal distribution
up.bound <- mu-sig/xi

# Get the profile likelihood cinfi
CIs.z.theta <- matrix(NA,ncol=2,nrow=length(z.theta))
for(i in 1:length(ps)){
  p <- ps[i]; print(ret.T[i])
  M <- fgev(year.max, prob = 1-p)
  prof.z <- profile(M, which = "quantile",mesh=0.01)
  prof.z$quantile[which.max(-prof.z$quantile[,2]),1]
  dev.thresh <- max(-prof.z$quantile[,2])-q
  ind.CI <- range(which(-prof.z$quantile[,2]>dev.thresh))
  ind.CI <- c(ind.CI[1]-1,ind.CI[2]+1)
  CIs.z.theta[i,] <- prof.z$quantile[ind.CI,1]
}

# Plot return levels and profile likelihood confidence intervals
par(mfrow=c(1,1),mar=c(6,6,1,0.2),mgp=c(4,1.6,0),pty="m")
plot(ret.T,z.theta,type="l",
     ylim=c(min(z.theta),1.02*max(CIs.z.theta)),
     xlab="Return period T (years)",
     ylab="Return level z_T")
abline(h=up.bound,lty="dotted")
lines(ret.T,CIs.z.theta[,1],lty="dashed")
lines(ret.T,CIs.z.theta[,2],lty="dashed")

##################
### POT Models ###
##################

### Declustering and parameter stability ###

# Define a sequence of thresholds u
us <- seq(0.90,0.9975,by=0.001)

# Specify a range of declustering run legnths
rs <- c(1,2,5,10)

# Define matrix to storeparameter estimates
pars <- matrix(NA,ncol=4,nrow=length(us)*length(rs))

cnt <- 1
# For each run length in rs, perform inner loop
for(j in 1:length(rs)){
  # Decluster and fit GP distribution for each threshold u
  for(i in 1:length(us)){
    u <- quantile(df.data$Y,us[i])
    clusts <- clusters(data=df.data$Y,u=u,r=rs[j])
    max.clust <- unlist(lapply(clusts,function(xx) max(xx)))
    pars[cnt,] <- c(fpot(unlist(max.clust),threshold=u)$estimate,length(max.clust),rs[j])
    cnt <- cnt +1
  }
}

par(mfrow=c(1,3),mar=c(6,6,1,0.2),mgp=c(4,1.6,0),pty="m")
# Plot sigma estimates as function of r and u
plot(us,pars[pars[,4]==1,1],type="l",ylim=range(pars[,1]),xlab="u",ylab="sigma")
for(r in rs[-1]){
  lines(us,pars[pars[,4]==r,1])
}
# Plot xi estimates as function of r and u
plot(us,pars[pars[,4]==1,2],type="l",xlab="u",ylab="xi")
for(r in rs[-1]){
  lines(us,pars[pars[,4]==r,2])
}
# Plot number of exceedances as function of r and u
plot(us,pars[pars[,4]==1,3],type="l",xlab="u",ylab="number exceedances")
for(r in rs[-1]){
  lines(us,pars[pars[,4]==r,3])
}

# Based on the stability plots above,
# fit model with exceedances of 0.94-quantile of Y and r=5
u <- quantile(df.data$Y,0.94)

# Decluster the data using a run length r=5
clusts <- clusters(data=df.data$Y,u=u,r=5)

# Find cluster maxima
max.clust <- lapply(clusts,function(xx) max(xx))
ind.max.clust <- lapply(clusts,function(xx) as.numeric(names(xx))[which.max(xx)])

# Find dates of cluster maxima
dt.max.clust <- df.data$dt[unlist(ind.max.clust)]
year.max <- sapply(unique(df.data$year),function(xx) max(df.data$Y[df.data$year==xx]))
dt.max <- sapply(unique(df.data$year),function(xx) df.data$dt[df.data$year==xx][which.max(df.data$Y[df.data$year==xx])])

# Find which observations correspond to both a cluster and an annual maxima
ind.both <- which(dt.max.clust %in% dt.max)
both.max <- max.clust[ind.both]
both.dt <- dt.max.clust[ind.both] # Dates of both

# Plot annual maxima
par(mfrow=c(1,1),mar=c(6,6,1,0.2),mgp=c(4.5,1.6,0),pty="m")
plot(c(as.Date(df.data$dt[df.data$Y>u]),as.Date(dt.max)),
     c(df.data$Y[df.data$Y>u],year.max),
     xlab = "Date", ylab= "\\widetilde{X}_t",ylim=c(0.95*u,1.02*max(df.data$Y)),
     cex = 1, pch = 20,col="white")
points(as.Date(dt.max),year.max,cex = 2.5, pch = 20,col="red")
points(as.Date(both.dt),both.max,cex = 2.5, pch = 20)

# Plot cluster maxima
par(mfrow=c(1,1),mar=c(6,6,1,0.2),mgp=c(4.5,1.6,0),pty="m")
plot(c(as.Date(df.data$dt[df.data$Y>u]),as.Date(dt.max)),
     c(df.data$Y[df.data$Y>u],year.max),
     xlab = "Date", ylab= "\\widetilde{X}_t",ylim=c(0.95*u,1.02*max(df.data$Y)),
     cex = 1, pch = 20,col="white")
points(as.Date(dt.max.clust),max.clust,cex = 2.5, pch = 20,col="red")
points(as.Date(both.dt),both.max,cex = 2.5, pch = 20)

##################
### GP fitting ###
##################

# Fit GP model to cluster maxima
fitted.gp <- fpot(unlist(max.clust),threshold=u)

# Compute model quantiles
U <- c(1:length(unlist(max.clust)))/(length(unlist(max.clust))+1)
Model.q <- qgpd(U,u,fitted.gp$estimate[1],fitted.gp$estimate[2])

# QQ plot for the fitted model
par(mfrow=c(1,1),mar=c(6,6,1,0.2),mgp=c(4,1.6,0),pty="s")
plot(sort(unlist(max.clust)),Model.q,cex = 0.8, pch = 20,col="grey35",
     xlab = "Empirical", ylab= "Model")
abline(a=0,b=1)
Ulow <- sapply(1:length(unlist(max.clust)), function(i) {
  qbeta(0.025, i, length(unlist(max.clust)) + 1 - i)
})
Uup <- sapply(1:length(unlist(max.clust)), function(i) {
  qbeta(0.975, i, length(unlist(max.clust)) + 1 - i)
})
lines(Model.q, qgpd(Ulow,u,fitted.gp$estimate[1],fitted.gp$estimate[2]), lty = 2)
lines(Model.q, qgpd(Uup,u,fitted.gp$estimate[1],fitted.gp$estimate[2]), lty = 2)

# Save profile likelihood conficence intervals
CI.list <- list()

# Values of declustering run lengths
rs <- c(1,2,5,10)

# Compute return levels and profile likelihood confidence intervals
for(j in 1:length(rs)){ # Loop over possible declustering run lengths values
  print(rs[j])

  # Identify non-maxima within cluster
  clusts <- clusters(data=df.data$Y,u=u,r=rs[j])
  non.max.clust.list <- lapply(clusts,function(xx) as.numeric(names(xx[-which.max(xx)])))
  non.max.clust <- unlist(non.max.clust.list)

  # Remove non-maxima from data
  declust.df.data <- df.data[-non.max.clust,]

  # Fit GP distribution to declustered data

  CIs.z.gp <- matrix(NA,ncol=3,nrow=length(ret.T))
  for(i in 1:length(ret.T)){
    print(ret.T[i])
    fitted.gp <- fpot(declust.df.data$Y,threshold=u, mper = ret.T[i],npp=(2*31+30))
    prof.z_gp <- profile(fitted.gp, which = "rlevel", conf=0.95, mesh=0.001)

    prof.z_gp$rlevel[which.max(-prof.z_gp$rlevel[,2]),1]
    dev.thresh <- max(-prof.z_gp$rlevel[,2])-q
    ind.CI <- range(which(-prof.z_gp$rlevel[,2]>dev.thresh))
    ind.CI <- c(ind.CI[1]-1,ind.CI[2]+1)
    CIs.z.gp[i,] <- c(fitted.gp$estimate[1],prof.z_gp$rlevel[ind.CI,1])
  }

  fitted.gp <- fpot(declust.df.data$Y,threshold=u)
  z_0 <- u-fitted.gp$estimate[1]/fitted.gp$estimate[2]
  CI.list[[j]] <- list(CIs.z.gp=CIs.z.gp,z_0=z_0)
}

# Plot return levels and profile likelihood confidence intervals
par(mfrow=c(1,1),mar=c(6,6,1,0.2),mgp=c(4,1.6,0),pty="m")
plot(ret.T,z.theta,type="l",col="red",
     ylim=c(min(z.theta),1.02*max(c(CIs.z.theta,unlist(CI.list)))),
     xlab="Return period T (years)",
     ylab="Return level z_T",lwd=2)
lines(ret.T,CIs.z.theta[,1],col="red",lty="dashed",lwd=2)
lines(ret.T,CIs.z.theta[,2],col="red",lty="dashed",lwd=2)

cols <- seq(0.8,0.2,len=length(rs))
for(i in 1:length(CI.list)){
  lines(ret.T,CI.list[[i]]$CIs.z.gp[,1],col=rgb(cols[i],cols[i],cols[i]),lwd=2)
  lines(ret.T,CI.list[[i]]$CIs.z.gp[,2],lty="dashed",col=rgb(cols[i],cols[i],cols[i]),lwd=2)
  lines(ret.T,CI.list[[i]]$CIs.z.gp[,3],lty="dashed",col=rgb(cols[i],cols[i],cols[i]),lwd=2)
  segments(-100,CI.list[[i]]$z_0,0,CI.list[[i]]$z_0,col=rgb(cols[i],cols[i],cols[i]),lwd=2)#lty="dotted",
}
segments(-100,up.bound,0,up.bound,col="red",lwd=2)


