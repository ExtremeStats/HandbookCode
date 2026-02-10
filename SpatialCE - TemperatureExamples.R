rm(list=ls())

set.seed(123)

##############################
### Load required packages ###
##############################
library(mvnfast)
library(mvtnorm)
library(ggplot2)
library(maptools)
library(maps)
library(geosphere)

########################################################################################
########## Load specific functions for spatial conditional extremes approach ###########
### These files can be found as supplementary material for Wadsworth and Tawn (2022) ###
##### Download before use and save in file called "Functions" in working directory #####
########################################################################################
source("Functions/SpatialCEMod3-Functions.r")
source("Functions/DeltaLaplaceFunctions.r")
source("Functions/SpatialCE-CommonFunctions.r")
source("SpatialCE - functions.R")

####################
### Read in data ###
####################
load("TemperatureData_Netherlands_Summer.RData")

####################################
### Transform to Laplace margins ###
####################################
temps.Laplace <- apply(max.temps.summer, 2, laplace.margins)

#####################################
### Plot of country and locations ###
#####################################
netherlands <- map_data('world')[map_data('world')$region == "Netherlands",]
ggplot() + 
  geom_polygon(data = map_data("world"), aes(x=long, y=lat, group = group), color = 1, fill = "white") +
  geom_polygon(data = netherlands, aes(x=long, y=lat, group = group), color = "grey60", fill = "grey60") +
  coord_fixed(1.2, xlim = c(3, 7.5), ylim = c(50.7,54)) +
  geom_point(aes(x=locations[,1], y=locations[,2]), size=3) +
  labs(x = "Longitude", y = "Latitude", size=2)

##############################################################################
### Pairwise distance and estimates of chi(0.95), chi(0.975) and chi(0.99) ###
##############################################################################
dists        <- matrix(NA, ncol=nrow(locations), nrow=nrow(locations))
chi.ests.95  <- matrix(NA, ncol=nrow(locations), nrow=nrow(locations))
chi.ests.975  <- matrix(NA, ncol=nrow(locations), nrow=nrow(locations))
chi.ests.99 <- matrix(NA, ncol=nrow(locations), nrow=nrow(locations))

for(i in 1:nrow(locations)){
  for(j in 1:nrow(locations)){
    # Calculate Haversine distance between locations i and j
    dists[i,j] <- distHaversine(locations[i,], locations[j,])/1000  ## Distances in km
    
    # Extract data for locations i and j
    Y1 <- temps.Laplace[,i]
    Y2 <- temps.Laplace[,j]
    
    # Empirical estimates of chi(u)
    chi.ests.95[i,j]   <- chi.u(Y1, Y2, 0.95)
    chi.ests.975[i,j]  <- chi.u(Y1, Y2, 0.975)
    chi.ests.99[i,j]   <- chi.u(Y1, Y2, 0.99)
    
  }
}
par(mfrow=c(1,3), mar=c(5,5,2,2))
plot(dists, chi.ests.95,  xlab="Distance (km)", ylab="",
     pch=16, col="grey60", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5)
plot(dists, chi.ests.975,  xlab="Distance (km)", ylab="",
     pch=16, col="grey60", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5)
plot(dists, chi.ests.99, xlab="Distance (km)", ylab="",
     pch=16, col="grey60", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5)


#####################
### Model fitting ###
#####################
dists <- dists/100

### Start with single conditioning site near the centre of the spatial domain ###
### Start with nlminb fit and then repeat optim until convergence is achieved ###
fit.init <- nlminb(spatialCEnllMod3.singlesite,start=c(1.6,1.3,0.2,3,0,0.5,1,1),
                   x=temps.Laplace, 
                   coord=locations, DM=dists, thresh=qdlaplace(0.95,mu=0,sigma=1,delta=1),
                   condsite=32,
                   diffGauss = F,usevgm=F, GA=F, print=T,control=list(iter.max=5000,rel.tol=1e-6))
fit.init <- optim(spatialCEnllMod3.singlesite, par=fit.init$par,
                  x=temps.Laplace, 
                  coord=locations, DM=dists, thresh=qdlaplace(0.95,mu=0,sigma=1,delta=1),
                  condsite=32,
                  diffGauss = F,usevgm=F, GA=F, print=T,control=list(maxit=5000,reltol=1e-6))

### Use single conditioning site results as starting point for composite fit ###
### Run until it has actually converged - sometimes optim stops too early ###
fit1 <- optim(spatialCEnllMod3, par=fit.init$par,
              x=temps.Laplace, 
              coord=locations, DM=dists, thresh=qdlaplace(0.95,mu=0,sigma=1,delta=1),
              diffGauss = F,usevgm=F, GA=F, print=T,control=list(maxit=5000,reltol=1e-6))
fit1 <- optim(spatialCEnllMod3, par=fit1$par,
              x=temps.Laplace, 
              coord=locations, DM=dists, thresh=qdlaplace(0.95,mu=0,sigma=1,delta=1),
              diffGauss = F,usevgm=F, GA=F, print=T,control=list(maxit=5000,reltol=1e-6))

#save(fit1, file = "compositeFit3.RData")

################################
### Comparison of chi values ###
################################
#load(file = "compositeFit3.RData")
par = fit1$par
newdata <- rCondSite.new3(n=20000, coord=locations, u=qdlaplace(0.95,mu=0,sigma=1,delta=1),  
                         kappa=par[1], lambda=par[2], b=par[3],
                         phi=par[4], mu=par[5], delta=par[6],
                         nu=par[7], sig=par[8],DM=dists, site=1)
chi.model1 <- apply(newdata, 2, function(x){mean(x>=qdlaplace(0.95,mu=0,sigma=1,delta=1))})
newdata <- rCondSite.new3(n=20000, coord=locations, u=qdlaplace(0.975,mu=0,sigma=1,delta=1),  
                         kappa=par[1], lambda=par[2], b=par[3],
                         phi=par[4], mu=par[5], delta=par[6],
                         nu=par[7], sig=par[8],DM=dists, site=1)
chi.model2 <- apply(newdata, 2, function(x){mean(x>=qdlaplace(0.975,mu=0,sigma=1,delta=1))})
newdata <- rCondSite.new3(n=20000, coord=locations, u=qdlaplace(0.99,mu=0,sigma=1,delta=1),  
                         kappa=par[1], lambda=par[2], b=par[3],
                         phi=par[4], mu=par[5], delta=par[6],
                         nu=par[7], sig=par[8],DM=dists, site=1)
chi.model3 <- apply(newdata, 2, function(x){mean(x>=qdlaplace(0.99,mu=0,sigma=1,delta=1))})

par(mfrow=c(1,3), mar=c(5,5,2,2))
plot(100*dists, chi.ests.95, xlab="Distance (km)", ylab="",
     pch=16, col="grey60", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5)
points(100*dists[sort(dists[,1], ind=T)$ix,1], chi.model1[sort(dists[,1], ind=T)$ix], 
       type="l", lwd=5)
plot(100*dists, chi.ests.975, xlab="Distance (km)", ylab="",
     pch=16, col="grey60", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5)
points(100*dists[sort(dists[,1], ind=T)$ix,1], chi.model2[sort(dists[,1], ind=T)$ix], 
       type="l", lwd=5)
plot(100*dists, chi.ests.99, xlab="Distance (km)", ylab="",
     pch=16, col="grey60", ylim=c(0,1), cex.lab=1.5, cex.axis=1.5)
points(100*dists[sort(dists[,1], ind=T)$ix,1], chi.model3[sort(dists[,1], ind=T)$ix], 
       type="l", lwd=5)

##########################
### Comparison of data ###
##########################
set.seed(123)
par(mfrow=c(1,2))
transect <- which(locations[,2]==52.375)
newdata <- rCondSite.new3(n=128, coord=locations, u=qdlaplace(0.95,mu=0,sigma=1,delta=1),  
                         kappa=par[1], lambda=par[2], b=par[3],
                         phi=par[4], mu=par[5], delta=par[6],
                         nu=par[7], sig=par[8],DM=dists, site=transect[1])

plot(locations[transect,1], newdata[1,transect], type="l", ylim=range(newdata), 
     xlab="Longitude", ylab="Y(s)",
     cex.lab=1, cex.axis=1)
for(i in 2:nrow(newdata)){
  points(locations[transect,1], newdata[i,transect], type="l")
}
ext <- which(temps.Laplace[,transect[1]]>qdlaplace(0.95,mu=0,sigma=1,delta=1))
plot(locations[transect,1], temps.Laplace[ext[1],transect], type="l", col="grey60", ylim=range(newdata),
     xlab="Longitude", ylab="Y(s)",
     cex.lab=1, cex.axis=1)
for(i in ext){
  points(locations[transect,1], temps.Laplace[i,transect], type="l", col="grey60")
}
