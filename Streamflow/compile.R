# This code compiles the results across the three models and makes plots and tables


rm(list=ls())
library(ggplot2)
library(tidyverse)

setwd("S:\\Documents\\My Papers\\EVA_Handbook\\BayesEVA\\Streamflow\\")
file <- "https://www4.stat.ncsu.edu/~bjreich/ST740/HCDN_annual_max.RData"
load(url(file))

y          <- log(Y+1)
ns         <- nrow(y)
nt         <- ncol(y)
x          <- (1:nt)/10 
x          <- x - mean(x)

load(paste0("model1.RData"))
  output1 <- output
load(paste0("model2.RData"))
  output2 <- output
load(paste0("model3.RData"))
  output3 <- output

# Tables
 print(output1$WAIC)
 print(output2$WAIC)
 print(output3$WAIC)

 tab <- output1$mu[,c(2,4,5)]
 tab <- cbind(tab,sqrt(output1$S1[c(1,6,11,16),c(2,4,5)]))
 tab <- cbind(tab,sqrt(output1$S2[c(1,6,11,16),c(2,4,5)]))
 print(round(tab,2))

# Plots

pdf("streamflow_results.pdf")

 r <- range(c(output1$beta0[,1],output3$beta0[,1]))
 plot(output1$beta0[,1],output3$beta0[,1],pch=19,cex=.5,
      xlim=r,ylim=r,cex.axis=1.5,cex.lab=1.5,
      xlab="Hierarchical model",ylab="Separate model")
 abline(0,1)

 r <- range(c(output1$beta1[,1],output3$beta1[,1]))
 plot(output1$beta1[,1],output3$beta1[,1],pch=19,cex=.5,
      xlim=r,ylim=r,cex.axis=1.5,cex.lab=1.5,
      xlab="Hierarchical model",ylab="Separate model")
 abline(0,1)

 r <- range(c(output1$lsig[,1],output3$lsig[,1]))
 plot(output1$lsig[,1],output3$lsig[,1],pch=19,cex=.5,
      xlim=r,ylim=r,cex.axis=1.5,cex.lab=1.5,
      xlab="Hierarchical model",ylab="Separate model")
 abline(0,1)

 r <- range(c(output1$xi[,1],output3$xi[,1]))
 plot(output1$xi[,1],output3$xi[,1],pch=19,cex=.5,
      xlim=r,ylim=r,cex.axis=1.5,cex.lab=1.5,
      xlab="Hierarchical model",ylab="Separate model")
 abline(0,1)


 sd_ratio_beta0 <- output1$beta0[,3]/output3$beta0[,3]
 sd_ratio_beta1 <- output1$beta1[,3]/output3$beta1[,3]
 sd_ratio_lsig  <- output1$lsig[,3]/output3$lsig[,3]
 sd_ratio_xi    <- output1$xi[,3]/output3$xi[,3]
 sd_ratio       <- cbind(sd_ratio_beta0,sd_ratio_beta1,sd_ratio_lsig,sd_ratio_xi)
 boxplot(sd_ratio,
         outline=FALSE,cex.lab=1.5,cex.axis=1.5,
         xlab="",ylab="Ratio of posterior SD",
         names=c(expression(beta[0]),expression(beta[1]),expression(log(sigma)),expression(xi)))
 abline(1,0)

est   <- output1$beta1[,2]
sig   <- (output1$beta1[,4]>0) | (output1$beta1[,5]<0)
state <- map_data("state")
out   <- data.frame(lon=s[,1],lat=s[,2],sig=as.factor(sig),est=est)
lim   <- max(abs(est))*c(-1,1)

gg <- ggplot() +
  geom_map(
    data = state, map = state,
    aes(long, lat, map_id = region),
    color = "gray", fill = "white", size = 0.01
  ) +
  geom_point(
    data = out,
    aes(lon, lat, color = est, shape = sig),
    alpha = 1,size=2) + 
  coord_fixed()+theme_bw() +
  scale_colour_gradient2(limits=lim)+
  xlab("Longitude")+
  ylab("Latitiude")+
  labs(colour = expression(beta[1]))+
  labs(shape = "Significant",size=2)
print(gg)



out1 <- out[sig,]
out2 <- out[!sig,]
gg <- ggplot() +
  geom_map(
    data = state, map = state,
    aes(long, lat, map_id = region),
    color = "gray", fill = "white", size = 0.01
  ) +
  geom_point(
    data = out2,
    aes(lon, lat, fill=est),
    alpha = 1,size=2,stroke=0.5,shape=21) + 
  geom_point(
    data = out1,
    aes(lon, lat, fill=est),
    alpha = 1,size=2,stroke=0.5,shape=24) + 
  coord_fixed()+theme_bw() +
  scale_colour_gradient2(limits=lim)+
  scale_fill_gradient2(limits=lim)+
  xlab("Longitude")+
  ylab("Latitiude")+
  labs(fill = expression(beta[1]))
print(gg)



dev.off()


