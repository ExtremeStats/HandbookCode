library(lubridate)
library(fGarch)
library(zoo)

#### Applying ARMA-GARCH filtering first 
load("Data/financial.RData")
data1=data1[,-c(2,12)]
missing=apply(data1,1,function(x) any(is.na(x)))
D=ncol(data1)-1
Returns=data.frame(as.Date(data1[-1,1]))
for(d in 1:D)
{
  Returns=cbind(Returns,as.numeric(-diff(log(data1[,1+d]))))
}

colnames(Returns)=colnames(data1)

N=nrow(Returns)
D=10
X0=U0=Z0=c()
par(mfrow=c(2,5),
    pty='s',mgp=c(2,1,0),omi=c(0,0,0,0)+0.1,mar=c(3,3,2,0),
    cex=1.2)
for(d in 1:D){
  fit3=  garchFit(~arma(1,0)+garch(1,1),Returns[,1+d],cond.dist='std') 
  nu=coef(fit3)['shape']
  se=fit3@fit$se.coef['shape']
  plot(fit3,which=13,main='',pch=20)
  text(-2.5,4,colnames(Returns)[d+1])
  text(-2.5,3,bquote(hat(nu)==.(round(nu,2))~'('~.(round(se,2))~')'))
  X0=cbind(X0,residuals(fit3))
  U0=cbind(U0,pt(residuals(fit3,standardize=TRUE),nu))
  Z0=cbind(Z0,qfrechet(pt(residuals(fit3,standardize=TRUE),nu)))
}
plot.zoo(X0)
plot.zoo(U0)
plot.zoo(log(Z0))

R0=rowSums(Z0)
r0=quantile(R0,0.96)
W0=Z0/rowSums(Z0)
I0=(R0>=r0)
sum(I0)
pairs(W0[I0,],pch=20,xlim=c(0,1),ylim=c(0,1))

par(mfrow=c(5,2),
    pty='m',mgp=c(1,0.5,0),omi=c(0,0,0.5,0)+0.1,mar=c(2,3,1,0),
    cex=1.2)
for( d in 1:D){
  plot(Returns[,1],log(Z0[,d]),
       col=alpha(8+I0,0.2+I0),
       type='h',xlab='',ylab='',main=colnames(Returns[1+d]))
  abline(h=0,lwd=0.4)
}

par(mfrow=c(1,1),pty='m',omi=c(0,0,0,0),mar=c(3,3,0.1,0.1),las=1)
plot(Returns[,1],Returns[,2],type='h',xlab='',ylab='Negative log-returns')
plot(Returns[,1],log(Z0[,2]),type='h',xlab='',ylab='GARCH-Rescaled negative log-returns (log)',col=alpha(I0+8,0.2+I0))

#### Alternatively using Rank-based rescaling directly 
#rm(list = ls())

data1 <- read.csv("data1.csv")
data1=data1[,-c(2,12)]
missing=apply(data1,1,function(x) any(is.na(x)))
D=ncol(data1)-1
Returns=data.frame(as.Date(data1[-1,1]))
for(d in 1:D)
{
  Returns=cbind(Returns,as.numeric(-diff(log(data1[,1+d]))))
}

colnames(Returns)=colnames(data1)


Z0=apply(Returns[,-1],2,
         function(x) qfrechet(rank(x)/(length(x)+1)))

R0=rowSums(Z0)
r0=quantile(R0,0.96)
W0=Z0/rowSums(Z0)
I0=(R0>=r0)
sum(I0)
pairs(W0[I0,],pch=20,xlim=c(0,1),ylim=c(0,1))


par(mfrow=c(1,1),pty='m',omi=c(0,0,0,0),mar=c(3,3,0.1,0.1),las=1)
plot(Returns[,1],Returns[,2],type='h',xlab='',ylab='Negative log-returns')
plot(Returns[,1],log(Z0[,2]),type='h',xlab='',ylab='Rank-Rescaled negative log-returns (log)',col=alpha(I0+8,0.2+I0))


par(mfrow=c(5,2),
    pty='m',mgp=c(1,0.5,0),omi=c(0,0,0.5,0)+0.1,mar=c(2,3,1,0),
    cex=1.2)
for( d in 1:D){
  plot(Returns[,1],log(Z0[,d]),
       col=alpha(8+I0,0.2+I0),
       type='h',xlab='',ylab='',main=colnames(Returns[1+d]))
  abline(h=0,lwd=0.4)
}


