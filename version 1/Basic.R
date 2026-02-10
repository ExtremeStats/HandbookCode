Lambda <- function (q, loc = 0, scale = 1, shape = 0) 
{
  if (min(scale) <= 0) 
    stop("invalid scale")
  if (length(shape) != 1) 
    stop("invalid shape")
  q <- (q - loc)/scale
  if (shape == 0) 
    p <- exp(-q)
  else p <- pmax(1 + shape * q, 0)^(-1/shape)
  p
}

lambda <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE) 
{
  if (min(scale) <= 0) 
    stop("invalid scale")
  if (length(shape) != 1) 
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0) 
    d <- log(1/scale) - x
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -Inf
  }
  if (!log) 
    d <- exp(d)
  d
}

nlogL0 <- function(th, z, y)
{ m <- ncol(y)+1 
ts1 <- -m*Lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3]) +
  lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3], log=TRUE) +log(m)
ts2 <- -Lambda(y, loc=th[1], scale=exp(th[2]), shape=th[3]) +
  lambda(y, loc=th[1], scale=exp(th[2]), shape=th[3], log=TRUE) 
ts3 <- (m-1)*Lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3])
-2*(sum(ts1) + sum(ts2) + sum(ts3))
}

nlogL <- function(th, z, y)
{ m <- ncol(y)+1 
  ts1 <- -m*Lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3]) +
            lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3], log=TRUE) + log(m)
  ts2 <- - Lambda(y, loc=th[4], scale=exp(th[5]), shape=th[6]) +
           lambda(y, loc=th[4], scale=exp(th[5]), shape=th[6], log=TRUE) 
  ts3 <- (m-1)*Lambda(z, loc=th[4], scale=exp(th[5]), shape=th[6])
  -2*(sum(ts1) + sum(ts2) + sum(ts3))
}

max.test <- function(y)
{
  # each row of y contains maxima
  
  fit2 <- fgev(c(y))  # original fit for comparison

  m <- ncol(y)
  n <- nrow(y)
  s <- y
  for (i in 1:nrow(y)) s[i,] <- sort(y[i,])
  z <- s[,m]
  y <- s[,-m]

#  fit1 <- fgev(z)
#  fit2 <- fgev(y)
#  print(c(fit0$deviance, fit1$deviance, fit2$deviance))
  init <- fit2$estimate
  init[2] <- log(init[2])
  fit0 <- optim(init, nlogL0, method="BFGS", z=z, y=y)
    
  init <- jitter(c(init, init))
  fit1 <- optim(init, nlogL, method="BFGS", z=z, y=y)
# fit1 <- optim(fit1$par, nlogL, method="CG", z=z, y=y)
#  L1 <- fit1$value
  list(L2 = fit2$deviance, L0=fit0$value, L1=fit1$value, W=fit0$value-fit1$value, 
       th2=fit2$estimate, th0=fit0$par, th1=fit1$par, 
       conv1=fit1$convergence, conv0=fit0$convergence)
}

# graphs to illustrate basic points

n <- 100; N <- 365*n; shape <- 2; scale <- 6
x <- matrix(rweibull(N,shape=shape,scale=scale),nrow=n)
y <- apply(x, MARGIN=1, FUN=max)

library(evd)
fit <- fgev(y)
u.p <- 0.9
u <- quantile(c(x),u.p)
fit1 <- fpot(c(x),threshold=u) 

amax <- 23
h <- hist(x,prob=T,nclass=200, plot=FALSE) 
a <- seq(from=0, to=amax, length=500)
da <- dweibull(a, shape=shape, scale=scale)

pdf(file="weibull.pdf",height=5,width=8)

par(mfrow=c(1,2),pty="s")

# density panel

plot(h$mids,h$density,type="s",xlab="Windspeed (m/s)",ylab="Density",log="y",ylim=c(10^-4,0.5),xlim=c(0,amax), 
     panel.first={ abline(v=u,col="grey",lty=2); lines(a,da,col="grey")},lwd=0.7)

# add EV approximations

lines(a,dgev(a,loc=fit$estimate[1],scale=fit$estimate[2],shape=fit$estimate[3])/365,lty=3,lwd=1.5,col="blue")
lines(a[a>u],fit1$pat*dgpd((a-u)[a>u],scale=fit1$estimate[1],shape=fit1$estimate[2]),lty=2,lwd=1.5,col="red")

# quantile panel

q <- qweibull(c(1:N)/(N+1),shape=shape,scale=scale)
z <- c(10,20,50,100,200,500,1000,5000)
plot(q,sort(c(x)),ylab="Empirical",xlab="Theoretical",pch=".",
     panel.first={abline(v=u,col="grey",lty=2); abline(0,1,col="grey"); 
       abline(v=qweibull(1-1/(365*z),shape=shape,scale=scale),col="grey",lwd=0.5)},xlim=c(0,amax),ylim=c(0,amax))
text(x=qweibull(1-1/(365*z),shape=shape,scale=scale),y=1*c(0:7),labels=z,cex=0.6,offset=0.5)

# add EV approximations

w <- seq(from=u, to=amax, length=500)
qw.gpd <- u+qgpd((pweibull(w,shape=shape,scale=scale)-u.p)/(1-u.p), scale=fit1$estimate[1], shape=fit1$estimate[2])
lines(w, qw.gpd, lty=2,lwd=1.5, col="red")

Np <- 5000
p <- c(1:Np)/(Np+1)
w.gev <- qweibull(p^(1/365), scale=scale, shape=shape)
# z <- seq(from=5,to=fit$estimate[1]-fit$estimate[2]/fit$estimate[3],length=100)
# gev.z <- pgev(z, loc=fit$estimate[1], scale=fit$estimate[2], shape=fit$estimate[3])
qw.gev <- qgev(p, loc=fit$estimate[1], scale=fit$estimate[2], shape=fit$estimate[3])
# qw.gev <- qgev((pweibull(w,shape=shape,scale=scale)-u.p)/(1-u.p), loc=fit$estimate[1], scale=fit$estimate[2], shape=fit$estimate[3])
lines(w.gev, qw.gev,lty=3,lwd=1.5,col="blue")

dev.off()

abisko.max <- matrix(NA, 102, 12)
year <- c(1913:2014)
mon <- c(1:12)
for (i in 1:102) for (j in 1:12) 
  { k <- (year(abisko$date)-1912==i & month(abisko$date)==j)
  abisko.max[i,j] <- max(abisko$precip[k]) }

# boxplots show many small values and seasonality, highest values in summer

boxplot(abisko$precip ~ month(abisko$date))

mon.max <- c(abisko.max)
mon.n <- length(mon.max)
year.max <- apply(abisko.max,1,max)
year.n <- length(year.max)
twoyear.max <- apply(matrix(year.max[-c(1,2)], ncol=2, byrow=T),1,max)
twoyear.n <- length(twoyear.max)
fiveyear.max <- apply(matrix(year.max[-c(1,2)], ncol=5, byrow=T),1,max)
fiveyear.n <- length(fiveyear.max)
tenyear.max <- apply(matrix(year.max[-c(1,2)], ncol=10, byrow=T),1,max)
tenyear.n <- length(tenyear.max)

library(evd)
pdf(file="abisko2.pdf",height=5,width=8)
par(pty="s", mfrow=c(1,2))
qqplot(qgumbel(c(1:mon.n)/(mon.n+1)), mon.max, pch="+",cex=0.5, 
       ylab="Ordered maxima", xlab="Gumbel plotting positions")
points(qgumbel(c(1:year.n)/(year.n+1)), sort(year.max), pch=16, cex=0.5)
points(qgumbel(c(1:twoyear.n)/(twoyear.n+1)), sort(twoyear.max), pch=15, cex=0.5)
points(qgumbel(c(1:fiveyear.n)/(fiveyear.n+1)), sort(fiveyear.max), pch=17, cex=0.5)
points(qgumbel(c(1:tenyear.n)/(tenyear.n+1)), sort(tenyear.max), pch=18, cex=0.5)

# analysis of annual maxima

library(mev)
library(evd)
# fit <- fgev(year.max)
# profile(fit)
out <- gev.pll(psi=seq(10,150,length=141),param='quant',dat=year.max, p=1/200)
# plot(out)

dev.off()

library(extremogram)
abisko.ts <- ts(abisko)

junk <- seq(from=min(abisko$date),to=max(abisko$date),by="day")
junk.y <- rep(0,length(junk))
junk.y[junk %in% abisko$date] <- abisko$precip


quant <- mean((junk.y>10))
n.u <- sum((junk.y>10))
q <- quant + 2*sqrt(quant*(1-quant)/n.u)

pdf(file="abisko4.pdf",height=4.5,width=8)
par(mfrow=c(1,2),pty="s")
exiplot(junk.y,c(0,30),ylab=expression(widehat(theta)[u]))
extremogram1(junk.y,quant=1-quant,maxlag=20,type=1)
abline(h=q,col="red")
abline(h=quant, col="grey")
dev.off()

# N <- 10
# T <- 10^seq(from=1, to=4, length=N) # return periods for extrapolation
# mid <- top <- bot <- rep(NA,N)
# for (i in 1:N) { 
#   out <- gev.pll(psi=seq(10,600,length=201),param='quant',dat=year.max, p=1/T[i])
#   mid[i] <- out$mle[2]
#   r <- out$psi[2*(out$maxpll-out$pll)<=qchisq(0.95,1)]
#   top[i] <- max(r)
#   bot[i] <- min(r)
# }
# plot(T,mid,ylim=range(mid,top,bot),type="l",log="x")

# POT analysis

library(mev)
library(evd)
pdf(file="/Users/davison/Desktop/abisko3.pdf",height=8,width=10)
par(mfrow=c(2,2))
NC.diag(xdat=c(abisko$precip),u=2*c(0:20),xlim=c(0,40)); abline(h=c(0.01,0.05),col="grey")
tcplot(data=c(abisko$precip), tlim=c(1,40), model="pp")
dev.off()

# create data with the added zero days; doesn't matter when, because of stationarity

y <- c(abisko$precip,rep(0,round(102*365.25)-length(abisko$precip)))

fpot(y,threshold=10)

library(ismev)
pp.fit(y, threshold=10, npy=365.25, method="BFGS", munit=20.4, siginit=5.85, shinit=0.05)

fpot(y, threshold=10, model="pp") #,std.err=F)

# Venice analysis

load("/Users/davison/Dropbox/Data/Venice/VeniceData1887To2023.R")

pdf(file="/Users/davison/Desktop/venice1.pdf",height=6,width=10)
# par(mfrow=c(2,2))
par(mar=c(3.1,3.1,1.1,1.1),mgp=c(1.5,0.5,0),mfrow=c(1,1)) 

plot(venice$year,venice$y[,1],ylim=c(50,200),xlab="",ylab="Sea level (cm)",pch=16,cex=0.6,col=colors()[17+7*venice$cens[,1]])
for (i in 2:10) points(venice$year,venice$y[,i],pch=16,cex=0.6,col=colors()[17+7*venice$cens[,i]])
lines(x+1900,fitted(fit.gam)[,1],col="slategrey")
dev.off()

pdf(file="/Users/davison/Desktop/venice2.pdf",height=8,width=10)
# par(mfrow=c(2,2))
par(mar=c(3.1,3.1,1.1,1.1),mgp=c(1.5,0.5,0),mfrow=c(2,2),pty="s") 

library(boot)
gum.sim <- boot(res,function(d) sort(rgumbel(length(d))), R=10000, sim="parametric")
gum.env <- envelope( gum.sim )
yl <- range(gum.sim$t,res)
qqplot(gum.pp, res,ylim=range(gum.sim$t,res),pch=16,cex=0.71,ylab="Ordered residuals",xlab="Gumbel plotting positions",panel.first={ abline(0,1,col="grey") })
lines(gum.pp, gum.env$point[1, ], lty = 4)
lines(gum.pp, gum.env$point[2, ], lty = 4)
lines(gum.pp, gum.env$overall[1, ], lty = 1)
lines(gum.pp, gum.env$overall[2, ], lty = 1)

fitr3 <- rlarg.fit(xdat=venice$y[venice$year<2020,1:3],ydat=X[venice$year<2020,],method="BFGS",mul=c(1))
junk <- make.exp.res(fitr3, venice$y[venice$year<2020,],3) 
dev.off()

load("/Users/davison/Dropbox/Data/Venice/VeniceData1887To2023.R")
summary(venice)
library(mgcv)
y <- venice$y[venice$year<2020,1]
x <- venice$year[venice$year<2020]-1900

fit.gam <- mgcv::gam(list(y~x,~1,~1),
                     family = gevlss(link = list("identity", "identity", "identity")))
fit.gam1 <- mgcv::gam(list(y~s(x),~1,~1),
                      family = gevlss(link = list("identity", "identity", "identity")))
lines(x+1900,fitted(fit.gam1)[,1],col="red")
lines(x+1900,fitted(fit.gam)[,1],col="blue")
