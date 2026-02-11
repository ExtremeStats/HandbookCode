rm(list=ls())

library(evd) 
library(evmix)
library(scoringRules)
library(forecast)
library(data.table)
library(murphydiagram)

## R-code to reproduce figures and tables in 
## "Evaluation of extreme forecasts and projections" 

## Run times are from running the code in RStudio on a regular laptop

#############################################################################
## Figure 1: Competing forecasts 
pdf(file="Benchmark.pdf", points=8, width=6, height=2)
par(mfrow=c(1,2), mex=0.75, mar=c(4.5,4,1,1)+0.1, mgp=c(3,1.5,0))

D <- 10000
delta <- 1
x <- seq(0,15,length.out=D)
y1 <- dexp(x, rate=delta)
y2 <- evd::dgpd(x, shape=0.25)
y3 <- dexp(x, rate=delta/2)

plot(x, y1, type="l", col="gray60", xlab="y", ylab="Density",
     main="", ylim=c(0,1), xlim=c(0,8))
polygon(c(x[1],x,x[D]),c(0,y1,0),col="gray90",border="gray90")
lines(x, y2, lty=1, col="grey30")
lines(x, y3, lty=2, col="grey30")

delta <- 0.5
y1 <- dexp(x, rate=delta)
y2 <- evd::dgpd(x, shape=0.25)
y3 <- dexp(x, rate=delta/2)
plot(x, y1, type="l", col="gray60", xlab="y", ylab="Density",
     main="", ylim=c(0,1), xlim=c(0,8))
polygon(c(x[1],x,x[D]),c(0,y1,0),col="gray90",border="gray90")
lines(x, y2, lty=1, col="grey30")
lines(x, y3, lty=2, col="grey30")

dev.off()
#############################################################################

#############################################################################
## Table 2: Average quantile scores 
gamma <- 0.25
nu1 <- 0.9
nu2 <- 2
n <- 100
q.vec <- c(0.90, 0.95, 0.99, 0.999, 0.9999, 0.99999)

set.seed(13)

m <- 1000
all.scores <- list()
all.exceed <- list()
ideal.all <- clim.all <- ext1.all <- ext2.all <- array(NA, dim=c(n*m, length(q.vec)))
Delta.vec <- rgamma(n*m, shape=1/gamma, rate=1/gamma)
Delta.all <- array(Delta.vec, dim=c(m,n))
y.vec <- rexp(n*m, rate=as.vector(Delta.all))
y.all <- array(y.vec, dim=c(m,n))

for(j in 1:4) {
  print(j)
  scores <- array(NA, dim=c(m,4))
  exceed <- array(NA, dim=c(m,4))
  q <- q.vec[j]
  q.ideal.all <- q.clim.all <- q.ext1.all <- q.ext2.all <- NULL
  for(i in 1:m) {
    Delta <- Delta.all[i,] ## rgamma(n, shape=1/gamma, rate=1/gamma)
    y <- y.all[i,] ## rexp(n, rate=Delta)
    
    q.ideal <- qexp(q, rate=Delta)
    q.clim <- qgpd(q, sigmau=1, xi=gamma)
    q.ext1 <- qexp(q, rate=Delta/nu1)
    q.ext2 <- qexp(q, rate=Delta/nu2)
    
    tmp.ideal <- 1*(y < q.ideal)
    tmp.clim <- 1*(y < q.clim)
    tmp.ext1 <- 1*(y < q.ext1)
    tmp.ext2 <- 1*(y < q.ext2)
    scores[i,1] <- mean((tmp.ideal - q) * (q.ideal - y))
    scores[i,2] <- mean((tmp.clim - q) * (q.clim - y))
    scores[i,3] <- mean((tmp.ext1 - q) * (q.ext1 - y))
    scores[i,4] <- mean((tmp.ext2 - q) * (q.ext2 - y))
    
    exceed[i,1] <- n - sum(tmp.ideal)
    exceed[i,2] <- n - sum(tmp.clim)
    exceed[i,3] <- n - sum(tmp.ext1)
    exceed[i,4] <- n - sum(tmp.ext2)
    
  }
  all.scores <- append(all.scores, list(scores))
  all.exceed <- append(all.exceed, list(exceed))
}

round(apply(all.scores[[1]], 2, mean)[c(1,2,4)],3)
round(apply(all.scores[[2]], 2, mean)[c(1,2,4)],3)
round(apply(all.scores[[3]], 2, mean)[c(1,2,4)],3)
round(apply(all.scores[[4]], 2, mean)[c(1,2,4)],3)
#############################################################################

#############################################################################
## Figure 2: Plot of quantile scores 
pdf(file="QSextremes.pdf", width=5, height=5, points=12)
cols <- c("gray70", "gray50", "black")
pchs <- c(16,15,1,1)
par(mar=c(4.5,4.5,2,2)+0.1, mex=0.75)
plot(all.scores[[4]][,4], all.scores[[4]][,2], col=cols[all.exceed[[4]][,2]+1],
     pch=pchs[all.exceed[[4]][,2]+1], cex=0.75,
     xlab="Quantile score for extremist forecast", ylab="Quantile score for marginal forecast")
abline(a=0, b=1)
dev.off()
#############################################################################

#############################################################################
## Figure 3: Score penalties 
pdf(file="CRPSvsLogS.pdf", points=10, width=6, height=4)
par(mfrow=c(2,2), mex=0.75, mar=c(4.5,4,1,1)+0.1, mgp=c(3,1.5,0))

## GEV distribution
x1 <- evd::rgev(1e6, loc=0, scale=1, shape=0.5)
x2 <- evd::rgev(1e6, loc=0, scale=1, shape=0.5)
e1 <- mean(abs(x1 - x2))
D <- 20000
x <- seq(-2,15,length.out=D)
evd::pgev(18, loc=0, scale=1, shape=0.5) 
e.xy.1 <- rep(NA,D)
## The for loop below takes about a minute to run 
for(i in 1:D) { 
  e.xy.1[i] <- mean(abs(x1-x[i]))
}
crps.1 <- e.xy.1 - 0.5 * e1
logs.1 <- -evd::dgev(x, loc=0, scale=1, shape=0.5, log=TRUE)
## Where do they intersect?
ind <- which.min(abs(logs.1-crps.1))
evd::pgev(x[ind], loc=0, scale=1, shape=0.5) 

y <- evd::dgev(x, loc=0, scale=1, shape=0.5)
plot(x, y, type="l", col="gray90", xlab="y", ylab="",
     axes=FALSE, main="", ylim=c(0,0.4))
axis(1)
polygon(c(x[1],x,x[D]),c(0,y,0),col="gray90",border="gray90")
box()
par(new = TRUE)
plot(x, crps.1, type="l", col="black", axes=FALSE, xlab="",
     ylab="Score", ylim=c(0,10), lty=3)
axis(2)
lines(x, logs.1, col="gray50", lty=1)
text(14.5, 9.6, "(a)", cex=1.2)

## Gumbel distribution
x1 <- evd::rgumbel(1e6, loc=0, scale=1)
x2 <- evd::rgumbel(1e6, loc=0, scale=1)
e1 <- mean(abs(x1 - x2))
D <- 20000
x <- seq(-4,13,length.out=D)
evd::pgumbel(13, loc=0, scale=1)
e.xy.1 <- rep(NA,D)
## The for loop below takes about a minute to run 
for(i in 1:D) {
  e.xy.1[i] <- mean(abs(x1-x[i]))
}
crps.1 <- e.xy.1 - 0.5 * e1
logs.1 <- -evd::dgumbel(x, loc=0, scale=1, log=TRUE)
## Where do they intersect?
ind <- which.min(abs(logs.1-crps.1))
evd::pgumbel(x[ind], loc=0, scale=1)

y <- evd::dgumbel(x, loc=0, scale=1)
plot(x, y, type="l", col="gray90", xlab="y", ylab="",
     axes=FALSE, main="", ylim=c(0,0.4))
axis(1)
polygon(c(x[1],x,x[D]),c(0,y,0),col="gray90",border="gray90")
box()
par(new = TRUE)
plot(x, crps.1, type="l", col="black", axes=FALSE, xlab="",
     ylab="Score", ylim=c(0,10), lty=3)
axis(2)
lines(x, logs.1, col="gray50", lty=1)
text(12.5, 9.6, "(b)", cex=1.2)

## GPD distribution
x1 <- evd::rgpd(1e6, loc=0, scale=1, shape=0.5)
x2 <- evd::rgpd(1e6, loc=0, scale=1, shape=0.5)
e1 <- mean(abs(x1 - x2))
D <- 20000
x <- seq(0,17,length.out=D)
evd::pgpd(17, loc=0, scale=1, shape=0.5)
e.xy.1 <- rep(NA,D)
## The for loop below takes about a minute to run 
for(i in 1:D) {
  e.xy.1[i] <- mean(abs(x1-x[i]))
}
crps.1 <- e.xy.1 - 0.5 * e1
logs.1 <- -evd::dgpd(x, loc=0, scale=1, shape=0.5, log=TRUE)
## Where do they intersect?
ind <- sort(abs(logs.1-crps.1), index=TRUE)$ix
evd::pgpd(x[ind[1:5]], loc=0, scale=1, shape=0.5) 

y <- evd::dgpd(x, loc=0, scale=1, shape=0.5)
plot(x, y, type="l", col="gray90", xlab="y", ylab="",
     axes=FALSE, main="", ylim=c(0,0.4))
axis(1)
polygon(c(x[1],x,x[D]),c(0,y,0),col="gray90",border="gray90")
box()
par(new = TRUE)
plot(x, crps.1, type="l", col="black", axes=FALSE, xlab="",
     ylab="Score", ylim=c(0,10), lty=3)
axis(2)
lines(x, logs.1, col="gray50", lty=1)
text(16.5, 9.6, "(c)", cex=1.2)

## Other GPD distribution
x1 <- evd::rgpd(1e6, loc=0, scale=1)
x2 <- evd::rgpd(1e6, loc=0, scale=1)
e1 <- mean(abs(x1 - x2))
D <- 20000
x <- seq(0,17,length.out=D)
evd::pgpd(13, loc=0, scale=1) 
e.xy.1 <- rep(NA,D)
## The for loop below takes about a minute to run 
for(i in 1:D) {
  e.xy.1[i] <- mean(abs(x1-x[i]))
}
crps.1 <- e.xy.1 - 0.5 * e1
logs.1 <- -evd::dgpd(x, loc=0, scale=1, log=TRUE)
## Where do they intersect?
ind <- which.min(abs(logs.1-crps.1))
evd::pgpd(x[ind], loc=0, scale=1) 

y <- evd::dgpd(x, loc=0, scale=1)
plot(x, y, type="l", col="gray90", xlab="y", ylab="",
     axes=FALSE, main="", ylim=c(0,0.4))
axis(1)
polygon(c(x[1],x,x[D]),c(0,y,0),col="gray90",border="gray90")
box()
par(new = TRUE)
plot(x, crps.1, type="l", col="black", axes=FALSE, xlab="",
     ylab="Score", ylim=c(0,10), lty=3)
axis(2)
lines(x, logs.1, col="gray50", lty=1)
text(16.5, 9.6, "(d)", cex=1.2)

dev.off()
#############################################################################

#############################################################################
## Figure 4: Threshold weighted and censored scores 
## Note: The code for this figure takes a while to run 
pdf(file="wScores.pdf", points=8, width=6, height=2)
par(mfrow=c(1,2), mex=0.75, mar=c(4.5,4,1,1)+0.1, mgp=c(3,1.5,0))

D <- 1000

## GEV distribution and weighted CRPS
x <- seq(-2,15,length.out=D)
sample_fc <- matrix(evd::rgev(5e8, loc=0, scale=1, shape=0.5), nrow=D)
u <- evd::qgev(0.9, loc=0, scale=1, shape=0.5)

## The next line of code takes about 3 minutes to run 
y1 <- twcrps_sample(y=x, dat=sample_fc, a=u)
## The next line of code takes about 10 minutes to run 
y2 <- owcrps_sample(y=x, dat=sample_fc, a=u)
## The next line of code takes about a minute to run 
y3 <- crps_sample(y=x, dat=sample_fc)
y4 <- 0.5 * y3 + 0.5 * y1

y <- evd::dgev(x, loc=0, scale=1, shape=0.5)
plot(x, y, type="l", col="gray90", xlab="y", ylab="",
     axes=FALSE, main="", ylim=c(0,0.4))
axis(1)
polygon(c(x[1],x,x[D]),c(0,y,0),col="gray90",border="gray90")
box()
par(new = TRUE)
plot(x, y1, type="l", col="gray50", axes=FALSE, xlab="",
     ylab="Score", ylim=c(0,10), lty=3)
axis(2)
lines(x, y2, col="black", lty=1)
lines(x, y3, col="black", lty=2)
lines(x, y4, col="gray50", lty=1)
text(14.5, 9.6, "(a)", cex=1.2)

## GEV distribution and weighted LogS
## The next line of code takes about 6 minutes to run 
z1 <- clogs_sample(y=x, dat=sample_fc, a=u, cens=FALSE)
## The next line of code takes about 7 minutes to run 
z2 <- clogs_sample(y=x, dat=sample_fc, a=u, cens=TRUE)
## The next line of code takes about 3 minutes to run 
z3 <- logs_sample(y=x, dat=sample_fc)

y <- evd::dgev(x, loc=0, scale=1, shape=0.5)
plot(x, y, type="l", col="gray90", xlab="y", ylab="",
     axes=FALSE, main="", ylim=c(0,0.4))
axis(1)
polygon(c(x[1],x,x[D]),c(0,y,0),col="gray90",border="gray90")
box()
par(new = TRUE)
plot(x, z1, type="l", col="gray50", axes=FALSE, xlab="",
     ylab="Score", ylim=c(0,10), lty=3)
axis(2)
lines(x, z2, col="black", lty=1)
lines(x, z3, col="black", lty=2)
text(14.5, 9.6, "(b)", cex=1.2)

dev.off()
#############################################################################

#############################################################################
## Table 3: Mean scores 
set.seed(11)

nu3 <- 0.8
nu4 <- 1.2

m <- 100 
n <- 1000
k <- 5e3

pit <- array(NA, dim=c(n,3,m))
crps <- array(NA, dim=c(n,3,m))
logS <- array(NA, dim=c(n,3,m))
cl <- array(NA,dim=c(n,3,m))
tw <- array(NA,dim=c(n,3,m))

for(j in 1:m) {
  print(j)
  
  Delta <- rgamma(n, shape=1/gamma, rate=1/gamma)
  y <- rexp(n, rate=Delta)
  
  pit[,1,j] <- pexp(y, rate=Delta)
  pit[,2,j] <- pexp(y, rate=Delta/nu3)
  pit[,3,j] <- pexp(y, rate=Delta/nu4)
  
  crps[,1,j] <- crps_exp(y, rate=Delta)
  crps[,2,j] <- crps_exp(y, rate=Delta/nu3)
  crps[,3,j] <- crps_exp(y, rate=Delta/nu4)
  
  logS[,1,j] <- -dexp(y, rate=Delta, log=TRUE)
  logS[,2,j] <- -dexp(y, rate=Delta/nu3, log=TRUE)
  logS[,3,j] <- -dexp(y, rate=Delta/nu4, log=TRUE)
  
  u <- 3.11
  
  cl[,1,j] <- -1 * (y > u) * log(dexp(y, rate=Delta, log=FALSE)/(1 - pexp(u, rate=Delta)))
  cl[,2,j] <- -1 * (y > u) * log(dexp(y, rate=Delta/nu3, log=FALSE)/(1 - pexp(u, rate=Delta/nu3)))
  cl[,3,j] <- -1 * (y > u) * log(dexp(y, rate=Delta/nu4, log=FALSE)/(1 - pexp(u, rate=Delta/nu4)))
  
  tw[,1,j] <- 1 * (y < u) * (1 - pexp(u, rate=Delta*2)) / (2*Delta) + 
    1 * (y > u) * ((y - u) - (2 * (pexp(y, rate=Delta) - pexp(u, rate=Delta)))/Delta 
                      + (1 - pexp(u, rate=Delta*2)) / (2*Delta) )
  tw[,2,j] <- 1 * (y < u) * (1 - pexp(u, rate=Delta*2/nu3)) / (2*Delta/nu3) + 
    1 * (y > u) * ((y - u) - (2 * (pexp(y, rate=Delta/nu3) - pexp(u, rate=Delta/nu3)) * nu3)/Delta 
                      + ((1 - pexp(u, rate=Delta*2/nu3)) * nu3 ) / (2*Delta) )
  tw[,3,j] <- 1 * (y < u) * (1 - pexp(u, rate=Delta*2/nu4)) / (2*Delta/nu4) + 
    1 * (y > u) * ((y - u) - (2 * (pexp(y, rate=Delta/nu4) - pexp(u, rate=Delta/nu4)) * nu4)/Delta 
                      + ((1 - pexp(u, rate=Delta*2/nu4)) * nu4 ) / (2*Delta) )
  
}

crps.m <- apply(crps, 2:3, mean)
logS.m <- apply(logS, 2:3, mean)
cl.m <- apply(cl, 2:3, mean)
tw.m <- apply(tw, 2:3, mean)

round(apply(crps.m, 1, mean), 3)
round(apply(logS.m, 1, mean), 3)
round(apply(tw.m, 1, mean), 3)
round(apply(cl.m, 1, mean), 3)

round(apply(crps.m, 1, sd), 3)
round(apply(logS.m, 1, sd), 3)
round(apply(tw.m, 1, sd), 3)
round(apply(cl.m, 1, sd), 3)
#############################################################################

#############################################################################
## Table 4: Testing predictive performance 
## This code is a continuation of the code for Table 3
perm.test <- function(s1, s2, n) {
  N <- length(s1)
  d <- s1 - s2
  s.perm <- rep(NA,n)
  for(i in 1:n) {
    si <- sample(c(-1,1), size=N, replace=TRUE)
    s.perm[i] <- mean(d * si)
  }
  s <- mean(d)
  return(ecdf(c(s, s.perm))(s))
}

dt <- data.table(score="CRPS", test="DM", m1=1, m2=2, p=0.05)
dt <- c("CRPS","DM",1,2,0.05)
N <- 1000 
## The for loop takes about 2 minutes to run 
for(i in 1:100) {
  p <- perm.test(crps[,1,i], crps[,2,i], n=5000)
  dt <- rbind(dt, c("CRPS", "perm", 1, 2, p))
  p <- perm.test(crps[,1,i], crps[,3,i], n=5000)
  dt <- rbind(dt, c("CRPS", "perm", 1, 3, p))
  p <- perm.test(crps[,2,i], crps[,3,i], n=5000)
  dt <- rbind(dt, c("CRPS", "perm", 2, 3, p))
  p <- perm.test(logS[,1,i], logS[,2,i], n=5000)
  dt <- rbind(dt, c("logS", "perm", 1, 2, p))
  p <- perm.test(logS[,1,i], logS[,3,i], n=5000)
  dt <- rbind(dt, c("logS", "perm", 1, 3, p))
  p <- perm.test(logS[,2,i], logS[,3,i], n=5000)
  dt <- rbind(dt, c("logS", "perm", 2, 3, p))
  p <- perm.test(cl[,1,i], cl[,2,i], n=5000)
  dt <- rbind(dt, c("CL", "perm", 1, 2, p))
  p <- perm.test(cl[,1,i], cl[,3,i], n=5000)
  dt <- rbind(dt, c("CL", "perm", 1, 3, p))
  p <- perm.test(cl[,2,i], cl[,3,i], n=5000)
  dt <- rbind(dt, c("CL", "perm", 2, 3, p))
  p <- perm.test(tw[,1,i], tw[,2,i], n=5000)
  dt <- rbind(dt, c("tw", "perm", 1, 2, p))
  p <- perm.test(tw[,1,i], tw[,3,i], n=5000)
  dt <- rbind(dt, c("tw", "perm", 1, 3, p))
  p <- perm.test(tw[,2,i], tw[,3,i], n=5000)
  dt <- rbind(dt, c("tw", "perm", 2, 3, p))
  p <- pnorm(sqrt(N) * mean(crps[,1,i]-crps[,2,i])/sd(crps[,1,i]-crps[,2,i])) 
  dt <- rbind(dt, c("CRPS", "DM", 1, 2, p))
  p <- pnorm(sqrt(N) * mean(crps[,1,i]-crps[,3,i])/sd(crps[,1,i]-crps[,3,i])) 
  dt <- rbind(dt, c("CRPS", "DM", 1, 3, p))
  p <- pnorm(sqrt(N) * mean(crps[,2,i]-crps[,3,i])/sd(crps[,2,i]-crps[,3,i])) 
  dt <- rbind(dt, c("CRPS", "DM", 2, 3, p))
  p <- pnorm(sqrt(N) * mean(logS[,1,i]-logS[,2,i])/sd(logS[,1,i]-logS[,2,i])) 
  dt <- rbind(dt, c("logS", "DM", 1, 2, p))
  p <- pnorm(sqrt(N) * mean(logS[,1,i]-logS[,3,i])/sd(logS[,1,i]-logS[,3,i])) 
  dt <- rbind(dt, c("logS", "DM", 1, 3, p))
  p <- pnorm(sqrt(N) * mean(logS[,2,i]-logS[,3,i])/sd(logS[,2,i]-logS[,3,i])) 
  dt <- rbind(dt, c("logS", "DM", 2, 3, p))
  p <- pnorm(sqrt(N) * mean(cl[,1,i]-cl[,2,i])/sd(cl[,1,i]-cl[,2,i]))
  dt <- rbind(dt, c("CL", "DM", 1, 2, p))
  p <- pnorm(sqrt(N) * mean(cl[,1,i]-cl[,3,i])/sd(cl[,1,i]-cl[,3,i]))
  dt <- rbind(dt, c("CL", "DM", 1, 3, p))
  p <- pnorm(sqrt(N) * mean(cl[,2,i]-cl[,3,i])/sd(cl[,2,i]-cl[,3,i]))
  dt <- rbind(dt, c("CL", "DM", 2, 3, p))
  p <- pnorm(sqrt(N) * mean(tw[,1,i]-tw[,2,i])/sd(tw[,1,i]-tw[,2,i]))
  dt <- rbind(dt, c("tw", "DM", 1, 2, p))
  p <- pnorm(sqrt(N) * mean(tw[,1,i]-tw[,3,i])/sd(tw[,1,i]-tw[,3,i]))
  dt <- rbind(dt, c("tw", "DM", 1, 3, p))
  p <- pnorm(sqrt(N) * mean(tw[,2,i]-tw[,3,i])/sd(tw[,2,i]-tw[,3,i]))
  dt <- rbind(dt, c("tw", "DM", 2, 3, p))
}
dt <- as.data.table(dt)
names(dt) <- c("score", "test", "m1", "m2", "p")
dt <- dt[-1,]
dt$p <- as.numeric(dt$p)

dt[,.(mean.p=mean(p), sign=(mean(p)<0.025 | mean(p)>0.975), signif=sum(p<0.025)+sum(p>0.975), sp=sum(p>0.975)), 
   by=.(m1,m2,score,test)]

## Testing for all cases combined. 
## Discussed in text below Table 4
perm.test(as.vector(logS[,1,]), as.vector(logS[,2,]), n=5000) 
perm.test(as.vector(logS[,1,]), as.vector(logS[,3,]), n=5000) 
perm.test(as.vector(logS[,2,]), as.vector(logS[,3,]), n=5000) 

s1 <- as.vector(logS[,2,])
s2 <- as.vector(logS[,3,])
N <- length(s1)
sqrt(N)*mean(s1-s2)/sd(s1-s2)

perm.test(as.vector(tw[,1,]), as.vector(tw[,2,]), n=5000) 
perm.test(as.vector(tw[,1,]), as.vector(tw[,3,]), n=5000) 
perm.test(as.vector(tw[,2,]), as.vector(tw[,3,]), n=5000) 

s1 <- as.vector(tw[,2,])
s2 <- as.vector(tw[,3,])
N <- length(s1)
sqrt(N)*mean(s1-s2)/sd(s1-s2)

## Investigating if model A or B rejected
m1 <- apply(tw[,1,],2,mean)
m2 <- apply(tw[,2,],2,mean)
m3 <- apply(tw[,3,],2,mean)
(m1-m3 < 0)
(m1-m2 < 0)
cbind(dt[score=="tw" & test=="DM" & m1==1 & m2==3]$p, (m1-m3 < 0))
#############################################################################

#############################################################################
## Figure 5: Calibration assessment 
## This code is a continuation of the code for Tables 3 and 4 
pdf(file="calibration.pdf", height=6, width=6, points=12)
par(mfrow=c(3,3), mex=0.75)
hist(pit[,1,1], nclass=10, ylim=c(0,140), xlab="PIT values", 
     main="Ideal forecast")
abline(h=100, lty=2, col="gray50")
hist(pit[,2,1], nclass=10, ylim=c(0,140), xlab="PIT values", 
     main="Light tails")
abline(h=100, lty=2, col="gray50")
hist(pit[,3,1], nclass=10, ylim=c(0,140), xlab="PIT values", 
     main="Heavy tails")
abline(h=100, lty=2, col="gray50")
par(mar=c(5,4,2,2)+0.1)
qqplot((1:n)/(n+1), pit[,1,1], col="gray70", pch=16, cex=0.6, 
       xlab="Expected", ylab="Observed") 
abline(a=0, b=1, col="gray30", lty=2)
qqplot((1:n)/(n+1), pit[,2,1], col="gray70", pch=16, cex=0.6, 
       xlab="Expected", ylab="Observed") 
abline(a=0, b=1, col="gray30", lty=2)
qqplot(punif((1:n)/(n+1)), pit[,3,1], col="gray70", pch=16, cex=0.6, 
       ylab="Observed", xlab="Expected") 
abline(a=0, b=1, col="gray30", lty=2)
qqplot(qgumbel(((1:n)-0.44)/(n-2*0.44+1)), qgumbel(pit[,1,1]), col="gray70", 
       pch=16, cex=0.6, 
       xlim=c(-2,8), ylim=c(-2,9), ylab="Observed (Gumbel-scale)",
       xlab="Expected") 
abline(a=0, b=1, col="gray30", lty=2)
qqplot(qgumbel(((1:n)-0.44)/(n-2*0.44+1)), qgumbel(pit[,2,1]), col="gray70", 
       pch=16, cex=0.6, 
       xlim=c(-2,8), ylim=c(-2,9), ylab="Observed (Gumbel-scale)",
       xlab="Expected") 
abline(a=0, b=1, col="gray30", lty=2)
qqplot(qgumbel(((1:n)-0.44)/(n-2*0.44+1)), qgumbel(pit[,3,1]), col="gray70", 
       pch=16, cex=0.6, 
       xlim=c(-2,8), ylim=c(-2,9), ylab="Observed (Gumbel-scale)",
       xlab="Expected") 
abline(a=0, b=1, col="gray30", lty=2)
dev.off()
#############################################################################

#############################################################################
## Figure 6: Murphy diagrams
## Note: This code takes a while to run 
K <- 1e5
set.seed <- 13
Delta <- rgamma(K, shape=1/gamma, rate=1/gamma)
y <- rexp(K, rate=Delta)
for(i in 1:4) {
  q <- q.vec[i]
  q.ideal <- qexp(q, rate=Delta)
  q.clim <- rep(evmix::qgpd(q, sigmau=1, xi=gamma), K)
  q.ext2 <- qexp(q, rate=Delta/nu2)
  if(i==1) {
    first.ext <- q.ext2 
    first.clim <- q.clim
  }
  if(i==4) {
    snd.ext <- q.ext2 
    snd.clim <- q.clim
  }
}

pdf(file="MurphyDiagram.pdf", width=10, height=5, points=10)
par(mfrow=c(1,2))
## The next line of code takes approximately 9 minutes to run 
murphydiagram(first.ext, first.clim, y, functional="quantile", alpha=0.9, 
              labels=c("Extremist", "Marginal"), colors=c("black", "gray50"))
## The next line of code takes approximately 9 minutes to run 
murphydiagram(snd.ext, snd.clim, y, functional="quantile", alpha=0.999, 
              labels=c("Extremist", "Marginal"), colors=c("black", "gray50"))
dev.off()
#############################################################################
