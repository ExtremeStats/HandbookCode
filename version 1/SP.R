library(evd)
library(extremogram)


data1 <- load("Data/SP500.RData")

sp <- ts(data1[,2],start=c(2002,28), frequency=252)
ret <- -100*diff(log(sp))

#### Home-made AR-GARCH filtering

garch <- function(x, m=30, y)
{
  n <- length(y)
  j <- (m+1):n
  sig <- rep(NA,n)
  sig[m] <- var(y[1:m])
  mu <- x[1]
  rho <- x[2]
  b0 <- x[3]
  b1 <- x[4]
  d <- x[5]
  for (i in (m+1):n) sig[i] <- b0+b1*(y[i-1]-mu)^2 + d*sig[i-1]
  -sum( dnorm(y[j], mu+rho*(y[j-1]-mu), sqrt(sig[j]), log=T) )
}

# first fit a model with normal errors to get initial values for t fit

init <- c(mean(ret),0.0,0.2,0.2,0.1)
fit.g30 <- nlm( garch, init, hessian=T, print.level=2, y=ret)

# t error fit

garch <- function(x, m=30, y)
{
  n <- length(y)
  j <- (m+1):n
  sig <- rep(NA,n)
  sig[m] <- var(y[1:m])
  mu <- x[1]
  rho <- x[2]
  b0 <- x[3]
  b1 <- x[4]
  d <- x[5]
  df <- x[6]
  for (i in (m+1):n) sig[i] <- (b0+b1*(y[i-1]-mu)^2 + d*sig[i-1])
  sig <- sig*(df-2)/df
  z <- (y[j]- mu-rho*(y[j-1]-mu))/sqrt(sig[j])
  -sum( dt(z, df, log=T) - log(sqrt(sig[j])) )
}


init <- c(fit.g30$estimate,20)
fit.g30.t <- nlm( garch, init, hessian=T, print.level=2, y=ret)
fit.g30.t <- nlm( garch, fit.g30.t$estimate, hessian=T, print.level=2, y=ret)

fit.g30.t$estimate
sqrt(diag(solve(fit.g30.t$hessian)))

# t fit

# > fit.g30.t <- nlm( garch, fit.g30.t$estimate, hessian=T, print.level=2, y=ret)
# iteration = 0
# Parameter:
#   [1] -0.07907895 -0.05415534  0.01423904  0.12213208  0.87342285  6.18729015
# Function Value
# [1] 7326.875
# Gradient:
#   [1] -3.596142e-03  1.825356e-03 -7.071321e-03 -4.938556e-04  4.997673e-03 -1.322946e-05
#
# Relative gradient close to zero.
# Current iterate is probably solution.
#
# > fit.g30.t$estimate
# [1] -0.07907895 -0.05415534  0.01423904  0.12213208  0.87342285  6.18729015
# > sqrt(diag(solve(fit.g30.t$hessian)))
# [1] 0.009234118 0.014045780 0.002976498 0.011249524 0.010592493 0.525872578


# make residuals from the t fit

garch.res <- function(fit, m=30, y)
{
  n <- length(y)
  j <- (m+1):n
  sig <- rep(NA,n)
  sig[m] <- var(y[1:m])
  x <- fit$estimate
  mu <- x[1]
  rho <- x[2]
  b0 <- x[3]
  b1 <- x[4]
  d <- x[5]
  for (i in (m+1):n) sig[i] <- b0+b1*(y[i-1]-mu)^2 + d*sig[i-1]
  (y[j]-( mu+rho*(y[j-1]-mu)))/sqrt(sig[j])
}

ret.gres <- ts(garch.res( fit=fit.g30.t, y=ret),start=c(2002,28),frequency=252)
par(mfrow=c(2,2))
acf(ret.gres); pacf(ret.gres); acf(ret.gres^2); pacf(ret.gres^2)
q <- 0.1

# ret.gres now has white noise residuals


pdf(file="/Users/davison/Desktop/sp1.pdf",height=6,width=10)
# par(mfrow=c(2,2))
par(mar=c(3.1,3.1,1.1,1.1),mgp=c(1.5,0.5,0),mfrow=c(2,2))

# plot(sp,ylab="S&P 500", panel.first = abline(v=c(2002:2025),col="grey",lwd=0.5))
plot(ret,ylab="Negative log returns (%)",type="n",xlab="")
abline(v=c(2002:2025),col="grey",lwd=0.5)
lines(ret)
plot(ret,ylab="Filtered negative log returns",type="n",xlab="")
abline(v=c(2002:2025),col="grey",lwd=0.5)
lines(ret.gres)

q <- 0.05

extremogram::extremogram1(ret,1-q,type=1, maxlag=250)
abline(h=q,col="grey")
abline(h=q + 2*sqrt(1/length(ret)),col="red")

extremogram::extremogram1(ret.gres,1-q,type=1, maxlag=250)
abline(h=q,col="grey")
abline(h=q + 2*sqrt(1/length(ret.gres)),col="red")

dev.off()

# > exi(ret,u=quantile(ret,0.95))
# [1] 0.8741007
# > exi(ret.gres,u=quantile(ret.gres,0.95))
# [1] 0.9492754


extremogram::extremogram1(abs(ret),1-q,type=1, maxlag=250)
abline(h=q,col="grey")
abline(h=q + 2*sqrt(1/length(ret)),col="red")

extremogram::extremogram1(abs(ret.gres),1-q,type=1, maxlag=250)
abline(h=q,col="grey")
abline(h=q + 2*sqrt(1/length(ret.gres)),col="red")

