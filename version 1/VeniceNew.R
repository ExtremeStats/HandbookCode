# Venice analysis

# useful functions first 

qqgumbel <- function (y, env=T, gum.env = NA, yl=NULL, ...) 
{ # qqplot of y against Gumbel order statistics, with simulation envelope if requested
  i <- !is.na(y) & is.finite(y) & !is.nan(y)
  y <- y[i]
  n <- length(y)
  x <- qgumbel(c(1:n)/(n + 1))
  if (is.null(yl)) yl <- range(y)
  if (env)
  {
    library(boot)
    if (is.na(gum.env)) { 
      gum.sim <- boot(y, function(d) sort(rgumbel(length(d))), R=10000, sim="parametric", mle=NA)
      gum.env <- envelope( gum.sim ) }
    if (is.null(yl)) yl <- range(gum.env$overall,y)
  }
  qqplot(x, y, xlab = "Gumbel plotting positions", ylim = yl,, pch=16, cex=0.7, 
         ylab = "Ordered spacings", panel.first=abline(0,1,col="grey"), ...)
  if (env)
  {
    lines(x, gum.env$point[1, ], lty = 4)
    lines(x, gum.env$point[2, ], lty = 4)
    lines(x, gum.env$overall[1, ], lty = 1)
    lines(x, gum.env$overall[2, ], lty = 1)
  }
  invisible()
}

make.exp.res <- function(fit, data, r, env=TRUE, yl=NULL)
{ # make residuals using r-largest fit and differences of exponentials on transformed scale
  y <- data[,1:r]
  eta <- fit$vals[,1]
  tau <- fit$vals[,2]
  xi <- fit$vals[,3]
  res <- (y-eta)/tau 
  pp.res <- res.exp <- (1 + xi*res)^(-1/xi)
  
  for (i in 2:r) res.exp[,i] <- pp.res[,i] - pp.res[,i-1]
  for (i in 1:r)
  {
    res <- -log(res.exp[,i])
    qqgumbel(res, yl=yl)
  }
  invisible(res.exp)
}

# load data up to 2023, but only use data before MOSE entered use (i.e., before 2020)
# because no need to worry about the censoring for this analysis

load("/Users/davison/Dropbox/Data/Venice/VeniceData1887To2023.R")
summary(venice)
library(mgcv)
y <- venice$y[venice$year<2020,1]
x <- venice$year[venice$year<2020]  # this is for ease of plotting; better to regress on x-1900 for interpretation
fit.gam <- mgcv::gam(list(y~x,~1,~1),  # simple linear least squares fit to maxima
                     family = gevlss(link = list("identity", "identity", "identity")))

# First Venice figure

pdf(file="/Users/davison/Desktop/venice1.pdf",height=6,width=10)
par(mar=c(3.1,3.1,1.1,1.1),mgp=c(1.5,0.5,0),mfrow=c(1,1)) 

plot(venice$year,venice$y[,1],ylim=c(50,200),xlab="",ylab="Sea level (cm)",pch=16,cex=0.6,
     col=colors()[17+7*venice$cens[,1]], 
     panel.first=abline(coef(fit.gam)[1:2],col="red",lwd=2))
for (i in 1:10) points(venice$year,venice$y[,i],pch=16,cex=0.6,col=colors()[17+7*venice$cens[,i]])
dev.off()

# fit to annual maxima only

library(ismev)
library(evd)

fit1 <- fgev(x=y,nsloc=x-1900) # straight-line regression
muhat <- fit1$estimate[1] + fit1$estimate[2]*(x-1900)
sigmahat <- fit1$estimate[3]
xihat <- fit1$estimate[4]
res <- (y-muhat)/sigmahat
res <- log( 1 + xihat*res )/xihat
n <- length(res)
gum.pp <- qgumbel(c(1:n)/(n+1))

# fit to top three order statistics each year

X <- cbind(venice$year-1900,(venice$year>1981),(venice$year-1900)*(venice$year>1981))
Xdat <- venice$y[venice$year<2020,1:3]
Ydat <- X[venice$year<2020,]

fitr3 <- rlarg.fit(xdat=Xdat, ydat=Ydat, method="BFGS", mul=c(1))
res.3 <- (1+fitr3$vals[,3]*(Xdat-fitr3$vals[,1])/fitr3$vals[,2])^(-1/fitr3$vals[,3])

# Second Venice figure

pdf(file="/Users/davison/Desktop/venice2.pdf",height=8,width=10)
par(mar=c(3.1,3.1,1.1,1.1),mgp=c(1.5,0.5,0),mfrow=c(2,2),pty="s") 

library(boot)
gum.sim <- boot(res,function(d) sort(rgumbel(length(d))), R=10000, sim="parametric")
gum.env <- envelope( gum.sim )
yl <- c(-2.5,12) #range(gum.env$overall,res)
qqplot(gum.pp, res,ylim=yl,pch=16,cex=0.71,ylab="Ordered residuals",
       xlab="Gumbel plotting positions",panel.first={ abline(0,1,col="grey") })
lines(gum.pp, gum.env$point[1, ], lty = 4)
lines(gum.pp, gum.env$point[2, ], lty = 4)
lines(gum.pp, gum.env$overall[1, ], lty = 1)
lines(gum.pp, gum.env$overall[2, ], lty = 1)

junk <- make.exp.res(fit=fitr3, data=venice$y[venice$year<2020,], r=3, yl=yl) 

dev.off()
