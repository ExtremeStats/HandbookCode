library("EnvStats")

set.seed(20240218)

n <- 10000
p <- c(0.8, 0.85, 0.9, 0.95, 0.99)

pdf("tdc_conceptual.pdf", width=8, height=4)
par(mar=c(4,4,2,9.5))

rho <- 0.3 
X1 <- rnorm(n)
X2 <- rho*X1 + sqrt(1-rho^2)*rnorm(n)

plot(p, 
     sapply(p, function(p) 
       sum( pmin(X1,X2) > qnorm(p) ) / sum(X2 > qnorm(p)) ),
    xlim=c(0.795,1), ylim=c(0,1), xaxs="i", yaxs="i", pch=19, type="b", 
    lty=3, lwd=2,
    xlab="quantile level", ylab="empirical tail dependence coefficients")

rho <- 0.7
X1 <- rnorm(n)
X2 <- rho*X1 + sqrt(1-rho^2)*rnorm(n)

lines(p, 
     sapply(p, function(p) 
       sum( pmin(X1,X2) > qnorm(p) ) / sum(X2 > qnorm(p)) ),
     pch=19, type="b", lty=1, lwd=2, col="red")

vario <- 1
X1 <- c(rep(1, times=n/2), exp(rnorm(n/2, sd=sqrt(vario)) - vario/2))
X2 <- c(exp(rnorm(n/2, sd=sqrt(vario)) - vario/2), rep(1, times=n/2))
Xsum <- X1 + X2
P  <- rpareto(n, location=1)
X1 <- P*X1/Xsum
X2 <- P*X2/Xsum

qval <- quantile(c(X1,X2), probs=p)

lines(p, 
      sapply(1:length(p), function(ind) 
        sum( pmin(X1,X2) > qval[ind] ) / sum(X2 > qval[ind]) ),
      pch=19, type="b", lty=2, lwd=2, col="blue")

par(xpd=T)

legend("right", inset=c(-0.3,0), 
       col=c("blue", "white", "red", "white", "black", "white"),
       lty=c(2,2,1,1,3,3), lwd=c(2,2,2,2,2,2),
       legend=c("Pareto", "process",  "Gaussian",  "(corr. 0.7)", 
                "Gaussian", "(corr. 0.3)"))
dev.off()