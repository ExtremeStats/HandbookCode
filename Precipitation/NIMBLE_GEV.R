rm(list = ls())
library(methods)  # otherwise new() not being found - weird
library(nimble)
library(extRemes)
library(coda)
library(postpack)

dgev <- nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0), xi = double(0), 
                   log = integer(0, default = 0)) {
        std <- (x - mu) / sigma
        if(1 + xi*std <= 0)
            logProb <- -Inf
        if(1 + xi*std > 0)
            logProb <- -log(sigma) - (1+1/xi) * log(1+xi*std) - (1+xi*std)^(-1/xi)
        returnType(double(0))
        if(log)  return(logProb)
        if(!log) return(exp(logProb))
    })

rgev <- nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0), xi = double(0)) {
        returnType(double(0))
        if(n != 1) print("rgev only allows n = 1; using n = 1.")
        u <- runif(1)
        if(xi == 0)
            x = -log(-log(u))
        if(xi != 0)
            x = ((-log(u))^(-xi)-1)/xi
        return (x*sigma + mu)
    })

data(Fort)
FortMax <- aggregate(Prec ~ year, data = Fort, FUN = max)
gev_model <- fevd(FortMax$Prec, type = 'GEV')
# gev_model$results$par[2] = log(gev_model$results$par[2])

code <- nimbleCode({
    # likelihood
    for(i in 1:n) {
        y[i] ~ dgev(mu, sig, xi)
    }
    # priors
    mu ~ dnorm(0, sd=10)
    log(sig) ~ dnorm(0, sd=2)
    xi ~ dnorm(0, sd = .5)
})

inits <- as.list(gev_model$results$par)
names(inits) <- c('mu', 'log_sig', 'xi')

output <- nimbleMCMC(code = code, data = list(y = FortMax$Prec),
                     constants = list(n = length(FortMax$Prec)),
                     inits = inits, niter = 11000, nburnin = 1000,
                     nchains = 2, summary = T, setSeed = 2:3)
output$summary
output$samples$chain1[,1] = exp(output$samples$chain1[,1])
output$samples$chain2[,1] = exp(output$samples$chain2[,1])
all.chains = do.call(rbind,output$samples)
dim(all.chains)
# all.chains[,1] = exp(all.chains[,1])
apply(all.chains,2,mean)
apply(all.chains,2,sd)
apply(all.chains,2,quantile,c(0.025,0.975))
# par(mfrow=c(1,3))
# plot(1:10000,output$samples[,1],'l')
# plot(1:10000,output$samples[,2],'l')
# plot(1:10000,output$samples[,3],'l')
# par(mfrow=c(1,1))


coda_samples = post_convert(output$samples)
plot(coda_samples)
geweke.diag(coda_samples)
gelman.diag(coda_samples)
effectiveSize(coda_samples)

par(mfrow=c(1,3))
plot(1:10000,all.chains[1:10000,2],'l',col='blue',xlab = 'Iteration',ylab = expression(mu),cex.axis=1.5,cex.lab=2)
lines(1:10000,all.chains[10001:20000,2],col='red')
plot(1:10000,all.chains[1:10000,1],'l',col='blue',xlab = 'Iteration',ylab = expression(sigma),cex.axis=1.5,cex.lab=2)
lines(1:10000,all.chains[10001:20000,1],col='red')
plot(1:10000,all.chains[1:10000,3],'l',col='blue',xlab = 'Iteration',ylab = expression(xi),cex.axis=1.5,cex.lab=2)
lines(1:10000,all.chains[10001:20000,3],col='red')
par(mfrow=c(1,1))


#
n   <- 100              # Number of observations
S   <- 2000             # Number of MCMC samples
Y   <- FortMax$Prec
indx = sample(1:20000,2000)
samps <- all.chains[indx,]

# Compute PPD summaries
tau <- seq(0.05,0.95,0.05)
nt  <- length(tau)
Qp     <- matrix(0,S,nt)
for(s in 1:S){
    Yp     <- evd::rgev(n,samps[s,2],samps[s,1],samps[s,3])
    Qp[s,] <- quantile(Yp,tau)
}

# Compute p-values
Q    <- quantile(Y,tau)
pval <- rep(0,nt)
for(j in 1:nt){pval[j]<-mean(Qp[,j]>Q[j])}

# I know you could do this in 5 mins, but anyways I coded up the computation of the PPD and Bayes p-value plot.  I'm thinking a side-by-side plot with the PPD on the left and pvals on the right.


# Plot PPD
pdf(file='plots/PPD_GEV.pdf',width=6,height = 5)
par(mar = c(5,5,2,1))
boxplot(t(Qp)~tau,outline=FALSE,ylim=c(0.45,5),
        cex.axis=1.5,cex.lab=1.5,
        xlab="Quantile level",ylab="Posterior predictive distribution")
points(Q,pch=19,col='blue',cex=0.5)
legend("topleft",c("PPD","Observed"),cex=1.5,pch=c(15,19),col=c("gray","blue"),bty="n")
par(mar = c(5,5,5,5))
dev.off()

# Plot p-values
plot(tau,pval,xlab="Quantile level",ylab="Bayesian p-value",
     cex.lab=1.5,cex.axis=1.5,ylim=0:1)

rlperiods = 2:100
rlvalues = matrix(NA,2000,99)
for(i in 1:2000){
    rlvalues[i,] = rlevd(rlperiods, loc=samps[i,2], scale=samps[i,1], shape=samps[i,3],type = 'GEV')    
}
rl_upper = apply(rlvalues, 2, quantile,0.975)
rl_lower = apply(rlvalues, 2, quantile,0.025)
rl_mean = apply(rlvalues, 2, mean)
pdf(file='plots/RL_GEV.pdf',width=6,height = 5)
par(mar = c(5,5,2,1))
plot(x=rlperiods,y=rl_mean,type = 'l',ylim = c(1.15,8), 
     xlab = 'Years (n)', ylab = 'Return level',cex.axis=1.5,cex.lab=1.5 )
lines(x=rlperiods, rl_lower,lty=2,col='grey')
lines(x=rlperiods, rl_upper,lty=2,col='grey')
par(mar = c(5,5,5,5))
dev.off()
