rm(list = ls())
library(nimble)
library(extRemes)

dgpd <- nimbleFunction(
run = function(x = double(0), mu = double(0), sigma = double(0), xi = double(0), 
               log = integer(0, default = 0)) {
    std <- (x - mu)/sigma
    r  <- log(1/sigma)
    id <- xi == 0
    good <- (std > 0 & ((1 + xi * std) > 0)) | is.na(std)
    if(id & good) 
        logProb <- r - std
    if(!id & good) 
        logProb <- r - (1/xi + 1) * log(1 + xi * std)
    if(!good) 
        logProb <- -Inf
    returnType(double(0))
    if(log)  return(logProb)
    if(!log) return(exp(logProb))
    })

rgpd <- nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0), xi = double(0)) {
        returnType(double(0))
        if(n != 1) print("rgpd only allows n = 1; using n = 1.")
        u <- runif(1)
        if (xi == 0) {
            x = mu - sigma * log(u)
        } else {
            x =  mu + sigma * (u^(-xi) - 1)/xi
        }
        return (x)
    })

assign('dgpd', dgpd, .GlobalEnv)
data(Fort)

## See monthly precip
Fortmonth <- aggregate(Prec ~ month, data = Fort, mean)
Fortmonth

# Subset the wet months
Fortwet = Fort[Fort$month %in% 4:9,]
head(Fortwet)

# Check autocorrelation
qq = quantile(Fortwet$Prec,0.95)
tmp1 <- tmp2 <- 0
n = nrow(Fortwet)
prec = Fortwet$Prec
for(i in 2:n){
    if(prec[i]>qq & prec[i-1]>qq)
        tmp1 = tmp1 + 1
    if(prec[i]>qq & prec[i-1]<=qq)
        tmp2 = tmp2 + 1
}
tmp1/sum(prec>qq)
tmp2/sum(prec<=qq)

gpd_model = fevd(x=Prec,data=Fortwet,threshold = qq,type = 'GP')

FortGPD = Fortwet$Prec[Fortwet$Prec>qq]
# acf(FortGPD)

code <- nimbleCode({
    # likelihood
    for(i in 1:n) {
        y[i] ~ dgpd(mu, sig, xi)
    }
    
    # priors
    log(sig) ~ dnorm(0, sd = 2)
    xi ~ dnorm(0, sd = .5)
})

inits <- as.list(gpd_model$results$par)
names(inits) <- c('log_sig', 'xi')
constants = list(n = length(FortGPD), mu=qq)

model <- nimbleModel(code, data = list(y = FortGPD), inits = inits,
                     constants = constants)
output <- nimbleMCMC(code = code, data = list(y = FortGPD),
                     constants = constants,
                     inits = inits, niter = 11000, nburnin = 1000,
                     nchains = 2, summary = T, setSeed = 2:3)
output$summary
output$samples$chain1[,1] = exp(output$samples$chain1[,1])
output$samples$chain2[,1] = exp(output$samples$chain2[,1])
all.chains = do.call(rbind,output$samples)
dim(all.chains)
apply(all.chains,2,mean)
apply(all.chains,2,sd)
apply(all.chains,2,quantile,c(0.025,0.975))

coda_samples = post_convert(output$samples)
plot(coda_samples)
geweke.diag(coda_samples)
gelman.diag(coda_samples)
effectiveSize(coda_samples)


#
n   <- 904              # Number of observations
S   <- 2000             # Number of MCMC samples
Y   <- FortGPD
indx = sample(1:20000,2000)
samps <- all.chains[indx,]

# Compute PPD summaries
tau <- seq(0.05,0.95,0.05)
nt  <- length(tau)
Qp     <- matrix(0,S,nt)
for(s in 1:S){
    Yp     <- evd::rgpd(n,qq,samps[s,1],samps[s,2])
    Qp[s,] <- quantile(Yp,tau)
}

# Compute p-values
Q    <- quantile(Y,tau)
pval <- rep(0,nt)
for(j in 1:nt){pval[j]<-mean(Qp[,j]>Q[j])}

# I know you could do this in 5 mins, but anyways I coded up the computation of the PPD and Bayes p-value plot.  I'm thinking a side-by-side plot with the PPD on the left and pvals on the right.


# Plot PPD
pdf(file='plots/PPD_GP.pdf',width=6,height = 5)
par(mar = c(5,5,2,1))
boxplot(t(Qp)~tau,outline=FALSE,ylim=range(Qp),
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
    rlvalues[i,] = rlevd(rlperiods, scale=samps[i,1], shape=samps[i,2], type = 'GP', threshold = qq)    
}
rl_upper = apply(rlvalues, 2, quantile,0.975)
rl_lower = apply(rlvalues, 2, quantile,0.025)
rl_mean = apply(rlvalues, 2, mean)
pdf(file='plots/RL_GP.pdf',width=6,height = 5)
par(mar = c(5,5,2,1))
plot(x=rlperiods,y=rl_mean,type = 'l',ylim = c(1.15,8), 
     xlab = 'Years (n)', ylab = 'Return level',cex.axis=1.5,cex.lab=1.5 )
lines(x=rlperiods, rl_lower,lty=2,col='grey')
lines(x=rlperiods, rl_upper,lty=2,col='grey')
par(mar = c(5,5,5,5))
dev.off()
