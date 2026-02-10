code <- nimbleCode({
   # likelihood
   for(s in 1:ns){for(t in 1:nt){
      y[s,t] ~ dgev(theta[s,1]+x[t]*theta[s,2], exp(theta[s,3]), theta[s,4])
   }}

   # Random effects
   for(s in 1:ns){
      theta[s,1:4] ~ dmnorm(RE[huc[s],1:4],Q2[1:4,1:4])
   }
   for(h in 1:H){
      RE[h,1:4] ~ dmnorm(mu[1:4],Q1[1:4,1:4])
   }
   # priors
   mu[1:4]     ~ dmnorm(A[1:4],B[1:4,1:4])
   Q1[1:4,1:4] ~ dwish(C[1:4, 1:4], D)
   Q2[1:4,1:4] ~ dwish(C[1:4, 1:4], D)
   S1[1:4,1:4] <- inverse(Q1[1:4,1:4])
   S2[1:4,1:4] <- inverse(Q2[1:4,1:4])
})

# Priors
 pri_sd <- c(10,10,10,0.5)
 pri_pc <- 1/pri_sd^2

# Prepare for nimble
 inits <- matrix(0,ns,4)
 for(s in 1:ns){
   temp      <- initialize(y[s,])
   inits[s,] <- c(temp[1],0,1+temp[2],temp[3])
 }
 inits <- list(theta=inits,mu=rep(0,4),Q1=0.1*diag(4),Q2=0.1*diag(4))
 data  <- list(y=y)
 const <- list(ns=ns,nt=nt,x=x,
               huc=huc,H=H,
               A = rep(0,4),       B = diag(pri_pc),
               C = diag(pri_pc)/5, D = 5)
# Fit the model
 tic <- proc.time()[3]
 fit <- nimbleMCMC(code=code, data = data, inits = inits,
                   constants=const,
                   monitors = c("mu","S1","S2","theta"), thin = thin,
                   niter = iters, nburnin = burn, nchains = 1,
                   summary = TRUE, WAIC = TRUE)
 toc   <- proc.time()[3]
 

# Summarize the output
 output       <- list()
 output$cpu   <- toc-tic
 output$WAIC  <- fit$WAIC
 output$S1    <- fit$summary[1:16,]
 output$S2    <- fit$summary[1:16+16,]
 output$mu    <- fit$summary[1:4+32,]
 output$beta0 <- fit$summary[1:ns+0*ns+36,]
 output$beta1 <- fit$summary[1:ns+1*ns+36,]
 output$lsig  <- fit$summary[1:ns+2*ns+36,]
 output$xi    <- fit$summary[1:ns+3*ns+36,]

 output$theta1_samps   <- fit$samples[,1   + ns*(0:3) + 36]
 output$theta10_samps  <- fit$samples[,10  + ns*(0:3) + 36]
 output$theta100_samps <- fit$samples[,100 + ns*(0:3) + 36]

 par(mfrow=c(2,2))
 plot(output$theta1_samps[,1],type="l")
 plot(output$theta1_samps[,2],type="l")
 plot(output$theta1_samps[,3],type="l")
 plot(output$theta1_samps[,4],type="l")

