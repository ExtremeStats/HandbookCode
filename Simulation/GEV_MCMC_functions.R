# Function to set initial value for the GEV parameters

initialize <- function(Y){
  loc   <- median(Y)+log(sqrt(6)*sd(Y)/pi)*log(log(2))
  scale <- sqrt(6)*sd(Y)/pi
  shape <- 0.001
return(c(loc,log(scale),shape))}


# Evaluate the log posterior PDF

log_post <- function(theta,Y,prior_mean,prior_sd){
   sum(dgev(Y,theta[1],exp(theta[2]),theta[3],log=TRUE))+
   sum(dnorm(theta,prior_mean,prior_sd,log=TRUE))
}

###########################################################################
# The main MCMC function to fit the model
#
#   Y[i] \sim GEV(theta[1],exp(theta[2]),theta[3])
#
#   Y        := Vector of observations
#   init     := Initial values of theta
#   prior_mn := Prior means of theta
#   prior_sd := Prior sds of theta
#   MH       := Initial values of the Metropolis candidate SDs
#   adapt    := Should MH be adapted during burn in?
#   iters    := Number of MCMC iterations
#   burn     := Length of burn-in period
#
###########################################################################

GEV_MCMC <- function(Y,init=NULL,
                     prior_mn=rep(0,3),prior_sd=rep(1,3),
                     MH=rep(0.3,3),adapt=TRUE,
                     iters=10000,burn=2000){
    library(evd)
    tick      <- proc.time()[3]

   # Set initial values 

    if(is.null(init)){init<-initialize(Y)}
    theta <- init

   # Storage    

    keepers           <- matrix(NA,iters,3)
    colnames(keepers) <- c("Location","Log scale","Shape")
    keepers[1,]       <- theta

   # Metropolis tuning

    att <- acc <- rep(0,3)
    curlp <- log_post(theta,Y,prior_mn,prior_sd) # current log posterior

   # MCMC sampling

    for(iter in 2:iters){

     # Parameter updates

      for(j in 1:3){
         att[j] <- att[j] + 1
         can    <- theta
         can[j] <- rnorm(1,theta[j],MH[j])
         canlp  <- log_post(can,Y,prior_mn,prior_sd)
         if(!is.na(canlp)){if(log(runif(1))<canlp-curlp){
            acc[j] <- acc[j] + 1
            theta  <- can
            curlp  <- canlp
         }}
      }   
   
     # Tuning the Metropolis candidate distribution

      if(adapt & iter<burn){for(j in 1:3){if(att[j]>50){
        if(acc[j]/att[j]<0.2){MH[j] <- MH[j]*0.8}
        if(acc[j]/att[j]>0.4){MH[j] <- MH[j]*1.2}
        acc[j] <- att[j] <- 0
      }}}

     # Store the MCMC samples

      keepers[iter,] <- theta  

    }
    tock <- proc.time()[3]
    out  <- list(samples = keepers,MH_acc = acc/att,time=tock-tick)
return(out)}
