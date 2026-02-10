# Function to set initial value for the GEV parameters

initialize_ref <- function(Y){
  loc   <- median(Y)+log(sqrt(6)*sd(Y)/pi)*log(log(2))
  scale <- sqrt(6)*sd(Y)/pi
  shape <- 0.001
return(c(loc,scale,shape))}

# Evaluate the log posterior PDF

log_post_ref <- function(theta,Y){
   sum(dgev(Y,theta[1],theta[2],theta[3],log=TRUE))+
   ref_log_prior(theta)
}

ref_log_prior <- function(theta){
   library(pracma)
   # (2.6) in https://www3.stat.sinica.edu.tw/ss_newpaper/SS-2021-0258_na.pdf
   mu  <- theta[1]
   sig <- theta[2]
   xi  <- theta[3]
   l   <- NA
   if((sig>0) & (xi> -0.5)){suppressWarnings({
     xi1 <- xi+1
     gam <- -digamma(1)
     der <- gamma(xi1)*psi(0,xi1)
     p   <- xi1*xi1*gamma(2*xi+1)
     q   <- xi*xi1*der + xi1*xi1*gamma(xi1)
     q   <- q/xi
     l   <- gamma(xi+2)/xi - q + 1 - gam
     l   <- pi*pi/6 + (1-gam)^2 - l*l/(1 + p - 2*gamma(xi+2))
     l   <- sqrt(l)/(sig*abs(xi))
     l   <- log(l)})
   }
return(l)}


###########################################################################
# The main MCMC function to fit the model
#
#   Y[i] \sim GEV(theta[1],exp(theta[2]),theta[3])
#
#   Y        := Vector of observations
#   init     := Initial values of theta
#   MH       := Initial values of the Metropolis candidate SDs
#   adapt    := Should MH be adapted during burn in?
#   iters    := Number of MCMC iterations
#   burn     := Length of burn-in period
#
###########################################################################

GEV_ref_MCMC <- function(Y,init=NULL,
                         MH=rep(0.3,3),adapt=TRUE,
                         iters=10000,burn=2000){
    library(evd)
    tick      <- proc.time()[3]

   # Set initial values 

    if(is.null(init)){init<-initialize_ref(Y)}
    theta <- init

   # Storage    

    keepers           <- matrix(NA,iters,3)
    colnames(keepers) <- c("Location","Scale","Shape")
    keepers[1,]       <- theta

   # Metropolis tuning

    att   <- acc <- rep(0,3)
    curlp <- log_post_ref(theta,Y) # current log posterior

   # MCMC sampling

    for(iter in 2:iters){

     # Parameter updates

      for(j in c(1,3)){
         att[j] <- att[j] + 1
         can    <- theta
         can[j] <- rnorm(1,theta[j],MH[j])
         canlp  <- log_post_ref(can,Y)
         if(!is.na(canlp)){
           if(log(runif(1))<canlp-curlp){
            acc[j] <- acc[j] + 1
            theta  <- can
            curlp  <- canlp
           }
         }
      }   

      att[2] <- att[2] + 1
      can    <- theta
      can[2] <- rlnorm(1,log(theta[2])-0.5*MH[2]^2,MH[2])
      canlp  <- log_post_ref(can,Y)
      if(!is.na(canlp)){
         R <- canlp-curlp +
              dlnorm(theta[2], log(can[2])-0.5*MH[2]^2, MH[2], log = TRUE)-
              dlnorm(can[2], log(theta[2])-0.5*MH[2]^2, MH[2], log = TRUE)
        if(log(runif(1))<R){
            acc[2] <- acc[2] + 1
            theta  <- can
            curlp  <- canlp
        }         
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


if(FALSE){
 mu  <- 3
 sig <- 2
 xi  <- .1
 Y   <- rgev(100,mu,sig,xi)
 fit <- GEV_ref_MCMC(Y)

 par(mfrow=c(3,1))
 plot(fit$samples[,1],type="l");abline(mu,0,col=2,lwd=2)
 plot(fit$samples[,2],type="l");abline(sig,0,col=2,lwd=2)
 plot(fit$samples[,3],type="l");abline(xi,0,col=2,lwd=2)
}