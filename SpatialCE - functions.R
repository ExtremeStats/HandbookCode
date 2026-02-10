library(VGAM)

###########################################################
### Function to apply rank transform to Laplace margins ###
###########################################################

laplace.margins <- function(data){
  
  ## Apply rank transform
  data.U <- rank(data)/(length(data)+1)
  
  ## Transform to standard Laplace
  data.L <- data
  data.L[data.U<0.5] <- log(2*data.U[data.U<0.5])
  data.L[data.U>=0.5] <- -log(2*(1-data.U[data.U>=0.5]))
  
  return(data.L)
  
}

############################################################################################
### Function for negative profile log-likelihood of bivariate conditional extremes model ###
############################################################################################

CE.nll.PL <- function(par, data){
  alpha <- par[1]; beta <- par[2]
  
  if(alpha < -1 | alpha > 1 | beta > 1){return(1e15)}
  
  y1 <- data[,1]; y2 <- data[,2]
  
  mu <- mean((y2-alpha*y1)/y1^beta)
  
  sigma <- sqrt(mean(((y2-alpha*y1-mu*y1^beta)/y1^beta)^2))
  
  nll <- -sum(dnorm(y2, mean=alpha*y1 + mu*y1^beta, 
                      sd=sigma*y1^beta, log=T))
  
  return(nll)
}

##################################################
### Function to carry out stationary bootstrap ###
##################################################
boot.stat <- function(data, mean.length){
  n <- nrow(data)
  block.start <- sample(n,n,rep=T)
  block.length <- 1 + rgeom(n, 1/mean.length)
  
  sample.index <- NULL
  iter <- 1
  while(length(sample.index) < n){
    sample.index <- c(sample.index, c(block.start[iter]:(block.start[iter]+block.length[iter]-1)))
    sample.index[sample.index > n] <- NA
    sample.index <- na.omit(sample.index)
    iter <- iter + 1
  }
  
  data.boot <- data[sample.index,]
  
  return(data.boot)
}

###############################################################################
### Function to calculate empirical estimate of chi(u) - on Laplace margins ###
###############################################################################
chi.u <- function(data1, data2, u){
  
  q <- qlaplace(u)

  #chi.u <- mean(data1>q & data2>q)/(1-u)
  chi.u <- mean(data2[data1>q]>q)
  
  return(chi.u)
}

###########################################################################################################
### Function to simulate observations (on Laplace scale), conditioning on a particular site being large ###
###########################################################################################################

### Note: adapted from the function "rCondSiteMod3" in the supplementary material of Wadsworth and Tawn for our specific modeling choices ###

rCondSite.new3 <- function(n, coord, u, lambda, b, phi, mu, nu, sig, kappa, delta, DM=dists, site){
  a<-alpha(c(DM[site,]),H=0,lambda=lambda,kappa=kappa)
  nsite<-dim(coord)[1]
  
  sigma<- sig^2 * exp(-(DM/phi)^nu)
  newMu<-rep(mu,nsite-1) - sigma[-site,site]%*%solve(sigma[site,site])%*%(mu-0)
  newSig<-sigma[-site,-site]-sigma[-site,site]%*%(solve(sigma[site,site])%*%sigma[site,-site])
  
  newZ<-rmvdlaplace(n=n,delta=delta,dim=dim(newSig)[1],mu=newMu,Sigma=newSig,sigmad=sqrt(diag(newSig)))
  if(dim(newSig)[1]==1){newZ<-t(newZ)} 
  
  newY1<-rexp(n)+u
  if(n>1){
    scaledZ<-(1+(diag(newY1)%*%matrix(rep(a[-site],n),nrow=n,byrow=T))^b) * newZ
    if(site==1){scaledZ<-cbind(rep(0,n),scaledZ)}
    else if(site==dim(coord)[1]){scaledZ<-cbind(scaledZ,rep(0,n))}
    else{scaledZ<-cbind(scaledZ[,1:(site-1)],rep(0,n),scaledZ[,site:(dim(coord)[1]-1)])}
    Data<-diag(newY1)%*%matrix(rep(a,n),nrow=n,byrow=T)+ scaledZ
  } else{
    scaledZ<-(1+((newY1)*matrix(rep(a[-site],n),nrow=n,byrow=T))^b)*newZ
    if(site==1){scaledZ<-c(rep(0,n),scaledZ)}
    else if(site==dim(coord)[1]){scaledZ<-c(scaledZ,rep(0,n))}
    else{scaledZ<-c(scaledZ[,1:(site-1)],rep(0,n),scaledZ[,site:(dim(coord)[1]-1)])}
    Data<-(newY1)*matrix(rep(a,n),nrow=n,byrow=T)+ scaledZ
  }
  return(Data)
}

