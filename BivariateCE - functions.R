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

####################################################################################
### Function for negative log-likelihood of bivariate conditional extremes model ###
####################################################################################

CE.nll <- function(par, data){  
  alpha <- par[1]; beta <- par[2]; mu <- par[3]; sigma <- par[4]
  
  if(alpha < -1 | alpha > 1 | beta > 1 | sigma < 0){return(1e15)}
  
  y1 <- data[,1]; y2 <- data[,2]
  
  nll <- -sum(dnorm(y2, mean=alpha*y1 + mu*y1^beta, 
                      sd=sigma*y1^beta, log=T))
  
  return(nll)
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


