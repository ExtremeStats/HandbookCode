install.packages("stats4")
install.packages("statmod")
install.packages("evd")
install.packages("cobs")

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function transforms data to unit-Frechet scale using rank transformation
library("evd")

trans2frechet <- function(Data, method, thd = 0.95)
{
  require("evd")
  len <- length(Data)
  trans.data <- double(len)
  trans.data <- 1/(rank(Data) / (len+1))
  return(1/log(trans.data))
}
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns of cobs object that contains the corrected estimator of the Pickands dependence 
# function based on the min-projection of
# Mhalla, L., Opitz, T. and Chavez-Demoulin, V. (2019) Exceedance-based nonlinear regression of tail dependence. Extremes, 22, 523–552.
### X is a bivariate vector
### quant is the quantile level at which the threshold for the min-projected r.v. is set
### n.simplex is the number of values at which the min-projection is to be performed

A.biv.corrected.thd <-function(X,quant,n.simplex=500){
  
  #Function to generate n observations from the simplex
  unit.simplex.custom <- function(p,n)
  {
    if(p<=0){stop("Impossible value for dimension p")}
    if(p==1){stop("Copulas not defined for univariate data")}
    if(p==2){
      n      <- n-2
      out    <- c(0,seq(0,0.2,length.out = (n/10)+1)[-1],seq(0.2,0.4,length.out = (n/5)+1)[-1],
                  seq(0.4,0.6,length.out = (2*n/5)+1)[-1],
                  seq(0.6,0.8,length.out = (n/5)+1)[-1],seq(0.8,1,length.out = (n/10)+1)[-1],1)
      output <- cbind(out,1 - out)
    }
    return(output)
  }
  
  MEV2Exp <- function(dataset,w){
    n.obs <- nrow(dataset)
    dim   <- ncol(dataset)
    if(length(w) != (dim-1)) stop("Dimensions of the data and the angular observation do not match")
    
    #matrix of uniform r.v.s
    U <- NULL
    U <- apply(dataset,c(1,2),function(x){1/x})
    
    #Projection of the uniform r.v.s to obtain the univariate exp r.v Y
    U <- t(t(U)/c(w,1-w)) # we consider dimension d = 2 only
    Y <- apply(U,1,min)
    as.vector(Y)
  }
  
  simplex.parallel.deficit.link <- function(dataa,w,quant){
    y         <- MEV2Exp(dataa,w)
    u         <- rep(as.numeric(quantile(y,1-quant)),length(y))
    n.deficit <- sum(y<u)
    n.excess  <- length(y)-n.deficit
    mle.A     <- n.deficit/((n.excess*unique(u))+sum(y[which(y<u)]))
    se.A      <- n.deficit/(mle.A^2)
    return(c(mle.A,sqrt(1/se.A)))
  }
  
  #correct the Pickands (estimated at different values in [0,1])
  library("cobs")
  Omega          <- unit.simplex.custom(2,n.simplex+2)[-c(1,n.simplex+2),]
  
  Pickands.w.est <- NULL
  for(i in 1:nrow(Omega)){
    Pickands.w.est[i] <- simplex.parallel.deficit.link(X,Omega[i,1],quant)[1]
  }
  
  Pickands.w.est <- c(1,Pickands.w.est,1)
  #---conditions on the lower bound to add in the cobs fit
  conditions.cobs <- matrix(NA,nrow=n.simplex+2,ncol=3)
  Omegaa <- rbind(c(0,1),Omega,c(1,0))
  for(i in which(Omegaa[,1] >= 0.5)){
    conditions.cobs[i,] <- c(1,Omegaa[i,1],Omegaa[i,1])
  }
  for(i in which(Omegaa[,1] <= 0.5)){
    conditions.cobs[i,] <- c(1,Omegaa[i,1],1-Omegaa[i,1])
  }
  #---------------------------------------------------------
  
  return(cobs(c(0,Omega[,1],1),Pickands.w.est,pointwise=rbind(c(0,0,1),c(0,1,1),conditions.cobs), lambda=-1, 
              constraint="convex",print.warn=FALSE,print.mesg = FALSE))
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the corrected estimate of the Pickands dependence function evaluated at a value w
# in the unit simplex
### w corresponds to the angular observation at which we evaluate the estimate of the Pickands function
### ret is the cobs objected containing the corrected estimator of the Pickands function

estimate.A <- function(w,ret){
  return(predict(ret,w)[2])
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the first derivative of the corrected estimate of the Pickands dependence function 
# evaluated at a value w in the unit simplex
### w corresponds to the angular observation at which we evaluate the derivative of A(w,1-w)
### ret is the cobs objected containing the corrected estimator of the Pickands function

estimate.A.prime <- function(w,ret){
  return(predict(ret,w,deriv=1)[2])
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the estimate of the extreme value copula C^{EV} as given by Eq.(4) in the paper
### u and v are the values in [0,1] at which C^{EV} is estimated
### ret is the cobs objected containing the corrected estimator of the Pickands function

estimate.C.EV <- function(u,v,ret){
  u1 <- log(u)
  u2 <- log(v)
  w  <- u1/(u1+u2)
  return(exp((u1+u2)*estimate.A(w,ret)))
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the derivative wrt u of the copula C^{EV} evaluated at (u,v)
### u and v are the values in [0,1] at which C^{EV} is estimated
### ret is the cobs objected containing the corrected estimator of the Pickands function

d.u.C.EV <- function(u,v,ret){
  u1         <- log(u)
  u2         <- log(v)
  w          <- u1/(u1+u2)
  C.ev       <- estimate.C.EV(u,v,ret)
  A.ev       <- estimate.A(w,ret)
  A.ev.prime <- estimate.A.prime(w,ret)
  ret        <- C.ev*((A.ev/u)+((1-w)*A.ev.prime/u))
  return(ret)
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the derivative wrt v of the copula C^{EV} evaluated at (u,v)
### u and v are the values in [0,1] at which C^{EV} is estimated
### ret is the cobs objected containing the corrected estimator of the Pickands function
d.v.C.EV <- function(u,v,ret){
  u1         <- log(u)
  u2         <- log(v)
  w          <- u1/(u1+u2)
  C.ev       <- estimate.C.EV(u,v,ret)
  A.ev       <- estimate.A(w,ret)
  A.ev.prime <- estimate.A.prime(w,ret)
  ret        <- C.ev*((A.ev/v)-(w*A.ev.prime/v))
  return(ret)
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the cobs objected containing the corrected estimator of the Pickands function
get.A <- function(X,method=c("gpd"),thd=0.95,n.simplex=500){
  ret <- A.biv.corrected.thd(X,thd,n.simplex)
  return(ret)
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function evaluates the conditional quantile Q_{Y^{ext}|X^{ext}=x}(tau) as given by Eq.(13) in the paper
### tau is the quantile level
### x is the value at which we condition
### ret is the cobs objected containing the corrected estimator of the Pickands function describing the tail dependence in Z
### Z is a bivariate random vector in the MDA of an EV copula and with unit-Frechet margins
### thd is the quantile level at which we fix the marginal threshold u_X

cond.quant.Y <- function(tau,x,Z,ret.A,method=c("gpd"),thd=0.95){
  library("stats4")
  
  if((method=="gpd")&(x<quantile(Z[,1],thd))) 
    stop("GPD method: The value on which to condition should be larger than the threshold")
  
  if(method=="gpd"){
    unif.y  <- NULL
    quant   <- quantile(Z[,1],thd)
    fit.gpd <- evd::fpot(Z[,1],quant,model=method, std.err = FALSE,method="BFGS")
    unif.x  <- evd::pgpd(x, quant, fit.gpd$estimate[1], fit.gpd$estimate[2], lower.tail = TRUE)
    
    fct.uniroot <- function(y){
      d.u.C.EV(unif.x,y,ret.A)-tau
    }
    unif.y <- uniroot(fct.uniroot,interval = c(0.001,1))$root
  }
  
  #need to transform unif.y on the Z[,2] scale
  ret <- NULL
  if(method=="gpd"){
    quant   <- quantile(Z[,2],thd)
    fit.gpd <- evd::fpot(Z[,2],quant,model=method, std.err = FALSE,method="BFGS")
    ret     <- evd::qgpd(unif.y, quant, fit.gpd$estimate[1], fit.gpd$estimate[2], lower.tail = TRUE)
  }
  return(ret)
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function evaluates the conditional quantile Q_{X^{ext}|Y^{ext}=y}(tau) as given by Eq.(12) in the paper
### tau is the quantile level
### x is the value at which we condition
### ret is the cobs objected containing the corrected estimator of the Pickands function describing the tail dependence in Z
### Z is a bivariate random vector in the MDA of an EV copula and with unit-Frechet margins
### thd is the quantile level at which we fix the marginal threshold u_Y

# The following function evaluates the conditional quantile F_{X|Y=y}^{-1}(tau)
cond.quant.X <- function(tau,y,Z,ret.A,method=c("gpd"),thd=0.95){
  library("stats4")
  
  if((method=="gpd")&(y<quantile(Z[,2],thd))) 
    stop("GPD method: The value on which to condition should be larger than the threshold")
  
  if(method=="gpd"){
    unif.x  <- NULL
    quant   <- quantile(Z[,2],thd)
    fit.gpd <- evd::fpot(Z[,2],quant,model=method, std.err = FALSE,method="BFGS")
    unif.y  <- evd::pgpd(y, quant, fit.gpd$estimate[1], fit.gpd$estimate[2], lower.tail = TRUE)
    
    fct.uniroot <- function(x){
      d.v.C.EV(x,unif.y,ret.A)-tau
    }
    
    unif.x <- uniroot(fct.uniroot,interval = c(0.001,1))$root
  }
  
  #need to transform unif.x on the Z[which(Z[,1]>quantile(Z[,1],thd.marginal)),1] scale
  ret <- NULL
  if(method=="gev"){
    fit.gev     <- fgev(Z[,1])
    estimates   <- fit.gev$estimate
    ret         <- evd::qgev(unif.x,loc=estimates[1],scale=estimates[2],shape=estimates[3])
  }
  if(method=="gpd"){
    quant   <- quantile(Z[,1],thd)
    fit.gpd <- evd::fpot(Z[,1],quant,model=method, std.err = FALSE,method="BFGS")
    ret     <- evd::qgpd(unif.x, quant, fit.gpd$estimate[1], fit.gpd$estimate[2], lower.tail = TRUE)
  }
  return(ret)
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the expected quantile score for a given quantile level
quantileScoring <- function(actual, pred, prob = 0.95) {
  mean((as.numeric(actual <= pred) - prob) * (pred - actual))
}

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# The following function returns the causal score of CausEV as given by Eq.(15) in the paper

# implementation details:
# if 'X -> Y'  output cd  = 1
# if 'Y -> X'  output cd = 0
# epsilon: confidence (or score)

### pair is the set of bivariate observations (X,Y)
### thd is the quantile level at which we fix the marginal threshold and the threshold for estimating the Pickands function
### n.simplex is the number of values at which the min-projection is to be performed
### n.sim.ext is the minimum number of joint exceedances required to compute the causal score
### n.tau is the number of quantiles used to compute the causal score

QCDD_extremes_upper_quad <- function(pair,method=c("gpd"),thd=0.95,n.simplex=500,n.sim.ext=100,n.tau=5){
  library("stats4")
  n <- length(pair[,1])
  X <- pair[,1]
  Y <- pair[,2]
  x.extreme <- X[which((X>quantile(X,thd+0.005))&(Y>quantile(Y,thd+0.005)))]
  y.extreme <- Y[which((X>quantile(X,thd+0.005))&(Y>quantile(Y,thd+0.005)))]
  
  if((length(x.extreme)<n.sim.ext)|(length(y.extreme)<n.sim.ext))
    stop("not enough joint extremes to do something meaningful!")
  
  #transform data on unit Frechet scale: needed to estimate the Pickands function and hence the EV copula
  X1          <- trans2frechet(X,method,thd)
  X2          <- trans2frechet(Y,method,thd)
  ret.A       <- get.A(cbind(X1,X2),method=method,thd=thd,n.simplex)
  
  
  library("stats4")
  library("statmod")
  tau <- gauss.quad.prob(n.tau,l=0,u=1)
  h   <- sapply(tau$nodes, function(uu){
    if(method=="gpd"){
      xp <- sapply(y.extreme,function(y){cond.quant.X(tau=uu,y=y,pair,ret.A,method=method,thd)})
      yp <- sapply(x.extreme,function(x){cond.quant.Y(tau=uu,x=x,pair,ret.A,method=method,thd)})
    }
    
    h1 <- quantileScoring(unlist(yp),y.extreme, uu)
    h2 <- quantileScoring(unlist(xp),x.extreme, uu)
    
    quant   <- quantile(X,thd)
    fit.gpd <- evd::fpot(X,quant,model=method, std.err = FALSE,method="BFGS")
    ret     <- evd::qgpd(uu,quant, fit.gpd$estimate[1], fit.gpd$estimate[2], lower.tail = TRUE)
    h11     <- quantileScoring(ret,x.extreme, uu)
    
    quant   <- quantile(Y,thd)
    fit.gpd <- evd::fpot(Y,quant,model=method, std.err = FALSE,method="BFGS")
    ret     <- evd::qgpd(uu, quant, fit.gpd$estimate[1], fit.gpd$estimate[2], lower.tail = TRUE)
    h22     <- quantileScoring(ret,y.extreme, uu)
    
    rel_sc_1 = (h22+h2)/(h11 + h22 + h1 + h2)
    rel_sc_2 = (h11+h1)/(h11 + h22 + h1 + h2)
    c(rel_sc_1,rel_sc_2)
  })
  
  r11 <- sum(tau$weights[!is.na(h[1,])]*h[1,!is.na(h[1,])])/sum(tau$weights[!is.na(h[1,])]) 
  r12 <- sum(tau$weights[!is.na(h[2,])]*h[2,!is.na(h[2,])])/sum(tau$weights[!is.na(h[2,])])   
  
  cd = ifelse(r11 > r12, 1, 0)
  return(list(cd = cd, epsilon = r11/(r11+r12)))
}