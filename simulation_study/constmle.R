constmle <- function(data, k, w, start=NULL, type="maxima", r0=NULL, q=NULL){
  ################################
  # START auxilary functions
  llik_pickands_min <- function(coef){
    # compute the difference of the coefficients:
    beta1 <- diff(coef); beta2 <- diff(beta1)
    
    # compute the Pickands and its derivatives:
    A <- c(bpb %*% coef)
    A1 <- c(k * (bpb1 %*% beta1))
    A2 <- c(k * (k-1) * (bpb2 %*% beta2))
    
    # GEV distribution
    lG <- -pseudo$z*A
    
    C <- pseudo$w*A1
    A1.s <- (A - C[,1])*(A + C[,2]) / pseudo$den[,2]
    A2.s <- A2*pseudo$w2[,1]/pseudo$r
    if(any((A1.s+A2.s)<0)){
      llik <- 1e10
    }else{
      llik <- lG + log(A1.s+A2.s) - log(pseudo$den[,1])
      llik <- -sum(llik)  
    }
    
    if(is.na(llik)) llik <- 1e10
    return(llik)
  }
  
  llik_h_beta <- function(coef){
    if (any(is.na(coef))) {
      llik <- Inf
      return(llik)
    }
    # define the eta coefficients
    #eta <- net(beta=coef,from='A')$eta
    eta <- 1/2 + k*diff(c(1, coef, 1))/2
    # compute the angular density:
    h <- (k-1)*c(bpb2 %*% diff(eta))
    # check if it is the same to dh in CODE_Simulation...
    # is it faster?
    if(any(h<0)){
      llik <- Inf
    }else{
      llik <- log(h)
      llik <- -mean(llik)      
    }
    
    if(is.na(llik)) llik <- Inf
    return(llik)
  }
  
  confun <- function(x) {
      return(constraint$r - constraint$R%*%c(1, x, 1))
  }
  
  
  # END auxilary functions
  ################################
  kp <- k+1
  km <- k-1
  # define the linear constraints
  constraint <- ExtremalDep:::constraints(k)
  # define the setup for the estimation
  pseudo <- list(z=rowSums(1/data), r=rowSums(data), w=data/rowSums(data), 
                 r2=(rowSums(data))^2, w2=(data/rowSums(data))^2, den=data^2)
  # define options
  #opts <- list(algorithm="NLOPT_LN_NELDERMEAD", xtol_rel=1e-04, maxeval=100)
  #opts <- list(algorithm="NLOPT_GN_ISRES", xtol_rel=1e-04, maxeval = 100000)
  opts <- list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-25, maxeval = 100000)
  # opts_local <- list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-4)
  # opts <- list(algorithm="NLOPT_LN_AUGLAG", maxeval = 100000, xtol_rel=1e-4,
  #              local_opts = opts_local, print_level = 0)
  #opts <- list(algorithm="NLOPT_LN_NEWUOA_BOUND", xtol_rel=1e-25, maxeval=100000)
  #opts <- list(algorithm="NLOPT_LN_SBPLX", xtol_rel=1e-08)
  # Compute starting values
  #if(is.null(start)){
  #  eta0 <- prior_eta_sampler(km)$eta
  #  beta0 <- net(eta=eta0, from='H')$beta}
  #else beta0 <- start
  if(is.null(start)) beta0 <- ExtremalDep:::rcoef(km)$beta
  if (abs(mean(start) - 1) < 10^{-4}) beta0 <- ExtremalDep:::rcoef(km)$beta
  else beta0 <- start
  # maximize the constraint log-likelihood function
  if(type=="maxima"){
    bpb <- bpb1 <- bpb2 <- NULL
    bpb <- ExtremalDep:::bp2d(pseudo$w, k)
    bpb1 <- ExtremalDep:::bp2d(pseudo$w, km)
    bpb2 <- ExtremalDep:::bp2d(pseudo$w, km - 1)
    fit <- nloptr::nloptr(beta0, eval_f=llik_pickands_min, eval_g_ineq=confun, 
                  lb=rep(0,kp), ub=rep(1,kp), opts=opts)
  }
  if(type=="rawdata"){
    bpb2 <- NULL 
    bpb2 <- ExtremalDep:::bp2d(pseudo$w[pseudo$r>r0, ], km-1)
    
    if (all(beta0 == 1)) beta0 <- as.matrix(ExtremalDep:::rcoef(km)$beta)
    
    fit <- list()
    fvalue <- c(Inf, numeric(14))
    fit[[1]] <- nloptr::nloptr(beta0[-c(1, kp)], eval_f=llik_h_beta, eval_g_ineq=confun,
                       lb=rep(0,kp - 2), ub=rep(1,kp - 2), opts=opts)
    for(i in 2:15){
      fit[[i]]<- nloptr::nloptr(fit[[i-1]]$solution[(kp-2):1], eval_f=llik_h_beta, eval_g_ineq=confun,
                        lb=rep(0,kp - 2), ub=rep(1,kp - 2), opts=opts)
      fvalue[i] <- fit[[i]]$objective}
    
    ind <- c(which(fvalue==min(fvalue)))[1]
    fit <- fit[[ind]]
  }
  if (type=="rawdata") beta <- c(1, fit$solution[(kp-2):1], 1)
  else beta <- fit$solution[kp:1]
  A <- ExtremalDep:::beed(data, cbind(w, 1-w), 2, 'md', 'emp', k, beta=beta, plot=FALSE)
  
  return(list(beta=beta, A=A$A, status=fit$status, start=beta0, iterations=fit$iterations, message=fit$message))
}