
#' Create blocks of observations following an exceedance of a threshold u
#'
#' @param df dataframe with columns 'dt' ('yyyymmdd' format), 'time.brk' indicating breaks in the series, 'X' data in standard Laplace margins
#' @param L number of lags after an exceedance to be considered in fitting procedure
#' @param u threshold to determine exceedances
#'
#' @return list of blocks of L observations
#' @export
#'
#' @examples
filterData <- function(df, L=20, u=-log(2*(1-0.95))){
  N <- nrow(df[df$X > u,]) # get number of exceedances
  out.lst <- vector("list", N)
  M <- 0  # counting variable

  if(inherits(unique(df$time.brk),"logical")){
    X   <- df$X
    exc <- which(X > u) # threshold exceedance times in year i
    for(j in 1:length(exc)){
      if(length(exc)!=0){
        M <- M + 1
        out.lst[[M]] <-  c(X[exc[j]],
                           X[(min(exc[j]+1, length(X))):(min(exc[j]+L, length(X)))])
      }
    }
  }else{
    t.b.range <- range(df$time.brk)
    for(i in t.b.range[1]:t.b.range[2]){
      time.brk.i.df <- df[df$time.brk==i,]#data  %>% filter(year==i)
      time.brk.i.X   <- time.brk.i.df$X
      exc <- which(time.brk.i.X > u) # threshold exceedance times in year i
      for(j in 1:length(exc)){
        if(length(exc)!=0){
          M <- M + 1
          out.lst[[M]] <-  c(time.brk.i.X[exc[j]],
                            time.brk.i.X[(min(exc[j]+1, nrow(time.brk.i.df))):(min(exc[j]+L, nrow(time.brk.i.df)))])
        }
      }
    }
  }
  out.lst
}


#' Fit time series autoregressive CE model of order 1
#'
#' @param pars initial value for the parameters theta and beta
#' @param data list returned by the filterData function
#' @param L number of lags after an exceedance to be considered in fitting procedure
#' @param lengthsList length of list 'data'
#' @param norming 'classic' or 'alt'
#'
#' @import stats
#' @return fitted model of order 1 with list of model estimates
#' @noRd
#'
#' @examples
Markov1 <- function(pars,
                    data,
                    L = 20,
                    lengthsList,
                    norming=c("classic", "alt")){
  norming <- match.arg(norming)
  alpha <- (2/pi)*atan(pars[1])
  beta <- 1 / (1 + exp(-pars[2]))
  alphas <- alpha^(1:L) # alpha parameters for lags 1 to L
  # now write function to get profile likelihood estimates of (mu_i, sigma_i), i=1...L
  prof.Mu.Sig <- function(no.Lags){
    list.indices <- which(lengthsList >= no.Lags+1) # + 1 is because initial exceedance is included in list
    lagged.data <- sapply(list.indices, function(i) data[[i]][no.Lags+1])
    data.exc <- sapply(list.indices, function(i) data[[i]][1])
    if(norming=="classic"){
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (data.exc^beta)
    } else {
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (1 + (alphas[no.Lags]*data.exc)^beta)
    }
    mu.z <- mean(z); sig.z <- sd(z)
    c(mu.z, sig.z)
  }
  # now sapply prof.Mu.Sig function to get all estimates of mu_i & sigma_i
  mus.sigs.z <- lapply(1:L, function(i) prof.Mu.Sig(i))
  mus.prof   <- sapply(1:L, function(i) mus.sigs.z[[i]][1])
  sigs.prof  <- sapply(1:L, function(i) mus.sigs.z[[i]][2])
  # now get log-lik contributions for each lag from 1 to L
  ll_cont <- rep(0, L) # to store likelihood contributions from different lags
  for(i in 1:L){
    list.indices <- which(lengthsList >= i+1)
    y <- sapply(list.indices, function(j) data[[j]][1])
    y.lag <- sapply(list.indices, function(j) data[[j]][i+1])
    if(norming=="classic"){
      mu.y.lag <- alphas[i]*y + (y^beta)*mus.prof[i]
      sd.y.lag <- (y^beta)*sigs.prof[i]
    }else{
      mu.y.lag <- alphas[i]*y + (1 + (alphas[i]*y)^beta)*mus.prof[i]
      sd.y.lag <- (1 + (alphas[i]*y)^beta)*sigs.prof[i]
    }
    ll_cont[i] <- sum(dnorm(x=y.lag, mean=mu.y.lag, sd=sd.y.lag, log=TRUE))
  }
  -sum(ll_cont)
}


#' Fit time series autoregressive CE model of order 2
#'
#' @param pars initial value for the parameters theta and beta
#' @param data list returned by the filterData function
#' @param L number of lags after an exceedance to be considered in fitting procedure
#' @param lengthsList length of list 'data'
#' @param norming 'classic' or 'alt'
#'
#' @import stats
#' @return fitted model of order 2 with list of model estimates
#' @noRd
#'
#' @examples
Markov2 <- function(pars,
                    data,
                    L = 20,
                    lengthsList,
                    norming=c("classic", "alt")){
  norming <- match.arg(norming)
  # pi1 and pi2 correspond to partial autocorrelations in autocorrelation function parametrization
  pi1 <- (2/pi)*atan(pars[1])
  pi2 <- (2/pi)*atan(pars[2])
  beta <- 1 / (1 + exp(-pars[3]))
  # theta1 and theta2 correspond to AR coefficients
  theta1 <- pi1*(1 - pi2)
  theta2 <- pi2
  alphas <- rep(0, L)
  alphas[1] <- theta1 / (1 - theta2)
  alphas[2] <- theta1*alphas[1] + theta2*1  # alpha0 = 1
  for(i in 3:L){
    alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2]
  }
  # now write function to get profile likelihood estimates of (mu_i, sigma_i), i=1...L
  prof.Mu.Sig <- function(no.Lags){
    list.indices <- which(lengthsList >= no.Lags+1) # + 1 is because initial exceedance is included in list
    lagged.data <- sapply(list.indices, function(i) data[[i]][no.Lags+1])
    data.exc <- sapply(list.indices, function(i) data[[i]][1])
    if(norming=="classic"){
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (data.exc^beta)
    } else {
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (1 + (alphas[no.Lags]*data.exc)^beta)
    }
    mu.z <- mean(z); sig.z <- sd(z)
    c(mu.z, sig.z)
  }
  # now sapply prof.Mu.Sig function to get all estimates of mu_i & sigma_i
  mus.sigs.z <- lapply(1:L, function(i) prof.Mu.Sig(i))
  mus.prof   <- sapply(1:L, function(i) mus.sigs.z[[i]][1])
  sigs.prof  <- sapply(1:L, function(i) mus.sigs.z[[i]][2])
  # now get log-lik contributions for each lag from 1 to L
  ll_cont <- rep(0, L) # to store likelihood contributions from different lags
  for(i in 1:L){
    list.indices <- which(lengthsList >= i+1)
    y <- sapply(list.indices, function(j) data[[j]][1])
    y.lag <- sapply(list.indices, function(j) data[[j]][i+1])
    if(norming=="classic"){
      mu.y.lag <- alphas[i]*y + (y^beta)*mus.prof[i]
      sd.y.lag <- (y^beta)*sigs.prof[i]
    }else{
      mu.y.lag <- alphas[i]*y + (1 + (alphas[i]*y)^beta)*mus.prof[i]
      sd.y.lag <- (1 + (alphas[i]*y)^beta)*sigs.prof[i]
    }
    ll_cont[i] <- sum(dnorm(x=y.lag, mean=mu.y.lag, sd=sd.y.lag, log=TRUE))
  }
  -sum(ll_cont)
}

#' Fit time series autoregressive CE model of order 3
#'
#' @param pars initial value for the parameters theta and beta
#' @param data list returned by the filterData function
#' @param L number of lags after an exceedance to be considered in fitting procedure
#' @param lengthsList length of list 'data'
#' @param norming 'classic' or 'alt'
#'
#' @import stats
#' @return fitted model of order 3 with list of model estimates
#' @noRd
#'
#' @examples
Markov3 <- function(pars,
                    data,
                    L = 20,
                    lengthsList,
                    norming=c("classic", "alt")){
  norming <- match.arg(norming)
  # pi1...pi3 correspond to partial autocorrelations
  pi1 <- (2/pi)*atan(pars[1])
  pi2 <- (2/pi)*atan(pars[2])
  pi3 <- (2/pi)*atan(pars[3])
  beta <-  1 / (1 + exp(-pars[4]))
  # theta1, theta2 & theta3 correspond to AR coefficients
  theta1 <- pi1 - pi1*pi2 - pi2*pi3
  theta2 <- pi2 - pi1*pi3 + pi1*pi2*pi3
  theta3 <- pi3
  alphas <- rep(0, L) # to store sequence of lagged alpha values
  alphas[1] <- (theta1 + theta2*theta3) / (1 - theta2 - theta1*theta3 - theta3^2)
  alphas[2] <- theta2 + (theta1 + theta3)*alphas[1]
  alphas[3] <- theta1*alphas[2] + theta2*alphas[1] + theta3*1  # alpha0 = 1
  for(i in 4:L){
    alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2] + theta3*alphas[i-3]
  }
  # now write function to get profile likelihood estimates of (mu_i, sigma_i), i=1...L
  prof.Mu.Sig <- function(no.Lags){
    list.indices <- which(lengthsList >= no.Lags+1) # + 1 is because initial exceedance is included in list
    lagged.data <- sapply(list.indices, function(i) data[[i]][no.Lags+1])
    data.exc <- sapply(list.indices, function(i) data[[i]][1])
    if(norming=="classic"){
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (data.exc^beta)
    } else {
      z <- (lagged.data - alphas[no.Lags]*data.exc) / (1 + (alphas[no.Lags]*data.exc)^beta)
    }
    mu.z <- mean(z); sig.z <- sd(z)
    c(mu.z, sig.z)
  }
  # now sapply prof.Mu.Sig function to get all estimates of mu_i & sigma_i
  mus.sigs.z <- lapply(1:L, function(i) prof.Mu.Sig(i))
  mus.prof   <- sapply(1:L, function(i) mus.sigs.z[[i]][1])
  sigs.prof  <- sapply(1:L, function(i) mus.sigs.z[[i]][2])
  # now get log-lik contributions for each lag from 1 to L
  ll_cont <- rep(0, L) # to store likelihood contributions from different lags
  for(i in 1:L){
    list.indices <- which(lengthsList >= i+1)
    y <- sapply(list.indices, function(j) data[[j]][1])
    y.lag <- sapply(list.indices, function(j) data[[j]][i+1])
    if(norming=="classic"){
      mu.y.lag <- alphas[i]*y + (y^beta)*mus.prof[i]
      sd.y.lag <- (y^beta)*sigs.prof[i]
    }else{
      mu.y.lag <- alphas[i]*y + (1 + (alphas[i]*y)^beta)*mus.prof[i]
      sd.y.lag <- (1 + (alphas[i]*y)^beta)*sigs.prof[i]
    }
    ll_cont[i] <- sum(dnorm(x=y.lag, mean=mu.y.lag, sd=sd.y.lag, log=TRUE))
  }
  -sum(ll_cont)
}


#' Fit time series autoregressive CE model of order 1, 2, or 3
#'
#' @param pars initial value for the parameters theta and beta
#' @param data list returned by the filterData function
#' @param L number of lags after an exceedance to be considered in fitting procedure
#' @param norming 'classic' or 'alt'
#' @param orderMarkov order of the markov process (1,2, or 3)
#' @param ...
#'
#' @return fitted model of order 1, 2, or 3 with list of model estimates
#' @export
#'
#' @examples
fit.Markov <- function(data,L = 20,pars,norming=c("classic", "alt"),orderMarkov,...){

  lengthsList <- sapply(data, function(i) length(i))
  # orderMarkov <- match.arg(orderMarkov)
  if(orderMarkov==1){
    fit <- optim(par=pars,fn=Markov1, data=data, L=L, lengthsList=lengthsList, norming=norming)
    raw.pars <- fit$par
    alpha <- (2/pi)*atan(raw.pars[1])
    beta <- 1 / (1 + exp(-raw.pars[2]))
    alphas <- alpha^(1:L)
    out.list <- list(alpha=alpha, beta=beta, alphas=alphas)
  } else if(orderMarkov==2){
    fit <- optim(par=pars,fn=Markov2, data=data, L=L, lengthsList=lengthsList, norming=norming)
    raw.pars <- fit$par
    pi1 <- (2/pi)*atan(raw.pars[1]) # pi1 & pi2 are partial autocorrelations
    pi2 <- (2/pi)*atan(raw.pars[2])
    beta <- 1 / (1 + exp(-raw.pars[3]))
    # theta1 and theta2 correspond to AR coefficients
    theta1 <- pi1*(1 - pi2)
    theta2 <- pi2
    alphas <- rep(0, L)
    alphas[1] <- theta1 / (1 - theta2)
    alphas[2] <- theta1*alphas[1] + theta2*1  # alpha0 = 1
    for(i in 3:L){
      alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2]
    }
    out.list <- list(pi1=pi1, pi2=pi2, theta1=theta1,theta2=theta2, beta=beta, alphas=alphas)
  }else if(orderMarkov==3){
    fit <- optim(par=pars,fn=Markov3, data=data, L=L, lengthsList=lengthsList, norming=norming)
    raw.pars <- fit$par
    # pi1...pi3 correspond to partial autocorrelations
    pi1 <- (2/pi)*atan(raw.pars[1])
    pi2 <- (2/pi)*atan(raw.pars[2])
    pi3 <- (2/pi)*atan(raw.pars[3])
    beta <-  1 / (1 + exp(-raw.pars[4]))
    # theta1, theta2 & theta3 correspond to AR coefficients
    theta1 <- pi1 - pi1*pi2 - pi2*pi3
    theta2 <- pi2 - pi1*pi3 + pi1*pi2*pi3
    theta3 <- pi3
    alphas <- rep(0, L)
    alphas[1] <- (theta1 + theta2*theta3) / (1 - theta2 - theta1*theta3 - theta3^2)
    alphas[2] <- theta2 + (theta1 + theta3)*alphas[1]
    alphas[3] <- theta1*alphas[2] + theta2*alphas[1] + theta3*1  # alpha0 = 1
    for(i in 4:L){
      alphas[i] <- theta1*alphas[i-1] + theta2*alphas[i-2] + theta3*alphas[i-3]
    }
    out.list <- list(pi1=pi1, pi2=pi2,pi3=pi3,
                     theta1=theta1,theta2=theta2,theta3=theta3,
                     beta=beta, alphas=alphas)
  }
  out.list
}

