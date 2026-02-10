
#' Get residual paths (or vector) from the observations and the fitted model
#'
#' @param fitted.mod fitted model returned by the 'fit.Markov' function
#' @param filt.data filtered data returen by the 'filterData' function
#'
#' @return residuals from the filt.data with respect to the fitted.mod
#' @export
#'
#' @examples
get.z <- function(fitted.mod,filt.data){

  # Get the fitted alphas and betas using exceedances of 0.95 quantile on Laplace margins
  alphas <- fitted.mod$alphas
  beta <- fitted.mod$beta
  lengthsList <- sapply(filt.data, function(i) length(i))
  L <- length(alphas)
  get.ind <- which(lengthsList >= L+1)    # get clusters in which there is a at least L subsequent values after exceedance i.e. vectors must be of length L+1
  z.list  <- lapply(get.ind, function(i) (filt.data[[i]][2:(L+1)] - alphas[1:L]*filt.data[[i]][1]) / filt.data[[i]][1]^beta)
  z.list
}

#' Forward sampling (Algorithm 1)
#'
#' @param fitted.mod fitted model returned by the 'fit.Markov' function
#' @param filt.data filtered data returen by the 'filterData' function
#' @param u.thresh threshold above which to sample from the fitted.mod
#' @param sim.t_0 '>' to sample above u.thresh, '=' to sample at u.thresh
#' @param n.sim number of forward simulations
#'
#' @return forward simulations from the fitted model
#' @export
#'
#' @examples
forwardSim <- function(fitted.mod,filt.data,u.thresh,sim.t_0=">",n.sim){

  # Get a list of residuals from the CEVT model fit
  z.list <- get.z(fitted.mod,filt.data)

  alphas <- fitted.mod$alphas
  beta <- fitted.mod$beta

  sim.t_0 <- ifelse(sim.t_0==">",1,0)
  single.sim <- function(sim.t_0){     # function for generating a single cluster
    if(sim.t_0==1){
      x.init <- u.thresh + rexp(1)
    }else{
      x.init <- u.thresh
    }
      # initial exceedance
    n <- length(z.list)
    which.z <- sample(1:n, size=1)
    forward.vals <- alphas*x.init + (x.init^beta)*z.list[[which.z]]
    c(x.init, forward.vals) # simulated cluster
  }
  lapply(1:n.sim, function(i) single.sim(sim.t_0))
}
