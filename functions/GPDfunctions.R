#' negative log-likelihood (nllh) of the GP distribution
#'
#' @param pars parameters to parameterise sigma and xi
#' @param x observed data
#' @param threshold value of threshold above which to return the nllh
#' @param xi.tol xi tolerance
#'
#' @return negative log-likelihood of the GP distribution with respect to parameters
#' @noRd
#'
#' @examples
negllGPD <- function(pars, x, threshold, xi.tol=0.01){
  x <- x[x > threshold] - threshold  # threshold excesses
  xi  <- 0.1*pars[1]
  sig <- pars[2] / (1 + xi)
  n <- length(x)
  if(sig < 0){
    ll <- -1e10
  } else {
    if(abs(xi) <= xi.tol){ # use Gumbel likelihood when xi is close to zero
      ll <- -n*log(sig) - (1/sig)*sum(x)
    } else {
      ll <- max(-1e10, -n*log(sig) - (1/xi + 1)*sum(log(pmax(1+xi*x/sig, 0)))) # need max as second term may eval to -Inf
    }
  }
  -ll
}

#' Fit the GP distribution to exceedances of a threshold
#'
#' @param pars parameters to parameterise sigma and xi
#' @param x observed data
#' @param threshold value of threshold above which to return the nllh
#' @param xi.tol xi tolerance
#' @param ...
#'
#' @return Maximum likelihood estimates of sigma and xi and rate of exceedance
#' @noRd
#'
#' @examples
my.fit.gpd <- function(pars, x, threshold, xi.tol=0.01,...){
  fit <- optim(par=pars, fn=negllGPD, x=x, threshold=threshold, xi.tol=xi.tol,...)
  fit.mle <- fit$par
  xi  <- 0.1*fit.mle[1]
  sig <- fit.mle[2] / (1 + xi)
  mles <- c(xi, sig)
  names(mles) <- c("shape", "scale")
  rate <- sum(x > threshold) / length(x) # proportion of data that exceed threshold
  list(mles=mles, rate=rate)
}


#' Transform a univariate vector of observations to have approximate standard Laplace distribution
#'
#' @param x observed data in original margins
#' @param threshold threshold above which to fit the GP distribution for the semi-parametric model
#' @param pars initial values for the optimiser of the GP parameters
#' @param ...
#'
#' @return list of estimated GP parameters and vector of observations in Laplace margins
#' @export
#'
#' @examples
ToLaplace <- function(x, threshold, pars=c(-1,2),...){
  # first step is to get PIT values for each data point using semi-parametric model;
  # store in vector called cdf
  ranks <- rank(x) # get ranks
  ecdf  <- ranks / (length(x) + 0.5) # empirical cdf
  exc   <- which(x > threshold)
  non.exc <- which(x <= threshold)
  fit1 <- my.fit.gpd(pars=pars, x=x, threshold=threshold,...)
  sig <- fit1[["mles"]]["scale"]
  xi  <- fit1[["mles"]]["shape"]
  exec.rate <- fit1[["rate"]]
  cdf <- NULL
  cdf[non.exc] <- ecdf[non.exc]
  cdf[exc] <- 1 - exec.rate*(pmax((1 + xi*(x[exc] - threshold)/sig), 0))^(-1/xi)
  lower.half <- which(cdf < 0.5)
  upper.half <- which(cdf >= 0.5)
  lap <- NULL # now transform to Laplace scale
  lap[lower.half] <- log(2*cdf[lower.half])
  lap[upper.half] <- -log(2*(1 - cdf[upper.half]))
  # return transformed x on Laplace scale and GPD parameter estimates & exceedance rate
  return(list(lap=lap,threshold=threshold, sig=sig, xi=xi,rate=exec.rate))
}

#' Transform observations X from the data dataframe to have approximate standard Laplace distribution
#'
#' @param data dataframe with columns 'dt' ('yyyymmdd' format), 'time.brk' indicating breaks in the series, 'X' data in original margins
#' @param thresh threshold above which to fit GP distribution for the semiparametric model
#'
#' @return dataframe with observations in Laplace margins and GP parameters for back transformation
#' @export
#'
#' @examples
ToLaplace_df <- function(df, thresh){
  X.lap <- ToLaplace(x=df$X, threshold=thresh) # returns a list
  df$X  <- X.lap[[1]]
  list(df=df, thresh=X.lap[["threshold"]], sig=X.lap[["sig"]], xi=X.lap[["xi"]], rate=X.lap[["rate"]])
}

#' Back transform observations from Laplace margins to the orginal scale
#'
#' @param x.lap samples in standard Laplace margins
#' @param y.star threshold value in original value for the GP transformation in function ToLaplace_df
#' @param rate probability of threshold exceedance
#' @param sig scale parameter of the fitted GP distribution
#' @param xi shape parameter of the fitted GP distribution
#' @param x.orig original observations
#'
#' @return samples converted back to original marginal distribution
#' @export
#'
#' @examples
BackTransform <- function(x.lap, y.star, rate, sig, xi, x.orig){
  single.transform <- function(x){ # first function to transform a single numeric
    if(x <= 0){
      p <- exp(x) / 2
    } else {
      p <- 1 - exp(-x) / 2
    }
    if(p <= (1 - rate)){
      q <- unname(quantile(x.orig, probs=p))
    } else {
      q <- y.star + (sig/xi)*(-1 + ((1 - p)/rate)^(-xi))
    }
    q
  } # end of subfunction
  # now apply subfunction to all elements of x.lap
  sapply(1:length(x.lap), function(i) single.transform(x=x.lap[i]))
}

