#' Log-Likelihood for Threshold Family Distributions
#'
#' This function calculates the log-likelihood of the mark-distribution process at the data level, i.e., the likelihood of density \code{A | theta_A, mu}. It supports both gamma and log-normal threshold families.
#'
#' @param cur_par A numeric vector representing the hyperparameters of the threshold family distribution. For gamma, this is the shape parameter \code{k}. For log-normal, this is the precision parameter \code{k}.
#' @param mu A numeric vector representing the log-median of the latent size process. Explicitly, \code{exp(mu)} is the median of the size distribution \code{A}.
#' @param A A numeric vector representing the size data.
#' @param thr.family A character string indicating the threshold family distribution. Can be either \code{"gamma"} or \code{"logNormal"}.
#' @param ind_zeros_counts A logical vector indicating which elements of \code{A} are zeros.
#' @param sum_dens A logical value indicating whether to calculate the sum of the log-likelihood values (\code{TRUE}) or return the log-density values individually (\code{FALSE}).
#'
#' @return If \code{sum_dens = TRUE}, returns a numeric value representing the sum of log-likelihood values. If \code{sum_dens = FALSE}, returns a numeric vector of log-density values for each element of \code{A}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' cur_par <- c(2) # Shape parameter for gamma distribution
#' mu <- rnorm(10)
#' A <- rgamma(10, shape = 2, scale = 1)
#' ind_zeros_counts <- A == 0
#' thr.family <- "gamma"
#' sum_dens <- TRUE
#' loglik_value <- log_lik_thr(cur_par, mu, A, thr.family, ind_zeros_counts, sum_dens)
#'
log_lik_thr<- function(
    cur_par,
    mu, 
    A, 
    thr.family,
    ind_zeros_counts,
    sum_dens=TRUE) {
  #browser()
  if(thr.family=="gamma"){
    k<- cur_par[1]
    if (sum_dens==TRUE){
      log.post_A<- sum(ifelse(ind_zeros_counts, 0, dgamma(x = A, shape = k, scale = exp(mu)/k , log = TRUE)))
    } else {
      log.post_A<- ifelse(ind_zeros_counts, 0, dgamma(x = A, shape = k, scale = exp(mu)/k , log = TRUE))
    }
  } 
  if(thr.family=="logNormal"){
    k<- cur_par[1]
    if (sum_dens==TRUE){
      log.post_A<- sum(ifelse(ind_zeros_counts, 0, dnorm(log(A), mean=mu, sd=sqrt(1/k), log = TRUE)))
    } else {
      log.post_A<- ifelse(ind_zeros_counts, 0, dnorm(log(A), mean=mu, sd=sqrt(1/k), log = TRUE))
    }
  }
  return(log.post_A)
}

