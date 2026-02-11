#' Log-liikelihood of the use mark-distribution
#'
#' This function computes the log-likelihood of the size distribution  given the parameters, latent size process, and data. It supports three types of mark distributions: extended Generalized Pareto Distribution (eGPD), mixture of beta-GPD (bGPD), and mixture of truncated gamma-GPD (tgGPD).
#'
#' @param cur_par A numeric vector representing the hyperparameters of the continuous mark distributions. For example, \code{k} and \code{xi} in the case of eGPD.
#' @param mu A numeric vector representing the log-median of the latent size process. Explicitly, \code{exp(mu)} is the median of the size distribution \code{A}.
#' @param A A numeric vector representing the landslide size data.
#' @param mark_dist A character string indicating the choice of mark distribution for the size process. Can be `"eGPD"`, `"bGPD"`, or `"tgGPD"`.
#' @param threshold A numeric vector representing the threshold values for the size data \code{A}.
#' @param sum_dens A logical value indicating whether to calculate the sum of the log-likelihood (\code{TRUE}) or return the log-density values individually (\code{FALSE}).
#'
#' @return A numeric value representing the sum of log-likelihood values if \code{sum_dens = TRUE}, or a numeric vector of log-density values if \code{sum_dens = FALSE}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' cur_par <- c(2, 0.1) # example parameters for eGPD
#' mu <- rnorm(10)
#' A <- rgamma(10, shape = 2, scale = 1)
#' mark_dist <- "eGPD"
#' threshold <- rep(1, 10)
#' sum_dens <- TRUE
#' loglik_value <- log_lik_shape_mixture(cur_par, mu, A, mark_dist, threshold, sum_dens)
#'
log_lik_shape_mixture <- function(cur_par,
                                  mu,
                                  A,
                                  mark_dist,
                                  threshold,
                                  sum_dens = TRUE) {
  ## browser()
  n1 <- length(mu)
  if (mark_dist == "eGPD") {
    k <- cur_par[1]
    xi <- cur_par[2]
    sigma <- exp(mu) / evd::qgpd(0.5 ^ (1 / k), scale = 1, shape = xi)
    
    if (sum_dens == TRUE) {
      log.post_A <- rep(0, n1)
      ind.postv <- A > 0
      log.post_A[ind.postv] <- dEGPD1(
        x = A[ind.postv],
        k = k,
        xi = xi,
        sigma = sigma[ind.postv],
        log = TRUE
      )
      log.post_A <- sum(log.post_A)
      
    } else {
      log.post_A <- rep(0, n1)
      ind.postv <- A > 0
      log.post_A[ind.postv] <- dEGPD1(
        x = A[ind.postv],
        k = k,
        xi = xi,
        sigma = sigma[ind.postv],
        log = TRUE
      )
      
    }
  } else if (mark_dist == "bGPD") {
    k <- cur_par[1]
    if (sum_dens == TRUE) {
      log.post_A <- rep(0, n1)
      ind <- (A < threshold) & (A != 0)
      log.post_A[ind] <- dbeta((A[ind] / threshold[ind]),
                               shape1 = k,
                               shape2 = (threshold[ind] / exp(mu[ind]) - 1) * k,
                               log = TRUE
      ) - log(threshold[ind])
      log.post_A <- sum(log.post_A)
      
    } else {
      log.post_A <- rep(0, n1)
      ind <- (A < threshold) & (A != 0)
      log.post_A[ind] <- dbeta((A[ind] / threshold[ind]),
                               shape1 = k,
                               shape2 = (threshold[ind] / exp(mu[ind]) - 1) * k,
                               log = TRUE
      ) - log(threshold[ind])
    }
  } else if (mark_dist == "tgGPD") {
    k <- cur_par[1]
    if (sum_dens == TRUE) {
      log.post_A <- rep(0, n1)
      ind <- (A < threshold) & (A != 0)
      log.post_A[ind] <-  dgamma(A[ind],
                                 shape = k,
                                 scale = exp(mu[ind]) / k,
                                 log = TRUE) - pgamma(
                                   threshold[ind],
                                   shape = k,
                                   scale =  exp(mu[ind]) / k,
                                   log.p = TRUE
                                 )
      log.post_A <- sum(log.post_A)
    } else {
      log.post_A <- rep(0, n1)
      ind <- (A < threshold) & (A != 0)
      log.post_A[ind] <-  dgamma(A[ind],
                                 shape = k,
                                 scale = exp(mu[ind]) / k,
                                 log = TRUE) -
        pgamma(
          threshold[ind],
          shape = k,
          scale =  exp(mu[ind]) / k,
          log.p = TRUE
        )
      
    }
  }
  return(log.post_A)
}




#' Log-Likelihood for Generalized Pareto Distribution (GPD) Parameters
#'
#' This function computes the log-likelihood of a Generalized Pareto Distribution (GPD) for given parameters and data, conditioned on exceeding a specified threshold.
#'
#' @param cur_par A numeric vector containing the GPD parameters: \code{sigma.GP} (scale parameter) and \code{xi} (shape parameter).
#' @param A A numeric vector representing the size data for the GPD model.
#' @param threshold A numeric vector representing the threshold values, above which the GPD is considered.
#'
#' @return A numeric value representing the sum of the log-likelihood for all elements of \code{A} that exceed the corresponding \code{threshold}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' cur_par <- c(1, 0.2) # sigma.GP = 1, xi = 0.2
#' A <- rnorm(10, mean = 5)
#' threshold <- rep(4, 10)
#' loglik_value <- log_lik_GPD_param(cur_par, A, threshold)
#'
log_lik_GPD_param <- function(cur_par, A, threshold) {
  sigma.GP <- cur_par[1]
  xi <- cur_par[2]
  thr.exceed.ind <- A > threshold
  log.post_A <- sum(evd::dgpd(
    x = A[thr.exceed.ind],
    loc = threshold[thr.exceed.ind],
    scale = sigma.GP,
    shape = xi,
    log = TRUE
  ))
  return(log.post_A)
}
