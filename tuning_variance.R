#' Adaptive Tuning of Variance Components in MH and MALA Algorithms
#'
#' This function adjusts the tuning parameter `sigma2_adapt` in the Metropolis-Hastings (MH) or MALA algorithm to reach a target acceptance rate.
#'
#' @param index_MCMC_iter Integer. The index of the current MCMC iteration.
#' @param sigma2_adapt Numeric. The current tuning parameter in the MH and MALA algorithm, which determines the proposal variance.
#' @param target_accept Numeric. The target acceptance probability (e.g., 0.224 for Random Walk MH, 0.57 for MALA).
#' @param rate_adapt Numeric. The acceptance rate calculated at fixed intervals over MCMC iterations.
#' @param burn_in1 Integer. The length of the first burn-in period in the MCMC algorithm.
#' @param burn_in2 Integer. The length of the second burn-in period in the MCMC algorithm.
#' @param adapt Integer. The sequence interval at which `sigma2_adapt` is updated.
#' @param adpat_param Numeric. A scaling parameter for the adaptive adjustment of `sigma2_adapt`.
#' @param adapt_seq Vector. A sequence indicating specific MCMC iterations for updating `sigma2_adapt`.
#' @param lower.acc Numeric. The lower bound for the acceptance rate threshold (e.g., 0.15 for Random Walk MH, 0.50 for MALA).
#' @param upper.acc Numeric. The upper bound for the acceptance rate threshold (e.g., 0.30 for Random Walk MH, 0.65 for MALA).
#'
#' @details
#' The function adapts the variance parameter `sigma2_adapt` based on the current acceptance rate relative to a target acceptance rate. During the first burn-in phase, the variance is adjusted to meet the target acceptance rate directly. During the second burn-in phase, updates occur only if the acceptance rate is outside of specified bounds (`lower.acc` and `upper.acc`).
#'
#' @return
#' A numeric value representing the updated tuning parameter `sigma2_adapt`.
#'
#' @export
#'
#' @examples
#' # Initialize parameters for adaptive tuning
#' index_MCMC_iter <- 100
#' sigma2_adapt <- 0.1
#' target_accept <- 0.224
#' rate_adapt <- 25
#' burn_in1 <- 500
#' burn_in2 <- 500
#' adapt <- 100
#' adpat_param <- 10
#' adapt_seq <- seq(100, 1000, by = 100)
#' lower.acc <- 0.15
#' upper.acc <- 0.30
#'
#' # Update the variance tuning parameter
#' sigma2_adapt <- adpative_function(index_MCMC_iter, sigma2_adapt, target_accept, rate_adapt, burn_in1, burn_in2, adapt, adpat_param, adapt_seq, lower.acc, upper.acc)
adpative_function <- function(index_MCMC_iter,
                              sigma2_adapt,
                              target_accept,
                              rate_adapt,
                              burn_in1,
                              burn_in2,
                              adapt,
                              adpat_param,
                              adapt_seq,
                              lower.acc,
                              upper.acc) {
  if (index_MCMC_iter < burn_in1 + burn_in2) {
    if (index_MCMC_iter %in% adapt_seq) {
      if (index_MCMC_iter < burn_in1) {
        sigma2_adapt <- exp(((rate_adapt / adapt) - target_accept) / adpat_param) * sigma2_adapt
      } else {
        sigma2_adapt <- ifelse((((
          rate_adapt / adapt
        ) > upper.acc) | ((
          rate_adapt / adapt
        ) < lower.acc)), exp(((rate_adapt / adapt) - target_accept
        ) / adpat_param) * sigma2_adapt, sigma2_adapt)
      }
    }
  } else {
    sigma2_adapt <- sigma2_adapt
  }
  return(sigma2_adapt)
}
