#' Log-Likelihood at the Data Level for Indicator Model
#'
#' This function computes the log-likelihood of the size data \(A\) given the latent size process parameterized by \(mu\). Specifically, it calculates the likelihood of density \(A | theta_A, mu\).
#'
#' @param mu Numeric value representing the log-median of the latent size process. More explicitly, \(exp(mu)\) is the median of the size distribution \(A\).
#' @param A Numeric vector representing the size data.
#' @param sum_dens Logical value indicating whether to calculate the log-density individually or the sum of log-likelihoods. Defaults to \code{FALSE}.
#' @param ind_zeros_counts Optional parameter for handling zero counts in the data. This can be used to specify which counts to include in the likelihood calculation.
#'
#' @return A numeric vector containing the log-likelihood for each observation in \(A\).
#' @export
#'
#' @examples
#' # Example size data
#' A <- c(0, 1, 0, 1, 1)
#' mu <- log(0.5)  # Assuming the median size is 0.5
#'
#' # Calculate the log-likelihood
#' log_likelihood <- log_lik_indicator_model(mu, A)
log_lik_indicator_model_probit<- function(mu, A) {
  log.post_A<-  A * log(logistic(mu)) + (1-A) * log(1- logistic(mu))  # dbinom(A, size=1, prob = exp(mu)/(1+exp(mu)), log = TRUE)
  return(log.post_A)
}





#' Update the Latent Parameter \code{mu} for Indicator Model
#'
#' This function generates samples for the latent parameter \code{mu} from a full conditional truncated normal distribution for for Indicator Model.
#'
#' @param mu_mean A numeric vector representing the mean values for the latent parameter \code{mu}.
#' @param kappa.mu A numeric value representing the precision parameter for the latent parameter \code{mu}.
#' @param ind_zero A logical vector indicating which entries of \code{mu} should be sampled from the left truncated normal distribution (i.e., mean ≤ 0).
#' @param ind.NA A logical vector indicating which entries of \code{mu} correspond to missing values (NA) that should be sampled from the normal distribution without truncation.
#' @param CV A character string indicating the type of cross-validation. If it is set to "OOS", the function will sample for the NA indices as well.
#'
#' @return A numeric vector of the same length as \code{mu_mean} containing the sampled values of the latent parameter \code{mu}.
#' @export
#'
#' @examples
#' # Example parameters
#' mu_mean <- c(-1, 0.5, 1, -0.5)
#' kappa.mu <- 2
#' ind_zero <- c(TRUE, FALSE, FALSE, TRUE)
#' ind.NA <- c(FALSE, FALSE, TRUE, FALSE)
#' CV <- "OOS"
#'
#' # Update mu values
#' updated_mu <- update_mu(mu_mean, kappa.mu, ind_zero, ind.NA, CV)

update_mu <- function(mu_mean, kappa.mu, ind_zero, ind.NA, CV) {
  #browser()
  sim <- rep(NA, length(A))
  sim[ind_zero] <- truncnorm::rtruncnorm(
    n = 1,
    a = -Inf,
    b = 0,
    mean = mu_mean[ind_zero],
    sd = sqrt(1 / kappa.mu)
  )
  sim[!ind_zero] <- truncnorm::rtruncnorm(
    n = 1,
    a = 0,
    b = Inf,
    mean = mu_mean[!ind_zero],
    sd = sqrt(1 / kappa.mu)
  )
  if (CV == "OOS") {
    sim[ind.NA] <- truncnorm::rtruncnorm(
      n = 1,
      a = -Inf,
      b = Inf,
      mean = mu_mean[ind.NA],
      sd = sqrt(1 / kappa.mu)
    )
  }
  return(sim)
}
