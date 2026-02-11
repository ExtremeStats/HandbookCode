#' Log-Likelihood for Latent (log-linear predictor) \code{mu} in the Indicator Model
#'
#' This function computes the log-likelihood for the latent variable \code{mu} in the context of an indicator model.
#'
#' @param mu A numeric vector representing the latent variable for which the log-likelihood is being calculated.
#' @param intercept2 A numeric value representing the intercept of the model.
#' @param W2 A numeric vector of random effects or additional covariates in \code{mu}.
#' @param beta2 A numeric vector of coefficients for the covariates in \code{Z2}.
#' @param kappa_mu A numeric value representing the precision parameter for \code{mu}.
#' @param Z2 A matrix of covariates in \code{mu}.
#'
#' @return A numeric vector representing the log-likelihood for each value of \code{mu}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' mu <- rnorm(10)
#' intercept2 <- 0.5
#' W2 <- rnorm(10)
#' beta2 <- c(0.3, -0.2)
#' kappa_mu <- 2
#' Z2 <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' loglik_values <- log_lik_latent_mu_thr_ind(mu, intercept2, W2, beta2, kappa_mu, Z2)
#'
log_lik_latent_mu_thr_ind_logit<- function(mu, intercept2, W2, beta2,  kappa_mu, Z2){
  loglik<- -0.5*kappa_mu* mu^2 + kappa_mu * mu *(intercept2 + Z2 %*% beta2 + W2)
  return(loglik)
  
}



#' Log-Likelihood for Latent \code{mu} in the Joint Model
#'
#' This function calculates the log-likelihood for the latent variable \code{mu} in the context of a joint model.
#'
#' @param mu A numeric vector representing the latent variable for which the log-likelihood is being calculated.
#' @param intercept2 A numeric value representing the intercept of the model.
#' @param W1 A numeric vector of random effects.
#' @param W2 A numeric vector of random effects.
#' @param beta A numeric value representing sharing parameter used as the coefficient of \code{W1} on \code{mu}.
#' @param beta2 A numeric vector of coefficients for the covariates in \code{Z2}.
#' @param kappa_mu A numeric value representing the precision or scale parameter for \code{mu}.
#' @param Z2 A matrix of covariates influencing \code{mu}.
#'
#' @return A numeric vector representing the log-likelihood for each value of \code{mu}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' mu <- rnorm(10)
#' intercept2 <- 0.5
#' W1 <- rnorm(10)
#' W2 <- rnorm(10)
#' beta <- 0.3
#' beta2 <- c(0.3, -0.2)
#' kappa_mu <- 2
#' Z2 <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' loglik_values <- log_lik_latent_mu_JM(mu, intercept2, W1, W2, beta, beta2, kappa_mu, Z2)
#'
log_lik_latent_mu_JM <- function(mu,
                                 intercept2,
                                 W1,
                                 W2,
                                 beta,
                                 beta2,
                                 kappa_mu,
                                 Z2) {
  loglik <- -0.5 * kappa_mu * mu ^ 2 + kappa_mu * mu * (intercept2 + Z2 %*% beta2 + beta * W1 + W2)
  return(loglik)
  
}


#' Log-Likelihood for Latent \code{eta} in Joint Model
#'
#' This function calculates the log-likelihood for the latent variable \code{eta} in the context of a joint model.
#'
#' @param eta A numeric vector representing the latent variable for which the log-likelihood is being calculated.
#' @param intercept1 A numeric value representing the intercept for the first part of the model.
#' @param beta1 A numeric vector of coefficients for the covariates in \code{Z1}.
#' @param kappa_eta A numeric value representing the precision or scale parameter for \code{eta}.
#' @param W1 A numeric vector of random effects or covariates influencing \code{eta}.
#' @param Z1 A matrix of covariates influencing \code{eta}.
#'
#' @return A numeric vector representing the log-likelihood for each value of \code{eta}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' eta <- rnorm(10)
#' intercept1 <- 0.5
#' W1 <- rnorm(10)
#' beta1 <- c(0.3, -0.2)
#' kappa_eta <- 2
#' Z1 <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' loglik_values <- log_lik_latent_eta_JM(eta, intercept1, beta1, kappa_eta, W1, Z1)
#'
log_lik_latent_eta_JM <- function(eta, intercept1, beta1, kappa_eta, W1, Z1) {
  loglik <- -0.5 * kappa_eta * eta ^ 2 + kappa_eta * eta  * (intercept1 + Z1 %*%
                                                               beta1 + W1)
  return(loglik)
}

#' Log-Likelihood for at the count data level in Joint Model
#'
#' This function calculates the log-likelihood for the count data.
#'
#' @param Y A numeric vector representing the response variable in the model.
#' @param eta A numeric vector representing the linear predictor or latent variable associated with the model.
#'
#' @return A numeric vector representing the log-likelihood for each value of \code{eta}.
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' Y <- rpois(10, lambda = 3)
#' eta <- rnorm(10)
#' loglik_values <- log_lik_eta_JM(Y, eta)
#'
log_lik_eta_JM <- function(Y, eta) {
  loglik <-  Y * eta - exp(eta)
  return(loglik)
}