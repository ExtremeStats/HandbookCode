#' Simulate the Precision Parameter of iCAR Random Effects for Threshold Model
#'
#' This function simulates the parameter \code{kappa_w2}, the precision parameter of the intrinsic Conditional Autoregressive (iCAR) model for a threshold model, from a full conditional Gamma distribution.
#'
#' @param W2 A numeric vector representing the latent variable (iCAR random effects) for the model.
#' @param node1 An integer vector indicating the first set of connected nodes (indices) involved in the computation.
#' @param node2 An integer vector indicating the second set of connected nodes (indices) involved in the computation.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution of \code{kappa_w2}, specifically the shape and rate parameters (\code{kappa_w2[1]} for shape and \code{kappa_w2[2]} for rate).
#'
#' @return A simulated value for \code{kappa_w2} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variable
#' W2 <- rnorm(100)
#' node1 <- sample(1:100, 50, replace = TRUE)
#' node2 <- sample(1:100, 50, replace = TRUE)
#' hyper_fixed <- list(kappa_w2 = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_w2 for threshold model
#' kappa_sim <- kappa_w2_sim_thr(W2, node1, node2, hyper_fixed)
#'
kappa_w2_sim_thr <- function(W2, node1, node2, hyper_fixed) {
  N <- length(W2)
  sim_kappa_w2 <- rgamma(
    1,
    shape = hyper_fixed$kappa_w2[1] + 0.5 * (N - 1),
    rate = hyper_fixed$kappa_w2[2] + 0.5 * sum((W2[node1] -
                                                  W2[node2]) ^ 2)
  )
  return(sim_kappa_w2)
}


#' Simulate the Precision Parameter of \code{mu} in Threshold Model
#'
#' This function simulates the parameter \code{kappa_mu}, representing the precision of \code{mu} in a threshold model, from a full conditional Gamma distribution.
#'
#' @param mu A numeric vector representing the current state of the latent mean parameters.
#' @param intercept2 A numeric value representing the intercept term in the model.
#' @param beta2 A numeric vector of coefficients corresponding to the covariate matrix \code{Z2}.
#' @param W2 A numeric vector representing the random effects for the model.
#' @param Z2 A matrix of covariates, where each row corresponds to an observation.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution of \code{kappa_mu}, specifically the shape and rate parameters (\code{kappa_mu[1]} for shape and \code{kappa_mu[2]} for rate).
#'
#' @return A simulated value for \code{kappa_mu} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample parameters
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta2 <- rnorm(5)
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(500), ncol=5)
#' hyper_fixed <- list(kappa_mu = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_mu for threshold model
#' kappa_sim <- kappa_mu_sim_thr(mu, intercept2, beta2, W2, Z2, hyper_fixed)
#'
kappa_mu_sim_thr <- function(mu, intercept2, beta2, W2, Z2, hyper_fixed) {
  n2 <- length(mu)
  sim_kappa_mu <- rgamma(1,
                         shape = hyper_fixed$kappa_mu[1] + 0.5 * n2,
                         rate = hyper_fixed$kappa_mu[2] + 0.5 * (sum((
                           mu - intercept2 - Z2 %*% beta2 - W2
                         ) ^ 2)))
  return(sim_kappa_mu)
}


##' Simulate the Intercept for a Threshold Model
#'
#' This function simulates the intercept for a threshold model based on provided parameters.
#'
#' @param mu A numeric vector of means.
#' @param beta2 A vector of regression coefficients associated with covariates in `Z2`.
#' @param beta A numeric vector of coefficients for random effects.
#' @param W2 A vector or matrix representing random effects.
#' @param kappa_mu A numeric value representing the precision for the mean `mu`.
#' @param Z2 A design matrix for covariates associated with `beta2`.
#' @param hyper_fixed A list containing hyperparameters for the model, particularly the fixed precision for the intercept.
#'
#' @return A simulated value of the intercept for the threshold model.
#' @export
#'
#' @examples
#' # Example usage:
#' mu <- rnorm(100, 0, 1)
#' beta2 <- rnorm(5, 0, 1)
#' beta <- rnorm(5, 0, 1)
#' W2 <- rnorm(100, 0, 1)
#' kappa_mu <- 0.5
#' Z2 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(intercept2 = 0.01)
#' intercept2_sim_thr(mu, beta2, beta, W2, kappa_mu, Z2, hyper_fixed)
intercept2_sim_thr <- function(mu,
                               beta2,
                               beta,
                               W2,
                               kappa_mu,
                               Z2,
                               hyper_fixed) {
  n2 <- length(mu)
  # prec.intercept2<-hyper_fixed[11]+n2*kappa_mu
  # mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)/prec.intercept2
  prec.intercept2 <- hyper_fixed$intercept2 + n2 * kappa_mu
  mean.intercept2 <- kappa_mu * sum(mu - Z2 %*% beta2 - W2) / prec.intercept2
  sim_intercept2 <- rnorm(n = 1,
                          mean = mean.intercept2,
                          sd = sqrt(1 / prec.intercept2))
  return(sim_intercept2)
}


##' Simulate the Intercept Parameter for Threshold Model
#'
#' This function simulates the parameter \code{intercept2}, which is the intercept term in a threshold model, based on the full conditional distribution derived from the data and model parameters.
#'
#' @param mu A numeric vector representing the current state of the latent mean parameters.
#' @param beta2 A numeric vector of coefficients for the covariate matrix \code{Z2}.
#' @param W2 A numeric vector representing the random effects for the model.
#' @param kappa_mu A numeric value representing the precision of the mean effects in the threshold model.
#' @param Z2 A matrix of covariates, where each row corresponds to an observation.
#' @param hyper_fixed A list containing the hyperparameter for \code{intercept2}, which includes its prior precision (\code{hyper_fixed$intercept2}).
#'
#' @return A simulated value for \code{intercept2} from its full conditional normal distribution.
#' @export
#'
#' @examples
#' # Sample parameters
#' mu <- rnorm(100)
#' beta2 <- rnorm(5)
#' W2 <- rnorm(100)
#' kappa_mu <- 2
#' Z2 <- matrix(rnorm(500), ncol=5)
#' hyper_fixed <- list(intercept2 = 0.1)
#'
#' # Simulate intercept2 for threshold model
#' intercept2_sim <- intercept2_sim_thr(mu, beta2, W2, kappa_mu, Z2, hyper_fixed)
#'
beta2_sim_thr <- function(mu,
                          intercept2,
                          beta,
                          kappa_mu,
                          W2,
                          Z2,
                          Z2.crossprd,
                          hyper_fixed) {
  q <- ncol(Z2)
  latent.cov.inv <- hyper_fixed$beta2 * diag(1, q) + kappa_mu * Z2.crossprd
  latent.mean.part <- kappa_mu * (t(Z2) %*% (mu - intercept2 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
}

#' Simulate the Precision Parameter \code{kappa_w2} for Threshold Indicator Models with probit link
#'
#' This function simulates the precision parameter \code{kappa_w2} for threshold indicator models based on a full conditional Gamma distribution.
#'
#' @param W2 A numeric vector representing the latent variable (iCAR) random effects for the model.
#' @param node1 A numeric vector indicating the first set of nodes (indices) involved in the computation.
#' @param node2 A numeric vector indicating the second set of nodes (indices) involved in the computation.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_w2}.
#'
#' @return A simulated value for \code{kappa_w2} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variable
#' W2 <- rnorm(100)
#' node1 <- sample(1:100, 50, replace = TRUE)
#' node2 <- sample(1:100, 50, replace = TRUE)
#' hyper_fixed <- list(kappa_w2 = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_w2 for threshold indicator model
#' kappa_sim <- kappa_w2_sim_thr_ind_probit(W2, node1, node2, hyper_fixed)
#'
kappa_w2_sim_thr_ind_probit <- function(W2, node1, node2, hyper_fixed) {
  N <- length(W2)
  sim_kappa_w2 <- rgamma(
    1,
    shape = hyper_fixed$kappa_w2[1] + 0.5 * (N - 1),
    rate = hyper_fixed$kappa_w2[2] + 0.5 * sum((W2[node1] -
                                                  W2[node2]) ^ 2)
  )
  return(sim_kappa_w2)
}

#' Simulate the Precision Parameter \code{kappa_w2} for Threshold Indicator Models with logit link
#'
#' This function simulates the precision parameter \code{kappa_w2} for threshold indicator models based on a full conditional Gamma distribution.
#'
#' @param W2 A numeric vector representing the latent variable (iCAR) random effects for the model.
#' @param node1 A numeric vector indicating the first set of nodes (indices) involved in the computation.
#' @param node2 A numeric vector indicating the second set of nodes (indices) involved in the computation.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_w2}.
#'
#' @return A simulated value for \code{kappa_w2} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variable
#' W2 <- rnorm(100)
#' node1 <- sample(1:100, 50, replace = TRUE)
#' node2 <- sample(1:100, 50, replace = TRUE)
#' hyper_fixed <- list(kappa_w2 = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_w2 for threshold indicator model
#' kappa_sim <- kappa_w2_sim_thr_ind_logit(W2, node1, node2, hyper_fixed)
#' 
kappa_w2_sim_thr_ind_logit<-function(W2, node1, node2, hyper_fixed){
  N<-length(W2)
  sim_kappa_w2<-rgamma(1, shape=hyper_fixed$kappa_w2[1] + 0.5*(N-1), 
                       rate= hyper_fixed$kappa_w2[2] + 0.5*sum((W2[node1]-W2[node2])^2))
  return(sim_kappa_w2)
}




#' Simulate the Precision Parameter \code{kappa_mu} for Threshold Indicator Models with probit link
#'
#' This function simulates the precision parameter \code{kappa_mu} for threshold indicator models from a full conditional Gamma distribution.
#'
#' @param mu A numeric vector representing the latent mean parameter of the model.
#' @param intercept2 A numeric value indicating the intercept term in the model.
#' @param beta2 A numeric vector representing the coefficients of covariates in the model.
#' @param W2 A numeric vector indicating the latent (iCAR) random effects.
#' @param Z2 A numeric matrix of covariates corresponding to the latent process \code{mu}.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_mu}.
#'
#' @return A simulated value for \code{kappa_mu} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variables and covariates
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta2 <- rnorm(10)
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(1000), ncol = 10)
#' hyper_fixed <- list(kappa_mu = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_mu for threshold indicator model
#' kappa_mu_sim <- kappa_mu_sim_thr_ind_probit(mu, intercept2, beta2, W2, Z2, hyper_fixed)
#'
kappa_mu_sim_thr_ind_probit <- function(mu, intercept2, beta2, W2, Z2, hyper_fixed) {
  n2 <- length(mu)
  sim_kappa_mu <- rgamma(1,
                         shape = hyper_fixed$kappa_mu[1] + 0.5 * n2,
                         rate = hyper_fixed$kappa_mu[2] + 0.5 * (sum((
                           mu - intercept2 - Z2 %*% beta2 - W2
                         ) ^ 2)))
  return(sim_kappa_mu)
}


#' Simulate the Precision Parameter \code{kappa_mu} for Threshold Indicator Models with logit link
#'
#' This function simulates the precision parameter \code{kappa_mu} for threshold indicator models from a full conditional Gamma distribution.
#'
#' @param mu A numeric vector representing the latent mean parameter of the model.
#' @param intercept2 A numeric value indicating the intercept term in the model.
#' @param beta2 A numeric vector representing the coefficients of covariates in the model.
#' @param W2 A numeric vector indicating the latent (iCAR) random effects.
#' @param Z2 A numeric matrix of covariates corresponding to the latent process \code{mu}.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_mu}.
#'
#' @return A simulated value for \code{kappa_mu} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variables and covariates
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta2 <- rnorm(10)
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(1000), ncol = 10)
#' hyper_fixed <- list(kappa_mu = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_mu for threshold indicator model
#' kappa_mu_sim <- kappa_mu_sim_thr_ind_logit(mu, intercept2, beta2, W2, Z2, hyper_fixed)
#'
kappa_mu_sim_thr_ind_logit<-function(mu, intercept2, beta2, W2, Z2, hyper_fixed){
  n2<-length(mu)
  sim_kappa_mu<- rgamma(1, shape=hyper_fixed$kappa_mu[1] + 0.5 * n2, 
                        rate= hyper_fixed$kappa_mu[2] + 0.5 * (sum((mu-intercept2- Z2%*%beta2 -W2)^2)))
  return(sim_kappa_mu)
}



#' Simulate the Intercept \code{intercept2} for Threshold Indicator Models with probit link
#'
#' This function simulates the intercept term \code{intercept2} for threshold indicator models using a normal full conditional distribution.
#'
#' @param mu A numeric vector representing the latent mean parameter of the model.
#' @param beta2 A numeric vector representing the coefficients of covariates in the model.
#' @param W2 A numeric vector indicating the latent (iCAR) random effects.
#' @param kappa_mu A numeric value representing the precision parameter for the latent process \code{mu}.
#' @param Z2 A numeric matrix of covariates corresponding to the latent process \code{mu}.
#' @param hyper_fixed A list containing hyperparameters for the intercept, specifically the precision of the normal prior for \code{intercept2}.
#'
#' @return A simulated value for \code{intercept2} from a normal distribution.
#' @export
#'
#' @examples
#' # Sample latent variables and covariates
#' mu <- rnorm(100)
#' beta2 <- rnorm(10)
#' W2 <- rnorm(100)
#' kappa_mu <- 1.5
#' Z2 <- matrix(rnorm(1000), ncol = 10)
#' hyper_fixed <- list(intercept2 = 0.1)
#'
#' # Simulate intercept2 for threshold indicator model
#' intercept2_sim <- intercept2_sim_thr_ind_probit(mu, beta2, W2, kappa_mu, Z2, hyper_fixed)
intercept2_sim_thr_ind_probit <- function(mu, beta2, W2, kappa_mu, Z2, hyper_fixed) {
  n2 <- length(mu)
  prec.intercept2 <- hyper_fixed$intercept2 + n2 * kappa_mu
  mean.intercept2 <- kappa_mu * sum(mu - Z2 %*% beta2 - W2) / prec.intercept2
  sim_intercept2 <- rnorm(n = 1,
                          mean = mean.intercept2,
                          sd = sqrt(1 / prec.intercept2))
  return(sim_intercept2)
}


#' Simulate the Intercept \code{intercept2} for Threshold Indicator Models with logit link
#'
#' This function simulates the intercept term \code{intercept2} for threshold indicator models using a normal full conditional distribution.
#'
#' @param mu A numeric vector representing the latent mean parameter of the model.
#' @param beta2 A numeric vector representing the coefficients of covariates in the model.
#' @param W2 A numeric vector indicating the latent (iCAR) random effects.
#' @param kappa_mu A numeric value representing the precision parameter for the latent process \code{mu}.
#' @param Z2 A numeric matrix of covariates corresponding to the latent process \code{mu}.
#' @param hyper_fixed A list containing hyperparameters for the intercept, specifically the precision of the normal prior for \code{intercept2}.
#'
#' @return A simulated value for \code{intercept2} from a normal distribution.
#' @export
#'
#' @examples
#' # Sample latent variables and covariates
#' mu <- rnorm(100)
#' beta2 <- rnorm(10)
#' W2 <- rnorm(100)
#' kappa_mu <- 1.5
#' Z2 <- matrix(rnorm(1000), ncol = 10)
#' hyper_fixed <- list(intercept2 = 0.1)
#'
#' # Simulate intercept2 for threshold indicator model
#' intercept2_sim <- intercept2_sim_thr_ind_logit(mu, beta2, W2, kappa_mu, Z2, hyper_fixed)
#' 
intercept2_sim_thr_ind_logit<-function(mu, beta2, beta,W2, kappa_mu, Z2, hyper_fixed){
  n2<-length(mu)
  # prec.intercept2<-hyper_fixed[11]+n2*kappa_mu
  # mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)/prec.intercept2
  prec.intercept2<-hyper_fixed$intercept2 +n2*kappa_mu
  mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2- W2)/prec.intercept2
  sim_intercept2<-rnorm(n=1, mean=mean.intercept2, sd=sqrt(1/prec.intercept2))
  return(sim_intercept2)
}


#' Simulate the Coefficients \code{beta2} for Threshold Indicator Models with probit link
#'
#' This function simulates the regression coefficients \code{beta2} for threshold indicator models using a normal full conditional distribution. It performs a precision-based update by considering both prior and data contributions.
#'
#' @param mu A numeric vector representing the latent mean parameter of the model.
#' @param intercept2 A numeric value representing the intercept term in the model.
#' @param beta A numeric vector of initial values for \code{beta2}.
#' @param kappa_mu A numeric value representing the precision parameter for the latent process \code{mu}.
#' @param W2 A numeric vector indicating the latent (iCAR) random effects.
#' @param Z2 A numeric matrix of covariates corresponding to the latent process \code{mu}.
#' @param Z2.crossprd A pre-computed cross-product matrix of \code{Z2} \(\code{t(Z2) Z2}\) for computational efficiency.
#' @param hyper_fixed A list containing hyperparameters, specifically the precision of the normal prior for \code{beta2}.
#'
#' @return A numeric vector of proposed values for \code{beta2}.
#' @export
#'
#' @examples
#' # Sample latent variables and covariates
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta <- rnorm(10)
#' kappa_mu <- 1.5
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(1000), ncol = 10)
#' Z2.crossprd <- crossprod(Z2) # pre-computed cross-product
#' hyper_fixed <- list(beta2 = 0.1)
#'
#' # Simulate beta2 for threshold indicator model
#' beta2_sim <- beta2_sim_thr_ind_probit(mu, intercept2, beta, kappa_mu, W2, Z2, Z2.crossprd, hyper_fixed)

beta2_sim_thr_ind_probit <- function(mu,
                              intercept2,
                              beta,
                              kappa_mu,
                              W2,
                              Z2,
                              Z2.crossprd,
                              hyper_fixed) {
  q <- ncol(Z2)
  latent.cov.inv <- hyper_fixed$beta2 * diag(1, q) + kappa_mu * Z2.crossprd
  latent.mean.part <- kappa_mu * (t(Z2) %*% (mu - intercept2 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
  
}


#' Simulate the Coefficients \code{beta2} for Threshold Indicator Models with logit link
#'
#' This function simulates the regression coefficients \code{beta2} for threshold indicator models using a normal full conditional distribution. It performs a precision-based update by considering both prior and data contributions.
#'
#' @param mu A numeric vector representing the latent mean parameter of the model.
#' @param intercept2 A numeric value representing the intercept term in the model.
#' @param beta A numeric vector of initial values for \code{beta2}.
#' @param kappa_mu A numeric value representing the precision parameter for the latent process \code{mu}.
#' @param W2 A numeric vector indicating the latent (iCAR) random effects.
#' @param Z2 A numeric matrix of covariates corresponding to the latent process \code{mu}.
#' @param Z2.crossprd A pre-computed cross-product matrix of \code{Z2} \(\code{t(Z2) Z2}\) for computational efficiency.
#' @param hyper_fixed A list containing hyperparameters, specifically the precision of the normal prior for \code{beta2}.
#'
#' @return A numeric vector of proposed values for \code{beta2}.
#' @export
#'
#' @examples
#' # Sample latent variables and covariates
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta <- rnorm(10)
#' kappa_mu <- 1.5
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(1000), ncol = 10)
#' Z2.crossprd <- crossprod(Z2) # pre-computed cross-product
#' hyper_fixed <- list(beta2 = 0.1)
#'
#' # Simulate beta2 for threshold indicator model
#' beta2_sim <- beta2_sim_thr_ind_probit(mu, intercept2, beta, kappa_mu, W2, Z2, Z2.crossprd, hyper_fixed)
beta2_sim_thr_ind_logit<-function(mu, intercept2, beta, kappa_mu,  W2, Z2, Z2.crossprd, hyper_fixed){
  q<-ncol(Z2)
  latent.cov.inv<- hyper_fixed$beta2 * diag(1,q) + kappa_mu * Z2.crossprd 
  latent.mean.part<- kappa_mu * (t(Z2)%*%(mu-intercept2 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
}




#' Simulate Precision Parameter \code{kappa_eta} of \code{eta} for Joint Model
#'
#' This function simulates the precision parameter \code{kappa_eta} for a joint model from a Gamma full conditional distribution.
#'
#' @param eta A numeric vector representing the latent variable for the model.
#' @param intercept1 A numeric value representing the intercept term in the joint model.
#' @param beta1 A numeric vector of regression coefficients associated with the covariates \code{Z1}.
#' @param W1 A numeric vector representing the latent random effects for the joint model.
#' @param Z1 A numeric matrix of covariates associated with \code{eta}.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_eta}.
#'
#' @return A simulated value for \code{kappa_eta} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variable and covariates
#' eta <- rnorm(100)
#' intercept1 <- 0.5
#' beta1 <- rnorm(10)
#' W1 <- rnorm(100)
#' Z1 <- matrix(rnorm(1000), ncol = 10)
#' hyper_fixed <- list(kappa_eta = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_eta for joint model
#' kappa_sim <- kappa_eta_sim_JM(eta, intercept1, beta1, W1, Z1, hyper_fixed)

kappa_eta_sim_JM <- function(eta, intercept1, beta1, W1, Z1, hyper_fixed) {
  n1 <- length(eta)
  sim_kappa_eta <- rgamma(1,
                          shape = hyper_fixed$kappa_eta[1] + 0.5 * n1,
                          rate = hyper_fixed$kappa_eta[2] + 0.5 * (sum((
                            eta - intercept1 - Z1 %*% beta1 - W1
                          ) ^ 2)))
  
  return(sim_kappa_eta)
}


#' Simulate Precision Parameter \code{kappa_w1} of iCAR \code{W_1} for Joint Model
#'
#' This function simulates the precision parameter \code{kappa_w1} for a joint model from a full conditional Gamma distribution.
#'
#' @param W1 A numeric vector representing the latent random effects for the joint model.
#' @param node1 A numeric vector indicating the first set of nodes (indices) involved in the computation.
#' @param node2 A numeric vector indicating the second set of nodes (indices) involved in the computation.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_w1}.
#'
#' @return A simulated value for \code{kappa_w1} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variable and nodes
#' W1 <- rnorm(100)
#' node1 <- sample(1:100, 50, replace = TRUE)
#' node2 <- sample(1:100, 50, replace = TRUE)
#' hyper_fixed <- list(kappa_w1 = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_w1 for joint model
#' kappa_sim <- kappa_w1_sim_JM(W1, node1, node2, hyper_fixed)
#'

kappa_w1_sim_JM <- function(W1, node1, node2, hyper_fixed) {
  # browser()
  N <- length(W1)
  sim_kappa_w1 <- rgamma(
    1,
    shape = hyper_fixed$kappa_w1[1] + 0.5 * (N - 1),
    rate = hyper_fixed$kappa_w1[2] + 0.5 * sum((W1[node1] -
                                                  W1[node2]) ^ 2)
  )
  return(sim_kappa_w1)
}


#' Simulate Precision Parameter \code{kappa_w2} of iCAR \code{W_2} for Joint Model
#'
#' This function simulates the precision parameter \code{kappa_w2} for a joint model from a full conditional Gamma distribution.
#'
#' @param W2 A numeric vector representing the latent random effects for the joint model.
#' @param node1 A numeric vector indicating the first set of nodes (indices) involved in the computation.
#' @param node2 A numeric vector indicating the second set of nodes (indices) involved in the computation.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_w2}.
#'
#' @return A simulated value for \code{kappa_w2} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Sample latent variable and nodes
#' W2 <- rnorm(100)
#' node1 <- sample(1:100, 50, replace = TRUE)
#' node2 <- sample(1:100, 50, replace = TRUE)
#' hyper_fixed <- list(kappa_w2 = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_w2 for joint model
#' kappa_sim <- kappa_w2_sim_JM(W2, node1, node2, hyper_fixed)
#'

kappa_w2_sim_JM <- function(W2, node1, node2, hyper_fixed) {
  N <- length(W2)
  sim_kappa_w2 <- rgamma(
    1,
    shape = hyper_fixed$kappa_w2[1] + 0.5 * (N - 1),
    rate = hyper_fixed$kappa_w2[2] + 0.5 * sum((W2[node1] -
                                                  W2[node2]) ^ 2)
  )
  return(sim_kappa_w2)
}


##' Simulate Precision Parameter \code{kappa_mu} of \code{mu} for Joint Model
#'
#' This function simulates the precision parameter \code{kappa_mu} for a joint model from a full conditional Gamma distribution.
#'
#' @param mu A numeric vector representing the response variable.
#' @param intercept2 A numeric value representing the intercept term in the model.
#' @param beta2 A numeric vector of coefficients for covariates in \code{Z2}.
#' @param beta A numeric value representing the coefficient associated with random effect \code{W1}.
#' @param W1 A numeric vector representing the first set of latent random effects in the model.
#' @param W2 A numeric vector representing the second set of latent random effects in the model.
#' @param Z2 A numeric matrix of covariates corresponding to the response variable.
#' @param hyper_fixed A list containing hyperparameters for the Gamma prior distribution, specifically the shape and rate parameters for \code{kappa_mu}.
#'
#' @return A simulated value for \code{kappa_mu} from the Gamma distribution.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta2 <- rnorm(5)
#' beta <- 1.2
#' W1 <- rnorm(100)
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(kappa_mu = c(shape = 2, rate = 1))
#'
#' # Simulate kappa_mu for joint model
#' kappa_sim <- kappa_mu_sim_JM(mu, intercept2, beta2, beta, W1, W2, Z2, hyper_fixed)

kappa_mu_sim_JM <- function(mu,
                            intercept2,
                            beta2,
                            beta,
                            W1,
                            W2,
                            Z2,
                            hyper_fixed) {
  n2 <- length(mu)
  sim_kappa_mu <- rgamma(1,
                         shape = hyper_fixed$kappa_mu[1] + 0.5 * n2,
                         rate = hyper_fixed$kappa_mu[2] + 0.5 * (sum((mu -
                                                                        intercept2 - Z2 %*% beta2 - beta * W1 - W2) ^ 2
                         )))
  return(sim_kappa_mu)
}



#' Simulate Sharing Coefficient \code{beta} of Joint Model
#'
#' This function simulates the sharing coefficient \code{beta} for a joint model using its full conditional distribution, which is normally distributed with precision based on covariates and latent variables.
#'
#' @param mu A numeric vector representing the response variable.
#' @param intercept2 A numeric value representing the intercept term in the model.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param beta2 A numeric vector of coefficients for covariates in \code{Z2}.
#' @param W1 A numeric vector representing the first set of latent random effects in the model.
#' @param W2 A numeric vector representing the second set of latent random effects in the model.
#' @param Z2 A numeric matrix of covariates corresponding to the response variable.
#' @param hyper_fixed A list containing hyperparameters for the Normal prior distribution, specifically the precision for \code{beta}.
#'
#' @return A simulated value for \code{beta} from its full conditional normal distribution.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' kappa_mu <- 2
#' beta2 <- rnorm(5)
#' W1 <- rnorm(100)
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(beta = 1.5)
#'
#' # Simulate beta for joint model
#' beta_sim <- beta_sim_JM(mu, intercept2, kappa_mu, beta2, W1, W2, Z2, hyper_fixed)

beta_sim_JM <- function(mu,
                        intercept2,
                        kappa_mu,
                        beta2,
                        W1,
                        W2,
                        Z2,
                        hyper_fixed) {
  #prec_beta<-kappa_mu * (t(W1) %*% t.A2.A2 %*% W1) + hyper_fixed[9]
  #mean_beta<-kappa_mu*(t(W1)%*% t(A2)%*% (mu-intercept2-Z2%*%beta2- W2))/prec_beta  ## may be we can simplify it later
  prec_beta <- kappa_mu * sum(W1 ^ 2) + hyper_fixed$beta
  mean_beta <- kappa_mu * (sum(W1 * (mu - intercept2 - Z2 %*% beta2 - W2))) /
    prec_beta  ## may be we can simplify it later
  sim_beta <- rnorm(n = 1,
                    mean = mean_beta,
                    sd = sqrt(1 / prec_beta))
  return(sim_beta)
}


#' Simulate Intercept \code{intercept1} of \code{eta} for Joint Model
#'
#' This function simulates the intercept parameter \code{intercept1} for a joint model based on its full conditional distribution, which is normally distributed with precision derived from the latent variables and covariates.
#'
#' @param eta A numeric vector representing the response variable for the first part of the model.
#' @param beta1 A numeric vector of coefficients for the covariates in \code{Z1}.
#' @param W1 A numeric vector representing the latent random effects for the first part of the model.
#' @param kappa_eta A numeric value representing the precision parameter associated with \code{eta}.
#' @param Z1 A numeric matrix of covariates corresponding to the response variable \code{eta}.
#' @param hyper_fixed A list containing hyperparameters, specifically the prior precision for \code{intercept1}.
#'
#' @return A simulated value for \code{intercept1} from its full conditional normal distribution.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' eta <- rnorm(100)
#' beta1 <- rnorm(5)
#' W1 <- rnorm(100)
#' kappa_eta <- 2
#' Z1 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(intercept1 = 1.5)
#'
#' # Simulate intercept1 for joint model
#' intercept1_sim <- intercept1_sim_JM(eta, beta1, W1, kappa_eta, Z1, hyper_fixed)

intercept1_sim_JM <- function(eta, beta1, W1, kappa_eta, Z1, hyper_fixed) {
  n1 <- length(eta)
  prec.intercept1 <- hyper_fixed$intercept1 + n1 * kappa_eta
  mean.intercept1 <- kappa_eta * sum((eta - Z1 %*% beta1 - W1)) / prec.intercept1
  sim_intercept1 <- rnorm(n = 1,
                          mean = mean.intercept1,
                          sd = sqrt(1 / prec.intercept1))
  return(sim_intercept1)
}


#' Simulate Intercept \code{intercept2} of \code{mu} for Joint Model
#'
#' This function simulates the intercept parameter \code{intercept2} for a joint model based on its full conditional distribution, which is normally distributed with precision derived from the latent variables, covariates, and random effects.
#'
#' @param mu A numeric vector representing the response variable for the second part of the model.
#' @param beta2 A numeric vector of coefficients for the covariates in \code{Z2}.
#' @param beta A numeric value representing the effect of the latent variable \code{W1} on \code{mu}.
#' @param W1 A numeric vector representing the latent random effects for the first part of the model.
#' @param W2 A numeric vector representing the latent random effects for the second part of the model.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param Z2 A numeric matrix of covariates corresponding to the response variable \code{mu}.
#' @param hyper_fixed A list containing hyperparameters, specifically the prior precision for \code{intercept2}.
#'
#' @return A simulated value for \code{intercept2} from its full conditional normal distribution.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' mu <- rnorm(100)
#' beta2 <- rnorm(5)
#' beta <- 0.5
#' W1 <- rnorm(100)
#' W2 <- rnorm(100)
#' kappa_mu <- 2
#' Z2 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(intercept1 = 1.5)
#'
#' # Simulate intercept2 for joint model
#' intercept2_sim <- intercept2_sim_JM(mu, beta2, beta, W1, W2, kappa_mu, Z2, hyper_fixed)

intercept2_sim_JM <- function(mu,
                              beta2,
                              beta,
                              W1,
                              W2,
                              kappa_mu,
                              Z2,
                              hyper_fixed) {
  n2 <- length(mu)
  # prec.intercept2<-hyper_fixed[11]+n2*kappa_mu
  # mean.intercept2<-kappa_mu*sum(mu-Z2%*%beta2-beta*(A2%*%W1)-A2%*%W2)/prec.intercept2
  prec.intercept2 <- hyper_fixed$intercept1 + n2 * kappa_mu
  mean.intercept2 <- kappa_mu * sum(mu - Z2 %*% beta2 - beta * W1 - W2) /
    prec.intercept2
  sim_intercept2 <- rnorm(n = 1,
                          mean = mean.intercept2,
                          sd = sqrt(1 / prec.intercept2))
  return(sim_intercept2)
}


#' Simulate Coefficients \code{beta1} of \code{eta} for Joint Model
#'
#' This function simulates the coefficients \code{beta1} for a joint model based on its full conditional distribution, which is derived from the latent variables and covariates.
#'
#' @param eta A numeric vector representing the response variable for the first part of the model.
#' @param intercept1 A numeric value representing the intercept for the first part of the model.
#' @param W1 A numeric vector representing the latent random effects for the first part of the model.
#' @param kappa_eta A numeric value representing the precision parameter associated with \code{eta}.
#' @param Z1 A numeric matrix of covariates corresponding to the response variable \code{eta}.
#' @param hyper_fixed A list containing hyperparameters, specifically the prior precision for \code{beta1}.
#'
#' @return A numeric vector of simulated values for \code{beta1} from its full conditional distribution.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' eta <- rnorm(100)
#' intercept1 <- 0.5
#' W1 <- rnorm(100)
#' kappa_eta <- 2
#' Z1 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(beta1 = 1)
#'
#' # Simulate beta1 for joint model
#' beta1_sim <- beta1_sim_JM(eta, intercept1, W1, kappa_eta, Z1, hyper_fixed)

beta1_sim_JM <- function(eta,
                         intercept1,
                         W1,
                         kappa_eta,
                         Z1,
                         hyper_fixed) {
  p <- ncol(Z1)
  latent.cov.inv <- hyper_fixed$beta1 * diag(1, p) + kappa_eta * (t(Z1) %*%
                                                                    Z1)
  latent.mean.part <- kappa_eta * (t(Z1) %*% (eta - intercept1 - W1))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(p)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
  
}

#' Simulate Coefficients \code{beta2} of \code{mu} for Joint Model
#'
#' This function simulates the coefficients \code{beta2} for a joint model based on its full conditional distribution, which incorporates the effects of covariates and random effects.
#'
#' @param mu A numeric vector representing the response variable for the second part of the model.
#' @param intercept2 A numeric value representing the intercept for the second part of the model.
#' @param beta A numeric value representing the coefficient for the first part of the model.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param W1 A numeric vector representing the latent random effects for the first part of the model.
#' @param W2 A numeric vector representing the latent random effects for the second part of the model.
#' @param Z2 A numeric matrix of covariates corresponding to the response variable \code{mu}.
#' @param hyper_fixed A list containing hyperparameters, specifically the prior precision for \code{beta2}.
#'
#' @return A numeric vector of simulated values for \code{beta2} from its full conditional distribution.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' mu <- rnorm(100)
#' intercept2 <- 0.5
#' beta <- 1
#' kappa_mu <- 2
#' W1 <- rnorm(100)
#' W2 <- rnorm(100)
#' Z2 <- matrix(rnorm(500), 100, 5)
#' hyper_fixed <- list(beta2 = 1)
#'
#' # Simulate beta2 for joint model
#' beta2_sim <- beta2_sim_JM(mu, intercept2, beta, kappa_mu, W1, W2, Z2, hyper_fixed)

beta2_sim_JM <- function(mu,
                         intercept2,
                         beta,
                         kappa_mu,
                         W1,
                         W2,
                         Z2,
                         hyper_fixed) {
  q <- ncol(Z2)
  latent.cov.inv <- hyper_fixed$beta2 * diag(1, q) + kappa_mu * (t(Z2) %*%
                                                                   Z2)
  latent.mean.part <- kappa_mu * (t(Z2) %*% (mu - intercept2 - beta * W1 - W2))
  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(q)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  proposals <- mm + vv
  return(proposals)
  
}
