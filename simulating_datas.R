#' Simulate from a Model with \eqn{F_A} = **eGPD** Distribution
#'
#' This function simulates data from a model where the distribution of the response variable follows an eGPD (extended Generalized Pareto Distribution). The function takes as input various parameters including precision matrices, hyperparameters, covariate matrices, and model type to generate the simulated response and latent variables.
#'
#' @param Q Precision matrix of dimensions \eqn{n \times n}, where \eqn{n} is the total number of slope units.
#' @param other.hyper A list containing hyperparameters of the model in the following order:
#'   \itemize{
#'     \item \code{kappa_w1}: Hyperparameter for \eqn{W_1}.
#'     \item \code{kappa_w2}: Hyperparameter for \eqn{W_2}.
#'     \item \code{kappa_eta}: Hyperparameter for \eqn{\eta}.
#'     \item \code{kappa_mu}: Hyperparameter for \eqn{\mu}.
#'     \item \code{intercept1}: Intercept parameter for \eqn{\eta}.
#'     \item \code{intercept2}: Intercept parameter for \eqn{\mu}.
#'     \item \code{beta}: Coefficient for \eqn{W_1} in \eqn{\mu}.
#'   }
#' @param beta1 Coefficient associated with covariates \eqn{Z_1}.
#' @param beta2 Coefficient associated with covariates \eqn{Z_2}.
#' @param Z1 Matrix of covariates for the predictor \eqn{\eta}.
#' @param Z2 Matrix of covariates for the predictor \eqn{\mu}.
#' @param mark_dist The type of distribution for the response variable; e.g., \code{"eGPD"}.
#' @param model_type A character string indicating the model type. If \code{"FE"}, fixed effects are assumed for the random effects or else \code{"jSp"} for joint spatial model
#' @param hyper.mu hyperparameters (\eqn{\kappa}, \eqn{\xi}) of eGPD
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{Y}}{Simulated count variable following the specified \code{mark_dist}.}
#'   \item{\code{mu}}{Log-latent predictor of the size variable.}
#'   \item{\code{A}}{Simulated size variable from the eGPD distribution.}
#'   \item{\code{W1}}{Simulated random effect \eqn{W_1}.}
#'   \item{\code{eta}}{Log-latent predictor of the count variable.}
#'   \item{\code{W2}}{Simulated random effect \eqn{W_2}.}
#' }
#'
#' @export
#'
#' @examples
#' # Example usage of sim_mark_function
#' library(JM.landslide)
#' data("Wenchuan_info_used_for_simulation") ## Load Wenchuan landslides data
#'
#' # Load necessary libraries
#' library(spdep)
#' library(INLA)
#'
#' # Create adjacency graph file
#' nb2INLA("adjgraph-sim.txt", poly2nb(shp_selected, queen = FALSE, row.names = shp_selected$SU_ID))
#' adjacensy <- inla.read.graph(filename = "adjgraph-sim.txt")
#'
#' # Create the precision matrix
#' N <- adjacensy$n
#' diag.Q <- diag(adjacensy$nnbs, N)
#' A.Q <- matrix(0, nrow = N, ncol = N)
#' for (i in 1:N) {
#'   A.Q[i, adjacensy$nbs[[i]]] <- 1
#' }
#' Q <- diag.Q - A.Q  # Precision matrix
#'
#' # Set fixed parameters
#' kappa_w1 <- 10
#' kappa_w2 <- 5
#' kappa_eta <- 10
#' kappa_mu <- 5
#' intercept1 <- 2
#' intercept2 <- 4
#' beta <- 1
#' other.hyper <- list(
#'   kappa_w1 = kappa_w1,
#'   kappa_w2 = kappa_w2,
#'   kappa_eta = kappa_eta,
#'   kappa_mu = kappa_mu,
#'   intercept1 = intercept1,
#'   intercept2 = intercept2,
#'   beta = beta
#' )
#'
#' # Simulate covariates
#' Z1 <- mvtnorm::rmvnorm(N, mean = rep(0, 3))
#' Z2 <- mvtnorm::rmvnorm(N, mean = rep(0, 2))
#' hyper.mu<- c(20,0.1)
#'
#' set.seed(1)
#' beta1 <- runif(ncol(Z1), -1, 1)
#' beta2 <- runif(ncol(Z2), -1, 1)
#'
#' model_type <- "FE"
#'
#' # Simulate data using the function
#' set.seed(1)
#' Sim_mark_data <- sim_mark_function(
#'   Q = Q,
#'   hyper.mu =  hyper.mu,
#'   other.hyper = other.hyper,
#'   beta1 = beta1, beta2 = beta2,
#'   Z1 = as.matrix(Z1), Z2 = as.matrix(Z2),
#'   mark_dist = "eGPD",
#'   model_type = model_type
#' )
#'
#' # Access results
#' Y <- Sim_mark_data$Y
#' A <- Sim_mark_data$A
#' mu <- Sim_mark_data$mu
#' eta <- Sim_mark_data$eta
#' W1 <- Sim_mark_data$W1
#' W2 <- Sim_mark_data$W2

sim_mark_function <- function(Q,
                              hyper.mu,
                              other.hyper,
                              beta1,
                              beta2,
                              Z1,
                              Z2,
                              mark_dist,
                              model_type) {
  N <- dim(Q)[1]
  n1 <- n2 <- N
  p <- ncol(Z1)
  q <- ncol(Z2)
  hyper.mu <- hyper.mu
  kappa_w1 <- other.hyper$kappa_w1
  kappa_w2 <- other.hyper$kappa_w2
  kappa_eta <- other.hyper$kappa_eta
  kappa_mu <-  other.hyper$kappa_mu
  intercept1 <- other.hyper$intercept1
  intercept2 <- other.hyper$intercept2
  beta <- other.hyper$beta
  
  beta1 <- beta1
  beta2 <- beta2
  if (model_type == "FE") {
    W1_s <- rep(0, N)
    W2_s <- rep(0, N)
  } else{
    W1_s <- rMVNormP_eigen(n = 1,
                           mu = rep(0, N),
                           Sigma = kappa_w1 * Q)
    W2_s <- rMVNormP_eigen(n = 1,
                           mu = rep(0, N),
                           Sigma = kappa_w2 * Q)
  }
  eta <- intercept1 + Z1 %*% beta1 +  W1_s + rnorm(n1, mean = 0, sd = sqrt(1 /
                                                                             kappa_eta))
  Y = rpois(n = n1, lambda = exp(eta))
  
  mu <- intercept2 + Z2 %*% beta2 + beta * W1_s + W2_s + rnorm(n2, mean =
                                                                 0, sd = sqrt(1 / kappa_mu))
  k <- hyper.mu[1]
  xi <-  hyper.mu[2]
  sigma <- exp(mu) / evd::qgpd(0.5 ^ (1 / k), scale = 1, shape = xi)
  A <- rEGPD1(n = n2,
              k = k,
              xi = xi,
              sigma = sigma)
  
  return(list(
    Y = Y,
    mu = mu,
    A = A,
    W1 = W1_s,
    eta = eta,
    W2 = W2_s
  ))
}
