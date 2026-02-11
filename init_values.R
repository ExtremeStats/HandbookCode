#' Initialize Parameters for Threshold Model
#'
#' This function generates initial values for the parameters of a threshold model based on the specified family.
#'
#' @param thr.family A character string specifying the family of the threshold model. Can be `"gamma"` or `"logNormal"`.
#' @param seed A numeric value used to set the seed for random number generation, ensuring reproducibility.
#'
#' @return A list containing:
#' \item{log.hyper.mu.init}{A numeric vector representing the initialized value of the model parameter.}
#' \item{log.hyper.mu.name}{An expression object with the name of the parameter (e.g., \code{expression(k)}).}
#'
#' @export
#'
#' @examplesIf FALSE
#' # The following is an example usage of the function and is not meant to be run directly.
#' thr.family <- "gamma"
#' seed <- 123
#' init_values <- init_fun_log.hyper.mu_thr_model(thr.family, seed)
#'
init_fun_log.hyper.mu_thr_model <- function(thr.family, seed) {
  set.seed(seed)
  if (thr.family == "gamma") {
    log.hyper.mu <- c(runif(1, 0, 1))
    log.hyper.mu.name <- c(expression(k))
  }
  if (thr.family == "logNormal") {
    log.hyper.mu <- c(runif(1, 0, 1))
    log.hyper.mu.name <- c(expression(k))
  }
  return(list(
    log.hyper.mu.init = log.hyper.mu,
    log.hyper.mu.name = log.hyper.mu.name
  ))
}


#' Initialize All Other Parameters for Threshold Model
#'
#' This function generates initial values for the parameters of a threshold model based on the specified family and additional covariates. It handles both `"gamma"` and `"logNormal"` families.
#'
#' @param Z2 A matrix of covariates for the model.
#' @param A A matrix (or data frame) representing additional covariates or transformations used in the model.
#' @param thr.family A character string specifying the family of the threshold model. Can be `"gamma"` or `"logNormal"`.
#' @param seed A numeric value used to set the seed for random number generation, ensuring reproducibility.
#' @param simulation A logical value indicating whether the function is being run in simulation mode.
#'
#' @return A list containing:
#' \item{init.all.other.param}{A numeric vector of initialized values for the model parameters, including kappas, intercepts, and covariate coefficients.}
#' \item{model.param.name.all.other.param}{A character vector or expression list representing the names of the parameters (e.g., \code{expression(kappa[w])}).}
#'
#' @export
#'
#' @examplesIf FALSE
#' # The following is an example usage of the function and is not meant to be run directly.
#' Z2 <- matrix(runif(10), nrow = 5, ncol = 2)
#' A <- matrix(runif(10), nrow = 5, ncol = 2)
#' thr.family <- "gamma"
#' seed <- 123
#' simulation <- TRUE
#' init_values <- init_fun_all_other_param_thr_model(Z2, A, thr.family, seed, simulation)
#'
init_fun_all_other_param_thr_model <- function(Z2, A, thr.family,  seed, simulation) {
  set.seed(seed)
  if (thr.family == "gamma") {
    kappa_w2_init <- runif(1, 1, 5)
    kappa_mu_init <- runif(1, 1, 5)
    intercept2.init <- runif(1, -1, 1)
    
    n1 <- nrow(Z2)
    
    beta2.init <- rep(runif(1, -1, 1), ncol(Z2))
    mu_init <- rep(runif(1, -0.1, 0.1), n1)
    w2_init <- rep(runif(1, -1, 1), times = n1)
    init1 <- c(kappa_w2_init,
               kappa_mu_init,
               intercept2.init,
               beta2.init,
               w2_init,
               mu_init)
    
    model.param.name <- c(
      expression(kappa[w]),
      expression(kappa[mu]),
      "intercept-e",
      paste0("beta2_", 1:ncol(Z2)),
      expression(mu[1]),
      expression(mu[2]),
      expression(W2[1]),
      expression(W2[2])
    )
  }
  if (thr.family == "logNormal") {
    kappa_w2_init <- runif(1, 1, 5)
    kappa_mu_init <- runif(1, 1, 5)
    intercept2.init <- runif(1, -1, 1)
    
    n1 <- nrow(Z2)
    
    beta2.init <- rep(runif(1, -1, 1), ncol(Z2))
    mu_init <- rep(runif(1, -0.1, 0.1), n1)
    w2_init <- rep(runif(1, -1, 1), times = n1)
    init1 <- c(kappa_w2_init,
               kappa_mu_init,
               intercept2.init,
               beta2.init,
               w2_init,
               mu_init)
    
    model.param.name <- c(
      expression(kappa[w]),
      expression(kappa[mu]),
      "intercept-e",
      paste0("beta2_", 1:ncol(Z2)),
      expression(mu[1]),
      expression(mu[2]),
      expression(W2[1]),
      expression(W2[2])
    )
  }
  return(
    list(
      init.all.other.param = init1,
      model.param.name.all.other.param = model.param.name
    )
  )
}

#' Initialize Parameters for Indicator Model with Probit Link
#'
#' This function generates initial values for the parameters of an indicator model, including coefficients for covariates and other model-specific parameters.
#'
#' @param Z2 A matrix of covariates for the model.
#' @param A A matrix (or data frame) representing additional covariates or transformations used in the model.
#' @param seed A numeric value used to set the seed for random number generation, ensuring reproducibility.
#'
#' @return A list containing:
#' \item{init.all.other.param}{A numeric vector of initialized values for model parameters, including kappas, intercepts, and covariate coefficients.}
#' \item{model.param.name.all.other.param}{A character vector or expression list representing the names of the parameters (e.g., \code{expression(kappa[w])}).}
#'
#' @export
#'
#' @examplesIf FALSE
#' # The following is an example usage of the function and is not meant to be run directly.
#' Z2 <- matrix(runif(10), nrow = 5, ncol = 2)
#' A <- matrix(runif(10), nrow = 5, ncol = 2)
#' seed <- 123
#' init_values <- init_fun_all_other_param_indicator_model(Z2, A, seed)
#'
init_fun_all_other_param_indicator_model_probit_link <- function(Z2, A, seed) {
  set.seed(seed)
  if(prob_link_ty)
  kappa_w2_init <- runif(1, 1, 5)
  kappa_mu_init <- runif(1, 1, 5)
  intercept2.init <- runif(1, -1, 1)
  
  n1 <- nrow(Z2)
  
  beta2.init <- rep(runif(1, -1, 1), ncol(Z2))
  mu_init <- rep(runif(1, -0.1, 0.1), n1)
  w2_init <- rep(runif(1, -1, 1), times = n1)
  init1 <- c(kappa_w2_init,
             kappa_mu_init,
             intercept2.init,
             beta2.init,
             w2_init,
             mu_init)
  
  model.param.name <- c(
    expression(kappa[w]),
    expression(kappa[mu]),
    "intercept-e",
    paste0("beta2_", 1:ncol(Z2)),
    expression(mu[1]),
    expression(mu[2]),
    expression(W2[1]),
    expression(W2[2])
  )
  
  return(
    list(
      init.all.other.param = init1,
      model.param.name.all.other.param = model.param.name
    )
  )
}



#' Initialize Parameters for Indicator Model with Logit Link
#'
#' This function generates initial values for the parameters of an indicator model, including coefficients for covariates and other model-specific parameters.
#'
#' @param Z2 A matrix of covariates for the model.
#' @param A A matrix (or data frame) representing additional covariates or transformations used in the model.
#' @param seed A numeric value used to set the seed for random number generation, ensuring reproducibility.
#'
#' @return A list containing:
#' \item{init.all.other.param}{A numeric vector of initialized values for model parameters, including kappas, intercepts, and covariate coefficients.}
#' \item{model.param.name.all.other.param}{A character vector or expression list representing the names of the parameters (e.g., \code{expression(kappa[w])}).}
#'
#' @export
#'
#' @examplesIf FALSE
#' # The following is an example usage of the function and is not meant to be run directly.
#' Z2 <- matrix(runif(10), nrow = 5, ncol = 2)
#' A <- matrix(runif(10), nrow = 5, ncol = 2)
#' seed <- 123
#' init_values <- init_fun_all_other_param_indicator_model(Z2, A, seed)
#'
init_fun_all_other_param_indicator_model_logit_link<-function(Z2, A, seed){
  set.seed(seed)
  kappa_w2_init<- runif(1, 1, 5)
  kappa_mu_init<- runif(1, 1, 5)
  intercept2.init<- runif(1, -1, 1)
  
  n1<- nrow(Z2)
  
  beta2.init<- rep(runif(1, -1, 1), ncol(Z2))
  mu_init<- rep(runif(1, -0.1, 0.1), n1) 
  w2_init<- rep(runif(1, -1, 1), times=n1)
  init1<-c( kappa_w2_init, kappa_mu_init, intercept2.init,
            beta2.init, w2_init, mu_init)
  
  model.param.name<- c(expression(kappa[w]), expression(kappa[mu]),
                       "intercept-e", paste0("beta2_",1:ncol(Z2)), 
                       expression(mu[1]), expression(mu[2]), 
                       expression(W2[1]), expression(W2[2]))
  
  return(list(init.all.other.param=init1, 
              model.param.name.all.other.param=model.param.name))
}









#' Initialize Parameters for Hyperparameter in the Bulk
#'
#' This function initializes the values and names for the hyperparameters based on the specified marginal distribution.
#'
#' @param mark_dist A character string indicating the type of marginal distribution. Possible values are:
#' `"eGPD"` for extended Generalized Pareto Distribution, `"bGPD"` for a mixture of beta-GPD, or `"tgGPD"` for a mixture of truncated gamma-GPD.
#' @param seed A numeric value to set the seed for random number generation to ensure reproducibility.
#'
#' @return A list containing:
#' \item{log.hyper.mu.init}{A named list of initial values for the hyperparameters, such as `k` and `xi`.}
#' \item{log.hyper.mu.name}{An expression vector with the names of the parameters for better display in plots or tables.}
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' mark_dist <- "eGPD"
#' seed <- 123
#' init_values <- init_fun_hyper.mu_JM(mark_dist, seed)
#'
init_fun_hyper.mu_JM <- function(mark_dist, seed) {
  set.seed(seed)
  if (mark_dist == "eGPD") {
    log.hyper.mu <- list("k" = runif(1, 0, 10), "xi" = runif(1, 0, 0.2))
    log.hyper.mu.name <- c(expression(k), expression(xi))
  } else if (mark_dist == "bGPD") {
    log.hyper.mu <- list("k" = runif(1, 0, 10))
    log.hyper.mu.name <- c(expression(k))
  } else if (mark_dist == "tgGPD") {
    log.hyper.mu <-  list("k" = runif(1, 0, 10))
    log.hyper.mu.name <- c(expression(k))
  }
  return(list(
    log.hyper.mu.init = log.hyper.mu,
    log.hyper.mu.name = log.hyper.mu.name
  ))
}


#' Initialize All Other Parameters for Joint Model
#'
#' This function generates initial values for the parameters of a joint model that are independent of the choice of marginal distributions. It handles initialization for fixed effects (FE) and other model types.
#'
#' @param Z1 A matrix of covariates for the first part of the model.
#' @param Z2 A matrix of covariates for the second part of the model.
#' @param A A matrix (or data frame) representing additional covariates or transformations used in the model.
#' @param Y A vector of response variables for the model.
#' @param seed A numeric value to set the seed for random number generation to ensure reproducibility.
#' @param simulation A logical value indicating whether the function is being run in simulation mode.
#' @param threshold A numeric value or vector specifying the threshold for handling certain parameters in the model.
#' @param mark_dist A character string indicating the type of marginal distribution. Accepts `"eGPD"` for extended Generalized Pareto Distribution, `"bGPD"` for a mixture of beta-GPD, or `"tgGPD"` for a mixture of truncated gamma-GPD.
#' @param model_type A character string indicating the type of model. Can be `"FE"` for fixed effects or other specified model types.
#' @param thr.acces.ind A logical vector indicating whether each observation exceeds the threshold.
#'
#' @return A list containing:
#' \item{init.all.other.param}{A named list of initialized values for the model parameters, including kappas, intercepts, coefficients for covariates, and other model-specific parameters.}
#' \item{model.param.name.all.other.param}{An expression vector or character vector representing the names of the parameters for better display in plots or tables.}
#'
#' @export
#'
#' @examplesIf FALSE
#' # Example usage (not meant to be run directly):
#' Z1 <- matrix(runif(10), nrow = 5, ncol = 2)
#' Z2 <- matrix(runif(10), nrow = 5, ncol = 2)
#' A <- matrix(runif(10), nrow = 5, ncol = 2)
#' Y <- rnorm(5)
#' seed <- 123
#' simulation <- TRUE
#' threshold <- c(0.5, 0.8)
#' mark_dist <- "eGPD"
#' model_type <- "FE"
#' thr.acces.ind <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
#' init_values <- init_fun_all_other_param_JM(Z1, Z2, A, Y, seed, simulation, threshold, mark_dist, model_type, thr.acces.ind)
#'
init_fun_hyper.GP_JM <- function(mark_dist, seed) {
  set.seed(seed)
  if (mark_dist == "eGPD") {
    log.hyper.mu <- NULL
    log.hyper.mu.name <- NULL
  } else if (mark_dist == "bGPD") {
    log.hyper.mu <- list("sigma.GP" = runif(1, 0, 10), "xi" = runif(1, 0, 0.2))
    log.hyper.mu.name <- c(expression(sigma[GP]), expression(xi))
  } else if (mark_dist == "tgGPD") {
    log.hyper.mu <-  list("sigma.GP" = runif(1, 0, 10), "xi" = runif(1, 0, 0.2))
    log.hyper.mu.name <- c(expression(sigma[GP]), expression(xi))
  }
  return(list(
    log.hyper.mu.init = log.hyper.mu,
    log.hyper.mu.name = log.hyper.mu.name
  ))
}

#' Initialize Values for Model Parameters in Joint Model
#'
#' Generates initial values for model parameters that are independent of the choice of mark distributions.
#'
#' @param Z1 A design matrix for covariates associated with the first set of parameters.
#' @param Z2 A design matrix for covariates associated with the second set of parameters.
#' @param A A matrix representing spatial relationships or adjacency information.
#' @param Y A vector or matrix of response variables (e.g., landslide counts or sizes).
#' @param seed An integer to set the random seed for reproducibility.
#' @param simulation A logical indicating if the function is being used in a simulation setting.
#' @param threshold A numeric value representing the threshold for the model.
#' @param thr.acces.ind A vector of indices indicating which data points exceed the threshold.
#'
#' @return A list containing initial values of parameters required for the model.
#' @export
#'
#' @examples
#' # Example usage:
#' Z1 <- matrix(rnorm(100), 20, 5)
#' Z2 <- matrix(rnorm(200), 20, 10)
#' A <- matrix(1, nrow = 20, ncol = 20)
#' Y <- rpois(20, lambda = 5)
#' seed <- 123
#' simulation <- TRUE
#' threshold <- 3
#' thr.acces.ind <- which(Y > threshold)
#' initialize_parameters(Z1, Z2, A, Y, seed, simulation, threshold, thr.acces.ind)
init_fun_all_other_param_JM <- function(Z1,
                                        Z2,
                                        A,
                                        Y,
                                        seed,
                                        simulation,
                                        threshold,
                                        mark_dist,
                                        model_type,
                                        thr.acces.ind) {
  #browser()
  set.seed(seed)
  kappa_w1_init <- runif(1, 1, 5)
  kappa_w2_init <- runif(1, 1, 5)
  kappa_eta_init <- runif(1, 1, 5)
  kappa_mu_init <- runif(1, 1, 5)
  intercept1.init <- runif(1, -1, 1)
  intercept2.init <- runif(1, -1, 1)
  beta.init <- runif(1, -1, 1)
  
  n.Q <- length(Y)
  
  beta1.init <- rep(runif(1, -1, 1), ncol(Z1))
  beta2.init <- rep(runif(1, -1, 1), ncol(Z2))
  eta_init <- rep(runif(1, -1, 1), times = length(Y))
  W1_init <- rep(runif(1, -1, 1), times = n.Q)
  if (mark_dist == "eGPD") {
    mu_init <- rep(runif(1, -0.1, 0.1), nrow(Z1))
  } else {
    mu_init <- ifelse(thr.acces.ind,
                      runif(1, -0.1, 0.1) + log(threshold),
                      runif(1, -0.1, 0.1))
  }
  w2_init <- rep(runif(1, -1, 1), times = n.Q)
  init <- list(
    "kappa_w1_init" = kappa_w1_init,
    "kappa_w2_init" = kappa_w2_init,
    "kappa_eta_init" = kappa_eta_init,
    "kappa_mu_init" = kappa_mu_init,
    "intercept1.init" = intercept1.init,
    "intercept2.init" = intercept2.init,
    "beta.init" = beta.init,
    "beta1.init" = beta1.init,
    "beta2.init" = beta2.init,
    "eta_init" = eta_init,
    "W1_init" = W1_init,
    "mu_init" = mu_init,
    "w2_init" = w2_init
  )
  if (model_type == "FE") {
    model.param.name <- c(
      expression(kappa[eta]),
      expression(kappa[mu]),
      "intercept1",
      "intercept2",
      paste0("beta1_", 1:ncol(Z1)),
      paste0("beta2_", 1:ncol(Z2)),
      expression(eta[1]),
      expression(eta[2]),
      expression(mu[1])
    )
  } else {
    model.param.name <- c(
      expression(kappa[w1]),
      expression(kappa[w2]),
      expression(kappa[eta]),
      expression(kappa[mu]),
      "intercept1",
      "intercept2",
      expression(beta),
      paste0("beta1_", 1:ncol(Z1)),
      paste0("beta2_", 1:ncol(Z2)),
      expression(eta[1]),
      expression(eta[2]),
      expression(W1[1]),
      expression(mu[1]),
      expression(W2[1])
    )
  }
  return(
    list(
      init.all.other.param = init,
      model.param.name.all.other.param = model.param.name
    )
  )
}
