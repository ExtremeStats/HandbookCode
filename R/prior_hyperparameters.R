#' Prior Distributions for Threshold Model's Distributions
#'
#' Computes the log prior density for the parameter of the distribution in a threshold model, supporting both gamma and log-normal families.
#'
#' @param cur_par A numeric vector representing the current parameter values of the mark distribution. For both gamma and log-normal families, it typically includes the shape parameter \code{k}.
#' @param thr.family A character string indicating the family of the threshold model. Can be either \code{"gamma"} or \code{"logNormal"}.
#' @param hyper.mu_fixed A numeric vector of length 2, providing the shape and rate hyperparameters of the gamma prior distribution for the mark distribution parameter.
#'
#' @return A numeric value representing the log-prior density of the parameter under the specified threshold family.
#'
#' @details
#' This function calculates the log-prior density for the shape parameter \code{k} of the mark distribution based on the specified threshold family:
#' - If \code{thr.family = "gamma"}, the prior for \code{k} is a gamma distribution with the specified hyperparameters.
#' - If \code{thr.family = "logNormal"}, the same gamma prior is applied to \code{k}.
#'
#' @export
#'
#' @examples
#' # Gamma threshold model
#' cur_par <- c(2) # Current parameter value (shape k)
#' thr.family <- "gamma"
#' hyper.mu_fixed <- c(1, 0.5) # Hyperparameters for the gamma prior
#' gamma_prior_markdist_thr(cur_par, thr.family, hyper.mu_fixed)
#'
#' # Log-normal threshold model
#' cur_par <- c(3) # Current parameter value (shape k)
#' thr.family <- "logNormal"
#' hyper.mu_fixed <- c(2, 1) # Hyperparameters for the gamma prior
#' gamma_prior_markdist_thr(cur_par, thr.family, hyper.mu_fixed)
gamma_prior_markdist_thr <- function(cur_par, thr.family, hyper.mu_fixed) {
  if (thr.family == "gamma") {
    k <- cur_par[1]
    logpriordens <- dgamma(k,
                           shape = hyper.mu_fixed[1],
                           rate = hyper.mu_fixed[2],
                           log = TRUE)
  }
  if (thr.family == "logNormal") {
    k <- cur_par[1]
    logpriordens <- dgamma(k,
                           shape = hyper.mu_fixed[1],
                           rate = hyper.mu_fixed[2],
                           log = TRUE)
  }
  
  return(logpriordens)
}





#' Prior Distributions for Joint Model's Landslide Size Distributions
#'
#' Computes the log-prior density for the parameters of the landslide size distribution in a joint model, supporting various distribution families like extended GPD (eGPD), beta-GPD (bGPD), and truncated gamma-GPD (tgGPD).
#'
#' @param cur_par A numeric vector representing the current parameter values of the mark distribution. For eGPD, includes both \code{k} and \code{xi}. For bGPD and tgGPD, only \code{k} is considered.
#' @param mark_dist A character string indicating the type of mark distribution. Possible values are \code{"eGPD"}, \code{"bGPD"}, and \code{"tgGPD"}.
#' @param hyper.mu_fixed A numeric vector of length 2, providing the shape and rate hyperparameters of the gamma prior distribution for the mark distribution parameter.
#'
#' @return A numeric value representing the log-prior density of the parameters under the specified mark distribution family.
#'
#' @details
#' This function calculates the log-prior density for the parameters of the mark distribution based on the specified family:
#' - **eGPD**: The prior for \code{k} is a gamma distribution, and the prior for \code{xi} is also a gamma distribution with shape 1 and rate 15.
#' - **bGPD**: Only the shape parameter \code{k} has a gamma prior distribution with the specified hyperparameters.
#' - **tgGPD**: Similar to bGPD, only \code{k} has a gamma prior distribution.
#'
#' @export
#'
#' @examples
#' # Prior for extended GPD (eGPD)
#' cur_par <- c(2, 0.5) # k and xi
#' mark_dist <- "eGPD"
#' hyper.mu_fixed <- c(1, 0.5)
#' gamma_prior_markdist_JM(cur_par, mark_dist, hyper.mu_fixed)
#'
#' # Prior for beta-GPD (bGPD)
#' cur_par <- c(3)
#' mark_dist <- "bGPD"
#' hyper.mu_fixed <- c(1, 0.5)
#' gamma_prior_markdist_JM(cur_par, mark_dist, hyper.mu_fixed)
#'
#' # Prior for truncated gamma-GPD (tgGPD)
#' cur_par <- c(4)
#' mark_dist <- "tgGPD"
#' hyper.mu_fixed <- c(2, 1)
#' gamma_prior_markdist_JM(cur_par, mark_dist, hyper.mu_fixed)
gamma_prior_markdist_JM <- function(cur_par, mark_dist, hyper.mu_fixed) {
  k <- cur_par[1]
  xi <- cur_par[2]
  if (mark_dist == "eGPD") {
    logpriordens <- dgamma(k,
                           shape = hyper.mu_fixed[1],
                           rate = hyper.mu_fixed[2],
                           log = TRUE) +
      dgamma(xi,
             shape = 1,
             rate = 15,
             log = TRUE)
  } else if (mark_dist == "bGPD") {
    logpriordens <- dgamma(k,
                           shape = hyper.mu_fixed[1],
                           rate = hyper.mu_fixed[2],
                           log = TRUE)
    
  } else if (mark_dist == "tgGPD") {
    logpriordens <- dgamma(k,
                           shape = hyper.mu_fixed[1],
                           rate = hyper.mu_fixed[2],
                           log = TRUE)
  }
  
  return(logpriordens)
}




#' Prior Distributions for GPD Parameters
#'
#' Computes the log-prior density for the parameters of a generalized Pareto distribution (GPD), specifically for the scale (\code{sigma}) and shape (\code{xi}) parameters.
#'
#' @param cur_par A numeric vector containing the current parameter values of the GPD: \code{cur_par[1]} is the scale parameter (\code{sigma}) and \code{cur_par[2]} is the shape parameter (\code{xi}).
#' @param hyper.mu_fixed A numeric vector of length 2 providing the shape and rate hyperparameters of the gamma prior distribution for the \code{sigma.GP} parameter.
#'
#' @return A numeric value representing the log-prior density of the GPD parameters.
#'
#' @details
#' The function computes the log-prior density for:
#' - **\code{sigma.GP}**: Assumed to have a gamma prior with shape and rate hyperparameters specified in \code{hyper.mu_fixed}.
#' - **\code{xi}**: Assumed to have a gamma prior with shape 1 and rate 15.
#'
#' @export
#'
#' @examples
#' # Example for calculating the log-prior density for GPD parameters
#' cur_par <- c(2, 0.5) # sigma.GP and xi
#' hyper.mu_fixed <- c(1, 0.5) # shape and rate hyperparameters for sigma.GP
#' gamma_prior_GPD_param(cur_par, hyper.mu_fixed)
gamma_prior_GPD_param <- function(cur_par, hyper.mu_fixed) {
  sigma.GP <- cur_par[1]
  xi <- cur_par[2]
  logpriordens <- dgamma(xi,
                         shape = 1,
                         rate = 15,
                         log = TRUE) +
    dgamma(sigma.GP,
           shape = hyper.mu_fixed[1],
           rate = hyper.mu_fixed[2],
           log = TRUE)
  return(logpriordens)
}






#' Lower and Upper Bounds for Landslide Size Distribution Parameters
#'
#' Provides the lower and upper bounds for the parameters of different mark distributions based on the specified distribution type.
#'
#' @param mark_dist A character string indicating the type of mark distribution. Supported values are:
#'   - `"eGPD"`: Extended generalized Pareto distribution.
#'   - `"bGPD"`: Mixture beta-GPD.
#'   - `"tgGPD"`: Mixture  truncated gamma-GPD.
#'
#' @return A list containing two elements:
#'   - \code{lb}: A numeric vector of the lower bounds for the parameters of the specified mark distribution.
#'   - \code{ub}: A numeric vector of the upper bounds for the parameters of the specified mark distribution.
#'
#' @details
#' The function returns parameter bounds depending on the type of mark distribution:
#'   - For `"eGPD"`, the lower bounds are \code{c(0, 0)}, and the upper bounds are \code{c(Inf, 1)}.
#'   - For `"bGPD"` and `"tgGPD"`, both have lower bounds of \code{c(0)} and upper bounds of \code{c(Inf)}.
#'
#' @export
#'
#' @examples
#' # Get parameter bounds for an extended GPD distribution
#' lb_ub_markdist("eGPD")
#'
#' # Get parameter bounds for a beta-GPD distribution
#' lb_ub_markdist("bGPD")
#'
#' # Get parameter bounds for a truncated GPD distribution
#' lb_ub_markdist("tgGPD")
lb_ub_markdist <- function(mark_dist) {
  if (mark_dist == "eGPD") {
    lb = c(0, 0)
    ub = c(Inf, 1)
  } else if (mark_dist == "bGPD") {
    lb = c(0)
    ub = c(Inf)
  } else if (mark_dist == "tgGPD") {
    lb = c(0)
    ub = c(Inf)
  }
  return(list("ub" = ub, "lb" = lb))
}





#' Lower and Upper Bounds for Generalized Pareto Distribution (GPD) Parameters
#'
#' Provides the lower and upper bounds for the parameters of certain Generalized Pareto Distributions based on the specified distribution type.
#'
#' @param mark_dist A character string indicating the type of GPD distribution. Supported values are:
#'   - `"bGPD"`: Beta-Generalized Pareto Distribution.
#'   - `"tgGPD"`: Truncated Generalized Pareto Distribution.
#'
#' @return A list containing two elements:
#'   - \code{lb}: A numeric vector of the lower bounds for the parameters of the specified GPD distribution.
#'   - \code{ub}: A numeric vector of the upper bounds for the parameters of the specified GPD distribution.
#'
#' @details
#' The function returns parameter bounds depending on the type of GPD distribution:
#'   - For both `"bGPD"` and `"tgGPD"`, the lower bounds are \code{c(0, 0)}, and the upper bounds are \code{c(Inf, 1)}.
#'
#' @export
#'
#' @examples
#' # Get parameter bounds for a beta-GPD distribution
#' lb_ub_GP("bGPD")
#'
#' # Get parameter bounds for a truncated GPD distribution
#' lb_ub_GP("tgGPD")
lb_ub_GP <- function(mark_dist) {
  if (mark_dist == "bGPD") {
    lb = c(0, 0)
    ub = c(Inf, 1)
  } else if (mark_dist == "tgGPD") {
    lb = c(0, 0)
    ub = c(Inf, 1)
  }
  return(list("ub" = ub, "lb" = lb))
}
