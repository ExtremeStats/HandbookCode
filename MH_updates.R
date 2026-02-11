#' Update Latent Parameters Using  Random Walk Metropolis (RWM) in Parallel
#'
#' This function updates latent parameters in parallel using the Random Walk Metropolis (RWM) algorithm.
#' The updates are performed through a Gibbs sampling step with random walk proposals for each latent parameter.
#'
#' @param par_cur A vector of length \code{n} containing the current values of the latent parameters.
#' @param par_name A character string representing the name of the latent parameter vector being updated.
#' @param loglik_data A function representing the log-likelihood at the data level, depending on the latent parameters.
#' @param loglik_latent A function representing the log-likelihood at the latent level, often acting as a prior.
#' @param n Integer. The total number of latent parameters.
#' @param var_prop A vector of length \code{n}, representing the variance for the normal random walk proposals.
#' @param ... Additional arguments to be passed to the log-likelihood functions.
#'
#' @return A vector of updated latent parameter values, with length \code{n}.
#'
#' @examples
#' # Example usage of the latent_MH_update_parallel function:
#' n <- 20
#' par_cur <- rnorm(n)
#' par_name <- "latent_parameter"
#' var_prop <- rep(0.1, n) # Constant variance for simplicity
#' loglik_data <- function(latent_parameter, ...) -0.5 * (latent_parameter - 1)^2
#' loglik_latent <- function(latent_parameter, ...) dnorm(latent_parameter, mean = 1, sd = 0.1, log = TRUE)
#' updated_params <- latent_MH_update_parallel(par_cur, par_name, loglik_data, loglik_latent, n, var_prop)
#'
#' @seealso
#' Other functions involved in the Gibbs sampling process.
#'
#' @export

latent_MH_update_parallel <-
  function(par_cur,
           par_name,
           loglik_data,
           loglik_latent,
           ns,
           nt,
           var_prop,
           ...) {
    #browser()
    ellipsis <- list(...)
    args_loglik <- formalArgs(loglik_data)
    args_prior <- formalArgs(loglik_latent)
    stopifnot(is.finite(par_cur),
              isTRUE(all(args_loglik == formalArgs(loglik_data))),
              isTRUE(all(
                args_prior == formalArgs(loglik_latent)
              )),
              dim(var_prop) == c(ns, nt))
    
    ## current parameter value
    loglik_args <- ellipsis[names(ellipsis) %in% args_loglik]
    prior_args <- ellipsis[names(ellipsis) %in% args_prior]
    
    loglik_args[[par_name]] <- par_cur
    prior_args[[par_name]] <- par_cur
    
    loglik_curr <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    ## proposing the parameter
    par_prop <- matrix(rnorm(
      n = ns * nt,
      mean = par_cur,
      sd = sqrt(var_prop)
    ),
    nrow = ns,
    ncol = nt)
    
    loglik_args[[par_name]] <- par_prop
    prior_args[[par_name]] <- par_prop
    
    loglik_prop <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    
    log_MH_ratio <-
      loglik_prop - loglik_curr
    
    par_new <- ifelse((matrix(
      log(runif(ns * nt)), nrow = ns, ncol = nt
    ) < as.matrix(log_MH_ratio)) & (!is.na(as.matrix(log_MH_ratio))), par_prop, par_cur)
    
    
    # par_new<- ifelse((matrix(log(runif(ns * nt)), nrow = ns, ncol = nt) < as.matrix(log_MH_ratio)) ,
    #                    par_prop, par_cur)
    
    return(par_new)
  }


#' Update Latent Parameters with Constraints Using Random Walk Metropolis (RWM) in Parallel
#'
#' This function updates latent parameters using a constrained Random Walk Metropolis (RWM) approach in parallel.
#' The updates are performed using random walk proposals on the transformed parameter scale, applying constraints defined by the bounds.
#'
#' @param par_cur A vector of length \code{n} containing the current values of the latent parameters.
#' @param par_name A character string representing the name of the latent parameter vector being updated.
#' @param loglik_data A function representing the log-likelihood at the data level, depending on the latent parameters.
#' @param loglik_latent A function representing the log-likelihood at the latent level, often acting as a prior.
#' @param ns Integer. The number of spatial locations.
#' @param nt Integer. The number of time points. Here nt=1, as only one temporal replicates
#' @param var_prop A vector of length \code{n}, representing the variance for the normal random walk proposals.
#' @param transform A logical vector of length \code{n} indicating whether the parameters should be transformed to apply constraints.
#' @param lb A vector of length \code{n}, specifying the lower bound constraints for the parameters.
#' @param ub A vector of length \code{n}, specifying the upper bound constraints for the parameters.
#' @param ... Additional arguments to be passed to the log-likelihood functions.
#'
#' @return A vector of updated latent parameter values of length \code{n}.
#'
#' @examples
#' # Example usage of the latent_MH_update_parallel_constrain function:
#' ns <- 5
#' nt <- 4
#' n <- ns * nt
#' par_cur <- rnorm(n)
#' par_name <- "latent_parameter"
#' var_prop <- rep(0.1, n) # Constant variance for simplicity
#' lb <- rep(-1, n)
#' ub <- rep(1, n)
#' transform <- rep(TRUE, n)
#' loglik_data <- function(latent_parameter, ...) -0.5 * (latent_parameter - 1)^2
#' loglik_latent <- function(latent_parameter, ...) dnorm(latent_parameter, mean = 1, sd = 0.1, log = TRUE)
#' updated_params <- latent_MH_update_parallel_constrain(par_cur, par_name, loglik_data, loglik_latent,
#'                                                      ns, nt, var_prop, transform, lb, ub)
#'
#' @seealso
#' Other functions involved in Gibbs sampling with constrained parameters.
#'
#' @export
latent_MH_update_parallel_constrain <-
  function(par_cur,
           par_name,
           loglik_data,
           loglik_latent,
           ns,
           nt,
           var_prop,
           transform,
           lb,
           ub,
           ...) {
    #browser()
    ellipsis <- list(...)
    args_loglik <- formalArgs(loglik_data)
    args_prior <- formalArgs(loglik_latent)
    stopifnot(is.finite(par_cur),
              isTRUE(all(args_loglik == formalArgs(loglik_data))),
              isTRUE(all(
                args_prior == formalArgs(loglik_latent)
              )),
              dim(var_prop) == c(ns, nt))
    
    ## current parameter value
    loglik_args <- ellipsis[names(ellipsis) %in% args_loglik]
    prior_args <- ellipsis[names(ellipsis) %in% args_prior]
    
    loglik_args[[par_name]] <- par_cur
    prior_args[[par_name]] <- par_cur
    
    loglik_curr <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    
    # transforming the parameters to constrain scales
    trpar_curr <- ifelse(transform, transfo(par = par_cur, lb = lb, ub = ub), par_cur)
    # Sample a proposal from the normal: proposing the parameter
    trpar_prop <- matrix(rnorm(
      n = ns * nt,
      mean = trpar_curr,
      sd = sqrt(var_prop)
    ),
    nrow = ns,
    ncol = nt)
    # Compute the difference of log-Jacobin
    # adj <- ifelse(transform, jac_inv_transfo(tpar = trpar_prop, lb = lb, ub = ub, log=TRUE) -
    #                 jac_inv_transfo(tpar = trpar_curr, lb = lb, ub = ub, log=TRUE), 0)
    adj <- jac_inv_transfo(
      tpar = trpar_prop,
      lb = lb,
      ub = ub,
      log = TRUE
    ) -
      jac_inv_transfo(
        tpar = trpar_curr,
        lb = lb,
        ub = ub,
        log = TRUE
      )
    ### transforming back to original scale
    par_prop <- ifelse(transform,
                       inv_transfo(tpar = trpar_prop, lb = lb, ub = ub),
                       trpar_prop)
    
    loglik_args[[par_name]] <- par_prop
    prior_args[[par_name]] <- par_prop
    
    loglik_prop <-
      do.call(what = loglik_data, args = loglik_args) +
      do.call(what = loglik_latent, args = prior_args)
    
    log_MH_ratio <-
      loglik_prop - loglik_curr + adj
    
    par_new <- ifelse((matrix(
      log(runif(ns * nt)), nrow = ns, ncol = nt
    ) < as.matrix(log_MH_ratio)) & (!is.na(as.matrix(log_MH_ratio))), par_prop, par_cur)
    # par_new<- ifelse((matrix(log(runif(ns * nt)), nrow = ns, ncol = nt) < as.matrix(log_MH_ratio)) ,
    #                    par_prop, par_cur)
    
    return(par_new)
  }


#' Block Random Walk Metropolis Update
#'
#' This function performs block updates of parameters using random walk proposals within
#' the bounds of the parameter support. It optionally applies transformations to ensure the proposals
#' lie within the appropriate bounds.
#'
#' @param par_curr A vector of the current parameter values.
#' @param par_name A character string representing the name of the parameter being updated. The log-likelihood
#'   and log-prior functions should use this name as their argument for the parameter of interest.
#' @param loglik A function that computes the log-likelihood, depending on the parameter of interest.
#' @param logprior A function that computes the log-prior for the parameter of interest.
#' @param var_markdist A vector of the variances for the normal random walk proposals.
#' @param lb A vector specifying the lower bounds of the parameter support.
#' @param ub A vector specifying the upper bounds of the parameter support.
#' @param transform Logical; if \code{TRUE}, transformations are applied to convert the current parameter space
#'   to the entire real line before making Gaussian proposals. If \code{FALSE}, truncated normal proposals are used
#'   within the bounds \code{(lb, ub)}.
#' @param ... Additional arguments to be passed to the \code{loglik} and \code{logprior} functions.
#'
#' @return A vector of updated parameter values.
#'
#' @details
#' The function supports both constrained and unconstrained random walk proposals. When \code{transform = TRUE},
#' the parameter is transformed to the real line before proposing new values, and the acceptance ratio is adjusted
#' by the Jacobian of the transformation. When \code{transform = FALSE}, a truncated normal distribution is used
#' for proposing new values within the given bounds \code{(lb, ub)}.
#'
#' @examples
#' # Example usage of the block_rw_update function:
#' set.seed(123)
#' par_curr <- rnorm(5)
#' par_name <- "param"
#' loglik <- function(param) sum(dnorm(param, mean = 1, sd = 2, log = TRUE))
#' logprior <- function(param) sum(dnorm(param, mean = 0, sd = 1, log = TRUE))
#' var_markdist <- rep(0.1, length(par_curr))
#' lb <- rep(-2, length(par_curr))
#' ub <- rep(2, length(par_curr))
#' updated_param <- block_rw_update(par_curr, par_name, loglik, logprior, var_markdist, lb, ub, transform = TRUE)
#'
#' @seealso
#' Functions for constrained MCMC sampling.
#'
#' @export
block_rw_update <-
  function(par_curr,
           par_name,
           loglik,
           logprior,
           var_markdist,
           lb = -Inf,
           ub = Inf,
           transform = FALSE,
           ...) {
    #browser()
    ellipsis <- list(...)
    args_loglik <- formalArgs(loglik)
    args_logprior <- formalArgs(logprior)
    
    stopifnot(is.finite(par_curr))
    # if(lb == -Inf & ub == Inf){
    #   transform <- FALSE
    # }
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args_loglik]
    logprior_args <- ellipsis[names(ellipsis) %in% args_logprior]
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    logprior_args[[par_name]] <- par_curr
    
    logpost_curr <-
      do.call(what = loglik, args = loglik_args) +
      do.call(what = logprior, args = logprior_args)
    if (!transform) {
      par_prop <- rtnorm(
        n = length(par_curr),
        a = lb,
        b = ub,
        mean = par_curr,
        sd = sqrt(var_markdist)
      )
      adj <- sum(dtnorm(
        par_curr,
        a = lb,
        b = ub,
        mean = par_prop,
        sd = sqrt(var_markdist),
        log = TRUE
      )) -
        sum(dtnorm(
          par_prop,
          a = lb,
          b = ub,
          mean = par_curr,
          sd = sqrt(var_markdist),
          log = TRUE
        ))
    } else{
      # Transformation is TRUE
      # Y = f(X)
      # Transform parameter to the real line
      trpar_curr <- transfo(par = par_curr, lb = lb, ub = ub)
      # Sample a proposal from the normal
      trpar_prop <- rnorm(
        n = length(par_curr),
        mean = trpar_curr,
        sd = sqrt(var_markdist)
      )
      # Compute the difference of log-Jacobin
      adj <- sum(jac_inv_transfo(
        tpar = trpar_prop,
        lb = lb,
        ub = ub,
        log = TRUE
      )) -
        sum(jac_inv_transfo(
          tpar = trpar_curr,
          lb = lb,
          ub = ub,
          log = TRUE
        ))
      ### tranforming back to origial scale
      par_prop <- inv_transfo(tpar = trpar_prop, lb = lb, ub = ub)
    }
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    logprior_args[[par_name]] <- par_prop
    
    logpost_prop <-
      do.call(what = loglik, args = loglik_args) +
      do.call(what = logprior, args = logprior_args)
    log_MH_ratio <-
      logpost_prop - logpost_curr + adj ### because adjusted in the Jacobean
    
    if ((log_MH_ratio > log(runif(1))) & (!is.na(log_MH_ratio))) {
      par_new <- par_prop
    }
    
    # if((log_MH_ratio > log(runif(1)))){
    #   par_new <- par_prop
    # }
    return(par_new)
  }















