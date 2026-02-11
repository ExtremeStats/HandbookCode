#' Impute NA Values for Landslide Counts
#'
#' This function imputes missing (NA) values in landslide count data using a Poisson distribution.
#'
#' @param ind_NA_Y A logical vector indicating which elements in the count data are missing (NA).
#' @param eta A numeric vector of linear predictors for the Poisson distribution.
#' @param CV A character string specifying the cross-validation method to use. Accepts either `"OOS"` for out-of-sample
#' cross-validation or `"WS"` for within-sample cross-validation.
#'
#' @return A vector of imputed values, replacing NA elements in the original count data.
#' If \code{CV = "OOS"}, only the NA values are imputed. If \code{CV = "WS"}, all values are replaced with new counts.
#'
#' @export
#'
#' @examples
#' ind_NA_Y <- c(TRUE, FALSE, TRUE, FALSE)
#' eta <- c(1.2, 0.5, 1.8, 0.3)
#' CV <- "OOS"
#' imputed_values <- impute.NA.Y(ind_NA_Y, eta, CV)
#'
impute.NA.Y <- function(ind_NA_Y, eta, CV) {
  n1 <- length(eta)
  if (CV == "OOS") {
    imputed_NA_Y <- rpois(sum(ind_NA_Y), lambda = exp(eta)[ind_NA_Y])
  } else if (CV == "WS") {
    imputed_NA_Y <- rpois(n1, lambda = exp(eta))
  }
  return(imputed_NA_Y)
}


#' Impute NA Values for Landslide sizes
#'
#' This function imputes missing (NA) values in landslide sizes data based on the specified marginal distribution.
#' It handles different distribution families, including eGPD (extended Generalized Pareto Distribution), bGPD
#' (mixture of beta-GPD), and tgGPD (mixture of truncated gamma-GPD).
#'
#' @param CV A character string specifying the cross-validation type. Can be `"WS"` for within-sample or `"OOS"` for
#' out-of-sample.
#' @param ind_NA_A A logical vector indicating which elements in the accumulation data are missing (NA).
#' @param ind_zeros_counts A logical vector indicating which elements in the count data are zeros.
#' @param mu A numeric vector of linear predictors, typically on the log scale.
#' @param thr.prob A numeric value indicating the threshold exceedance probability.
#' @param cur_par A numeric vector of hyper parameters for the respective distribution (e.g., `k`, `xi`).
#' @param mark_dist A character string indicating the marginal distribution family. Accepts `"eGPD"`, `"bGPD"`, or `"tgGPD"`.
#' @param threshold A numeric vector representing the threshold for the generalized Pareto distribution.
#' @param thr.exceed.ind (Optional) A logical vector indicating if the threshold is exceeded. Default is `NULL`.
#'
#' @return A list containing:
#' \item{thr.exceed.ind}{A logical vector indicating whether the threshold is exceeded.}
#' \item{imputed_NA_A}{A numeric vector of imputed values for the missing accumulations.}
#'
#' @export
#'
#' @examplesIf FALSE
#' # The following is an example of the function usage and is not meant to be run directly.
#' CV <- "WS"
#' ind_NA_A <- c(TRUE, FALSE, TRUE)
#' ind_zeros_counts <- c(FALSE, TRUE, FALSE)
#' mu <- c(1.5, 0.8, 1.1)
#' thr.prob <- 0.9
#' cur_par <- c(0.5, 1.2, 0.7)
#' mark_dist <- "eGPD"
#' threshold <- c(0.6, 1.0, 0.8)
#' impute.NA.A(CV, ind_NA_A, ind_zeros_counts, mu, thr.prob, cur_par, mark_dist, threshold)
impute.NA.A <- function(CV,
                        ind_NA_A,
                        ind_zeros_counts,
                        mu,
                        thr.prob,
                        cur_par,
                        mark_dist,
                        threshold,
                        thr.exceed.ind = NULL) {
  #browser()
  n1 <- length(mu)
  if (mark_dist == "eGPD") {
    k <- cur_par[1]
    xi <- cur_par[2]
    if (CV == "WS") {
      mu <- mu #mu[!ind_zeros_counts]
      sigma <- exp(mu) / evd::qgpd(0.5 ^ (1 / k), scale = 1, shape = xi)
      imputed_NA_A <- rEGPD1(
        n = length(mu),
        k = k,
        xi = xi,
        sigma = sigma
      )
    } else{
      mu <- mu[ind_NA_A]
      sigma <- exp(mu) / evd::qgpd(0.5 ^ (1 / k), scale = 1, shape = xi)
      imputed_NA_A <- rEGPD1(
        n = length(mu),
        k = k,
        xi = xi,
        sigma = sigma
      )
    }
    
  } else if (mark_dist == "bGPD") {
    ### I assume parsimonious gamma distribution whenever we have missing observations
    k <- cur_par[1]
    sigma.GP <- cur_par[2]
    xi <- cur_par[3]
    
    if (CV == "WS") {
      mu <- mu #mu[!ind_zeros_counts]
      threshold <- threshold # threshold[!ind_zeros_counts]
      thr.exceed.ind <- exp(mu) > threshold
      sigma <- ifelse(thr.exceed.ind, sigma.GP, (threshold / exp(mu) - 1) * k)
      
      imputed_NA_A <-  ifelse(
        thr.exceed.ind,
        threshold + evd::rgpd(
          n = 1,
          loc = 0,
          scale = sigma.GP,
          shape = xi
        ),
        threshold * rbeta(
          n = 1,
          shape1 = k,
          shape2 = sigma
        )
      )
    } else {
      mu <- mu[ind_NA_A]
      threshold <- threshold[ind_NA_A]
      thr.exceed.ind <- exp(mu) > threshold
      sigma <- ifelse(thr.exceed.ind, sigma.GP, (threshold / exp(mu) - 1) * k)
      
      imputed_NA_A <-  ifelse(
        thr.exceed.ind,
        threshold + evd::rgpd(
          n = 1,
          loc = 0,
          scale = sigma.GP,
          shape = xi
        ),
        threshold * rbeta(
          n = 1,
          shape1 = k,
          shape2 = sigma
        )
      )
    }
    
  } else if (mark_dist == "tgGPD") {
    k <- cur_par[1]
    sigma.GP <- cur_par[2]
    xi <- cur_par[3]
    
    
    if (CV == "WS") {
      mu <- mu #mu[!ind_zeros_counts]
      threshold <- threshold #threshold[!ind_zeros_counts]
      thr.exceed.ind <- exp(mu) > threshold
      sigma <- ifelse(thr.exceed.ind, sigma.GP, exp(mu) / k)
      
      imputed_NA_A <-  ifelse(
        thr.exceed.ind,
        threshold + evd::rgpd(
          n = 1,
          loc = 0,
          scale = sigma.GP,
          shape = xi
        ),
        rgammat(
          n = 1,
          upper.bound  = threshold,
          shape = k,
          scale = sigma
        )
      )
      
    } else {
      mu <- mu[ind_NA_A]
      threshold <- threshold[ind_NA_A]
      thr.exceed.ind <- exp(mu) > threshold
      sigma <- ifelse(thr.exceed.ind, sigma.GP, exp(mu) / k)
      
      imputed_NA_A <- ifelse(
        thr.exceed.ind,
        threshold + evd::rgpd(
          n = 1,
          loc = 0,
          scale = sigma.GP,
          shape = xi
        ),
        rgammat(
          n = 1,
          upper.bound  = threshold,
          shape = k,
          scale = sigma
        )
      )
    }
  }
  return(list ("thr.exceed.ind" = thr.exceed.ind, "imputed_NA_A" = imputed_NA_A))
}
