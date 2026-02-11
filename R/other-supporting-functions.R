#' CDF of Truncated Gamma Distribution
#'
#' This function calculates the cumulative distribution function (CDF) of a Gamma distribution
#' truncated to an upper bound. It returns either the log of the CDF or the CDF itself, depending
#' on the input arguments.
#'
#' @param x A numeric vector of quantiles.
#' @param upper.bound A numeric value indicating the upper bound of the truncation for the Gamma distribution.
#' @param shape A positive numeric value representing the shape parameter of the Gamma distribution.
#' @param scale A positive numeric value representing the scale parameter of the Gamma distribution.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as \code{log(p)}. Defaults to \code{TRUE}.
#'
#' @return A numeric vector of the CDF values of the truncated Gamma distribution at the specified quantiles \code{x}.
#'
#' @details
#' The function computes the CDF of a Gamma distribution truncated at the given upper bound.
#' If \code{log.p = TRUE}, the log of the CDF is returned, otherwise the CDF itself is returned.
#' The CDF is normalized by the probability that the Gamma variable falls within the truncated range.
#'
#' @examples
#' # Example usage of the pgammat function:
#' x <- seq(0, 5, by = 0.1)
#' shape <- 2
#' scale <- 1
#' upper.bound <- 5
#' cdf_values <- pgammat(x, upper.bound, shape, scale, log.p = FALSE)
#'
#' @seealso
#' \code{\link{pgamma}} for the Gamma distribution function without truncation.
#'
#' @export
pgammat <- function(x, upper.bound, shape, scale, log.p = TRUE) {
  F.a <- 0
  F.b <- pgamma(upper.bound, shape = shape, scale = scale)
  if (log.p) {
    cdf <- pgamma(x,
                  shape = shape,
                  scale = scale,
                  log.p = TRUE) - log(F.b - F.a)
  } else{
    cdf <- pgamma(x, shape = shape, scale = scale) / (F.b - F.a)
  }
  return(cdf)
}



#' Density of the Truncated Gamma Distribution
#'
#' This function calculates the density of a Gamma distribution truncated to an upper bound.
#' It provides the log-density or the density itself based on the user's input.
#'
#' @param x A numeric vector of quantiles at which the density is evaluated.
#' @param upper.bound A numeric value indicating the upper truncation bound for the Gamma distribution.
#' @param shape A positive numeric value representing the shape parameter of the Gamma distribution.
#' @param scale A positive numeric value representing the scale parameter of the Gamma distribution.
#' @param log Logical; if \code{TRUE}, the log-density is returned. If \code{FALSE}, the density is returned. Defaults to \code{TRUE}.
#'
#' @return A numeric vector of density values for the truncated Gamma distribution at the specified quantiles \code{x}.
#'
#' @details
#' The Gamma distribution is truncated at an upper bound specified by \code{upper.bound}.
#' The function normalizes the density by the probability that the Gamma variable lies within the range from zero to \code{upper.bound}.
#' If \code{log = TRUE}, the function returns the log-density; otherwise, it returns the density.
#'
#' @export
#'
#' @examples
#' # Calculate the density for a sequence of quantiles
#' x <- seq(0, 5, by = 0.1)
#' shape <- 2
#' scale <- 1
#' upper.bound <- 5
#' density_values <- dgammat(x, upper.bound, shape, scale, log = FALSE)
#' log_density_values <- dgammat(x, upper.bound, shape, scale, log = TRUE)
#'
#' @seealso
#' \code{\link{dgamma}} for the standard Gamma density function.
dgammat <- function(x, upper.bound, shape, scale, log = TRUE) {
  F.a <- 0
  F.b <- pgamma(upper.bound, shape = shape, scale = scale)
  if (log) {
    density <- dgamma(x,
                      shape = shape,
                      scale = scale,
                      log = TRUE) - log((F.b - F.a))
  } else{
    density <- dgamma(x, shape = shape, scale = scale) / (F.b - F.a)
  }
  return(density)
}


#' Simulation from Truncated Gamma Distribution
#'
#' This function generates random samples from a Gamma distribution that is truncated at a specified upper bound.
#'
#' @param n An integer specifying the number of observations to generate.
#' @param upper.bound A numeric value indicating the upper truncation bound for the Gamma distribution.
#' @param shape A positive numeric value representing the shape parameter of the Gamma distribution.
#' @param scale A positive numeric value representing the scale parameter of the Gamma distribution. Defaults to 1.
#'
#' @return A numeric vector of length \code{n} containing random samples drawn from the truncated Gamma distribution.
#'
#' @details
#' The Gamma distribution is truncated at an upper bound specified by \code{upper.bound}. The function generates samples from
#' the truncated distribution by first sampling from a uniform distribution over the cumulative distribution function (CDF)
#' range and then transforming these samples using the inverse Gamma CDF.
#'
#' @export
#'
#' @examples
#' # Generate 100 samples from a truncated Gamma distribution
#' n <- 100
#' shape <- 2
#' scale <- 1
#' upper.bound <- 5
#' samples <- rgammat(n, upper.bound, shape, scale)
#' hist(samples, main = "Histogram of Truncated Gamma Samples", xlab = "Value")
#'
#' @seealso
#' \code{\link{rgamma}} for sampling from the standard Gamma distribution.
rgammat <- function(n, upper.bound, shape, scale = 1) {
  F.a <- 0
  F.b <- pgamma(upper.bound, shape = shape, scale = scale)
  u <- runif(n, min = F.a, max = F.b)
  x <- qgamma(u, shape = shape, scale = scale)
  return(x)
}



#' Density Function of Scaled Beta Distribution
#'
#' This function computes the density of a Beta distribution scaled to have its support in the interval \eqn{(0, u.thr)}.
#'
#' @param x A numeric vector of quantiles for which the density is to be evaluated.
#' @param u.thr A positive numeric value specifying the upper bound of the support of the Beta distribution.
#' @param shape1 A positive numeric value representing the first shape parameter of the Beta distribution.
#' @param shape2 A positive numeric value representing the second shape parameter of the Beta distribution.
#' @param log Logical; if \code{TRUE}, the log-density is returned. Defaults to \code{TRUE}.
#'
#' @return A numeric vector of density values corresponding to the quantiles in \code{x}.
#'
#' @details
#' The standard Beta distribution has support on the interval \eqn{(0, 1)}. This function rescales the distribution so that
#' its support lies in \eqn{(0, u.thr)}. The density is computed by applying a change of variables to the standard Beta
#' density function.
#'
#' If \code{log} is \code{TRUE}, the log of the density is returned.
#'
#' @export
#'
#' @examples
#' # Evaluate the density of the scaled Beta distribution
#' x <- seq(0.01, 5, length.out = 100)
#' u.thr <- 5
#' shape1 <- 2
#' shape2 <- 3
#' density_values <- dbetat(x, u.thr, shape1, shape2, log = FALSE)
#' plot(x, density_values, type = "l", main = "Scaled Beta Density", ylab = "Density", xlab = "x")
#'
#' @seealso
#' \code{\link{dbeta}} for the standard Beta distribution density.
dbetat <- function(x, u.thr, shape1, shape2, log = TRUE) {
  # browser()
  if (log) {
    dens <- dbeta(x / u.thr,
                  shape1 = shape1,
                  shape2 = shape2,
                  log = TRUE) - log(u.thr)
  } else {
    dens <- dbeta(x / u.thr, shape1 = shape1, shape2 = shape2) / u.thr
  }
  return(dens)
}


#' Distribution Function of Scaled Beta Distribution
#'
#' This function computes the cumulative distribution function (CDF) of a Beta distribution scaled to have its support in the interval \eqn{(0, u.thr)}.
#'
#' @param x A numeric vector of quantiles for which the CDF is to be evaluated.
#' @param u.thr A positive numeric value specifying the upper bound of the support of the Beta distribution.
#' @param shape1 A positive numeric value representing the first shape parameter of the Beta distribution.
#' @param shape2 A positive numeric value representing the second shape parameter of the Beta distribution.
#' @param log.p Logical; if \code{TRUE}, the log of the CDF is returned. Defaults to \code{TRUE}.
#'
#' @return A numeric vector of CDF values corresponding to the quantiles in \code{x}.
#'
#' @details
#' The standard Beta distribution has support on the interval \eqn{(0, 1)}. This function rescales the distribution so that
#' its support lies in \eqn{(0, u.thr)}. The CDF is computed by applying a change of variables to the standard Beta
#' distribution function.
#'
#' If \code{log.p} is \code{TRUE}, the log of the CDF is returned.
#'
#' @export
#'
#' @examples
#' # Evaluate the CDF of the scaled Beta distribution
#' x <- seq(0.01, 5, length.out = 100)
#' u.thr <- 5
#' shape1 <- 2
#' shape2 <- 3
#' cdf_values <- pbetat(x, u.thr, shape1, shape2, log.p = FALSE)
#' plot(x, cdf_values, type = "l", main = "Scaled Beta CDF", ylab = "CDF", xlab = "x")
#'
#' @seealso
#' \code{\link{pbeta}} for the standard Beta distribution CDF.
pbetat <- function(x, u.thr, shape1, shape2, log.p = TRUE) {
  if (log.p) {
    dens <- pbeta(x / u.thr,
                  shape1 = shape1,
                  shape2 = shape2,
                  log.p = TRUE)
  } else {
    dens <- pbeta(x / u.thr, shape1 = shape1, shape2 = shape2)
  }
  return(dens)
}


#' Simulate from Scaled Beta Distribution
#'
#' This function generates random variates from a Beta distribution scaled to have its support in the interval \eqn{(0, u.thr)}.
#'
#' @param n An integer specifying the number of random variates to generate.
#' @param u.thr A positive numeric value specifying the upper bound of the support of the Beta distribution.
#' @param shape1 A positive numeric value representing the first shape parameter of the Beta distribution.
#' @param shape2 A positive numeric value representing the second shape parameter of the Beta distribution.
#'
#' @return A numeric vector of length \code{n} containing random variates from the scaled Beta distribution.
#'
#' @details
#' The standard Beta distribution has support on the interval \eqn{(0, 1)}. This function rescales the distribution
#' so that its support lies in \eqn{(0, u.thr)}. The random variates are generated by sampling from a standard Beta
#' distribution and then scaling the results by \code{u.thr}.
#'
#' @export
#'
#' @examples
#' # Generate 100 random variates from a scaled Beta distribution
#' set.seed(123)
#' n <- 100
#' u.thr <- 5
#' shape1 <- 2
#' shape2 <- 3
#' random_variates <- rbetat(n, u.thr, shape1, shape2)
#' hist(random_variates, main = "Scaled Beta Distribution", xlab = "Value", breaks = 20)
#'
#' @seealso
#' \code{\link{rbeta}} for the standard Beta distribution random variate generation.
rbetat <- function(n, u.thr, shape1, shape2) {
  x <- u.thr * rbeta(n = n,
                     shape1 = shape1,
                     shape2 = shape2)
  return(x)
}


#' Density Function of Extended GPD (Type 1)
#'
#' This function computes the density of the extended Generalized Pareto Distribution (GPD) of type 1.
#' The extended GPD modifies the standard GPD by raising its cumulative distribution function (CDF) to
#' the power \eqn{k}, providing a flexible model for tail behavior.
#'
#' @param x A numeric vector of quantiles at which to evaluate the density.
#' @param k A positive numeric value representing the extension parameter, which modifies the GPD.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A positive numeric value representing the scale parameter of the GPD.
#' @param log A logical value indicating whether to return the log-density (default is \code{TRUE}).
#'
#' @return A numeric vector of the density values of the extended GPD (Type 1) evaluated at \code{x}.
#'
#' @details
#' The extended GPD (Type 1) modifies the standard GPD by raising its cumulative distribution function (CDF)
#' to the power \eqn{k}. The density function is given by:
#' \deqn{f(x; k, \xi, \sigma) = k \cdot f_{\text{GPD}}(x; \xi, \sigma) \cdot F_{\text{GPD}}(x; \xi, \sigma)^{k-1}}
#' where \eqn{f_{\text{GPD}}} and \eqn{F_{\text{GPD}}} are the density and CDF of the standard GPD, respectively.
#'
#' @seealso
#' \code{\link[evd]{dgpd}} for the standard GPD density, and \code{\link[evd]{pgpd}} for the standard GPD cumulative distribution function.
#'
#' @export
#'
#' @examples
#' # Evaluate the density of the extended GPD (Type 1)
#' x <- seq(0, 10, length.out = 100)
#' k <- 2
#' xi <- 0.5
#' sigma <- 1
#' density <- dEGPD1(x, k, xi, sigma, log = FALSE)
#' plot(x, density, type = "l", main = "Extended GPD (Type 1) Density", xlab = "x", ylab = "Density")
#'
#' @importFrom evd dgpd pgpd
dEGPD1 <- function(x, k, xi, sigma, log = TRUE) {
  # browser()
  if (log) {
    dens <- log(k) + (evd::dgpd(
      x = x,
      loc = 0,
      scale = sigma,
      shape = xi,
      log = TRUE
    )) +
      ((k - 1) * log(
        evd::pgpd(
          q = x,
          loc = 0,
          scale = sigma,
          shape = xi,
          lower.tail = TRUE
        )
      ))
  } else{
    dens <- k * (evd::dgpd(
      x = x,
      loc = 0,
      scale = sigma,
      shape = xi
    )) * (evd::pgpd(
      q = x,
      loc = 0,
      scale = sigma,
      shape = xi,
      lower.tail = TRUE
    )) ^ (k - 1)
  }
  return(dens)
}


#' Cumulative Distribution Function of Extended GPD (Type 1)
#'
#' This function computes the cumulative distribution function (CDF) of the extended Generalized Pareto Distribution (GPD) of type 1.
#' The extended GPD modifies the standard GPD by raising its CDF to the power \eqn{k}, providing a flexible model for tail behavior.
#'
#' @param x A numeric vector of quantiles at which to evaluate the CDF.
#' @param k A positive numeric value representing the extension parameter, which modifies the GPD.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A positive numeric value representing the scale parameter of the GPD.
#' @param log A logical value indicating whether to return the log-CDF (default is \code{TRUE}).
#'
#' @return A numeric vector of the CDF values of the extended GPD (Type 1) evaluated at \code{x}.
#'
#' @details
#' The extended GPD (Type 1) modifies the standard GPD by raising its cumulative distribution function (CDF)
#' to the power \eqn{k}. The CDF is given by:
#' \deqn{F(x; k, \xi, \sigma) = F_{\text{GPD}}(x; \xi, \sigma)^{k}}
#' where \eqn{F_{\text{GPD}}} is the CDF of the standard GPD.
#'
#' @seealso
#' \code{\link[evd]{pgpd}} for the standard GPD cumulative distribution function.
#'
#' @export
#'
#' @examples
#' # Evaluate the CDF of the extended GPD (Type 1)
#' x <- seq(0, 10, length.out = 100)
#' k <- 2
#' xi <- 0.5
#' sigma <- 1
#' cdf <- pEGPD1(x, k, xi, sigma, log = FALSE)
#' plot(x, cdf, type = "l", main = "Extended GPD (Type 1) CDF", xlab = "x", ylab = "CDF")
#'
#' @importFrom evd pgpd
pEGPD1 <- function(x, k, xi, sigma, log = TRUE) {
  if (log) {
    cdf <- (k) * log(evd::pgpd(
      q = x,
      loc = 0,
      scale = sigma,
      shape = xi,
      lower.tail = TRUE
    ))
  } else{
    cdf <-  (evd::pgpd(
      q = x,
      loc = 0,
      scale = sigma,
      shape = xi,
      lower.tail = TRUE
    )) ^ (k)
  }
  return(cdf)
}


#' Random Number Generation from Extended GPD (Type 1)
#'
#' This function generates random variates from the extended Generalized Pareto Distribution (GPD) of type 1.
#' The extended GPD modifies the standard GPD by raising its CDF to the power \eqn{k}, allowing for more flexibility in tail behavior.
#'
#' @param n An integer specifying the number of random variates to generate.
#' @param k A positive numeric value representing the extension parameter that modifies the GPD.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A positive numeric value representing the scale parameter of the GPD.
#'
#' @return A numeric vector of length \code{n} containing random variates from the extended GPD (Type 1).
#'
#' @details
#' The extended GPD (Type 1) modifies the standard GPD by raising its cumulative distribution function (CDF)
#' to the power \eqn{k}. The random variates are generated using the inverse transform sampling method:
#' \deqn{X = \frac{\sigma}{\xi} \left( (1 - U^{1/k})^{-\xi} - 1 \right)}
#' where \eqn{U} is a uniform random variable on (0, 1).
#'
#' @seealso
#' \code{\link[evd]{rgpd}} for generating random variates from the standard GPD.
#'
#' @export
#'
#' @examples
#' # Generate random variates from the extended GPD (Type 1)
#' set.seed(123)
#' n <- 1000
#' k <- 2
#' xi <- 0.5
#' sigma <- 1
#' samples <- rEGPD1(n, k, xi, sigma)
#' hist(samples, breaks = 30, main = "Histogram of Extended GPD (Type 1) Samples", xlab = "Value")
rEGPD1 <- function(n, k, xi, sigma) {
  X <- (sigma / xi) * (((1 - (runif(
    n = n
  )) ^ (1 / k)) ^ (-xi)) - 1)
  return(X)
}




#' Density Function of Mixture of Truncated Gamma and GPD Distribution
#'
#' This function computes the density of a distribution that is a mixture of a truncated Gamma distribution and a Generalized Pareto Distribution (GPD).
#' The mixture model has a probability mass defined by \code{thr.prob} at the threshold \code{u.thr}, below which the distribution follows a truncated Gamma, and above which it follows a GPD.
#'
#' @param x A numeric vector of quantiles at which the density function is evaluated.
#' @param u.thr A numeric value specifying the threshold. Values of \code{x} less than \code{u.thr} follow the truncated Gamma distribution, and those greater than \code{u.thr} follow the GPD.
#' @param k A positive numeric value representing the shape parameter of the truncated Gamma distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector or value specifying the scale parameter of the truncated Gamma distribution.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param thr.prob A numeric value in (0, 1) specifying the probability mass of the GPD part in the mixture.
#' @param log Logical; if \code{TRUE}, the logarithm of the density is returned.
#'
#' @return A numeric vector of the same length as \code{x}, containing the (log) density values of the mixture distribution at the specified quantiles.
#'
#' @details
#' The mixture distribution is defined as follows:
#' - For \eqn{x < u.thr}, the density follows a truncated Gamma distribution with shape parameter \eqn{k}, scale parameter \eqn{\sigma}, and support in \eqn{(0, u.thr)}.
#' - For \eqn{x \geq u.thr}, the density follows a Generalized Pareto Distribution (GPD) with shape parameter \eqn{\xi}, scale parameter \eqn{\sigma.GP}, and location \eqn{u.thr}.
#' - The mixing probability \eqn{thr.prob} controls the weight of the GPD part in the overall distribution.
#'
#' @seealso
#' \code{\link[evd]{dgpd}} for the density of the standard GPD,
#' \code{\link{dgammat}} for the density of the truncated Gamma distribution.
#'
#' @export
#'
#' @examples
#' # Example usage of the dtgGPD function
#' set.seed(123)
#' x <- seq(0, 10, by = 0.1)
#' u.thr <- 5
#' k <- 2
#' xi <- 0.5
#' sigma <- 1
#' sigma.GP <- 0.8
#' thr.prob <- 0.3
#' dens <- dtgGPD(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log = FALSE)
#' plot(x, dens, type = "l", main = "Density of Mixture of Truncated Gamma and GPD", ylab = "Density", xlab = "x")
dtgGPD <- function(x,
                   u.thr,
                   k,
                   xi,
                   sigma,
                   sigma.GP,
                   thr.prob,
                   log = TRUE) {
  if (log) {
    dens <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    dens[ind] <- log(1 - thr.prob) + dgammat(
      x = x[ind],
      upper.bound = u.thr ,
      shape = k,
      scale = sigma[ind],
      log = TRUE
    )
    ind.e <- x > u.thr
    dens[ind.e] <- log(thr.prob) + evd::dgpd(
      x = x[ind.e],
      loc = u.thr,
      scale = sigma.GP,
      shape = xi,
      log = TRUE
    )
  } else{
    dens <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    dens[ind] <- log(1 - thr.prob) * dgammat(
      x = x[ind],
      upper.bound = u.thr,
      shape = k,
      scale = sigma[ind]
    )
    ind.e <- x > u.thr
    dens[ind.e] <- (thr.prob) * evd::dgpd(
      x = x[ind.e],
      loc = u.thr,
      scale = sigma.GP,
      shape = xi
    )
  }
  return(dens)
}


#' Distribution Function of Mixture of Truncated Gamma and GPD Distribution
#'
#' This function computes the cumulative distribution function (CDF) of a distribution that is a mixture of a truncated Gamma distribution and a Generalized Pareto Distribution (GPD).
#' The mixture model is characterized by a threshold \code{u.thr}, below which the CDF follows a truncated Gamma distribution and above which it follows a GPD.
#'
#' @param x A numeric vector of quantiles at which the CDF is evaluated.
#' @param u.thr A numeric value specifying the threshold. Values of \code{x} less than \code{u.thr} follow the truncated Gamma distribution, and those greater than \code{u.thr} follow the GPD.
#' @param k A positive numeric value representing the shape parameter of the truncated Gamma distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector or value specifying the scale parameter of the truncated Gamma distribution.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param thr.prob A numeric value in (0, 1) specifying the probability mass of the GPD part in the mixture.
#' @param log.p Logical; if \code{TRUE}, the logarithm of the CDF is returned.
#'
#' @return A numeric vector of the same length as \code{x}, containing the (log) CDF values of the mixture distribution at the specified quantiles.
#'
#' @details
#' The mixture distribution is defined as follows:
#' - For \eqn{x < u.thr}, the CDF follows a truncated Gamma distribution with shape parameter \eqn{k}, scale parameter \eqn{\sigma}, and support in \eqn{(0, u.thr)}.
#' - For \eqn{x \geq u.thr}, the CDF is a mixture of a truncated Gamma and a Generalized Pareto Distribution (GPD), controlled by the mixing probability \eqn{thr.prob}.
#'
#' @seealso
#' \code{\link[evd]{pgpd}} for the CDF of the standard GPD,
#' \code{\link{pgammat}} for the CDF of the truncated Gamma distribution.
#'
#' @export
#'
#' @examples
#' # Example usage of the ptgGPD function
#' set.seed(123)
#' x <- seq(0, 10, by = 0.1)
#' u.thr <- 5
#' k <- 2
#' xi <- 0.5
#' sigma <- 1
#' sigma.GP <- 0.8
#' thr.prob <- 0.3
#' cdf_values <- ptgGPD(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log.p = FALSE)
#' plot(x, cdf_values, type = "l", main = "CDF of Mixture of Truncated Gamma and GPD", ylab = "CDF", xlab = "x")
ptgGPD <- function(x,
                   u.thr,
                   k,
                   xi,
                   sigma,
                   sigma.GP,
                   thr.prob,
                   log.p = TRUE) {
  if (log.p) {
    cdf <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    #cdf[ind]<- log(1 - thr.prob) + pgammat(x = x[ind], range = c(0, u.thr), shape = k, scale = sigma[ind], log.p=TRUE)
    cdf[ind] <- log(1 - thr.prob) + pgammat(
      x = x[ind],
      upper.bound = u.thr,
      shape = k,
      scale = sigma[ind],
      log.p = TRUE
    )
    ind.e <- x > u.thr
    cdf[ind.e] <- log((1 - thr.prob) + thr.prob  * evd::pgpd(
      x = x[ind.e],
      loc = u.thr,
      shape = xi,
      scale = sigma.GP
    ))
  } else{
    cdf <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    cdf[ind] <- (1 - thr.prob) * pgammat(
      x = x[ind],
      upper.bound = u.thr,
      shape = k,
      scale = sigma[ind]
    )
    ind.e <- x > u.thr
    cdf[ind.e] <- (1 - thr.prob)  + (thr.prob) * evd::pgpd(
      x = x[ind.e],
      loc = u.thr,
      scale = sigma.GP,
      shape = xi
    )
  }
  return(cdf)
}


#' Conditional Cumulative Distribution Function of Mixture of Truncated Gamma and GPD
#'
#' This function calculates the conditional cumulative distribution function (CDF) for a mixture of truncated Gamma and Generalized Pareto Distributions (GPD), given the parameter vector \code{mu}.
#' The CDF is computed based on whether the input value \code{x} is greater than or less than the threshold \code{u.thr}.
#'
#' @param x A numeric value specifying the quantile at which the CDF is evaluated.
#' @param mu A numeric vector representing the location parameter values for the conditional CDF.
#' @param u.thr A numeric value specifying the threshold. Values of \code{x} greater than \code{u.thr} are treated as GPD, while those less than \code{u.thr} follow the truncated Gamma distribution.
#' @param k A positive numeric value representing the shape parameter of the truncated Gamma distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector specifying the scale parameter of the truncated Gamma distribution.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param thr.prob A numeric value in (0, 1) specifying the probability mass of the GPD component in the mixture.
#'
#' @return A numeric vector of the same length as \code{mu}, containing the conditional CDF values of the mixture distribution at the specified quantile \code{x}.
#'
#' @details
#' The mixture distribution is modeled conditionally based on the threshold \code{u.thr}. Specifically:
#' - If \code{x > u.thr}, the conditional CDF is computed using the GPD with a certain probability mass, while the truncated Gamma distribution contributes when \code{exp(mu)} is less than \code{u.thr}.
#' - If \code{x <= u.thr}, the conditional CDF is calculated using only the truncated Gamma distribution.
#'
#' @seealso
#' \code{\link[evd]{pgpd}} for the CDF of the standard GPD,
#' \code{\link{pgammat}} for the CDF of the truncated Gamma distribution.
#'
#' @export
#'
#' @examples
#' # Example usage of cond.ptgGPD function
#' set.seed(123)
#' x <- 5
#' mu <- rnorm(10)
#' u.thr <- 4
#' k <- 2
#' xi <- 0.5
#' sigma <- rep(1, length(mu))
#' sigma.GP <- 0.8
#' thr.prob <- 0.3
#' cdf_values <- cond.ptgGPD(x, mu, u.thr, k, xi, sigma, sigma.GP, thr.prob)
#' print(cdf_values)
cond.ptgGPD <- function(x,
                        mu,
                        u.thr,
                        k,
                        xi,
                        sigma,
                        sigma.GP,
                        thr.prob) {
  #browser()
  cdf <- rep(NA, length(mu))
  if (x > u.thr) {
    ind <- exp(as.numeric(mu)) < u.thr
    cdf[ind] <- 1
    cdf[!ind] <-  1 - thr.prob + thr.prob * evd::pgpd(
      x = x,
      loc = u.thr,
      scale = sigma.GP,
      shape = xi
    )
  } else{
    ind <- exp(as.numeric(mu)) > u.thr
    cdf[ind] <- 0
    cdf[!ind] <-  (1 - thr.prob) * pgammat(
      x = x,
      upper.bound = u.thr,
      shape = k,
      scale = sigma[!ind]
    )
  }
  return(cdf)
  
}


#' Simulate Random Numbers from a Mixture of Truncated Gamma and GPD Distribution
#'
#' This function generates random samples from a mixture of a truncated Gamma and a Generalized Pareto Distribution (GPD). The samples are drawn based on an indicator vector, where elements specify which distribution (Gamma or GPD) to sample from.
#'
#' @param n An integer specifying the total number of random samples to generate.
#' @param u.thr A numeric value specifying the threshold that separates the Gamma and GPD components.
#' @param k A positive numeric value representing the shape parameter of the truncated Gamma distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector specifying the scale parameter for the truncated Gamma distribution.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param ind A logical vector of length \code{n} indicating which samples should be generated from the GPD (\code{TRUE}) and which should be generated from the truncated Gamma distribution (\code{FALSE}).
#'
#' @return A numeric vector of length \code{n}, containing the random samples drawn from the mixture of truncated Gamma and GPD distributions.
#'
#' @details
#' The function generates \code{n} random samples from a mixture distribution where:
#' - Samples corresponding to \code{ind = TRUE} are drawn from the GPD with threshold \code{u.thr}, scale parameter \code{sigma.GP}, and shape parameter \code{xi}.
#' - Samples corresponding to \code{ind = FALSE} are drawn from a truncated Gamma distribution with shape \code{k}, scale \code{sigma}, and truncated at \code{u.thr}.
#'
#' @seealso
#' \code{\link[evd]{rgpd}} for sampling from the Generalized Pareto Distribution,
#' \code{\link{rgammat}} for sampling from the truncated Gamma distribution.
#'
#' @export
#'
#' @examples
#' # Example usage of rtgGPD function
#' set.seed(123)
#' n <- 100
#' u.thr <- 5
#' k <- 2
#' xi <- 0.3
#' sigma <- 1
#' sigma.GP <- 0.8
#' ind <- runif(n) > 0.7  # 70% Gamma, 30% GPD
#' samples <- rtgGPD(n, u.thr, k, xi, sigma, sigma.GP, ind)
#' hist(samples, breaks = 30, main = "Samples from Mixture of Truncated Gamma and GPD")
rtgGPD <- function(n, u.thr, k, xi, sigma, sigma.GP, ind) {
  n1 <- sum(ind)
  x <- ifelse(
    ind,
    evd::rgpd(
      n = n1,
      loc = u.thr,
      scale = sigma.GP,
      shape = xi
    ),
    rgammat(
      n = n - n1,
      upper.bound = u.thr,
      shape = k,
      scale = sigma
    )
  )
  return(x)
}


#' Density Function of Mixture of Truncated Beta and GPD
#'
#' This function computes the density of a mixture distribution consisting of a truncated beta distribution on the interval \code{(0, u.thr)} and a Generalized Pareto Distribution (GPD) for values exceeding the threshold \code{u.thr}.
#'
#' @param x A numeric vector of observations at which the density needs to be evaluated.
#' @param u.thr A positive numeric value specifying the threshold that separates the truncated beta and GPD components.
#' @param k A positive numeric value representing the shape1 parameter of the truncated beta distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector specifying the shape2 parameters of the truncated beta distribution; must be of the same length as \code{x}.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param thr.prob A numeric value between 0 and 1 indicating the mixture probability for the GPD component (i.e., the probability that a sample exceeds the threshold \code{u.thr}).
#' @param log A logical value; if \code{TRUE}, the log density is returned. If \code{FALSE}, the density is returned.
#'
#' @return A numeric vector of the same length as \code{x}, containing the density (or log density) values of the mixture distribution evaluated at each element of \code{x}.
#'
#' @details
#' The function computes the density of a mixture distribution where:
#' - For values of \code{x} in the interval \code{(0, u.thr)}, the density is derived from a truncated beta distribution with parameters \code{k} and \code{sigma}.
#' - For values of \code{x} greater than \code{u.thr}, the density is derived from a GPD with location \code{u.thr}, scale \code{sigma.GP}, and shape \code{xi}.
#' The mixture is determined by \code{thr.prob}, which specifies the probability of being in the GPD component.
#'
#' @seealso
#' \code{\link{dbetat}} for the truncated beta distribution density function,
#' \code{\link[evd]{dgpd}} for the Generalized Pareto Distribution density function.
#'
#' @export
#'
#' @examples
#' # Example usage of dbGPD function
#' set.seed(123)
#' x <- seq(0, 10, length.out = 100)
#' u.thr <- 5
#' k <- 2
#' xi <- 0.3
#' sigma <- rep(1, length(x))
#' sigma.GP <- 0.8
#' thr.prob <- 0.3
#' dens <- dbGPD(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log = FALSE)
#' plot(x, dens, type = "l", main = "Density of Mixture of Truncated Beta and GPD")
dbGPD <- function(x,
                  u.thr,
                  k,
                  xi,
                  sigma,
                  sigma.GP,
                  thr.prob,
                  log = TRUE) {
  #browser()
  if (log) {
    dens <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    dens[ind] <- log(1 - thr.prob) + dbetat(
      x = x[ind],
      u.thr = u.thr,
      shape1 = k,
      shape2 = sigma[ind],
      log = TRUE
    )
    ind.e <- x > u.thr
    dens[ind.e] <- log(thr.prob) + evd::dgpd(
      x = x[ind.e],
      loc = u.thr,
      scale = sigma.GP,
      shape = xi,
      log = TRUE
    )
  } else{
    dens <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    dens[ind] <- (1 - thr.prob) * dbetat(
      x = x[ind],
      u.thr = u.thr,
      shape1 = k,
      shape2 = sigma[ind]
    )
    ind.e <- x > u.thr
    dens[ind.e] <- (thr.prob) * evd::dgpd(
      x = x[ind.e],
      loc = u.thr,
      scale = sigma.GP,
      shape = xi
    )
  }
  return(dens)
  
}

#' Cumulative Distribution Function of Mixture of Truncated Beta and GPD
#'
#' This function computes the cumulative distribution function (CDF) of a mixture distribution that consists of a truncated beta distribution on the interval \code{(0, u.thr)} and a Generalized Pareto Distribution (GPD) for values exceeding the threshold \code{u.thr}.
#'
#' @param x A numeric vector of observations at which the CDF needs to be evaluated.
#' @param u.thr A positive numeric value specifying the threshold that separates the truncated beta and GPD components.
#' @param k A positive numeric value representing the shape1 parameter of the truncated beta distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector specifying the shape2 parameters of the truncated beta distribution; must be of the same length as \code{x}.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param thr.prob A numeric value between 0 and 1 indicating the mixture probability for the GPD component (i.e., the probability that a sample exceeds the threshold \code{u.thr}).
#' @param log.p A logical value; if \code{TRUE}, the log of the CDF is returned. If \code{FALSE}, the CDF is returned.
#'
#' @return A numeric vector of the same length as \code{x}, containing the CDF (or log CDF) values of the mixture distribution evaluated at each element of \code{x}.
#'
#' @details
#' The function computes the CDF of a mixture distribution where:
#' - For values of \code{x} in the interval \code{(0, u.thr)}, the CDF is derived from a truncated beta distribution with parameters \code{k} and \code{sigma}.
#' - For values of \code{x} greater than \code{u.thr}, the CDF is derived from a GPD with location \code{u.thr}, scale \code{sigma.GP}, and shape \code{xi}.
#' The mixture is determined by \code{thr.prob}, which specifies the probability of being in the GPD component.
#'
#' @seealso
#' \code{\link{pbetat}} for the truncated beta distribution CDF function,
#' \code{\link[evd]{pgpd}} for the Generalized Pareto Distribution CDF function.
#'
#' @export
#'
#' @examples
#' # Example usage of pbGPD function
#' set.seed(123)
#' x <- seq(0, 10, length.out = 100)
#' u.thr <- 5
#' k <- 2
#' xi <- 0.3
#' sigma <- rep(1, length(x))
#' sigma.GP <- 0.8
#' thr.prob <- 0.3
#' cdf <- pbGPD(x, u.thr, k, xi, sigma, sigma.GP, thr.prob, log.p = FALSE)
#' plot(x, cdf, type = "l", main = "CDF of Mixture of Truncated Beta and GPD")
pbGPD <- function(x,
                  u.thr,
                  k,
                  xi,
                  sigma,
                  sigma.GP,
                  thr.prob,
                  log.p = TRUE) {
  if (log.p) {
    cdf <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    cdf[ind] <- log(1 - thr.prob) + pbetat(
      x = x[ind],
      u.thr = u.thr,
      shape1 = k,
      shape2 = sigma[ind],
      log.p = TRUE
    )
    ind.e <- x > u.thr
    cdf[ind.e] <- log((1 - thr.prob) + thr.prob * evd::pgpd(
      x = x[ind.e],
      loc = u.thr,
      shape = xi,
      scale = sigma.GP,
      lower.tail = TRUE
    )
    )
  } else{
    cdf <- rep(0, length(x))
    ind <- (x < u.thr) & (x > 0)
    cdf[ind] <- (1 - thr.prob) * dbetat(
      x = x[ind],
      u.thr = u.thr,
      shape1 = k,
      shape2 = sigma[ind]
    )
    ind.e <- x > u.thr
    cdf[ind.e] <-  (1 - thr.prob) + (thr.prob) * evd::pgpd(
      x = x[ind.e],
      loc = u.thr,
      shape = xi,
      scale = sigma.GP
    )
  }
  return(cdf)
  
}


#' Conditional Cumulative Distribution Function of Mixture of Truncated Beta and GPD
#'
#' This function computes the conditional cumulative distribution function (CDF) for a mixture distribution, which includes a truncated beta distribution for values below the threshold \code{u.thr} and a Generalized Pareto Distribution (GPD) for values above \code{u.thr}.
#'
#' @param x A numeric value representing the high quantile at which the CDF is evaluated. It should be greater than \code{u.thr}.
#' @param mu A numeric vector representing the location parameter of the mixture distribution.
#' @param u.thr A positive numeric value specifying the threshold separating the truncated beta and GPD components.
#' @param k A positive numeric value representing the shape1 parameter of the truncated beta distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector specifying the scale parameter for the truncated beta distribution; must be of the same length as \code{mu}.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param thr.prob A numeric value between 0 and 1 indicating the mixture probability for the GPD component (i.e., the probability that a sample exceeds the threshold \code{u.thr}).
#' @param log.p A logical value; if \code{TRUE}, the log of the CDF is returned. If \code{FALSE}, the CDF is returned.
#'
#' @return A numeric vector of the same length as \code{mu}, containing the conditional CDF values (or log CDF values) evaluated at \code{x}.
#'
#' @details
#' This function provides a conditional CDF for a mixture distribution:
#' - For \code{x} greater than \code{u.thr}, the CDF is derived from the GPD component.
#' - For \code{x} less than or equal to \code{u.thr}, the CDF is derived from the truncated beta distribution.
#' The value of \code{mu} is used to conditionally determine which distribution is applicable, and \code{thr.prob} specifies the probability of being in the GPD component.
#'
#' @seealso
#' \code{\link{pbetat}} for the truncated beta distribution CDF function,
#' \code{\link[evd]{pgpd}} for the Generalized Pareto Distribution CDF function.
#'
#' @export
#'
#' @examples
#' # Example usage of cond.pbGPD function
#' x <- 8
#' mu <- seq(log(1), log(20), length.out = 100)
#' u.thr <- 5
#' k <- 2
#' xi <- 0.5
#' sigma <- rep(1, length(mu))
#' sigma.GP <- 1.5
#' thr.prob <- 0.4
#' cdf <- cond.pbGPD(x, mu, u.thr, k, xi, sigma, sigma.GP, thr.prob)
#' plot(mu, cdf, type = "l", main = "Conditional CDF of Mixture Distribution")
cond.pbGPD <- function(x,
                       mu,
                       u.thr,
                       k,
                       xi,
                       sigma,
                       sigma.GP,
                       thr.prob) {
  #browser()
  cdf <- rep(NA, length(mu))
  if (x > u.thr) {
    ind <- exp(as.numeric(mu)) < u.thr
    cdf[ind] <- 1
    cdf[!ind] <-  1 - thr.prob + thr.prob * evd::pgpd(
      q = x,
      loc = u.thr,
      shape = xi,
      scale = sigma.GP ,
      lower.tail = TRUE
    )
  } else{
    ind <- exp(as.numeric(mu)) > u.thr
    cdf[ind] <- 0
    cdf[!ind] <-  (1 - thr.prob) * pbetat(
      x = x,
      u.thr = u.thr,
      shape1 = k,
      shape2 = sigma[!ind],
      log.p = TRUE
    )
  }
  
  
  cdf <- ifelse(
    x > u.thr,
    ifelse(
      exp(as.numeric(mu)) < u.thr,
      1,
      1 - thr.prob + thr.prob * evd::pgpd(
        q = x,
        loc = u.thr,
        shape = xi,
        scale = sigma.GP ,
        lower.tail = TRUE
      )
    ),
    ifelse(
      exp(as.numeric(mu)) > u.thr,
      0,
      (1 - thr.prob) * pbetat(
        x = x,
        u.thr = u.thr,
        shape1 = k,
        shape2 = sigma[!ind],
        log.p = TRUE
      )
    )
  )
  
  return(cdf)
}


#' Simulate from Mixture of Beta and GPD Distributions
#'
#' This function simulates values from a mixture distribution, composed of a beta distribution truncated to (0, \code{u.thr}) and a Generalized Pareto Distribution (GPD) for values greater than \code{u.thr}.
#'
#' @param n An integer specifying the number of observations to generate.
#' @param u.thr A positive numeric value indicating the threshold separating the truncated beta and GPD components.
#' @param k A positive numeric value representing the shape1 parameter of the truncated beta distribution.
#' @param xi A numeric value representing the shape parameter of the GPD.
#' @param sigma A numeric vector of scale parameters for the truncated beta distribution, with length equal to \code{n}.
#' @param sigma.GP A positive numeric value representing the scale parameter of the GPD.
#' @param ind A logical vector of length \code{n}, indicating which observations should be generated from the GPD component (\code{TRUE}) and which from the truncated beta component (\code{FALSE}).
#'
#' @return A numeric vector of length \code{n} containing simulated values from the mixture distribution.
#'
#' @details
#' This function generates values from a mixture of two distributions:
#' - The truncated beta distribution is used for values below the threshold \code{u.thr}.
#' - The GPD is used for values greater than \code{u.thr}.
#' The argument \code{ind} is used to determine which distribution to sample from for each observation.
#'
#' @seealso
#' \code{\link{rbetat}} for simulating from the truncated beta distribution,
#' \code{\link[evd]{rgpd}} for simulating from the Generalized Pareto Distribution.
#'
#' @export
#'
#' @examples
#' # Example usage of rbGPD function
#' n <- 100
#' u.thr <- 5
#' k <- 2
#' xi <- 0.5
#' sigma <- rep(1, n)
#' sigma.GP <- 1.5
#' ind <- runif(n) < 0.3  # 30% of the samples will be from GPD
#' simulated_values <- rbGPD(n, u.thr, k, xi, sigma, sigma.GP, ind)
#' hist(simulated_values, breaks = 30, main = "Simulated Mixture of Beta and GPD")
rbGPD <- function(n, u.thr, k, xi, sigma, sigma.GP, ind) {
  n1 <- sum(ind)
  x <- ifelse(
    ind,
    evd::rgpd(
      n = n1,
      loc = u.thr,
      scale = sigma.GP,
      shape = xi
    ),
    rbetat(
      n = n - n1,
      u.thr = u.thr,
      shape1 = k,
      shape2 = sigma
    )
  )
  return(x)
}


#' Simulate from Univariate Truncated Normal Distribution
#'
#' This function generates random samples from a univariate normal distribution truncated to a specified interval. The truncation bounds, mean, and standard deviation of the distribution are provided as parameters. The implementation is suitable for cases where \code{a < 37} or \code{b > -37}.
#'
#' @param n An integer specifying the number of random samples to generate. Default is 1.
#' @param mean A numeric scalar representing the location parameter (mean) of the normal distribution. Default is 0.
#' @param sd A positive numeric scalar representing the scale parameter (standard deviation) of the normal distribution. Default is 1.
#' @param a A numeric value specifying the lower bound of truncation. Default is \code{-Inf}.
#' @param b A numeric value specifying the upper bound of truncation. Default is \code{Inf}.
#'
#' @return A numeric vector of length \code{n} containing random samples from the truncated normal distribution.
#'
#' @details
#' If the standardized bounds are extreme (\code{b_std < -37} or \code{a_std > 37}), the function issues a warning and uses the \code{TruncatedNormal} package to sample from the distribution.
#'
#' For \code{a_std < 0}, the function computes the truncated distribution's cumulative distribution function (CDF) values at \code{a} and \code{b}, then generates samples within the interval using the inverse CDF (quantile function). For \code{a_std >= 0}, it uses a complementary CDF transformation.
#'
#' @seealso
#' \code{\link[TruncatedNormal]{rtnorm}} for more advanced truncated normal sampling, especially for extreme tail events.
#'
#' @export
#'
#' @examples
#' # Simulate 1000 samples from a truncated normal distribution
#' set.seed(123)
#' samples <- rtnorm(n = 1000, mean = 0, sd = 1, a = -2, b = 2)
#' hist(samples, breaks = 30, main = "Histogram of Truncated Normal Samples")
#'
#' # Simulate with extreme bounds, triggering a fallback to the TruncatedNormal package
#' set.seed(123)
#' extreme_samples <- rtnorm(n = 1000, mean = 0, sd = 1, a = 38, b = Inf)
#' hist(extreme_samples, breaks = 30, main = "Histogram of Extreme Truncated Normal Samples")
rtnorm <- function(n = 1,
                   mean = 0,
                   sd = 1,
                   a = -Inf,
                   b = Inf) {
  stopifnot(
    length(a) == 1L,
    length(b) == 1L,
    length(mean) == 1L,
    length(sd) == 1L,
    isTRUE(a < b),
    isTRUE(sd > 0)
  )
  a_std <- (a - mean) / sd
  b_std <- (b - mean) / sd
  if (b_std < -37 | a_std > 37) {
    warning(
      "Interval requested is beyond numerical tolerance.\nUse \"TruncatedNormal\" package \"rtnorm\" instead for rare events simulation."
    )
    return(TruncatedNormal::rtnorm(
      n = n,
      mu = mean,
      sd = sd,
      lb = a,
      ub = b
    ))
  }
  if (a_std < 0) {
    Fa <- pnorm(a_std)
    Fb <- pnorm(b_std)
    mean + sd * qnorm(Fa + runif(n) * (Fb - Fa))
  } else {
    Fa <- pnorm(a_std, lower.tail = FALSE)
    Fb <- pnorm(b_std, lower.tail = FALSE)
    mean + sd * (-qnorm(Fa - (Fa - Fb) * runif(n)))
  }
}

#' Density Function of a Univariate Truncated Normal Distribution
#'
#' This function calculates the density (or log density) of a univariate normal distribution that has been truncated to a specified interval. The density is adjusted for the truncation bounds, and the mean and standard deviation of the distribution are provided as parameters.
#'
#' @inheritParams rtnorm
#' @param x A numeric vector of observations at which the density is to be evaluated.
#' @param log A logical value indicating whether to return the log density. Default is \code{FALSE}.
#'
#' @return A numeric vector of density (or log density) values for the truncated normal distribution evaluated at \code{x}.
#'
#' @details
#' The function utilizes the \code{TruncatedNormal} package to efficiently compute the normalizing constant, ensuring accurate density values even for intervals with extreme truncation.
#'
#' If both \code{a} and \code{b} are finite, the function computes the truncated normal density adjusted by the log of the normalizing constant. In cases where one bound is infinite, it defaults to the standard normal distribution behavior.
#'
#' @seealso
#' \code{\link{rtnorm}} for simulating from a truncated normal distribution.
#'
#' @export
#'
#' @examples
#' # Evaluate the density of a truncated normal distribution
#' set.seed(123)
#' x <- seq(-2, 2, length.out = 100)
#' density <- dtnorm(x, mean = 0, sd = 1, a = -1, b = 1)
#' plot(x, density, type = "l", main = "Density of Truncated Normal Distribution")
#'
#' # Evaluate the log density of the truncated normal distribution
#' log_density <- dtnorm(x, mean = 0, sd = 1, a = -1, b = 1, log = TRUE)
#' plot(x, log_density, type = "l", main = "Log Density of Truncated Normal Distribution")
dtnorm <- function(x,
                   mean = 0,
                   sd = 1,
                   a = -Inf,
                   b = Inf,
                   log = FALSE) {
  stopifnot(
    length(a) == 1L,
    length(b) == 1L,
    length(mean) == 1L,
    length(sd) == 1L,
    isTRUE(a < b),
    isTRUE(sd > 0)
  )
  dens <- dnorm(x, mean = mean, sd = sd, log = TRUE)
  dens <- dens  - TruncatedNormal::lnNpr(a = (a - mean) / sd, b = (b - mean) /
                                           sd)
  return(dens)
}

#' Transform Parameter to an Unconstrained Scale
#'
#' This function transforms a parameter constrained by specified lower and upper bounds to an unconstrained scale. The transformation depends on whether the bounds are finite or infinite.
#'
#' @param par A scalar parameter to be transformed.
#' @param lb A scalar specifying the lower bound of the parameter.
#' @param ub A scalar specifying the upper bound of the parameter.
#'
#' @return The transformed parameter on an unconstrained scale.
#'
#' @details
#' The transformation is determined by the bounds provided:
#' - If both \code{lb} and \code{ub} are infinite, the parameter remains unchanged.
#' - If only \code{lb} is finite, a log transformation is applied: \code{log(par - lb)}.
#' - If only \code{ub} is finite, a log transformation is applied: \code{log(ub - par)}.
#' - If both \code{lb} and \code{ub} are finite, a logistic transformation is applied: \code{qlogis((par - lb) / (ub - lb))}.
#'
#' @seealso
#' Section 56 of the Stan Reference Manual, version 2.9, provides more details on transformations for unconstrained parameter spaces. \url{https://github.com/stan-dev/stan/releases/download/v2.9.0/stan-reference-2.9.0.pdf}
#'
#' @export
#'
#' @examples
#' # Transform parameters with different bounds
#' transfo(par = 0.5, lb = -Inf, ub = Inf)  # No transformation
#' transfo(par = 3, lb = 1, ub = Inf)       # Log transformation with lower bound
#' transfo(par = 7, lb = -Inf, ub = 10)     # Log transformation with upper bound
#' transfo(par = 2, lb = 1, ub = 5)         # Logistic transformation between bounds
#'
#' # Vectorized transformation
#' params <- c(2, 3, 4)
#' transfo(par = params, lb = 1, ub = 5)
transfo <- Vectorize(function(par, lb, ub) {
  stopifnot(length(par) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if (lb == -Inf & ub == Inf) {
    return(par)
  } else if (lb > -Inf & ub == Inf) {
    return(log(par - lb))
  } else if (lb == -Inf & ub < Inf) {
    return(log(ub - par))
  } else if (lb > -Inf & ub < Inf) {
    return(qlogis((par - lb) / (ub - lb)))
  }
}, vectorize.args = c("par", "lb", "ub"))



#' Jacobian of the Transformations from \code{transfo}
#'
#' This function computes the Jacobian (or its logarithm) of the transformation applied to a parameter by \code{transfo}. The Jacobian is useful for reversing the transformation and ensuring the correct change of variable during inference.
#'
#' @param tpar A scalar representing the transformed parameter.
#' @param lb A scalar specifying the lower bound of the parameter.
#' @param ub A scalar specifying the upper bound of the parameter.
#' @param log A logical value indicating whether to return the log-Jacobian. If \code{TRUE}, the log-Jacobian is returned. If \code{FALSE}, the Jacobian is returned.
#'
#' @return The (log-)Jacobian of the transformation. If \code{log = TRUE}, returns the log-Jacobian. Otherwise, returns the Jacobian.
#'
#' @details
#' The function calculates the Jacobian based on the type of bounds provided for the parameter:
#' - If both \code{lb} and \code{ub} are infinite, the Jacobian is 1 (log-Jacobian is 0).
#' - If only \code{lb} is finite, the Jacobian is \code{exp(tpar)}.
#' - If only \code{ub} is finite, the Jacobian is \code{-exp(tpar)}.
#' - If both \code{lb} and \code{ub} are finite, a logistic transformation is applied, and the log-Jacobian includes the terms \code{log(ub - lb)} and the log of the sigmoid function.
#'
#' @export
#'
#' @examples
#' # Compute the log-Jacobian of various transformations
#' jac_inv_transfo(tpar = 0.5, lb = -Inf, ub = Inf, log = TRUE)  # No transformation
#' jac_inv_transfo(tpar = 1, lb = 0, ub = Inf, log = TRUE)       # Log-transformation with lower bound
#' jac_inv_transfo(tpar = -1, lb = -Inf, ub = 5, log = TRUE)     # Log-transformation with upper bound
#' jac_inv_transfo(tpar = 0.2, lb = 1, ub = 5, log = TRUE)       # Logistic transformation between bounds
#'
#' # Compute the Jacobian (not in log-scale)
#' jac_inv_transfo(tpar = 0.5, lb = -Inf, ub = Inf, log = FALSE)
jac_inv_transfo <- Vectorize(function(tpar, lb, ub, log = FALSE) {
  if (lb == -Inf & ub == Inf) {
    ljac <- 0  ## log(tpar) should be zero, why there is log(tpar) we need to check
  }
  if (lb > -Inf & ub == Inf) {
    ljac <- tpar
  } else if (lb == -Inf & ub < Inf) {
    ljac <- tpar
  } else if (lb > -Inf & ub < Inf) {
    ljac <- log(ub - lb) + plogis(tpar, log.p = TRUE) + plogis(tpar, log.p = TRUE, lower.tail = FALSE)
  }
  if (log) {
    return(ljac)
  } else{
    return(exp(ljac))
  }
}, vectorize.args = c("tpar", "lb", "ub"))


#' Gradient of the Log-Jacobian of Transformations
#'
#' This function computes the gradient of the log-Jacobian of the parameter transformation applied by \code{transfo}. The gradient is essential for adjusting parameter transformations during inference, particularly when using gradient-based methods.
#'
#' @param tpar A scalar representing the transformed parameter.
#' @param lb A scalar specifying the lower bound of the parameter.
#' @param ub A scalar specifying the upper bound of the parameter.
#'
#' @return The gradient of the log-Jacobian of the transformation, evaluated at \code{tpar}.
#'
#' @details
#' The function calculates the gradient based on the type of bounds provided for the parameter:
#' - If both \code{lb} and \code{ub} are infinite, the gradient is 0, as no transformation is applied.
#' - If only \code{lb} is finite, the gradient is 1, indicating a shift to unconstrained space.
#' - If only \code{ub} is finite, the gradient is also 1.
#' - If both \code{lb} and \code{ub} are finite, the gradient is computed as \code{-1 + 2 * plogis(-tpar)}, which accounts for the logistic transformation applied to map the parameter between the bounds.
#'
#' @export
#'
#' @examples
#' # Compute the gradient of the log-Jacobian of various transformations
#' dlogjac_inv_transfo(tpar = 0.5, lb = -Inf, ub = Inf)   # No transformation
#' dlogjac_inv_transfo(tpar = 1, lb = 0, ub = Inf)        # Lower-bound transformation
#' dlogjac_inv_transfo(tpar = -1, lb = -Inf, ub = 5)      # Upper-bound transformation
#' dlogjac_inv_transfo(tpar = 0.2, lb = 1, ub = 5)        # Logistic transformation between bounds
dlogjac_inv_transfo <- Vectorize(function(tpar, lb, ub) {
  if (lb == -Inf & ub == Inf) {
    return(0)
  }
  if (lb > -Inf & ub == Inf) {
    return(1)
  } else if (lb == -Inf & ub < Inf) {
    return(1)
  } else{
    -1 + 2 * plogis(-tpar)
  }
}, vectorize.args = c("tpar", "lb", "ub"))


#' Inverse Transformation to Original Parameter Scale
#'
#' This function reverses the transformations applied to parameters, converting them back to their original scale. It is used in conjunction with parameter transformations that constrain parameters within specific bounds.
#'
#' @param tpar A scalar representing the transformed parameter, typically on the unconstrained scale.
#' @param lb A scalar specifying the lower bound of the original parameter's support.
#' @param ub A scalar specifying the upper bound of the original parameter's support.
#'
#' @return The parameter transformed back to its original scale.
#'
#' @details
#' The function handles different types of transformations based on the bounds provided:
#' - If both \code{lb} and \code{ub} are infinite, the transformation is the identity, and \code{tpar} is returned unchanged.
#' - If only \code{lb} is finite, the transformation is \code{exp(tpar) + lb}, which shifts and scales \code{tpar} to enforce the lower bound.
#' - If only \code{ub} is finite, the transformation is \code{ub - exp(tpar)}, shifting and scaling \code{tpar} to enforce the upper bound.
#' - If both \code{lb} and \code{ub} are finite, the transformation uses the logistic function: \code{lb + (ub - lb) * plogis(tpar)}. This maps \code{tpar} smoothly between \code{lb} and \code{ub}.
#'
#' @export
#'
#' @examples
#' # Inverse transformations for various bounds
#' inv_transfo(tpar = 1, lb = -Inf, ub = Inf)   # No transformation
#' inv_transfo(tpar = 0.5, lb = 0, ub = Inf)    # Lower-bound transformation
#' inv_transfo(tpar = -1, lb = -Inf, ub = 5)    # Upper-bound transformation
#' inv_transfo(tpar = 0.2, lb = 1, ub = 5)      # Transformation between finite bounds

inv_transfo <- Vectorize(function(tpar, lb, ub) {
  stopifnot(length(tpar) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if (lb == -Inf & ub == Inf) {
    return(tpar)
  } else if (lb > -Inf & ub == Inf) {
    return(exp(tpar) + lb)
  } else if (lb == -Inf & ub < Inf) {
    return(ub - exp(tpar))
  } else{
    return(lb + (ub - lb) * plogis(tpar))
  }
}, vectorize.args = c("tpar", "lb", "ub"))



##' Truncated Eigen Square Root (or Inverse Square Root)
#'
#' Computes the square root or inverse square root of a given matrix, designed for use with covariance or precision matrices. The function handles small negative eigenvalues due to numerical errors by setting them to zero and throws an error if any significant negative eigenvalues are detected.
#'
#' @param Sigma A p x p matrix for which to compute the truncated square root. Typically, this will be a covariance or precision matrix.
#' @param inv A logical value indicating whether to compute the inverse square root (\code{TRUE}) or the square root (\code{FALSE}) of the matrix.
#'
#' @return A p x k matrix (where k is the rank of \code{Sigma}) representing the truncated square root or inverse square root of \code{Sigma}.
#'
#' @details
#' This function performs the following steps:
#' - Computes the eigen decomposition of \code{Sigma}.
#' - Small eigenvalues within a tolerance (\code{1e-12}) are set to zero to handle potential numerical errors.
#' - Throws an error if any large negative eigenvalues are present, indicating that something may be wrong with \code{Sigma}.
#' - Computes the square root or inverse square root of \code{Sigma} based on the eigenvalues and returns the truncated matrix.
#'
#' @note
#' The function is particularly useful for covariance or precision matrices, where all eigenvalues should be non-negative. Any large negative eigenvalues could indicate issues with the input matrix.
#'
#' @export
#'
#' @examples
#' # Example covariance matrix
#' Sigma <- matrix(c(4, 2, 2, 3), nrow = 2)
#' # Compute square root
#' trunc_eigen_sqrt(Sigma, inv = FALSE)
#' # Compute inverse square root
#' trunc_eigen_sqrt(Sigma, inv = TRUE)
trunc_eigen_sqrt <- function(Sigma, inv) {
  es <- eigen(Sigma)
  
  # Small negative eigenvalues can occur just due to numerical error and should
  # be set back to zero.
  es$values[(es$values < 1e-12) & (es$values > -1e-12)] <- 0
  
  # If any eigen values are large negative throw error (something wrong)
  if (any(es$values < -1e-12))
    stop("Non-trivial negative eigenvalues present")
  
  # calculate square root and reveal rank (k)
  k <- sum(es$values > 0)
  if (!inv) {
    L <- es$vectors %*% diag(sqrt(es$values))
  } else if (inv) {
    L <- es$vectors %*% diag(1 / sqrt(es$values))
  }
  return(L[, 1:k, drop = F])
}

#' Covariance Parameterization using Eigen Decomposition
#'
#' Generates samples from a multivariate normal distribution with a specified mean vector and covariance matrix, utilizing the eigen decomposition for parameterization.
#'
#' @param n An integer specifying the number of samples to draw.
#' @param mu A numeric vector of length \code{p} representing the mean of the distribution.
#' @param Sigma A \code{p} x \code{p} covariance matrix.
#'
#' @return A \code{p} x \code{n} matrix of samples drawn from the multivariate normal distribution with mean vector \code{mu} and covariance matrix \code{Sigma}.
#'
#' @details
#' This function uses the eigen decomposition of the covariance matrix \code{Sigma} to generate samples from the corresponding multivariate normal distribution. The process is as follows:
#' - Perform the truncated eigen decomposition of \code{Sigma}.
#' - Draw samples from a standard normal distribution.
#' - Transform these samples to have the specified mean \code{mu} and covariance structure \code{Sigma}.
#'
#' @export
#'
#' @examples
#' # Example parameters
#' mu <- c(0, 0) # Mean vector
#' Sigma <- matrix(c(2, 0.5, 0.5, 1), nrow = 2) # Covariance matrix
#' # Generate 100 samples
#' samples <- rMVNormC_eigen(n = 100, mu = mu, Sigma = Sigma)
rMVNormC_eigen <- function(n, mu, Sigma) {
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv = FALSE)
  k <- ncol(L)
  Z <- matrix(rnorm(k * n), k, n)
  X <- L %*% Z
  X <- sweep(X, 1, mu, FUN = `+`)
}

#' Precision Parameterization using Eigen Decomposition
#'
#' Generates samples from a multivariate normal distribution with a specified mean vector and precision matrix, utilizing the eigen decomposition for parameterization.
#'
#' @param n An integer specifying the number of samples to draw.
#' @param mu A numeric vector of length \code{p} representing the mean of the distribution.
#' @param Sigma A \code{p} x \code{p} precision matrix (the inverse of the covariance matrix).
#'
#' @return A \code{p} x \code{n} matrix of samples drawn from the multivariate normal distribution with mean vector \code{mu} and precision matrix \code{Sigma}.
#'
#' @details
#' This function uses the eigen decomposition of the precision matrix \code{Sigma} to generate samples from the corresponding multivariate normal distribution. The process is as follows:
#' - Perform the truncated eigen decomposition of the precision matrix \code{Sigma}.
#' - Draw samples from a standard normal distribution.
#' - Transform these samples to have the specified mean \code{mu} and precision structure \code{Sigma}.
#'
#' Note that the precision matrix is the inverse of the covariance matrix. The function inverts the square root transformation to properly adjust for the precision parameterization.
#'
#' @export
#'
#' @examples
#' # Example parameters
#' mu <- c(0, 0) # Mean vector
#' Precision <- matrix(c(2, 0.5, 0.5, 1), nrow = 2) # Precision matrix (inverse covariance)
#' # Generate 100 samples
#' samples <- rMVNormP_eigen(n = 100, mu = mu, Sigma = Precision
rMVNormP_eigen <- function(n, mu, Sigma) {
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv = TRUE)
  k <- ncol(L)
  Z <- matrix(rnorm(k * n), k, n)
  X <- L %*% Z
  X <- sweep(X, 1, mu, FUN = `+`)
}






#' Logistic (Sigmoid) function
#'
#' This function computes the logistic function, which transforms the input \code{x} into a value between 0 and 1. The logistic function is widely used to model probabilities and is defined as \code{1 / (1 + exp(-x))}.
#'
#' @param x A numeric vector of values.
#'
#' @return A numeric vector of the same length as \code{x}, with values between 0 and 1.
#'
#' @details The logistic function is used to map real-valued numbers to the (0, 1) interval, making it useful for modeling probabilities in various statistical models, such as binary classification or logistic regression.
#'
logistic <- function(x) {
  1 / (1 + exp(-x))
}






