library(mvtnorm)
library(extraDistr)
library(mvQuad)
library(docstring)

dtbt <- function(x, df, rho) {
  if(NCOL(x) != 2) stop("x must be a matrix with two columns.")
  
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  ind <- apply(x, 1, function(z) any(z < 0))
  dens_x <- dmvt(x = x, 
                 sigma = Sigma,
                 df = df,
                 log = F) /
    pmvt(lower = c(0, 0), sigma = Sigma, df = df, keepAttr = FALSE)
  dens_x[ind] <- 0
  return(dens_x)
}

ptbt <- function(x, df, rho) {
  if(NCOL(x) != 2) stop("x must be a matrix with two columns.")
  
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  ind <- apply(x, 1, function(z) any(z < 0))
  prob_x <- apply(x, 
                  1, 
                  function(z) pmvt(lower = c(0, 0),
                                   upper = z,
                                   sigma = Sigma,
                                   df = df,
                                   keepAttr = FALSE)) /
    pmvt(lower = c(0, 0), sigma = Sigma, df = df, keepAttr = FALSE)
  prob_x[ind] <- 0
  return(prob_x)
}

qtbt <- function(p, 
                 df, 
                 rho, 
                 type = c("or", "and", "density"),
                 x1 = NULL,
                 m = 1000,
                 tol = .Machine$double.eps^0.25) {
  #' @title 
  #' Extremal regions for a truncated bivariate t-distribution
  #' 
  #' @description
  #' Returns a set of points representing the extremal regions 
  #' \eqn{Q^{(1)}}, \eqn{Q^{(2)}} or \eqn{Q^{(3)}} 
  #' of a truncated bivariate t-distribution.
  #' 
  #' @param p desired tail quantile level.
  #' @param df degree of freedom of the bivariate t-distribution.
  #' @param rho correlation of the underlying multivariate Gaussian.
  #' @param type type of the extremal region. See details.
  #' @param x1 vector of points of the first coordinate for the representative set.
  #' If \code{x1} is \code{NULL}, 2 * (\code{m} - 1) points in a suitable range will be used.
  #' @param m specifies the number of points in the representative set. 
  #' Will be ignored if \code{x1} is not \code{NULL}.
  #' @param tol tolerance level for \code{uniroot}.
  #' 
  #' @return
  #' \code{Q}: the set of points representing the extremal region.
  #' 
  #' @details
  #' The density of truncated bivariate t-distribution is given by
  #' \eqn{f(x_1, x_2) = \frac{t(x_1, x_2)}{T(0, 0)}, (x_1, x_2) \in \mathbb{R}_{+}^2} where
  #' \eqn{t} and \eqn{T} are the pdf and cdf of a t-distribution.
  #' The scale matrix  \eqn{\Sigma} is of the form \eqn{\Sigma_{ii} = 1}, \eqn{\Sigma_{ij} = \rho}.
  #' The three extremal regions are given by
  #' \itemize{
  #'    \item \eqn{Q^{(1)} = \{ (x_1, x_2) : P(X_1 > x_1 \cup X_2 > x_2) = p\}}
  #'    \item \eqn{Q^{(2)} = \{ (x_1, x_2) : P(X_1 > x_1 \cap X_2 > x_2) = p\}}
  #'    \item \eqn{Q^{(3)} = \{ (x_1, x_2) : f(x_1, x_2) \leq \beta \}} s.t. \eqn{P(Q^{(3)}) = p}
  #' }
  #' The function \code{uniroot} is used to find the representative points for the
  #' extremal sets.

  if (!is.null(x1)) {
    if (any(x1 < 0)) stop("x1 has to be non-negative.")
  }
  type = match.arg(type)
  
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  p_0_0 <- pmvt(lower = c(0, 0), sigma = Sigma, df = df, keepAttr = FALSE)

  
  
  
  if (type == "or") {
    p_star <- (1 - p) * p_0_0
    
    lb <- 0
    ub <- min(10^6, qt(1 - .Machine$double.eps, df = df)^2)
    interval <- c(lb, ub)
    
    rbt <- function(x1, x2) {
      ret <- pmvt(lower = c(0, 0), upper = c(x1, x2), sigma = Sigma, df = df, keepAttr = FALSE) - p_star
      return(ret)
    }
    
    x1_min <- uniroot(function(x) rbt(x, Inf), interval, extendInt = "upX", tol = tol)$root + 10 * tol
    
    ub <- uniroot(function(x) rbt(x1_min, x), interval, extendInt = "upX", tol = tol)$root
    
    if (is.null(x1)) {
      x_star <- uniroot(function(x) rbt(x, x), interval, extendInt = "upX", tol = tol)$root
      
      x1_grid <-  seq(x_star, x1_min, length.out = m)
      x2_grid <- c(x_star, rep(NA, m - 2), ub)
      
      if (m > 2) {
        for (i in seq(2, m - 1)) {
          lb <- x2_grid[i - 1]
          interval <- c(lb, ub)
          x2_grid[i] <- uniroot(function(x) rbt(x1_grid[i], x), interval, extendInt = "upX", tol = tol)$root
        }
      }
      
      ind <- order(x1_grid)
      
      x_grid <- c(x1_grid[ind], x2_grid[-1])
      Q <- cbind(x_grid, rev(x_grid))
      colnames(Q) <- c("x1_hat", "x2_hat")
      return(Q)
      
    } else {
      order_x1 <- order(x1, decreasing = T)
      x1_s <- x1[order_x1]
      ind_min <- (x1_s >= x1_min)
      x1_grid <- x1_s[ind_min]
      m <- sum(ind_min)
      x2_grid <- rep(Inf, length(x1))
      
      for (i in seq_len(m)) {
        if (i == 1) lb <- x1_min
        else lb <- x2_grid[i - 1]
        interval <- c(lb, ub)
        x2_grid[i] <- uniroot(function(x) rbt(x1_grid[i], x), interval, extendInt = "upX", tol = tol)$root
      }
      
      Q <- cbind(x1, x2_grid[order(order_x1)])
      colnames(Q) <- c("x1_hat", "x2_hat")
      return(Q)
      
    }
  } else if (type == "and") {
    p_star <- p * p_0_0
    lb <- 0
    ub <- qt(p_star, df = df, lower.tail = F)
    interval <- c(lb, ub)
    rbt <- function(x1, x2) {
      ret <- pmvt(lower = c(x1, x2), sigma = Sigma, df = df, keepAttr = FALSE) - p_star
      return(ret)
    }
    
    x1_max <- uniroot(function(x) rbt(x, 0), interval, extendInt = "downX", tol = tol)$root
    ub <- x1_max
    interval <- c(lb, ub)
    
    if (is.null(x1)) {
      x_star <- uniroot(function(x) rbt(x, x), interval, extendInt = "downX", tol = tol)$root
      
      x1_grid <-  seq(x_star, ub, length.out = m)
      x2_grid <- c(ub, rep(NA, m - 2), 0)
      
      if (m > 2) {
        for (i in seq(2, m - 1)) {
          ub <- x2_grid[i - 1]
          interval <- c(lb, ub)
          x2_grid[i] <- uniroot(function(x) rbt(x1_grid[i], x), interval, extendInt = "downX", tol = tol)$root
        }
      }
      ind <- order(x2_grid)
      
      x_grid <- c(x2_grid[ind][-m], x1_grid)
      Q <- cbind(x_grid, rev(x_grid))
      colnames(Q) <- c("x1_hat", "x2_hat")
      return(Q)
    } else {
      order_x1 <- order(x1, decreasing = F)
      x1_s <- x1[order_x1]
      ind_max <- (x1_s <= x1_max)
      x1_grid <- x1_s[ind_max]
      m <- sum(ind_max)
      x2_grid <- rep(NA, length(x1))
      
      for (i in seq_len(m)) {
        if (i == 1) ub <- x1_max
        else ub <- x2_grid[i - 1]
        interval <- c(lb, ub)
        x2_grid[i] <- uniroot(function(x) rbt(x1_grid[i], x), interval, extendInt = "downX", tol = tol)$root
      }
      
      Q <- cbind(x1, x2_grid[order(order_x1)])
      colnames(Q) <- c("x1_hat", "x2_hat")
      return(Q)
      }
    } else {
    Q <- ellipse(sigma = Sigma, df = df, prob = 1 - p, npoints = m, pos = TRUE)
    colnames(Q) <- c("x1_hat", "x2_hat")
    return(Q)
  }
}

rtbt <- function(n, df, rho) {
  #' @title 
  #' Pseudo-random samples of a truncated bivariate t-distribution
  #' 
  #' @description
  #' Samples truncated bivariate t-distribution using rejection sampling.
  #' 
  #' @param n number of samples.
  #' @param df degree of freedom of the bivariate t-distribution.
  #' @param rho correlation of the underlying multivariate Gaussian.
  #' 
  #' @return
  #' A matrix of dimension (n x 2) pseudo-random samples of a truncated bivariate t-distribution
  #' 
  #' @details
  #' The density of truncated bivariate t-distribution is given by
  #' \eqn{f(x_1, x_2) = \frac{t(x_1, x_2)}{T(0, 0)}, (x_1, x_2) \in \mathbb{R}_{+}^2} where
  #' \eqn{t} and \eqn{T} are the pdf and cdf of a t-distribution.
  #' The scale matrix  \eqn{\Sigma} is of the form \eqn{\Sigma_{ii} = 1}, \eqn{\Sigma_{ij} = \rho}.
  #' The package \code{mvtnorm} is used to generate the samples.
  
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  p_0_0 <- pmvt(lower = c(0, 0), sigma = Sigma, df = df, keepAttr = FALSE)
  N <- ceiling(n / p_0_0)
  
  Z <- rmvt(N, sigma = Sigma, df = df)
  ind <- apply(Z, 1, FUN = function(z) (z[1] > 0) & (z[2] > 0))
  X <- Z[ind, , drop = FALSE]
  
  n_eff <- nrow(X)
  if (n_eff >= n) return(X[1:n, , drop = FALSE])
  else {
    return(rbind(X, rtbt(n - n_eff, df, rho)))
  }
}
