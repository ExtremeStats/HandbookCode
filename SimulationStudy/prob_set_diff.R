source("truncated_bivariate_t_distribution.R")

prob_set_diff <- function(Q_hat, 
                          p,
                          rho, 
                          df,
                          type = c("or", "and", "density"),
                          quadrature_rule = c("cNC1", "cNC2", "cNC3", "cNC4", "cNC5", "cNC6",
                                              "oNC0", "oNC1", "oNC2", "oNC3",
                                              "GLe", "GKr", "nLe", "Leja"),
                          level = 1000,
                          tol = .Machine$double.eps^0.25) {
  #' @title 
  #' Set difference probability.
  #' 
  #' @description
  #' Approximates the probability of the difference between the estimated extremal region and the
  #' true extremal region. 
  #' 
  #' @param Q_hat the set of points representing the estimated extremal region
  #' @param p nominal tail quantile level
  #' @param df degree of freedom of the bivariate t-distribution.
  #' @param rho correlation of the underlying multivariate Gaussian.
  #' @param type type of the extremal region. See details.
  #' @param quadrature_rule rule of the quadrature approximation. See details.
  #' @param level accuracy level for the quadrature rule.
  #' @param tol tolerance level for \code{uniroot}.
  #' 
  #' @return
  #' The approximated set difference probability
  #' 
  #' @details
  #' The density of truncated bivariate t-distribution is given by
  #' \eqn{f(x_1, x_2) = \frac{t(x_1, x_2)}{T(0, 0)}, (x_1, x_2) \in \mathbb{R}_{+}^2} where
  #' \eqn{t} and \eqn{T} are the pdf and cdf of a t-distribution.
  #' The scale matrix  \eqn{\Sigma} is of the form \eqn{\Sigma_{ii} = 1}, \eqn{\Sigma_{ij} = \rho}.
  #' The package \code{mvtnorm} is used to generate the samples.
  #' A quadrature approximation implemented in the package \code{mvQuad} is used to estimate
  #' the set difference probability. Linear intrapolation is used to obtain the
  #' estimated extremal region using the set of representative points in \code{Q_hat}.
  
  type = match.arg(type)
  quadrature_rule <- match.arg(quadrature_rule )
  
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  p_0_0 <- pmvt(lower = c(0, 0), sigma = Sigma, df = df, keepAttr = FALSE)
  
  Q_hat <- Q_hat[complete.cases(Q_hat), ]
  ind_x1_pos <- (Q_hat[, 1] >= 0)
  ind_x2_pos <- (Q_hat[, 2] >= 0)
  ind_pos <- ind_x1_pos & ind_x2_pos
  Q_hat <- Q_hat[ind_pos, ]
  
  if (type == "or") {
    p_star <- (1 - p) * p_0_0
    
    lb <- 0
    ub1 <- qt(1 - p * 10^(-3), df = df)
    ub2 <- qt(1 - .Machine$double.eps, df = df)^2
    interval <- c(lb, ub2)
    interval1 <- c(lb, ub1)
    interval2 <- c(ub1, ub2)
    
    rbt <- function(x1, x2) {
      ret <- pmvt(lower = c(0, 0), upper = c(x1, x2), sigma = Sigma, df = df, keepAttr = FALSE) - p_star
      return(ret)
    }
    
    x1_min <- uniroot(function(x) rbt(x, Inf), interval, extendInt = "upX", tol = tol)$root + 10 * tol
    
    
    interval1 <- c(max(0 ,min(c(Q_hat[, 1], x1_min))), ub1)
    interval2 <- c(ub1, ub2)
    
    nw <- createNIGrid(dim = 1, type = quadrature_rule, level = level)
    rescale(nw, interval1)
    
    f <- function(x1) {
      l <- qtbt(p = p, df = df, rho = rho, type = type, x1 = x1, tol = tol)
      y_hat <- sapply(seq_len(length(x1)), 
                       FUN = function(i) {
                         ind_b <- x1[i] < Q_hat[, 1]
                         ind_l <- x1[i] >= Q_hat[, 1]
                          if(all(ind_b)) {
                            x_l <- 0
                            y_l <- Inf
                          } else {
                            x_l <- max(Q_hat[ind_l, 1])
                            y_l <- Q_hat[which(Q_hat[, 1] == x_l), 2]
                          }
                         ind_b <- x1[i] > Q_hat[, 1]
                         ind_u <- x1[i] <= Q_hat[, 1]
                         if(all(ind_b)) {
                           x_u <- x1[i]
                           y_u <- Q_hat[which(Q_hat[, 1] == max(Q_hat[, 1])), 2]
                         } else {
                           x_u <- min(Q_hat[ind_u, 1])
                           y_u <- Q_hat[which(Q_hat[, 1] == x_u), 2]
                         }
                         w <- (x1[i] - x_l) / (x_u - x_l)
                         y_hat <- (1 - w) * y_l + w * y_u
                         return(y_hat)
      })
      y_true <- l[, 2]
      Q_diff <- t(apply(cbind(y_true, y_hat), 1, sort))
      cond_pt(Q_diff, x1, Sigma, df, tol = tol) * dt(x1, df = df)
    }
    
    A1 <- quadrature(f, nw)
    
    nw <- createNIGrid(dim = 1, type = quadrature_rule, level = level)
    rescale(nw, interval2)
    A2 <- quadrature(f, nw)
    
    return((A1 + A2) / p_0_0)
  } else if (type == "and") {
    p_star <- p * p_0_0
    
    lb <- 0
    ub <- qt(p_star, df = df, lower.tail = F)
    
    interval <- c(0, max(c(Q_hat[, 1], ub)))
    nw <- createNIGrid(dim = 1, type = quadrature_rule, level = level)
    rescale(nw, interval)
    f <- function(x1) {
      l <- qtbt(p = p, df = df, rho = rho, type = type, x1 = x1, tol = tol)
      y_hat <- sapply(seq_len(length(x1)), 
                      FUN = function(i) {
                        ind_b <- x1[i] < Q_hat[, 1]
                        ind_l <- x1[i] >= Q_hat[, 1]
                        if(all(ind_b)) {
                          x_l <- x1[i]
                          y_l <- min(Q_hat[which(Q_hat[, 1] == min(Q_hat[, 1])), 2])
                        } else {
                          x_l <- max(Q_hat[ind_l, 1])
                          y_l <- min(Q_hat[which(Q_hat[, 1] == x_l), 2])
                        }
                        ind_b <- x1[i] > Q_hat[, 1]
                        ind_u <- x1[i] <= Q_hat[, 1]
                        if(all(ind_b)) {
                          x_u <- x1[i]
                          y_u <- 0
                        } else {
                          x_u <- min(Q_hat[ind_u, 1])
                          y_u <- max(Q_hat[which(Q_hat[, 1] == x_u), 2])
                        }
                        w <- (x1[i] - x_l) / (x_u - x_l)
                        y_hat <- (1 - w) * y_l + w * y_u
                        return(y_hat)
                        })
      y_true <- l[, 2]
      y_true[is.na(l[, 2])] <- 0
      Q_diff <- t(apply(cbind(y_true, y_hat), 1, sort))
      cond_pt(Q_diff, x1, Sigma, df, tol = tol) * dt(x1, df = df)
    }
    
    A <- quadrature(f, nw)
    return(A / p_0_0)
  } else {
    Q_hat <- Q_hat[complete.cases(Q_hat), ]
    
    Q_true <- ellipse(sigma = Sigma, df = df, prob = 1 - p, npoints = level, pos = TRUE)

    interval <- c(0, max(Q_true[, 1], Q_hat[, 1]))
    
    nw <- createNIGrid(dim = 1, type = quadrature_rule, level = level)
    rescale(nw, interval)
    
    W_hat <- Q_hat[, 1] / rowSums(Q_hat)
    W_hat <- W_hat[!is.na(W_hat)]
    ind_w <- (W_hat > 0.5)
    Q_hat_lower <- Q_hat[ind_w, ]
    Q_hat_upper <- Q_hat[!ind_w, ]
    
    
    
    f <- function(x1) {
      y_hat <- t(sapply(seq_len(length(x1)),
                      FUN = function(i) {
                        #lower half
                        ind_u <- x1[i] < Q_hat_lower[, 1]
                        ind_l <- x1[i] >= Q_hat_lower[, 1]
                        
                        if (all(ind_u)) {
                          y_hat_lower <- 0
                        } else if (all(ind_l)) {
                          y_hat_lower <- mean(max(Q_hat_lower[, 2]), min(Q_hat_upper[, 2]))
                        } else {
                          x_l <- max(Q_hat_lower[ind_l, 1])
                          x_u <- min(Q_hat_lower[ind_u, 1])
                          
                          a <- (x1[i] - x_l) / (x_u - x_l)
                          pos_x_l <- min(which(Q_hat_lower[, 1] == x_l))
                          pos_x_u <- min(which(Q_hat_lower[, 1] == x_u))
                          y_hat_lower <- (1 - a)  * Q_hat_lower[pos_x_l, 2] +  
                            a * Q_hat_lower[pos_x_u, 2]
                        }
                        
                        #upper half
                        
                        ind_u <- x1[i] < Q_hat_upper[, 1]
                        ind_l <- x1[i] >= Q_hat_upper[, 1]
                        
                        if (all(ind_u)) {
                          y_hat_upper <- Q_hat_upper[min(which(Q_hat_upper[, 1] == min(Q_hat_upper[, 1]))), 2]
                        } else if (all(ind_l)) {
                          y_hat_upper <- mean(max(Q_hat_lower[, 2]), min(Q_hat_upper[, 2]))
                        } else {
                          x_l <- max(Q_hat_upper[ind_l, 1])
                          x_u <- min(Q_hat_upper[ind_u, 1])
                          
                          a <- (x1[i] - x_l) / (x_u - x_l)
                          pos_x_l <- min(which(Q_hat_upper[, 1] == x_l))
                          pos_x_u <- min(which(Q_hat_upper[, 1] == x_u))
                          y_hat_upper <- (1 - a)  * Q_hat_upper[pos_x_l, 2] +  
                            a * Q_hat_upper[pos_x_u, 2]
                        }
                        return(c(y_hat_lower, y_hat_upper))
                      }))
      y_hat_lower <- y_hat[, 1]
      y_hat_upper <- y_hat[, 2]
      
      
      W_true <- Q_true[, 1] / rowSums(Q_true)
      ind_w <- W_true > 0.5
      
      Q_true_lower <- Q_true[ind_w, ]
      Q_true_upper <- Q_true[!ind_w, ]
      
      y_true <- t(sapply(seq_len(length(x1)),
                        FUN = function(i) {
                          #lower half
                          ind_u <- x1[i] < Q_true_lower[, 1]
                          ind_l <- x1[i] >= Q_true_lower[, 1]
                          
                          
                          if (all(ind_u)) {
                            y_true_lower <- 0
                          } else if (all(ind_l)) {
                            y_true_lower <- mean(max(Q_true_lower[, 2]), min(Q_true_upper[, 2]))
                          } else {
                            x_l <- max(Q_true_lower[ind_l, 1])
                            x_u <- min(Q_true_lower[ind_u, 1])
                            
                            a <- (x1[i] - x_l) / (x_u - x_l)
                            pos_x_l <- min(which(Q_true_lower[, 1] == x_l))
                            pos_x_u <- min(which(Q_true_lower[, 1] == x_u))
                            y_true_lower <- (1 - a)  * Q_true_lower[pos_x_l, 2] +  
                              a * Q_true_lower[pos_x_u, 2]
                          }
                          
                          #upper half
                          
                          ind_u <- x1[i] < Q_true_upper[, 1]
                          ind_l <- x1[i] >= Q_true_upper[, 1]
                          
                          if (all(ind_u)) {
                            y_true_upper <- Q_true_upper[min(which(Q_true_upper[, 1] == min(Q_true_upper[, 1]))), 2]
                          } else if (all(ind_l)) {
                            y_true_upper <- mean(max(Q_true_lower[, 2]), min(Q_true_upper[, 2]))
                          } else {
                            x_l <- max(Q_true_upper[ind_l, 1])
                            x_u <- min(Q_true_upper[ind_u, 1])
                            
                            a <- (x1[i] - x_l) / (x_u - x_l)
                            pos_x_l <- min(which(Q_true_upper[, 1] == x_l))
                            pos_x_u <- min(which(Q_true_upper[, 1] == x_u))
                            y_true_upper <- (1 - a)  * Q_true_upper[pos_x_l, 2] +  
                              a * Q_true_upper[pos_x_u, 2]
                          }
                          return(c(y_true_lower, y_true_upper))
                        }))
      y_true_lower <- y_true[, 1]
      y_true_upper <- y_true[, 2]
      
      
      N <- length(x1)
      Q_diff1 <- cbind(sapply(seq_len(N), FUN = function(i) max(c(y_true_upper[i], y_hat_lower[i]))),
                       sapply(seq_len(N), FUN = function(i) max(c(y_true_upper[i], y_hat_upper[i]))))
      Q_diff2 <- cbind(sapply(seq_len(N), FUN = function(i) max(c(y_hat_upper[i], y_true_lower[i]))),
                       sapply(seq_len(N), FUN = function(i) max(c(y_hat_upper[i], y_true_upper[i]))))
      Q_diff3 <- cbind(sapply(seq_len(N), FUN = function(i) min(c(y_true_lower[i], y_hat_lower[i]))),
                       sapply(seq_len(N), FUN = function(i) min(c(y_true_lower[i], y_hat_upper[i]))))
      Q_diff4 <- cbind(sapply(seq_len(N), FUN = function(i) min(c(y_hat_lower[i], y_true_lower[i]))),
                       sapply(seq_len(N), FUN = function(i) min(c(y_hat_lower[i], y_true_upper[i]))))
      (cond_pt(Q_diff1, x1, Sigma, df, tol = tol) + cond_pt(Q_diff2, x1, Sigma, df, tol = tol) +
          cond_pt(Q_diff3, x1, Sigma, df, tol = tol) + cond_pt(Q_diff4, x1, Sigma, df, tol = tol)) * 
        dt(x1, df = df)
    }
    A <- quadrature(f, nw)
    return(A / p_0_0)
  }
}
