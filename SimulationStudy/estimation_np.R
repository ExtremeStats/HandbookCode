# A_np <- function(w, X, k, method = c("c", "b")) {
#   method <- match.arg(method)
#   R_X1 <- order(X[, 1])
#   R_X2 <- order(X[, 2])
#   
#   n <- nrow(X)
#   Z1 <- n / (n + 1 - R_X1)
#   Z2 <- n / (n + 1 - R_X2)
#   
#   S <- Z1 + Z2
#   W <- Z1 / S
#   
#   ind <- (S >= n/k)
#   N <- sum(ind)
#   U <- W[ind]
#   
#   if (method == "c") {
#     U_bar <- mean(U)
#     S <- var(U) * (N - 1) / N
#     weights <- (1 - (U_bar - 0.5) / S * (U - U_bar)) / N
#   } else {
#     weights <- 1 / N
#   }
# 
# 
#   A_w <- sapply(w, FUN = function(x) 2 * sum(weights * sapply(U, FUN = function(u) max(c(u * x, (1 - u) * (1 - x))))))
#   #sapply(A_w, FUN = function(z) min(c(1, max(c(0.5, z)))))
#   return(A_w)
# }

A_np <- function(w, X, k, smoothed = FALSE, v = exp(seq(1, 10, by = 0.5))) {
  R_X1 <- rank(X[, 1])
  R_X2 <- rank(X[, 2])
  
  n <- nrow(X)
  Z1 <- n / (n + 1 - R_X1)
  Z2 <- n / (n + 1 - R_X2)
  
  S <- Z1 + Z2
  W <- Z1 / S
  
  ind <- (S >= n/k)
  N <- sum(ind)
  U <- W[ind]
  
  U_bar <- mean(U)
  S <- var(U) * (N - 1) / N
  weights <- (1 - (U_bar - 0.5) / S * (U - U_bar)) / N
  if (smoothed) {
    if (length(v) > 1) v <- cv_h(v, U)
    
    f <- function(z) sapply(U, FUN = function(u) pbeta(z, u * v, (1 - u) * v)) %*% weights
    
    m <- length(w)
    ind <- order(w)
    A <- rep(0, m)
    for (i in seq_len(m)) {
      if (i == 1) {
        A[i] <- w[ind[i]] + 
          2 * integrate(f, 0, 1 - w[ind[i]], abs.tol = 10^(-4))$value
      } else {
        A[i] <- 
          A[i - 1] + (w[ind[i]] - w[ind[i - 1]]) - 
          2 * integrate(f, 1 - w[ind[i]], 1 - w[ind[i - 1]], abs.tol = 10^(-4))$value
      }
    }
  } else {
    A <- 
      sapply(w, FUN = function(x) 2 * 
               sum(weights * sapply(U, FUN = function(u) max(c(u * x, (1 - u) * (1 - x))))))
  }
  return(A)
}

H_np <- function(w, X, k, 
                 smoothed = FALSE,
                 v = exp(seq(1, 10, by = 0.5))) {
  R_X1 <- rank(X[, 1])
  R_X2 <- rank(X[, 2])
  
  n <- nrow(X)
  Z1 <- n / (n + 1 - R_X1)
  Z2 <- n / (n + 1 - R_X2)
  
  S <- Z1 + Z2
  W <- Z1 / S
  
  ind <- (S >= n/k)
  N <- sum(ind)
  U <- W[ind]
  
  U_bar <- mean(U)
  S <- var(U) * (N - 1) / N
  weights <- (1 - (U_bar - 0.5) / S * (U - U_bar)) / N
  
  if (smoothed) {
    if (length(v) > 1) v <- cv_h(v, U)
    sapply(U, FUN = function(u) pbeta(w, u * v, (1 - u) * v)) %*% weights
  } else {
    sapply(w, FUN = function(x) sum(weights * (U <= x)))
  }
}

h_np <- function(w, X, k,
                 v = exp(seq(1, 10, by = 0.5))) {
  R_X1 <- rank(X[, 1])
  R_X2 <- rank(X[, 2])
  
  n <- nrow(X)
  Z1 <- n / (n + 1 - R_X1)
  Z2 <- n / (n + 1 - R_X2)
  
  S <- Z1 + Z2
  W <- Z1 / S
  
  ind <- (S >= n/k)
  N <- sum(ind)
  U <- W[ind]
  
  U_bar <- mean(U)
  S <- var(U) * (N - 1) / N
  weights <- (1 - (U_bar - 0.5) / S * (U - U_bar)) / N
  
  if (length(v) > 1) v <- cv_h(v, U)
  
  sapply(w, FUN = function(x) weights %*% sapply(U, FUN = function(u) dbeta(x, u * v, (1 - u) * v)))
}

cv_h <- function(v_grid, U, k = 10) {
  h_np0 <- function(w, U, weights, v) 
    sapply(w, FUN = function(x) weights %*% sapply(U, FUN = function(u) dbeta(x, u * v, (1 - u) * v)))
  L <- length(U)
  ind <- sample(seq_len(L))
  M <- ceiling(L / k)
  h_dens <- rep(0, length(v_grid))
  
  for (i in seq_len(k)) {
    ind_fold <- ind[((i - 1) * M + 1) : min((i * M), L)]
    N <- L - length(ind_fold)
    U_bar_cv <- mean(U[-ind_fold])
    S_cv <- var(U[-ind_fold]) * (N - 1) / N
    weights_cv <- (1 - (U_bar_cv - 0.5) / S_cv * (U[-ind_fold] - U_bar_cv)) / N
    h_dens <- h_dens + sapply(v_grid, FUN = function(v) sum(log(h_np0(U[ind_fold], U[-ind_fold], weights_cv, v))))
  }
  v <- min(v_grid[which(h_dens == max(h_dens))])
  return(v)
}
