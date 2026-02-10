cv_h_lik <- function(Y, k_grid = c(4: 10), nfolds = 10, q = 0.95) {
  rdata <- rowSums(Y)
  u <- quantile(rdata, probs = q)
  extdata <- Y[rdata >= u, ]
  
  L <- nrow(extdata)
  ind <- sample(seq_len(L))
  M <- ceiling(L / nfolds)
  llik <- rep(0, length(k_grid))
  
  for (i in seq_len(nfolds)) {
    ind_fold <- ind[((i - 1) * M + 1) : min((i * M), L)]
    
    f <- function(k) fExtDep.np(method = "Frequentist",
                                data = extdata[-ind_fold, ],
                                u = 0,
                                mar.fit = FALSE,
                                type = "rawdata", 
                                k0 = k)$Ahat$beta
    w <- extdata[ind_fold, 1] / rowSums(extdata[ind_fold, ])
    for (j in seq_len(length(k_grid))) {
      beta <- f(k_grid[j])
      lik <- sapply(1:length(w), function(i)  ExtremalDep:::dh(w[i], beta))
      lik[lik < 0] <- 0
      llik[j] <- llik[j] + sum(log(lik))
    }
  }
  
  k <- min(k_grid[which(llik == max(llik))])
  return(k)
}

cv_h_emp <- function(Y, nu_grid, nfolds = 10, tau = 0.95, 
                  method = "euclidean", raw = TRUE) {
  if (any(nu_grid <= 0))
    stop("nu_grid must be postive")
  if (tau >= 1) 
    stop("tau must be smaller than 1")
  if (dim(Y)[2] > 2) 
    stop("this function only applies to bivariate extremes")
  
  
  if (raw == TRUE) {
    n <- dim(Y)[1]
    FX <- ecdf(Y[, 1])
    FY <- ecdf(Y[, 2])
    x <- -1/log(n/(n + 1) * FX(Y[, 1]))
    y <- -1/log(n/(n + 1) * FY(Y[, 2]))
  }
  if (raw == FALSE) {
    x <- Y[, 1]
    y <- Y[, 2]
  }
  
  w <- x/(x + y)
  u <- quantile(x + y, tau)
  Y0 <- cbind(x, y)[which(x + y > u), ]
  w0 <- w[which(x + y > u)]
  k0 <- length(w0)
  
  ind <- sample(seq_len(k0))
  M <- ceiling(k0 / nfolds)
  llik <- rep(0, length(nu_grid))
  
  for (i in seq_len(nfolds)) {
    ind_fold <- ind[((i - 1) * M + 1) : min((i * M), k0)]
    Y <- rbind(rep(u, 2)/4, Y0[-ind_fold, ])
    llik <- llik + 
      sapply(nu_grid, 
             FUN = function(nu) sum(log(angdensity(Y, 
                                                   tau = 0, 
                                                   nu = nu,
                                                   grid = w0[ind_fold],
                                                   method = method, 
                                                   raw = FALSE)$h)))
  }
  v <- min(nu_grid[which(llik == max(llik))])
  return(v)
}