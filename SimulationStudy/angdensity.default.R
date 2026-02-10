angdensity.default <- function (Y, tau = 0.95, nu, grid = seq(0.01, 0.99, length = 2^8), 
                                method = "euclidean", raw = TRUE) 
{
  if (tau >= 1) 
    stop("tau must be smaller than 1")
  if (nu < 0) 
    stop("nu must be positive")
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
  w <- w[which(x + y > u)]
  k <- length(w)
  if (method == "euclidean") {
    #v <- var(w) * (k/(k - 1))
    v <- var(w) * ((k - 1) / k)
    #p <- 1/k * (1 - (mean(w) - 1/2) * v^{-1} * (w - 1/2))
    p <- 1/k * (1 - (mean(w) - 1/2) * v^{-1} * (w - mean(w)))
  }
  if (method == "empirical") {
    lambda <- emplik:::el.test(w, 1/2)$lambda
    p <- 1/k * (1 + (w - 1/2) * lambda)^-1
  }
  smooth <- function(W) sum(p * dbeta(W, w * nu, (1 - w) * 
                                        nu))
  h <- sapply(grid, smooth)
  outputs <- list(h = h, w = w, nu = nu, grid = grid, Y = Y, p = p)
  class(outputs) <- "angdensity"
  return(outputs)
}
