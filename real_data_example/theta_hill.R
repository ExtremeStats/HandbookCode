theta_hill <- function(X, k) {
  n <- length(X)
  Z <- sort(X)
  xi_hat <- 1 / k * sum(log(Z[(n - k + 1) : n]) -log(Z[n - k]))
  b_hat <- quantile(Z, probs = 1 - k / n)
  a_hat <- xi_hat * b_hat
  theta_hat <- c(a_hat, b_hat,  xi_hat)
  return(theta_hat)
}