bbeed.thresh.mar <- function (data, cov1 = as.matrix(rep(1, nrow(data))), 
                              cov2 = as.matrix(rep(1, nrow(data))), U, mar = TRUE, 
                              par10, par20, sig10, sig20, 
                              param0 = NULL, k0 = NULL, pm0 = NULL, prior.k = "nbinom", 
                              prior.pm = "unif", nk = 70, 
                              hyperparam = list(mu.nbinom = 3.2, var.nbinom = 4.48), 
                              nsim = NULL, warn = FALSE) 
{
  if (nrow(cov1) != nrow(data) || nrow(cov2) != nrow(data)) {
    stop("Wrong dimensions of the covariates")
  }
  if (length(par10) != (ncol(cov1) + 2)) {
    stop("Wrong length of initial first marginal parameter vector")
  }
  if (length(par20) != (ncol(cov2) + 2)) {
    stop("Wrong length of initial second marginal parameter vector")
  }
  kn <- c(sum(data[, 1] > U[1]), sum(data[, 2] > U[2]))/nrow(data)
  if (mar) {
    message("\n Preliminary on margin 1 \n")
    mar1.fit <- ExtremalDep:::fGEV(data = data[, 1], par.start = par10, 
                                   u = U[1], method = "Bayesian", cov = cov1, sig0 = sig10, 
                                   nsim = 50000)
    par10 <- apply(mar1.fit$param_post[-c(1:30000), ], 2, 
                   mean)
    sig10 <- tail(mar1.fit$sig.vec, 1)
    message("\n Preliminary on margin 2 \n")
    mar2.fit <- ExtremalDep:::fGEV(data = data[, 2], par.start = par20, 
                                   u = U[2], method = "Bayesian", cov = cov2, sig0 = sig20, 
                                   nsim = 50000)
    par20 <- apply(mar2.fit$param_post[-c(1:30000), ], 2, 
                   mean)
    sig20 <- tail(mar2.fit$sig.vec, 1)
  }
  message("\n Estimation of the extremal dependence and margins \n")
  mcmc <- ExtremalDep:::cens.bbeed.mar(data = data, cov1 = cov1, cov2 = cov2, 
                                       pm0 = pm0, param0 = param0, k0 = k0, par10 = par10, 
                                       par20 = par20, sig10 = sig10, sig20 = sig20, kn = kn, 
                                       prior.k = prior.k, prior.pm = prior.pm, hyperparam = hyperparam, 
                                       nsim = nsim, U = U, p.star = 0.234, warn = warn,
                                       nk = nk)
  return(mcmc)
}