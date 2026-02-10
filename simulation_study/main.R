# This script runs the simulations performed in subsection "Simulation study".
# WARNING: This will take a long time. To run it on Windows, the variable
#  no.cores at line 44 has to be set to 1.

source("functions_simulation_study.R")

assignInNamespace("angdensity.default", angdensity.default, pos="package:extremis")
assignInNamespace("constmle", constmle, pos="package:ExtremalDep")
assignInNamespace("cens.bbeed.mar", cens.bbeed.mar, pos="package:ExtremalDep")

#### Parameters of the simulation setting ####
###### number of simulation runs
nsim <- 1500
###### number of samples to generate
n <- 1500
###### correlation of the underlying bivariate gaussian distribution
rho <- 0
###### degrees of freedom
df <- 1
###### extremal quantiles to estimate
q_n <- c(2 / n, 1 /n, 1 / (2 * n))
###### quantile specifying the threshold to use
q <- 0.9

m <- 100
w <- seq(1/m, 1 - 1/m, length.out = m)

parameters_simulation_setting <- list(
  nsim = nsim,
  n = n,
  rho = rho,
  df = df,
  q_n = q_n,
  q = q,
  w = w
)

### set seeds for reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(2)
seeds <- sample(c(1:10000), nsim, replace = FALSE)

### set number of cores you want to use for parallel computing
no.cores <- 75



k <- (1 - q) * n


##### Empirical approach #####

results_emp <- mclapply(seq_len(nsim), function(i) {
  set.seed(seeds[i])
  cat(".")
  if(i %% 50 == 0) cat("iter:", i, "/", nsim, "\n")
  
  # Generates samples from a bivariate truncated t-distribution.
  X <- rtbt(n = n, rho = rho, df = df)

  # Estimates the marginal parameters using the Hill estimator.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  theta_emp <- cbind(theta_X1_hat, theta_X2_hat)

  # Selects a concentration parameter for the beta density mixtures using k-fold cv.
  nu <- cv_h_emp(X, nu_grid = exp(seq(1, 7, by = 0.5)), tau = q, nfolds = 5)

  # Fits the angular density.
  fit_emp <- angdensity(X, tau = q, nu = nu, grid = w)

  # Estimated angular density.
  h_emp <- function(w) {
    sapply(w,
           FUN = function(x) fit_emp$p %*% sapply(fit_emp$w,
                                                  FUN = function(u) dbeta(x, u * nu, (1 - u) * nu)))
  }
  
  # Estimated Pickands dependence function.
  A_emp <- function(w) {
    f <- function(z) sapply(fit_emp$w,
                            FUN = function(u) pbeta(z, u * nu, (1 - u) * nu)) %*% fit_emp$p

    A <- 1 - w + 2 * sapply(w, function(t) integrate(f, 0, t, abs.tol = 10^(-5))$value)
    return(A)
  }

  # Estimation of the extremal regions and calculation of the symmetric set difference probabilities.
  Q1_emp_q1 <- Q_hat(w = w, p = q_n[1], theta = theta_emp,
                  A = A_emp,
                  t = n / k, type = "or")
  p_Q1_diff_emp_q1 <- prob_set_diff(Q_hat = Q1_emp_q1, p = q_n[1],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  Q2_emp_q1 <- Q_hat(w = w, p = q_n[1], theta = theta_emp,
                  A = A_emp,
                  t = n / k, type = "and")
  p_Q2_diff_emp_q1 <- prob_set_diff(Q_hat = Q2_emp_q1, p = q_n[1],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  Q3_emp_q1 <- Q_hat(w = w, p = q_n[1], theta = theta_emp,
                  h = h_emp,
                  t = n / k, type = "density")
  p_Q3_diff_emp_q1 <- prob_set_diff(Q_hat = Q3_emp_q1, p = q_n[1],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)

  Q1_emp_q2 <- Q_hat(w = w, p = q_n[2], theta = theta_emp,
                     A = A_emp,
                     t = n / k, type = "or")
  p_Q1_diff_emp_q2 <- prob_set_diff(Q_hat = Q1_emp_q2, p = q_n[2],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  Q2_emp_q2 <- Q_hat(w = w, p = q_n[2], theta = theta_emp,
                     A = A_emp,
                     t = n / k, type = "and")
  p_Q2_diff_emp_q2 <- prob_set_diff(Q_hat = Q2_emp_q2, p = q_n[2],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  Q3_emp_q2 <- Q_hat(w = w, p = q_n[2], theta = theta_emp,
                     h = h_emp,
                     t = n / k, type = "density")
  p_Q3_diff_emp_q2 <- prob_set_diff(Q_hat = Q3_emp_q2, p = q_n[2],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)

  Q1_emp_q3 <- Q_hat(w = w, p = q_n[3], theta = theta_emp,
                     A = A_emp,
                     t = n / k, type = "or")
  p_Q1_diff_emp_q3 <- prob_set_diff(Q_hat = Q1_emp_q3, p = q_n[3],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  Q2_emp_q3 <- Q_hat(w = w, p = q_n[3], theta = theta_emp,
                     A = A_emp,
                     t = n / k, type = "and")
  p_Q2_diff_emp_q3 <- prob_set_diff(Q_hat = Q2_emp_q3, p = q_n[3],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  Q3_emp_q3 <- Q_hat(w = w, p = q_n[3], theta = theta_emp,
                     h = h_emp,
                     t = n / k, type = "density")
  p_Q3_diff_emp_q3 <- prob_set_diff(Q_hat = Q3_emp_q3, p = q_n[3],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)


  return(list(w = w,
              theta_hat = theta_emp,
              h_hat = h_emp(w),
              A_hat = A_emp(w),
              Q1_hat_q1 = Q1_emp_q1,
              Q2_hat_q1 = Q2_emp_q1,
              Q3_hat_q1 = Q3_emp_q1,
              Q1_hat_q2 = Q1_emp_q2,
              Q2_hat_q2 = Q2_emp_q2,
              Q3_hat_q2 = Q3_emp_q2,
              Q1_hat_q3 = Q1_emp_q3,
              Q2_hat_q3 = Q2_emp_q3,
              Q3_hat_q3 = Q3_emp_q3,
              set_diff_Q1_q1 = p_Q1_diff_emp_q1,
              set_diff_Q2_q1 = p_Q2_diff_emp_q1,
              set_diff_Q3_q1 = p_Q3_diff_emp_q1,
              set_diff_Q1_q2 = p_Q1_diff_emp_q2,
              set_diff_Q2_q2 = p_Q2_diff_emp_q2,
              set_diff_Q3_q2 = p_Q3_diff_emp_q2,
              set_diff_Q1_q3 = p_Q1_diff_emp_q3,
              set_diff_Q2_q3 = p_Q2_diff_emp_q3,
              set_diff_Q3_q3 = p_Q3_diff_emp_q3,
              seed = seeds[i]))

}, mc.cores = no.cores)


now <- Sys.time()
print(now)
now.date <- format(now, "%d_%B_%Y")
now.time <- format(now, "%H_%M_%S")
cidr <- getwd()
dirname <- file.path(cidr, paste("/results__", now.date, "__", now.time, sep = ""))

# creates a directory in which to store the results
if (!dir.exists(dirname)) {
  dir.create(dirname, recursive = T)
}
setwd(dirname)
filename = paste("Output_Empirical_Approach", ".rds",sep="")
saveRDS(results_emp, filename)
saveRDS(parameters_simulation_setting, "parameters_simulation_setting.rds")
setwd(cidr)

##### Likelihood-Frequentist approach #####

results_lik <- mclapply(seq_len(nsim), function(i) {
  set.seed(seeds[i])
  cat(".")
  if(i %% 50 == 0) cat("iter:", i, "/", nsim, "\n")

  # Generates samples from a bivariate truncated t-distribution.
  X <- rtbt(n = n, rho = rho, df = df)

  # Estimates the marginal parameters using the Hill estimator.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)

  theta_lik <- cbind(theta_X1_hat, theta_X2_hat)

  Y <- trans2UFrechet(X)
  
  # Selects a degree for the Bernstein polynomials using k-fold cv.
  k0 <- cv_h_lik(Y, k_grid = c(4 : 8), q = q, nfolds = 5)

  u <- c(quantile(Y[, 1], probs = q), quantile(Y[, 2], probs = q))

  fit_freq <- fExtDep.np(method = "Frequentist",
                         data = Y,
                         u = u[1],
                         mar.fit = FALSE,
                         type = "rawdata",
                         k0 = k0)

  beta <- fit_freq$Ahat$beta

  p0 <- fit_freq$p0
  p1 <- fit_freq$p1

  # Estimated angular density.
  h_lik <- function(w) {
    K <- length(beta) - 1
    eta <- K / 2 * (beta[-1] - beta[-(K + 1)] + 1 / K)
    diff_eta <- eta[-1] - eta[-K]
    h <- rowSums(sapply(seq(0, K - 2),
                        FUN = function(j) dbeta(w, j + 1, K - j - 1) * diff_eta[j + 1]))
    return(h)
  }
  
  # Estimated Pickands dependence function.
  A_lik <- function(w) {
    K <- length(beta) - 1
    A <- sapply(c(0 : K),
                FUN = function(k) sapply(1 - w,
                                         FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }

  # Estimation of the extremal regions and calculation of the symmetric set difference probabilities.
  Q1_lik_q1 <- Q_hat(w = w, p = q_n[1], theta = theta_lik,
                     A = A_lik,
                     t = n / k, type = "or")
  p_Q1_diff_lik_q1 <- prob_set_diff(Q_hat = Q1_lik_q1, p = q_n[1],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  Q2_lik_q1 <- Q_hat(w = w, p = q_n[1], theta = theta_lik,
                     A = A_lik,
                     t = n / k, type = "and")
  p_Q2_diff_lik_q1 <- prob_set_diff(Q_hat = Q2_lik_q1, p = q_n[1],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  Q3_lik_q1 <- Q_hat(w = w, p = q_n[1], theta = theta_lik,
                     h = h_lik,
                     t = n / k, type = "density")
  p_Q3_diff_lik_q1 <- prob_set_diff(Q_hat = Q3_lik_q1, p = q_n[1],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)

  Q1_lik_q2 <- Q_hat(w = w, p = q_n[2], theta = theta_lik,
                     A = A_lik,
                     t = n / k, type = "or")
  p_Q1_diff_lik_q2 <- prob_set_diff(Q_hat = Q1_lik_q2, p = q_n[2],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  Q2_lik_q2 <- Q_hat(w = w, p = q_n[2], theta = theta_lik,
                     A = A_lik,
                     t = n / k, type = "and")
  p_Q2_diff_lik_q2 <- prob_set_diff(Q_hat = Q2_lik_q2, p = q_n[2],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  Q3_lik_q2 <- Q_hat(w = w, p = q_n[2], theta = theta_lik,
                     h = h_lik,
                     t = n / k, type = "density")
  p_Q3_diff_lik_q2 <- prob_set_diff(Q_hat = Q3_lik_q2, p = q_n[2],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)

  Q1_lik_q3 <- Q_hat(w = w, p = q_n[3], theta = theta_lik,
                     A = A_lik,
                     t = n / k, type = "or")
  p_Q1_diff_lik_q3 <- prob_set_diff(Q_hat = Q1_lik_q3, p = q_n[3],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  Q2_lik_q3 <- Q_hat(w = w, p = q_n[3], theta = theta_lik,
                     A = A_lik,
                     t = n / k, type = "and")
  p_Q2_diff_lik_q3 <- prob_set_diff(Q_hat = Q2_lik_q3, p = q_n[3],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  Q3_lik_q3 <- Q_hat(w = w, p = q_n[3], theta = theta_lik,
                     h = h_lik,
                     t = n / k, type = "density")
  p_Q3_diff_lik_q3 <- prob_set_diff(Q_hat = Q3_lik_q3, p = q_n[3],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)


  return(list(w = w,
              theta_hat = theta_lik,
              h_hat = h_lik(w),
              A_hat = A_lik(w),
              Q1_hat_q1 = Q1_lik_q1,
              Q2_hat_q1 = Q2_lik_q1,
              Q3_hat_q1 = Q3_lik_q1,
              Q1_hat_q2 = Q1_lik_q2,
              Q2_hat_q2 = Q2_lik_q2,
              Q3_hat_q2 = Q3_lik_q2,
              Q1_hat_q3 = Q1_lik_q3,
              Q2_hat_q3 = Q2_lik_q3,
              Q3_hat_q3 = Q3_lik_q3,
              set_diff_Q1_q1 = p_Q1_diff_lik_q1,
              set_diff_Q2_q1 = p_Q2_diff_lik_q1,
              set_diff_Q3_q1 = p_Q3_diff_lik_q1,
              set_diff_Q1_q2 = p_Q1_diff_lik_q2,
              set_diff_Q2_q2 = p_Q2_diff_lik_q2,
              set_diff_Q3_q2 = p_Q3_diff_lik_q2,
              set_diff_Q1_q3 = p_Q1_diff_lik_q3,
              set_diff_Q2_q3 = p_Q2_diff_lik_q3,
              set_diff_Q3_q3 = p_Q3_diff_lik_q3,
              beta_hat = beta,
              p0_hat = p0,
              p1_hat = p1,
              seed = seeds[i]))

}, mc.cores = no.cores)


now <- Sys.time()
print(now)
now.date <- format(now, "%d_%B_%Y")
now.time <- format(now, "%H_%M_%S")
cidr <- getwd()
dirname <- file.path(cidr, paste("/results__", now.date, "__", now.time, sep = ""))

# creates a directory in which to store the results
if (!dir.exists(dirname)) {
  dir.create(dirname, recursive = T)
}
setwd(dirname)
filename = paste("Output_Likelihood_Frequentist_Approach", ".rds",sep="")
saveRDS(results_lik, filename)
saveRDS(parameters_simulation_setting, "parameters_simulation_setting.rds")
setwd(cidr)

##### Likelihood-Bayesian approach #####

hyperparam <- list(a.unif = 0, b.unif = 0.1,
                   mu.nbinom = 3.2, var.nbinom = 4.48)

mcmc_samples <- 2 * 10^4
burn_in <- c(1:floor(mcmc_samples/2))
ind <- seq(1, ceil(mcmc_samples/2), by = 5)

k0 <- 5

results_bay <- mclapply(seq_len(nsim), function(i) {
  set.seed(seeds[i])
  cat(".")
  if(i %% 50 == 0) cat("iter:", i, "/", nsim, "\n")

  # Generates samples from a bivariate truncated t-distribution.
  X <- rtbt(n = n, rho = rho, df = df)

  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))

  # Fits first the frequentist model for the initialization of the parameters.
  fit0 <- fExtDep.np(method = "Frequentist",
                     data = X,
                     u = u[1],
                     mar.fit = TRUE,
                     type = "rawdata", k0 = k0)

  beta0 <- fit0$Ahat$beta
  eta0 <- eta <- k0 / 2 * (diff(beta0) + 1 / k0)
  param0 <- list(eta = eta0, beta = beta0)
  pm0 <- list(p0 = eta0[1], p1 = 1 - eta0[k0])

  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X1_hat <- fGEV(X[, 1],
                       par.start = theta_X1_hat[c(2, 1, 3)],
                       method = "Frequentist",
                       u = quantile(X[, 1], probs = q),
                       optim.method = "Nelder-Mead")$est
  
  theta_X2_hat <- theta_hill(X[, 2], k)
  theta_X2_hat <- fGEV(X[, 2],
                       par.start = theta_X2_hat[c(2, 1, 3)],
                       method = "Frequentist",
                       u = quantile(X[, 2], probs = q),
                       optim.method = "Nelder-Mead")$est
  
  if (df == 1) sig <- 0.5 else sig <- 0.25
  
  cov1 <- as.matrix(rep(1, nrow(X)))
  cov2 <- as.matrix(rep(1, nrow(X)))
  
  sig10 <- sig20 <- sig
  par10 <- theta_X1_hat
  par20 <- theta_X2_hat
  
  fit_bayes <- fExtDep.np(method = "Bayesian", data = X,
                          mar.fit = TRUE, mar.prelim = FALSE,
                          cov1 = cov1, 
                          cov2 = cov2,
                          u = u,
                          par10 = par10, par20 = par20,
                          sig10 = sig10, sig20 = sig20,
                          param0 = param0,
                          k0 = k0 - 1, pm0 = pm0,
                          prior.k = "nbinom", prior.pm = "unif",
                          nk = 20, lik = TRUE, hyperparam = hyperparam,
                          nsim = mcmc_samples, type = "rawdata")

  # Discards burn_in phase and stores only every fifth posterior sample to speed up computations,
  n_iter <- length(ind)
  theta_hat1_bayesian <- fit_bayes$mar1[-burn_in, c(2, 1, 3)]
  theta_hat2_bayesian <- fit_bayes$mar2[-burn_in, c(2, 1, 3)]
  k_bayesian <- fit_bayes$k[-burn_in]
  eta_bayesian <- fit_bayes$eta[-burn_in, seq_len(21)]
  pm_bayesian <- fit_bayes$pm[-burn_in, ]
  accepted <- fit_bayes$accepted[-burn_in]
  accepted.mar1 <- fit_bayes$accepted.mar1[-burn_in]
  accepted.mar2 <- fit_bayes$accepted.mar2[-burn_in]
  m <- length(w)


  h_bay_array <- array(NA, dim = c(m, n_iter))
  A_bay_array <- array(NA, dim = c(m, n_iter))

  Q1_bay_q1_array <- array(NA, dim = c(m, 2, n_iter))
  Q1_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q1_bay_q3_array <- array(NA, dim = c(m, 2, n_iter))

  Q2_bay_q1_array <- array(NA, dim = c(m, 2, n_iter))
  Q2_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q2_bay_q3_array <- array(NA, dim = c(m, 2, n_iter))

  Q3_bay_q1_array <- array(NA, dim = c(m, 2, n_iter))
  Q3_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q3_bay_q3_array <- array(NA, dim = c(m, 2, n_iter))

  for (i in seq_len(n_iter)) {
    K <- k_bayesian[ind[i]] + 1
    eta <- eta_bayesian[ind[i] , seq_len(K)]
    diff_eta <- diff(eta)
    beta <- c(1, 1 / K * (2 * cumsum(eta) + K - c(1:K)))
    h_bay <- function(w) {
      h <- rowSums(sapply(seq(0, K - 2), FUN = function(j) dbeta(w, j + 1, K - j - 1) * diff_eta[j + 1]))
      return(h)
    }
    A_bay <- function(w) {
      A <- sapply(c(0 : K),
                  FUN = function(k) sapply(w,
                                           FUN = function(z) bernsteinb(k, K, z))) %*% beta
      return(drop(A))
    }

    theta_hat_b <- cbind(theta_hat1_bayesian[ind[i], ],
                         theta_hat2_bayesian[ind[i], ])

    h_bay_array[, i] <- h_bay(w)

    A_bay_array[, i] <- A_bay(w)

    Q1_bay_q1_array[, , i] <- Q_hat(w = w, p = q_n[1],
                                theta = theta_hat_b,
                                A = A_bay,
                                t = n / k, type = "or")
    Q1_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                theta = theta_hat_b,
                                A = A_bay,
                                t = n / k, type = "or")
    Q1_bay_q3_array[, , i] <- Q_hat(w = w, p = q_n[3],
                                theta = theta_hat_b,
                                A = A_bay,
                                t = n / k, type = "or")

    Q2_bay_q1_array[, , i] <- Q_hat(w = w, p = q_n[1],
                                theta = theta_hat_b,
                                A = A_bay,
                                t = n / k, type = "and")
    Q2_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                theta = theta_hat_b,
                                A = A_bay,
                                t = n / k, type = "and")
    Q2_bay_q3_array[, , i] <- Q_hat(w = w, p = q_n[3],
                                theta = theta_hat_b,
                                A = A_bay,
                                t = n / k, type = "and")

    Q3_bay_q1_array[, , i] <- Q_hat(w = w, p = q_n[1],
                                theta = theta_hat_b,
                                h = h_bay,
                                t = n / k, type = "density")
    Q3_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                theta = theta_hat_b,
                                h = h_bay,
                                t = n / k, type = "density")
    Q3_bay_q3_array[, , i] <- Q_hat(w = w, p = q_n[3],
                                theta = theta_hat_b,
                                h = h_bay,
                                t = n / k, type = "density")

  }

  theta_X1_bay_mean <- colMeans(theta_hat1_bayesian[ind, ])
  theta_X2_bay_mean <- colMeans(theta_hat2_bayesian[ind, ])

  p0_mean <- mean(pm_bayesian[ind, 1])
  p1_mean <- mean(pm_bayesian[ind, 2])

  k_mean <- mean(k_bayesian[ind]) + 1

  theta_bay <- cbind(theta_X1_bay_mean, theta_X2_bay_mean)

  h_bay <- rowMeans(h_bay_array, na.rm = T)

  A_bay <- rowMeans(A_bay_array, na.rm = T)

  Q1_bay_q1 <- t(apply(Q1_bay_q1_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_bay_q2 <- t(apply(Q1_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_bay_q3 <- t(apply(Q1_bay_q3_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))

  Q2_bay_q1 <- t(apply(Q2_bay_q1_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_bay_q2 <- t(apply(Q2_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_bay_q3 <- t(apply(Q2_bay_q3_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))

  Q3_bay_q1 <- t(apply(Q3_bay_q1_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_bay_q2 <- t(apply(Q3_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_bay_q3 <- t(apply(Q3_bay_q3_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))


  p_Q1_diff_bay_q1 <- prob_set_diff(Q_hat = Q1_bay_q1, p = q_n[1],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  p_Q2_diff_bay_q1 <- prob_set_diff(Q_hat = Q2_bay_q1, p = q_n[1],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  p_Q3_diff_bay_q1 <- prob_set_diff(Q_hat = Q3_bay_q1, p = q_n[1],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)

  p_Q1_diff_bay_q2 <- prob_set_diff(Q_hat = Q1_bay_q2, p = q_n[2],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  p_Q2_diff_bay_q2 <- prob_set_diff(Q_hat = Q2_bay_q2, p = q_n[2],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  p_Q3_diff_bay_q2 <- prob_set_diff(Q_hat = Q3_bay_q2, p = q_n[2],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)

  p_Q1_diff_bay_q3 <- prob_set_diff(Q_hat = Q1_bay_q3, p = q_n[3],
                                    rho = rho, df = df, type = "or",
                                    quadrature_rule = "GLe", level = 100)
  p_Q2_diff_bay_q3 <- prob_set_diff(Q_hat = Q2_bay_q3, p = q_n[3],
                                    rho = rho, df = df, type = "and",
                                    quadrature_rule = "GLe", level = 100)
  p_Q3_diff_bay_q3 <- prob_set_diff(Q_hat = Q3_bay_q3, p = q_n[3],
                                    rho = rho, df = df, type = "density",
                                    quadrature_rule = "GLe", level = 100)


  return(list(w = w,
              theta_hat = theta_bay,
              h_hat = h_bay,
              A_hat = A_bay,
              Q1_hat_q1 = Q1_bay_q1,
              Q2_hat_q1 = Q2_bay_q1,
              Q3_hat_q1 = Q3_bay_q1,
              Q1_hat_q2 = Q1_bay_q2,
              Q2_hat_q2 = Q2_bay_q2,
              Q3_hat_q2 = Q3_bay_q2,
              Q1_hat_q3 = Q1_bay_q3,
              Q2_hat_q3 = Q2_bay_q3,
              Q3_hat_q3 = Q3_bay_q3,
              set_diff_Q1_q1 = p_Q1_diff_bay_q1,
              set_diff_Q2_q1 = p_Q2_diff_bay_q1,
              set_diff_Q3_q1 = p_Q3_diff_bay_q1,
              set_diff_Q1_q2 = p_Q1_diff_bay_q2,
              set_diff_Q2_q2 = p_Q2_diff_bay_q2,
              set_diff_Q3_q2 = p_Q3_diff_bay_q2,
              set_diff_Q1_q3 = p_Q1_diff_bay_q3,
              set_diff_Q2_q3 = p_Q2_diff_bay_q3,
              set_diff_Q3_q3 = p_Q3_diff_bay_q3,
              p0_hat = p0_mean,
              p1_hat = p1_mean,
              k_hat = k_mean,
              accepted = accepted,
              accepted.mar1 = accepted.mar1,
              accepted.mar2 = accepted.mar2,
              theta_X1_hat = theta_hat1_bayesian[ind, ],
              theta_X2_hat = theta_hat2_bayesian[ind, ],
              eta_hat = eta_bayesian[ind, ],
              k_bayesian = k_bayesian[ind],
              ind = ind,
              seed = seeds[i]))

}, mc.cores = no.cores)


now <- Sys.time()
print(now)
now.date <- format(now, "%d_%B_%Y")
now.time <- format(now, "%H_%M_%S")
cidr <- getwd()
dirname <- file.path(cidr, paste("/results__", now.date, "__", now.time, sep = ""))

# creates a directory in which to store the results
if (!dir.exists(dirname)) {
  dir.create(dirname, recursive = T)
}
setwd(dirname)
filename = paste("Output_Likelihood_Bayesian_Approach", ".rds",sep="")
saveRDS(results_bay, filename)
saveRDS(parameters_simulation_setting, "parameters_simulation_setting.rds")
setwd(cidr)


