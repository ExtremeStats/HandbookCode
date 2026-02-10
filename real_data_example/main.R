# This script contains all the code used in subsection "Real data example".
# In line 15, it's possible to specify whether the models should be fitted
#  (this will be very slow) or pre-loaded. 


#### Import of packages and functions used in the analysis. ####
source("functions_real_data_analysis.R")

#### Specification of graphical parameters. ####
save_plots <- TRUE
if (save_plots) pdf("real_data_example.pdf", width = 6, height = 6)
par(mai=c(.7,1.2,.7,.2), mgp=c(3,0.8,0))

#### Specification whether models should be fitted or loaded. ####
fit_models <- FALSE
save_models <- FALSE
if (!fit_models) load("posterior_samples.RData")

#### Import of the data. ####
folder <- paste0(dirname(getwd()), "/data/")
filename <- "NegLogReturns.Rdata"
load(paste0(folder, filename))

nlr <- subset(nlr, select = c("Date",
                              "Bitcoin", "Ether", "Crypto_Index",
                              "Dow_Jones", "SPUSBI", "Gold_ETF"))



#### Estimation of the tail indices. ####
assets <- names(nlr)
assets[which(assets == "Crypto_Index")] <- "Crypto Index"
assets[which(assets == "Dow_Jones")] <- "Dow Jones"
assets[which(assets == "SPUSBI")] <- "S&P Bond Index"
assets[which(assets == "Gold_ETF")] <- "Gold"
k <- round(seq(150, 10, by = -1))
for (j in seq(2, ncol(nlr))) {
  xi_hill <- rep(NA, length(k))
  sd_hill <- rep(NA, length(k))
  for (i in seq_len(length(k))) {
    est_hill <- theta_hill(nlr[, j], k[i])
    xi_hill[i] <- est_hill[3]
    sd_hill[i] <- sqrt(xi_hill[i] / k[i])
    
  }
  plot(k, xi_hill, 
       ylim = c(-.1,1.),
       xlim = c(0, 150),
       main = assets[j], type = "l", xlab = "", ylab = "",
       col = "red", lwd = 4, cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2)
  title(xlab = "k", cex.lab = 2, mgp=c(2,1,0))
  lines(k, xi_hill - 1.64 * sd_hill, col = "red", lty = "dashed", lwd = 5)
  lines(k, xi_hill + 1.64 * sd_hill, col = "red", lty = "dashed", lwd = 5)
  abline(h=0, col="black", lty=4, lwd=5)
  legend("top", col="red", lwd=4, bty="n", cex=2, lty=1,
         legend=expression(hat(gamma)[n]~"HILL"))
}

# Based on the plots of the estimated extremal indices, we choose k = 50.
k <- 50 
# q will be the quantile used to determine the threshold for the censored likelihood.
q <- 1 - k / nrow(nlr)

#### Specification of the hyperparameters for the Bayesian approach. ####
if (fit_models) {
  prior.k <- "nbinom"
  prior.pm <- "unif"
  hyperparam <- list(a.unif = 0, b.unif = 0.5,
                     mu.nbinom = 10, var.nbinom = 100)
  k0 <- 5
  sig10 <- sig20 <- 0.1
  mcmc_samples <-  2 * 10^5
}
burn_in <- seq_len(1 * 10^5)

w <- seq(10^(-7), 1 - 10^(-7), by = 0.01)

#### Estimation of the extremal dependence between Dow Jones and Bitcoin. ####
X <- cbind(nlr$Dow_Jones, nlr$Bitcoin) * 100


t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_dow_btc <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                  mar.fit = TRUE, mar.prelim = FALSE,
                                  par10 = theta_X1_hat[c(2, 1, 3)],
                                  par20 = theta_X2_hat[c(2, 1, 3)],
                                  sig10 = sig10, sig20 = sig20,
                                  k0 = k0,
                                  prior.k = prior.k, prior.pm = prior.pm,
                                  nk = 1000, hyperparam = hyperparam,
                                  nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_dow_btc, file = "fit_bayes_dow_btc.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_dow_btc$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_dow_btc$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_dow_btc$k[-burn_in]
eta_bayesian <- fit_bayes_dow_btc$eta[-burn_in, ]
p0_bayesian <- fit_bayes_dow_btc$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_dow_btc$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

Q2_bay_array <- array(NA, dim = c(m, 2, n_iter))

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
  Q2_bay_array[, , i] <- 1 - exp(-1 / 100 * Q_hat(w = w, p = 1 / 1000,
                                        theta = theta_hat_b,
                                        A = A_bay,
                                        t = t, type = "and"))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_dow_btc <- c(extremal_coefficient_bay_mean, 
                      extremal_coefficient_bay_lb, 
                      extremal_coefficient_bay_ub)
names(ext_coef_dow_btc) <- c("mean", "lb", "ub")
point_mass_dow_btc <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                            c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_dow_btc) <- c("mean", "lb", "ub")
rownames(point_mass_dow_btc) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 3), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Dow Jones vs. Bitcoin",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

Q2_bay_mean <- t(apply(Q2_bay_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
Q2_bay_lb <- t(apply(Q2_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
Q2_bay_ub <- t(apply(Q2_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))

plot(100 * Q2_bay_mean, type = "l",
     lwd = 4, ylab = "Bitcoin", xlab = "",
     xlim = c(1, 18), ylim = c(5, 45),
     col = "forestgreen",
     main = expression(hat(Q)[n]^(2)),
     cex.lab=2, cex.axis=2, cex.main=2,
     las = 1)
title(xlab = "Dow Jones", cex.lab = 2, mgp=c(2,1,0))
lines(100 *Q2_bay_lb, col = "forestgreen", 
      lwd = 5, lty = "dashed")
lines(100 * Q2_bay_ub, col = "forestgreen", 
      lwd = 5, lty = "dashed")
legend("topright", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between Dow Jones and Ether. ####
X <- cbind(nlr$Dow_Jones, nlr$Ether) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_dow_eth <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                  mar.fit = TRUE, mar.prelim = FALSE,
                                  par10 = theta_X1_hat[c(2, 1, 3)],
                                  par20 = theta_X2_hat[c(2, 1, 3)],
                                  sig10 = sig10, sig20 = sig20,
                                  k0 = k0,
                                  prior.k = prior.k, prior.pm = prior.pm,
                                  nk = 1000, hyperparam = hyperparam,
                                  nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_dow_eth, file = "fit_bayes_dow_eth.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_dow_eth$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_dow_eth$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_dow_eth$k[-burn_in]
eta_bayesian <- fit_bayes_dow_eth$eta[-burn_in, ]
p0_bayesian <- fit_bayes_dow_eth$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_dow_eth$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_dow_eth <- c(extremal_coefficient_bay_mean, 
                      extremal_coefficient_bay_lb, 
                      extremal_coefficient_bay_ub)
names(ext_coef_dow_eth) <- c("mean", "lb", "ub")
point_mass_dow_eth <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                            c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_dow_eth) <- c("mean", "lb", "ub")
rownames(point_mass_dow_eth) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 3), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Dow Jones vs. Ether",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between Dow Jones and gold. ####
X <- cbind(nlr$Dow_Jones, nlr$Gold_ETF) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_dow_gold <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                  mar.fit = TRUE, mar.prelim = FALSE,
                                  par10 = theta_X1_hat[c(2, 1, 3)],
                                  par20 = theta_X2_hat[c(2, 1, 3)],
                                  sig10 = sig10, sig20 = sig20,
                                  k0 = k0,
                                  prior.k = prior.k, prior.pm = prior.pm,
                                  nk = 1000, hyperparam = hyperparam,
                                  nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_dow_gold, file = "fit_bayes_dow_gold.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_dow_gold$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_dow_gold$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_dow_gold$k[-burn_in]
eta_bayesian <- fit_bayes_dow_gold$eta[-burn_in, ]
p0_bayesian <- fit_bayes_dow_gold$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_dow_gold$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_dow_gold <- c(extremal_coefficient_bay_mean, 
                      extremal_coefficient_bay_lb, 
                      extremal_coefficient_bay_ub)
names(ext_coef_dow_gold) <- c("mean", "lb", "ub")
point_mass_dow_gold <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                            c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_dow_gold) <- c("mean", "lb", "ub")
rownames(point_mass_dow_gold) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 3), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Dow Jones vs. Gold",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between SP US aggregated bond index and Bitcoin. ####
X <- cbind(nlr$SPUSBI, nlr$Bitcoin) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_bonds_btc <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                   mar.fit = TRUE, mar.prelim = FALSE,
                                   par10 = theta_X1_hat[c(2, 1, 3)],
                                   par20 = theta_X2_hat[c(2, 1, 3)],
                                   sig10 = sig10, sig20 = sig20,
                                   k0 = k0,
                                   prior.k = prior.k, prior.pm = prior.pm,
                                   nk = 1000, hyperparam = hyperparam,
                                   nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_bonds_btc, file = "fit_bayes_bonds_btc.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_bonds_btc$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_bonds_btc$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_bonds_btc$k[-burn_in]
eta_bayesian <- fit_bayes_bonds_btc$eta[-burn_in, ]
p0_bayesian <- fit_bayes_bonds_btc$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_bonds_btc$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_bonds_btc <- c(extremal_coefficient_bay_mean, 
                      extremal_coefficient_bay_lb, 
                      extremal_coefficient_bay_ub)
names(ext_coef_bonds_btc) <- c("mean", "lb", "ub")
point_mass_bonds_btc <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                            c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_bonds_btc) <- c("mean", "lb", "ub")
rownames(point_mass_bonds_btc) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 3), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "S&P Bond Index vs. Bitcoin",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between SP US aggregated bond index and Ether. ####
X <- cbind(nlr$SPUSBI, nlr$Ether) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_bonds_eth <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                    mar.fit = TRUE, mar.prelim = FALSE,
                                    par10 = theta_X1_hat[c(2, 1, 3)],
                                    par20 = theta_X2_hat[c(2, 1, 3)],
                                    sig10 = sig10, sig20 = sig20,
                                    k0 = k0,
                                    prior.k = prior.k, prior.pm = prior.pm,
                                    nk = 1000, hyperparam = hyperparam,
                                    nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_bonds_eth, file = "fit_bayes_bonds_eth.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_bonds_eth$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_bonds_eth$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_bonds_eth$k[-burn_in]
eta_bayesian <- fit_bayes_bonds_eth$eta[-burn_in, ]
p0_bayesian <- fit_bayes_bonds_eth$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_bonds_eth$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_bonds_eth <- c(extremal_coefficient_bay_mean, 
                        extremal_coefficient_bay_lb, 
                        extremal_coefficient_bay_ub)
names(ext_coef_bonds_eth) <- c("mean", "lb", "ub")
point_mass_bonds_eth <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                              c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_bonds_eth) <- c("mean", "lb", "ub")
rownames(point_mass_bonds_eth) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 3), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "S&P Bond Index vs. Ether",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between SP US aggregated bond index and gold. ####
X <- cbind(nlr$SPUSBI, nlr$Gold_ETF) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_bonds_gold <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                  mar.fit = TRUE, mar.prelim = FALSE,
                                  par10 = theta_X1_hat[c(2, 1, 3)],
                                  par20 = theta_X2_hat[c(2, 1, 3)],
                                  sig10 = sig10, sig20 = sig20,
                                  k0 = k0,
                                  prior.k = prior.k, prior.pm = prior.pm,
                                  nk = 1000, hyperparam = hyperparam,
                                  nsim = mcmc_samples, type = "rawdata")
  
  if (save_models) {
    saveRDS(fit_bayes_bonds_gold, file = "fit_bayes_bonds_gold.rds")
  }
}


# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_bonds_gold$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_bonds_gold$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_bonds_gold$k[-burn_in]
eta_bayesian <- fit_bayes_bonds_gold$eta[-burn_in, ]
p0_bayesian <- fit_bayes_bonds_gold$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_bonds_gold$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

Q2_bay_array <- array(NA, dim = c(m, 2, n_iter))

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
  Q2_bay_array[, , i] <- 1 - exp(-1 / 100 * Q_hat(w = w, p = 1 / 1000,
                                        theta = theta_hat_b,
                                        A = A_bay,
                                        t = t, type = "and"))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_bonds_gold <- c(extremal_coefficient_bay_mean, 
                        extremal_coefficient_bay_lb, 
                        extremal_coefficient_bay_ub)
names(ext_coef_bonds_gold) <- c("mean", "lb", "ub")
point_mass_bonds_gold <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                              c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_bonds_gold) <- c("mean", "lb", "ub")
rownames(point_mass_bonds_gold) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 3), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "S&P Bond Index vs. Gold",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

Q2_bay_mean <- t(apply(Q2_bay_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
Q2_bay_lb <- t(apply(Q2_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
Q2_bay_ub <- t(apply(Q2_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))

plot(100 * Q2_bay_mean, type = "l",
     lwd = 4, ylab = "Gold", xlab = "",
     xlim = c(0.1, 3), ylim = c(1.5, 8),
     col = "forestgreen",
     main = expression(hat(Q)[n]^(2)),
     cex.lab=2, cex.axis=2, cex.main=2,
     las = 1)
title(xlab = "S&P Bond Index", cex.lab = 2, mgp=c(2,1,0))
lines(100 *Q2_bay_lb, col = "forestgreen", 
      lwd = 5, lty = "dashed")
lines(100 * Q2_bay_ub, col = "forestgreen", 
      lwd = 5, lty = "dashed")
legend("topright", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between Crypto Index and Bitcoin. ####
X <- cbind(nlr$Crypto_Index, nlr$Bitcoin) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_crypto_btc <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                    mar.fit = TRUE, mar.prelim = FALSE,
                                    par10 = theta_X1_hat[c(2, 1, 3)],
                                    par20 = theta_X2_hat[c(2, 1, 3)],
                                    sig10 = sig10, sig20 = sig20,
                                    k0 = k0,
                                    prior.k = prior.k, prior.pm = prior.pm,
                                    nk = 1000, hyperparam = hyperparam,
                                    nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_crypto_btc, file = "fit_bayes_crypto_btc.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_crypto_btc$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_crypto_btc$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_crypto_btc$k[-burn_in]
eta_bayesian <- fit_bayes_crypto_btc$eta[-burn_in, ]
p0_bayesian <- fit_bayes_crypto_btc$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_crypto_btc$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_crypto_btc <- c(extremal_coefficient_bay_mean, 
                        extremal_coefficient_bay_lb, 
                        extremal_coefficient_bay_ub)
names(ext_coef_crypto_btc) <- c("mean", "lb", "ub")
point_mass_crypto_btc <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                              c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_crypto_btc) <- c("mean", "lb", "ub")
rownames(point_mass_crypto_btc) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 2), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Crypto Index vs. Bitcoin",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between Crypto Index and Ether. ####
X <- cbind(nlr$Crypto_Index, nlr$Ether) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_crypto_eth <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                     mar.fit = TRUE, mar.prelim = FALSE,
                                     par10 = theta_X1_hat[c(2, 1, 3)],
                                     par20 = theta_X2_hat[c(2, 1, 3)],
                                     sig10 = sig10, sig20 = sig20,
                                     k0 = k0,
                                     prior.k = prior.k, prior.pm = prior.pm,
                                     nk = 1000, hyperparam = hyperparam,
                                     nsim = mcmc_samples, type = "rawdata")
  if (save_models) {
    saveRDS(fit_bayes_crypto_eth, file = "fit_bayes_crypto_eth.rds")
  }
}


# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_crypto_eth$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_crypto_eth$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_crypto_eth$k[-burn_in]
eta_bayesian <- fit_bayes_crypto_eth$eta[-burn_in, ]
p0_bayesian <- fit_bayes_crypto_eth$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_crypto_eth$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_crypto_eth <- c(extremal_coefficient_bay_mean, 
                         extremal_coefficient_bay_lb, 
                         extremal_coefficient_bay_ub)
names(ext_coef_crypto_eth) <- c("mean", "lb", "ub")
point_mass_crypto_eth <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                               c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_crypto_eth) <- c("mean", "lb", "ub")
rownames(point_mass_crypto_eth) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 2), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Crypto Index vs. Ether",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the extremal dependence between Bitcoin and Ether. ####
X <- cbind(nlr$Bitcoin, nlr$Ether) * 100

t <- nrow(X) / k

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_btc_eth <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                  mar.fit = TRUE, mar.prelim = FALSE,
                                  par10 = theta_X1_hat[c(2, 1, 3)],
                                  par20 = theta_X2_hat[c(2, 1, 3)],
                                  sig10 = sig10, sig20 = sig20,
                                  k0 = k0,
                                  prior.k = prior.k, prior.pm = prior.pm,
                                  nk = 1000, hyperparam = hyperparam,
                                  nsim = mcmc_samples, type = "rawdata")
  
  if (save_models) {
    saveRDS(fit_bayes_btc_eth, file = "fit_bayes_btc_eth.rds")
  }
}


# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_btc_eth$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_btc_eth$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_btc_eth$k[-burn_in]
eta_bayesian <- fit_bayes_btc_eth$eta[-burn_in, ]
p0_bayesian <- fit_bayes_btc_eth$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_btc_eth$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

Q1_bay_array <- array(NA, dim = c(m, 2, n_iter))
Q2_bay_array <- array(NA, dim = c(m, 2, n_iter))

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
  Q1_bay_array[, , i] <- 1 - exp(-1 / 100 * Q_hat(w = w, p = 1 / 1000,
                                        theta = theta_hat_b,
                                        A = A_bay,
                                        t = t, type = "or"))
  
  Q2_bay_array[, , i] <- 1 - exp(-1 / 100 * Q_hat(w = w, p = 1 / 1000,
                                        theta = theta_hat_b,
                                        A = A_bay,
                                        t = t, type = "and"))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_btc_eth <- c(extremal_coefficient_bay_mean, 
                         extremal_coefficient_bay_lb, 
                         extremal_coefficient_bay_ub)
names(ext_coef_btc_eth) <- c("mean", "lb", "ub")
point_mass_btc_eth <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                               c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_btc_eth) <- c("mean", "lb", "ub")
rownames(point_mass_btc_eth) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 2), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Bitcoin vs. Ether",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

Q1_bay_mean <- t(apply(Q1_bay_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
Q1_bay_lb <- t(apply(Q1_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
Q1_bay_lb <- rbind(Q1_bay_lb, 
                   c(10^6, min(Q1_bay_lb[, 2]) - 10^(-16)),
                   c(min(Q1_bay_lb[, 1]) - 10^(-16), 10^6))
Q1_bay_lb <- Q1_bay_lb[order(Q1_bay_lb[, 1]), ]
Q1_bay_ub <- t(apply(Q1_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))

plot(100 * Q1_bay_mean, type = "l",
     lwd = 4, ylab = "Ether", xlab = "",
     xlim = c(15, 50), ylim = c(15, 55),
     col = "forestgreen",
     main = expression(hat(Q)[n]^(1)),
     cex.lab=2, cex.axis=2, cex.main=2,
     las = 1)
title(xlab = "Bitcoin", cex.lab = 2, mgp=c(2,1,0))
lines(100 *Q1_bay_lb, col = "forestgreen", 
      lwd = 5, lty = "dashed")
lines(100 * Q1_bay_ub, col = "forestgreen", 
      lwd = 5, lty = "dashed")
legend("topright", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)


Q2_bay_mean <- t(apply(Q2_bay_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
Q2_bay_lb <- t(apply(Q2_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
Q2_bay_ub <- t(apply(Q2_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))

plot(100 * Q2_bay_mean, type = "l",
     lwd = 4, ylab = "Ether", xlab = "",
     xlim = c(5, 30), ylim = c(10, 40),
     col = "forestgreen",
     main = expression(hat(Q)[n]^(2)),
     cex.lab=2, cex.axis=2, cex.main=2,
     las = 1)
title(xlab = "Bitcoin", cex.lab = 2, mgp=c(2,1,0))
lines(100 *Q2_bay_lb, col = "forestgreen", 
      lwd = 5, lty = "dashed")
lines(100 * Q2_bay_ub, col = "forestgreen", 
      lwd = 5, lty = "dashed")
legend("topright", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### Estimation of the conditional extremal dependence between Bitcoin and Ether. ####
X <- cbind(nlr$Bitcoin, nlr$Ether) * 100
X <- X[which((X[, 1] >= 0) & (X[, 2] >= 0)), ]

t <- nrow(X) / k
q <- 1 - 1 / t

if (fit_models) {
  # Threshold for the censored likelihood.
  u <- c(quantile(X[, 1], probs = q), quantile(X[, 2], probs = q))
  # We use the Hill estimates as initial values for the marginal parameters go get a warm start.
  theta_X1_hat <- theta_hill(X[, 1], k)
  theta_X2_hat <- theta_hill(X[, 2], k)
  
  set.seed(1)
  fit_bayes_btc_eth_cond <- fExtDep.np(method = "Bayesian", data = X, u = u,
                                  mar.fit = TRUE, mar.prelim = FALSE,
                                  par10 = theta_X1_hat[c(2, 1, 3)],
                                  par20 = theta_X2_hat[c(2, 1, 3)],
                                  sig10 = sig10, sig20 = sig20,
                                  k0 = k0,
                                  prior.k = prior.k, prior.pm = prior.pm,
                                  nk = 1000, hyperparam = hyperparam,
                                  nsim = mcmc_samples, type = "rawdata")
  
  if (save_models) {
    saveRDS(fit_bayes_btc_eth_cond, file = "fit_bayes_btc_eth_cond.rds")
  }
}



# We discard the first posterior samples and keep the remaining ones.
theta_X1_hat_bayesian <- fit_bayes_btc_eth_cond$mar1[-burn_in, c(2, 1, 3)]
theta_X2_hat_bayesian <- fit_bayes_btc_eth_cond$mar2[-burn_in, c(2, 1, 3)]
k_bayesian <- fit_bayes_btc_eth_cond$k[-burn_in]
eta_bayesian <- fit_bayes_btc_eth_cond$eta[-burn_in, ]
p0_bayesian <- fit_bayes_btc_eth_cond$pm[-burn_in, 1]
p1_bayesian <- fit_bayes_btc_eth_cond$pm[-burn_in, 2]

# To decrease running time we'll only consider every 10th posterior sample.
ind <- seq(1, length(k_bayesian), by = 10)
n_iter <- length(ind)

m <- length(w)
h_bay_array <- array(NA, dim = c(m, n_iter))
A_bay_array <- array(NA, dim = c(m, n_iter))
extremal_coefficient_bay_array <- rep(NA, n_iter)

Q3_bay_array <- array(NA, dim = c(m, 2, n_iter))

# Here, we estimate the angular density, the Pickands dependence function,
# the marginal parameters and the extremal coefficient for every posterior sample
# considered. 
# A more convenient way to obtain them is to use the function summmary_ExtDep.

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
                FUN = function(k) sapply(w, FUN = function(z) bernsteinb(k, K, z))) %*% beta
    return(drop(A))
  }
  
  theta_hat_b <- cbind(theta_X1_hat_bayesian[ind[i], ],
                       theta_X2_hat_bayesian[ind[i], ])
  
  h_bay_array[, i] <- h_bay(w)
  
  A_bay_array[, i] <- A_bay(w)
  f_bay <- function(w) {
    z <- sapply(w, FUN = function(t) max(t, 1 - t))
    z * h_bay(w)
  }
  extremal_coefficient_bay_array[i] <- 
    2 * (integrate(f_bay, 10^-16, 1 - 10^-16)$value + eta[1] + (1 - eta[K]))
  
  Q3_bay_array[, , i] <- 1 - exp(-1 / 100 * Q_hat(w = w, p = 1 / 1000,
                                        theta = theta_hat_b,
                                        h = h_bay,
                                        t = t, type = "density"))
  
}

h_bay_mean <- rowMeans(h_bay_array, na.rm = T)
h_bay_lb <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.025))
h_bay_ub <- apply(h_bay_array, 1, 
                  FUN = function(z) quantile(z, probs = 0.975))
p0_bay_mean <- mean(p0_bayesian[ind])
p0_bay_lb <- quantile(p0_bayesian[ind], probs = 0.025)
p0_bay_ub <- quantile(p0_bayesian[ind], probs = 0.975)
p1_bay_mean <- mean(p1_bayesian[ind])
p1_bay_lb <- quantile(p1_bayesian[ind], probs = 0.025)
p1_bay_ub <- quantile(p1_bayesian[ind], probs = 0.975)

extremal_coefficient_bay_mean <- mean(extremal_coefficient_bay_array)
extremal_coefficient_bay_lb <- quantile(extremal_coefficient_bay_array, probs = 0.025)
extremal_coefficient_bay_ub <- quantile(extremal_coefficient_bay_array, probs = 0.975)

ext_coef_btc_eth_cond <- c(extremal_coefficient_bay_mean, 
                      extremal_coefficient_bay_lb, 
                      extremal_coefficient_bay_ub)
names(ext_coef_btc_eth_cond) <- c("mean", "lb", "ub")
point_mass_btc_eth_cond <- rbind(c(p0_bay_mean, p0_bay_lb, p0_bay_ub),
                            c(p1_bay_mean, p1_bay_lb, p1_bay_ub))
colnames(point_mass_btc_eth_cond) <- c("mean", "lb", "ub")
rownames(point_mass_btc_eth_cond) <- c("p0", "p1")

plot(w, h_bay_mean, type = "l", ylim = c(0, 2), xlim = c(0, 1),
     ylab = expression(hat(h)(w)), xlab = "", main =  "Bitcoin vs. Ether (cond. distr.)",
     lwd = 4, cex.lab=2, cex.axis=2, cex.main =2, cex.sub =2,
     las = 1, col = "forestgreen")
title(xlab = "w", cex.lab = 2, mgp=c(2,1,0))
lines(w, h_bay_lb, lty = "dashed", lwd = 5, col = "forestgreen")
lines(w, h_bay_ub, lty = "dashed", lwd = 5, col = "forestgreen")
points(c(0, 1), c(p0_bay_mean, p1_bay_mean), pch = 16, col = "forestgreen", cex = 2)
arrows(x0 = 0, y0 = p0_bay_lb, x1 = 0, y1 = p0_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
arrows(x0 = 1, y0 = p1_bay_lb, x1 = 1, y1 = p1_bay_ub, 
       code = 3, angle = 90, length = 0.1, lwd = 4, col = "forestgreen")
legend("top", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

Q3_bay_mean <- t(apply(Q3_bay_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
Q3_bay_lb <- t(apply(Q3_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
Q3_bay_ub <- t(apply(Q3_bay_array, 1, FUN = function(x) 
  apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))

plot(100 * Q3_bay_mean, type = "l",
     lwd = 4, ylab = "Ether", xlab = "",
     xlim = c(10, 50), ylim = c(10, 80),
     col = "forestgreen",
     main = expression(hat(Q)[n]^(3)),
     cex.lab=2, cex.axis=2, cex.main=2,
     las = 1)
title(xlab = "Bitcoin", cex.lab = 2, mgp=c(2,1,0))
lines(100 *Q3_bay_lb, col = "forestgreen", 
      lwd = 5, lty = "dashed")
lines(100 * Q3_bay_ub, col = "forestgreen", 
      lwd = 5, lty = "dashed")
legend("topright", 
       legend = c("Lik.-Bayes."),
       lty = 1, col = c("forestgreen"),
       bty = "n", lwd = 4, cex=2, seg.len=1)

#### save plots ####
if (save_plots) dev.off()



#### diagnostic plots ####
if (save_plots) pdf("mcmc_real_data_example.pdf", width = 6, height = 6)
par(mai=c(.7,1.2,.7,.2), mgp=c(3,0.8,0))
plot(fit_bayes_dow_btc$mar1[, 3], 
     ylim = c(min(fit_bayes_dow_btc$mar1[, 3], fit_bayes_dow_btc$mar2[, 3]),
              max(fit_bayes_dow_btc$mar1[, 3], fit_bayes_dow_btc$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Dow Jones vs. Bitcoin", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_dow_btc$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Dow Jones"), expression(hat(xi)~"Bitcoin")))


plot(fit_bayes_dow_eth$mar1[, 3], 
     ylim = c(min(fit_bayes_dow_eth$mar1[, 3], fit_bayes_dow_eth$mar2[, 3]),
              max(fit_bayes_dow_eth$mar1[, 3], fit_bayes_dow_eth$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Dow Jones vs. Ether", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_dow_eth$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Dow Jones"), expression(hat(xi)~"Ether")))

plot(fit_bayes_dow_gold$mar1[, 3], 
     ylim = c(min(fit_bayes_dow_gold$mar1[, 3], fit_bayes_dow_gold$mar2[, 3]),
              max(fit_bayes_dow_gold$mar1[, 3], fit_bayes_dow_gold$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Dow Jones vs. Gold", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_dow_gold$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Dow Jones"), expression(hat(xi)~"Gold")))

plot(fit_bayes_bonds_btc$mar1[, 3], 
     ylim = c(min(fit_bayes_bonds_btc$mar1[, 3], fit_bayes_bonds_btc$mar2[, 3]),
              max(fit_bayes_bonds_btc$mar1[, 3], fit_bayes_bonds_btc$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Bond Index vs. Bitcoin", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_bonds_btc$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Bond Index"), expression(hat(xi)~"Bitcoin")))

plot(fit_bayes_bonds_eth$mar1[, 3], 
     ylim = c(min(fit_bayes_bonds_eth$mar1[, 3], fit_bayes_bonds_eth$mar2[, 3]),
              max(fit_bayes_bonds_eth$mar1[, 3], fit_bayes_bonds_eth$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Bond Index vs. Gold", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_bonds_eth$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Bond Index"), expression(hat(xi)~"Ether")))

plot(fit_bayes_bonds_gold$mar1[, 3], 
     ylim = c(min(fit_bayes_bonds_gold$mar1[, 3], fit_bayes_bonds_gold$mar2[, 3]),
              max(fit_bayes_bonds_gold$mar1[, 3], fit_bayes_bonds_gold$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Bond Index vs. Gold", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_bonds_gold$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Bond Index"), expression(hat(xi)~"Gold")))

plot(fit_bayes_crypto_btc$mar1[, 3], 
     ylim = c(min(fit_bayes_crypto_btc$mar1[, 3], fit_bayes_crypto_btc$mar2[, 3]),
              max(fit_bayes_crypto_btc$mar1[, 3], fit_bayes_crypto_btc$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Crypto Index vs. Bitcoin", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_crypto_btc$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Crypto Index"), expression(hat(xi)~"Bitcoin")))

plot(fit_bayes_crypto_eth$mar1[, 3], 
     ylim = c(min(fit_bayes_crypto_eth$mar1[, 3], fit_bayes_crypto_eth$mar2[, 3]),
              max(fit_bayes_crypto_eth$mar1[, 3], fit_bayes_crypto_eth$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Crypto Index vs. Ether", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_crypto_eth$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Crypto Index"), expression(hat(xi)~"Ether")))

plot(fit_bayes_btc_eth$mar1[, 3], 
     ylim = c(min(fit_bayes_btc_eth$mar1[, 3], fit_bayes_btc_eth$mar2[, 3]),
              max(fit_bayes_btc_eth$mar1[, 3], fit_bayes_btc_eth$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Bitcoin vs. Ether", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_btc_eth$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Bitcoin"), expression(hat(xi)~"Ether")))

plot(fit_bayes_btc_eth_cond$mar1[, 3], 
     ylim = c(min(fit_bayes_btc_eth_cond$mar1[, 3], fit_bayes_btc_eth_cond$mar2[, 3]),
              max(fit_bayes_btc_eth_cond$mar1[, 3], fit_bayes_btc_eth_cond$mar2[, 3]) + 0.2),
     xlim = c(1, max(seq_len(mcmc_samples))),
     main = "Bitcoin vs. Ether (cond.)", type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2, cex.lab=2, cex.axis=2, cex.main = 2,
     las = 1, cex.sub = 2)
lines(fit_bayes_btc_eth_cond$mar2[, 3], col = "blue", lwd = 2)
abline(v = max(burn_in), lwd = 3, lty = "dashed")
legend("topleft", col=c("red", "blue"), lwd=2, bty="n", cex=1.5, lty=1,
       legend=c(expression(hat(xi)~"Bitcoin"), expression(hat(xi)~"Ether")))

dev.off()