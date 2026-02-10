# This script creates the plots from the subsection "Simulation_study" 
#  after running "main.R".

source("functions_simulation_study.R")

##### simulation setting 1 #####
sim_setting <- "xxx/"

folder_emp <- "xxx/"
filename_output_emp <- "Output_Empirical_Approach.rds"
results_emp1 <- readRDS(paste0(sim_setting, folder_emp, filename_output_emp))


folder_lik <- "xxx/"
filename_output_lik <- "Output_Likelihood_Frequentist_Approach.rds"
results_lik1 <- readRDS(paste0(sim_setting, folder_lik, filename_output_lik))

folder_bay <- "xxx/"
filename_output_bay <- "Output_Likelihood_Bayesian_Approach.rds"
results_bay1 <- readRDS(paste0(sim_setting, folder_bay, filename_output_bay))


filename_param <- "parameters_simulation_setting.rds"
params_emp1 <- readRDS(paste0(sim_setting, folder_emp, filename_param))
params_lik1 <- readRDS(paste0(sim_setting, folder_lik, filename_param))
params_bay1 <- readRDS(paste0(sim_setting, folder_bay, filename_param))


check_params_1 <- all(sapply(seq_len(length(params_emp1)), 
                             FUN = function(i) identical(params_emp1[[i]], params_lik1[[i]])))
check_params_2 <- all(sapply(seq_len(length(params_emp1)), 
                             FUN = function(i) identical(params_emp1[[i]], params_bay1[[i]])))

if (!check_params_1 | !check_params_2) warning("Simulation parameters do not match!")

params1 <- params_emp1

#### Plots full run ####
plot_full1 <- function() {
  nsim <- params1$nsim
  n <- params1$n
  rho <- params1$rho
  df <- params1$df
  q_n <- params1$q_n
  q <- params1$q
  w <- params1$w
  
  k <- (1 - q) * n
  
  ind_stuck <- rep(FALSE, nsim)
  for (i in seq_len(nsim)) {
    accepted <- results_bay1[[i]]$accepted
    acc_diff <- c(1, accepted) - c(accepted, 1)
    pos_0 <- which(acc_diff == -1)
    pos_1 <- which(acc_diff == 1)
    pos_diff <- pos_0 - pos_1
    if (any(pos_diff >= 300)) ind_stuck[i] <- TRUE
  }
  
  ind_keep <- which(!ind_stuck)[seq_len(1000)]
  
  Q1_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "or", m = 100)$Q
  
  Q1_hat_emp_q2 <- sapply(ind_keep, FUN = function(i) results_emp1[[i]]$Q1_hat_q2, simplify = "array")
  Q1_hat_lik_q2 <- sapply(ind_keep, FUN = function(i) results_lik1[[i]]$Q1_hat_q2, simplify = "array")
  Q1_hat_bay_q2 <- sapply(ind_keep, FUN = function(i) results_bay1[[i]]$Q1_hat_q2, simplify = "array")
  
  Q1_emp_mean_q2 <- t(apply(Q1_hat_emp_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_emp_mean_q2 <- rbind(Q1_emp_mean_q2, 
                          c(min(Q1_emp_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_emp_mean_q2[, 2])))
  Q1_emp_mean_q2 <- Q1_emp_mean_q2[order(Q1_emp_mean_q2[, 1]), ]
  Q1_lik_mean_q2 <- t(apply(Q1_hat_lik_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_lik_mean_q2 <- rbind(Q1_lik_mean_q2, 
                          c(min(Q1_lik_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_lik_mean_q2[, 2])))
  Q1_lik_mean_q2 <- Q1_lik_mean_q2[order(Q1_lik_mean_q2[, 1]), ]
  Q1_bay_mean_q2 <- t(apply(Q1_hat_bay_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_bay_mean_q2 <- rbind(Q1_bay_mean_q2, 
                          c(min(Q1_bay_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_mean_q2[, 2])))
  Q1_bay_mean_q2 <- Q1_bay_mean_q2[order(Q1_emp_mean_q2[, 1]), ]
  Q1_emp_lb_q2 <- t(apply(Q1_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q1_emp_lb_q2 <- rbind(Q1_emp_lb_q2, 
                          c(min(Q1_emp_lb_q2[, 1]), 10^10),
                          c(10^10, min(Q1_emp_lb_q2[, 2])))
  Q1_emp_lb_q2 <- Q1_emp_lb_q2[order(Q1_emp_lb_q2[, 1]), ]
  Q1_emp_ub_q2 <- t(apply(Q1_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q1_emp_ub_q2 <- rbind(Q1_emp_ub_q2, 
                          c(min(Q1_emp_ub_q2[, 1]), 10^10),
                          c(10^10, min(Q1_emp_ub_q2[, 2])))
  Q1_emp_ub_q2 <- Q1_emp_ub_q2[order(Q1_emp_ub_q2[, 1]), ]
  Q1_lik_lb_q2 <- t(apply(Q1_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q1_lik_lb_q2 <- rbind(Q1_lik_lb_q2, 
                          c(min(Q1_lik_lb_q2[, 1]), 10^10),
                          c(10^10, min(Q1_lik_lb_q2[, 2])))
  Q1_lik_lb_q2 <- Q1_lik_lb_q2[order(Q1_lik_lb_q2[, 1]), ]
  Q1_lik_ub_q2 <- t(apply(Q1_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q1_lik_ub_q2 <- rbind(Q1_lik_ub_q2, 
                          c(min(Q1_lik_ub_q2[, 1]), 10^10),
                          c(10^10, min(Q1_lik_ub_q2[, 2])))
  Q1_lik_ub_q2 <- Q1_lik_ub_q2[order(Q1_lik_ub_q2[, 1]), ]
  Q1_bay_lb_q2 <- t(apply(Q1_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q1_bay_lb_q2 <- rbind(Q1_bay_lb_q2, 
                          c(min(Q1_bay_lb_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_lb_q2[, 2])))
  Q1_bay_lb_q2 <- Q1_bay_lb_q2[order(Q1_bay_lb_q2[, 1]), ]
  Q1_bay_ub_q2 <- t(apply(Q1_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q1_bay_ub_q2 <- rbind(Q1_bay_ub_q2, 
                          c(min(Q1_bay_ub_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_ub_q2[, 2])))
  Q1_bay_ub_q2 <- Q1_bay_ub_q2[order(Q1_bay_ub_q2[, 1]), ]
  
  Q2_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "and", m = 100)$Q
  
  Q2_hat_emp_q2 <- sapply(ind_keep, FUN = function(i) results_emp1[[i]]$Q2_hat_q2, simplify = "array")
  Q2_hat_lik_q2 <- sapply(ind_keep, FUN = function(i) results_lik1[[i]]$Q2_hat_q2, simplify = "array")
  Q2_hat_bay_q2 <- sapply(ind_keep, FUN = function(i) results_bay1[[i]]$Q2_hat_q2, simplify = "array")
  
  Q2_emp_mean_q2 <- t(apply(Q2_hat_emp_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_lik_mean_q2 <- t(apply(Q2_hat_lik_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_bay_mean_q2 <- t(apply(Q2_hat_bay_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_emp_lb_q2 <- t(apply(Q2_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q2_emp_ub_q2 <- t(apply(Q2_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q2_lik_lb_q2 <- t(apply(Q2_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q2_lik_ub_q2 <- t(apply(Q2_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q2_bay_lb_q2 <- t(apply(Q2_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q2_bay_ub_q2 <- t(apply(Q2_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  
  Q3_true_q2 <- ellipse(sigma = matrix(c(1, rho, rho, 1), nrow = 2), 
                        df = df, prob = 1 - q_n[2],
                        npoints = 100, pos = TRUE)
  
  Q3_hat_emp_q2 <- sapply(ind_keep, FUN = function(i) results_emp1[[i]]$Q3_hat_q2, simplify = "array")
  Q3_hat_lik_q2 <- sapply(ind_keep, FUN = function(i) results_lik1[[i]]$Q3_hat_q2, simplify = "array")
  Q3_hat_bay_q2 <- sapply(ind_keep, FUN = function(i) results_bay1[[i]]$Q3_hat_q2, simplify = "array")
  
  Q3_emp_mean_q2 <- t(apply(Q3_hat_emp_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_lik_mean_q2 <- t(apply(Q3_hat_lik_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_bay_mean_q2 <- t(apply(Q3_hat_bay_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_emp_lb_q2 <- t(apply(Q3_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q3_emp_ub_q2 <- t(apply(Q3_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q3_lik_lb_q2 <- t(apply(Q3_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q3_lik_ub_q2 <- t(apply(Q3_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q3_bay_lb_q2 <- t(apply(Q3_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q3_bay_ub_q2 <- t(apply(Q3_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  
  plot(Q1_true_q2,
       xlim = c(0, 8),
       ylim = c(0, 8), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(1)))
  lines(Q1_emp_mean_q2 / 1000, col = "red", 
        lwd = 4)
  lines(Q1_lik_mean_q2/ 1000, col = "blue", 
        lwd = 4)
  lines(Q1_bay_mean_q2/ 1000, col = "forestgreen", 
        lwd = 4)
  lines(Q1_emp_lb_q2/ 1000, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q1_emp_ub_q2/ 1000, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q1_lik_lb_q2/ 1000, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q1_lik_ub_q2/ 1000, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q1_bay_lb_q2/ 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q1_bay_ub_q2/ 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q2_true_q2 / 1000, 
       xlim = c(0, 2.5),
       ylim = c(0, 2.5), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(2)))
  lines(Q2_emp_mean_q2 / 1000, col = "red", 
        lwd = 4)
  lines(Q2_lik_mean_q2 / 1000, col = "blue", 
        lwd = 4)
  lines(Q2_bay_mean_q2 / 1000, col = "forestgreen", 
        lwd = 4)
  lines(Q2_emp_lb_q2 / 1000, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q2_emp_ub_q2 / 1000, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q2_lik_lb_q2 / 1000, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q2_lik_ub_q2 / 1000, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q2_bay_lb_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q2_bay_ub_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q3_true_q2 / 1000, 
       xlim = c(0, 5),
       ylim = c(0, 5), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(3)))
  lines(Q3_emp_mean_q2 / 1000, col = "red", 
        lwd = 4)
  lines(Q3_lik_mean_q2 / 1000, col = "blue", 
        lwd = 4)
  lines(Q3_bay_mean_q2 / 1000, col = "forestgreen", 
        lwd = 4)
  lines(Q3_emp_lb_q2 / 1000, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q3_emp_ub_q2 / 1000, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q3_lik_lb_q2 / 1000, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q3_lik_ub_q2 / 1000, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q3_bay_lb_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q3_bay_ub_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
}

#### Plots single run ####
plot_single_run1 <- function(j = 1) {
  nsim <- params1$nsim
  n <- params1$n
  rho <- params1$rho
  df <- params1$df
  q_n <- params1$q_n
  q <- params1$q
  w <- params1$w
  
  k <- (1 - q) * n
  
  ind_stuck <- rep(FALSE, nsim)
  for (i in seq_len(nsim)) {
    accepted <- results_bay1[[i]]$accepted
    acc_diff <- c(1, accepted) - c(accepted, 1)
    pos_0 <- which(acc_diff == -1)
    pos_1 <- which(acc_diff == 1)
    pos_diff <- pos_0 - pos_1
    if (any(pos_diff >= 300)) ind_stuck[i] <- TRUE
  }
  
  ind_keep <- which(!ind_stuck)[seq_len(1000)]
  
  Q1_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "or", m = 100)$Q
  Q2_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "and", m = 100)$Q
  Q3_true_q2 <- ellipse(sigma = matrix(c(1, rho, rho, 1), nrow = 2), 
                        df = df, prob = 1 - q_n[2],
                        npoints = 100, pos = TRUE)
  
  Q1_emp_mean_q2 <- results_emp1[[ind_keep[j]]]$Q1_hat_q2
  Q1_emp_mean_q2 <- rbind(Q1_emp_mean_q2, 
                          c(min(Q1_emp_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_emp_mean_q2[, 2])))
  Q1_emp_mean_q2 <- Q1_emp_mean_q2[order(Q1_emp_mean_q2[, 1]), ]
  Q1_lik_mean_q2 <- results_lik1[[ind_keep[j]]]$Q1_hat_q2
  Q1_lik_mean_q2 <- rbind(Q1_lik_mean_q2, 
                          c(min(Q1_lik_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_lik_mean_q2[, 2])))
  Q1_lik_mean_q2 <- Q1_lik_mean_q2[order(Q1_lik_mean_q2[, 1]), ]
  Q2_emp_mean_q2 <- results_emp1[[ind_keep[j]]]$Q2_hat_q2
  Q2_lik_mean_q2 <- results_lik1[[ind_keep[j]]]$Q2_hat_q2
  Q3_emp_mean_q2 <- results_emp1[[ind_keep[j]]]$Q3_hat_q2
  Q3_lik_mean_q2 <- results_lik1[[ind_keep[j]]]$Q3_hat_q2
  
  
  theta_hat1_bayesian <- results_bay1[[ind_keep[j]]]$theta_X1_hat
  theta_hat2_bayesian <- results_bay1[[ind_keep[j]]]$theta_X2_hat
  k_bayesian <- results_bay1[[ind_keep[j]]]$k_bayesian
  eta_bayesian <- results_bay1[[ind_keep[j]]]$eta_hat
  
  m <- length(results_bay1[[ind_keep[j]]]$w)
  n_iter <- length(k_bayesian)
  Q1_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q2_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q3_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  
  for (i in seq_len(n_iter)) {
    K <- k_bayesian[i] + 1
    eta <- eta_bayesian[i , seq_len(K)]
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
    
    theta_hat_b <- cbind(theta_hat1_bayesian[i, ],
                         theta_hat2_bayesian[i, ])
    
    Q1_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                    theta = theta_hat_b,
                                    A = A_bay,
                                    t = n / k, type = "or")
    Q2_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                    theta = theta_hat_b,
                                    A = A_bay,
                                    t = n / k, type = "and")
    Q3_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                    theta = theta_hat_b,
                                    h = h_bay,
                                    t = n / k, type = "density")
  }
  
  Q1_bay_mean_q2 <- 
    t(apply(Q1_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_bay_mean_q2 <- rbind(Q1_bay_mean_q2, 
                          c(min(Q1_bay_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_mean_q2[, 2])))
  Q1_bay_mean_q2 <- Q1_bay_mean_q2[order(Q1_bay_mean_q2[, 1]), ]
  Q2_bay_mean_q2 <- 
    t(apply(Q2_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_bay_mean_q2 <- 
    t(apply(Q3_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  
  Q1_bay_lb_q2 <- 
    t(apply(Q1_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
  Q1_bay_lb_q2 <- rbind(Q1_bay_lb_q2, 
                          c(min(Q1_bay_lb_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_lb_q2[, 2])))
  Q1_bay_lb_q2 <- Q1_bay_lb_q2[order(Q1_bay_lb_q2[, 1]), ]
  Q2_bay_lb_q2 <- 
    t(apply(Q2_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
  Q3_bay_lb_q2 <- 
    t(apply(Q3_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
  
  Q1_bay_ub_q2 <- 
    t(apply(Q1_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))
  Q1_bay_ub_q2 <- rbind(Q1_bay_ub_q2, 
                          c(min(Q1_bay_ub_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_ub_q2[, 2])))
  Q1_bay_ub_q2 <- Q1_bay_ub_q2[order(Q1_bay_ub_q2[, 1]), ]
  Q2_bay_ub_q2 <- 
    t(apply(Q2_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))
  Q3_bay_ub_q2 <- 
    t(apply(Q3_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))
  
  plot(Q1_true_q2,
       xlim = c(0, 8),
       ylim = c(0, 8), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(1)))
  lines(Q1_emp_mean_q2 / 1000, col = "red", 
        lwd = 4)
  lines(Q1_lik_mean_q2/ 1000, col = "blue", 
        lwd = 4)
  lines(Q1_bay_mean_q2/ 1000, col = "forestgreen", 
        lwd = 4)
  lines(Q1_bay_lb_q2/ 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q1_bay_ub_q2/ 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q2_true_q2 / 1000, 
       xlim = c(0, 3),
       ylim = c(0, 2.5), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(2)))
  lines(Q2_emp_mean_q2 / 1000, col = "red", 
        lwd = 4)
  lines(Q2_lik_mean_q2 / 1000, col = "blue", 
        lwd = 4)
  lines(Q2_bay_mean_q2 / 1000, col = "forestgreen", 
        lwd = 4)
  lines(Q2_bay_lb_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q2_bay_ub_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q3_true_q2 / 1000, 
       xlim = c(0, 4.5),
       ylim = c(0, 4.5), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(3)))
  lines(Q3_emp_mean_q2 / 1000, col = "red", 
        lwd = 4)
  lines(Q3_lik_mean_q2 / 1000, col = "blue", 
        lwd = 4)
  lines(Q3_bay_mean_q2 / 1000, col = "forestgreen", 
        lwd = 4)
  lines(Q3_bay_lb_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q3_bay_ub_q2 / 1000, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
}

##### simulation setting 2 #####
sim_setting <- "xxx/"

folder_emp <- "xxx/"
filename_output_emp <- "Output_Empirical_Approach.rds"
results_emp2 <- readRDS(paste0(sim_setting, folder_emp, filename_output_emp))


folder_lik <- "xxx/"
filename_output_lik <- "Output_Likelihood_Frequentist_Approach.rds"
results_lik2 <- readRDS(paste0(sim_setting, folder_lik, filename_output_lik))

folder_bay <- "xxx/"
filename_output_bay <- "Output_Likelihood_Bayesian_Approach.rds"
results_bay2 <- readRDS(paste0(sim_setting, folder_bay, filename_output_bay))


filename_param <- "parameters_simulation_setting.rds"
params_emp2 <- readRDS(paste0(sim_setting, folder_emp, filename_param))
params_lik2 <- readRDS(paste0(sim_setting, folder_lik, filename_param))
params_bay2 <- readRDS(paste0(sim_setting, folder_bay, filename_param))


check_params_1 <- all(sapply(seq_len(length(params_emp2)), 
                             FUN = function(i) identical(params_emp2[[i]], params_lik2[[i]])))
check_params_2 <- all(sapply(seq_len(length(params_emp2)), 
                             FUN = function(i) identical(params_emp2[[i]], params_bay2[[i]])))

if (!check_params_1 | !check_params_2) warning("Simulation parameters do not match!")

params2 <- params_emp2

#### Plots full run ####
plot_full2 <- function() {
  nsim <- params2$nsim
  n <- params2$n
  rho <- params2$rho
  df <- params2$df
  q_n <- params2$q_n
  q <- params2$q
  w <- params2$w
  
  k <- (1 - q) * n
  
  ind_stuck <- rep(FALSE, nsim)
  for (i in seq_len(nsim)) {
    accepted <- results_bay2[[i]]$accepted
    acc_diff <- c(1, accepted) - c(accepted, 1)
    pos_0 <- which(acc_diff == -1)
    pos_1 <- which(acc_diff == 1)
    pos_diff <- pos_0 - pos_1
    if (any(pos_diff >= 300)) ind_stuck[i] <- TRUE
  }
  
  ind_keep <- which(!ind_stuck)[seq_len(1000)]
  
  Q1_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "or", m = 100)$Q
  
  Q1_hat_emp_q2 <- sapply(ind_keep, FUN = function(i) results_emp2[[i]]$Q1_hat_q2, simplify = "array")
  Q1_hat_lik_q2 <- sapply(ind_keep, FUN = function(i) results_lik2[[i]]$Q1_hat_q2, simplify = "array")
  Q1_hat_bay_q2 <- sapply(ind_keep, FUN = function(i) results_bay2[[i]]$Q1_hat_q2, simplify = "array")
  
  Q1_emp_mean_q2 <- t(apply(Q1_hat_emp_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_emp_mean_q2 <- rbind(Q1_emp_mean_q2, 
                          c(min(Q1_emp_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_emp_mean_q2[, 2])))
  Q1_emp_mean_q2 <- Q1_emp_mean_q2[order(Q1_emp_mean_q2[, 1]), ]
  Q1_lik_mean_q2 <- t(apply(Q1_hat_lik_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_lik_mean_q2 <- rbind(Q1_lik_mean_q2, 
                          c(min(Q1_lik_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_lik_mean_q2[, 2])))
  Q1_lik_mean_q2 <- Q1_lik_mean_q2[order(Q1_lik_mean_q2[, 1]), ]
  Q1_bay_mean_q2 <- t(apply(Q1_hat_bay_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_bay_mean_q2 <- rbind(Q1_bay_mean_q2, 
                          c(min(Q1_bay_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_mean_q2[, 2])))
  Q1_bay_mean_q2 <- Q1_bay_mean_q2[order(Q1_emp_mean_q2[, 1]), ]
  Q1_emp_lb_q2 <- t(apply(Q1_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q1_emp_lb_q2 <- rbind(Q1_emp_lb_q2, 
                        c(min(Q1_emp_lb_q2[, 1]), 10^10),
                        c(10^10, min(Q1_emp_lb_q2[, 2])))
  Q1_emp_lb_q2 <- Q1_emp_lb_q2[order(Q1_emp_lb_q2[, 1]), ]
  Q1_emp_ub_q2 <- t(apply(Q1_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q1_emp_ub_q2 <- rbind(Q1_emp_ub_q2, 
                        c(min(Q1_emp_ub_q2[, 1]), 10^10),
                        c(10^10, min(Q1_emp_ub_q2[, 2])))
  Q1_emp_ub_q2 <- Q1_emp_ub_q2[order(Q1_emp_ub_q2[, 1]), ]
  Q1_lik_lb_q2 <- t(apply(Q1_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q1_lik_lb_q2 <- rbind(Q1_lik_lb_q2, 
                        c(min(Q1_lik_lb_q2[, 1]), 10^10),
                        c(10^10, min(Q1_lik_lb_q2[, 2])))
  Q1_lik_lb_q2 <- Q1_lik_lb_q2[order(Q1_lik_lb_q2[, 1]), ]
  Q1_lik_ub_q2 <- t(apply(Q1_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q1_lik_ub_q2 <- rbind(Q1_lik_ub_q2, 
                        c(min(Q1_lik_ub_q2[, 1]), 10^10),
                        c(10^10, min(Q1_lik_ub_q2[, 2])))
  Q1_lik_ub_q2 <- Q1_lik_ub_q2[order(Q1_lik_ub_q2[, 1]), ]
  Q1_bay_lb_q2 <- t(apply(Q1_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q1_bay_lb_q2 <- rbind(Q1_bay_lb_q2, 
                        c(min(Q1_bay_lb_q2[, 1]), 10^10),
                        c(10^10, min(Q1_bay_lb_q2[, 2])))
  Q1_bay_lb_q2 <- Q1_bay_lb_q2[order(Q1_bay_lb_q2[, 1]), ]
  Q1_bay_ub_q2 <- t(apply(Q1_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q1_bay_ub_q2 <- rbind(Q1_bay_ub_q2, 
                        c(min(Q1_bay_ub_q2[, 1]), 10^10),
                        c(10^10, min(Q1_bay_ub_q2[, 2])))
  Q1_bay_ub_q2 <- Q1_bay_ub_q2[order(Q1_bay_ub_q2[, 1]), ]
  
  Q2_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "and", m = 100)$Q
  
  Q2_hat_emp_q2 <- sapply(ind_keep, FUN = function(i) results_emp2[[i]]$Q2_hat_q2, simplify = "array")
  Q2_hat_lik_q2 <- sapply(ind_keep, FUN = function(i) results_lik2[[i]]$Q2_hat_q2, simplify = "array")
  Q2_hat_bay_q2 <- sapply(ind_keep, FUN = function(i) results_bay2[[i]]$Q2_hat_q2, simplify = "array")
  
  Q2_emp_mean_q2 <- t(apply(Q2_hat_emp_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_lik_mean_q2 <- t(apply(Q2_hat_lik_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_bay_mean_q2 <- t(apply(Q2_hat_bay_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q2_emp_lb_q2 <- t(apply(Q2_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q2_emp_ub_q2 <- t(apply(Q2_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q2_lik_lb_q2 <- t(apply(Q2_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q2_lik_ub_q2 <- t(apply(Q2_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q2_bay_lb_q2 <- t(apply(Q2_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q2_bay_ub_q2 <- t(apply(Q2_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  
  Q3_true_q2 <- ellipse(sigma = matrix(c(1, rho, rho, 1), nrow = 2), 
                        df = df, prob = 1 - q_n[2],
                        npoints = 100, pos = TRUE)
  
  Q3_hat_emp_q2 <- sapply(ind_keep, FUN = function(i) results_emp2[[i]]$Q3_hat_q2, simplify = "array")
  Q3_hat_lik_q2 <- sapply(ind_keep, FUN = function(i) results_lik2[[i]]$Q3_hat_q2, simplify = "array")
  Q3_hat_bay_q2 <- sapply(ind_keep, FUN = function(i) results_bay2[[i]]$Q3_hat_q2, simplify = "array")
  
  Q3_emp_mean_q2 <- t(apply(Q3_hat_emp_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_lik_mean_q2 <- t(apply(Q3_hat_lik_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_bay_mean_q2 <- t(apply(Q3_hat_bay_q2, 1, 
                            FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_emp_lb_q2 <- t(apply(Q3_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q3_emp_ub_q2 <- t(apply(Q3_hat_emp_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q3_lik_lb_q2 <- t(apply(Q3_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q3_lik_ub_q2 <- t(apply(Q3_hat_lik_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  Q3_bay_lb_q2 <- t(apply(Q3_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.025, na.rm = T))))
  Q3_bay_ub_q2 <- t(apply(Q3_hat_bay_q2, 1, 
                          FUN = function(A) 
                            apply(A, 1, 
                                  FUN = function(x) quantile(x, probs = 0.975, na.rm = T))))
  
  plot(Q1_true_q2,
       xlim = c(0, 2.5),
       ylim = c(0, 2.5), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(1)))
  lines(Q1_emp_mean_q2 / 100, col = "red", 
        lwd = 4)
  lines(Q1_lik_mean_q2/ 100, col = "blue", 
        lwd = 4)
  lines(Q1_bay_mean_q2/ 100, col = "forestgreen", 
        lwd = 4)
  lines(Q1_emp_lb_q2/ 100, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q1_emp_ub_q2/ 100, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q1_lik_lb_q2/ 100, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q1_lik_ub_q2/ 100, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q1_bay_lb_q2/ 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q1_bay_ub_q2/ 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q2_true_q2 / 100, 
       xlim = c(0, 1),
       ylim = c(0, 1), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(2)))
  lines(Q2_emp_mean_q2 / 100, col = "red", 
        lwd = 4)
  lines(Q2_lik_mean_q2 / 100, col = "blue", 
        lwd = 4)
  lines(Q2_bay_mean_q2 / 100, col = "forestgreen", 
        lwd = 4)
  lines(Q2_emp_lb_q2 / 100, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q2_emp_ub_q2 / 100, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q2_lik_lb_q2 / 100, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q2_lik_ub_q2 / 100, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q2_bay_lb_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q2_bay_ub_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q3_true_q2 / 100, 
       xlim = c(0, 1.5),
       ylim = c(0, 1.5), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(3)))
  lines(Q3_emp_mean_q2 / 100, col = "red", 
        lwd = 4)
  lines(Q3_lik_mean_q2 / 100, col = "blue", 
        lwd = 4)
  lines(Q3_bay_mean_q2 / 100, col = "forestgreen", 
        lwd = 4)
  lines(Q3_emp_lb_q2 / 100, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q3_emp_ub_q2 / 100, col = "red", 
        lwd = 5, lty = "dashed")
  lines(Q3_lik_lb_q2 / 100, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q3_lik_ub_q2 / 100, col = "blue", 
        lwd = 5, lty = "dashed")
  lines(Q3_bay_lb_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q3_bay_ub_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
}

#### Plots single run ####
plot_single_run2 <- function(j = 1) {
  nsim <- params2$nsim
  n <- params2$n
  rho <- params2$rho
  df <- params2$df
  q_n <- params2$q_n
  q <- params2$q
  w <- params2$w
  
  k <- (1 - q) * n
  
  ind_stuck <- rep(FALSE, nsim)
  for (i in seq_len(nsim)) {
    accepted <- results_bay2[[i]]$accepted
    acc_diff <- c(1, accepted) - c(accepted, 1)
    pos_0 <- which(acc_diff == -1)
    pos_1 <- which(acc_diff == 1)
    pos_diff <- pos_0 - pos_1
    if (any(pos_diff >= 300)) ind_stuck[i] <- TRUE
  }
  
  ind_keep <- which(!ind_stuck)[seq_len(1000)]
  
  Q1_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "or", m = 100)$Q
  Q2_true_q2 <- qtbt(p = q_n[2], df = df, rho = rho, type = "and", m = 100)$Q
  Q3_true_q2 <- ellipse(sigma = matrix(c(1, rho, rho, 1), nrow = 2), 
                        df = df, prob = 1 - q_n[2],
                        npoints = 100, pos = TRUE)
  
  Q1_emp_mean_q2 <- results_emp2[[ind_keep[j]]]$Q1_hat_q2
  Q1_emp_mean_q2 <- rbind(Q1_emp_mean_q2, 
                          c(min(Q1_emp_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_emp_mean_q2[, 2])))
  Q1_emp_mean_q2 <- Q1_emp_mean_q2[order(Q1_emp_mean_q2[, 1]), ]
  Q1_lik_mean_q2 <- results_lik2[[ind_keep[j]]]$Q1_hat_q2
  Q1_lik_mean_q2 <- rbind(Q1_lik_mean_q2, 
                          c(min(Q1_lik_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_lik_mean_q2[, 2])))
  Q1_lik_mean_q2 <- Q1_lik_mean_q2[order(Q1_lik_mean_q2[, 1]), ]
  Q2_emp_mean_q2 <- results_emp2[[ind_keep[j]]]$Q2_hat_q2
  Q2_lik_mean_q2 <- results_lik2[[ind_keep[j]]]$Q2_hat_q2
  Q3_emp_mean_q2 <- results_emp2[[ind_keep[j]]]$Q3_hat_q2
  Q3_lik_mean_q2 <- results_lik2[[ind_keep[j]]]$Q3_hat_q2
  
  theta_hat1_bayesian <- results_bay2[[ind_keep[j]]]$theta_X1_hat
  theta_hat2_bayesian <- results_bay2[[ind_keep[j]]]$theta_X2_hat
  k_bayesian <- results_bay2[[ind_keep[j]]]$k_bayesian
  eta_bayesian <- results_bay2[[ind_keep[j]]]$eta_hat
  
  m <- length(results_bay2[[ind_keep[j]]]$w)
  n_iter <- length(k_bayesian)
  Q1_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q2_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  Q3_bay_q2_array <- array(NA, dim = c(m, 2, n_iter))
  
  for (i in seq_len(n_iter)) {
    K <- k_bayesian[i] + 1
    eta <- eta_bayesian[i , seq_len(K)]
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
    
    theta_hat_b <- cbind(theta_hat1_bayesian[i, ],
                         theta_hat2_bayesian[i, ])
    
    Q1_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                    theta = theta_hat_b,
                                    A = A_bay,
                                    t = n / k, type = "or")
    Q2_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                    theta = theta_hat_b,
                                    A = A_bay,
                                    t = n / k, type = "and")
    Q3_bay_q2_array[, , i] <- Q_hat(w = w, p = q_n[2],
                                    theta = theta_hat_b,
                                    h = h_bay,
                                    t = n / k, type = "density")
  }
  
  Q1_bay_mean_q2 <- 
    t(apply(Q1_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q1_bay_mean_q2 <- rbind(Q1_bay_mean_q2, 
                          c(min(Q1_bay_mean_q2[, 1]), 10^10),
                          c(10^10, min(Q1_bay_mean_q2[, 2])))
  Q1_bay_mean_q2 <- Q1_bay_mean_q2[order(Q1_bay_mean_q2[, 1]), ]
  Q2_bay_mean_q2 <- 
    t(apply(Q2_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  Q3_bay_mean_q2 <- 
    t(apply(Q3_bay_q2_array, 1, FUN = function(x) rowMeans(x, na.rm = T)))
  
  Q1_bay_lb_q2 <- 
    t(apply(Q1_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
  Q1_bay_lb_q2 <- rbind(Q1_bay_lb_q2, 
                        c(min(Q1_bay_lb_q2[, 1]), 10^10),
                        c(10^10, min(Q1_bay_lb_q2[, 2])))
  Q1_bay_lb_q2 <- Q1_bay_lb_q2[order(Q1_bay_lb_q2[, 1]), ]
  Q2_bay_lb_q2 <- 
    t(apply(Q2_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
  Q3_bay_lb_q2 <- 
    t(apply(Q3_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.025))))
  
  Q1_bay_ub_q2 <- 
    t(apply(Q1_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))
  Q1_bay_ub_q2 <- rbind(Q1_bay_ub_q2, 
                        c(min(Q1_bay_ub_q2[, 1]), 10^10),
                        c(10^10, min(Q1_bay_ub_q2[, 2])))
  Q1_bay_ub_q2 <- Q1_bay_ub_q2[order(Q1_bay_ub_q2[, 1]), ]
  Q2_bay_ub_q2 <- 
    t(apply(Q2_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))
  Q3_bay_ub_q2 <- 
    t(apply(Q3_bay_q2_array, 1, 
            FUN = function(x) apply(x, 1, FUN = function(z) quantile(z, probs = 0.975))))
  
  plot(Q1_true_q2,
       xlim = c(0, 3),
       ylim = c(0, 2), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(1)))
  lines(Q1_emp_mean_q2 / 100, col = "red", 
        lwd = 4)
  lines(Q1_lik_mean_q2/ 100, col = "blue", 
        lwd = 4)
  lines(Q1_bay_mean_q2/ 100, col = "forestgreen", 
        lwd = 4)
  lines(Q1_bay_lb_q2/ 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q1_bay_ub_q2/ 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q2_true_q2 / 100, 
       xlim = c(0, 1.6),
       ylim = c(0, 0.8), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(2)))
  lines(Q2_emp_mean_q2 / 100, col = "red", 
        lwd = 4)
  lines(Q2_lik_mean_q2 / 100, col = "blue", 
        lwd = 4)
  lines(Q2_bay_mean_q2 / 100, col = "forestgreen", 
        lwd = 4)
  lines(Q2_bay_lb_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q2_bay_ub_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
  
  plot(Q3_true_q2 / 100, 
       xlim = c(0, 2.5),
       ylim = c(0, 1.2), 
       type = "l", lwd = 4, ylab = "", xlab = "",
       cex.lab=2, cex.axis=2, cex.main = 2,
       las = 1, cex.sub = 2,
       main = expression(hat(Q)[n]^(3)))
  lines(Q3_emp_mean_q2 / 100, col = "red", 
        lwd = 4)
  lines(Q3_lik_mean_q2 / 100, col = "blue", 
        lwd = 4)
  lines(Q3_bay_mean_q2 / 100, col = "forestgreen", 
        lwd = 4)
  lines(Q3_bay_lb_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  lines(Q3_bay_ub_q2 / 100, col = "forestgreen", 
        lwd = 5, lty = "dashed")
  legend("topright", 
         legend = c("Empirical", 
                    "Lik.-Freq.",
                    "Lik.-Bayes."),
         lty = 1, col = c("red", "blue", "forestgreen"),
         bty = "n", lwd = 4, cex=2, seg.len=1)
}

#### Generate plots ####
pdf("simulation_study.pdf", width=6, height=6)
par(mai=c(.7,1.2,.7,.2), mgp=c(3,0.8,0))
plot_single_run1()
plot_single_run2()

plot_full1()
plot_full2()
dev.off()
