#  ==========================================================================  #
#  This program is free software; you can redistribute it and/or modify        #
#  it under the terms of the GNU General Public License as published by        #
#  the Free Software Foundation; either version 2 of the License, or           #
#  (at your option) any later version.                                         #
#                                                                              #
#  This program is distributed in the hope that it will be useful,             #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of              #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
#  GNU General Public License for more details.                                #
#                                                                              #
#  You should have received a copy of the GNU General Public License           #
#  along with this program; if not, a copy is available at                     #
#  http://www.r-project.org/Licenses/                                          #
#  --------------------------------------------------------------------------  #
#                                                                              #  
#  This code was primarily designed for instructional purposes. Effective      #
#  sample sizes can be improved as noted in de Carvalho et al (2022). For      #
#  details see the references below.                                           #
#                                                                              #
#  References                                                                  #
#  - de Carvalho, M., Palacios Ramirez, V., Henriques-Rodrigues, L. & Lee, J.  #  
#   (2026). In: Handbook on Statistics of Extremes. Chapman & Hall/CRC. Boca   #
#   Raton, FL.                                                                 #
#                                                                              #
#  - de Carvalho, M., Pereira, S., Pereira, S. & de Zea Bermudez, P. (2022).   #
#   "An Extreme Value Bayesian Lasso for the Conditional Left and Right Tails" #
#    Journal of Agricultural, Biological and Environmental Statistics, 27,     # 
#    222–39.                                                                   # 
#                                                                              #
#  ==========================================================================  #

## CLEAN WORKSPACE & LOAD PACKAGES
rm(list = ls())
packages <- c("coda", "DATAstudio", "dplyr", "ismev2", "MASS", "lubridate", "ggplot2")
sapply(packages, require, character.only = TRUE)

## LOAD & PREPARE DATA
data(madeira)
attach(madeira)
q <- dim(madeira)[2] - 1
mad <- filter(madeira, prec > 0)
n <- dim(mad)[1]
y <- mad$prec
x <- as.matrix(cbind(ones = rep(1, n), mad[, 3:(q + 1)]))
x[, -1] <- scale(x[, -1])

## VISUALIZE RAW DATA
mad$yearmonth <- ymd(paste0(mad$yearmonth, "01"))
pdata <- ggplot(data = data.frame(x = mad$yearmonth, y = y)) +
  geom_point(data = mad, aes(x = yearmonth, y = y, color  = y), 
             alpha = 0.3, shape = 16, size = 3, show.legend = FALSE) +
  scale_color_gradient(low = "steelblue1", high = "red") + 
  xlab('Time (years)') +
  ylab('Total Monthly Precipitation') + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16)) 

ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/madeira.pdf", width = 10, height = 5)  

## LOG LIKELIHOOD & MLE
l <- function(beta) {
  betas <- beta[1:q]; betax <- beta[(q+1):(2*q)]; betak <- beta[(2*q+1):(3*q)]
  h <- H <- numeric()
  sig <- exp(x %*% betas)
  sh <- x %*% betax
  # if outside the support apply penalty  
  if (any(1 + sh * y / sig  <= 0)) {
    ans <- -10^10
  } else {
  # standard likelihood  
    for (i in 1:n) {
      h[i] <- ismev2::dgpd(y[i], sig = sig[i], sh = sh[i], log = TRUE)
      H[i] <- log(ismev2::pgpd(y[i], sig = sig[i], sh = sh[i]))
    }
    ans <- sum(h + x %*% betak + (exp(x %*% betak) - 1) * H)
  }
  return(ans)
}

mle <- optim(rep(0, 3 * q), l, method = "L-BFGS-B", hessian = TRUE, 
             lower = -30, upper = 30, control = list(fnscale = -1))
mle$convergence
betahat <- mle$par
Iinv <- solve(-mle$hessian)

## MCMC SETTINGS & PRIOR
burn <- 1000
T <- 20000
betas <- matrix(0, nrow = T, ncol = 3 * q)
betas[1, ] <- betahat
c <- 10^2
set.seed(1)

## LOG POSTERIOR
p <- function(beta) 
  l(beta) - sum(dnorm(beta, sd = c, log = TRUE))

## INDEPENDENCE SAMPLER FOR EGP REGRESSION
for (i in 1:(T - 1)) {
  # generate from proposal
  betastar <- mvrnorm(1, betahat, Iinv)
  # evaluate log posterior of proposal vs current
  lpstar <- p(betastar) 
  lp <- p(betas[i, ])
  alpha <- exp(lpstar - lp)
  # accept or reject
  if (alpha > runif(1)) {
    betas[i + 1, ] <- betastar
  } else {
    betas[i + 1, ] <- betas[i, ]
  } 
}

## EXAMINE MCMC OUTPUTS
# mixing can be improved as noted above
par(mfrow = c(2, 2))
plot(betas[, 1], type = "l")
abline(h = betahat[1])
plot(betas[, 2], type = "l")
abline(h = betahat[2])
plot(betas[, 3], type = "l")
abline(h = betahat[3])
plot(betas[, 4], type = "l")
abline(h = betahat[4])
effectiveSize(betas[burn: T, ])

## VISUALIZATIONS
packages <- c("ggplot2", "gridExtra")
sapply(packages, require, character.only = TRUE)

# fitted regression coefficients
l <- apply(betas, 2, quantile, 0.025)
m <- apply(betas, 2, median)
u <- apply(betas, 2, quantile, 0.975)

# regression coefficients for kappa
a <- ggplot(data.frame(l = l[(2 * q + 2):(3 * q)], 
                       m = m[(2 * q + 2):(3 * q)], 
                       u = u[(2 * q + 2):(3 * q)], group = factor(2:q)), 
            aes(group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1.4) +
  geom_point(aes(x = group, y = m), size = 4.5) + 
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, size = 1.2) + 
  labs(x = "Regression Coefficients", y = "") + 
  coord_cartesian(ylim=c(-0.8, 0.6))+
  ggtitle(expression(paste("", kappa, "(x)"))) + 
  scale_x_discrete(labels = c("AMO","ENSO", "NP", "PDO", "SOI", "NAO")) + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 31)) 

# regression coefficients for sigma
b <- ggplot(data.frame(l = l[2:q], m = m[2:q], u = u[2:q], group = factor(2:q)), 
            aes(group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.4) +
  geom_point(aes(x = group, y = m), size = 4.5) + 
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth = 1.2) + 
  labs(x = "Regression Coefficients", y = "") + 
  coord_cartesian(ylim=c(-0.8, 0.6))+
  ggtitle(expression(paste("", sigma, "(x)"))) + 
  scale_x_discrete(labels = c("AMO","ENSO", "NP", "PDO", "SOI", "NAO")) + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 31)) 

# regression coefficients for xi
c <- ggplot(data.frame(l = l[(q + 2):(2 * q)], 
                       m = m[(q + 2):(2 * q)], 
                       u = u[(q + 2):(2 * q)], group = factor(2:q)), 
            aes(group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1.4) +
  geom_point(aes(x = group, y = m), size = 4.5) + 
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, size = 1.2) + 
  labs(x = "Regression Coefficients", y = "") + 
  coord_cartesian(ylim=c(-0.8, 0.6))+
  ggtitle(expression(paste("", xi, "(x)"))) + 
  scale_x_discrete(labels = c("AMO","ENSO", "NP", "PDO", "SOI", "NAO")) + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 31)) 

# grid.arrange(a, b, c, ncol = 3)

ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/a.pdf", a)
ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/b.pdf", b)
ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/c.pdf", c)

## QQ-plot
betasig <- betas[1:T, 1:q] 
betash <- betas[1:T, (q+1):(2*q)] 
betak <- betas[1:T, (2*q+1):(3*q)]
etraj <- matrix(nrow = T - burn, ncol = n)
for(i in 1:(T - burn))
  for(j in 1:n)
    etraj[i, j] <- qnorm(pegpd(y[j], exp(x[j, ] %*% betak[i, ]), 
                               exp(x[j, ] %*% betasig[i, ]), x[j, ] %*% betash[i, ]))
for(i in 1:(T - burn))
  etraj[i, ] <- sort(etraj[i, ])

elband <- apply(etraj, 2, quantile, probs = 0.025)
emedian <- apply(etraj, 2, median)
euband <- apply(etraj, 2, quantile, probs = 0.975)

r <- ggplot() +
  geom_qq_line(data = as.data.frame(emedian), aes(sample = emedian), 
               colour = "gray", linewidth = 2) + 
  geom_ribbon(aes(x = qnorm(1:n / (n + 1)), ymin = elband, ymax = euband), 
              fill = 'steelblue', alpha = 0.2) +
  geom_point(data = data.frame(x = qnorm(1:n / (n + 1)), y = sort(emedian)), 
             aes(x = x, y = y), color = "steelblue", alpha = 0.3, size = 0.3) + 
  theme_minimal() +
  xlim(-3.4, 3.4) + 
  ylim(-3.4, 3.4) + 
  labs(title = "c) QQ Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") + 
  theme(legend.position = "none") 

# qqboxplot of residuals evaluated at posterior median 
packages <- c("ggplot2", "gridExtra", "qqboxplot")
sapply(packages, require, character.only = TRUE)

ggplot(data = data.frame(y = emedian), aes(y = y)) + 
  geom_qqboxplot(notch = TRUE, varwidth = TRUE, reference_dist = "norm") + 
  ylab('Normal Scale') + 
  ggtitle('Madeira Rainfall') + 
  theme_minimal() + 
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  ylim(-3, 3) + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 23)) 
ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/qqboxegp.pdf", width = 6)

## Density 
# Histogram of y
p <- ggplot() +
  geom_histogram(aes(x = y, y = ..density..), bins = 30, fill = "steelblue", 
                 color = "gray", alpha = 0.2) +
  xlab('Precipitation') + 
  ggtitle("Histogram") + 
  ylab('Density') + 
  ylim(0, 0.5) + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 23)) 

# Density conditional on mean covariate 
xx <- colMeans(x)
XX <- c(xx, xx, xx)
dens <- function(y) {
  denstraj <- numeric()
  for (i in 1:dim(betas)[1])
    denstraj[i] <- (ismev2::dgpd(y, exp(xx %*% betasig[i,]), xx %*% betash[i,])) * exp(xx %*% betak[i,]) * 
      (ismev2::pgpd(y, exp(xx %*% betasig[i,]), xx %*% betash[i,]))^(exp(xx %*% betak[i,]) - 1)
  return(denstraj)
}

ygrid <- seq(0.05, 16, by = 0.1)
lgrid <- length(ygrid)
dtraj <- matrix(0, nrow = lgrid, ncol = length(dens(1)))
for(i in 1:lgrid)
  dtraj[i, ] <- dens(ygrid[i])
L <- apply(dtraj, 1, quantile, 0.025)
M <- apply(dtraj, 1, median)
U <- apply(dtraj, 1, quantile, 0.975)
data <- data.frame(
  ygrid = ygrid,
  L = L,
  M = M,
  U = U
)

q <- ggplot() +
  geom_ribbon(data = data, aes(x = ygrid, ymin = L, ymax = U), fill = 'steelblue', alpha = 0.2) +
  geom_line(data = data, aes(x = ygrid, y = M), color = "steelblue", size = 1.2) +
  labs(
    x = "y",
    y = "Density"
  ) +
  ylim(0, 0.5) + 
  geom_line(data = data, aes(x = ygrid, y = M), color = "steelblue", size = 1.2) +
  labs(
    x = "y",
    y = "Density",
    title = "Given Mean Covariate"
  ) +
  xlab('Precipitation') + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 23)) 

ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/histmadeira.pdf", p)
ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/dens1madeira.pdf", q)

# Density conditional on maximal covariate 
xx <- x[70,]
XX <- c(xx, xx, xx)
dens <- function(y) {
  denstraj <- numeric()
  for (i in 1:dim(betas)[1])
    denstraj[i] <- (ismev2::dgpd(y, exp(xx %*% betasig[i,]), xx %*% betash[i,])) * exp(xx %*% betak[i,]) * 
      (ismev2::pgpd(y, exp(xx %*% betasig[i,]), xx %*% betash[i,]))^(exp(xx %*% betak[i,]) - 1)
  return(denstraj)
}

ygrid <- seq(0.05, 16, by = 0.1)
lgrid <- length(ygrid)
dtraj <- matrix(0, nrow = lgrid, ncol = length(dens(1)))
for(i in 1:lgrid)
  dtraj[i, ] <- dens(ygrid[i])
L <- apply(dtraj, 1, quantile, 0.025)
M <- apply(dtraj, 1, median)
U <- apply(dtraj, 1, quantile, 0.975)
data <- data.frame(
  ygrid = ygrid,
  L = L,
  M = M,
  U = U
)

# Credible interval (ribbon)
r <- ggplot() +
  geom_ribbon(data = data, aes(x = ygrid, ymin = L, ymax = U), fill = "steelblue", alpha = 0.2) +
  geom_line(data = data, aes(x = ygrid, y = M), color = "steelblue", size = 1.2) +
  labs(
    x = "y",
    y = "Density"
  ) +
  xlab('Precipitation') + 
  ylim(0, 0.5) + 
  geom_line(data = data, aes(x = ygrid, y = M), color = "steelblue", size = 1.2) +
  # Labels and theme
  labs(
    x = "Precipitation",
    y = "Density",
    title = "Given Maximal Covariate"
  ) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 23)) 

ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/dens2madeira.pdf", r)

## COMPUTE WAIC
require(loo)
# (L is log-likelihood matrix with columns for observations and rows for iterations)
L <- matrix(NA, nrow = T - burn, ncol = n)
for(t in 1:(T - burn))
  for(i in 1:n)
    l[t, i] <- log(degpd(y[i], exp(x[i, ] %*% betak[t, ]), 
                         exp(x[i, ] %*% betas[t, ]), x[i, ] %*% betax[t, ]))
model.eval <- waic(L)
model.eval$estimates
WAIC <- -2 * (model.eval$estimates[1, 1])