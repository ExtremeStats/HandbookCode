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
#   (2025). In: Handbook on Statistics of Extremes. Chapman & Hall/CRC. Boca   #
#   Raton, FL.                                                                 #
#                                                                              #
#  ==========================================================================  #
  
## CLEAN WORKSPACE & LOAD PACKAGES
rm(list = ls())
packages <- c("coda", "dplyr", "ggplot2", "gridExtra", "MASS", "splines",
              "ismev2")
sapply(packages, require, character.only = TRUE)

data(california)
attach(california)

# LOAD & PREPARE DATA
y <- Acres; n <- length(Acres)
t <- as.numeric(california$Date) / 10000
u <- quantile(y, 0.99)
t <- t[y > u]
y <- y[y > u] 
x <- cbind(rep(1, length(y)), bs(t))
q <- ncol(x)

# LOG LIKELIHOOD & MLE
l <- function(beta)
  -sum(exp(-x %*% beta) * log(y / u) + x %*% beta) 
mle <- optim(rep(1, q), l, hessian = TRUE, control = list(fnscale = -1))
betahat <- mle$par
I <- solve(-mle$hessian)

## SETUP MH 
burn <- 1000
T <- 10000
betas <- matrix(1, nrow = T, ncol = 4)
c <- 10^2
set.seed(1)

## LOG POSTERIOR 
p <- function(beta) {
  -sum(exp(-x %*% beta) * log(y / u) + x %*% beta) + 
   sum(dnorm(beta, sd = c, log = TRUE))
}

## METROPOLIS-HASTINGS FOR PT REGRESSION
for (i in 1:(T - 1)) {
  # generate from proposal
  betastar <- mvrnorm(1, betas[i, ], I)
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

# examine MCMC outputs
par(mfrow = c(2, 2))
plot(betas[, 1], type = "l")
abline(h = betahat[1])
plot(betas[, 2], type = "l")
abline(h = betahat[2])
plot(betas[, 3], type = "l")
abline(h = betahat[3])
plot(betas[, 4], type = "l")
abline(h = betahat[4])
effectiveSize(betas[, 4])

xitraj <- exp(betas[(burn + 1): T, ] %*% t(x))
lband <- apply(xitraj, 2, quantile, probs = 0.025)
xibayes <- apply(xitraj, 2, median)
ximle <- exp(x %*% betahat)
uband <- apply(xitraj, 2, quantile, probs = 0.975)

## VISUALIZATIONS
# plot fit and data
df <- data.frame(date = as.Date(t * 10000, origin = "1970-01-01"), 
                 y = xibayes)

thresholded_data <- california %>% 
  filter(Acres > u)

# Load the necessary libraries
library(zoo)      # For rolling windows
data("california")
# function to compute the Hill estimator for a given window
#hill_estimator <- function(data) {
#  u <- quantile(data, 0.95)
#  iexc <- which(data > u)
#  exc <- data[iexc]
#  k <- length(exc)
#  1 / k * sum(log(exc / u))
#}
#window_size <- 500
#rolling_hill <- rollapply(california$Acres, window_size, hill_estimator, by.column = FALSE, fill = NA)
# Plot the results

a <- ggplot(df, aes(x = date, y = y)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_ribbon(aes(ymin = lband, ymax = uband), fill = 'steelblue', alpha = 0.2) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10)) + 
  #geom_line(data = data.frame(x = california$Date, y = rolling_hill), aes(x = x, y = y))
  labs(x = "Time (years)", y = expression(xi[t]), 
       title = "a) Conditional Extreme Value Index") 

b <- ggplot(thresholded_data, aes(x = Date, y = log(Acres), color = log(Acres))) +
  geom_point(alpha = 0.5) + 
  scale_color_gradient(low = "steelblue1", high = "red") + 
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(size = 10)) + 
  labs(x = "Time (years)", y = "Acres (log)", 
       title = "b) Acres Above Threshold (log scale)")

## QQ-plot
k <- length(y)
etraj <- matrix(nrow = T - burn, ncol = k)
for(i in 1:(T - burn))
  for(j in 1:k)
    etraj[i, j] <- qnorm(1 - 1/(y[j]/u)^(1 / xitraj[i, j]))
for(i in 1:(T - burn))
  etraj[i, ] <- sort(etraj[i, ])

elband <- apply(etraj, 2, quantile, probs = 0.025)
emedian <- apply(etraj, 2, median)
euband <- apply(etraj, 2, quantile, probs = 0.975)

r <- ggplot() +
  geom_qq_line(data = as.data.frame(emedian), aes(sample = emedian), 
               colour = "gray", linewidth = 2) + 
  geom_ribbon(aes(x = qnorm(1:k / (k + 1)), ymin = elband, ymax = euband), 
              fill = 'steelblue', alpha = 0.2) +
  geom_point(data = data.frame(x = qnorm(1:k / (k + 1)), y = sort(emedian)), 
             aes(x = x, y = y), color = "steelblue", alpha = 0.3, size = 0.3) + 
  theme_minimal() +
  xlim(-3, 3) + 
  ylim(-3, 3) + 
  labs(title = "c) QQ-Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") + 
  theme(legend.position = "none", plot.title = element_text(size = 10)) 

## Intensity function 
s <- ggplot(data = data.frame(x = df$date), aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), fill = "steelblue", 
                 color = "gray", alpha = .2) +
  #scale_fill_gradient(low = "lightblue", high = "red") +  
  geom_rug(color = "gray", alpha = 0.3) + 
  labs(title = "d) Intensity Function") + 
  xlab('Time (years)') + 
  ylab('Density') + 
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(size = 10)) 

gg <- grid.arrange(a, b, r, s, ncol = 2)

ggsave("/Users/mdecarvalho/Dropbox/Apps/Overleaf/Handbook\ on\ Statistics\ of\ Extremes\ -\ Ch6/Ch6/LaTeX/figures/PTR_california.pdf", 
       gg, width = 7, height = 7)

# qqboxplot of residuals evaluated at posterior median 
packages <- c("ggplot2", "gridExtra", "qqboxplot")
sapply(packages, require, character.only = TRUE)

ggplot(data = data.frame(y = emedian), aes(y = y)) + 
  geom_qqboxplot(notch = TRUE, varwidth = TRUE, reference_dist = "norm") + 
  ylab('Normal Scale') + 
  ggtitle('California Wildfires') + 
  theme_minimal() + 
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  ylim(-3, 3) + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 23)) 
ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/qqboxtir.pdf", width = 6)