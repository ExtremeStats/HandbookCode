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
#  This code is mainly written with instructional purposes in mind. For        #
#  details see the references below.                                           #
#                                                                              #
#  References                                                                  #
#  - de Carvalho, M., Palacios Ramirez, V., Henriques-Rodrigues, L. & Lee, J.  #  
#   (2025). In: Handbook on Statistics of Extremes. Chapman & Hall/CRC. Boca   #
#   Raton, FL.                                                                 #
#                                                                              #
#  ==========================================================================  #
  
## CLEAN WORKSPACE & LOAD PACKAGES
#rm(list = ls())
packages <- c("DATAstudio", "dplyr", "evd", "ggplot2", "gridExtra", "ismev2", 
              "MASS", "lubridate", "splines")
sapply(packages, require, character.only = TRUE)

## LOAD & PREPARE DATA
hblockmax <- hongkong %>%
  mutate(month = floor_date(date, "year")) %>%  
  group_by(month) %>%
  summarize(max_value = max(value))
attach(hblockmax)

## LOG LIKELIHOOD & MLE
y <- max_value; n <- length(max_value)
t <- as.numeric(month) / 10000
x <- cbind(rep(1, length(y)), bs(t, degree = 3))
q <- ncol(x)
l <- function(beta) {
    betamu <- beta[1:q]; betasig <- beta[(q+1):(2*q)]; betash <- beta[(2*q+1):(3*q)]
    mu <- x %*% betamu
    sig <- exp(x %*% betasig)
    sh <- x %*% betash
    ll <- as.numeric()
    # if outside the support apply penalty  
    if (any(1 + sh * (y - mu) / sig  <= 0) || any(sig <= 0)) {
      ans <- -10^20
    } else {
    # standard likelihood  
      for (i in 1:n) {
        ll[i] <- dgev(y[i], mu[i], sig[i], sh[i], log = TRUE)
      }  
      ans <- sum(ll)  
    }
    return(ans)
}
init <- gev.fit(y)$mle
start <- c(init[1], rep(0, q - 1), log(init[2]), rep(0, q - 1), init[3], rep(0, q - 1))
mle <- optim(start, l, hessian = TRUE, method = "BFGS", 
             control = list(fnscale = -1))
betahat <- mle$par
Iinv <- solve(-mle$hessian, tol = 10^(-35))
print(mle$par)

## PROJECT INVERSE FISHER INFORMATION INTO THE SPACE OF POSITIVE DEFINITE MATRICES
packages <- c("matrixcalc", "Matrix")
sapply(packages, require, character.only = TRUE)
Iinv <- nearPD(Iinv)$mat
test <- is.positive.definite(as.matrix(Iinv))
print(test)

## MCMC SETTINGS & PRIOR
burn <- 1000
T <- 10000
betas <- matrix(0, nrow = T, ncol = 3 * q)
betas[1, ] <- betahat
c <- 10^2
set.seed(1)

## LOG POSTERIOR
p <- function(beta) {
  betamu <- beta[1:q]; betasig <- beta[(q+1):(2*q)]; betash <- beta[(2*q+1):(3*q)]
  return(l(beta) + sum(dnorm(betamu, sd = c, log = TRUE)) + 
         sum(dnorm(betasig, sd = c, log = TRUE)) +
         sum(dnorm(betash, sd = c, log = TRUE)))
}

## INDEPENDENCE SAMPLER FOR GEV REGRESSION
for (i in 1:(T - 1)) {
  # generate from proposal
  parsstar <- mvrnorm(1, betahat, Iinv)
  # evaluate log posterior of proposal vs current
  lpstar <- p(parsstar) 
  lp <- p(betas[i, ])
  alpha <- exp(lpstar - lp)
  # accept or reject
  if (alpha > runif(1)) {
    betas[i + 1, ] <- parsstar
  } else {
    betas[i + 1, ] <- betas[i, ]
  } 
}

## EXAMINE MCMC OUTPUTS
require(coda)
mcmc <- mcmc(betas)
effectiveSize(mcmc)
## regression coefficients for mu
par(mfrow = c(2, 2))
plot(betas[, 1], type = "l")
abline(h = mle$par[1])
plot(betas[, 2], type = "l")
abline(h = mle$par[2])
plot(betas[, 3], type = "l")
abline(h = mle$par[3])
plot(betas[, 4], type = "l")
abline(h = mle$par[4])
## regression coefficients for sigma
plot(betas[, 5], type = "l")
abline(h = mle$par[5])
plot(betas[, 6], type = "l")
abline(h = mle$par[6])
plot(betas[,7], type = "l")
abline(h = mle$par[7])
plot(betas[, 8], type = "l")
abline(h = mle$par[8])
## regression coefficients for xi
plot(betas[, 9], type = "l")
abline(h = mle$par[9])
plot(betas[, 10], type = "l")
abline(h = mle$par[10])
plot(betas[, 11], type = "l")
abline(h = mle$par[11])
plot(betas[, 12], type = "l")
abline(h = mle$par[12])

# fitted regression coefficients and CI
betamu <-  betas[(burn + 1):T, 1:q]
btrajmu <- betamu %*% t(x)
lbandmu <- apply(btrajmu, 2, quantile, probs = 0.025)
bmedianmu <- apply(btrajmu, 2, median)
ubandmu <- apply(btrajmu, 2, quantile, probs = 0.975)
# plot fit and data
dfmu <- data.frame(date = as.Date(t * 10000, origin = "1970-01-01"), 
                   y = bmedianmu)

pmu <- ggplot(dfmu, aes(x = date, y = y)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_ribbon(aes(ymin = lbandmu, ymax = ubandmu), fill = 'steelblue', 
              alpha = 0.2) +
  theme_minimal() +
  labs(x = "Time (years)", y = expression(mu[t]), title = "")+
  scale_y_continuous(name = expression(mu[t]), breaks = seq(0, 50, by = 2), 
                     limits = c(30, 38))

betasig <- betas[(burn+1):T, (q+1):(2*q)]
btrajsig <- exp(betasig %*% t(x))
lbandsig <- apply(btrajsig, 2, quantile, probs = 0.025)
bmediansig <- apply(btrajsig, 2, median)
ubandsig <- apply(btrajsig, 2, quantile, probs = 0.975)
# plot fit and data
dfsig <- data.frame(date = as.Date(t * 10000, origin = "1970-01-01"), 
                   y = bmediansig)

psig <- ggplot(dfsig, aes(x = date, y = y)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_ribbon(aes(ymin = lbandsig, ymax = ubandsig), fill = 'steelblue', 
              alpha = 0.2) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank()) + 
  labs(x = "", y = expression(sigma[t]), title = "") + 
  ggtitle("b) Scale and Shape")

betash <- betas[(burn + 1):T, (2*q+1):(3*q)]
btrajsh <- betash %*% t(x)
lbandsh <- apply(btrajsh, 2, quantile, probs = 0.025)
bmediansh <- apply(btrajsh, 2, median)
ubandsh <- apply(btrajsh, 2, quantile, probs = 0.975)
# plot fit and data
dfsh <- data.frame(date = as.Date(t * 10000, origin = "1970-01-01"), 
                    y = bmediansh)

psh <- ggplot(dfsh, aes(x = date, y = y)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_ribbon(aes(ymin = lbandsh, ymax = ubandsh), fill = 'steelblue', 
              alpha = 0.2) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12)) + 
  labs(x = "Time (years)", y = expression(xi[t]), title = "") + 
  scale_y_continuous(name = expression(xi[t]), breaks = seq(-0.7, 1, by = 0.4), 
                     limits = c(-0.7, 1)) 

## plot data, mu, and 99% time-varying quantile 
traj.tvq <- matrix(nrow = T - burn, ncol = n)
for(i in 1:(T - burn))
  for(j in 1:n)
    traj.tvq[i, j] <-  qgev(0.95, x[j, ] %*% betamu[i, ], 
                           exp(x[j, ] %*% betasig[i, ]), x[j, ] %*% betash[i, ])
lband.tvq <- apply(traj.tvq, 2, quantile, probs = 0.025)
median.tvq <- apply(traj.tvq, 2, quantile, probs = 0.5)
uband.tvq <- apply(traj.tvq, 2, quantile, probs = 0.975)
  
pdata <- ggplot() +
  geom_point(data = hblockmax, aes(x = month, y = max_value, color  = exp(y)), 
             alpha = 0.3, shape = 16, size = 3, show.legend = FALSE) +
  scale_color_gradient(low = "steelblue1", high = "red") + 
  xlab('Time (years)') +
  ylab('Maximum temperature (ºC)') + 
  geom_ribbon(data = dfmu, aes(x = date, y = y, ymin = lbandmu, ymax = ubandmu), 
              fill = 'steelblue', alpha = 0.2) +
  geom_line(data = dfmu, aes(x = date, y = y),colour = "steelblue", linewidth = 1) +
  geom_line(data = data.frame(x = dfmu$date, y = median.tvq), aes(x = x, y = y),colour = "red", linetype = "dashed") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12)) + 
  ggtitle("a) Location")

p <- grid.arrange(pdata,arrangeGrob(psig, psh), ncol = 2) 

ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/fitgev.pdf",p)  

## QQ-plot
etraj <- matrix(nrow = T - burn, ncol = n)
for(i in 1:(T - burn))
  for(j in 1:n)
    etraj[i, j] <- qnorm(pgev(y[j], x[j, ] %*% betamu[i, ], 
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
  xlim(-3.9, 3.9) + 
  ylim(-3.9, 3.9) + 
  labs(title = "c) QQ Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") + 
  theme(legend.position = "none") 

# qqboxplot of residuals evaluated at posterior median 
packages <- c("dplyr", "ggplot2", "gridExtra","qqboxplot")
sapply(packages, require, character.only = TRUE)
ggplot(data = data.frame(y = emedian), aes(y = y)) + 
  geom_qqboxplot(notch = TRUE, varwidth = TRUE, reference_dist = "norm") + 
  ylab('Normal Scale') + 
  ggtitle('Hong Kong Data') + 
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) + 
  ylim(-3, 3) + 
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(size = 23)) 
ggsave("~/Dropbox/Apps/Overleaf/Handbook on Statistics of Extremes - Ch6/Ch6/LaTeX/figures/qqboxgev.pdf", width = 6)  

## COMPUTE WAIC
require(loo)
# (L is log-likelihood matrix with columns for observations and rows for iterations)
L3 <- matrix(NA, nrow = T - burn, ncol = n)
for(t in 1:(T - burn))
  for(i in 1:n)
    L3[t, i] <- log(dgev(y[i], x[i, ] %*% betamu[t, ], 
                         exp(x[i, ] %*% betasig[t, ]), x[i, ] %*% betash[t, ]))
model.eval3 <- waic(L3)
model.eval3$estimates
WAIC <- -2 * (model.eval3$estimates[1, 1])

comp <- loo_compare(model.eval1, model.eval2, model.eval3)
print(comp, digits = 1)
print(model.eval1, digits = 2)
print(model.eval2, digits = 2)
print(model.eval3, digits = 2)
