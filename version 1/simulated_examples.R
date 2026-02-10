library(evd)
library(mev)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(ggdist)
library(data.table)
library(pracma)
library(graphics)
library(EnvStats)
library(gridExtra)
library(ggpubr)
source("Code/utils.R", local=TRUE)

set.seed(60)
block_size <- 50
nb_blocks <- 100
n <- block_size * nb_blocks
quants <- 1 - 1 / logspace(0.01, 1.75, 25)

# Gamma
sample_gamma <- rgamma(n, 1)
theo_params_gamma <- cbind(log(block_size), 1, 0)

# Uniform
sample_unif <- runif(n)
theo_params_unif <- cbind(block_size / (1 + block_size), 1 / (block_size + 1) * sqrt(block_size / (block_size + 2)), -1)

# Pareto
alpha <- 2
sample_pareto <- rpareto(n, 1, 2)
c <- (block_size)**(1 / alpha)
xi <- 1 / alpha
theo_params_pareto <- cbind(c, c * xi, xi)

plot_density <- function(sample, theo_params) {
  params <- theo_params
  blocks <- split(sample, f = rep(1:(length(sample) %/% block_size), each = block_size))
  block_max <- sapply(blocks, FUN = max)
  fig <- ggplot() +
    geom_histogram(aes(x = block_max, y = after_stat(density)), bins = 20, color = "red", fill = "red", alpha = 0.3) +
    geom_line(aes(x = block_max, y = dgev(block_max, loc = params[1], scale = params[2], shape = params[3])), color = "grey") +
    labs(title = "b) Density function", x = "x", y = "pdf(x)") +
    theme_classic(base_size=16)+
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  return(fig)
}

plot_survival <- function(theo, theo_origin, empirical) {
  fig <- ggplot() +
    geom_line(aes(x = theo, y = 1 - quants), color = "grey", alpha = 0.7) +
    geom_line(aes(x = theo, y = theo_origin), color = "black", alpha = 0.7) +
    geom_point(aes(x = theo, y = empirical), color = "red", pch = 'x', cex = 3) +
    scale_y_log10() +
    scale_x_log10() +
    labs(title = "c) Survival function", x = "x", y = "Pr(M>x)") +
    theme_classic(base_size=16)+
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  return(fig)
}

plot_survival_gamma <- function() {
  params <- theo_params_gamma
  blocks <- split(sample_gamma, f = rep(1:(length(sample_gamma) %/% block_size), each = block_size))
  block_max <- sapply(blocks, FUN = max)
  theo <- qgev(quants, loc = params[1], scale = params[2], shape = params[3])
  theo_origin <- 1 - pgamma(theo, shape = 1)
  ecdf_block <- ecdf(block_max)
  empirical <- sapply(theo, FUN = function(x) 1 - ecdf_block(x))
  fig <- plot_survival(theo, theo_origin, empirical)
  return(fig)
}

plot_survival_pareto <- function() {
  params <- theo_params_pareto
  blocks <- split(sample_pareto, f = rep(1:(length(sample_pareto) %/% block_size), each = block_size))
  block_max <- sapply(blocks, FUN = max)
  theo <- qgev(quants, loc = params[1], scale = params[2], shape = params[3])
  theo_origin <- 1 - ppareto(theo, shape = 2, location = 1)
  ecdf_block <- ecdf(block_max)
  empirical <- sapply(theo, FUN = function(x) 1 - ecdf_block(x))
  fig <- plot_survival(theo, theo_origin, empirical)
  return(fig)
}

plot_survival_uniform <- function() {
  params <- theo_params_unif
  blocks <- split(sample_unif, f = rep(1:(length(sample_unif) %/% block_size), each = block_size))
  block_max <- sapply(blocks, FUN = max)
  theo <- qgev(quants, loc = params[1], scale = params[2], shape = params[3])
  theo_origin <- 1 - punif(theo)
  ecdf_block <- ecdf(block_max)
  empirical <- sapply(theo, FUN = function(x) 1 - ecdf_block(x))
  fig <- plot_survival(theo, theo_origin, empirical)
  return(fig)
}

plot_block_maxima <- function(sample) {
  blocks <- split(sample, rep(1:(length(sample) %/% block_size), each = block_size))
  block_max <- sapply(blocks, max)
  indices <- sapply(1:nb_blocks, function(i) {
    v <- blocks[[i]]
    return(which.max(v) + (i - 1) * block_size)
  })
  full_sample <- as.data.frame(sample)
  full_sample$id <- 1:nrow(full_sample)
  fig <- ggplot() +
    geom_point(aes(x = full_sample$id, y = full_sample$sample), color = "black", pch = 'x', cex = 3) +
    geom_point(aes(x = indices, y = block_max), color = "red", pch = 'x', cex = 3) +
    labs(title = "a) Block maxima (r=50)", x = "", y = "") +
    theme_classic(base_size=16)+
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  return(fig)
}

plot_survival_pot_pareto <- function(threshold = 7) {
  sigma_u <- c * xi + xi * (c - threshold)
  theo_params <- cbind(threshold, sigma_u, xi)
  full_sample <- as.data.frame(sample_pareto)
  full_sample$id <- 1:nrow(full_sample)
  pot <- full_sample$sample[full_sample$sample > threshold]
  theo <- qgpd(quants, loc = theo_params[1], scale = theo_params[2], shape = theo_params[3])
  theo_origin <- 1 - ppareto(theo, shape = 2, location = 1)
  ecdf_block <- ecdf(pot)
  empirical <- sapply(theo, FUN = function(x) 1 - ecdf_block(x))
  fig <- plot_survival(theo, theo_origin, empirical)
  return(fig)
}

# Pareto
fig1 <- plot_block_maxima(sample_pareto)
fig2 <- plot_density(sample_pareto, theo_params_pareto)
fig3 <- plot_survival_pareto()
pareto <- grid.arrange(fig1, fig2, fig3, nrow = 1, top = text_grob("Pareto distribution with shape=2", size=20))
ggsave('figures/example-xi-pos.png', pareto, width = 3600, height = 1080, unit = 'px')

# POT Pareto
fig1 <- plot_mean_excess(sample_pareto)
fig2 <- plot_pot(sample_pareto)
fig3 <- plot_survival_pot_pareto()
paretopot <- grid.arrange(fig1, fig2, fig3, nrow = 1, top = text_grob("Pareto distribution with shape=2", size=20))
ggsave('figures/example-xi-pos-pot.png', paretopot, width = 3600, height = 1080, unit = 'px')

# Gamma
fig1 <- plot_block_maxima(sample_gamma)
fig2 <- plot_density(sample_gamma, theo_params_gamma)
fig3 <- plot_survival_gamma()
gamma <- grid.arrange(fig1, fig2, fig3, nrow = 1, top = text_grob("Exponential distribution with shape=1", size=20))
ggsave('figures/example-xi-null.png', gamma, width = 3600, height = 1080, unit = 'px')

# Uniform
fig1 <- plot_block_maxima(sample_unif)
fig2 <- plot_density(sample_unif, theo_params_unif)
fig3 <- plot_survival_uniform()
unif <- grid.arrange(fig1, fig2, fig3, nrow = 1, top = text_grob("Standard uniform distribution", size=20))
ggsave('figures/example-xi-neg.png', unif, width = 3600, height = 1080, unit = 'px')