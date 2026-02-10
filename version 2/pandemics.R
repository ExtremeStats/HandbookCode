library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(ggdist)
library(data.table)
library(pracma)
library(graphics)
library(gridExtra)
library(ggpubr)
library(extRemes)
source("Code/utils.R", local = TRUE)
load("Data/pandemics.RData")
quants <- seq(1e-4, 1 - 1e-4, length.out = 100)

plot_pot_and_survival <- function() {
  u <- 2.4
  y <- log(1+10^5 * pandemic_deaths_year["prop_population"])
  meplot <- plot_mean_excess(y, threshold = u)
  potplot <- plot_pot(y$prop_population, threshold = u, x=pandemic_deaths_year$year)
  fit <- fevd(y$prop_population, threshold = u, type = 'GP')
  bootstrap_samples <- ci(fit, return.period = 100, R=5000, method = 'boot', return.samples = TRUE)[, 1:2]
  rl_list <- list()
  shapes <- list()
  upbound <- list()
  for (i in 1:size(bootstrap_samples)[1]) {
    scale <- bootstrap_samples[i, 1]
    shape <- bootstrap_samples[i, 2]
    rls <- rlevd(1 / (1 - quants), type = "GP", threshold = u, scale = scale, shape = shape, rate = 1 / 365.25)
    df <- data.frame(value = rls)
    names(df) <- i
    rl_list[[i]] <- df
    shapes[[i]] <- shape
    if (shape<0){
    upbound[[i]] <- u-scale/shape} else {upbound[[i]] <- NA}
  }
  d<-data.frame(upbound)
  j <- d[!is.na(d)] # then get percentiles
  print(quantile((exp(j)/10^5), 0.3))
  print(median(exp(j)/10^5))
  print(quantile((exp(j)/10^5), 0.6))
  total <- data.frame(t(bind_cols(rl_list)))
  cis <- setDT(total)[, lapply(.SD, quantile, probs = c(0.05, 0.95), na.rm = TRUE), .SDcols = names(total)]
  theo <- extRemes::qevd(quants, type = "GP", threshold = u, scale = fit$results$par["scale"], shape = fit$results$par["shape"])
  df <- data.frame(lower = as.numeric(cis[1,]), upper = as.numeric(cis[2,]), rl = theo, p = 1 - quants)
  excesses <- y$prop_population[y$prop_population >= u]
  ecdf_y <- ecdf(excesses)
  empirical <- sapply(theo, FUN = function(x) 1 - ecdf_y(x))
  survival <- ggplot() +
    geom_line(aes(x = theo, y = 1 - quants), color = "slategrey") +
    geom_line(data = df, aes(x = upper, y = p), color = "slategrey", linetype = "dashed") +
    geom_line(data = df, aes(x = lower, y = p), color = "slategrey", linetype = "dashed") +
    geom_point(aes(x = theo, y = empirical), color = "black", pch = 16, cex = 0.7) +
    labs(x = "x", y = "Pr(M>x)") +
    xlim(NA, 15)+
    xscale("log")+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 17),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  plot <- grid.arrange(meplot, potplot, survival, nrow = 1)
  ggsave(filename = "figures/pandemic-data-pot.png", plot = plot, bg = "white", width = 3800, height = 1080, unit = 'px')
  return(plot)
}

plot_pot_and_survival()