library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(graphics)
library(data.table)
library(POT)
library(extremogram)
library(readr)
library(tidyverse)
library(lubridate)
library(ggdist)
load("Data/Venice.RData")
venice_hourly_over80_1982_2023$time <- as.POSIXct(venice_hourly_over80_1982_2023$time)

plot_tenmax <- function() {
  fig <- ggplot() +
    geom_point(data = tidyr::gather(ten_yearly_max_venice_1887_2023, max_n, level_cm, -year), aes(x = year, y = level_cm), color = "black", pch = 'x', cex = 3) +
    labs(x = "year", y = "level(cm)", title = "a) Ten largest recorded tides per year between 1887 and 2023") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  print(fig)
  ggsave("figures/tenmax-venice.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900,)
  return(fig)
}

plot_number_exceedances_above_80cm <- function() {
  nb <- dplyr::rename(venice_hourly_over80_1982_2023, level_cm = data)
  nb <- nb %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::mutate(rolling_mean = zoo::rollmean(count, k = 5, fill = NA))
  fig <- ggplot(nb, aes(x = year, y = rolling_mean)) +
    geom_line(color = "black") +
    geom_point(color = "black", , pch = 'x', cex = 3) +
    labs(x = "year", y = "", title = "c) Number of exceedances above 80cm each year") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  print(fig)
  ggsave("figures/nb-ex-80cm-venice.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900,)
  return(fig)
}

plot_exceedances_above_80cm <- function() {
  venice_hourly_over80_1982_2023$MOSE <- as.integer(venice_hourly_over80_1982_2023$MOSE)
  fig <- ggplot(data = venice_hourly_over80_1982_2023, mapping = aes(x = time, y = data, colour = factor(MOSE))) +
    geom_point(pch = "x", cex = 3) +
    scale_color_manual(name = "MOSE", guide = "legend", values = c("0" = "black", "1" = "red"), labels = c(NA, NA)) +
    labs(x = "time", y = "level(cm)", title = "b) Recorded exceedances above 80cm between 1982 and 2023") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  print(fig)
  ggsave("figures/ex-80cm-venice.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900,)
  return(fig)
}

plot_extremogram <- function(threshold = 80) {
  data <- venice_hourly_over80_1982_2023 %>%
    mutate(timedelta = time - time[1]) %>%
    mutate(days_since_start = as.numeric(timedelta) / (24 * 60 * 60)) %>%
    arrange(days_since_start)
  data_extremogram <- data %>%
    filter(data >= threshold) %>%
    dplyr::select(days_since_start, data)
  data_extremogram <- data_extremogram %>%
    mutate(iat = days_since_start - lag(days_since_start, default = 0))
  h_range <- seq(1, 100, by = 2)
  indices_exceedances <- as.integer(data_extremogram$days_since_start)
  extremogram_realized <- extremogram_loop(h_range, indices_exceedances)
  extremogram_realized <- as.vector(extremogram_realized)
  df <- data.frame(h = h_range, extremogram = extremogram_realized)
  fig <- ggplot(data = df, aes(x = h, y = extremogram)) +
    geom_line(color = "black") +
    geom_point(color = "black", pch = 'x', cex = 3) +
    ggtitle("a) Extremogram for tide levels above u=80cm") +
    xlab("h") +
    ylab(expression(pi[h](u))) +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  print(fig)
  ggsave("figures/extremogram-venice.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900,)
  return(fig)
}

extremogram_loop <- function(h_range, indices_exceedances) {
  counts <- vector()
  for (h in h_range) {
    counts <- append(counts, length(intersect(indices_exceedances + h, indices_exceedances)) / length(indices_exceedances))
  }
  return(counts)
}

plot_cumulative_nb_exceedances <- function() {
  sorted <- venice_hourly_over80_1982_2023 %>% arrange(time)
  days_since_start <- difftime(sorted$time, sorted$time[1], units = 'days')
  iat <- c(0, diff(days_since_start))
  sorted$bool <- 1
  cum_nb_ex <- data.frame(nb_ex = cumsum(sorted$bool), time = sorted$time,
                          label = "Observed", ls = "value", name = "obs")
  exp <- data.frame(nb_ex = days_since_start * (1 / mean(iat)),
                    time = sorted$time, label = "Exponential IAT", ls = "value", name = "exp")
  params <- list()
  for (i in 1:1000) {
    s <- sample(iat, length(iat), replace = TRUE)
    params[[i]] <- 1 / mean(s)
  }
  lower <- data.frame(nb_ex = quantile(as.numeric(params), 0.05) * days_since_start,
                      time = sorted$time, label = "Exponential IAT", ls = "95% CI", name = "lower")
  upper <- data.frame(nb_ex = quantile(as.numeric(params), 0.95) * days_since_start,
                      time = sorted$time, label = "Exponential IAT", ls = "95% CI", name = "upper")
  data_cne <- rbind(cum_nb_ex, exp, lower, upper)
  data_cne <- data_cne %>% arrange(time)
  fig <- ggplot(data = data_cne, aes(x = time, y = nb_ex, colour = label, linetype = ls, group = name)) +
    geom_line() +
    scale_color_manual(name = c("source"), guide = "legend", values = c("Observed" = "black", "Exponential IAT" = "red")) +
    scale_linetype_manual(name = c("type"), guide = "legend", values = c("value" = "solid", "95% CI" = "dashed")) +
    labs(x = "time", y = "number of exceedances", title = "b) Cumulative number of exceedances above 80cm") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  print(fig)
  ggsave("figures/cum-nb-exceedances-venice.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900)
  return(fig)
}

plot_mean_cluster_size <- function() {
  sorted <- venice_hourly_over80_1982_2023 %>% arrange(time)
  data <- sorted %>%
    mutate(timedelta = time - time[1]) %>%
    mutate(days_since_start = as.numeric(timedelta) / (24 * 60 * 60)) %>%
    mutate(time_norm = days_since_start / max(days_since_start))
  data_extremal_index <- data %>%
    dplyr::select(days_since_start, data, time_norm) %>%
    mutate(iat_normed = c(0, diff(time_norm))) %>%
    mutate(iat = c(0, diff(days_since_start)))
  iat <- as.numeric(data_extremal_index$iat)
  block_sizes <- seq(1, 50, by = 2)
  empirical_mean_cluster_size <- loop_extremal_index(data_extremal_index$days_since_start, block_sizes)
  empirical_mean_cluster_size <- data.frame(bs = block_sizes, mcs = empirical_mean_cluster_size, source = "Observed", type = "value", name = "obs")
  params <- list()
  for (i in 1:1000) {
    s <- sample(iat, length(iat), replace = TRUE)
    exp_param <- 1 / mean(s)
    exp_sample <- cumsum(rexp(n = length(data_extremal_index$days_since_start), rate = exp_param))
    to_concat <- data.frame(loop_extremal_index(exp_sample, block_sizes))
    names(to_concat) = i
    params[[i]] <- to_concat
  }
  total <- data.frame(t(bind_cols(params)))
  mean_exp <- setDT(total)[, lapply(.SD, mean, na.rm = TRUE)]
  quants <- setDT(total)[, lapply(.SD, quantile, probs = c(.05, .95), na.rm = TRUE)]
  lower <- data.frame(bs = block_sizes, mcs = as.numeric(t(quants[1,])), source = "Exponential IAT", type = "95% CI", name = "lower")
  upper <- data.frame(bs = block_sizes, mcs = as.numeric(t(quants[2,])), source = "Exponential IAT", type = "95% CI", name = "upper")
  exp <- data.frame(bs = block_sizes, mcs = as.numeric(t(mean_exp)), source = "Exponential IAT", type = "value", name = "mean")
  data_mcs <- rbind(empirical_mean_cluster_size, exp, lower, upper)
  fig <- ggplot(data = data_mcs, aes(x = bs, y = mcs, colour = source, linetype = type, group = name)) +
    geom_line() +
    geom_point(pch = 'x', cex = 2) +
    scale_color_manual(guide = "legend", values = c("Observed" = "black", "Exponential IAT" = "red")) +
    scale_linetype_manual(guide = "legend", values = c("value" = "solid", "95% CI" = "dashed")) +
    labs(x = "block size r (days)", title = "c) Mean cluster size") +
    ylab(expression(theta[r]^-1 ~ (u))) +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  print(fig)
  ggsave("figures/mean-cluster-size-venice.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900,)
  return(fig)
}

loop_extremal_index <- function(indices_exceedances, block_sizes) {
  local_count <- c()
  for (bs in block_sizes) {
    hp <- hist(indices_exceedances,
               breaks = seq(0, bs * ceiling(max(indices_exceedances) / bs), by = bs),
               plot = FALSE)$counts
    hp <- hp[hp > 0]
    local_count <- c(local_count, mean(hp))
  }
  local_count <- set_names(local_count, block_sizes)
  return(local_count)
}

plot_tenmax()
plot_exceedances_above_80cm()
plot_number_exceedances_above_80cm()
plot_extremogram()
plot_cumulative_nb_exceedances()
plot_mean_cluster_size()