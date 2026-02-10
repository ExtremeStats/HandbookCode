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
padova<-read.csv("Data/Padova_precipitation.csv")
padova <- padova |> mutate(date = make_date(YEAR, MONTH, DAY))
padova$date <- as.POSIXct(padova$date)


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

plot_exceedances_above_80cm()
plot_number_exceedances_above_80cm()
plot_extremogram()
plot_cumulative_nb_exceedances()
plot_mean_cluster_size()