library(sf)
library(EnvStats)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(ggdist)
library(data.table)
library(pracma)
library(graphics)
library(zoo)
library(tidyverse)
library(ncdf4)
library(lattice)
library(RColorBrewer)
library(extRemes)
library(stars)
library(raster)
library(rnaturalearth)
library(mgcv)
library(gamlss)
load("Data/heatwaves.RData")

# 2023 heatwave
plot_regional_cut <- function() {
  fig <- ggplot(north_hemisphere_polygons) +
    geom_sf(aes(fill = region)) +
    labs(title = "Regional division of the North Hemisphere") +
    theme_classic(base_size=16) +
    scale_fill_manual(values = brewer.pal(9, "Paired")) +
    theme(plot.title = element_text(size = 20), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  ggsave(fig, filename = 'figures/regional-cut.png', width = 3000, height = 1200, unit = 'px', dpi = 200, bg = "white")
}

plot_regional_time_series <- function() {
  regional <- era5_maxtemp_1950_2023 %>% mutate(across(everything(), as.numeric))
  plot_list <- list()
  for (reg in names(regional)[-1]) {
    dat <- data.frame(year = regional$year, ts=regional[reg][,1], source="Observed")
    b <- gam(list(ts~s(year), ~1, ~1), data=dat, family=gevlss)
    mu <- fitted(b)[,1]
    rho <- fitted(b)[,2]
    xi <- fitted(b)[,3]
    fv <- mu + exp(rho)*(gamma(1-xi)-1)/xi
    shape <- round(xi[1],2)
    fitted <- paste("GEV ", "\U03Be", "=", shape, sep="")
    gev <- data.frame(year = regional$year, ts=fv, source=fitted)
    dd <- rbind(dat, gev)
    plot <- ggplot(data = dd, aes(x = year, colour=source, y=ts)) +
      geom_step()+
      scale_color_manual(guide = "legend", name="", values = setNames(c("black", "red"), c("Observed", fitted))) +
      #geom_line(data = regional %>% mutate(rolling_mean = zoo::rollmean(!!sym(reg), k = 10, fill = NA, align = "center")),
               # aes(y = rolling_mean), color = "red") +
      labs(title = gsub(".", " ", reg, fixed = TRUE), x = "year", y = "TxN (\u00B0C)") +
      theme_classic(base_size=16) +
      theme(plot.title = element_text(size = 20), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))+
      theme(legend.position="bottom")
    plot_list[[reg]] <- plot
  }
  plot_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)
  ggsave(filename = "figures/regional-time-series-gam.png", width = 4000, height = 1600, unit = 'px', dpi = 200, bg = "white", plot = plot_grid)
}

plot_examples_reference <- function() {
  plot_list <- list()
  datalist <- c(era5_us_midjuly2023, era5_china_midjuly2023, era5_europe_midjuly2023)
  Nlist <- c(18, 14, 7) # number of days for the rolling mean
  titlelist <- c("USA/Mexico", "China", "Southern Europe") # region names
  coastlines <- ne_coastline(returnclass = 'sf')
  countries <- ne_countries(returnclass = 'sf')
  for (i in 1:length(datalist)) {
    stars_obj <- datalist[[i]]$data
    lat <- datalist[[i]]$lat
    lon <- datalist[[i]]$lon
    title <- titlelist[i]
    n <- Nlist[[i]]
    countries$color <- "blanchedalmond"
    plot <- ggplot() +
      geom_sf(data = coastlines) +
      geom_sf(data = countries, aes(fill = factor(1)), fill = "blanchedalmond") +
      coord_sf(ylim = c(min(lat), max(lat)), xlim = c(min(lon), max(lon))) +
      geom_stars(data = stars_obj, na.action = na.omit) +
      labs(title = title, x = "longitude", y = "latitude") +
      scale_fill_gradientn(colors = (brewer.pal(9, "OrRd")),
                           guide = guide_colourbar(barwidth = 0.5, barheight = 6),
                           name = paste0(paste0("Tx", n), " (\u00B0C)")) +
      theme(panel.background = element_rect("lightskyblue"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), plot.title = element_text(size = 20))
    plot_list[[i]] <- plot
  }
  top_grid <- cowplot::plot_grid(plotlist = plot_list[1:2],
                                 labels = c("", ""),
                                 ncol = 2)
  plot <- cowplot::plot_grid(top_grid, plot_list[[3]],
                             labels = c("TxN on the second week of July 2023", ""),
                             label_size = 20,
                             label_fontface = "plain",
                             ncol = 1,
                             hjust = -0.5,
                             rel_heights = c(1.1, 1),
                             rel_widths = c(1.3, 1))
  ggsave(filename = "figures/examples-july-2023.png", plot = plot, bg = "white", width = 2200, height = 2000, unit = 'px', dpi = 200)
}

polynomial <- function(x, a, b, c) {
  return(a + b * x + c * x**2)
}

plot_regional_return_levels <- function(region = "Southern.Europe") {
  y <- era5_maxtemp_1950_2023[[region]]
  time <- era5_maxtemp_1950_2023[["year"]]
  time <- (time - min(time)) / (max(time) - min(time))
  logrange <- 10^seq(log10(1 + 1e-2), log10(10000), length.out = 50)
  # nonstationary fit
  gev_fit <- fevd(y, data.frame(y = y, time = time, time2 = time^2),
                  location.fun = ~time + time2)
  pars <- gev_fit$results$par
  data <- data.frame(value = y,
                     time = era5_maxtemp_1950_2023[["year"]],
                     location = polynomial(time, pars["mu0"], pars["mu1"], pars["mu2"]))
  # detrended series for the years of interest
  y_ref <- y - data$location + data$location[41]
  gev_fit_ref <- fevd(y_ref)
  ci_loc_ref <- ci(gev_fit_ref, type = "parameter", which.par = 1, method = "proflik", xrange = c(34, 36, 0.01),
                   nint = 200)
  y_now <- y - data$location + data$location[length(y)]
  gev_fit_now <- fevd(y_now)
  ci_loc_now <- ci(gev_fit_now, type = "parameter", which.par = 1, method = "proflik", xrange = c(36, 38, 0.01),
                   nint = 200)
  sorted_yref <- sort(y_ref)
  sorted_ynow <- sort(y_now)
  n <- length(sorted_ynow)
  real_rt <- 1 / (1 - seq(1, n) / n)
  # profile likelihood CI don't work so well with this package, so doing bootstrap-based in R
  ci_ref <- ci(gev_fit_ref, type = "return.level", return.period = logrange, do.ci = TRUE, method = "boot")
  ci_now <- ci(gev_fit_now, type = "return.level", return.period = logrange, do.ci = TRUE, method = "boot")
  df_ref <- data.frame(lower = as.numeric(ci_ref[, 1]), upper = as.numeric(ci_ref[, 3]), rl = as.numeric(ci_ref[, 2]), year = 1990, rp = logrange)
  df_now <- data.frame(lower = as.numeric(ci_now[, 1]), upper = as.numeric(ci_now[, 3]), rl = as.numeric(ci_now[, 2]), year = 2023, rp = logrange)
  data_rls <- rbind(df_ref, df_now)
  fig1 <- ggplot(data = data, mapping = aes(x = time)) +
    geom_point(aes(y = value), color = "black", pch = 'x', cex = 1) +
    geom_line(aes(y = location), color = "red") +
    geom_segment(aes(x = 1990, xend = 1990, y = ci_loc_ref[1], yend = ci_loc_ref[3]), color = "red") +
    geom_segment(aes(x = 2023, xend = 2023, y = ci_loc_now[1], yend = ci_loc_now[3]), color = "red") +
    labs(title = paste0("a) TxN 1950-2023 ", gsub(".", " ", region, fixed = TRUE)), x = "year", y = "TxN (\u00B0C)") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  fig2 <- ggplot() +
    geom_line(data = data_rls, aes(x = rp, y = rl, colour = factor(year))) +
    geom_line(data = data_rls, aes(x = rp, y = lower, colour = factor(year)), linetype = "dashed") +
    geom_line(data = data_rls, aes(x = rp, y = upper, colour = factor(year)), linetype = "dashed") +
    geom_point(aes(x = real_rt, y = sorted_yref), color = "royalblue", pch = 'x', cex = 3) +
    geom_point(aes(x = real_rt, y = sorted_ynow), color = "red", pch = 'x', cex = 3) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_manual(name = NULL, guide = "legend", values = c("1990" = "royalblue", "2023" = "red")) +
    labs(title = paste0("b) Estimated return levels in ", gsub(".", " ", region, fixed = TRUE)), x = "return period (years)", y = "TxN (\u00B0C)") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  plot <- cowplot::plot_grid(fig1, fig2, nrow = 2)
  ggsave(filename = "figures/southern-europe-rls.png", plot = plot, width = 2400, height = 1600, unit = 'px', bg = "white", dpi = 200)
  return(plot)
}

# Phalodi
jodhpur_name <- "IN019180500"
jod_maxtemp <- drop_na(phalodi_maxtemp[1:46, c(1, 3)])

plot_phalodi_rls <- function() {
  y <- jod_maxtemp
  l <- length(y[[jodhpur_name]])
  yinc <- y[[jodhpur_name]]
  yex <- yinc[1:l - 1]
  logrange <- 10^seq(log10(1 + 1e-2), log10(1000000), length.out = 50)
  level <- 48.8
  sorted_y <- sort(yinc)
  n <- length(sorted_y)
  real_rt <- 1 / (1 - seq(1, n) / n)
  gev_fit <- fevd(yinc)
  gev_fit_ex <- fevd(yex)
  cis <- ci(gev_fit, type = "return.level", return.period = logrange, do.ci = TRUE, method = "boot")
  cis_ex <- ci(gev_fit_ex, type = "return.level", return.period = logrange, do.ci = TRUE, method = "boot")
  df <- data.frame(lower = as.numeric(cis[, 1]), upper = as.numeric(cis[, 3]), rl = as.numeric(cis[, 2]), rp = logrange, method = "Including")
  df_ex <- data.frame(lower = as.numeric(cis_ex[, 1]), upper = as.numeric(cis_ex[, 3]), rl = as.numeric(cis_ex[, 2]), rp = logrange, method = "Excluding")
  data_rls <- rbind(df, df_ex)
  fig1 <- ggplot(data = y, mapping = aes(x = YEAR, y = yinc)) +
    geom_point(color = "black", pch = 'x', cex = 3) +
    labs(title = "a) TXx Phalodi May-June 1944-2016", x = "year", y = "TXx (\u00B0C)") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  fig2 <- ggplot() +
    geom_line(data = data_rls, aes(x = rp, y = rl, colour = factor(method))) +
    geom_line(data = data_rls, aes(x = rp, y = lower, colour = factor(method)), linetype = "dashed") +
    geom_line(data = data_rls, aes(x = rp, y = upper, colour = factor(method)), linetype = "dashed") +
    geom_hline(yintercept = level, color = "goldenrod") +
    geom_point(aes(x = real_rt, y = sorted_y), color = "black", pch = 'x', cex = 3) +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(ylim = c(43, 53)) +
    scale_color_manual(name = NULL, guide = "legend", values = c("Excluding" = "royalblue", "Including" = "red")) +
    labs(title = paste0("b) Estimated return levels in Phalodi"), x = "return period (years)", y = "TXx (\u00B0C)") +
    theme_classic(base_size=16) +
    theme(plot.title = element_text(size = 18), panel.background = element_rect(fill = "#eeeeee",
                                colour = "#eeeeee",
                                size = 0.5, linetype = "blank"))
  plot <- cowplot::plot_grid(fig1, fig2, ncol = 2, rel_widths = c(1, 1.4))
  ggsave(filename = "figures/phalodi-time-series.png", plot = plot, width = 3136, height = 1080, unit = 'px', bg = "white")
  return(plot)
}

plot_phalodi_rls()
plot_examples_reference()
plot_regional_time_series()
plot_regional_return_levels()
plot_regional_cut()