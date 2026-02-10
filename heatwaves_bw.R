library(sf)
library(EnvStats)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggpattern)
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
load("../Data/heatwaves.RData")

# 2023 heatwave
plot_regional_cut <- function() {
  regions <- st_make_valid(north_hemisphere_polygons, )
  valid_regions <- filter(regions, st_is_valid(regions))
  invalid <- filter(regions, !st_is_valid(regions))
  x<-filter(invalid, region %in% c("Middle East", "North-Eastern Asia", "South-Eastern Asia"))
  merged_region <- valid_regions %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry))
  df <- bind_rows(merged_region, x)
  fig <- ggplot(df, aes(fill=region, pattern = region, pattern_angle = region)) +
  geom_sf_pattern(pattern_color="black",
                  pattern_density=0.01, pattern_spacing = 0.025) +
    theme_classic(base_size=11)+ 
    scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "wave", "circle", "crosshatch", "stripe", "weave", "wave")) +
  scale_fill_manual(values = brewer.pal(9, "Greys")[2:10]) +
    theme(axis.text = element_text(size = 10), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  ggsave(fig, filename = 'figures/regional-cut.png', width = 2200, height = 900, unit = 'px', dpi = 200, bg = "white")
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
    fitted <- paste0("\u03be", "=", "\u2212",-shape, sep="") #\u03be+0302
    gev <- data.frame(year = regional$year, ts=fv, source=fitted)
    dd <- rbind(dat, gev)
    plot <- ggplot(data = dd, aes(x = year, colour=source, y=ts)) +
       scale_color_manual(guide = "legend", name="", values = setNames(c(NA,"black"), c(NA, fitted))) +
      geom_step()+
      labs(title = gsub(".", " ", reg, fixed = TRUE), x = "", y = "TxN (\u00B0C)") +
      theme_classic(base_size=11) +
      theme(axis.text = element_text(size = 10),
            plot.title = element_text(size = 11),
            panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "blank"),
            legend.position=c(0.8,0.2))
    plot_list[[reg]] <- plot
  }
  plot_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)
  ggsave(filename = "figures/regional-time-series-gam.png", width = 2000, height = 1000, unit = 'px', dpi = 200, bg = "white", plot = plot_grid)
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
    n <- Nlist[[i]]
    countries$color <- "blanchedalmond"
    plot <- ggplot() +
      geom_sf(data = coastlines) +
      geom_sf(data = countries, aes(fill = factor(1)), fill = "blanchedalmond") +
      coord_sf(ylim = c(min(lat), max(lat)), xlim = c(min(lon), max(lon))) +
      geom_stars(data = stars_obj, na.action = na.omit) +
      labs(x = "longitude", y = "latitude") +
      scale_fill_gradientn(colors = (brewer.pal(9, "OrRd")),
                           guide = guide_colourbar(barwidth = 0.5, barheight = 6),
                           name = paste0(paste0("Tx", n), " (\u00B0C)")) +
      theme_classic(base_size=11)+ 
      theme(axis.text = element_text(size = 10),
        panel.background = element_rect("lightskyblue"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
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
  ggsave(filename = "figures/examples-july-2023.png", plot = plot, bg = "white", width = 2000, height = 1000, unit = 'px', dpi = 200)
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
    geom_line(aes(y = location), color = "slategrey") +
    annotate("segment",x = 1990, xend = 1990, 
             y = ci_loc_ref[1], yend = ci_loc_ref[3], color = "slategrey") + 
    annotate("segment",x = 2023, xend = 2024, 
           y = ci_loc_now[1], yend = ci_loc_now[3], color = "slategrey") +
    geom_point(aes(y = value), color = "black", pch = 16, cex = 2) +
    labs(x = "", y = "TxN (\u00B0C)") +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  fig2 <- ggplot() +
    geom_line(data = data_rls, aes(x = rp, y = rl, colour = factor(year))) +
    geom_line(data = data_rls, aes(x = rp, y = lower, colour = factor(year)), linetype = "dashed") +
    geom_line(data = data_rls, aes(x = rp, y = upper, colour = factor(year)), linetype = "dashed") +
    geom_point(aes(x = real_rt, y = sorted_yref), color = "black",pch = 16, cex = 2) +
    geom_point(aes(x = real_rt, y = sorted_ynow), color = "slategrey", pch = 16, cex = 2) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_manual(name = NULL, guide = "legend", values = c("1990" = "black", "2023" = "slategrey")) +
    labs(x = "Return period (years)", y = "TxN (\u00B0C)") +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"),
          legend.position = c(0.8,0.2))
  plot <- cowplot::plot_grid(fig1, fig2, ncol = 2)
  ggsave(filename = "figures/southern-europe-rls.png", plot = plot, width = 2000, height = 1000, unit = 'px', bg = "white", dpi = 200)
  return(plot)
}

# Phalodi
jodhpur_name <- "IN019180500"
jod_maxtemp <- drop_na(phalodi_maxtemp[1:46, c(1, 3)])

plot_phalodi_rls <- function() {
  y <- jod_maxtemp
  l <- length(y[[jodhpur_name]])
  yinc <- y[[jodhpur_name]]
  yex <- yinc[1:(l - 1)]
  logrange <- 10^seq(log10(1 + 1e-2), log10(1000), length.out = 50)
  level <- 48.8
  sorted_y <- sort(yinc)
  n <- length(sorted_y)
  real_rt <- 1 / (1 - (seq(1, n)-0) / (n+0))
  gev_fit <- fevd(yinc)
  gev_fit_ex <- fevd(yex)
  cis <- ci(gev_fit, type = "return.level", return.period = logrange, do.ci = TRUE, method = "boot", R=1000)
  cis_ex <- ci(gev_fit_ex, type = "return.level", return.period = logrange, do.ci = TRUE, method = "boot", R=1000)
  df <- data.frame(lower = as.numeric(cis[, 1]), upper = as.numeric(cis[, 3]), rl = as.numeric(cis[, 2]), rp = logrange, method = "Including")
  df_ex <- data.frame(lower = as.numeric(cis_ex[, 1]), upper = as.numeric(cis_ex[, 3]), rl = as.numeric(cis_ex[, 2]), rp = logrange, method = "Excluding")
  data_rls <- rbind(df, df_ex)
  fig1 <- ggplot(data = y, mapping = aes(x = YEAR, y = yinc)) +
    theme_classic(base_size=11) +
    geom_point(color = "black", pch = 16, cex = 2) +
    labs(x = "", y = "TXx (\u00B0C)") +
    scale_y_continuous(breaks=c(44,46,48,50,52)) +
    theme(axis.text = element_text(size = 10), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  fig2 <- ggplot() +
    geom_line(data = data_rls, aes(x = rp, y = rl, colour = factor(method))) +
    geom_line(data = data_rls, aes(x = rp, y = lower, colour = factor(method)), linetype = "dashed") +
    geom_line(data = data_rls, aes(x = rp, y = upper, colour = factor(method)), linetype = "dashed") +
    geom_hline(yintercept = level, color = "grey", linetype="dotted") +
    geom_point(aes(x = real_rt, y = sorted_y), color = "black", pch = 16, cex = 2) +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(ylim = c(43, 53)) +
    scale_y_continuous(breaks=c(44,46,48,50,52)) +
    scale_color_manual(name = NULL, guide = "legend", values = c("Excluding" = "black", "Including" = "slategrey")) +
    labs(x = "Return period (years)", y = "TXx (\u00B0C)") +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10), legend.position=c(0.8,0.2))
  plot <- cowplot::plot_grid(fig1, fig2, ncol = 2)
  ggsave(filename = "figures/phalodi-time-series.png", plot = plot, width = 2800, height = 1100, unit = 'px', bg = "white")
  return(plot)
}

plot_regional_time_series()
plot_regional_cut()
plot_phalodi_rls()
plot_regional_return_levels()
