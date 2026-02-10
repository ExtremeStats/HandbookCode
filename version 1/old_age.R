library(ggplot2)
library(ggdist)
library(data.table)
library(pracma)
library(gridExtra)
library(graphics)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(bshazard)
library(flexsurv)
load("Data/old_age.RData")

plot_Lexis <- function() {
  sample_size <- 200
  c1_list <- c(as.Date("2003-01-01"), as.Date("1968-01-01")) # from IDL metadata
  c2_list <- c(as.Date('2014-12-31'), as.Date('2020-12-31'))
  title_list <- c("Austria", "United Kingdom")
  data_list <- list(idl_aut, idl_uk)
  plot_list <- list()
  for (i in 1:2) {
    title <- title_list[[i]]
    c1 <- c1_list[[i]]
    c2 <- c2_list[[i]]
    d <- data_list[[i]]
    midyear <- c1 + (c2 - c1) / 2
    nb_ind_max_before_midyear <- min(nrow(d[d$deathday < midyear,]), sample_size / 2)
    nb_ind_max_after_midyear <- min(nrow(d[d$deathday >= midyear,]), sample_size - nb_ind_max_before_midyear)
    rds <- dplyr::bind_rows(d[d$deathday < midyear,] %>% dplyr::sample_n(nb_ind_max_before_midyear),
                            d[d$deathday >= midyear,] %>% dplyr::sample_n(nb_ind_max_after_midyear))
    rds$date_u0 <- as.Date(rds$date_u0)
    rds$deathday <- as.Date(rds$deathday)
    plot <- ggplot(data = rds) +
      geom_segment(aes(x = date_u0, xend = deathday, y = u0, yend = excess_over_u0), color = "black") +
      geom_vline(xintercept = c1, color = "red", linetype = "dashed") +
      annotate("text", x = c1, y = 0, angle = 90, label = "c1\n", color = "red") +
      geom_vline(xintercept = c2, color = "red", linetype = "dashed") +
      annotate("text", x = c2, y = 0, angle = 90, label = "c2\n", color = "red") +
      labs(title = title, x = "time", y = "excess lifetime over 105 years") +
      theme_classic(base_size = 16) +
      scale_x_date(limit = c(c1, c2)) +
      ylim(0, NA) +
      theme(plot.title = element_text(size = 18),
            panel.background = element_rect(fill = "#eeeeee",
                                            colour = "#eeeeee",
                                            size = 0.5, linetype = "blank"))
    plot_list[[i]] <- plot
  }
  plot_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
  title <- text_grob("Lexis diagram", size = 20)
  plot_grid_with_title <- cowplot::plot_grid(title, plot_grid, nrow = 2, rel_heights = c(0.1, 1))
  ggsave(filename = "figures/Lexis-idl-samples.png", width = 3000, height = 1200, unit = 'px', dpi = 200, bg = "white", plot = plot_grid_with_title)
}

plot_survival <- function() {
  title_list <- c("Austria", "United Kingdom")
  data_list <- list(idl_aut, idl_uk)
  plot_list <- list()
  for (i in 1:2) {
    title <- title_list[[i]]
    d <- data_list[[i]]
    fit <- survfit(Surv(time = AGEDAYS, event = as.logical(Status)) ~ 1, data = d)
    plot <- ggsurvplot(fit, lw = 1, data = d,
                       conf.int = TRUE, palette = "black",
                       ggtheme = theme_classic(base_size = 16) + theme(,
                         panel.background = element_rect(fill = "#eeeeee",
                                                         colour = "#eeeeee",
                                                         size = 0.5, linetype = "blank")),
                       legend.labs = NA,
                       break.time.by = 365.25,
                       xscale = 365.25, ylab = expression(hat(S)(t)),
                       xlab = "time t (years)",
                       linetype = c("solid"),
                       conf.int.style = "step",
                       font.title = c(size = 18),
                       xlim = c(min(d$AGEDAYS) + 100, max(d$AGEDAYS)),
                       title = title) # empirical
    plot_list[[i]] <- plot
  }
  plot_grid <- arrange_ggsurvplots(plot_list, nrow = 2, ncol = 1,
                                   title = text_grob("Estimated probability of survival", size = 20))
  ggsave(filename = "figures/survival-idl-samples.png", width = 1600, height = 1600, unit = 'px', dpi = 200, bg = "white", plot = plot_grid)
}

# same as previous function for cumulative hazard instead of survival
plot_cum_hazard <- function() {
  title_list <- c("Austria", "United Kingdom")
  data_list <- list(idl_aut, idl_uk)
  plot_list <- list()
  for (i in 1:2) {
    title <- title_list[[i]]
    d <- data_list[[i]]
    fit <- survfit(Surv(time = AGEDAYS, event = as.logical(Status)) ~ 1, data = d)
    plot <- ggsurvplot(fit, fun = "cumhaz", lw = 1, data = d,
                       conf.int = TRUE, palette = "black",
                       yscale = "log10",
                       ggtheme = theme_classic(base_size = 16) + theme(panel.background = element_rect(fill = "#eeeeee",
                                                                                                       colour = "#eeeeee",
                                                                                                       size = 0.5, linetype = "blank")),
                       legend.labs = NA,
                       break.time.by = 365.25,
                       xscale = 365.25, ylab = expression(hat(H)(t)),
                       xlab = "time t (years)",
                       linetype = c("solid"),
                       conf.int.style = "step",
                       font.title = c(size = 18),
                       xlim = c(min(d$AGEDAYS) + 100, max(d$AGEDAYS)),
                       title = title) # empirical
    plot_list[[i]] <- plot
  }
  plot_grid <- arrange_ggsurvplots(plot_list, nrow = 2, ncol = 1, title = text_grob("Estimated cumulative hazard", size = 20))
  ggsave(filename = "figures/hazard-idl-samples.png", width = 1600, height = 1600, unit = 'px', dpi = 200, bg = "white", plot = plot_grid)
}

plot_hazard <- function() {
  title_list <- c("Austria", "United Kingdom")
  data_list <- list(idl_aut, idl_uk)
  fit <- list()
  for (i in 1:2) {
    title <- title_list[[i]]
    d <- data_list[[i]]
    fit[[i]] <- bshazard(Surv(time = AGEDAYS, event = as.logical(Status)) ~ 1, data = d)
  }
  ausplot <- ggplot(mapping = aes(x = fit[[1]]$time / 365.25)) +
    geom_line(aes(y = fit[[1]]$hazard)) +
    geom_line(aes(y = fit[[1]]$lower.ci), linetype = 'dashed') +
    geom_line(aes(y = fit[[1]]$upper.ci), linetype = 'dashed') +
    #geom_ribbon(aes(ymin = fit[[1]]$lower.ci, ymax = fit[[1]]$upper.ci), fill = "gray", alpha = 0.5) +
    labs(title = title_list[[1]], x = "time t (years)", y = expression(hat(h)(t))) +
    theme_classic(base_size = 16) +
    coord_cartesian(ylim = c(0, max(fit[[1]]$hazard))) +
    theme(plot.title = element_text(size = 18), ,
          panel.background = element_rect(fill = "#eeeeee",
                                          colour = "#eeeeee",
                                          size = 0.5, linetype = "blank"))
  ukplot <- ggplot(mapping = aes(x = fit[[2]]$time / 365.25)) +
    geom_line(aes(y = fit[[2]]$hazard)) +
    labs(title = title_list[[2]], x = "time t (years)", y = expression(hat(h)(t))) +
    theme_classic(base_size = 16) +
    coord_cartesian(ylim = c(0, max(fit[[2]]$hazard))) +
    geom_line(aes(y = fit[[2]]$lower.ci), linetype = 'dashed') +
    geom_line(aes(y = fit[[2]]$upper.ci), linetype = 'dashed') +
    #geom_ribbon(aes(ymin = fit[[2]]$lower.ci, ymax = fit[[2]]$upper.ci), fill = "gray", alpha = 0.5) +
    theme(plot.title = element_text(size = 18), ,
          panel.background = element_rect(fill = "#eeeeee",
                                          colour = "#eeeeee",
                                          size = 0.5, linetype = "blank"))
  plot_grid <- grid.arrange(ausplot, ukplot,
                            nrow = 2,
                            top = text_grob("Estimated force of mortality", size = 20))
  ggsave(filename = "figures/inst-hazard-idl-samples.png", width = 1600, height = 1600, unit = 'px', dpi = 200, bg = "white", plot = plot_grid)
}

plot_Lexis()
plot_survival()
plot_cum_hazard()
plot_hazard()


