library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(graphics)
library(data.table)
library(POT)
library(extremogram)
library(readr)
library(mgcv)
library(tidyverse)
library(lubridate)
library(ggdist)
library(boot)
library(ismev)
library(boot)
library(evd)
load("Data/Venice.RData")
venice_hourly_over80_1982_2023$time <- as.POSIXct(venice_hourly_over80_1982_2023$time)
venice <- ten_yearly_max_venice_1887_2023 %>%
  tidyr::pivot_longer(cols = -year, names_to = "level", values_to = "y")
y <-  venice$y[venice$level=="y1"&venice$year<2020]
x <- venice$year[venice$level=="y1"&venice$year<2020]-1900
fit.gam <- mgcv::gam(list(y~x,~1,~1),
                     family = gevlss(link = list("identity", "identity", "identity")))

plot_tenmax <- function() {
  fig <- ggplot() +
    geom_point(data = tidyr::gather(ten_yearly_max_venice_1887_2023, max_n, level_cm, -year), aes(x = year, y = level_cm), color = "black", pch = 16, cex = 0.7) +
    labs(x = "", y = "sea level(cm)") +
    geom_line(aes(x+1900,fitted(fit.gam)[,1]),col="slategrey")+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18), panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  ggsave("figures/venice1.png", bg = "white", dpi = 200, unit = 'px', width = 1600, height = 900,plot=fig)
}

# useful functions first

qqgumbel <- function (y, env=T, gum.env = NA, yl=NULL, ...)
{ # qqplot of y against Gumbel order statistics, with simulation envelope if requested
  i <- !is.na(y) & is.finite(y) & !is.nan(y)
  y <- y[i]
  n <- length(y)
  x <- qgumbel(c(1:n)/(n + 1))
  if (is.null(yl)) yl <- range(y)
  if (env)
  {
    library(boot)
    if (is.na(gum.env)) {
      gum.sim <- boot(y, function(d) sort(rgumbel(length(d))), R=10000, sim="parametric", mle=NA)
      gum.env <- envelope( gum.sim ) }
    if (is.null(yl)) yl <- range(gum.env$overall,y)
  }
  qq <- qqplot(x, y)
  dev.off()
  plot <- ggplot(data=data.frame(qq), aes(x=x))+
        geom_point(aes(y=y), pch=16,cex=0.7)+
    geom_abline(aes(slope=1, intercept=0), color="slategrey")+
    labs(x="Gumbel plotting positions", y="ordered spacings")+
    ylim(yl)+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  if (env)
  {
    plot <- plot +
      geom_line(aes(x, gum.env$point[1, ]), lty = 4)+
    geom_line(aes(x, gum.env$point[2, ]), lty = 4)+
    geom_line(aes(x, gum.env$overall[1, ]), lty = 1)+
    geom_line(aes(x, gum.env$overall[2, ]), lty = 1)
  }
  return(plot)
}

make.exp.res <- function(fit, data, r, env=TRUE, yl=NULL)
{ # make residuals using r-largest fit and differences of exponentials on transformed scale
  j <- r+1
  y <- data[2:j]
  eta <- fit$vals[,1]
  tau <- fit$vals[,2]
  xi <- fit$vals[,3]
  res <- (y-eta)/tau
  pp.res <- res.exp <- (1 + xi*res)^(-1/xi)

  for (i in 2:r) res.exp[i] <- pp.res[i] - pp.res[i-1]
  plot_list <- list()
  for (i in 1:r)
  {
    res <- -log(res.exp[i])
    plot_list[[i]] <- qqgumbel(res[[paste0("y", i)]], yl=yl, env=env)
  }
  return(plot_list)
}

plot_gam_fit <- function (){
  fit1 <- fgev(x=y,nsloc=x-1900) # straight-line regression
  muhat <- fit1$estimate[1] + fit1$estimate[2]*(x-1900)
  sigmahat <- fit1$estimate[3]
  xihat <- fit1$estimate[4]
  res <- (y-muhat)/sigmahat
  res <- log( 1 + xihat*res )/xihat
  n <- length(res)
  gum.pp <- qgumbel(c(1:n)/(n+1))

  # fit to top three order statistics each year
  Xdat<-ten_yearly_max_venice_1887_2023[ten_yearly_max_venice_1887_2023$year<2020, 2:4]
  X <- cbind(ten_yearly_max_venice_1887_2023$year-1900,(ten_yearly_max_venice_1887_2023$year>1981),(ten_yearly_max_venice_1887_2023$year-1900)*(ten_yearly_max_venice_1887_2023$year>1981))
  Ydat <- X[ten_yearly_max_venice_1887_2023$year<2020,]
  fitr3 <- rlarg.fit(xdat=Xdat, ydat=Ydat, method="BFGS", mul=c(1))
  gum.sim <- boot(res,function(d) sort(rgumbel(length(d))), R=10000, sim="parametric")
  gum.env <- envelope( gum.sim )
  yl <- c(-2.5,12)
  qq <- qqplot(gum.pp, res)
  dev.off()
  plot <- ggplot(data=data.frame(qq), aes(x=x))+
        geom_point(aes(y=y), pch=16,cex=0.71)+
    geom_abline(aes(slope=1, intercept=0), color="slategrey")+
    geom_line(aes(gum.pp, gum.env$point[1, ]), lty = 4)+
    geom_line(aes(gum.pp, gum.env$point[2, ]), lty = 4)+
    geom_line(aes(gum.pp, gum.env$overall[1, ]), lty = 1)+
    geom_line(aes(gum.pp, gum.env$overall[2, ]), lty = 1)+
    labs(x="Gumbel plotting positions", y="ordered residuals")+
    ylim(yl)+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  junk <- make.exp.res(fit=fitr3, data=ten_yearly_max_venice_1887_2023[ten_yearly_max_venice_1887_2023$year<2020, ], r=3, yl=yl)
  plot_list <- c(list(plot), junk)
  total <- cowplot::plot_grid(plotlist = plot_list,
                                 labels = c("", ""),
                                 ncol = 2, nrow=2)
  ggsave("figures/venice2.png", bg = "white", dpi = 200, unit = 'px', width = 2400, height = 1400,plot=total)

}
plot_tenmax()
plot_gam_fit()