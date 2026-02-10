library(evd)
library(extremogram)
library(evd)
library(extremogram)
library(evd)
library(mev)
library(lubridate)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(data.table)

load("Data/SP500.RData")

sp <- ts(sp500[,2],start=c(2002,28), frequency=252)
ret <- -100*diff(log(sp))

#### Home-made AR-GARCH filtering
garch <- function(x, m=30, y)
{
  n <- length(y)
  j <- (m+1):n
  sig <- rep(NA,n)
  sig[m] <- var(y[1:m])
  mu <- x[1]
  rho <- x[2]
  b0 <- x[3]
  b1 <- x[4]
  d <- x[5]
  for (i in (m+1):n) sig[i] <- b0+b1*(y[i-1]-mu)^2 + d*sig[i-1]
  -sum( dnorm(y[j], mu+rho*(y[j-1]-mu), sqrt(sig[j]), log=T) )
}

# first fit a model with normal errors to get initial values for t fit
init <- c(mean(ret),0.0,0.2,0.2,0.1)
fit.g30 <- nlm( garch, init, hessian=T, print.level=2, y=ret)

# t error fit
garch <- function(x, m=30, y)
{
  n <- length(y)
  j <- (m+1):n
  sig <- rep(NA,n)
  sig[m] <- var(y[1:m])
  mu <- x[1]
  rho <- x[2]
  b0 <- x[3]
  b1 <- x[4]
  d <- x[5]
  df <- x[6]
  for (i in (m+1):n) sig[i] <- (b0+b1*(y[i-1]-mu)^2 + d*sig[i-1])
  sig <- sig*(df-2)/df
  z <- (y[j]- mu-rho*(y[j-1]-mu))/sqrt(sig[j])
  -sum( dt(z, df, log=T) - log(sqrt(sig[j])) )
}


init <- c(fit.g30$estimate,20)
fit.g30.t <- nlm( garch, init, hessian=T, print.level=2, y=ret)
fit.g30.t <- nlm( garch, fit.g30.t$estimate, hessian=T, print.level=2, y=ret)

fit.g30.t$estimate
sqrt(diag(solve(fit.g30.t$hessian)))

# make residuals from the t fit
garch.res <- function(fit, m=30, y)
{
  n <- length(y)
  j <- (m+1):n
  sig <- rep(NA,n)
  sig[m] <- var(y[1:m])
  x <- fit$estimate
  mu <- x[1]
  rho <- x[2]
  b0 <- x[3]
  b1 <- x[4]
  d <- x[5]
  for (i in (m+1):n) sig[i] <- b0+b1*(y[i-1]-mu)^2 + d*sig[i-1]
  (y[j]-( mu+rho*(y[j-1]-mu)))/sqrt(sig[j])
}

ret.gres <- ts(garch.res( fit=fit.g30.t, y=ret),start=c(2002,28),frequency=252)
#par(mfrow=c(2,2))
#acf(ret.gres); pacf(ret.gres); acf(ret.gres^2); pacf(ret.gres^2)

plot_sp1 <- function (){
  df <- data.frame(time = time(ret), value = as.numeric(ret))
  plot1 <- ggplot(data=df, aes(x=time, y=value))+
    geom_vline(xintercept = c(2002:2025),col="grey",lwd=0.5)+
    geom_line()+
    labs(y="negative log returns(%)")+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18), panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "blank"))

  df_filtered <- data.frame(time = time(ret.gres), value = as.numeric(ret.gres))
  plot2 <- ggplot(data=df_filtered, aes(x=time, y=value))+
    geom_vline(xintercept = c(2002:2025),col="grey",lwd=0.5)+
    geom_line()+
    labs(y="filtered negative log returns(%)")+
    ylim(c(min(ret), max(ret)))+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18), panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "blank"))
  q <- 0.05
  maxlag<-250
  ex1<-extremogram::extremogram1(ret,1-q,type=1, maxlag=maxlag)
  dev.off()
  plot3 <- ggplot()+
    geom_segment(aes(x=0 : (maxlag - 1), xend=0 : (maxlag - 1), y=0, yend=ex1), lwd=0.3)+
    geom_hline(yintercept = q,col="slategrey",lwd=0.5)+
    geom_hline(yintercept=q + 2*sqrt(1/length(ret)),col="slategrey",linetype="dotdash")+
    labs(y="extremogram", x="lag")+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18), panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "blank"))

  ex2<-extremogram::extremogram1(ret.gres,1-q,type=1, maxlag=250)
  dev.off()
  plot4 <- ggplot()+
    geom_segment(aes(x=0 : (maxlag - 1), xend=0 : (maxlag - 1), y=0, yend=ex2), lwd=0.3)+
    geom_hline(yintercept = q,col="slategrey",lwd=0.5)+
    geom_hline(yintercept=q + 2*sqrt(1/length(ret.gres)),col="slategrey", linetype="dotdash")+
    labs(y="extremogram", x="lag")+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18), panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "blank"))
  total <- cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4),
                                   labels = c("", "", "", ""),
                                   ncol = 2, nrow=2)
  ggsave("figures/sp1.png", bg = "white", dpi = 200, unit = 'px', width = 2400, height = 1700,plot=total)
}

plot_sp1()