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

load("../Data/SP500.RData")

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

# first fit a model with normal errors to get initial values for t fit
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

plot_sp1 <- function (){
  df <- data.frame(time = time(ret), value = as.numeric(ret))
  plot1 <- ggplot(data=df, aes(x=time, y=value))+
    geom_vline(xintercept = c(2002:2024),col="grey",linewidth=0.3)+
    geom_line()+
    labs(y="Daily losses (%)",x="")+
    scale_x_continuous(breaks=seq(2000,2024,4))+
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10), panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "blank"))

  df_filtered <- data.frame(time = time(ret.gres), value = as.numeric(ret.gres))
  plot2 <- ggplot(data=df_filtered, aes(x=time, y=value))+
    geom_vline(xintercept = c(2002:2024),col="grey",linewidth=0.3)+
    scale_x_continuous(breaks=seq(2000,2024,4))+
    geom_line()+
    labs(y="Filtered daily losses (%)",x="")+
    ylim(c(min(ret), max(ret)))+
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10), panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "blank"))
  q <- 0.05
  qtop <- q + 2*sqrt(1/length(ret))
  maxlag<-251
  ex1<-extremogram::extremogram1(ret,1-q,type=1, maxlag=maxlag)
  junk <- ex1[-1]
  n <- length(junk)
  df <- data.frame(h=c(1:n),th=junk)
  plot3 <- ggplot(df,aes(h,th)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(0,0.2) +
    labs(y=expression(widehat(pi)[u](h)), x=expression(h)) +
    geom_linerange(x=df$h,ymin=rep(0,n),ymax=df$th,color="grey") +
    annotate("segment",x=min(df$h),xend=max(df$h), y=q) +
    annotate("segment",x=min(df$h),xend=max(df$h), y=qtop, linetype=2) 
  
  ex2<-extremogram::extremogram1(ret.gres,1-q,type=1, maxlag=250)
  junk <- ex2[-1]
  n <- length(junk)
  df <- data.frame(h=c(1:n),th=junk)
  plot4 <- ggplot(df,aes(h,th)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(0,0.2) +
    labs(y=expression(widehat(pi)[u](h)), x=expression(h)) +
    geom_linerange(x=df$h,ymin=rep(0,n),ymax=df$th,color="grey") +
    annotate("segment",x=min(df$h),xend=max(df$h), y=q) +
    annotate("segment",x=min(df$h),xend=max(df$h), y=qtop, linetype=2) 
  
  library(ggpubr)
  top <- ggarrange(plot1, plot2, ncol=2,nrow=1)
  bottom <- ggarrange(plot3, plot4, ncol=2,nrow=1)
  all <- ggarrange(top,bottom,ncol=1,nrow=2,heights=c(1,0.7))
  ggsave(filename = "figures/sp1.png", plot = all, bg = "white", 
         width = 2000, height = 1600, unit = 'px', dpi=250)
}

plot_sp1()

library(evd)
(fit1 <- fit.pp(ret,threshold=quantile(ret,0.95),npp=252))
(fit2 <- fit.pp(ret.gres,threshold=quantile(ret.gres,0.95),npp=252))