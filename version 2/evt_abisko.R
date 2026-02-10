library(evd)
library(mev)
library(scales)
library(lubridate)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(ggdist)
library(data.table)
library(extremogram)

plot_windspeed <- function (){
  Lambda <- function (q, loc = 0, scale = 1, shape = 0)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  q <- (q - loc)/scale
  if (shape == 0)
    p <- exp(-q)
  else p <- pmax(1 + shape * q, 0)^(-1/shape)
  p
}

lambda <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0)
    d <- log(1/scale) - x
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -Inf
  }
  if (!log)
    d <- exp(d)
  d
}

nlogL0 <- function(th, z, y)
{ m <- ncol(y)+1
ts1 <- -m*Lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3]) +
  lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3], log=TRUE) +log(m)
ts2 <- -Lambda(y, loc=th[1], scale=exp(th[2]), shape=th[3]) +
  lambda(y, loc=th[1], scale=exp(th[2]), shape=th[3], log=TRUE)
ts3 <- (m-1)*Lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3])
-2*(sum(ts1) + sum(ts2) + sum(ts3))
}

nlogL <- function(th, z, y)
{ m <- ncol(y)+1
  ts1 <- -m*Lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3]) +
            lambda(z, loc=th[1], scale=exp(th[2]), shape=th[3], log=TRUE) + log(m)
  ts2 <- - Lambda(y, loc=th[4], scale=exp(th[5]), shape=th[6]) +
           lambda(y, loc=th[4], scale=exp(th[5]), shape=th[6], log=TRUE)
  ts3 <- (m-1)*Lambda(z, loc=th[4], scale=exp(th[5]), shape=th[6])
  -2*(sum(ts1) + sum(ts2) + sum(ts3))
}
max.test <- function(y)
{
  # each row of y contains maxima
  fit2 <- fgev(c(y))  # original fit for comparison
  m <- ncol(y)
  s <- y
  for (i in 1:nrow(y)) s[i,] <- sort(y[i,])
  z <- s[,m]
  y <- s[,-m]
  init <- fit2$estimate
  init[2] <- log(init[2])
  fit0 <- optim(init, nlogL0, method="BFGS", z=z, y=y)
  init <- jitter(c(init, init))
  fit1 <- optim(init, nlogL, method="BFGS", z=z, y=y)
  list(L2 = fit2$deviance, L0=fit0$value, L1=fit1$value, W=fit0$value-fit1$value,
       th2=fit2$estimate, th0=fit0$par, th1=fit1$par,
       conv1=fit1$convergence, conv0=fit0$convergence)
}
# graphs to illustrate basic points
n <- 100; N <- 365*n; shape <- 2; scale <- 6
x <- matrix(rweibull(N,shape=shape,scale=scale),nrow=n)
y <- apply(x, MARGIN=1, FUN=max)
fit <- fgev(y)
u.p <- 0.9
u <- quantile(c(x),u.p)
fit1 <- fpot(c(x),threshold=u)
amax <- 23
h <- hist(x,prob=T,nclass=200, plot=FALSE)
a <- seq(from=0, to=amax, length=500)
da <- dweibull(a, shape=shape, scale=scale)
  par(mfrow=c(1,2),pty="s")
# density panel
  plot1 <- ggplot()+
    geom_step(aes(h$mids,h$density))+
    geom_line(aes(a[-1],da[-1]),color="slategrey",lwd=0.7)+
    geom_vline(xintercept = u, color="slategrey", linetype=2)+
    geom_line(aes(a,dgev(a,loc=fit$estimate[1],scale=fit$estimate[2],shape=fit$estimate[3])/365),lty=3,lwd=0.7,color="blue")+
    geom_line(aes(a[a>u],fit1$pat*dgpd((a-u)[a>u],scale=fit1$estimate[1],shape=fit1$estimate[2])),lty=2,lwd=0.7,col="red")+
    labs(x="windspeed (m/s)", y="density")+
    xlim(0,amax)+
    scale_y_log10(limits=c(10^-4,0.5), expand=c(0,0))+
    theme_classic(base_size=18)
  # quantile panel
  q <- qweibull(c(1:N)/(N+1),shape=shape,scale=scale)
  z <- c(10,20,50,100,200,500,1000,5000)
  w <- seq(from=u, to=amax, length=500)
  qw.gpd <- u+qgpd((pweibull(w,shape=shape,scale=scale)-(u.p-0.01))/(1-(u.p-0.01)), scale=fit1$estimate[1], shape=fit1$estimate[2])
  Np <- 5000
  p <- c(1:Np)/(Np+1)
  w.gev <- qweibull(p^(1/365), scale=scale, shape=shape)
  qw.gev <- qgev(p, loc=fit$estimate[1], scale=fit$estimate[2], shape=fit$estimate[3])
  plot2 <- ggplot()+
    geom_point(aes(q,sort(c(x))), pch=".", cex=2)+
    geom_line(aes(w, qw.gpd), lty=2,lwd=0.7, col="red")+
    geom_line(aes(w.gev, qw.gev),lty=3,lwd=0.7, col="blue")+
    geom_vline(xintercept = u, col="slategrey",lty=2)+
    geom_vline(xintercept=qweibull(1-1/(365*z),shape=shape,scale=scale),color="slategrey",lwd=0.5)+
    geom_abline(aes(slope=1,intercept=0), col="slategrey")+
    geom_text(aes(x=qweibull(1-1/(365*z),shape=shape,scale=scale),y=1*c(0:7),label=z, cex=0.4), show.legend = FALSE)+
    xlim(0, amax)+
    ylim(0, amax)+
    labs(y="empirical", x="theoretical")+
    theme_classic(base_size=18)
  pl <- cowplot::plot_grid(plotlist = list(plot1, plot2),
                           #rel_widths = c(1, 1.2),
                                 labels = c("", ""),
                                 ncol = 2)
  ggsave(filename = "figures/weibull.png", plot = pl,
         bg = "white", width = 3200, height = 1080, unit = 'px', dpi=250)
}

plot_gev <- function (){
n <- 1001
y <- seq(from=-6,to=10,length=n)
d1 <- dgev(y,shape=-0.5)
d2 <- dgev(y,shape=0)
d3 <- dgev(y,shape=0.5)
d <- list(d1, d2, d3)
shapes <- c(-0.5, 0, 0.5)
df <- list()
for (i in 1:3){
  di <- d[[i]]
  xi <- shapes[[i]]
  df[[i]] <- data.frame(y = y[di>0], pdf = di[di>0], shape=xi)
}
df_t <- dplyr::bind_rows(df)
plot <- ggplot(df_t, mapping=aes(x=y, y=pdf, linetype=as.factor(shape)))+
  geom_line()+
  geom_point(x=-2,y=0,pch=16)+
  geom_point(x=2,y=0,pch=15)+
  labs(x = expression(y),
         y = "density function") +
    ylim(0, max(d1,d2,d3))+
    theme_classic(base_size=18) +
    scale_linetype_manual(name = "\U03Be", guide = "legend", values=c("solid", "dotted", "dashed"), labels =  sort(shapes))+
  theme(axis.text = element_text(size = 18),
        panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
ggsave(filename = "figures/gev-pdfs.png", plot = plot, bg = "white", width = 2500, height = 1080, unit = 'px')
}


plot_gpd <- function (){
n <- 500
y <- seq(from=10^(-6),to=10,length=n)
d1 <- dgpd(y,shape=0)
d2 <- dgpd(y,shape=-0.3)
d3 <- dgpd(y,shape=0.3)
d4 <- dgpd(y,shape=-0.7,scale=2)
d5 <- dgpd(y,shape=-1,scale=1/0.6)
d6 <- dgpd(y,shape=-2,scale=6)
d <- list(d1, d2, d3, d4, d5, d6)
shapes <- c(0, 0.3, -0.3, -0.7, -1, -2)
scales <- c(1, 1, 1, 2, 1/0.6, 6)
df <- list()
for (i in 2:6){
  di <- d[[i]]
  xi <- shapes[[i]]
  sigma <- scales[[i]]
  df[[i]] <- data.frame(y = y[di>0], pdf = di[di>0], scale=sigma, shape=xi)
}
df_t <- dplyr::bind_rows(df)
  df_t$scale <- round(df_t$scale, 2)
df_t$shape <- round(df_t$shape, 2)
  # Create a new variable for the combination of scale and shape
df_t$scale_shape <- interaction(df_t$scale, df_t$shape, sep='/')

# Define the unique linetype and color values for each combination
unique_combinations <- unique(df_t$scale_shape)
linetypes <- c("solid", "dotted", "solid", "dotdash", "dotted")
colors <- c("grey", "grey", "black", "black", "black")

# Ensure we have enough linetypes and colors
all_linetypes <- rep(linetypes, length.out = length(unique_combinations))
all_colors <- rep(colors, length.out = length(unique_combinations))

# Create a named vector for linetypes and colors
linetype_mapping <- setNames(all_linetypes, unique_combinations)
color_mapping <- setNames(all_colors, unique_combinations)

# Create the plot
plot <-  ggplot(df_t, aes(x = y, y = pdf, linetype = scale_shape, color = scale_shape)) +
  geom_line() +
  labs(x = expression(y),
       y = "density function",
       linetype = "\u03C3/\U03BE",  # Combined legend title
       color = "\u03C3/\U03BE") +
  ylim(0, 1) +
  theme_classic(base_size = 18) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = linetype_mapping) +
  theme(axis.text = element_text(size = 18),
        legend.position = c(0.8, 0.8),  # Position legend inside the graph
        legend.background = element_rect(fill = alpha('white', 0.6)),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "blank"))
ggsave(filename = "figures/gpd-pdfs.png", plot = plot, bg = "white", width = 2400, height = 1400, unit = 'px')
}

# abisko
plot_abisko1 <- function () {
dates <- as.Date(c("1940-01-01","1979-12-31"))
plot <- ggplot(abisko, mapping=aes(x=date, y=precip))+
  geom_vline(xintercept=dates, linetype=2)+
  geom_segment(aes(x = as.Date("1940-01-01"), xend = as.Date("1979-12-31"), y = 10, yend = 10), linetype=2)+
  geom_point(pch=16,cex=0.7)+
  ylim(0, NA)+
  labs(x="", y = "daily precipitation (mm)")+
  scale_y_continuous(limits = c(0,79), expand = c(0, 0))+
  scale_x_date(limits = as.Date(c('1913-01-01','2015-01-01')), , expand = c(0, 0))+
    theme_classic(base_size=18)+ theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
ggsave(filename = "figures/abisko1.png", plot = plot, bg = "white", width = 2200, height = 1080, unit = 'px')
}

plot_abisko2 <- function (){
  abisko.max <- matrix(NA, 102, 12)
  year <- c(1913:2014)
  for (i in 1:102) for (j in 1:12)
    { k <- (year(abisko$date)-1912==i & month(abisko$date)==j)
    abisko.max[i,j] <- max(abisko$precip[k]) }
  mon.max <- c(abisko.max)
  mon.n <- length(mon.max)
  year.max <- apply(abisko.max,1,max)
  year.n <- length(year.max)
  twoyear.max <- apply(matrix(year.max[-c(1,2)], ncol=2, byrow=T),1,max)
  twoyear.n <- length(twoyear.max)
  fiveyear.max <- apply(matrix(year.max[-c(1,2)], ncol=5, byrow=T),1,max)
  fiveyear.n <- length(fiveyear.max)
  tenyear.max <- apply(matrix(year.max[-c(1,2)], ncol=10, byrow=T),1,max)
  tenyear.n <- length(tenyear.max)
  plot1 <- ggplot()+
    geom_point(aes(x=qgumbel(c(1:mon.n)/(mon.n+1)), y=sort(mon.max)), pch="+",cex=0.5)+
    geom_point(aes(x=qgumbel(c(1:year.n)/(year.n+1)), y=sort(year.max)), pch=16, cex=0.5)+
    geom_point(aes(x=qgumbel(c(1:twoyear.n)/(twoyear.n+1)), y=sort(twoyear.max)), pch=15, cex=0.5)+
    geom_point(aes(x=qgumbel(c(1:fiveyear.n)/(fiveyear.n+1)), y=sort(fiveyear.max)), pch=17, cex=0.5)+
    geom_point(aes(x=qgumbel(c(1:tenyear.n)/(tenyear.n+1)), y=sort(tenyear.max)), pch=18, cex=0.5)+
    labs(y="ordered maxima(mm)", x="Gumbel plotting positions")+
    theme_classic(base_size=18)
  pll <- gev.pll(psi=seq(10,150,length=141),param='quant',dat=year.max, p=1/200, PLOT=FALSE)
  dev.off()
  x<- pll$psi
  y<-pll$pll-pll$maxpll
  c <- confint(pll)
  d <- confint(pll, level=0.99)
  inf <- c[2]
  inf_ <- d[2]
  plot2 <- ggplot()+
    geom_line(aes(x=x[x>42], y=y[x>42]))+
    geom_vline(xintercept=pll$psi.max)+
    geom_hline(yintercept=y[x>=inf][1], color="slategrey")+
    geom_hline(yintercept=y[x>=inf_][1], color="slategrey")+
    xlim(10,150)+
    scale_y_continuous(expand = c(0,0))+
    labs(x=expression(psi), y="profile log-likelihood")+
    theme_classic(base_size=18) +
    theme(axis.text = element_text(size = 18),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  pl <- cowplot::plot_grid(plotlist = list(plot1, plot2),
                                 labels = c("", ""),
                                 ncol = 2)
  ggsave(filename = "figures/abisko2.png", plot = pl, bg = "white", width = 2500, height = 1080, unit = 'px')
}

plot_abisko3 <- function (){
  png("figures/abisko3.png", width = 1800, height = 1100, unit = 'px')
  par(mfrow=c(2,2), cex=3, mar = c(4, 4, 2, 0.5), las=1)
  NC.diag(xdat=c(abisko$precip),u=2*c(0:20),
               xlim=c(0,40))
  abline(h=c(0.01,0.05),col="grey")
  tcplot(data=c(abisko$precip), tlim=c(1,40), model="pp",
              ylab=c(expression(mu),expression(sigma), expression(xi)), xlab="threshold")
  dev.off()
}

plot_abisko4 <- function (){
  junk <- seq(from=min(abisko$date),to=max(abisko$date),by="day")
  junk.y <- rep(0,length(junk))
  junk.y[junk %in% abisko$date] <- abisko$precip
  quant <- mean((junk.y>10))
  n.u <- sum((junk.y>10))
  q <- quant + 2*sqrt(quant*(1-quant)/n.u)
  png(file="figures/abisko4.png", width = 1600, height = 500, unit = 'px')
  par(mfrow=c(1,2),mar = c(4, 4.5, 0.5, 0.5), cex=3)
  exiplot(junk.y,c(0,30),ylab=expression(widehat(theta)[u]),ylim=c(0.4,1), xlab = "threshold", lwd=2)
  extremogram1(junk.y,quant=1-quant,maxlag=20,type=1)
  abline(h=q,col="black", lty=2, lwd=1.5)
  abline(h=quant, col="black", lwd=1.5)
  dev.off()
}

# plot_windspeed()
# plot_gpd()
# plot_gev()
# plot_abisko1()
# plot_abisko2()
# plot_abisko3()
plot_abisko4()