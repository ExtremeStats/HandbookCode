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
library(ggpubr)
library(exdex)

# simulated wind speed data  

plot_windspeed <- function()
  {
  set.seed(22); 
  n <- 100; N <- 365*n; shape <- 2; scale <- 6
  daily <- data.frame(day = as.Date("1899-12-31") + c(1:N),  value = rweibull(N,shape=shape,scale=scale))

u <- 12.5
x <- matrix(daily$value,nrow=n,byrow=T)
y <- apply(x, MARGIN=1, FUN=max)
f <- daily$value %in% y
g <- (daily$value > u)

panel1 <- ggplot(daily, aes(x=day, y=value)) +
  theme_classic(base_size=11) + 
  labs(x="", y="Wind speed (m/s)")+
  geom_line(color="grey") + 
  annotate(geom="segment",x=as.Date("1900-01-01"),xend=as.Date("1999-12-31"),y=12.5) +
  annotate(geom="point",x=daily$day[f],y=daily$value[f]) +
  annotate(geom="segment",x=daily$day[g],y=u,yend=daily$value[g],lwd=0.5) +
  theme(axis.text = element_text(size = 10), legend.position ="none", 
        panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "blank"))

# panel for annual maxima
fit.max <- fgev(daily$value[f])
x <- seq(from=10,to=23,length=1000)
y1 <- 365*dweibull(x,shape,scale)*pweibull(x,shape,scale)^364
y2 <- dgev(x, loc=fit.max$estimate[1],scale=fit.max$estimate[2], shape=fit.max$estimate[3])
panel2 <- ggplot(daily[f,], aes(x=value)) + 
  theme_classic(base_size=11) + 
  ylab("Density") +
  xlab("Annual maximum wind speed (m/s)") +
  geom_histogram(bins=25,fill="lightgrey",aes(y=after_stat(density)))  +
  geom_rug() +
  annotate("line",x=x,y=y1) + 
  annotate("line",x=x,y=y2,col="red",linetype=2) + 
  theme(axis.text = element_text(size = 10), legend.position ="none", 
        panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "blank"))

# panel for exceedances
fit.ex <- fpot(daily$value[g],threshold=u)
x <- seq(from=u,to=23,length=1000)
y1 <- dweibull(x,shape,scale)/(1-pweibull(u,shape,scale))
y2 <- dgpd(x, loc=u,scale=fit.ex$estimate[1], shape=fit.ex$estimate[2])
panel3 <- ggplot(daily[g,], aes(x=value)) + 
  theme_classic(base_size=11) + 
  ylab("Density") +
  xlab("Wind speed (m/s)") +
  geom_histogram(bins=25,fill="lightgrey",aes(y=after_stat(density)))  +
  geom_rug() +
  annotate("line",x=x,y=y1) + 
  annotate("line",x=x[-1],y=y2[-1],col="red",linetype=2)+ 
  theme(axis.text = element_text(size = 10), legend.position ="none", 
        panel.background = element_rect(fill = "white",colour = "white",linewidth = 0.5, linetype = "blank"))

all <- ggarrange(panel1,ggarrange(panel2,panel3,ncol=2),nrow=2)

ggsave(filename = "figures/weibull.png", plot = all,
         bg = "white", width = 2000, height = 2000, unit = 'px', dpi=250)
}

# some GEV densities

plot_gev <- function (){
n <- 10001
y <- seq(from=-6,to=10,length=n)
d1 <- dgev(y,shape=0)
d2 <- dgev(y,shape=-0.5)
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
  geom_point(x=2,y=0,pch=16)+
  labs(x = expression(y), y = "Density") +
    ylim(0, max(d1,d2,d3))+
    theme_classic(base_size=11) +
  theme(axis.text = element_text(size = 10), legend.position ="none", 
        panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
ggsave(filename = "figures/gev-pdfs.png", plot = plot, bg = "white", width = 2000, height = 1000, unit = 'px', dpi=250)
}

# some GPD densities

plot_gpd <- function (){
  x <- seq(from=10^(-10),to=10,length=1000)
  plot(x,dgpd(x,shape=0),type="l",xlab="x",ylab="PDF", ylim=c(0,1),xlim=c(0,10))
  lines(x,dgpd(x,shape=1/2),col="red")
  lines(x[1:500],dgpd(x[1:500],shape=-1/5),col="blue")
  points(5,0,pch=16,col="blue")
  
  x <- seq(from=10^(-10),to=1-10^(-10),length=1000)
  plot(x,dgpd(x,shape=-1),type="l",xlab="x",ylab="PDF", ylim=c(0,5),xlim=c(0,1))
  lines(x,dgpd(x,shape=-2, scale=2),col="red")
  lines(x,dgpd(x,shape=-1/2,scale=1/2),col="blue")
  lines(x,dgpd(x,shape=-1/1.25,scale=1/1.25),col="green")
  
n <- 500
y1 <- seq(from=10^(-10),to=10,length=n)
y2 <- seq(from=10^(-10),to=1-10^(-10),length=n)
shapes1 <- c(0, 0.5, -0.2)
scales1 <- c(1,1,1)
shapes2 <- c(-1,-2,-1/2, -1/1.25)
scales2 <- c(1, 2, 1/2, 1/1.25)
d1 <- dgpd(y1,shape=shapes1[1])
d2 <- dgpd(y1,shape=shapes1[2])
d3 <- dgpd(y1,shape=shapes1[3])
d4 <- dgpd(y2,shape=shapes2[1],scale=scales2[1])
d5 <- dgpd(y2,shape=shapes2[2],scale=scales2[2])
d6 <- dgpd(y2,shape=shapes2[3],scale=scales2[3])
d7 <- dgpd(y2,shape=shapes2[4],scale=scales2[4])
D1 <- list(d1, d2, d3)
D2 <- list(d4, d5, d6, d7)
df <- list()
for (i in 1:3){
  di <- D1[[i]]
  xi <- shapes1[[i]]
  sigma <- scales1[[i]]
  df[[i]] <- data.frame(y = y1[di>0], pdf = di[di>0], scale=sigma, shape=i)
}
df_t1 <- dplyr::bind_rows(df)
df <- list()
for (i in 1:4){
  di <- D2[[i]]
  xi <- shapes2[[i]]
  sigma <- scales2[[i]]
  df[[i]] <- data.frame(y = y2[di>0], pdf = di[di>0], scale=sigma, shape=i)
}
df_t2 <- dplyr::bind_rows(df)

# Create the first plot
p1 <-  ggplot(df_t1, aes(x = y, y = pdf, linetype = as.factor(shape))) +
  geom_line() +
  geom_point(x=5,y=0,pch=16) + 
  labs(x = expression(y),  y = "Density") +
  ylim(0, 1) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(size = 10),
        legend.position = "none",
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5, linetype = "blank"))

p2 <-  ggplot(df_t2, aes(x = y, y = pdf, linetype = as.factor(shape))) +
  geom_line() +
  labs(x = expression(y),  y = "Density") +
  ylim(0, 5) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(size = 10),
        legend.position = "none",
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5, linetype = "blank"))

pl <- cowplot::plot_grid(plotlist = list(p1, p2),
                         labels = c("", ""),
                         ncol = 2)

ggsave(filename = "figures/gpd-pdfs.png", plot = pl, bg = "white", width = 2000, height = 1000, unit = 'px', dpi=250)
}

# plots for Abisko precipitation analysis
# these data don't have zeros, so exceedance probabilities will be incorrect, so 
# we add zeros on days with no rain

data(abisko) 
junk <- seq(from=min(abisko$date),to=max(abisko$date),by="day")
junk.y <- rep(0,length(junk))
junk.y[junk %in% abisko$date] <- abisko$precip
abisko <- data.frame(date=as.Date(junk),precip=junk.y)

# abisko, first plot of data

plot_abisko1 <- function () {
dates <- as.Date(c("1940-01-01","1979-12-31"))
plot <- ggplot(abisko, mapping=aes(x=date, y=precip))+
  annotate("rect", xmin = as.Date("1940-01-01"), xmax = as.Date("1979-12-31"), ymin = 10, ymax = 79, linetype=0, fill="grey")+
  geom_point(pch=16,cex=0.7)+
  labs(x="", y = "Precipitation (mm)")+
  scale_y_continuous(limits = c(0,79), expand = c(0, 0))+
  scale_x_date(limits = as.Date(c('1913-01-01','2015-01-01')), expand = c(0, 0))+
    theme_classic(base_size=11)+ 
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
ggsave(filename = "figures/abisko1.png", plot = plot, bg = "white", width = 2000, height = 1000, unit = 'px', dpi=250)
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
    geom_point(aes(x=qgumbel(c(1:mon.n)/(mon.n+1)), y=sort(mon.max)), pch=16,cex=0.7)+
    geom_point(aes(x=qgumbel(c(1:year.n)/(year.n+1)), y=sort(year.max)), pch=16, cex=0.7)+
    geom_point(aes(x=qgumbel(c(1:twoyear.n)/(twoyear.n+1)), y=sort(twoyear.max)), pch=15, cex=0.7)+
    geom_point(aes(x=qgumbel(c(1:fiveyear.n)/(fiveyear.n+1)), y=sort(fiveyear.max)), pch=17, cex=0.7)+
    geom_point(aes(x=qgumbel(c(1:tenyear.n)/(tenyear.n+1)), y=sort(tenyear.max)), pch=18, cex=1)+
    labs(y="Ordered maxima (mm)", x="Gumbel plotting positions")+
    theme_classic(base_size=11)+
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5, linetype = "blank"))
  
  pll <- gev.pll(psi=seq(10,150,length=141),param='quant',dat=year.max, p=1/200, PLOT=FALSE)
  dev.off()
  x<- pll$psi
  y<-pll$pll-pll$maxpll
  subs <- (y>=-6.5)
  c <- confint(pll)
  d <- confint(pll, level=0.99)
  inf <- c[2]
  inf_ <- d[2]
  plot2 <- ggplot()+
    geom_vline(xintercept=pll$psi.max, color="grey")+
    geom_hline(yintercept=y[x>=inf][1], color="slategrey")+
    geom_hline(yintercept=y[x>=inf_][1], color="slategrey")+
    geom_line(aes(x=x[subs], y=y[subs]))+
    xlim(40,140)+
    labs(x=expression(psi), y="Profile log-likelihood")+
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  pl <- cowplot::plot_grid(plotlist = list(plot1, plot2),
                                 labels = c("", ""),
                                 ncol = 2)
  ggsave(filename = "figures/abisko2.png", plot = pl, bg = "white", width = 2000, height = 900, unit = 'px', dpi=250)
}

plot_abisko3 <- function (){
 # Northrop-Coleman plot 
  abisko.NC <- NC.diag(xdat=c(abisko$precip),u=2*c(0:20),plot=FALSE) #),cex.lab=2) #,mar=c(3,5,2,2))
  df <- data.frame(u=abisko.NC$u,p=abisko.NC$e.p.values)
  panel1 <- ggplot(df,aes(u,p)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(0,1) + 
    xlim(0,40) +
    labs(y="P-value", x=expression(u)) +
    annotate("segment",x=0,xend=40,y=0.05,color="grey")+
    geom_point()+ geom_line()
  
  # mean residual life plot
  junk <-  mrlplot(data=c(abisko$precip),tlim=c(0,40),nt=21)
  df <- data.frame(u=junk$x,junk$y)
  panel2 <- ggplot(df,aes(u,upper)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(0,14) + 
    xlim(0,40) +
    labs(y="Mean excess", x=expression(u)) +
    annotate("line",x=df$u,y=df$mrl) +  
    annotate("line",x=df$u,y=df$lower,linetype=2) +  
    annotate("line",x=df$u,y=df$upper,linetype=2)
  
  # threshold choice plots
  junk <- tcplot(data=c(abisko$precip), tlim=c(0,40), nt=21, model="pp",which=1)
  df <- data.frame(u=2*c(0:20),junk$locs)
  panel3 <- ggplot(df,aes(u,upper)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(range(df[,-1])) + 
    xlim(0,40) +
    labs(y="Location", x=expression(u)) +
    geom_linerange(x=df$u,ymin=df[,2],ymax=df[,4],col="grey") +  
    annotate("line",x=df$u,y=df[,3]) + annotate("point",x=df$u,y=df[,3]) 
   
  junk <- tcplot(data=c(abisko$precip), tlim=c(0,40), nt=21, model="pp",which=2)
  df <- data.frame(u=2*c(0:20),junk$scales)
  panel4 <- ggplot(df,aes(u,upper)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(range(df[,-1])) + 
    xlim(0,40) +
    labs(y="Scale", x=expression(u)) +
    geom_linerange(x=df$u,ymin=df[,2],ymax=df[,4],col="grey") +  
    annotate("line",x=df$u,y=df[,3]) + annotate("point",x=df$u,y=df[,3]) 
  
  junk <- tcplot(data=c(abisko$precip), tlim=c(0,40), nt=21, model="pp",which=3)
  df <- data.frame(u=2*c(0:20),junk$shapes)
  panel5 <- ggplot(df,aes(u,upper)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(range(df[,-1])) + 
    xlim(0,40) +
    labs(y="Shape", x=expression(u)) +
    geom_linerange(x=df$u,ymin=df[,2],ymax=df[,4],col="grey") +  
    annotate("line",x=df$u,y=df[,3]) + annotate("point",x=df$u,y=df[,3]) 
  
  # final layout
  right <- ggarrange(panel2,panel1,ncol=1,nrow=2) 
  left <- ggarrange(panel3,panel4,panel5,ncol=1,nrow=3)
  all <- ggarrange(left, right, ncol=2,nrow=1)
  
  ggsave(filename = "figures/abisko3.png", plot = all, bg = "white", 
         width = 2000, height = 1800, unit = 'px', dpi=250)
}

plot_abisko4 <- function (){
  # extremal index
  junk <- seq(from=min(abisko$date),to=max(abisko$date),by="day")
  junk.y <- rep(0,length(junk))
  junk.y[junk %in% abisko$date] <- abisko$precip
  quant <- mean((junk.y>12))
  n.u <- sum((junk.y>12))
  q <- quant + 2*sqrt(quant*(1-quant)/n.u)
  u <- th.u <- seq(from=0, to=20, length=101)
  lims.u <- cbind(u,u)
  for (i in 1:length(u)) 
  { 
    th <- exdex::kgaps(junk.y, u=u[i], k=1)
    th.u[i] <- th$theta
    lims.u[i,] <- th$theta + th$se*1.96*c(-1,1)
  }
  df <- data.frame(u=u,th=th.u,llim=lims.u[,1],ulim=lims.u[,2])
  
  # extremal index plot
  panel1 <- ggplot(df,aes(u,th)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(0.4,1) +
    xlim(0,20) +
    labs(y=expression(widehat(theta)[u]), x=expression(u)) +
    annotate("line", x=df$u, y=df$llim, col="grey") +
    annotate("line", x=df$u, y=df$ulim, col="grey") +    
    geom_line() 

  # extremogram
  junk <- extremogram::extremogram1(junk.y,quant=1-quant,maxlag=21,plot=FALSE,type=1)
  junk <- junk[-1]
  n <- length(junk)
  df <- data.frame(h=c(1:n),th=junk)
  panel2 <- ggplot(df,aes(h,th)) +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5, linetype = "blank")) +
    ylim(0,0.2) +
    labs(y=expression(widehat(pi)[u](h)), x=expression(h)) +
    annotate("segment",x=1,xend=20, y=q, linetype=2,col="grey") +
    annotate("segment",x=1,xend=20, y=quant, col="grey") +
    geom_linerange(x=df$h,ymin=rep(0,n),ymax=df$th)
  
  # final layout
  all <- ggarrange(panel1, panel2, ncol=2,nrow=1)
  
  ggsave(filename = "figures/abisko4.png", plot = all, bg = "white", 
         width = 2000, height = 900, unit = 'px', dpi=250)
  
}

# make plots

plot_windspeed()
plot_gpd()
plot_gev()
plot_abisko1()
plot_abisko2()
plot_abisko3()
plot_abisko4()

# analysis of abisko data

fpot(abisko$precip,threshold=12,"pp", # point process fit 
             start=list(loc=20,scale=6,shape=0.08),npp=365.25)

fpot(abisko$precip,threshold=12,"gpd") # GP fit

# make annual maxima

abisko.max <- matrix(NA, 102, 12)
year <- c(1913:2014)
for (i in 1:102) for (j in 1:12)
{ k <- (year(abisko$date)-1912==i & month(abisko$date)==j)
abisko.max[i,j] <- max(abisko$precip[k]) }
year.max <- apply(abisko.max,1,max)

fgev(year.max)  # GEV fit to annual maxima
