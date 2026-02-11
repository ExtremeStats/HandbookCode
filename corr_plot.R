library("ggplot2")

x.seq <- seq(-1,1,0.01)
y.seq <- seq(-1,1,0.01)
t.seq <- 0:2

full.grid <- expand.grid(x=x.seq, y=y.seq, t=t.seq)

breaks <- 0.999*exp(-1/6*(0:15))

## separable corrrelation
alpha.s <- 0.8
alpha.t <- 1
range.s <- 2
range.t <- 3

cov.sep <- function(x,y,t) exp(- ((x/range.s)^2+(y/range.s)^2)^(alpha.s/2)  - 
                                 abs(t/range.t)^alpha.t )

cov.sep.res <- data.frame(full.grid,
                      t.text= paste("t=", full.grid[,3], sep=""),
                      corr=apply(full.grid, 1, 
                             function(x) cov.sep(x[1], x[2], x[3])))

ggplot(cov.sep.res, aes(x, y, z=corr)) +  geom_raster(aes(fill = corr)) +
  ggtitle("(a) Separable Space-Time Correlation Function") +
  scale_fill_gradientn(colours=rev(hcl.colors(50, palette="viridis"))[1:50], limits=c(0,1)) +
  geom_contour(colour="white", breaks=breaks)  + facet_grid(~t.text)
ggsave(filename="cov_sep.pdf",  plot = last_plot(), width=10, height=3.8)

## correlation gneiting
alpha <- 0.7
beta  <- 0.75
delta <- 0.8
a <- 0.3
c <- 0.5

cov.nsst <- function(x,y,t) 1/(a*abs(t)^(2*alpha)+1)^beta*exp(-c*sqrt(x^2+y^2)^delta/(a*abs(t)^(2*alpha)+1)^beta)

cov.nsst.res <- data.frame(full.grid,
                      t.text= paste("t=", full.grid[,3], sep=""),
                      corr=apply(full.grid, 1, 
                                 function(x) cov.nsst(x[1], x[2], x[3])))

ggplot(cov.nsst.res, aes(x, y, z=corr)) +  geom_raster(aes(fill = corr)) +
  ggtitle("(b) Non-Separable Space-Time Correlation Function from Gneiting's Class") +
  scale_fill_gradientn(colours=rev(hcl.colors(50, palette="viridis"))[1:50], limits=c(0,1)) +
  geom_contour(colour="white", breaks=breaks)  + facet_grid(~t.text)
ggsave(filename="cov_nsst.pdf",  plot = last_plot(), width=10, height=3.8)

## Cox-Isham Model
rho <- 0.5
beta <- 2
mu <- c(0.3,0.3)

cov.ci <- function(x,y,t) {
    
    ci.mat <- matrix(c(1+abs(t)^beta, abs(t)^beta*rho, 
                       abs(t)^beta*rho, 1+abs(t)^beta), nr=2, nc=2)
    shift.vec <- c(x-t*mu[1], y-t*mu[2])
    return(1/sqrt(det(ci.mat))*exp(-sqrt(t(shift.vec%*%solve(ci.mat,shift.vec)))/range.s))
}

cov.ci.res <- data.frame(full.grid,
                      t.text= paste("t=", full.grid[,3], sep=""),
                      corr=apply(full.grid, 1, 
                                 function(x) cov.ci(x[1], x[2], x[3])))

ggplot(cov.ci.res, aes(x, y, z=corr)) +  geom_raster(aes(fill = corr)) +
  ggtitle("(c) Non-Separable Space-Time Correlation Function via Lagrangian Approach") +
  scale_fill_gradientn(colours=rev(hcl.colors(50, palette="viridis"))[1:50], limits=c(0,1)) +
  geom_contour(colour="white", breaks=breaks)  + facet_grid(~t.text)
ggsave(filename="cov_ci.pdf",  plot = last_plot(), width=10, height=3.8)

## tdc extremal-t
nu <- 1
tdc.et <- 2*(1-pt(sqrt((nu+1)*(1-cov.nsst.res$corr)/(1+cov.nsst.res$corr)), df=nu+1))
tdc.et.res <- data.frame(full.grid,
                     t.text= paste("t=", full.grid[,3], sep=""),
                     tdc=tdc.et)
ggplot(tdc.et.res, aes(x, y, z=tdc)) +  geom_raster(aes(fill = tdc)) +
  ggtitle("(a) Extremal-t Process with Non-Separable Space-Time Correlation Function from Gneiting's Class") +
  scale_fill_gradientn(colours=rev(hcl.colors(50, palette="viridis"))[1:50], limits=c(0,1)) +
  geom_contour(colour="white", breaks=breaks)  + facet_grid(~t.text)
ggsave(filename="tdc_et.pdf",  plot = last_plot(), width=10, height=3.8)

## tdc Brown-resnick
vario.sep <- function(x,y,t) ((x/range.s)^2+(y/range.s)^2)^(alpha.s/2) + abs(t/range.t)^alpha.t
tdc.br.res <- data.frame(full.grid,
                         t.text= paste("t=", full.grid[,3], sep=""),
                         tdc=apply(full.grid, 1, 
                                   function(x) 2*(1-pnorm(vario.sep(x[1], x[2], x[3])))))
ggplot(tdc.br.res, aes(x, y, z=tdc)) +  geom_raster(aes(fill = tdc)) +
  ggtitle("(b) Brown-Resnick Process with Separable Variogram of Power Type") +
  scale_fill_gradientn(colours=rev(hcl.colors(50, palette="viridis"))[1:50], limits=c(0,1)) +
  geom_contour(colour="white", breaks=breaks)  + facet_grid(~t.text)
ggsave(filename="tdc_br.pdf",  plot = last_plot(), width=10, height=3.8)