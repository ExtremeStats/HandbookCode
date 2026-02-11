


##################
#### Section 1.3.1
##################

library(devtools); install_github("kohrrelation/EVTxgboost")
library(EVTxgboost); library(mev)
library(scoringRules); library(plotrix)
data(GB_DF)

head(colnames(GB_DF)); DF <- GB_DF

indx.train <- which(DF$year<2011)
indx.test <- which(DF$year>=2011)
seqs <- seq(quantile(DF$BA[indx.train],0.8),
            quantile(DF$BA[indx.train],0.99), length.out=35)
par_est <- st_est <- numeric()
for (i in 1:length(seqs)){
  par_est[i] <- gp.fit(DF$BA,seqs[i])$estimate[2]
  st_est[i] <- gp.fit(DF$BA,seqs[i])$std.err[2]
}

pdf(file='ba_dist.pdf', width=5.5, height=2.5)
par(mfrow=c(1,3))
hist(log10(DF$BA[indx.train]+1), breaks=150, xlab='log(BA+1)')
s <- NC.diag(x=DF$BA[indx.train], u=seqs, my.xlab='BA (ac)')
u <- seqs[which(s$e.p.values>0.15)[1]]
abline(v=u, lty='dashed')
plotCI(x=seqs, y=par_est, li=par_est-1.96*st_est,
       ui=par_est+ 1.96*st_est, pch=19, cex=0.5, ylab='xi')
abline(v=u, lty='dashed')
dev.off()
u
quantile(DF$BA[indx.train],0.98)

fit_gpd <- mev::gp.fit(DF$BA[indx.train], u)$estimate
init.sigma <- fit_gpd[1]
y <- DF$BA[indx.train] - u; which.excess <- which(y - u>0)
y <- y[which.excess]; X <- DF[indx.train,3:37][which.excess,]

pdf(file='cv_xi.pdf', width=6, height=4)
model_cv <- GPDxgb.cv(y, X, xi=0.4, init.sigma=init.sigma,
                      stepsize=0.02, nrounds=200, simultaneous=FALSE,
                      xis=seq(0.2,0.6,by=0.05), cv.nfold = 5)
dev.off()
c(model_cv$indx.1se, model_cv$chosen_xi, model_cv$indx.min)

y_test <- DF$BA[indx.test] - u
indx.excess.test <- which(y_test - u>0)
X_test <- DF[indx.test,3:37][indx.excess.test,]
model_fit <- GPDxgb.train(y, X, xi=model_cv$chosen_xi,
                          init.sigma=init.sigma, stepsize=0.02,
                          nrounds=model_cv$indx.1se, simultaneous=FALSE,
                          X_test=X_test)


model_cv_2 <- GPDxgb.cv(y, X, xi=0.4, init.sigma=init.sigma,
                        stepsize=0.01, stepsize_xi=0.01, nrounds=200,
                        simultaneous=TRUE, cv.nfold = 5)
model_fit_2 <- GPDxgb.train(y, X, xi=0.4, init.sigma=init.sigma,
                            stepsize=0.01, stepsize_xi=0.01,
                            nrounds=round(model_cv_2$indx.1se/2),
                            simultaneous = TRUE, X_test=X_test)

test_excess <- y_test[indx.excess.test]
n.exc <- length(test_excess)
crps_basic <- crps_gpd(test_excess, scale=rep(fit_gpd[1], n.exc),
                       shape=rep(fit_gpd[2], n.exc))
crps_model <- crps_gpd(test_excess, scale=model_fit$sigma.pred,
                       shape=rep(model_fit$xi, n.exc))
crps_model_2 <- crps_gpd(test_excess,
                         scale=model_fit_2$sigma.pred,
                         shape=(model_fit_2$xi.pred))
c(mean(crps_basic), mean(crps_model), mean(crps_model_2))






##################
#### Section 1.3.2
##################



library(INLA); library(EVTxgboost); data(France)
hyper.share.gp <- list(prior='gaussian', param=c(0,1/20))

formulae <- y ~ 0 + alpha.bin +
  f(fwi_bin, model=spde1_fwi) +
  f(fa_bin, model=spde1_fa) +
  f(ii_gp_global_2, copy='i_gp_global',
    hyper=list(theta=hyper.share.gp)) +
  f(year_time_indx_bin, model='rw1', diagonal=10^(-4)) +
  alpha.gp + f(fwi, model=spde1_fwi) +
  f(fa, model=spde1_fa) + f(i_gp_global, model=spde) +
  f(year_time_indx_gp, model='rw1', diagonal=10^(-4))

a <- inla.nonconvex.hull(points=loc.data, convex=-0.005)

mesh <- inla.mesh.2d(boundary=a, max.edge=c(15,100))
plot(mesh); nv <- mesh$n

spde <- inla.spde2.pcmatern(mesh=mesh,
                            prior.range=c(20, 0.5),
                            prior.sigma=c(1, 0.5))
A_global <- inla.spde.make.A(mesh=mesh,
                             loc=cbind(dats_gp$loc_x, dats_gp$loc_y))

mesh1d_fwi <- inla.mesh.1d(seq(min(dats_gp$FWI),
                               max(dats_gp$FWI), length.out=4),
                           degree=2, boundary=c('dirichlet','free'))
A_fwi <- inla.spde.make.A(mesh1d_fwi, dats_gp$FWI)
mesh1d_fa <- inla.mesh.1d(seq(-2.1, 1.5, length.out=4), degree=2)
A_fa <- inla.spde.make.A(mesh1d_fa, dats_gp$FA)
spde1_fwi <- inla.spde2.matern(mesh1d_fwi)
spde1_fwi.idx <- inla.spde.make.index("fwi",
                                      n.spde=spde1_fwi$n.spde)
spde1_fwi.idx.b <- inla.spde.make.index("fwi_bin",
                                        n.spde=spde1_fwi$n.spde)
spde1_fa <- inla.spde2.matern(mesh1d_fa, constr=TRUE)
spde1_fa.idx <- inla.spde.make.index("fa",
                                     n.spde=spde1_fa$n.spde)
spde1_fa.idx.b <- inla.spde.make.index("fa_bin",
                                       n.spde=spde1_fa$n.spde)
hyper.gp <- list(theta=list(prior="loggamma", param=c(1,1)))


stk.ygp <- inla.stack(data=list(y=cbind(dats_gp$excess, NA)),
                      A=list(1, A_global, A_fwi, A_fa, 1),
                      effects=list(alpha.gp=rep(1, length(dats_gp$excess)),
                                   i_gp_global=1:nv, spde1_fwi.idx, spde1_fa.idx,
                                   year_time_indx_gp=dats_gp$year_time_indx))
stk.ybin <- inla.stack(data=list(y=cbind(NA,dats_gp$event)),
                       A=list(1, A_global, A_global, A_fwi, A_fa,1),
                       effects=list(alpha.bin=rep(1, length(dats_gp$event)),
                                    i_bin_global=1:nv, ii_gp_global_2=1:nv, spde1_fwi.idx.b,
                                    spde1_fa.idx.b, year_time_indx_bin=dats_gp$year_time_indx))
c.stk <- inla.stack(stk.ygp, stk.ybin)

result <- inla(formulae, data=inla.stack.data(c.stk),
               family=c('gp','binomial'), Ntrials=1,
               control.inla=list(strategy='adaptive', int.strategy='eb'),
               control.family=list(list(control.link=list(quantile=0.5),
                                        hyper=hyper.gp), list()), only.hyperparam=FALSE,
               control.predictor=list(A=inla.stack.A(c.stk),
                                      compute=FALSE))




## PLOTS



# BASIC

xx <- seq((1-mean_FWI)/scale_FWI, (100-mean_FWI)/scale_FWI, length.out=100)
A.xx <- inla.spde.make.A(mesh1d_fwi , xx)
xx2 <- 0
A.xx2 <- inla.spde.make.A(mesh1d_fa, rep(xx2,100) )
mean_FA

stk.ygp.pred <- inla.stack(
  data = list(y=cbind(rep(NA, length(xx)),NA)), #cbind one of each likelihood
  A = list(1, A.xx, A.xx2),
  effects = list(alpha.gp = rep(1, length(xx)),
                 spde1_fwi.idx, spde1_fa.idx), tag='gp.pred' )

stk.ybin.pred <- inla.stack(data = list(y = cbind(rep(NA, length(xx)), NA)),
                            A = list(1, A.xx, A.xx2),
                            effects = list(alpha.bin = rep(1, length(xx)),
                                           spde1_fwi.idx.b, spde1_fa.idx.b), tag='bin.pred' )

pred.stack <- inla.stack(stk.ygp.pred, stk.ybin.pred)

## ------------------------------------------------------------------------
joint.stack <- inla.stack(c.stk , pred.stack)


m.spde <- inla(formulae, data = inla.stack.data(joint.stack) ,
               family =c('gp', 'binomial'),
               control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE),
               control.mode = list(theta=result$mode$theta , restart = FALSE, fixed=TRUE),
               control.family = list(list(control.link=list(quantile=0.5), hyper=hyper.gp),
                                     list()),
               control.inla = list(strategy='adaptive', int.strategy = 'eb'),
               control.compute = list(), verbose=FALSE)


## ------------------------------------------------------------------------

## GPD - FWI relationship
idx <- inla.stack.index(joint.stack, 'gp.pred')$data
tab.spde1 <- data.frame(x = xx*scale_FWI + mean_FWI ,
                        y = exp(m.spde$summary.fitted.values[idx, "mean"]),
                        ll95 = exp(m.spde$summary.fitted.values[idx, "0.025quant"]),
                        ul95 = exp(m.spde$summary.fitted.values[idx, "0.975quant"])
)
pdf(file='GPD_FWI.pdf', width=4.5, height=3)
ggplot() +
  geom_point() +
  # ggtitle("GPD-FWI Relationship") +
  geom_line(aes(x = x, y = y), linetype = 1, data = tab.spde1) +
  geom_ribbon(aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = tab.spde1) +
  xlab('FWI') +ylab('Median BA excesses (ha)')+theme_bw()
#  geom_point(data=data.frame(FWI=dats_gp$FWI[is.na(dats_gp$excess)]*scale_FWI + mean_FWI), aes(y=3, x=FWI), size=0.1, col='blue' ) +
#  geom_point(data=data.frame(FWI=dats_gp$FWI[!is.na(dats_gp$excess)]*scale_FWI + mean_FWI), aes(y=2.95, x=FWI), size=0.1, col='red') +
dev.off()

ilogit <- function(x){
  exp(x)/(1+exp(x))
}


## BIN - FWI relationship
idx <- inla.stack.index(joint.stack, 'bin.pred')$data
tab.spde1 <- data.frame(x = xx*scale_FWI + mean_FWI ,
                        y = ilogit(m.spde$summary.fitted.values[idx, "mean"]),
                        ll95 = ilogit(m.spde$summary.fitted.values[idx, "0.025quant"]),
                        ul95 = ilogit(m.spde$summary.fitted.values[idx, "0.975quant"])
)
pdf(file='BIN_FWI.pdf', width=4.5, height=3)
ggplot() +
  geom_point() +
  # ggtitle("BIN-FWI Relationship") +
  geom_line(aes(x = x, y = y), linetype = 1, data = tab.spde1) +
  geom_ribbon(aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = tab.spde1) +
  xlab('FWI') +ylab('Probability of exceedance')+theme_bw()
dev.off()




###FA



## ------------------------------------------------------------------------
# Predict at a finer grid


## ------------------------------------------------------------------------
# Predict at a finer grid
xx <- seq((0- mean_FA)/ scale_FA , (1- mean_FA)/ scale_FA, length.out=100)
A.xx2 <- inla.spde.make.A(mesh1d_fa, xx)
xx2 <- (mean_FWI-mean_FWI)/scale_FWI
A.xx <- inla.spde.make.A(mesh1d_fwi, rep(xx2,100) )
mean_FWI

stk.ygp.pred <- inla.stack(
  data = list(y=cbind(rep(NA, length(xx)),NA)), #cbind one of each likelihood
  A = list(1, A.xx, A.xx2),
  effects = list(alpha.gp = rep(1, length(xx)),
                 spde1_fwi.idx, spde1_fa.idx), tag='gp.pred' )

stk.ybin.pred <- inla.stack(data = list(y = cbind(rep(NA, length(xx)), NA)),
                            A = list(1, A.xx, A.xx2),
                            effects = list(alpha.bin = rep(1, length(xx)),
                                           spde1_fwi.idx.b, spde1_fa.idx.b), tag='bin.pred' )

pred.stack <- inla.stack(stk.ygp.pred, stk.ybin.pred)

## ------------------------------------------------------------------------
joint.stack <- inla.stack(c.stk , pred.stack)


m.spde <- inla(formulae, data = inla.stack.data(joint.stack) ,
               family =c('gp', 'binomial'),
               control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE),
               control.mode = list(theta=result$mode$theta , restart = FALSE, fixed=TRUE),
               control.family = list(list(control.link=list(quantile=0.5), hyper=hyper.gp),
                                     list()),
               control.inla = list(strategy='adaptive', int.strategy = 'eb'),
               control.compute = list(), verbose=FALSE)


## GPD - FA relationship
idx <- inla.stack.index(joint.stack, 'gp.pred')$data
tab.spde1 <- data.frame(x = xx*scale_FA + mean_FA ,
                        y = exp(m.spde$summary.fitted.values[idx, "mean"]),
                        ll95 = exp(m.spde$summary.fitted.values[idx, "0.025quant"]),
                        ul95 = exp(m.spde$summary.fitted.values[idx, "0.975quant"])
)
pdf(file='GPD_FA.pdf', width=4.5, height=3)
ggplot() +
  geom_point() +
  # ggtitle("GPD-FA Relationship") +
  geom_line(aes(x = x, y = y), linetype = 1, data = tab.spde1) +
  geom_ribbon(aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = tab.spde1) +
  xlab('FA') +ylab('Median BA excesses (ha)') +xlim(0,1)+theme_bw()
dev.off()


ilogit <- function(x){
  exp(x)/(1+exp(x))
}

## BIN - FA relationship
idx <- inla.stack.index(joint.stack, 'bin.pred')$data
tab.spde1 <- data.frame(x = xx*scale_FA + mean_FA ,
                        y = ilogit(m.spde$summary.fitted.values[idx, "mean"]),
                        ll95 = ilogit(m.spde$summary.fitted.values[idx, "0.025quant"]),
                        ul95 = ilogit(m.spde$summary.fitted.values[idx, "0.975quant"])
)
pdf(file='BIN_FA.pdf', width=4.5, height=3)
ggplot() +
  geom_point() +
  # ggtitle("BIN-FA Relationship") +
  geom_line(aes(x = x, y = y), linetype = 1, data = tab.spde1) +
  geom_ribbon(aes(x = x, y = y, ymin = ll95, ymax = ul95), alpha = 0.2, data = tab.spde1) +
  xlab('FA') +ylab('Probability of exceedance')+xlim(0,1)+theme_bw()
dev.off()







