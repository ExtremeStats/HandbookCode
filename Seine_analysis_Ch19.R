wd <- "path/to/working/directory/Seine_analysis_Ch19/"
setwd(wd)
set.seed(22)

library(devtools)
devtools::install_github("nicolagnecco/causalXtreme")
library(causalXtreme)
library(sna)
install_version("SID",version="1.0")
library(lubridate)

### load data
load("Data/data_seine_Ch19.RData")
# daily average water levels for each station

################################################################################
############################# EASE method ######################################
################################################################################
reach_EASE <- causal_discovery(data_seine[,-c(1,5)],set_arg= list("both_tails" = FALSE))$est_g

################################################################################
############################# RMLM-based method ################################
################################################################################
source("Code/extreme_DAG_Krali_scaling.R")
order_data2 <- causalorder(frechet2(t(data_seine[,-c(1,5)])), floor(sqrt(nrow(data_seine))) , .0001, .0001, a=1.01)
MLmat_data2 <- t(apply(Bmatrix(pmax(frechet2(t(data_seine[,-c(1,5)])[order_data2,])-
                                      apply(frechet2(t(data_seine[,-c(1,5)])[order_data2,]),1,mean),0),
                               floor(sqrt(nrow(data_seine)))),1,function(x)x/sum(x^2)^.5))
reach_RMLM  <- t(apply(MLmat_data2,c(1,2), function(x) as.numeric(x>0)))
diag(reach_RMLM) <- rep(0,4)
reach_RMLM[order(order_data2, decreasing = FALSE),order(order_data2, decreasing = FALSE)]

################################################################################
############################# CausEV method ####################################
################################################################################
source("Code/QCDD.R")

B        <- 300
toremove <- which(year(as.Date(data_seine[,1],"%d/%m/%Y",origin="01/01/1970")) %in% 2008)
dd       <- data_seine[-toremove,]
dd_order <- order(as.Date(dd[,1],"%d/%m/%Y",origin="01/01/1970"))
tmp      <- year(as.Date(dd[,1],"%d/%m/%Y",origin="01/01/1970"))

QCDD21   <- QCDD31 <- QCDD41 <- QCDD32 <- QCDD42 <- QCDD43 <- NULL
# DO NOT RUN (takes about 6 hours). Could be parallelized (takes about 1h30 on 20 cores)
for (i in 1:B)
{
  yy <- as.numeric(sample(names(table(tmp)),replace=T))
  new.dd <- NULL
  for (k in 1:length(yy))
  {
    indefix <- which(year(as.Date(dd[,1],"%d/%m/%Y",origin="01/01/1970")) == yy[k])
    new.dd  <- rbind(new.dd, dd[dd_order[indefix],])
  }
  
  thd <- 0.9
  n.ext.min   <- length(which((new.dd[,2+1]>quantile(new.dd[,2+1],thd+0.005))&
                                (new.dd[,1+1]>quantile(new.dd[,1+1],thd+0.005))))
  QCDD21[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,2+1],new.dd[,1+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,3+1]>quantile(new.dd[,3+1],thd+0.005))&
                                (new.dd[,1+1]>quantile(new.dd[,1+1],thd+0.005))))
  QCDD31[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,3+1],new.dd[,1+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,4+1]>quantile(new.dd[,4+1],thd+0.005))&
                                (new.dd[,1+1]>quantile(new.dd[,1+1],thd+0.005))))
  QCDD41[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,4+1],new.dd[,1+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,3+1]>quantile(new.dd[,3+1],thd+0.005))&
                                (new.dd[,2+1]>quantile(new.dd[,2+1],thd+0.005))))
  QCDD32[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,3+1],new.dd[,2+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  n.ext.min   <- length(which((new.dd[,4+1]>quantile(new.dd[,4+1],thd+0.005))&
                                (new.dd[,2+1]>quantile(new.dd[,2+1],thd+0.005))))
  
  QCDD42[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,4+1],new.dd[,2+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,4+1]>quantile(new.dd[,4+1],thd+0.005))&
                                (new.dd[,3+1]>quantile(new.dd[,3+1],thd+0.005))))
  QCDD43[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,4+1],new.dd[,3+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
}
QCDD <- cbind(QCDD21,QCDD31,QCDD41,QCDD32,QCDD42,QCDD43)

# construct the reachability matrix for CausEV
mat <- matrix(0,ncol=4,nrow=4)

if((quantile(QCDD[,1], 0.025)>0.5)&(quantile(QCDD[,1],0.975)>0.5))
  mat[2,1] <- 1
if((quantile(QCDD[,1], 0.025)<0.5)&(quantile(QCDD[,1],0.975)<0.5))
  mat[1,2] <- 1

if((quantile(QCDD[,2], 0.025)>0.5)&(quantile(QCDD[,2],0.975)>0.5))
  mat[3,1] <- 1
if((quantile(QCDD[,2], 0.025)<0.5)&(quantile(QCDD[,2],0.975)<0.5))
  mat[1,3] <- 1

if((quantile(QCDD[,3], 0.025)>0.5)&(quantile(QCDD[,3],0.975)>0.5))
  mat[4,1] <- 1
if((quantile(QCDD[,3], 0.025)<0.5)&(quantile(QCDD[,3],0.975)<0.5))
  mat[1,4] <- 1

if((quantile(QCDD[,4], 0.025)>0.5)&(quantile(QCDD[,4],0.975)>0.5))
  mat[3,2] <- 1
if((quantile(QCDD[,4], 0.025)<0.5)&(quantile(QCDD[,4],0.975)<0.5))
  mat[2,3] <- 1

if((quantile(QCDD[,5], 0.025)>0.5)&(quantile(QCDD[,5],0.975)>0.5))
  mat[4,2] <- 1
if((quantile(QCDD[,5], 0.025)<0.5)&(quantile(QCDD[,5],0.975)<0.5))
  mat[2,4] <- 1

if((quantile(QCDD[,6], 0.025)>0.5)&(quantile(QCDD[,6],0.975)>0.5))
  mat[4,3] <- 1
if((quantile(QCDD[,6], 0.025)<0.5)&(quantile(QCDD[,6],0.975)<0.5))
  mat[3,4] <- 1

reach_CausEv       <- reachability(mat)
diag(reach_CausEv) <- rep(0,4)


################################################################################
########################### Bootstrap for SID ##################################
################################################################################

######## bootstrap for CausalEx #########

aux.yy   <- year(as.Date(data_seine[,1],"%d/%m/%Y",origin="01/01/1970"))
table(aux.yy) # shows that there is only one value in 2008. We remove it:
toremove <- which(year(as.Date(data_seine[,1],"%d/%m/%Y",origin="01/01/1970")) %in% 2008)
dd       <- data_seine[-toremove,]
dd_order <- order(as.Date(dd[,1],"%d/%m/%Y",origin="01/01/1970"))
aux.yy   <- year(as.Date(dd[,1],"%d/%m/%Y",origin="01/01/1970"))

# ###################
# true reachability matrix
dag <- matrix(0,ncol=4,nrow=4)
dag[2,] <- c(1,0,0,0)
dag[3,] <- c(1,0,0,0)
dag[4,] <- c(1,0,1,0)

B <- 500
QCDD21 <- numeric(B)
QCDD31 <- numeric(B)
QCDD41 <- numeric(B)
QCDD32 <- numeric(B)
QCDD42 <- numeric(B)
QCDD43 <- numeric(B)

SIDgnecco <- SID_Mario <- numeric(B)

# Krali's method
source("Code/extreme_DAG_Krali_scaling.R")
#source functions from causalXtreme package
source("https://raw.githubusercontent.com/nicolagnecco/causalXtreme/refs/heads/master/R/graph_theory_functions.R")

# DO NOT RUN (takes about 10 hours). Could be parallelized
for (i in 1:B){
  yy     <- as.numeric(sample(names(table(aux.yy)),replace=T))
  new.dd <- NULL
  for (k in 1:length(yy))
  {
    indefix <- which(year(as.Date(dd[,1],"%d/%m/%Y",origin="01/01/1970")) == yy[k])
    new.dd  <- rbind(new.dd, dd[dd_order[indefix],])
  }
  
  thd <- 0.9
  
  SIDgnecco[i] <- compute_str_int_distance(dag,causal_discovery(new.dd[,-c(1,5)],set_arg= list("both_tails" = FALSE))$est_g)
  
  order_data2    <- causalorder(frechet2(t(new.dd[,-c(1,5)])), floor(sqrt(nrow(new.dd))) , .0001, .0001, a=1.01)
  MLmat_data2    <- t(apply(Bmatrix(pmax(frechet2(t(new.dd[,-c(1,5)])[order_data2,])-apply(frechet2(t(new.dd[,-c(1,5)])[order_data2,]),1,mean),0),floor(sqrt(nrow(new.dd)))),1,function(x)x/sum(x^2)^.5))
  dag_RMLM       <- t(apply(MLmat_data2,c(1,2), function(x) as.numeric(x>0)))
  diag(dag_RMLM) <- rep(0,4)
  SID_Mario[i]   <- compute_str_int_distance(dag[order_data2,order_data2], dag_RMLM)
  
  ### QCDD
  n.ext.min   <- length(which((new.dd[,2+1]>quantile(new.dd[,2+1],thd+0.005))&
                                (new.dd[,1+1]>quantile(new.dd[,1+1],thd+0.005))))
  QCDD21[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,2+1],new.dd[,1+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,3+1]>quantile(new.dd[,3+1],thd+0.005))&
                                (new.dd[,1+1]>quantile(new.dd[,1+1],thd+0.005))))
  QCDD31[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,3+1],new.dd[,1+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,4+1]>quantile(new.dd[,4+1],thd+0.005))&
                                (new.dd[,1+1]>quantile(new.dd[,1+1],thd+0.005))))
  QCDD41[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,4+1],new.dd[,1+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,3+1]>quantile(new.dd[,3+1],thd+0.005))&
                                (new.dd[,2+1]>quantile(new.dd[,2+1],thd+0.005))))
  QCDD32[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,3+1],new.dd[,2+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  n.ext.min   <- length(which((new.dd[,4+1]>quantile(new.dd[,4+1],thd+0.005))&
                                (new.dd[,2+1]>quantile(new.dd[,2+1],thd+0.005))))
  
  QCDD42[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,4+1],new.dd[,2+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
  
  n.ext.min   <- length(which((new.dd[,4+1]>quantile(new.dd[,4+1],thd+0.005))&
                                (new.dd[,3+1]>quantile(new.dd[,3+1],thd+0.005))))
  QCDD43[i] <- QCDD_extremes_upper_quad(data.frame(new.dd[,4+1],new.dd[,3+1]),thd =thd,n.sim.ext = n.ext.min)$epsilon
}
QCDD <- cbind(QCDD21,QCDD31,QCDD41,QCDD32,QCDD42,QCDD43)

### compute_str_int_distance returns a normalized SID which is useful when comparing performances on dags of different sizes. Here, we have a unique dag with four nodes so we work with the non-normalized SID
SIDgnecco.l <- quantile(SIDgnecco,0.025)*(nrow(dag)*(nrow(dag)-1))#0
SIDgnecco.u <- quantile(SIDgnecco,0.975)*(nrow(dag)*(nrow(dag)-1))#5.525

SID_Mario.l <- quantile(SID_Mario,0.025)*(nrow(dag_RMLM)*(nrow(dag_RMLM)-1))#0
SID_Mario.u <- quantile(SID_Mario,0.975)*(nrow(dag_RMLM)*(nrow(dag_RMLM)-1))#6

SIDcausEv <- numeric()
dag <- matrix(0,ncol=4,nrow=4)
dag[2,] <- c(1,0,0,0)
dag[3,] <- c(1,0,0,0)
dag[4,] <- c(1,0,1,0)

for (i in 1:nrow(QCDD))
{
  mat <- matrix(0,ncol=4,nrow=4)
  if(QCDD[i,1] > 0.5 ) {mat[2,1] <- 1} else mat[1,2] <- 1
  if(QCDD[i,2] > 0.5 ) {mat[3,1] <- 1} else mat[1,3] <- 1
  if(QCDD[i,3] > 0.5 ) {mat[4,1] <- 1} else mat[1,4] <- 1
  if(QCDD[i,4] > 0.5 ) {mat[3,2] <- 1} else mat[2,3] <- 1
  if(QCDD[i,5] > 0.5 ) {mat[4,2] <- 1} else mat[2,4] <- 1
  if(QCDD[i,6] > 0.5 ) {mat[4,3] <- 1} else mat[3,4] <- 1
  mat_reach <- reachability(mat)
  diag(mat_reach) <- rep(0,4)
  SIDcausEv[i] <- compute_str_int_distance(dag,mat_reach)
}

SIDcausev.l <- quantile(SIDcausEv,0.025)*(nrow(dag)*(nrow(dag)-1))#0
SIDcausev.u <- quantile(SIDcausEv,0.975)*(nrow(dag)*(nrow(dag)-1))#4


