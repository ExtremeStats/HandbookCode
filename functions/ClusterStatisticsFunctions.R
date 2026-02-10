#' Estimation of p(v, d, r)
#'
#' @param frwrd.sim forward simulations from the 'forwardSim' function
#' @param u threshold above which to consider clusters
#' @param L number of lags after an exceedance to be considered in fitting procedure
#' @param plot.pmf boolean; true to plot the pmf, false otherwise
#'
#' @return pmf of cluster sizes
#' @export
#'
#' @examples
cluster.size.pmf <- function(frwrd.sim, u, L){
  # Function to compute p(u, L, j)
  pmf.clust.length.scalar <- function(frwrd.sim,u,L,j){
    res <-apply(frwrd.sim[,1:L],1,function(x){
      as.numeric(sum(x>u)==j)}
    )
    mean(res)
  }

  sapply(c(1:L),function(j) pmf.clust.length.scalar(frwrd.sim,u,L,j))
}


#' Estimation of the maximum path following a threshold exceedance
#'
#' @param frwrd.sim forward simulations from the 'forwardSim' function
#' @param plot.max boolean; true to plot, false otherwise
#'
#' @import graphics
#' @return path-wise maximum of forward simulations
#' @export
#'
#' @examples
cluster.max <- function(frwrd.sim,plot.max=FALSE){
  cluster.m <- apply(frwrd.sim,2,max)
  if(plot.max==TRUE){
    par(mfrow=c(1,1),mar=c(4.1,5.1,0.5,1.1),mgp=c(2.6,0.8,0))
    plot(frwrd.sim)
    matplot(t(frwrd.sim),xlab="L",ylab=expression(X[L]),
            type="l",col="black",lty="solid")
    lines(cluster.m,col="red")
  }
  cluster.m
}

#' Estimation of the mean path following a threshold exceedance
#'
#' @param frwrd.sim
#' @param plot.mean
#'
#' @import graphics
#' @return
#' @export
#'
#' @examples
cluster.mean <- function(frwrd.sim,plot.mean=FALSE){
  cluster.m <- apply(frwrd.sim,2,mean)
  if(plot.mean==TRUE){
    par(mfrow=c(1,1),mar=c(4.1,5.1,0.5,1.1),mgp=c(2.6,0.8,0))
    plot(frwrd.sim)
    matplot(t(frwrd.sim),xlab="L",ylab=expression(X[L]),
            type="l",col="black",lty="solid")
    lines(cluster.m,col="red")
  }
  cluster.m
}


#' Simultaneous prediction intervals for a path following an exceedance
#'
#' @param frwrd.sim forward simulations from the 'forwardSim' function
#'
#' @import excursions
#' @import graphics
#' @return
#' @export
#'
#' @examples
cluster.simPredict <- function(frwrd.sim,alpha=0.05,plot.sim=FALSE){
  require(excursions)
  excurs <- simconf.mc(samples = t(frwrd.sim),alpha = alpha)
  sim.conf <- cbind(lower=excurs$a,upper=excurs$b)
  if(plot.sim==TRUE){
    matplot(t(frwrd.sim),xlab="L",ylab=expression(X[L]),
            type="l",col="black",lty="solid")
    lines(sim.conf[,1],col="red")
    lines(sim.conf[,2],col="red")
  }
  sim.conf
}

#' Estimation of the extremal index using the intervals estimator of Ferro/Seger
#'
#' @param X numeric vector of observations
#' @param u numeric, threshold for identifying exdeedances
#'
#' @return estimate of the extremal index
#' @export
#'
#' @examples
ExtrIndex_Intervals <- function(X, u){
  exc <- which(X > u)  # exceedance times
  T.exc <- diff(exc)  # interexceedance times
  if(max(T.exc) <= 2){
    n <- length(T)
    Estimate <- (2*(sum(T.exc))^2)/(n*sum(T.exc^2))
    theta <- min(Estimate,1)
    return(theta)
  }else{
    n <- length(T.exc)
    Num <- 2*(sum(T.exc-1))^2
    Den <- n*sum((T.exc-1)*(T.exc-2))
    Estimate <- Num/Den
    theta <- min(Estimate,1)
    return(theta)
  }
}
