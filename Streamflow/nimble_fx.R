# This code has a few functions need for NIMBLE

# https://github.com/paciorek/lbl-extremes-2019/blob/master/3-nimble-eva.Rmd

library(methods)
library(nimble)

dgev <- nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0), xi = double(0), 
        log = integer(0, default = 0)) {
        
        returnType(double(0))
        std <- (x - mu) / sigma
        if(1 + xi*std <= 0) {
             logProb <- -Inf
        } else logProb <- -log(sigma) - (1+1/xi) * log(1 + xi * std) - (1 + xi * std)^(-1/xi)
        if(log) return(logProb)
        else return(exp(logProb))
    })

rgev <- nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0), xi = double(0)) {
        returnType(double(0))
        if(n != 1) print("rgev only allows n = 1; using n = 1.")
        u <- runif(1)
        if (xi == 0) {
           x = -log(-log(u))
        } else {
           x = ((-log(u))^(-xi) - 1.0)/xi
        }
        return (x*sigma + mu)
    })

assign('dgev', dgev, .GlobalEnv)
assign('rgev', rgev, .GlobalEnv)

# Method of moments for initial values
initialize <- function(Y){
  scale <- sqrt(6)*sd(Y,na.rm=TRUE)/pi
  loc   <- median(Y,na.rm=TRUE)+scale*log(log(2))
  shape <- 0.001
return(c(loc,log(scale),shape))}
