# For plotting:
library(ggplot2)
library(scales)

# For times series models
library(rugarch)
library(tseries)

# For the extremal index 
library(exdex)

### Load and plot the daily negative log-returns
load("sp500.RData")
n <- nrow(data)
head(data)

ggplot(data, aes(x = Date, y = SP500)) + 
  geom_hline(yintercept = 0, color = "grey") +
  geom_line(aes(y =  SP500), lwd = 0.5) +
  coord_cartesian(ylim=c(-21.5,21.5)) +  
  theme_bw() + xlab("Time (years)") + ylab("Negative log-returns (%)") + 
  theme(axis.title.x = element_text(size=19), 
        axis.title.y = element_text(size=19), 
        axis.text.x = element_text(size=16),  
        axis.text.y = element_text(size=16)) 

### Funtion to make the Pareto quantile plots
parquantplot <- function(data, k){
  topdata <- sort(data, decreasing = T)[1:k]
  x <- -log((1:k)/(k+1))
  y <- log(topdata)
  datapar <- data.frame('x' = x, 'y' = y)
  return(datapar)
}

datapar1 <- parquantplot(data$SP500, k = 50)
datapar2 <- parquantplot(data$SP500, k = 500)
regr1 <- lm(datapar1$y ~ datapar1$x)
regr2 <- lm(datapar2$y ~ datapar2$x)
### Initial estimates for the tail index alpha:
1/coef(regr1)[2]
1/coef(regr2)[2]


ggplot(datapar1, aes(x = x, y = y)) + 
  geom_point(aes(x = x, y =  y)) +
  geom_line(aes(x = x, y = fitted(regr1)), col = "blue") + 
  xlab(expression(-log(i/(k+1)))) + 
  ylab(expression(log(X[(n-i+1)]))) + 
  theme_bw() + ggtitle("k = 50") + 
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 

ggplot(datapar2, aes(x = x, y = y)) + 
  geom_point(aes(x = x, y =  y)) +
  geom_line(aes(x = x, y = fitted(regr2)), col = "blue") + 
  xlab(expression(-log(i/(k+1)))) + 
  ylab(expression(log(X[(n-i+1)]))) + 
  theme_bw() + ggtitle("k = 500") + 
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 
dev.off()

### Autocorrelation plots for the raw data

bacf <- acf(data$SP500, plot = FALSE, lag.max = 20)
bacfdf <- with(bacf, data.frame(lag, acf))[-1,]
ciline <- qnorm((1 - 0.95)/2)/sqrt(n)
ggplot(data = bacfdf, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + ylab("ACF") + 
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'blue', lwd = 0.5) + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'blue', lwd = 0.5) +
  theme_bw() + ggtitle("Daily negative log-returns of the S&P500") + 
  theme(plot.title = element_text(size=19, hjust=0.5),
        axis.title.x =element_text(size=17), 
        axis.title.y =element_text(size=17),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 

bacf <- acf(data$SP500^2, plot = FALSE, lag.max = 20)
bacfdf <- with(bacf, data.frame(lag, acf))[-1,]
ggplot(data = bacfdf, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + ylab("ACF") + 
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'blue', lwd = 0.5) + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'blue', lwd = 0.5) +
  theme_bw() + ggtitle("Squared daily negative log-returns of the S&P500") + 
  theme(plot.title = element_text(size=19, hjust=0.5),
        axis.title.x =element_text(size=17), 
        axis.title.y =element_text(size=17),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 


### Fit an AR(1)-GARCH(1,1) model and extract the residuals
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)))
gfit <- ugarchfit(model, data$SP500) 
eps <- gfit@fit$z

### Autocorrelation plots for the residuals

bacf <- acf(eps, plot = FALSE, lag.max = 20)
bacfdf <- with(bacf, data.frame(lag, acf))[-1,]
ggplot(data = bacfdf, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + ylab("ACF") + 
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'blue', lwd = 0.5) + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'blue', lwd = 0.5) +
  theme_bw() + ggtitle("GARCH(1,1) residuals") + 
  theme(plot.title = element_text(size=19, hjust=0.5),
        axis.title.x =element_text(size=17), 
        axis.title.y = element_text(size=17), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 

bacf <- acf(eps^2, plot = FALSE, lag.max = 20)
bacfdf <- with(bacf, data.frame(lag, acf))[-1,]
ggplot(data = bacfdf, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + ylab("ACF") + 
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'blue', lwd = 0.5) + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'blue', lwd = 0.5) +
  theme_bw() + ggtitle("Squared GARCH(1,1) residuals") + 
  theme(plot.title = element_text(size=19, hjust=0.5),
        axis.title.x =element_text(size=17), 
        axis.title.y = element_text(size=17), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 


### Estimate the extremal index (takes one minute)
bval <- seq(20, 800, by = 10) # block sizes
res <- choose_b(data$SP500, bval, interval_type = "lik") 
theta <- data.frame('b' = bval, 'thest' = res$theta_sl[,2])

ggplot(theta, aes(x = b, y = thest)) + 
  geom_line(aes(y =  thest), lwd = 0.9) +
  geom_ribbon(aes(ymin= res$lower_sl[,2], ymax = res$upper_sl[,2], color = "Black"), show.legend = FALSE, alpha = 0.125, linetype = 'dotted') +
  coord_cartesian(ylim=c(0,0.7)) + theme_bw() + 
  xlab("Block size") + ylab("Extremal index") + 
  scale_color_manual(name = "Legend", values =  c('Black' = 'black')) +
  theme(axis.title.x =element_text(size=18), 
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))

# Give point estimate and confidence interval for the extremal index
# We use the "BB2018 sliding blocks" method
theta <- spm(data$SP500, b = 500)
confint(theta, interval_type = "lik")

### Declustering: take every fifth observation
dataMo <- data[weekdays(data$Date) == "Monday",]
dataTu <- data[weekdays(data$Date) == "Tuesday",]
dataWe <- data[weekdays(data$Date) == "Wednesday",]
dataTh <- data[weekdays(data$Date) == "Thursday",]
dataFr <- data[weekdays(data$Date) == "Friday",]

thetaMo <- spm(dataMo$SP500, b = 200)
cithetaMo <- confint(thetaMo, interval_type = "lik")
thetaTu <- spm(dataTu$SP500, b = 200)
cithetaTu <- confint(thetaTu, interval_type = "lik")
thetaWe <- spm(dataWe$SP500, b = 200)
cithetaWe <- confint(thetaWe, interval_type = "lik")
thetaTh <- spm(dataTh$SP500, b = 200)
cithetaTh <- confint(thetaTh, interval_type = "lik")
thetaFr <- spm(dataFr$SP500, b = 200)
cithetaFr <- confint(thetaFr, interval_type = "lik")

# Gives the point estimates and confidence intervals:
rbind(c(thetaMo$theta_sl[2], cithetaMo$ciBB$prof_CI),
      c(thetaTu$theta_sl[2], cithetaTu$ciBB$prof_CI),
      c(thetaWe$theta_sl[2], cithetaWe$ciBB$prof_CI),
      c(thetaTh$theta_sl[2], cithetaTh$ciBB$prof_CI),
      c(thetaFr$theta_sl[2], cithetaFr$ciBB$prof_CI))

### Declustering: second procedure
declFun <- function(data, range = 2){ 
  n <- nrow(data)
  dataord <- data[order(data$SP500, decreasing = T),]
  indzero <- which(dataord$SP500 <= 0)[1]
  indneg <- which(dataord$SP500 < 0)[1]
  temp1 <- as.numeric(rownames(dataord[1:(indzero-1),]))
  temp2 <- as.numeric(rownames(dataord[indzero:(indneg-1),]))
  temp3 <- rev(as.numeric(rownames(dataord[indneg:n,])))

  i <- 2
  while(i <= length(temp1)){
    check <- any(sapply(c(1:(i-1)), function(j) temp1[i] %in% (temp1[j]-range):(temp1[j]+range)))
    if(check){
      temp1 <- temp1[-i]
    } else{
      i <- (i+1)
    }
  }
  i <- 2
  while(i <= length(temp3)){
    check <- any(sapply(c(1:(i-1)), function(j) temp3[i] %in% (temp3[j]-range):(temp3[j]+range)))
    if(check){
      temp3 <- temp3[-i]
    } else{
      i <- (i+1)
    }
  }
  datadecl <- data[sort(c(temp1,temp2,temp3)),] 
  return(datadecl)
}

# Takes about a minute
datadecl <- declFun(data)
thetadecl <- spm(datadecl$SP500, b = 250)
confint(thetadecl, interval_type = "lik")

# Takes less than a minute
datadecl <- declFun(data, range = 9)
thetadecl <- spm(datadecl$SP500, b = 150)
confint(thetadecl, interval_type = "lik")

# Takes less than a minute
datadecl <- declFun(data, range = 20)
thetadecl <- spm(datadecl$SP500, b = 100)
confint(thetadecl, interval_type = "lik")

### Simulation study
fixed.p <- as.list(coef(gfit))
varModel <- list(model = "sGARCH", garchOrder = c(1,1))
spec <- ugarchspec(varModel, mean.model = list(armaOrder = c(1,0)), fixed.pars = fixed.p) 

#takes about two hours; results are stored in SimStudy.RDS
#res <- matrix(0, nrow = 100, ncol = 4)
#for(i in 1:100){ 
#  x <- ugarchpath(spec, n.sim = n, m.sim = 1, rseed = i) 
#  X <- fitted(x)[,1] 
#  theta <- spm(X, b = 500)
#  res[i,1] <- theta$theta_sl[2]
  
#  dataMo <- X[weekdays(data$Date) == "Monday"]
#  thetaMo <- spm(dataMo, b = 200)
#  res[i,2] <- thetaMo$theta_sl[2]
  
#  dataX <- data
#  dataX$SP500 <- X
#  datadecl <- declFun(dataX)
#  thetadecl <- spm(datadecl$SP500, b = 250)
#  res[i,3] <- thetadecl$theta_sl[2]
  
#  datadecl2 <- declFun(dataX, range = 9)
#  thetadecl2 <- spm(datadecl2$SP500, b = 150)
#  res[i,4] <- thetadecl2$theta_sl[2]
#}

res <- readRDS("SimStudy.RDS")
apply(res, 2, mean)
