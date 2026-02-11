# For plotting
library(ggplot2)

# For the block bootstrap
library(boot)

### Functions needed for the different Hill estimators
Hill <- function(datas, k, a){ 
  temp <- lapply(a, function(i) log(datas[1:k]/datas[k+1])^i)
  return(sapply(temp, sum)/k)
}

HillcorrRho <- function(data, krmax = NULL, fixrho = FALSE){
  if(is.null(krmax)){
    m <- length(which(data > 0))
    krmax <- floor(min(m-1,(2*m)/(log(log(m)))))
  }
  krseq <- c(1:krmax)
  datas <- rev(sort(data))
  temp <- sapply(krseq, function(k) Hill(datas, k, c(1:4)))
  if(fixrho){
    rho <- -1
  } else{
    numS2 <- 3*(temp[4,] - 24*temp[1,]^4)*(temp[2,] - 2*temp[1,]^2)
    denomS2 <- 4*(temp[3,] - 6*temp[1,]^3)^2 #missing square in paper
    S2 <- numS2/denomS2
    indx <- which(S2 >= 2/3 & S2 <= 3/4)
    rho <- (-4 + 6*S2[indx] + sqrt(3*S2[indx] - 2))/(4*S2[indx] -3)
  }
  return(cbind(indx,rho))
}

HillBoot <- function(data, kseq){
  datas <- rev(sort(data))
  hill <- sapply(kseq, function(k) Hill(datas, k, 1))
  return(1/hill)
}

HillcorrBoot <- function(data, kseq, rho){
  datas <- rev(sort(data))
  temp <- sapply(kseq, function(k) Hill(datas, k, c(1:2)))
  rhofunc <- rho/(1 - rho)
  hillcorr <- temp[1,] - (temp[2,] - 2*temp[1,]^2)/(2*temp[1,]*rhofunc)
  return(1/hillcorr)
}

### Load the daily negative log-returns
load("sp500.RData")
n <- nrow(data)

rho <- -1
### Motivation for choosing rho = -1
chooseRho <- HillcorrRho(data$SP500, krmax = 5000)
# all values
plot(chooseRho, ylim = c(-10,1), type = "l") 
# zoom in:
plot(chooseRho[chooseRho[,1] >= 600,], ylim = c(-4,0.5), type = "l") #stable region for krho in [750,1250]
# gives a result very close to -1: 
chooseRho[chooseRho[,1] == 1250,2]


### Takes a couple of minutes 
kseq <- c(1:2500)
hillboottemp <- tsboot(data$SP500, statistic = HillBoot, R = 500, kseq = kseq,
                   sim = "geom", l = 200)
hillboot <- cbind(hillboottemp$t0, apply(hillboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                  apply(hillboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))
hillcorrboottemp <- tsboot(data$SP500, statistic = HillcorrBoot, rho = rho, R = 500, kseq = kseq,
                       sim = "geom", l = 200)
hillcorrboot <- cbind(hillcorrboottemp$t0, apply(hillcorrboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                  apply(hillcorrboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))

k <- c(50:2500)
alpha <- data.frame('k' = k, 'hill' = hillboot[k,1], 'hillL' = hillboot[k,2],
                    'hillU' = hillboot[k,3],'hillcorr' = hillcorrboot[k,1],
                    'hillcorrL' = hillcorrboot[k,2],'hillcorrU' = hillcorrboot[k,3])

ggplot(alpha, aes(x = k, y = hill)) + 
  geom_line(aes(y =  hill, color = 'Standard Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= hillL, ymax = hillU, color = 'Standard Hill'), show.legend = FALSE, alpha = 0.125, fill = 'blue', linetype = 'dotted') +
  coord_cartesian(ylim=c(1.5,4.25)) + xlab(expression(k[alpha])) + 
  ylab("Tail index") + theme_bw() + 
  theme(axis.title.x =element_text(size=18), 
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))  +
  geom_line(aes(y =  hillcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= hillcorrL, ymax = hillcorrU, color = 'Corrected Hill'),  show.legend = FALSE, linetype = 'dotted', alpha = 0.125) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.025,0.15), legend.key.width = unit(3, "line"))

### obtain estimates of alpha
alphaHill <- alpha[alpha[,1] == 250,2]
alphaHillCorr <- alpha[alpha[,1] == 1000,5]

### Weissman estimator
quantEstim <- function(data, alpha, kseq, nsize, p = 0.99){
  Xnk <- sort(data)[nsize-kseq]
  est <- Xnk*(kseq/(nsize*(1-p)))^(1/alpha)
  return(est)
}

kseq <- c(1:1000)
quantboottemp <- tsboot(data$SP500, statistic = quantEstim, R = 500, kseq = kseq,
                       alpha = alphaHill, nsize = n, sim = "geom", l = 200)
quantboot <- cbind(quantboottemp$t0, apply(quantboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                  apply(quantboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))
quantcorrboottemp <- tsboot(data$SP500, statistic = quantEstim, R = 500, kseq = kseq,
                            alpha = alphaHillCorr, nsize = n, sim = "geom", l = 200)
quantcorrboot <- cbind(quantcorrboottemp$t0, apply(quantcorrboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                      apply(quantcorrboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))

kq <- c(75:1000)
quanttot <- data.frame('k' = kq, 'quant' = quantboot[kq,1], 'quantL' = quantboot[kq,2],
                    'quantU' = quantboot[kq,3],'quantcorr' = quantcorrboot[kq,1],
                    'quantcorrL' = quantcorrboot[kq,2],'quantcorrU' = quantcorrboot[kq,3])
quantemp <- quantile(data$SP500,0.99)

ggplot(quanttot, aes(x = kq, y = quant)) + theme_bw() + 
  geom_line(aes(y =  quant, color = 'Standard Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= quantL, ymax = quantU, color = 'Standard Hill'), show.legend = FALSE, alpha = 0.125, fill = 'blue', linetype = 'dotted') +
  coord_cartesian(ylim=c(2.35,3.4)) + xlab("k") +
  theme_bw() + ylab("99% quantile") +
  theme(axis.title.x =element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))  +
  geom_hline(aes(yintercept = quantemp, color = 'Empirical'), lwd = 0.9) +
  geom_line(aes(y =  quantcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= quantcorrL, ymax = quantcorrU, color = 'Corrected Hill'),  show.legend = FALSE, linetype = 'dotted', alpha = 0.125) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Corrected Hill' = 'black',
                                                  'Empirical' = 'forestgreen')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.5,0.99), legend.key.width = unit(3, "line"))

### Rolling windows
k1 <- 50
k2 <- 200
kquant <- 50
yearly <- split(data, format(data$Date, "%Y"))
m <- length(yearly)
alpha <- matrix(0, nrow = (m-8), ncol = 2)
estim <- excess <- matrix(0, nrow = (m-8), ncol = 3)
excess8year <- matrix(0, nrow = (m-15), ncol = 3)

for(i in 1:(m-8)){
  datat <- do.call(rbind, yearly[i:(i+7)])
  alpha[i,1] <- HillBoot(datat$SP500, k1)  
  alpha[i,2] <- HillcorrBoot(datat$SP500, k2, rho)  
  
  nsize <- length(datat$SP500)
  estim[i,1] <- quantEstim(datat$SP500, alpha = alpha[i,1], kquant, nsize)
  estim[i,2] <- quantEstim(datat$SP500, alpha = alpha[i,2], kquant, nsize)
  estim[i,3] <- quantile(datat$SP500, 0.99)
  excess[i,1] <- length(which(yearly[[i+8]]$SP500 > estim[i,1])) 
  excess[i,2] <- length(which(yearly[[i+8]]$SP500 > estim[i,2]))
  excess[i,3] <- length(which(yearly[[i+8]]$SP500 > estim[i,3]))
  if(i <= m-15){
    datat8year <- do.call(rbind, yearly[(i+8):(i+15)])
    excess8year[i,1] <- length(which(datat8year$SP500 > estim[i,1])) 
    excess8year[i,2] <- length(which(datat8year$SP500 > estim[i,2])) 
    excess8year[i,3] <- length(which(datat8year$SP500 > estim[i,3])) 
  }
}

quant <- data.frame('Year' = c(1969:2022), 'quantHill' = estim[,1],
                    'quantHillcorr' = estim[,2], 'quantEmp' = estim[,3])

ggplot(quant, aes(x = Year, y = quantHill)) + 
  geom_line(aes(y =  quantHill, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(1.5,5)) + ylab("99% quantile") +
  theme_bw() + ggtitle("Rolling window estimates") + xlab("Time (years)") + 
  theme(plot.title = element_text(size=17, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  quantHillcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  quantEmp, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen',
                                                  'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=16), legend.position = c(0.03,0.98), legend.key.width = unit(3, "line"))


excessmat <- data.frame('Year' = c(1969:2022), 'excessHill' = excess[,1],
                        'excessHillcorr' = excess[,2], 'excessEmp' = excess[,3])

ggplot(excessmat, aes(x = Year, y = excessHill)) + 
  geom_line(aes(y =  excessHill, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(0,30)) + xlab("Time (years)") + 
  theme_bw() + ggtitle("Next year") + ylab("Number of exceedances") + 
  theme(plot.title = element_text(size=17, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  excessHillcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  excessEmp, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen',
                                                  'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=16), legend.position = c(0.2,0.975), legend.key.width = unit(3, "line"))

excessmat8 <- data.frame('Year' = c(1969:2015), 'excessHill' = excess8year[,1],
                         'excessHillcorr' = excess8year[,2], 'excessEmp' = excess8year[,3])

ggplot(excessmat8, aes(x = Year, y = excessHill)) + 
  geom_line(aes(y =  excessHill, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(0,120)) + xlab("Time (years)") +
  theme_bw() + ggtitle("Next 8 years") + ylab("Number of exceedances") +
  theme(plot.title = element_text(size=17, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  excessHillcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  excessEmp, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen',
                                                  'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=16), legend.position = c(0.1,0.975), legend.key.width = unit(3, "line"))

apply(excess, 2, mean)
apply(excess8year, 2, mean)

