# For plotting:
library(ggplot2)

# For times series models
library(rugarch)

# For the extremal index 
library(exdex)

# For the block bootstrap
library(boot)

### Functions needed for the different Hill estimators
Hill <- function(datas, k, a){ 
  temp <- lapply(a, function(i) log(datas[1:k]/datas[k+1])^i)
  return(sapply(temp, sum)/k)
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

### Load the daily negative log-returns and fit a AR(1)-GARCH(1,1) model
load("sp500.RData")
n <- nrow(data)
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)))
gfit <- ugarchfit(model, data$SP500) 
eps <- gfit@fit$z 

### Estimating the extremal index (takes about a minute)
bval <- seq(15, 300, by = 5)
res <- choose_b(eps, bval, interval_type = "lik")
theta <- data.frame('b' = bval, 'thest' = res$theta_sl[,2])

ggplot(theta, aes(x = b, y = thest)) + 
  geom_line(aes(y =  thest), lwd = 0.75) +
  geom_ribbon(aes(ymin= res$lower_sl[,2], ymax = res$upper_sl[,2]), color = 'black', show.legend = FALSE, linetype = 'dotted', alpha = 0.125, fill = 'black') +
  coord_cartesian(ylim=c(0.6,1)) + xlab("Block size") + 
  theme_bw() + ylab("Extremal index") +
  theme(axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) 

theta <- spm(eps, b = 75)
confint(theta, interval_type = "lik")

### Takes a couple of minutes
kseq <- c(1:2500)
rho <- -1
hillboottemp <- tsboot(eps, statistic = HillBoot, R = 500, kseq = kseq,
                       sim = "geom", l = 200)
hillboot <- cbind(hillboottemp$t0, apply(hillboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                  apply(hillboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))
hillcorrboottemp <- tsboot(eps, statistic = HillcorrBoot, rho = rho, R = 500, kseq = kseq,
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
  coord_cartesian(ylim=c(2,6)) + xlab(expression(k[alpha])) + 
  theme_bw() + ylab("Tail index") + 
  theme(axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))  +
  geom_line(aes(y =  hillcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= hillcorrL, ymax = hillcorrU, color = 'Corrected Hill'),  show.legend = FALSE, linetype = 'dotted', alpha = 0.125) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.49,0.99), legend.key.width = unit(3, "line"))

### obtain estimates of alpha
alphaHill <- alpha[alpha[,1] == 250,2]
alphaHillCorr <- alpha[alpha[,1] == 1500,5]

### Weissman estimator
quantEstim <- function(data, alpha, kseq, nsize, p = 0.99){
  Xnk <- sort(data)[nsize-kseq]
  est <- Xnk*(kseq/(nsize*(1-p)))^(1/alpha)
  return(est)
}

kseq <- c(1:1000)
quantboottemp <- tsboot(eps, statistic = quantEstim, R = 500, kseq = kseq,
                        alpha = alphaHill, nsize = n, sim = "geom", l = 200)
quantboot <- cbind(quantboottemp$t0, apply(quantboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                   apply(quantboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))
quantcorrboottemp <- tsboot(eps, statistic = quantEstim, R = 500, kseq = kseq,
                            alpha = alphaHillCorr, nsize = n, sim = "geom", l = 200)
quantcorrboot <- cbind(quantcorrboottemp$t0, apply(quantcorrboottemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                       apply(quantcorrboottemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))

kq <- c(50:750)
quanttot <- data.frame('k' = kq, 'quant' = quantboot[kq,1], 'quantL' = quantboot[kq,2],
                       'quantU' = quantboot[kq,3],'quantcorr' = quantcorrboot[kq,1],
                       'quantcorrL' = quantcorrboot[kq,2],'quantcorrU' = quantcorrboot[kq,3])
quantemp <- quantile(eps,0.99)

ggplot(quanttot, aes(x = kq, y = quant)) + theme_bw() + 
  geom_line(aes(y =  quant, color = 'Standard Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= quantL, ymax = quantU, color = 'Standard Hill'), show.legend = FALSE, alpha = 0.125, fill = 'blue', linetype = 'dotted') +
  coord_cartesian(ylim=c(2.35,2.8)) + xlab('k') +
  theme_bw() + ylab("99% quantile") +
  theme(axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))  +
  geom_hline(aes(yintercept = quantemp, color = 'Empirical'), lwd = 0.9) +
  geom_line(aes(y =  quantcorr, color = 'Corrected Hill'), lwd = 0.9) +
  geom_ribbon(aes(ymin= quantcorrL, ymax = quantcorrU, color = 'Corrected Hill'),  show.legend = FALSE, linetype = 'dotted', alpha = 0.125) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Corrected Hill' = 'black',
                                                  'Empirical' = 'forestgreen')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.49,0.99), legend.key.width = unit(3, "line"))

### Rolling windows
k1 <- 50
k2 <- 200
kquant <- 50

window_size <- 2000
test_window <- 250
test_window8year <- 2000
alpha <- alphaG <- matrix(0, nrow = (n - window_size), ncol = 2)
estim <- estimG <- matrix(0, nrow = (n - window_size), ncol = 3)

### this function is a bit faster that the "boot" versions
Hillcorr <- function(datasorted, k, rho){
  temp <- Hill(datasorted, k, c(1:2))
  rhofunc <- rho/(1 - rho)
  hillcorr <- temp[1] - (temp[2] - 2*temp[1]^2)/(2*temp[1]*rhofunc)
  return(1/hillcorr)
}

#### Takes 30 minutes; results are stored in estim.RDS and estimG.RDS

# for(i in 1:(n - window_size)){
#   datat <- data$SP500[i:(i+window_size-1)]
#   gfit <- ugarchfit(model, datat, solver = 'hybrid') 
#   eps <- gfit@fit$z 
#   frcst <- ugarchforecast(gfit, n.ahead = 1)
#   mu <- as.numeric(fitted(frcst))
#   sig <- as.numeric(sigma(frcst))
#   
#   datasortedG <- sort(eps, decreasing = T)
#   alphaG[i,1] <- 1/Hill(datasortedG, k1, 1)  
#   alphaG[i,2] <- Hillcorr(datasortedG, k2, rho)  
#   
#   datasorted <- sort(datat, decreasing = T)
#   alpha[i,1] <- 1/Hill(datasorted, k1, 1)  
#   alpha[i,2] <- Hillcorr(datasorted, k2, rho)  
#   
#   estimG[i,1] <- mu + sig*quantEstim(eps, alpha = alphaG[i,1], kquant, window_size)
#   estimG[i,2] <- mu + sig*quantEstim(eps, alpha = alphaG[i,2], kquant, window_size)
#   estimG[i,3] <- mu + sig*quantile(eps, 0.99)
#   estim[i,1] <- quantEstim(datat, alpha = alpha[i,1], kquant, window_size)
#   estim[i,2] <- quantEstim(datat, alpha = alpha[i,2], kquant, window_size)
#   estim[i,3] <- quantile(datat, 0.99)
# }

estimG <- readRDS("quantCond.RDS")
estim <- readRDS("quantUncond.RDS")

garchuc_test <- garchcc_test <- uc_test <- cc_test <- matrix(1, nrow = (n - test_window - window_size+1), ncol = 3)
excess <- excessG <- matrix(NA, nrow = (n - test_window - window_size+1), ncol = 3)
for (i in 1:(n- test_window - window_size+1)){
  datat <- data$SP500[(i+window_size):(i+window_size+test_window-1)]
  VaRt <- estim[i:(i+test_window-1), ]
  VaRtG <- estimG[i:(i+test_window-1), ]
  excess[i,] <- apply(VaRt, 2, function(j) sum(datat > j))
  excessG[i,] <- apply(VaRtG, 2, function(j) sum(datat > j))
  for(j in 1:3){
    if(excess[i,j]>1){  
        res <- VaRTest(alpha = 0.01, -datat, -VaRt[,j])
        uc_test[i,j] <- res$uc.LRp
        cc_test[i,j] <- res$cc.LRp
    }
    if(excessG[i,j]>1){  
      resG <- VaRTest(alpha = 0.01, -datat, -VaRtG[,j])
      garchuc_test[i,j] <- resG$uc.LRp
      garchcc_test[i,j] <- resG$cc.LRp
    }
  }
}  

# Takes a couple of minutes
garchuc_test8 <- garchcc_test8 <- uc_test8 <- cc_test8 <- matrix(1, nrow = (n - test_window8year - window_size+1), ncol = 3)
excess8 <- excessG8 <- matrix(NA, nrow = (n - test_window8year - window_size+1), ncol = 3)
for (i in 1:(n- test_window8year - window_size+1)){
  datat <- data$SP500[(i+window_size):(i+window_size+test_window8year-1)]
  VaRt <- estim[i:(i+test_window8year-1), ]
  VaRtG <- estimG[i:(i+test_window8year-1), ]
  excess8[i,] <- apply(VaRt, 2, function(j) sum(datat > j))
  excessG8[i,] <- apply(VaRtG, 2, function(j) sum(datat > j))
  for(j in 1:3){
    if(excess8[i,j]>1){  
      res <- VaRTest(alpha = 0.01, -datat, -VaRt[,j])
      uc_test8[i,j] <- res$uc.LRp
      cc_test8[i,j] <- res$cc.LRp
    }
    if(excessG8[i,j]>1){  
      resG <- VaRTest(alpha = 0.01, -datat, -VaRtG[,j])
      garchuc_test8[i,j] <- resG$uc.LRp
      garchcc_test8[i,j] <- resG$cc.LRp
    }
  }
}  

exceedances <- data.frame('Date' = data$Date[(test_window+window_size):n], 
                          'excess1' = excess[,1],'excess2' = excess[,2],'excess3' = excess[,3],
                          'excessG1' = excessG[,1],'excessG2' = excessG[,2], 'excessG3' = excessG[,3])
apply(exceedances[,-1],2,mean)
apply(exceedances[,-1],2,max)

ggplot(exceedances, aes(x = Date, y = excess1)) + 
  geom_line(aes(y =  excess1, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(0,30)) + xlab("Time (years)") +
  ylab("Number of exceedances") +
  theme_bw() + ggtitle("Unconditional risk analysis") +  
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  excess2, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  excess3, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen', 'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.2,0.975), legend.key.width = unit(3, "line"))

ggplot(exceedances, aes(x = Date, y = excessG1)) + 
  geom_line(aes(y =  excessG1, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(0,15)) + xlab("Time (years)") +
  ylab("Number of exceedances") +
  theme_bw() + ggtitle("Conditional risk analysis") +  
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  excessG2, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  excessG3, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen', 'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.2,0.975), legend.key.width = unit(3, "line"))



exceedances8 <- data.frame('Date' = data$Date[(test_window8year+window_size):n], 
                           'excess1' = excess8[,1],'excess2' = excess8[,2],'excess3' = excess8[,3],
                           'excessG1' = excessG8[,1],'excessG2' = excessG8[,2], 'excessG3' = excessG8[,3])
apply(exceedances8[,-1],2,mean)
apply(exceedances8[,-1],2,max)

ggplot(exceedances8, aes(x = Date, y = excess1)) + 
  geom_line(aes(y =  excess1, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(0,60)) +
  theme_bw() + ggtitle("Unconditional risk analysis") +  
  ylab("Number of exceedances") + xlab("Time (years)") + 
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  excess2, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  excess3, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen', 'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.03,0.99), legend.key.width = unit(3, "line"))


ggplot(exceedances8, aes(x = Date, y = excessG1)) + 
  geom_line(aes(y =  excessG1, color = 'Standard Hill'), lwd = 0.9) +
  coord_cartesian(ylim=c(0,40)) +
  theme_bw() + ggtitle("Conditional risk analysis") +  
  ylab("Number of exceedances") + xlab("Time (years)") + 
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y =element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15)) +
  geom_line(aes(y =  excessG2, color = 'Corrected Hill'), lwd = 0.9) +
  geom_line(aes(y =  excessG3, color = 'Empirical'), lwd = 0.9) +
  scale_color_manual(name = "Legend", values =  c('Standard Hill' = 'blue', 'Empirical' = 'forestgreen', 'Corrected Hill' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=20), legend.position = c(0.03,0.24), legend.key.width = unit(3, "line"))


### Section 5: backtesting

uc_tests <- data.frame('Date' = data$Date[(test_window+window_size):n], 
                       'uc1' = uc_test[,1],'uc2' = uc_test[,2], 'uc3' = uc_test[,3],
                       'ucG1' = garchuc_test[,1],'ucG2' = garchuc_test[,2], 'ucG3' = garchuc_test[,3])

round(apply(uc_tests[,-1], 2, function(i) length(which(i < 0.05)))/nrow(uc_tests),3)


cc_tests <- data.frame('Date' = data$Date[(test_window+window_size):n], 
                       'cc1' = cc_test[,1],'cc2' = cc_test[,2], 'cc3' = cc_test[,3],
                       'ccG1' = garchcc_test[,1],'ccG2' = garchcc_test[,2], 'ccG3' = garchcc_test[,3])

round(apply(cc_tests[,-1], 2, function(i) length(which(i < 0.05)))/nrow(cc_tests),3)

###### 8 year

uc_tests8 <- data.frame('Date' = data$Date[(test_window8year+window_size):n], 
                        'uc1' = uc_test8[,1],'uc2' = uc_test8[,2], 'uc3' = uc_test8[,3],
                        'ucG1' = garchuc_test8[,1],'ucG2' = garchuc_test8[,2], 'ucG3' = garchuc_test8[,3])

round(apply(uc_tests8[,-1], 2, function(i) length(which(i < 0.05)))/nrow(uc_tests8),3)

cc_tests8 <- data.frame('Date' = data$Date[(test_window8year+window_size):n], 
                        'cc1' = cc_test8[,1],'cc2' = cc_test8[,2], 'cc3' = cc_test8[,3],
                        'ccG1' = garchcc_test8[,1],'ccG2' = garchcc_test8[,2], 'ccG3' = garchcc_test8[,3])

round(apply(cc_tests8[,-1], 2, function(i) length(which(i < 0.05)))/nrow(cc_tests8),3)


