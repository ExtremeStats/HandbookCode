# For plotting
library(ggplot2)

# For times series models
library(rugarch)

# For the block bootstrap
library(boot)

#### data_bivar1: SP500 and Dow Jones Industrial Average. Only the dates for which both prices are available are retained
#### data_bivar2: SP500 and FTSE100. Only the dates for which both prices are available are retained

load("data_bivar1.RData")
n1 <- nrow(data1)
load("data_bivar2.RData")
n2 <- nrow(data2)


ggplot(data1, aes(x = SP500, y = DJIA)) + 
  geom_point() + theme_bw() + 
  coord_cartesian(ylim=c(-15,15), xlim=c(-15,15)) +  
  xlab("S&P500") + ylab("DJIA") + 
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y = element_text(size=15)) 

ggplot(data2, aes(x = SP500, y = FTSE)) + 
  geom_point() + theme_bw() + 
  coord_cartesian(ylim=c(-15,15), xlim=c(-15,15)) +  
  xlab("S&P500") + ylab("FTSE") + 
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y = element_text(size=15)) 

### Fitting the AR(1)-GARCH(1,1) model and extracting the residuals
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1, 0)))
gfit1 <- ugarchfit(model, data1$SP500) 
gfit2 <- ugarchfit(model, data1$DJIA) 
gfit3 <- ugarchfit(model, data2$SP500) 
gfit4 <- ugarchfit(model, data2$FTSE) 
eps1 <- gfit1@fit$z
eps2 <- gfit2@fit$z 
eps3 <- gfit3@fit$z 
eps4 <- gfit4@fit$z 

resid1 <- data.frame('SP500' =  eps1,'DJIA' = eps2)
resid2 <- data.frame('SP500' =  eps3,'FTSE' = eps4)

ggplot(resid1, aes(x = SP500, y = DJIA)) + 
  geom_point() + theme_bw() + 
  coord_cartesian(ylim=c(-8,8), xlim=c(-8,8)) +  
  xlab("S&P500") + ylab("DJIA") + 
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y = element_text(size=15)) 

ggplot(resid2, aes(x = SP500, y = FTSE)) + 
  geom_point() + theme_bw() + 
  coord_cartesian(ylim=c(-8,8), xlim=c(-8,8)) +  
  xlab("S&P500") + ylab("FTSE") + 
  theme(axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        axis.text.x = element_text(size=15),  
        axis.text.y = element_text(size=15)) 

# Function to calculate an estimate of chi
chiEmp <-function(data, kseq){
  n <- nrow(data)
  ranks <- apply(data, 2, rank)
  minranks <- apply(ranks, 1, min)
  res <- sapply(kseq, function(i) sum(minranks > n - i)/i)
  return(res)
}

kseq <- seq(20,1500, by = 5)
uncondtemp <- tsboot(data1, statistic = chiEmp, R = 500, kseq = kseq,
                  sim = "geom", l = 200)
uncond <- cbind(uncondtemp$t0, apply(uncondtemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                 apply(uncondtemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))
condtemp <- tsboot(resid1, statistic = chiEmp, R = 500, kseq = kseq,
                     sim = "geom", l = 200)
cond <- cbind(condtemp$t0, apply(condtemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                apply(condtemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))

chires1 <- data.frame('k' = kseq, 'uncond' = uncond[,1], 'uncondL' = uncond[,2],
                    'uncondU' = uncond[,3], 'cond' = cond[,1], 'condL' = cond[,2],
                    'condU' = cond[,3])

uncondtemp <- tsboot(data2, statistic = chiEmp, R = 500, kseq = kseq,
                     sim = "geom", l = 200)
uncond <- cbind(uncondtemp$t0, apply(uncondtemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
                apply(uncondtemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))
condtemp <- tsboot(resid2, statistic = chiEmp, R = 500, kseq = kseq,
                   sim = "geom", l = 200)
cond <- cbind(condtemp$t0, apply(condtemp$t, 2, function(i) quantile(i, 0.05, na.rm =T)),
              apply(condtemp$t, 2, function(i) quantile(i, 0.95, na.rm =T)))

chires2 <- data.frame('k' = kseq, 'uncond' = uncond[,1], 'uncondL' = uncond[,2],
                      'uncondU' = uncond[,3], 'cond' = cond[,1], 'condL' = cond[,2],
                      'condU' = cond[,3])

subset(chires1, k == 500)
subset(chires2, k == 500)

ggplot(chires1, aes(x = k, y = uncond)) + 
  geom_line(aes(y =  uncond, color = 'Raw data'), lwd = 0.9) +
  geom_ribbon(aes(ymin= uncondL, ymax = uncondU, color = 'Raw data'), show.legend = FALSE, alpha = 0.125, fill = 'blue', linetype = 'dotted') +
  coord_cartesian(ylim=c(0.64,0.96)) + ylab("Tail dependence coefficient") + theme_bw() + 
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))  +
  ggtitle("S&P500 vs DJIA") + 
  geom_line(aes(y =  cond, color = 'Residuals'), lwd = 0.9) +
  geom_ribbon(aes(ymin= condL, ymax = condU, color = 'Residuals'),  show.legend = FALSE, linetype = 'dotted', alpha = 0.125) +
  scale_color_manual(name = "Legend", values =  c('Raw data' = 'blue', 'Residuals' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.55,0.2), legend.key.width = unit(3, "line"))

ggplot(chires2, aes(x = k, y = uncond)) + 
  geom_line(aes(y =  uncond, color = 'Raw data'), lwd = 0.9) +
  geom_ribbon(aes(ymin= uncondL, ymax = uncondU, color = 'Raw data'), show.legend = FALSE, alpha = 0.125, fill = 'blue', linetype = 'dotted') +
  coord_cartesian(ylim=c(0.01,0.56)) + ylab("Tail dependence coefficient") + theme_bw() + 
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x =element_text(size=18), 
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=15),  
        axis.text.y=element_text(size=15))  +
  ggtitle("S&P500 vs FTSE100") + 
  geom_line(aes(y =  cond, color = 'Residuals'), lwd = 0.9) +
  geom_ribbon(aes(ymin= condL, ymax = condU, color = 'Residuals'),  show.legend = FALSE, linetype = 'dotted', alpha = 0.125) +
  scale_color_manual(name = "Legend", values =  c('Raw data' = 'blue', 'Residuals' = 'black')) +
  theme(legend.justification = c(0,1), legend.title = element_blank(), 
        legend.text = element_text(size=18), legend.position = c(0.55,0.2), legend.key.width = unit(3, "line"))

