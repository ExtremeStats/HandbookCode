## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(JM.landslide)  # Contains functions for joint modeling of landslide counts and sizes data.
library(spdep)         # Provides functions for spatial dependence analysis.
library(INLA)          # Used here for creating the neighborhood structures.
library(ggplot2)
library(gridExtra)
library(viridis)

# The dataset contains the necessary spatial information, such as the shapefile with spatial locations and coordinates.
data("Wenchuan_info_used_for_simulation")

## -----------------------------------------------------------------------------
nb2INLA("adjgraph-sim.txt", poly2nb(shp_selected, queen = FALSE, row.names = shp_selected$SU_ID))
adjacensy <- inla.read.graph(filename = "adjgraph-sim.txt")
N <- adjacensy$n
diag.Q <- diag(adjacensy$nnbs, N)
A.Q <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  A.Q[i, adjacensy$nbs[[i]]] <- 1
}
Q <- diag.Q - A.Q  # Precision matrix

## -----------------------------------------------------------------------------
# Set fixed parameters
kappa_w1 <- 10
kappa_w2 <- 5
kappa_eta <- 10
kappa_mu <- 5
intercept1 <- 2
intercept2 <- 4
beta <- 1
other.hyper <- list(
  kappa_w1 = kappa_w1,
  kappa_w2 = kappa_w2,
  kappa_eta = kappa_eta,
  kappa_mu = kappa_mu,
  intercept1 = intercept1,
  intercept2 = intercept2,
  beta = beta
)
hyper.mu<- c(20,0.1)

set.seed(1)
Z1 <- mvtnorm::rmvnorm(N, mean = rep(0, 3))
Z2 <- mvtnorm::rmvnorm(N, mean = rep(0, 2))
beta1 <- runif(ncol(Z1), -1, 1)
beta2 <- runif(ncol(Z2), -1, 1)

## ----warning=FALSE------------------------------------------------------------
set.seed(1)
Sim_mark_data <- JM.landslide::sim_mark_function(
  Q = Q,
  hyper.mu =  hyper.mu,
  other.hyper = other.hyper,
  beta1 = beta1, beta2 = beta2,
  Z1 = as.matrix(Z1), Z2 = as.matrix(Z2),
  mark_dist = "eGPD",
  model_type = "jSp"  
)
# Access results
Y <- Sim_mark_data$Y
A <- Sim_mark_data$A
mu <- Sim_mark_data$mu
eta <- Sim_mark_data$eta
W1 <- Sim_mark_data$W1
W2 <- Sim_mark_data$W2

## ----warning=FALSE, fig.width = 7, fig.height = 6-----------------------------
generate_histogram <- function(data, var_name, binwidth = 30, title = NULL, fill_color = "blue", y_label = "Frequency") {
  ggplot(data, aes_string(x = var_name)) +
    geom_histogram(binwidth = binwidth, fill = fill_color, alpha = 0.7, color = "black") +
    labs(title = title, y = y_label) +  # Add custom y-axis label
    theme_minimal()
}
Sim_data <- data.frame(Y = Sim_mark_data$Y, 
                            A = Sim_mark_data$A, 
                            mu = Sim_mark_data$mu,
                            eta = Sim_mark_data$eta, 
                            W1 = Sim_mark_data$W1, 
                            W2 = Sim_mark_data$W2)
p1 <- generate_histogram(Sim_data, var_name = "Y", binwidth = 30, title = "Histogram of Y", fill_color = "blue", y_label = "Frequency")
p2 <- generate_histogram(Sim_data, var_name = "A", binwidth = 30, title = "Histogram of A", fill_color = "green", y_label = "Frequency")
p3 <- generate_histogram(Sim_data, var_name = "eta", binwidth = 0.5, title = "Histogram of eta", fill_color = "red", y_label = "Frequency")
p4 <- generate_histogram(Sim_data, var_name = "mu", binwidth = 0.5, title = "Histogram of mu", fill_color = "purple", y_label = "Frequency")
p5 <- generate_histogram(Sim_data, var_name = "W1", binwidth = 0.5, title = "Histogram of W1", fill_color = "orange", y_label = "Frequency")
p6 <- generate_histogram(Sim_data, var_name = "W2", binwidth = 0.5, title = "Histogram of W2", fill_color = "yellow", y_label = "Frequency")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)


## ----warning=FALSE------------------------------------------------------------
outputs <- JM.landslide::mcmc_sampler(
  Y = Y,               
  A = A,               
  Z1 = Z1,                           
  Z2 = Z2,                           
  CV = "WS",                         
  mark_dist = "eGPD",                
  thr.family = "gamma",              
  model_type = "jSp",                
  adjacensy = adjacensy,            
  no.rm.obs = 200,                 
  N.MCMC = 2000,                     
  samples.store = 250,              
  print.result = FALSE            
)

## ----warning=FALSE------------------------------------------------------------
outputs$JM.info$summry_hyper

## ----warning=FALSE------------------------------------------------------------
outputs$JM.info$summ_fixed_effects

## -----------------------------------------------------------------------------
coords<- cbind(shp_selected$POINT_X, shp_selected$POINT_Y)
WS_results.Y<- data.frame(true.Y = outputs$JM.info$qqplots_with_CIs$WS$true.Y, 
                          est.Y = outputs$JM.info$qqplots_with_CIs$WS$est.Y,
                          lci.Y = outputs$JM.info$qqplots_with_CIs$WS$lci.Y,
                          uci.Y = outputs$JM.info$qqplots_with_CIs$WS$uci.Y
)

## ----fig.width=7, fig.height=4------------------------------------------------
range<- c(min(WS_results.Y$true.Y, WS_results.Y$est.Y, WS_results.Y$lci.Y, WS_results.Y$uci.Y),
          max(WS_results.Y$true.Y, WS_results.Y$est.Y, WS_results.Y$lci.Y, WS_results.Y$uci.Y))

p.WS.Y <- ggplot(WS_results.Y, aes(x = est.Y, y = true.Y)) +
  geom_point(size = 0.5) +  # Smaller points for estimates
  geom_ribbon(aes(ymin = lci.Y, ymax = uci.Y), fill = "blue", alpha = 0.2) +  # Shaded area for confidence intervals
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +  # Diagonal line for reference
  labs(title = "count (WS)", x = "Estimated Quantile", y = "Data Quantile") +
  xlim(range) + ylim(range) +
  #theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))
## For sizes
WS_results.A<- data.frame(true.A = outputs$JM.info$qqplots_with_CIs$WS$true.A[outputs$JM.info$qqplots_with_CIs$WS$true.A>0], 
                          est.A =  outputs$JM.info$qqplots_with_CIs$WS$est.A,
                          lci.A =  outputs$JM.info$qqplots_with_CIs$WS$lci.A,
                          uci.A = outputs$JM.info$qqplots_with_CIs$WS$uci.A
)

range<- c(min(WS_results.A$true.A, WS_results.A$est.A, WS_results.A$lci.A, WS_results.A$uci.A),
          max(WS_results.A$true.A, WS_results.A$est.A, WS_results.A$lci.A, WS_results.A$uci.A))
p.WS.A <- ggplot(WS_results.A, aes(x = est.A, y = true.A)) +
  geom_point(size = 0.5) +  # Smaller points for estimates
  geom_ribbon(aes(ymin = lci.A, ymax = uci.A), fill = "blue", alpha = 0.2) +  # Shaded area for confidence intervals
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +  # Diagonal line for reference
  labs(title = "size (WS)", x = "Estimated Quantile", y = "Data Quantile") +
  xlim(range) + ylim(range) +
  #theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))

grid.arrange(p.WS.Y, p.WS.A, ncol = 2)

## -----------------------------------------------------------------------------
post_mean_w1<- outputs$JM.info$est_random_effects$W1$mean 
post_mean_w2<- outputs$JM.info$est_random_effects$W2$mean 
post_std_w1<- outputs$JM.info$est_random_effects$W1$sd
post_std_w2<-  outputs$JM.info$est_random_effects$W2$sd

## ----fig.width=7, fig.height=4------------------------------------------------
# Data frame for plotting
df <- data.frame(
  lat = coords[,1], 
  lon = coords[,2], 
  post_mean_w1 = post_mean_w1,
  post_mean_w2 = post_mean_w2,
  post_std_w1 = post_std_w1,
  post_std_w2 = post_std_w2
)

# W1 mean
p1 <- ggplot(df, aes(x = lat, y = lon, col = post_mean_w1)) +
  geom_point() + 
  ggtitle("W1 (mean)") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = 'mean') +
  scale_color_viridis()

#  W1 standard deviation
p2 <- ggplot(df, aes(x = lat, y = lon, col = post_std_w1)) +
  geom_point() + 
  ggtitle("W1 (sd)") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = 'sd') +
  scale_color_viridis()

# W2 mean
p3 <- ggplot(df, aes(x = lat, y = lon, col = post_mean_w2)) +
  geom_point() + 
  ggtitle("W2 (mean)") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = 'mean') +
  scale_color_viridis()

# W2 standard deviation
p4 <- ggplot(df, aes(x = lat, y = lon, col = post_std_w2)) +
  geom_point() + 
  ggtitle("W2 (sd)") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = 'sd') +
  scale_color_viridis()

# Arrange the four plots into a 2x2 grid
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)


## -----------------------------------------------------------------------------
est.counts<- outputs$JM.info$est_counts_and_sizes$WS$post.mean.Y
est.sizes<- rep(0, length(est.counts))
est.sizes[outputs$data$size.full==0]<- NA
est.sizes[!(outputs$data$size.full==0)]<- outputs$JM.info$est_counts_and_sizes$WS$post.mean.A 
est.log.hazard<- outputs$JM.info$est_hazards$mean
est.suscept<- outputs$JM.info$est_susceptibility$mean

## ----fig.width=7, fig.height=4------------------------------------------------
# Load the data
df <- data.frame(
  lat = coords[,1], 
  lon = coords[,2], 
  log.counts = log(est.counts), 
  log.areas = log(est.sizes),
  suscep = est.suscept,
  log.hazard = est.log.hazard
)

# Plot 1: estimated log-count
p1 <- ggplot(df, aes(x = lat, y = lon, col = log.counts)) +
  geom_point() + 
  ggtitle("Estimated Log-Count") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = 'count') +
  scale_color_viridis()

# Plot 2: estimated log-size
p2 <- ggplot(df, aes(x = lat, y = lon, col = log.areas)) +
  geom_point() + 
  ggtitle("Estimated Log-Size") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = 'size') +
  scale_color_viridis()

# Plot 3: susceptibility
p3 <- ggplot(df, aes(x = lat, y = lon, col = suscep)) +
  geom_point() + 
  ggtitle("Susceptibility") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = "prob") +
  scale_color_viridis()

# Plot 4: log-hazard
p4 <- ggplot(df, aes(x = lat, y = lon, col = log.hazard)) +
  geom_point() + 
  ggtitle("Log-Hazard") +
  xlab(expression("Lon ["*degree*"]")) + 
  ylab(expression("Lat ["*degree*"]")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  labs(color = "hazard") +
  scale_color_viridis()
grid.arrange(p1, p2, p3, p4, ncol = 2)


