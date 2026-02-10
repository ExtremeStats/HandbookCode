# This is the main file that calls MCMC and stores the results

rm(list=ls())

setwd("S:\\Documents\\My Papers\\EVA_Handbook\\BayesEVA\\Streamflow\\")
source("nimble_fx.R")

# Load the streamflow data
file <- "https://www4.stat.ncsu.edu/~bjreich/ST740/HCDN_annual_max.RData"
load(url(file))

y          <- log(Y+1)    # Transform the data
ns         <- nrow(y)     # Number of stations
nt         <- ncol(y)     # Number of years  
x          <- (1:nt)/10   # Time covariates
x          <- x - mean(x)

model      <- 3           # 1 = full model; 2 = model without REs; 3 = separate fits by site 
iters      <- 25000       # MCMC iterations
burn       <- 5000        # MCMC burn-in
thin       <- 10          # MCMC thinning
HUC        <- as.factor(HUC02) # HUC indicators
huc        <- as.numeric(HUC)
H          <- max(huc)

# Run the model

 print(paste("model =",model))
 source(paste0("model",model,".R"))
 print(output$WAIC)

# Save the output

 save(output,code,pri_sd,file=paste0("model",model,".RData"))
