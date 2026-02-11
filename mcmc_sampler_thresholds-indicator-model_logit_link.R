#' MCMC Sampler for Threshold Indicator Model Using Logit Link
#'
#' This function implements a Markov Chain Monte Carlo (MCMC) sampler for a threshold indicator model.
#' It handles Gibbs sampling and Metropolis-Hastings steps for parameter estimation, focusing on
#' modeling binary indicators derived from size data.
#'
#' @param N.MCMC Integer. Number of MCMC iterations.
#' @param A Numeric vector. Size data of length `n2`.
#' @param ind.NA Logical vector. Indicator for missing values in `A`.
#' @param CV Character. Cross-validation type: `"WS"` for within-sample or `"OOS"` for out-of-sample.
#' @param Z2 Matrix. Covariates for the size data (dimensions `n2 x q`).
#' @param thin Integer. Thinning interval for MCMC sampling.
#' @param adapt Integer. Adaptation interval for tuning parameter updates.
#' @param burn_in1 Integer. First burn-in period for MCMC sampling.
#' @param burn_in2 Integer. Second burn-in period for MCMC sampling.
#' @param ind_zero Logical vector. Indicator for zero entries in `A`.
#' @param hyper_fixed List. Fixed hyperparameters for the model.
#' @param print.result Logical. If `TRUE`, prints progress during MCMC.
#' @param traceplot Logical. If `TRUE`, generates traceplots for parameter diagnostics.
#' @param true.values Numeric vector. True parameter values for validating simulation experiments.
#' @param simulation Logical. If `TRUE`, runs the function as a simulation experiment.
#' @param nbd_info Matrix. Information on the adjacency structure of spatial units.
#' @param no_of_nbd Integer. Number of neighbors for each spatial unit.
#' @param node.set Matrix. Node connections used in spatial modeling.
#' @param hyper.mu_adapt_seq2 Numeric vector. Adaptation sequence for `mu` parameters.
#' @param mu_adapt_seq2 Numeric vector. Adaptation sequence for `mu`.
#' @param eta_adapt_seq2 Numeric vector. Adaptation sequence for `eta`.
#' @param init.seed Integer. Seed for reproducibility of the random number generation.
#'
#' @return A list containing the following elements:
#' \item{samples}{Matrix. MCMC samples for the model parameters, including hyperparameters and latent variables.}
#' \item{imputed.A.WSD}{Numeric vector. Imputed values of `A` in within-sample diagnostics.}
#' \item{imputed.A.squre}{Numeric vector. Sum of squared imputed values of `A`.}
#' \item{post.sum.mean.mu}{Numeric vector. Posterior sum of `mu` parameter means.}
#' \item{post.sum.squre.mu}{Numeric vector. Posterior sum of squared `mu` parameters.}
#' \item{post.sum.mean.w2}{Numeric vector. Posterior sum of `w2` parameter means.}
#' \item{post.sum.squre.w2}{Numeric vector. Posterior sum of squared `w2` parameters.}
#'
#' @export
#'
#' @examples
#' # Example of how to run the MCMC sampler
#' result <- mcmc_sampler_indicator_model(N.MCMC = 1000, A = size_data, ind.NA = is.na(size_data),
#'                                       CV = "WS", Z2 = covariate_matrix,
#'                                       thin = 10, adapt = 50, burn_in1 = 100,
#'                                       burn_in2 = 200, ind_zero = zero_indicator,
#'                                       hyper_fixed = list(), print.result = TRUE,
#'                                       traceplot = FALSE, true.values = NULL,
#'                                       simulation = TRUE, nbd_info = adjacency_matrix,
#'                                       no_of_nbd = 4, node.set = node_connections,
#'                                       hyper.mu_adapt_seq2 = seq(0.05, 0.5, length.out = 100),
#'                                       mu_adapt_seq2 = seq(0.05, 0.5, length.out = 100),
#'                                       eta_adapt_seq2 = seq(0.05, 0.5, length.out = 100),
#'                                       init.seed = 123)
mcmc_sampler_indicator_model_logit<-function(
    N.MCMC,
    A,
    ind.NA,
    CV,
    Z2,
    thin,
    adapt,
    burn_in1,
    burn_in2,
    tun.mu,
    ind_zero,
    hyper_fixed,
    print.result = TRUE,
    traceplot = FALSE,
    true.values = NULL,
    simulation,
    nbd_info,
    no_of_nbd,
    node.set,
    hyper.mu_adapt_seq2,
    mu_adapt_seq2,
    eta_adapt_seq2,
    init.seed)
{ 
  #browser()
  n2<- nrow(Z2)
  q<-ncol(Z2)
  n1<- nrow(Z2)
  n.Q<- n1
  
  tun.mu<- rep(tun.mu, n1)
  
  Z2.crossprd<- (t(Z2)%*%Z2)
  ####### Imputation in case of within sample diagnostics
  post.sum.mean.mu <- rep(0, times = n1)
  post.sum.squre.mu <- rep(0, times = n1)
  post.sum.mean.w2 <- rep(0, times = n1)
  post.sum.squre.w2 <- rep(0, times = n1)
  
  imputed.A.sum <- rep(0, times = n1)
  imputed.A.sum.squre <- rep(0, times = n1)
  
  
  ## Storing the tuning parameteres
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),ncol=2,byrow=TRUE) ## storing the adaptive scale parameters
  sigma.matrix[1,]<- c(tun.mu[1], tun.mu[2])
  
  init.all.other.param_and.name<- init_fun_all_other_param_indicator_model_logit_link(Z2=Z2, A=A, seed=init.seed)
  
  init.other.param<- init.all.other.param_and.name$init.all.other.param
  cur.samples.kappa_w2<-init.other.param[1] #current samples for \kappa_w2
  cur.samples.kappa_mu<-init.other.param[2] #current samples for \kappa_w2
  cur.samples.intercept2<-init.other.param[3]
  cur.samples.beta2<-init.other.param[(3+1):(4-1+q)]    #current samples for beta2
  cur.samples.mu<-init.other.param[((4+q)):((4-1+q+n1))] #current samples for mu
  cur.samples.w2<-init.other.param[((4+q+n1)):((4-1+q+n1+n.Q))]  #current samples for w2
  cur.samples.w2<- cur.samples.w2-mean(cur.samples.w2) ## making sure that sum of W2 is zero
  
  ### saving the Chain for all the hyperparameters and some latent parmeters
  samples <- matrix(nrow=floor(N.MCMC/thin),ncol=3+q+4,byrow=TRUE) # storing the samples
  samples[1,]<- c(init.other.param[c(1:3, 4:(3+q), (4+q):(3+q+2), (4+q+n1):(3+q+n1+2))]) ## saving only few samples
  
  
  j<-1
  l<-1
  m<-1
  k<-1
  for (i in 1:(N.MCMC-1)){
    # if((i%%(adapt))-1==0 & i< (burn_in1+burn_in2+2)){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
    #   rate.latent<-0
    #   rate.hyper<-0
    # }
    
    if(((i%%(adapt))-1==0) & (i< (burn_in1+burn_in2+2))){ #to calculate the acceptance rate based on only current samples, burning+2 to calculate the acceptance rate after the burning samples
      rate.mu<- rep(0, n1)
    }
    
    if((i%%thin)-1==0){  ## by is the thinning steps and adapt is the number of iterations after which i update the variance of MALA and random walk algorithms
      if(simulation==TRUE){
        par(mfrow=c(4,4),oma=c(0,0,2,0),mar=c(4,5,1,1))
      } else {
        par(mfrow=c(5,5),oma=c(0,0,2,0),mar=c(4,5,1,1))
      }
      if(print.result==TRUE){
        if(i< (burn_in1+burn_in2+2)){
          cat(paste0(" THRESHOLD INDICATOR MODEL (Logit Link):", "\n",
                     " Iteration: ",i, "\n",
                      " Accep rate mu[1] = ", round(rate.mu[1]/((i%%(adapt))+1), digits = 3),
                      " | sigma mu[1] = ", round(tun.mu[1], digits = 5), "\n",
                      " Accep rate mu[2] = ",round(rate.mu[2]/((i%%(adapt))+1), digits = 3),
                      " | sigma mu[2] = ", round(tun.mu[2], digits = 5), "\n",
                     "-------------------------------------------------------------------",  "\n",
                      sep=""))
        } else{
          cat(paste0(" THRESHOLD INDICATOR MODEL (Logit Link):", "\n",
                     " Iteration: ",i, "\n",
                     " Accep rate mu[1] = ", round(rate.mu[1]/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma mu[1] = ", round(tun.mu[1], digits = 5), "\n",
                     " Accep rate mu[2]= ",round(rate.mu[2]/(i-(burn_in1+burn_in2+2)), digits = 3),
                     " | sigma mu[2] = ", round(tun.mu[2], digits = 5), "\n",
                     "-------------------------------------------------------------------",  "\n",
                     sep=""))
        }
      }
      
      if(traceplot){
        model.param.name<-c(init.all.other.param_and.name$model.param.name.all.other.param)
        ## Plotting the traceplots
        for (ll  in 1: length(model.param.name)) {
          if(simulation==TRUE){
            plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            abline(h=true.values[ll], col=2)
          }  else {
            plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=model.param.name[ll]) # Plot for alpha.tilde
            
          }
        }
      }
      
      if((i%%adapt)-1==0){
        #plot(sigma.matrix[1:k,1],sigma.matrix[1:k,2],xlab="sigma.latent",ylab="sigma.hyper2") #plot for the scale parameters chosen adaptively
        k<-k+1
      }
      
      l<-l+1
    }
    ################################ imputations ############################################################
    if(CV=="OOS"){
      imputed_A<- rbinom(n=n1, size = 1, prob=exp(cur.samples.mu)/(1+exp(cur.samples.mu)))
      A[ind.NA]<- imputed_A[ind.NA]
    } else if(CV=="WS"){
      imputed_A<- rbinom(n=n1, size = 1, prob=exp(cur.samples.mu)/(1+exp(cur.samples.mu)))
    }
    
    if(i>burn_in1+burn_in2){
      imputed.A.sum<- imputed.A.sum+imputed_A
      imputed.A.sum.squre<- imputed.A.sum.squre + imputed_A^2
      
      post.sum.mean.mu<- post.sum.mean.mu + exp(cur.samples.mu) /(1+exp(cur.samples.mu)) ### estimated mean 
      post.sum.squre.mu<- post.sum.squre.mu + (exp(cur.samples.mu) /(1+exp(cur.samples.mu)))^2 ### estimated mean counts 
      post.sum.mean.w2<- post.sum.mean.w2 + cur.samples.w2
      post.sum.squre.w2 <- post.sum.squre.w2 + cur.samples.w2^2
    }
    
    ############ updating all the model parameters #######
    ###### Proposing new parameters for kappa_w2 (Gibbs steps) #####
    cur.samples.kappa_w2<- kappa_w2_sim_thr_ind_logit(W2=cur.samples.w2, 
                                                      node1 = node.set[,1], 
                                                      node2 = node.set[,2],
                                                      hyper_fixed=hyper_fixed)  # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_a<- kappa_a  # Proposed hyperparameters using uniform random walks
    #### Proposing new parameters for kappa_a (Gibbs steps) #####
    cur.samples.kappa_mu<- kappa_mu_sim_thr_ind_logit(mu=cur.samples.mu, intercept2 = cur.samples.intercept2, beta2 = cur.samples.beta2,
                                        W2=cur.samples.w2, Z2=Z2, hyper_fixed = hyper_fixed) # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_mu<- true.values[2]  # Proposed hyperparameters using uniform random walks
    
    #cur.samples.kappa_mu<- 1000
    
    ###### Proposing new parameters for intercept2 (Gibbs steps) ######
    cur.samples.intercept2<- intercept2_sim_thr_ind_probit(mu=cur.samples.mu, beta2=cur.samples.beta2,
                                            W2=cur.samples.w2, kappa_mu = cur.samples.kappa_mu, Z2=Z2, 
                                            hyper_fixed = hyper_fixed)
    
    ###### updating beta2 using Gibbs #########
    cur.samples.beta2<- beta2_sim_thr_ind_logit(mu=cur.samples.mu, intercept2=cur.samples.intercept2, kappa_mu=cur.samples.kappa_mu, 
                                  W2=cur.samples.w2, Z2=Z2, Z2.crossprd=Z2.crossprd, hyper_fixed=hyper_fixed)
    
    ###### updating mu using MALA #########
    prop_mu<- latent_MH_update_parallel(par_cur = cur.samples.mu,
                                        par_name = "mu",
                                        loglik_data = log_lik_indicator_model_probit,
                                        loglik_latent = log_lik_latent_mu_thr_ind_logit,
                                        var_prop = tun.mu,
                                        A = A,
                                        intercept2 = cur.samples.intercept2,
                                        W2=cur.samples.w2,
                                        beta2=cur.samples.beta2,
                                        kappa_mu = cur.samples.kappa_mu,
                                        Z2=Z2,
                                        ns=n1,
                                        nt =1
    )
    
    rate.mu<- ifelse(cur.samples.mu==prop_mu, rate.mu, rate.mu+1)
    tun.mu<-adpative_function(index_MCMC_iter=i, sigma2_adapt=tun.mu, target_accept=0.40,
                                 rate_adapt=rate.mu, burn_in1=burn_in1, burn_in2=burn_in2,
                                 adapt=adapt, adpat_param=1, adapt_seq=mu_adapt_seq2, lower.acc=0.30, upper.acc=0.50)
    
    cur.samples.mu<- prop_mu
    #cur.samples.mu<- mu
    ###### updating W2 using Gibbs #########
    mean_mu_latent_vec<- cur.samples.mu - cur.samples.intercept2 -  Z2 %*% cur.samples.beta2 
    for (ind in 1:n1) {
      cur.samples.w2[ind]<- W2_sim_Gibbs_componentwise_thr_ind_logit(node_index = ind,  
                                                       W2 = cur.samples.w2, 
                                                       mean_mu_latent = mean_mu_latent_vec[ind] ,
                                                       kappa_w2 = cur.samples.kappa_w2, 
                                                       kappa_mu = cur.samples.kappa_mu , 
                                                       nbd_info = nbd_info, 
                                                       no_of_nbd = no_of_nbd)
    }
    #cur.samples.w2<- W2
    cur.samples.w2<- cur.samples.w2- mean(cur.samples.w2) ## making sure that sum of W2 is zero
    
    ### Saving the samples after thinning the samples at every by iterations .  
    if((i%%thin)-1==0){
      samples[j,]<- c(cur.samples.kappa_w2,
                      cur.samples.kappa_mu,
                      cur.samples.intercept2,
                      cur.samples.beta2,
                      cur.samples.mu[1:2],
                      cur.samples.w2[1:2])
      
      j=j+1
    }
    ### storing the adaptive tuning parameters
    if((i%%adapt)-1==0){ # to save allexp the scale parameter of the MALA
      sigma.matrix[m,]<- c(tun.mu[1],tun.mu[2])
      m=m+1
    }
  }
  
  burnin <- burn_in1 + burn_in2
  exceed_prob.mean <- imputed.A.sum / (N.MCMC - burnin)
  exceed_prob.sd <- sqrt(imputed.A.sum.squre / (N.MCMC - burnin) - exceed_prob.mean^2)
  
  post.sum.mean.mu<- imputed.A.sum / (N.MCMC - burnin)
  post.sum.sd.mu<- sqrt(imputed.A.sum.squre / (N.MCMC - burnin) - post.sum.mean.mu^2)
  
  post.sum.mean.W2<- post.sum.mean.w2 / (N.MCMC - burnin)
  post.sum.sd.w2<- sqrt(post.sum.squre.w2 / (N.MCMC - burnin) - post.sum.mean.W2^2)
  
  
  return(
    list(
      "samples" = samples,
      "exceed_prob_est" = list( mean = exceed_prob.mean, sd = exceed_prob.sd),
      "est_random_effects" =  list(mean = post.sum.mean.W2, sd = post.sum.sd.w2),
      "est_linear_predictor" = list(mean = post.sum.mean.mu, sd=  post.sum.sd.mu)
    )
  )

}

