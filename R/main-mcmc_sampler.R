#' MCMC Sampler for Joint Spatial Modeling of Landslide Counts and Sizes
#'
#' @description
#' Implements an MCMC sampler for the joint modeling of landslide count data and landslide size data, as described in the book chapter by Yadav et al. (2024), "Statistics of Extremes for Natural Hazards: Landslides and Earthquakes." The function allows for various covariate inputs, cross-validation types, and distributional assumptions, offering flexibility for both real data analysis and simulation experiments.
#'
#' @param Y A vector of length `n` containing landslide count data.
#' @param A A vector of length `n` containing landslide size data corresponding to the landslide counts.
#' @param Z1 A design matrix of dimension `n x p` consisting of covariates used in the log-linear predictor of counts.
#' @param Z2 A design matrix of dimension `n x q` consisting of covariates used in the log-linear predictor of sizes.
#' @param CV A character string specifying the type of cross-validation. Two choices:
#' \itemize{
#'   \item `"OOS"`: Out-of-sample cross-validation.
#'   \item `"WS"`: Within-sample assessment.
#' }
#' @param mark_dist A character string specifying the distribution type for landslide sizes. Three choices:
#' \itemize{
#'   \item `"eGPD"`: Extended Generalized Pareto Distribution.
#'   \item `"bGPD"`: Beta-GPD mixtures.
#'   \item `"tgGPD"`: Truncated gamma-GPD mixtures.
#' }
#' @param thr.family A character string specifying the distribution used to estimate the non-stationary thresholds of landslide sizes. Only used when `mark_dist` is `"bGPD"` or `"tgGPD"`. Two choices:
#' \itemize{
#'   \item `"gamma"`: Gamma distribution.
#'   \item `"logNormal"`: Log-normal distribution.
#' }
#' @param q.probs A numeric vector of quantiles at which to calculate the landslide risk. Default is `quantile(sqrt(A), probs = seq(0.50, 0.99, 0.05))`.
#' @param q.probs.thr A numeric vector specifying the threshold quantiles to extract threshold exceedances. Default is `0.83`.
#' @param no.rm.obs The number of observations randomly removed when performing the out-of-sample (OOS) experiment.
#' @param N.MCMC An integer specifying the number of MCMC samples. Default is `2e4`.
#' @param tun.hyper.mu A numeric value for the tuning parameter of the normal Metropolis-Hastings steps when updating the hyperparameters for `mark_dist`.
#' @param tun.hyper.GP A numeric value for the tuning parameter of the normal Metropolis-Hastings steps when updating the hyperparameters for the generalized Pareto distribution.
#' @param tun.eta A numeric value for the tuning parameter of the normal Metropolis-Hastings steps when updating the log-linear predictor `eta`.
#' @param tun.mu A numeric value for the tuning parameter of the normal Metropolis-Hastings steps when updating the log-linear predictor `mu`.
#' @param thin An integer specifying the thinning interval used in the MCMC. Default is `10`.
#' @param adapt An integer specifying the number of samples after which to perform adaptive MCMC. Default is `100`.
#' @param burn_in1 The first burn-in period of MCMC samples. If not specified, it defaults to `N.MCMC / 4`.
#' @param burn_in2 The second burn-in period of MCMC samples. If not specified, it defaults to `N.MCMC / 2`.
#' @param true.values A named list or vector of the true parameter values, used only for simulation experiments.
#' @param simulation A logical value. If `TRUE`, a simulation experiment is conducted; otherwise, real data experiments are performed. Default is `FALSE`.
#' @param print.result A logical value. If `TRUE`, results are printed at fixed numbers of iterations. Default is `TRUE`.
#' @param traceplot A logical value. If `TRUE`, visual traceplots are shown as the MCMC progresses. Default is `FALSE`.
#' @param init.seed initial seed (leave it to default NULL Set it NULL all the time)
#' @param model_type A character string specifying the type of model. Two choices:
#' \itemize{
#'   \item `"FE"`: Fixed effects model.
#'   \item `"jSp"`: Joint spatial model.
#' }
#' @param adjacensy A matrix or list specifying the adjacency structure of the slope units for the spatial model.
#' @param model.base A logical value. If `TRUE`, sets `beta = 0`, indicating that counts and sizes are modeled as independent of each other. Default is `FALSE`.
#' @param samples.store An integer specifying the number of samples stored for calculating summary statistics (e.g., QQ-plots, uncertainty quantification).
#' @param fit_thr_model_only A logical value. If `TRUE` fit only threshold and threshold indicator models only.
#' @param ind_model_link_fun A character string specifying the type of link function used in indictaor model. Two choices:
#' \itemize{
#'   \item `"probit"`: for probit regression.
#'   \item `"logit"`: for logistic regression. 
#' }
#'
#' @return A list containing results of the MCMC sampling, including estimated parameters, risk quantiles, and any diagnostics.
#' @export
#'
#' @references
#' 1. Yadav, R., Lombardo, L., & Huser, R. (2024). *Statistics of Extremes for Natural Hazards: Landslides and Earthquakes*. arXiv preprint. [arXiv:2404.09156](https://arxiv.org/abs/2404.09156).
#' 2. Yadav, R., Huser, R., Opitz, T., & Lombardo, L. (2024). *Joint modeling of landslide counts and sizes using spatial marked point processes with sub-asymptotic mark distributions*. Journal of the Royal Statistical Society Series C (JRSSC): Applied Statistics, qlad077. [https://doi.org/10.1093/jrsssc/qlad077](https://doi.org/10.1093/jrsssc/qlad077)
#' @examples
#' \dontrun{
#' # Example usage
#' mcmc_sampler(Y = count_data, A = size_data, Z1 = covariates_counts, Z2 = covariates_sizes,
#'              CV = "OOS", mark_dist = "eGPD", thr.family = "gamma", model_type = "jSp",
#'              adjacensy = adj_matrix, q.probs = seq(0.50, 0.99, by = 0.05),
#'              N.MCMC = 20000, simulation = TRUE, init.seed = 123)
#' }


mcmc_sampler <- function(Y,
                         A,
                         Z1,
                         Z2,
                         CV,
                         mark_dist,
                         thr.family = "gamma",
                         model_type,
                         adjacensy,
                         ind_model_link_fun = "logit",
                         q.probs = as.numeric(quantile(sqrt(A), probs = seq(0.50, 0.99, 0.05))),
                         q.probs.thr = 0.83,
                         no.rm.obs = 2000,
                         N.MCMC = 2e4,
                         model.base =  FALSE,
                         tun.hyper.mu = 0.01,
                         tun.hyper.GP = 0.01,
                         tun.eta = 1,
                         tun.mu = 1,
                         thin = 10,
                         adapt = 100,
                         burn_in1 = NULL,
                         burn_in2 = NULL,
                         true.values = NULL,
                         simulation = FALSE,
                         print.result = TRUE,
                         traceplot = FALSE,
                         samples.store,
                         fit_thr_model_only = FALSE,
                         init.seed = NULL) {
  if (is.null(burn_in1) & is.null(burn_in2)) {
    burn_in1 <- floor(N.MCMC / 2)
    burn_in2 <- floor(N.MCMC / 4)
  }
  
  if (samples.store > (N.MCMC - burn_in1 - burn_in2)) {
    stop(
      paste0(
        "Error: `samples.store` (",
        samples.store,
        ") should be less than the number of post burn-in samples `N.MCMC - burn_in1 - burn_in2` (",
        N.MCMC - burn_in1 - burn_in2,
        ")."
      )
    )
  }
  
  if (no.rm.obs > length(A)) {
    stop(
      paste0(
        "Error: `no.rm.obs` (",
        no.rm.obs,
        ") should be less than the number of spatial locations `length(A)` (",
        length(A),
        ")."
      )
    )
  }
  
  
  nbd_info = adjacensy$nbs
  no_of_nbd = adjacensy$nnbs
  m.bar <- mean(no_of_nbd)
  
  N <- adjacensy$n
  node.set <- c()
  for (i in 1:length(adjacensy$nbs)) {
    for (j in 1:length(adjacensy$nbs[[i]])) {
      if (i < adjacensy$nbs[[i]][j]) {
        node.set <- rbind(node.set, c(i, adjacensy$nbs[[i]][j]))
      }
    }
  }
  
  
  eta_adapt_seq2 <- seq(from = adapt, to = N.MCMC, by = adapt)
  mu_adapt_seq2 <- seq(from = adapt, to = N.MCMC, by = adapt)
  hyper.mu_adapt_seq2 <- seq(from = adapt, to = N.MCMC, by = adapt)
  
  
  if (CV == "WS") {
    ind_zeros_counts <- A == 0
    A_data_frame_size <- data.frame(original_A = A, A_with_NA = A)
    ind_miss <- NULL
  } else if (CV == "OOS") {
    A_data_frame_size <- data.frame(original_A = A, A_with_NA = A)
    set.seed(1)
    ind_miss <- sample(1:length(A), size = no.rm.obs, replace = FALSE)
    A_data_frame_size[ind_miss, 2] <- NA
    ind_zeros_counts <- A_data_frame_size$A_with_NA == 0
    ind_zeros_counts[ind_miss] <- FALSE
  }
  
  
  hyper.mu_fixed <- c(0.25, 0.25)  ## hyperparameters in the hyper.mu parameters
  hyper_fixed <- list(
    "kappa_eta" = c(0.75, 3),
    # kappa_eta ~ Gamma(shape=hyper_fixed[1], rate=hyper_fixed[2])
    "kappa_w1" = c(0.75, 3 / (m.bar * 0.7 ^ 2)),
    # kappa_w1 ~ Gamma(shape=hyper_fixed[3], rate=hyper_fixed[4])
    "kappa_w2"  = c(0.75, 3 / (m.bar * 0.7 ^ 2)),
    # kappa_w2 ~ Gamma(shape=hyper_fixed[5], rate=hyper_fixed[6])
    "kappa_mu" = c(0.75, 3),
    # kappa_mu ~ Gamma(shape=hyper_fixed[7], rate=hyper_fixed[8])
    "beta" = 0.01,
    # beta ~ norm(mean=0, precision=hyper_fixed[13]))
    "beta1" = 0.01,
    # beta1 ~ mvtnorm_p(mean=0, precision=hyper_fixed[14] I_p))
    "beta2" = 0.01,
    # beta2 ~ mvtnorm_q(mean=0, precision=hyper_fixed[15] I_q))
    "intercept1" = 0.01,
    # intercept1~ rnorm(1,mean=0, precision=hyper_fixed[16])
    "intercept2" = 0.01
  )     # intercept2~ rnorm(1,mean=0, precision=hyper_fixed[17])
  if (mark_dist == "bGPD" | mark_dist == "tgGPD") {
    threshold_model_output <- mcmc_sampler_threhshold_model(
      N.MCMC = N.MCMC,
      A = A_data_frame_size$A_with_NA,
      ind.NA = is.na(A_data_frame_size$A_with_NA),
      Z2 = as.matrix(Z2),
      thin = thin,
      adapt = adapt,
      burn_in1 = burn_in1,
      burn_in2 = burn_in2,
      tun.mu = tun.mu,
      tun.hyper.mu = tun.hyper.mu,
      hyper_fixed = hyper_fixed,
      print.result = print.result,
      traceplot = traceplot,
      CV = CV,
      true.values = true.values,
      simulation = simulation,
      nbd_info = nbd_info,
      no_of_nbd = no_of_nbd,
      node.set = node.set,
      thr.family = thr.family,
      ind_zeros_counts = ind_zeros_counts,
      q.probs = q.probs.thr,
      hyper.mu_adapt_seq2 = hyper.mu_adapt_seq2,
      mu_adapt_seq2 = mu_adapt_seq2,
      eta_adapt_seq2 = eta_adapt_seq2,
      ind_miss = ind_miss,
      init.seed = init.seed
    )
    if(ind_model_link_fun == "probit"){
    thresholds_indcator_model_output <-  mcmc_sampler_indicator_model_probit(
      N.MCMC = N.MCMC,
      A = threshold_model_output$est_thr_indicator,
      ind.NA = is.na(
        threshold_model_output$est_thr_indicator
      ),
      CV = CV,
      Z2 = as.matrix(Z2),
      thin = thin,
      adapt = adapt,
      burn_in1 = burn_in1,
      burn_in2 = burn_in2,
      ind_zero = threshold_model_output$est_thr_indicator==0,
      hyper_fixed = hyper_fixed,
      print.result = print.result,
      traceplot = traceplot,
      true.values = true.values,
      simulation = simulation,
      nbd_info = nbd_info,
      no_of_nbd = no_of_nbd,
      node.set = node.set,
      hyper.mu_adapt_seq2 =
        hyper.mu_adapt_seq2,
      mu_adapt_seq2 =
        mu_adapt_seq2,
      eta_adapt_seq2 = eta_adapt_seq2,
      init.seed = init.seed
    )
    } else {
      thresholds_indcator_model_output <-  mcmc_sampler_indicator_model_logit(
        N.MCMC = N.MCMC,
        A = threshold_model_output$est_thr_indicator,
        ind.NA = is.na(
          threshold_model_output$est_thr_indicator
        ),
        CV = CV,
        Z2 = as.matrix(Z2),
        thin = thin,
        adapt = adapt,
        burn_in1 = burn_in1,
        burn_in2 = burn_in2,
        tun.mu =  tun.mu,
        ind_zero = threshold_model_output$est_thr_indicator==0,
        hyper_fixed = hyper_fixed,
        print.result = print.result,
        traceplot = traceplot,
        true.values = true.values,
        simulation = simulation,
        nbd_info = nbd_info,
        no_of_nbd = no_of_nbd,
        node.set = node.set,
        hyper.mu_adapt_seq2 =
          hyper.mu_adapt_seq2,
        mu_adapt_seq2 =
          mu_adapt_seq2,
        eta_adapt_seq2 = eta_adapt_seq2,
        init.seed = init.seed
      )
    }
    
    if (CV == "WS") {
      ######### Divide the datasets for the OSD setting
      Y_data_frame_count <- data.frame(original_Y = Y, Y_with_NA = Y)
      threshold <-  threshold_model_output$threshold
      thr.acces.ind <- threshold_model_output$est_thr_indicator  # A_data_frame_size$original_A > threshold
      
    } else if (CV == "OOS") {
      ### missingness will be fixed once the thresholds and the threshold probability is fixed
      Y_data_frame_count <- data.frame(original_Y = Y, Y_with_NA = Y)
      Y_data_frame_count[ind_miss, 2] <- NA
      
      threshold <- threshold_model_output$threshold
      thr.acces.ind<- threshold_model_output$est_thr_indicator
      # thr.acces.ind <- A_data_frame_size$A_with_NA > threshold
      # thr.acces.ind[ind_miss] <- FALSE
    }
    
  } else{
    exceed_prob_miss <- NULL
    threshold <- NULL
    thr.acces.ind <- NULL
    thresholds_indcator_model_output <- NULL
    threshold_model_output <- NULL
    
    
    if (CV == "WS") {
      ######### Divide the datasets for the OSD setting
      Y_data_frame_count <- data.frame(original_Y = Y, Y_with_NA = Y)
    } else if (CV == "OOS") {
      ### missingens will be fixed once the thresholds and the threshold probability is fixed
      Y_data_frame_count <- data.frame(original_Y = Y, Y_with_NA = Y)
      Y_data_frame_count[ind_miss, 2] <- NA
    }
  }
  if (fit_thr_model_only) {
    results <- list(
      "thr.info" =  threshold_model_output,
      "thr.indicator.info" = thresholds_indcator_model_output,
      "JM.info" =  NULL,
      "mark_dist"  = mark_dist,
      "thr.family" = thr.family,
      "model_type" =  model_type,
      "data" = list(
        "count.full" = Y_data_frame_count$Y_with_NA,
        "count.with.NA" = Y_data_frame_count$original_Y,
        "size.full" =  A_data_frame_size$A_with_NA,
        "size.with.NA" =  A_data_frame_size$original_A
      ),
      "CV" = CV
    )
  } else {
    if (model_type == "FE") {
      results_JM <- mcmc_sampler_joint_model_FE(
        N.MCMC = N.MCMC,
        Y = Y_data_frame_count$Y_with_NA,
        ind_NA_Y = is.na(Y_data_frame_count$Y_with_NA),
        A = A_data_frame_size$A_with_NA,
        ind_NA_A = is.na(A_data_frame_size$A_with_NA),
        Z1 = as.matrix(Z1),
        Z2 = as.matrix(Z2),
        model_type = model_type,
        thin = thin,
        adapt = adapt,
        burn_in1 = burn_in1,
        burn_in2 = burn_in2,
        tun.eta = tun.eta,
        tun.mu = tun.mu,
        tun.hyper.mu = tun.hyper.mu,
        tun.hyper.GP = tun.hyper.GP,
        mark_dist = mark_dist,
        hyper_fixed = hyper_fixed,
        print.result = print.result,
        traceplot = traceplot,
        model.base = model.base,
        CV = CV,
        true.values = true.values,
        simulation = simulation,
        nbd_info = nbd_info,
        no_of_nbd = no_of_nbd,
        ind_zeros_counts = ind_zeros_counts,
        threshold = threshold,
        thr.acces.ind = thr.acces.ind,
        thr.prob = thr.prob,
        q.probs =  q.probs,
        hyper.mu_adapt_seq2 = hyper.mu_adapt_seq2,
        mu_adapt_seq2 = mu_adapt_seq2,
        eta_adapt_seq2 = eta_adapt_seq2,
        samples.store = samples.store,
        init.seed = init.seed
      )
      
    } else{
      results_JM <- mcmc_sampler_joint_model_jSp(
        N.MCMC = N.MCMC,
        Y = Y_data_frame_count$Y_with_NA,
        ind_NA_Y = is.na(Y_data_frame_count$Y_with_NA),
        A = A_data_frame_size$A_with_NA,
        ind_NA_A = is.na(A_data_frame_size$A_with_NA),
        Z1 = as.matrix(Z1),
        Z2 = as.matrix(Z2),
        model_type = model_type,
        Sim_data = Sim_mark_data,
        thin = thin,
        adapt = adapt,
        burn_in1 = burn_in1,
        burn_in2 = burn_in2,
        tun.eta = tun.eta,
        tun.mu = tun.mu,
        tun.hyper.mu = tun.hyper.mu,
        tun.hyper.GP = tun.hyper.GP,
        mark_dist = mark_dist,
        hyper_fixed = hyper_fixed,
        print.result = print.result,
        traceplot = traceplot,
        model.base = model.base,
        CV = CV,
        true.values = true.values,
        simulation = simulation,
        nbd_info = nbd_info,
        no_of_nbd = no_of_nbd,
        node.set = node.set,
        ind_zeros_counts = ind_zeros_counts,
        threshold = threshold,
        thr.acces.ind = thr.acces.ind,
        thr.prob = thr.prob,
        q.probs =  q.probs,
        hyper.mu_adapt_seq2 = hyper.mu_adapt_seq2,
        mu_adapt_seq2 = mu_adapt_seq2,
        eta_adapt_seq2 = eta_adapt_seq2,
        samples.store = samples.store,
        init.seed = init.seed
      )
    }
    #### save all the required results
    results <- list(
      "thr.info" =  threshold_model_output,
      "thr.indicator.info" = thresholds_indcator_model_output,
      "JM.info" =  results_JM,
      "mark_dist"  = mark_dist,
      "thr.family" = thr.family,
      "model_type" =  model_type,
      "data" = list(
        "count.full" = Y_data_frame_count$original_Y,
        "count.with.NA" = Y_data_frame_count$Y_with_NA,
        "size.full" =  A_data_frame_size$original_A,
        "size.with.NA" =  A_data_frame_size$A_with_NA
      ),
      "CV" = CV
    )
  }
  
  return(results)
}
