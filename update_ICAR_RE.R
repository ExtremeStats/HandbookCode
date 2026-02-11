#' Simulate iCAR Random Effects \code{W2} Using Gibbs Sampling for Threshold Model
#'
#' This function simulates the iCAR random effects \code{W2} for the threshold model using a component-wise Gibbs sampling approach. It updates \code{W2} based on the local neighborhood information and the latent mean.
#'
#' @param node_index An integer index indicating the specific node for which \code{W2} is being simulated.
#' @param W2 A numeric vector representing the current values of the latent variable \code{W2}.
#' @param mean_mu_latent A numeric value representing the mean of the latent variable \code{mu}.
#' @param kappa_w2 A numeric value representing the precision parameter associated with \code{W2}.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param nbd_info A list containing neighborhood information for each node, indicating which nodes are neighbors.
#' @param no_of_nbd A numeric vector representing the number of neighbors for each node.
#'
#' @return A numeric value representing the simulated value for \code{W2} at the specified \code{node_index}.
#' @export
#'
#' @examples
#' # Example data for the threshold model
#' W2 <- rnorm(100)
#' mean_mu_latent <- 0.5
#' kappa_w2 <- 1.0
#' kappa_mu <- 2.0
#' nbd_info <- list(1 = c(2, 3), 2 = c(1, 3), 3 = c(1, 2))
#' no_of_nbd <- c(2, 2, 2)
#'
#' # Simulate W2 for a specific node index
#' simulated_W2 <- W2_sim_Gibbs_componentwise_thr(1, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd)

W2_sim_Gibbs_componentwise_thr <- function(node_index,
                                           W2,
                                           mean_mu_latent,
                                           kappa_w2,
                                           kappa_mu,
                                           nbd_info,
                                           no_of_nbd) {
  mu_W2_index <- sum(W2[nbd_info[[node_index]]]) / no_of_nbd[node_index]
  var_w2 <- 1 / (no_of_nbd[node_index] * kappa_w2)
  var_mu <- 1 / kappa_mu
  
  var_sim <- (var_mu * var_w2 / (var_mu + var_w2))
  mean_sim <- var_sim * (mu_W2_index / var_w2 + mean_mu_latent / var_mu)
  sim <- rnorm(1, mean = mean_sim, sd = sqrt(var_sim))
  return(sim)
}




#' Simulate iCAR Random Effects \code{W2} Using Gibbs Sampling for Indicator Models with probit link
#'
#' This function simulates the liCAR random effects \code{W2} for the indicator model using a component-wise Gibbs sampling approach. It updates \code{W2} based on the local neighborhood information and the latent mean.
#'
#' @param node_index An integer index indicating the specific node for which \code{W2} is being simulated.
#' @param W2 A numeric vector representing the current values of the latent variable \code{W2}.
#' @param mean_mu_latent A numeric value representing the mean of the latent variable \code{mu}.
#' @param kappa_w2 A numeric value representing the precision parameter associated with \code{W2}.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param nbd_info A list containing neighborhood information for each node, indicating which nodes are neighbors.
#' @param no_of_nbd A numeric vector representing the number of neighbors for each node.
#'
#' @return A numeric value representing the simulated value for \code{W2} at the specified \code{node_index}.
#' @export
#'
#' @examples
#' # Example data for the indicator model
#' W2 <- rnorm(100)
#' mean_mu_latent <- 0.5
#' kappa_w2 <- 1.0
#' kappa_mu <- 2.0
#' nbd_info <- list(1 = c(2, 3), 2 = c(1, 3), 3 = c(1, 2))
#' no_of_nbd <- c(2, 2, 2)
#'
#' # Simulate W2 for a specific node index
#' simulated_W2 <- W2_sim_Gibbs_componentwise_thr_ind_probit(1, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd)

W2_sim_Gibbs_componentwise_thr_ind_probit <- function(node_index,
                                               W2,
                                               mean_mu_latent,
                                               kappa_w2,
                                               kappa_mu,
                                               nbd_info,
                                               no_of_nbd) {
  mu_W2_index <- sum(W2[nbd_info[[node_index]]]) / no_of_nbd[node_index]
  var_w2 <- 1 / (no_of_nbd[node_index] * kappa_w2)
  var_mu <- 1 / kappa_mu
  
  var_sim <- (var_mu * var_w2 / (var_mu + var_w2))
  mean_sim <- var_sim * (mu_W2_index / var_w2 + mean_mu_latent / var_mu)
  sim <- rnorm(1, mean = mean_sim, sd = sqrt(var_sim))
  return(sim)
}



#' Simulate iCAR Random Effects \code{W2} Using Gibbs Sampling for Indicator Models with logit link
#'
#' This function simulates the liCAR random effects \code{W2} for the indicator model using a component-wise Gibbs sampling approach. It updates \code{W2} based on the local neighborhood information and the latent mean.
#'
#' @param node_index An integer index indicating the specific node for which \code{W2} is being simulated.
#' @param W2 A numeric vector representing the current values of the latent variable \code{W2}.
#' @param mean_mu_latent A numeric value representing the mean of the latent variable \code{mu}.
#' @param kappa_w2 A numeric value representing the precision parameter associated with \code{W2}.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param nbd_info A list containing neighborhood information for each node, indicating which nodes are neighbors.
#' @param no_of_nbd A numeric vector representing the number of neighbors for each node.
#'
#' @return A numeric value representing the simulated value for \code{W2} at the specified \code{node_index}.
#' @export
#'
#' @examples
#' # Example data for the indicator model
#' W2 <- rnorm(100)
#' mean_mu_latent <- 0.5
#' kappa_w2 <- 1.0
#' kappa_mu <- 2.0
#' nbd_info <- list(1 = c(2, 3), 2 = c(1, 3), 3 = c(1, 2))
#' no_of_nbd <- c(2, 2, 2)
#'
#' # Simulate W2 for a specific node index
#' simulated_W2 <- W2_sim_Gibbs_componentwise_thr_ind_logit(1, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd)

W2_sim_Gibbs_componentwise_thr_ind_logit<- function(node_index, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd){
  mu_W2_index<- sum(W2[nbd_info[[node_index]]])/no_of_nbd[node_index]
  var_w2<- 1/(no_of_nbd[node_index] * kappa_w2)
  var_mu<- 1/kappa_mu
  
  var_sim<- (var_mu * var_w2 /(var_mu + var_w2))
  mean_sim<- var_sim * (mu_W2_index / var_w2 + mean_mu_latent/var_mu)
  sim<- rnorm(1, mean= mean_sim, sd= sqrt(var_sim))
  return(sim)
}



#' Simulate iCAR Random Effects \code{W1} Using Gibbs Sampling for Joint Models
#'
#' This function simulates the iCAR random effects \code{W1} for joint models using a component-wise Gibbs sampling approach. It updates \code{W1} based on the local neighborhood information and the latent means.
#'
#' @param node_index An integer index indicating the specific node for which \code{W1} is being simulated.
#' @param W1 A numeric vector representing the current values of the latent variable \code{W1}.
#' @param mean_eta_latent A numeric value representing the mean of the latent variable \code{eta}.
#' @param mean_mu_latent A numeric value representing the mean of the latent variable \code{mu}.
#' @param kappa_w1 A numeric value representing the precision parameter associated with \code{W1}.
#' @param kappa_eta A numeric value representing the precision parameter associated with \code{eta}.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param beta A numeric value representing the effect size of the latent variable \code{mu} on \code{W1}.
#' @param nbd_info A list containing neighborhood information for each node, indicating which nodes are neighbors.
#' @param no_of_nbd A numeric vector representing the number of neighbors for each node.
#'
#' @return A numeric value representing the simulated value for \code{W1} at the specified \code{node_index}.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' W1 <- rnorm(100)
#' mean_eta_latent <- 0.5
#' mean_mu_latent <- 1.0
#' kappa_w1 <- 1.0
#' kappa_eta <- 1.0
#' kappa_mu <- 2.0
#' beta <- 0.8
#' nbd_info <- list(1 = c(2, 3), 2 = c(1, 3), 3 = c(1, 2))
#' no_of_nbd <- c(2, 2, 2)
#'
#' # Simulate W1 for a specific node index
#' simulated_W1 <- W1_sim_Gibbs_componentwise_JM(1, W1, mean_eta_latent, mean_mu_latent, kappa_w1, kappa_eta, kappa_mu, beta, nbd_info, no_of_nbd)

W1_sim_Gibbs_componentwise_JM <- function(node_index,
                                          W1,
                                          mean_eta_latent,
                                          mean_mu_latent,
                                          kappa_w1,
                                          kappa_eta,
                                          kappa_mu,
                                          beta,
                                          nbd_info,
                                          no_of_nbd) {
  #browser()
  mu_W1_index <- sum(W1[nbd_info[[node_index]]]) / no_of_nbd[node_index]
  var_w1 <- 1 / (no_of_nbd[node_index] * kappa_w1)
  var_eta <- 1 / kappa_eta
  var_mu <- 1 / kappa_mu
  
  var_sim <- (1 / (beta ^ 2 / var_mu +  1 / var_eta + 1 / var_w1))
  mean_sim <- var_sim * (mu_W1_index / var_w1 +  mean_eta_latent / var_eta + beta * mean_mu_latent /
                           var_mu)
  sim <- rnorm(1, mean = mean_sim, sd = sqrt(var_sim))
  return(sim)
}



#' Simulate iCAR Random Effects \code{W2} Using Gibbs Sampling for Joint Models
#'
#' This function simulates the iCAR random effects \code{W2} for joint models using a component-wise Gibbs sampling approach. It updates \code{W2} based on the local neighborhood information and the latent mean.
#'
#' @param node_index An integer index indicating the specific node for which \code{W2} is being simulated.
#' @param W2 A numeric vector representing the current values of the latent variable \code{W2}.
#' @param mean_mu_latent A numeric value representing the mean of the latent variable \code{mu}.
#' @param kappa_w2 A numeric value representing the precision parameter associated with \code{W2}.
#' @param kappa_mu A numeric value representing the precision parameter associated with \code{mu}.
#' @param nbd_info A list containing neighborhood information for each node, indicating which nodes are neighbors.
#' @param no_of_nbd A numeric vector representing the number of neighbors for each node.
#'
#' @return A numeric value representing the simulated value for \code{W2} at the specified \code{node_index}.
#' @export
#'
#' @examples
#' # Example data for the joint model
#' W2 <- rnorm(100)
#' mean_mu_latent <- 0.5
#' kappa_w2 <- 1.0
#' kappa_mu <- 1.0
#' nbd_info <- list(1 = c(2, 3), 2 = c(1, 3), 3 = c(1, 2))
#' no_of_nbd <- c(2, 2, 2)
#'
#' # Simulate W2 for a specific node index
#' simulated_W2 <- W2_sim_Gibbs_componentwise_JM(1, W2, mean_mu_latent, kappa_w2, kappa_mu, nbd_info, no_of_nbd)

W2_sim_Gibbs_componentwise_JM <- function(node_index,
                                          W2,
                                          mean_mu_latent,
                                          kappa_w2,
                                          kappa_mu,
                                          nbd_info,
                                          no_of_nbd) {
  mu_W2_index <- sum(W2[nbd_info[[node_index]]]) / no_of_nbd[node_index]
  var_w2 <- 1 / (no_of_nbd[node_index] * kappa_w2)
  var_mu <- 1 / kappa_mu
  
  var_sim <- (var_mu * var_w2 / (var_mu + var_w2))
  mean_sim <- var_sim * (mu_W2_index / var_w2 + mean_mu_latent / var_mu)
  sim <- rnorm(1, mean = mean_sim, sd = sqrt(var_sim))
  return(sim)
}
