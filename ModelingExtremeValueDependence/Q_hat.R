Q_hat <- function(w, p, t, theta, h = NULL, A = NULL, type = c("or", "and", "density")) {
  type <- match.arg(type)
  
  if ((type == "or") |  (type == "and")) {
    if (is.null(A) & is.null(h)) stop("At least one of A or h must be provided.")
    if (is.null(A)) A <- 
        sapply(w, FUN = function(t) integrate(
          function(u) sapply(u, FUN = function(v) max(v * t, (1 - v) * (t - v)) * h(v)), 0, t))
  }
  if (type == "or") {
    Y1_hat <- A(w) / (t * p * (1 - w))
    Y2_hat <- A(w) / (t * p * w)
  } else if (type == "and") {
    Y1_hat <- (1 - A(w)) / (t * p * (1 - w))
    Y2_hat <- (1 - A(w)) / (t * p * w)
    
    ind_w0 <- (w == 0)
    ind_w1 <- (w == 1)
    Y1_hat[ind_w1] <- 1 / (t * p)
    Y2_hat[ind_w0] <- 1 / (t * p)
  } else {
    if (is.null(h)) stop("h must be provided.")
    xi1 <- theta[3, 1]
    xi2 <- theta[3, 2]
    # Following the notation of Beranger, Padoan, Sisson (2020), we have the following:
    # h(w) = 1 / 2 * g(w, 1 - w)
    # q(w, 1 - w) = 1 / (xi1 * xi2) * 2 * w^(1 - xi1) * (1 - w)^(1 - xi2) * h(w)
    # q_*(w) = q(w, 1 - w)^(-1 / (1 + xi1 + xi2))
    # Here q_* is denoted by s.
    s <- function(w) (1 / (xi1 * xi2) * w^(1 - xi1) * (1 - w)^(1 - xi2) * 2 * h(w))^(-1 / (xi1 + xi2 + 1))
    f <- function(w) 2 * s(w) * h(w)
    v_S_hat <- integrate(f, 10^(-6), 1 - 10^(-6))$value
    
    r = 1 / s(w)
    
    # For fixed w, (r * w, r  * (1 - w)) is a boundary point of the set S
    Y1_hat <- v_S_hat * r * w / (t * p)
    Y2_hat <- v_S_hat * r * (1 - w) / (t * p)
  }
  
  
  
  X1_hat <- theta[2, 1] + theta[1, 1] * (Y1_hat^theta[3, 1] - 1) / theta[3, 1]
  X2_hat <- theta[2, 2] + theta[1, 2] * (Y2_hat^theta[3, 2] - 1) / theta[3, 2]
  
  return(cbind(X1_hat, X2_hat))
}