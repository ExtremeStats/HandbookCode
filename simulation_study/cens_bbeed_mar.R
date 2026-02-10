cens.bbeed.mar <- function (data, cov1, cov2, pm0, param0, k0, par10, par20, sig10, 
                            sig20, kn, prior.k = c("pois", "nbinom"), 
                            prior.pm = c("unif", "beta"), 
                            hyperparam, nk = 70, nsim, U = NULL, p.star, 
                            warn = FALSE) 
{
  if (nrow(data) != nrow(cov1) || nrow(data) != nrow(cov2)) {
    stop("Different number of observations/covariates")
  }
  if (length(par10) != (ncol(cov1) + 2)) {
    stop("Wrong parameter length for margin 1")
  }
  if (length(par20) != (ncol(cov2) + 2)) {
    stop("Wrong parameter length for margin 2")
  }
  if (ncol(cov1) == 1 && all(cov1 == 1)) {
    Par10 <- t(as.matrix(par10))
  }
  else {
    Par10 <- cbind(cov1 %*% par10[1:(ncol(cov1))], rep(par10[ncol(cov1) + 
                                                               1], nrow(cov1)), rep(par10[ncol(cov1) + 2], nrow(cov1)))
  }
  if (ncol(cov2) == 1 && all(cov2 == 1)) {
    Par20 <- t(as.matrix(par20))
  }
  else {
    Par20 <- cbind(cov2 %*% par20[1:(ncol(cov2))], rep(par20[ncol(cov2) + 
                                                               1], nrow(cov2)), rep(par20[ncol(cov2) + 2], nrow(cov2)))
  }
  if (is.null(U)) {
    U <- apply(data, 2, function(x) quantile(x, probs = 0.9, 
                                             type = 3))
    message("U set to 90% quantile by default on both margins \n")
  }
  data.cens <- data
  data.cens[data.cens[, 1] <= U[1], 1] <- U[1]
  data.cens[data.cens[, 2] <= U[2], 2] <- U[2]
  if (ncol(cov1) == 1 && all(cov1 == 1) && ncol(cov2) == 1 && 
      all(cov2 == 1)) {
    data.censored <- apply(data.cens, 1, function(x) ((x[1] == 
                                                         U[1]) && (x[2] == U[2])))
    censored <- which(data.censored == 1)[1]
    n.censored <- sum(data.censored)
    uncensored <- which(data.censored == 0)
    data.cens <- cbind(data.cens[c(censored, uncensored), 
    ], c(n.censored, rep(1, length(uncensored))))
    xy <- as.matrix(data[c(censored, uncensored), ])
  }
  else {
    xy <- as.matrix(data)
  }
  if (ncol(cov1) == 1 && all(cov1 == 1)) {
    data.cens.uv <- t(apply(data.cens, 1, function(val) c(ExtremalDep:::marg(val[1], 
                                                               kn = kn[1], par = Par10, der = FALSE), ExtremalDep:::marg(val[2], 
                                                                                                           kn = kn[2], par = Par20, der = FALSE))))
  }
  else {
    data.cens.uv <- t(apply(cbind(data.cens, Par10, Par20), 
                            1, function(val) c(ExtremalDep:::marg(val[1], par = val[3:5], 
                                                    kn = kn[1], der = FALSE), ExtremalDep:::marg(val[2], par = val[6:8], 
                                                                                   kn = kn[2], der = FALSE))))
  }
  r.cens.uv <- r.cens.uv_new <- rowSums(data.cens.uv)
  w.cens.uv <- w.cens.uv_new <- data.cens.uv/r.cens.uv
  data0 <- list(r.cens.uv = r.cens.uv, w.cens.uv = w.cens.uv, 
                xy = as.matrix(xy), data.cens = as.matrix(data.cens), 
                data.cens.uv = as.matrix(data.cens.uv))
  bpb <- ExtremalDep:::bp(cbind(w.cens.uv[, 2], w.cens.uv[, 1]), k0 + 1)
  bpb1 <- ExtremalDep:::bp(cbind(w.cens.uv[, 2], w.cens.uv[, 1]), k0)
  bpb2 <- ExtremalDep:::bp(cbind(w.cens.uv[, 2], w.cens.uv[, 1]), k0 - 1)
  Check <- ExtremalDep:::check.bayes(prior.k = prior.k, prior.pm = prior.pm, 
                       hyperparam = hyperparam, pm0 = pm0, param0 = param0, 
                       k0 = k0, warn = warn)
  prior.k <- Check$prior.k
  prior.pm <- Check$prior.pm
  hyperparam <- Check$hyperparam
  pm0 <- Check$pm0
  param0 <- Check$param0
  k0 <- Check$k0
  a <- Check$a
  b <- Check$b
  mu.pois <- Check$mu.pois
  pnb <- Check$pnb
  rnb <- Check$rnb
  n0 = round(5/(p.star * (1 - p.star)))
  iMax = 100
  Numbig1 <- Numbig2 <- 0
  Numsmall1 <- Numsmall2 <- 0
  n <- nrow(data)
  acc.vec.mar1 <- acc.vec.mar2 <- acc.vec <- rep(NA, nsim)
  accepted.mar1 <- accepted.mar2 <- accepted <- rep(0, nsim)
  straight.reject1 <- straight.reject2 <- rep(0, nsim)
  smar1 <- par1 <- par10
  smar2 <- par2 <- par20
  Par1 <- Par10
  Par2 <- Par20
  pm <- pm0
  param <- param0
  spm <- c(pm$p0, pm$p1)
  sk <- k <- k_new <- k0
  seta <- array(0, dim = c(nsim + 1, nk + 2))
  seta[1, 1:(k + 1)] <- param$eta
  if (k == 3) 
    q <- 0.5
  else q <- 1
  alpha <- -qnorm(p.star/2)
  d <- length(par1)
  sig1 <- sig10
  sig2 <- sig20
  sig1.vec <- sig1
  sig2.vec <- sig2
  sig1Mat <- sig2Mat <- diag(d)
  sig1.start <- sig1
  sig2.start <- sig2
  sig1.restart <- sig1
  sig2.restart <- sig2
  bpb_mat <- NULL
  pb = txtProgressBar(min = 0, max = nsim, initial = 0, style = 3)
  print.i <- seq(0, nsim, by = 100)
  for (i in 1:nsim) {
    par1_new <- as.vector(rmvnorm(1, mean = par1, sigma = sig1^2 * 
                                    sig1Mat))
    if (ncol(cov1) == 1 && all(cov1 == 1)) {
      Par1_new <- t(as.matrix(par1_new))
    }
    else {
      Par1_new <- cbind(cov1 %*% par1_new[1:(ncol(cov1))], 
                        rep(par1_new[ncol(cov1) + 1], nrow(cov1)), rep(par1_new[ncol(cov1) + 
                                                                                  2], nrow(cov1)))
    }
    if (any(U[1] < (Par1_new[, 1] - Par1_new[, 2]/Par1_new[, 
                                                           3])) || any(Par1_new[, 2] < 0) || any(Par1_new[, 
                                                                                                          3] < 0)) {
      straight.reject1[i] <- 1
      acc.vec.mar1[i] <- 0
    }
    else {
      if (ncol(cov1) == 1 && all(cov1 == 1)) {
        data.cens.uv_new <- cbind(sapply(data.cens[, 
                                                   1], function(val) ExtremalDep:::marg(val, par = Par1_new, 
                                                                          kn = kn[1], der = FALSE)), data0$data.cens.uv[, 
                                                                                                                        2])
      }
      else {
        data.cens.uv_new <- cbind(apply(cbind(data.cens[, 
                                                        1], Par1_new), 1, function(val) ExtremalDep:::marg(val[1], 
                                                                                             par = val[2:4], kn = kn[1], der = FALSE)), 
                                  data0$data.cens.uv[, 2])
      }
      r.cens.uv_new <- rowSums(data.cens.uv_new)
      w.cens.uv_new <- data.cens.uv_new/r.cens.uv_new
      data_new <- list(xy = xy, data.cens = as.matrix(data.cens), 
                       data.cens.uv = data.cens.uv_new, w.cens.uv = w.cens.uv_new, 
                       r.cens.uv = r.cens.uv_new)
      bpb_new <- ExtremalDep:::bp(cbind(w.cens.uv_new[, 2], w.cens.uv_new[, 
                                                            1]), k + 1)
      bpb1_new <- ExtremalDep:::bp(cbind(w.cens.uv_new[, 2], w.cens.uv_new[, 
                                                             1]), k)
      bpb2_new <- ExtremalDep:::bp(cbind(w.cens.uv_new[, 2], w.cens.uv_new[, 
                                                             1]), k - 1)
      ratio.mar1 <- min(exp(ExtremalDep:::cens.llik.marg(data = data_new, 
                                           coef = param$eta, bpb = bpb_new, bpb1 = bpb1_new, 
                                           bpb2 = bpb2_new, thresh = U, par1 = Par1_new, 
                                           par2 = Par2, kn = kn) - ExtremalDep:::cens.llik.marg(data = data0, 
                                                                                  coef = param$eta, bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, 
                                                                                  thresh = U, par1 = Par1, par2 = Par2, kn = kn) + 
                              log(Par1[1, 2]) - log(Par1_new[1, 2])), 1)
      acc.vec.mar1[i] <- ratio.mar1
      u.mar1 <- runif(1)
      if (u.mar1 < ratio.mar1) {
        data0 <- data_new
        bpb <- bpb_new
        bpb1 <- bpb1_new
        bpb2 <- bpb2_new
        par1 <- par1_new
        Par1 <- Par1_new
        accepted.mar1[i] <- 1
      }
    }
    smar1 <- rbind(smar1, par1)
    par2_new <- as.vector(rmvnorm(1, mean = par2, sigma = sig2^2 * 
                                    sig2Mat))
    if (ncol(cov2) == 1 && all(cov2 == 1)) {
      Par2_new <- t(as.matrix(par2_new))
    }
    else {
      Par2_new <- cbind(cov2 %*% par2_new[1:(ncol(cov2))], 
                        rep(par2_new[ncol(cov2) + 1], nrow(cov2)), rep(par2_new[ncol(cov2) + 
                                                                                  2], nrow(cov2)))
    }
    if (any(U[2] < (Par2_new[, 1] - Par2_new[, 2]/Par2_new[, 
                                                           3])) || any(Par2_new[, 2] < 0) || any(Par2_new[, 
                                                                                                          3] < 0)) {
      straight.reject2[i] <- 1
      acc.vec.mar2[i] <- 0
    }
    else {
      if (ncol(cov2) == 1 && all(cov2 == 1)) {
        data.cens.uv_new <- cbind(data0$data.cens.uv[, 
                                                     1], sapply(data.cens[, 2], function(val) ExtremalDep:::marg(val[1], 
                                                                                                   par = Par2_new, kn = kn[2], der = FALSE)))
      }
      else {
        data.cens.uv_new <- cbind(data0$data.cens.uv[, 
                                                     1], apply(cbind(data.cens[, 2], Par2_new), 
                                                               1, function(val) ExtremalDep:::marg(val[1], par = val[2:4], 
                                                                                     kn = kn[2], der = FALSE)))
      }
      r.cens.uv_new <- rowSums(data.cens.uv_new)
      w.cens.uv_new <- data.cens.uv_new/r.cens.uv_new
      data_new <- list(xy = xy, data.cens = as.matrix(data.cens), 
                       data.cens.uv = data.cens.uv_new, w.cens.uv = w.cens.uv_new, 
                       r.cens.uv = r.cens.uv_new)
      bpb_new <- ExtremalDep:::bp(cbind(w.cens.uv_new[, 2], w.cens.uv_new[, 
                                                            1]), k + 1)
      bpb1_new <- ExtremalDep:::bp(cbind(w.cens.uv_new[, 2], w.cens.uv_new[, 
                                                             1]), k)
      bpb2_new <- ExtremalDep:::bp(cbind(w.cens.uv_new[, 2], w.cens.uv_new[, 
                                                             1]), k - 1)
      ratio.mar2 <- min(exp(ExtremalDep:::cens.llik.marg(data = data_new, 
                                           coef = param$eta, bpb = bpb_new, bpb1 = bpb1_new, 
                                           bpb2 = bpb2_new, thresh = U, par1 = Par1, par2 = Par2_new, 
                                           kn = kn) - ExtremalDep:::cens.llik.marg(data = data0, coef = param$eta, 
                                                                     bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, 
                                                                     par1 = Par1, par2 = Par2, kn = kn) + log(Par2[1, 
                                                                                                                   2]) - log(Par2_new[1, 2])), 1)
      acc.vec.mar2[i] <- ratio.mar2
      u.mar2 <- runif(1)
      if (u.mar2 < ratio.mar2) {
        data0 <- data_new
        bpb <- bpb_new
        bpb1 <- bpb1_new
        bpb2 <- bpb2_new
        par2 <- par2_new
        Par2 <- Par2_new
        accepted.mar2[i] <- 1
      }
    }
    smar2 <- rbind(smar2, par2)
    k_new <- ExtremalDep:::prior_k_sampler(k)
    bpb_new <- ExtremalDep:::bp(cbind(data0$w.cens.uv[, 2], data0$w.cens.uv[, 
                                                              1]), k_new + 1)
    bpb1_new <- ExtremalDep:::bp(cbind(data0$w.cens.uv[, 2], data0$w.cens.uv[, 
                                                               1]), k_new)
    bpb2_new <- ExtremalDep:::bp(cbind(data0$w.cens.uv[, 2], data0$w.cens.uv[, 
                                                               1]), k_new - 1)
    if (prior.k == "pois") {
      p <- stats:::dpois(k_new - 3, mu.pois)/stats:::dpois(k - 3, mu.pois)
    }
    if (prior.k == "nbinom") {
      p <- stats:::dnbinom(k_new - 3, prob = pnb, size = rnb)/stats:::dnbinom(k - 
                                                                3, prob = pnb, size = rnb)
    }
    if (k_new == nk) {
      message("maximum value k reached")
      break
    }
    pm_new <- ExtremalDep:::prior_p_sampler(a = a, b = b, prior.pm = prior.pm)
    while (ExtremalDep:::check.p(pm_new, k_new)) pm_new <- ExtremalDep:::prior_p_sampler(a = a, 
                                                             b = b, prior.pm = prior.pm)
    param_new <- ExtremalDep:::rcoef(k = k_new, pm = pm_new)
    ratio <- min(exp(ExtremalDep:::cens.llik.marg(data = data0, coef = param_new$eta, 
                                    bpb = bpb_new, bpb1 = bpb1_new, bpb2 = bpb2_new, 
                                    thresh = U, par1 = Par1, par2 = Par2, kn = kn) + 
                       log(q) - ExtremalDep:::cens.llik.marg(data = data0, coef = param$eta, 
                                               bpb = bpb, bpb1 = bpb1, bpb2 = bpb2, thresh = U, 
                                               par1 = Par1, par2 = Par2, kn = kn) + log(p)), 1)
    acc.vec[i] <- ratio
    u <- runif(1)
    if (u < ratio) {
      pm <- pm_new
      param <- param_new
      bpb <- bpb_new
      bpb1 <- bpb1_new
      bpb2 <- bpb2_new
      k <- k_new
      if (k == 3) 
        q <- 0.5
      else q <- 1
      accepted[i] <- 1
    }
    sk <- c(sk, k)
    spm <- rbind(spm, c(pm$p0, pm$p1))
    seta[i + 1, 1:(k + 1)] <- param$eta
    if (i %in% print.i) 
      setTxtProgressBar(pb, i)
    if (i > 100) {
      if (i == 101) {
        sig1Mat <- cov(smar1)
        thetaM1 <- apply(smar1, 2, mean)
        sig2Mat <- cov(smar2)
        thetaM2 <- apply(smar2, 2, mean)
      }
      else {
        tmp1 <- ExtremalDep:::update.cov(sigMat = sig1Mat, i = i, 
                           thetaM = thetaM1, theta = par1, d = d)
        sig1Mat <- tmp1$sigMat
        thetaM1 <- tmp1$thetaM
        tmp2 <- ExtremalDep:::update.cov(sigMat = sig2Mat, i = i, 
                           thetaM = thetaM2, theta = par2, d = d)
        sig2Mat <- tmp2$sigMat
        thetaM2 <- tmp2$thetaM
      }
    }
    if (i > n0) {
      sig1 <- ExtremalDep:::update.sig(sig = sig1, acc = ratio.mar1, 
                         d = d, p = p.star, alpha = alpha, i = i)
      sig1.vec <- c(sig1.vec, sig1)
      if ((i <= (iMax + n0)) && (Numbig1 < 5 || Numsmall1 < 
                                 5)) {
        Toobig1 <- (sig1 > (3 * sig1.start))
        Toosmall1 <- (sig1 < (sig1.start/3))
        if (Toobig1 || Toosmall1) {
          sig1.restart <- c(sig1.restart, sig1)
          Numbig1 <- Numbig1 + Toobig1
          Numsmall1 <- Numsmall1 + Toosmall1
          i <- n0
          sig1.start <- sig1
        }
      }
      sig2 <- ExtremalDep:::update.sig(sig = sig2, acc = ratio.mar2, 
                         d = d, p = p.star, alpha = alpha, i = i)
      sig2.vec <- c(sig2.vec, sig2)
      if ((i <= (iMax + n0)) && (Numbig2 < 5 || Numsmall2 < 
                                 5)) {
        Toobig2 <- (sig2 > (3 * sig2.start))
        Toosmall2 <- (sig2 < (sig2.start/3))
        if (Toobig2 || Toosmall2) {
          sig2.restart <- c(sig2.restart, sig2)
          Numbig2 <- Numbig2 + Toobig2
          Numsmall2 <- Numsmall2 + Toosmall2
          i <- n0
          sig2.start <- sig2
        }
      }
    }
  }
  close(pb)
  return(list(method = "Bayesian", mar.fit = TRUE, type = "rawdata", 
              data = data, cov1 = cov1, cov2 = cov2, pm = spm, eta = seta, 
              k = sk, mar1 = smar1, mar2 = smar2, accepted.mar1 = accepted.mar1, 
              accepted.mar2 = accepted.mar2, accepted = accepted, 
              straight.reject1 = straight.reject1, straight.reject2 = straight.reject2, 
              acc.vec = acc.vec, acc.vec.mar1 = acc.vec.mar1, acc.vec.mar2 = acc.vec.mar2, 
              nsim = nsim, sig1.vec = sig1.vec, sig2.vec = sig2.vec, 
              threshold = U, kn = kn, prior = list(hyperparam = hyperparam, 
                                                   k = prior.k, pm = prior.pm)))
}