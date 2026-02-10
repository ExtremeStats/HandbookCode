library(evd)
library(mev)
library(survival)
library(EnvStats)


mrl <- function(data, tlim, nt = max(100, length(data)), conf = 0.95) {
  data <- sort(data[!is.na(data)])
  nn <- length(data)
  if (nn <= 5) stop("`data\' has too few non-missing values")
  if (missing(tlim)) {
    tlim <- c(data[1], data[nn - 4])
    tlim <- tlim - .Machine$double.eps^0.5
  }
  if (all(data <= tlim[2])) stop("upper limit for threshold is too high")
  u <- seq(tlim[1], tlim[2], length = nt)
  x <- matrix(NA, nrow = nt, ncol = 3, dimnames = list(NULL, c("lower", "mrl", "upper")))
  for (i in 1:nt) {
    data <- data[data > u[i]]
    x[i, 2] <- mean(data - u[i])
    sdev <- sqrt(var(data))
    sdev <- (qnorm((1 + conf) / 2) * sdev) / sqrt(length(data))
    x[i, 1] <- x[i, 2] - sdev
    x[i, 3] <- x[i, 2] + sdev
  }
  return(list(x = u, y = x))
}

plot_mean_excess <- function(sample, threshold = 7) {
  b <- mrl(sample)
  df <- data.frame(x = b$x, b$y)[-1,]
  interc <- df$mrl[df$x >= threshold][1]
  x <- b$x
  y <- b$y
  fig <- ggplot(data = df, aes(x=x)) +
    geom_point(aes(y = mrl), alpha = 0.5, pch = 16, cex = 0.7) +
    geom_line(aes(y=lower), linetype='dashed')+
    geom_line(aes(y=upper), linetype='dashed')+
    #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5)+
    annotate("segment",x = threshold, y = -Inf, xend = threshold, yend = interc, 
             linetype="dashed", color = "slategrey") +
    labs(x = "u", y = "Mean excess over u") +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  return(fig)
}

plot_pot <- function(sample, threshold = 7, x=NA) {
  full_sample <- as.data.frame(sample)
  if (all(is.na(x))){
    full_sample$id <- 1:nrow(full_sample)
  } else {
    full_sample$id <- x
  }
  pot <- full_sample$sample[full_sample$sample > threshold]
  indices <- full_sample$id[full_sample$sample > threshold]
  fig <- ggplot() +
    geom_point(aes(x = full_sample$id, y = full_sample$sample), color = "black", pch = 16, cex = 0.7) +
    geom_point(aes(x = indices, y = pot), color = "slategrey", pch = 16, cex = 0.7) +
    geom_hline(yintercept = threshold, color = "slategrey", linetype = "dashed") +
    labs(x = "", y = "") +
    theme(legend.position = "best") +
    theme_classic(base_size=11) +
    theme(axis.text = element_text(size = 10),
          panel.background = element_rect(fill = "white",
                                colour = "white",
                                size = 0.5, linetype = "blank"))
  return(fig)
}