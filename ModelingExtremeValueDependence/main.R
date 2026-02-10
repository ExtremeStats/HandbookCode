#### This file contains all codes to generate the plots in the section "Modeling Extreme Value Dependence"
source("functions_examples.R")

scaleFUN <- function(x) sprintf("%.1f", x)

#### The exponent function, the angular density and Pickand's dependence function#####

V <- function(y1, y2, a, rho) {
  tau <- sqrt((a + 1) / (1 - rho^2))
  c <- pt(-rho * tau, df = a + 1)
  t_1 <- pt(tau * (y2 / y1 - rho), df = a + 1)
  t_2 <- pt(tau * (y1 / y2 - rho), df = a + 1)
  return(1 / (1 - c) * (1/y1 * (t_1 - c) + 1/y2 * (t_2 - c)))
}

h <- function(w, a, rho) {
  tau <- sqrt((a + 1) / (1 - rho^2))
  1 / (2 * a * w^3) * tau * ((1 - w) / w)^((1 - a) / a) *
    dt( tau * (((1 - w) / w)^(1 / a) - rho), df = a + 1) /
    (1 - pt(-rho * tau, df = a + 1))
}

A <- function(w, a, rho) {
  1 - w + 
    2 * sapply(w, function(u) 
      integrate(function(t) (u - t) * h(t, a, rho), 0, u)$value)
}


#### The summary measures \vartheta and \chi ####
a_grid <- seq(1, 10, by = 0.1)
rho_grid <- seq(0, 0.9, by = 0.01)
param_grid <- expand.grid(a = a_grid, rho = rho_grid)
param_grid$theta <- V(1, 1, param_grid$a, param_grid$rho)
param_grid$chi <- 2 - param_grid


p_theta <- ggplot(param_grid, aes(a, rho, fill = theta)) +
  geom_tile() +
  scale_fill_gradient(low="#F8766D", high="#00BFC4") +
  theme_minimal() + 
  labs(x = "a", y = expression(rho), fill = expression(vartheta)) + 
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  theme(text = element_text(size = 28)) 

#### The angular density function ####
a_grid <- c(1, 2)
rho_grid <- c(0, 0.5)
w_grid <- seq(0.01, 0.99, by = 0.001)
param_grid <- expand.grid(a = a_grid, rho = rho_grid, w = w_grid)
param_grid$h <- apply(param_grid, 1, FUN = function(z) h(z[3], z[1], z[2]))

chi_grid <- expand.grid(a = a_grid, rho = rho_grid)
chi_grid$chi <- 2 - V(1, 1, chi_grid$a, chi_grid$rho)
chi_grid$param <- interaction(chi_grid$a, chi_grid$rho)
chi_grid$param_names <- apply(chi_grid, 1, FUN = function(z) paste0("(", z[1], ", ",  z[2], ")"))
param_levels  <- levels(chi_grid$param)[order(chi_grid$chi, decreasing = T)]
param_labels  <- chi_grid$param_names[order(chi_grid$chi, decreasing = T)]
param_grid$col <- "#F8766D"
param_grid$col[(param_grid$a == 1) & (param_grid$rho == 0.5)] <- "#CD9600"
param_grid$col[(param_grid$a == 2) & (param_grid$rho == 0)] <- "#228B22"
param_grid$col[(param_grid$a == 2) & (param_grid$rho == 0.5)] <- "#00BFC4"

p_h  <- ggplot(param_grid, aes(w, h, colour = col)) +
  geom_line(linewidth = 2, alpha = 1) +
  theme_minimal() + 
  labs(x = "w", y = "h(w)") + 
  scale_colour_manual(name = expression(paste("(a, ", rho, ")")),
                      values = c("#F8766D" = "#F8766D",
                                 "#CD9600" = "#CD9600",
                                 "#228B22" = "#228B22",
                                 "#00BFC4" = "#00BFC4"),
                      limits = c("#CD9600","#F8766D", "#00BFC4", "#228B22"),
                      labels = param_labels) +  
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  theme(text = element_text(size = 28), 
        legend.position=c(.5,.85),
        legend.text.align = 0) +
  guides(colour = guide_legend(nrow = 2)) 

#### The extremal sets ####
nu_1 <- 1
rho_1 <- 0

nu_2 <- 2
rho_2 <- 0.5

n <- 1500
k <- 75
p_n <- 1 /n
m <- n / k

gamma_1 <- 1 / nu_1
b_m_1 <- uniroot(function(y) ptbt(matrix(c(y, Inf), ncol = 2), nu_1, rho_1) - 
                   (1 - 1/m), interval = c(0, m))$root
a_m_1 <- gamma_1 * b_m_1
theta_1 <- cbind(c(a_m_1, b_m_1, gamma_1), c(a_m_1, b_m_1, gamma_1))

gamma_2 <- 1 / nu_2
b_m_2 <- uniroot(function(y) ptbt(matrix(c(y, Inf), ncol = 2), nu_2, rho_2) - 
                   (1 - 1/m), interval = c(0, m))$root
a_m_2 <- gamma_2 * b_m_2
theta_2 <- cbind(c(a_m_2, b_m_2, gamma_2), c(a_m_2, b_m_2, gamma_2))

w <- seq(0.001, 0.999, by = 0.001)

Q_n_1 <- Q_hat(w, p_n, m, 
               theta = theta_1, 
               function(w) h(w, nu_1, rho_1), function(w) A(w, nu_1, rho_1), 
               type = "or")
Q_n_1_1 <- Q_n_1[, 1]
Q_n_1_2 <- Q_n_1[, 2]

Q1_true <- qtbt(p_n, nu_1, rho_1, type = "or")
Q_1_1 <- Q1_true[, 1]
Q_1_2 <- Q1_true[, 2]

df_Q1_1 <- data.frame("x1" = c(Q_n_1_1, Q_1_1) / median(Q_1_1),
                    "x2" = c(Q_n_1_2, Q_1_2) / median(Q_1_1),
                    "type" = c(rep("approx", length(Q_n_1_1)), 
                               rep("exact", length(Q_1_1))),
                    "a" = nu_1,
                    "rho" = rho_1)

Q_n_1 <- Q_hat(w, p_n, m, 
               theta = theta_2, 
               function(w) h(w, nu_2, rho_2), function(w) A(w, nu_2, rho_2), 
               type = "or")
Q_n_1_1 <- Q_n_1[, 1]
Q_n_1_2 <- Q_n_1[, 2]

Q1_true <- qtbt(p_n, nu_2, rho_2, type = "or")
Q_1_1 <- Q1_true[, 1]
Q_1_2 <- Q1_true[, 2]

df_Q1_2 <- data.frame("x1" = c(Q_n_1_1, Q_1_1) / median(Q_1_1),
                      "x2" = c(Q_n_1_2, Q_1_2) / median(Q_1_1),
                      "type" = c(rep("approx", length(Q_n_1_1)), 
                                 rep("exact", length(Q_1_1))),
                      "a" = nu_2,
                      "rho" = rho_2)

df_Q1 <- rbind(df_Q1_1, df_Q1_2)
df_Q1$param <- as.factor(apply(df_Q1, 1, 
                               FUN = function(z) paste0("(", z[4], ", ",  z[5], ")")))
df_Q1$type <- as.factor(df_Q1$type)


col <- c("exact" = "black", "approximation" = "#F8766D")
p_Q1 <- ggplot(df_Q1, aes(x1, x2, 
                          colour = param,
                          linetype = type)) +
  geom_line(linewidth = 2, alpha = 0.5) +
  theme_minimal() + 
  labs(x = element_blank(), 
       y = element_blank(),
       colour = expression(paste("(a, ", rho, ")"))) + 
  coord_cartesian(xlim = c(0.6, 1.6), ylim = c(0.6, 1.6)) +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  scale_linetype_manual(name = element_blank(),
                        values = c('approx' = 'dotted', 'exact' = 'solid'),
                        limits = c('approx', 'exact'),
                        labels = c(expression(Q[n]^(1)), expression(Q^(1)))
                        ) +
  theme(text = element_text(size = 28), 
        legend.position=c(.75,.75),
        legend.text.align = 0,
        legend.key.width = unit(2, "line"))
p_Q1

Q_n_2 <- Q_hat(w, p_n, m, 
               theta = theta_1, 
               function(w) h(w, nu_1, rho_1), function(w) A(w, nu_1, rho_1), 
               type = "and")
Q_n_2_1 <- Q_n_2[, 1]
Q_n_2_2 <- Q_n_2[, 2]

Q2_true <- qtbt(p_n, nu_1, rho_1, type = "and")
Q_2_1 <- Q2_true[, 1]
Q_2_2 <- Q2_true[, 2]

df_Q2_1 <- data.frame("x1" = c(Q_n_2_1, Q_2_1) / median(Q_2_1),
                      "x2" = c(Q_n_2_2, Q_2_2) / median(Q_2_1),
                      "type" = c(rep("approx", length(Q_n_2_1)), 
                                 rep("exact", length(Q_2_1))),
                      "a" = nu_1,
                      "rho" = rho_1)

Q_n_2 <- Q_hat(w, p_n, m, 
               theta = theta_2, 
               function(w) h(w, nu_2, rho_2), function(w) A(w, nu_2, rho_2), 
               type = "and")
Q_n_2_1 <- Q_n_2[, 1]
Q_n_2_2 <- Q_n_2[, 2]

Q2_true <- qtbt(p_n, nu_2, rho_2, type = "and")
Q_2_1 <- Q2_true[, 1]
Q_2_2 <- Q2_true[, 2]

df_Q2_2 <- data.frame("x1" = c(Q_n_2_1, Q_2_1) / median(Q_2_1),
                      "x2" = c(Q_n_2_2, Q_2_2) / median(Q_2_1),
                      "type" = c(rep("approx", length(Q_n_2_1)), 
                                 rep("exact", length(Q_2_1))),
                      "a" = nu_2,
                      "rho" = rho_2)

df_Q2 <- rbind(df_Q2_1, df_Q2_2)
df_Q2$param <- as.factor(apply(df_Q2, 1, 
                               FUN = function(z) paste0("(", z[4], ", ",  z[5], ")")))
df_Q2$type <- as.factor(df_Q2$type)


col <- c("exact" = "black", "approximation" = "#F8766D")
p_Q2 <- ggplot(df_Q2, aes(x1, x2, 
                          colour = param,
                          linetype = type)) +
  geom_line(linewidth = 2, alpha = 0.5) +
  theme_minimal() + 
  labs(x = element_blank(), 
       y = element_blank(),
       colour = expression(paste("(a, ", rho, ")"))) + 
  coord_cartesian(xlim = c(0.6, 1.6), ylim = c(0.6, 1.6)) +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  scale_linetype_manual(name = element_blank(),
                        values = c('approx' = 'dotted', 'exact' = 'solid'),
                        limits = c('approx', 'exact'),
                        labels = c(expression(Q[n]^(2)), expression(Q^(2)))
  ) +
  theme(text = element_text(size = 28), 
        legend.position=c(.75,.75),
        legend.text.align = 0,
        legend.key.width = unit(2, "line"))
p_Q2


Q_n_3 <- Q_hat(w, p_n, m, 
               theta = theta_1, 
               function(w) h(w, nu_1, rho_1), function(w) A(w, nu_1, rho_1), 
               type = "density")
Q_n_3_1 <- Q_n_3[, 1]
Q_n_3_2 <- Q_n_3[, 2]

Q3_true <- qtbt(p_n, nu_1, rho_1, type = "density")
Q_3_1 <- Q3_true[, 1]
Q_3_2 <- Q3_true[, 2]

df_Q3_1 <- data.frame("x1" = c(Q_n_3_1, Q_3_1) / median(Q_3_1),
                      "x2" = c(Q_n_3_2, Q_3_2) / median(Q_3_1),
                      "type" = c(rep("approx", length(Q_n_3_1)), 
                                 rep("exact", length(Q_3_1))),
                      "a" = nu_1,
                      "rho" = rho_1)

Q_n_3 <- Q_hat(w, p_n, m, 
               theta = theta_2, 
               function(w) h(w, nu_2, rho_2), function(w) A(w, nu_2, rho_2), 
               type = "density")
Q_n_3_1 <- Q_n_3[, 1]
Q_n_3_2 <- Q_n_3[, 2]

Q3_true <- qtbt(p_n, nu_2, rho_2, type = "density")
Q_3_1 <- Q3_true[, 1]
Q_3_2 <- Q3_true[, 2]

df_Q3_2 <- data.frame("x1" = c(Q_n_3_1, Q_3_1) / median(Q_3_1),
                      "x2" = c(Q_n_3_2, Q_3_2) / median(Q_3_1),
                      "type" = c(rep("approx", length(Q_n_3_1)), 
                                 rep("exact", length(Q_3_1))),
                      "a" = nu_2,
                      "rho" = rho_2)

df_Q3 <- rbind(df_Q3_1, df_Q3_2)
df_Q3$param <- as.factor(apply(df_Q3, 1, 
                               FUN = function(z) paste0("(", z[4], ", ",  z[5], ")")))
df_Q3$type <- as.factor(df_Q3$type)


col <- c("exact" = "black", "approximation" = "#F8766D")
p_Q3 <- ggplot(df_Q3, aes(x1, x2, 
                          colour = param,
                          linetype = type)) +
  geom_path(linewidth = 2, alpha = 0.5) +
  theme_minimal() + 
  labs(x = element_blank(), 
       y = element_blank(),
       colour = expression(paste("(a, ", rho, ")"))) + 
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) + 
  scale_linetype_manual(name = element_blank(),
                        values = c('approx' = 'dotted', 'exact' = 'solid'),
                        limits = c('approx', 'exact'),
                        labels = c(expression(Q[n]^(3)), expression(Q^(3)))
  ) +
  theme(text = element_text(size = 28), 
        legend.position=c(.75,.75),
        legend.text.align = 0,
        legend.key.width = unit(2, "line"))
p_Q3


#### Save plots ####

pdf("examples_t_distr.pdf", width=6, height=6)
par(mai=c(.7,1.2,.7,.2), mgp=c(3,0.8,0))
p_theta
p_h
p_Q1
p_Q2
p_Q3
dev.off()

