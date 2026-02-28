library(rjags)
library(coda)

model_string <- "
model{
  # Priors
  Se1 ~ dbeta(1,1)
  Sp1 ~ dbeta(1,1)
  Se2 ~ dbeta(1,1)
  Sp2 ~ dbeta(1,1)
  pi1 ~ dbeta(1,1)
  pi2 ~ dbeta(1,1)

  # Population 1 cell probabilities: (Y1,Y2) = (1,1),(1,0),(0,1),(0,0)
  p1[1] <- pi1*Se1*Se2 + (1-pi1)*(1-Sp1)*(1-Sp2)
  p1[2] <- pi1*Se1*(1-Se2) + (1-pi1)*(1-Sp1)*Sp2
  p1[3] <- pi1*(1-Se1)*Se2 + (1-pi1)*Sp1*(1-Sp2)
  p1[4] <- pi1*(1-Se1)*(1-Se2) + (1-pi1)*Sp1*Sp2

  # Population 2 (same Se/Sp, different prevalence)
  p2[1] <- pi2*Se1*Se2 + (1-pi2)*(1-Sp1)*(1-Sp2)
  p2[2] <- pi2*Se1*(1-Se2) + (1-pi2)*(1-Sp1)*Sp2
  p2[3] <- pi2*(1-Se1)*Se2 + (1-pi2)*Sp1*(1-Sp2)
  p2[4] <- pi2*(1-Se1)*(1-Se2) + (1-pi2)*Sp1*Sp2

  # Likelihood (2×2 tables per population)
  n1[1:4] ~ dmulti(p1[1:4], N1)
  n2[1:4] ~ dmulti(p2[1:4], N2)
}
"

cell_probs <- function(pi, Se1, Sp1, Se2, Sp2) {
  c(
    pi*Se1*Se2 + (1-pi)*(1-Sp1)*(1-Sp2),             # (1,1)
    pi*Se1*(1-Se2) + (1-pi)*(1-Sp1)*Sp2,             # (1,0)
    pi*(1-Se1)*Se2 + (1-pi)*Sp1*(1-Sp2),             # (0,1)
    pi*(1-Se1)*(1-Se2) + (1-pi)*Sp1*Sp2              # (0,0)
  )
}

sim_hw <- function(N1, N2, Se1, Sp11, Sp12, Se2, Sp2, pi1, pi2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p1 <- cell_probs(pi1, Se1, Sp11, Se2, Sp2)
  p2 <- cell_probs(pi2, Se1, Sp12, Se2, Sp2)
  n1 <- as.vector(rmultinom(1, N1, p1))
  n2 <- as.vector(rmultinom(1, N2, p2))
  list(n1 = n1, n2 = n2, p1 = p1, p2 = p2)
}

# True values
truth <- list(Se11 = 0.5, Se1 = 0.7, Se2 = 0.75, Sp11 = 0.7, Sp12 = 0.9, Sp2 = 0.95, pi1 = 0.1, pi2 = 0.9)
N1 <- 5000; N2 <- 5000

# One Simulated dataset
data <- sim_hw(N1 = N1, N2 = N2, 
               Se1 = truth$Se1, 
               Sp11 = truth$Sp11, Sp12 = truth$Sp12, 
               Se2 = truth$Se2, Sp2 = truth$Sp2, 
               pi1 = truth$pi1, pi2 = truth$pi2, 
               seed = NULL)

# JAGS data list (order (11,10,01,00)) 
jdata <- list(
  n1 = data$n1, n2 = data$n2,
  N1 = N1, N2 = N2
)

# Initial values 
inits_fun <- function() list(
  Se1 = runif(1, 0.5, 0.999),
  Sp1 = runif(1, 0.5, 0.999),
  Se2 = runif(1, 0.5, 0.999),
  Sp2 = runif(1, 0.5, 0.999),
  pi1 = runif(1, 0.05, 0.95),
  pi2 = runif(1, 0.05, 0.95)
)


jm <- jags.model(textConnection(model_string),
                 data = jdata, inits = inits_fun, n.chains = 1, n.adapt = 2000)

# Burn-in
update(jm, 5000) 

samples <- coda.samples(jm,
                        variable.names = c("Se1","Sp1","Se2","Sp2", "pi1", "pi2"),
                        n.iter = 30000)

s_summary <- summary(samples)
post_medians <- s_summary$quantiles[,"50%"] #Alternative median(samples[[1]][,1]) etc
post_q025  <- s_summary$quantiles[,"2.5%"]
post_q975  <- s_summary$quantiles[,"97.5%"]

post_summary <- data.frame(
  Median = post_medians,
  Q2.5 = post_q025,
  Q97.5 = post_q975
)

round(post_summary, 3)

#post_all <- as.matrix(samples)

#par(mfrow = c(2,3))
#plot(density(post_all[,"Se1"]),
#     main = "Posterior density of Se1",
#     xlab = "Se1", lwd = 2, col = "blue",
#     xlim = c(0,1))
#abline(v = truth$Se1, col = "red", lwd = 2, lty = 2)

#plot(density(post_all[,"Sp1"]),
#     main = "Posterior density of Sp1",
#     xlab = "Sp1", lwd = 2, col = "blue",
#     xlim = c(0,1))
#abline(v = truth$Sp1, col = "red", lwd = 2, lty = 2)

#plot(density(post_all[,"Se2"]),
#     main = "Posterior density of Se2",
#     xlab = "Se2", lwd = 2, col = "blue",
#     xlim = c(0,1))
#abline(v = truth$Se2, col = "red", lwd = 2, lty = 2)

#plot(density(post_all[,"Sp2"]),
#     main = "Posterior density of Sp2",
#     xlab = "Sp2", lwd = 2, col = "blue",
#     xlim = c(0,1))
#abline(v = truth$Sp2, col = "red", lwd = 2, lty = 2)

#plot(density(post_all[,"pi1"]),
#     main = "Posterior density of Pi1",
#     xlab = "Pi1", lwd = 2, col = "blue",
#     xlim = c(0,1))
#abline(v = truth$pi1, col = "red", lwd = 2, lty = 2)

#plot(density(post_all[,"pi2"]),
#     main = "Posterior density of Pi2",
#     xlab = "Pi2", lwd = 2, col = "blue",
#     xlim = c(0,1))
#abline(v = truth$pi2, col = "red", lwd = 2, lty = 2)