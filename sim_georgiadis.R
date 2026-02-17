cell_probs_georgiadis <- function(pi, Se1, Se21, Se20, Sp1, Sp21, Sp20) {
  c(
    pi*Se1*Se21 + (1-pi)*(1-Sp1)*(1-Sp21),             # (1,1)
    pi*Se1*(1-Se21) + (1-pi)*(1-Sp1)*Sp21,             # (1,0)
    pi*(1-Se1)*Se20 + (1-pi)*Sp1*(1-Sp20),             # (0,1)
    pi*(1-Se1)*(1-Se20) + (1-pi)*Sp1*Sp20              # (0,0)
  )
}

sim_georgiadis <- function(N1, N2, Se1, Se21, Se20, Sp1, Sp21, Sp20, pi1, pi2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p1 <- cell_probs_georgiadis(pi1, Se1, Se21, Se20, Sp1, Sp21, Sp20)
  p2 <- cell_probs_georgiadis(pi2, Se1, Se21, Se20, Sp1, Sp21, Sp20)
  n1 <- as.vector(rmultinom(1, N1, p1))
  n2 <- as.vector(rmultinom(1, N2, p2))
  list(n1 = n1, n2 = n2, p1 = p1, p2 = p2)
}

# True values
rhoD <- 0.5
rhonD <- 0.6
Se1 <- 0.85
Se2 <- 0.95
Sp1 <- 0.8
Sp2 <- 0.9

Se11 <- rhoD*sqrt(Se1*(1-Se1)*Se2*(1-Se2)) + Se1*Se2
Se21 <- Se11/Se1
Se01 <- Se2 - Se11
Se20 <- Se01/(1-Se1)

Sp00 <- rhonD*sqrt(Sp1*(1-Sp1)*Sp2*(1-Sp2)) + Sp1*Sp2
Sp20 <- Sp00/Sp1
Sp10 <- Sp2 - Sp00
Sp21 <- Sp10/(1-Sp1)

truth <- list(Se1 = Se1, Se21 = Se21, Se20 = Se20, 
              Sp1 = Sp1 , Sp21 = Sp21, Sp20 = Sp20, 
              pi1 = 0.01, pi2 = 0.99)
N1 <- 5000; N2 <- 5000

require(prevalence)
parSe1 <- betaExpert(best = 0.85, lower = 0.5, p = 0.95, method = "mean")
parSp1 <- betaExpert(best = 0.8, lower = 0.45, p = 0.95, method = "mean")

#parSe1$alpha
#parSe1$beta

#parSp1$alpha
#parSp1$beta

library(rjags)
model_string <- "
model{
  # Priors
  Se1 ~ dbeta(3.18, 0.56)
  Sp1 ~ dbeta(3.5, 0.87)
  Se21 ~ dbeta(1, 1)
  Se20 ~ dbeta(1, 1)
  Sp21 ~ dbeta(1, 1)
  Sp20 ~ dbeta(1, 1)
  #pi1 ~ dbeta(1, 1)
  #pi2 ~ dbeta(1, 1)
  pi1 ~ dunif(0, 0.3)
  pi2 ~ dunif(0.7, 1)
  
  # Likelihood (2×2 tables per population)
  n1[1:4] ~ dmulti(p1[1:4], N1)
  n2[1:4] ~ dmulti(p2[1:4], N2)

  # Population 1 cell probabilities: (Y1,Y2) = (1,1),(1,0),(0,1),(0,0)
  p1[1] <- pi1*Se1*Se21 + (1-pi1)*(1-Sp1)*(1-Sp21)
  p1[2] <- pi1*Se1*(1-Se21) + (1-pi1)*(1-Sp1)*Sp21
  p1[3] <- pi1*(1-Se1)*Se20 + (1-pi1)*Sp1*(1-Sp20)
  p1[4] <- pi1*(1-Se1)*(1-Se20) + (1-pi1)*Sp1*Sp20

  # Population 2 (same Se/Sp, different prevalence)
  p2[1] <- pi2*Se1*Se21 + (1-pi2)*(1-Sp1)*(1-Sp21)
  p2[2] <- pi2*Se1*(1-Se21) + (1-pi2)*(1-Sp1)*Sp21
  p2[3] <- pi2*(1-Se1)*Se20 + (1-pi2)*Sp1*(1-Sp20)
  p2[4] <- pi2*(1-Se1)*(1-Se20) + (1-pi2)*Sp1*Sp20

  Se2 <- Se1*Se21 + Se20*(1-Se1)
  Sp2  <- Sp21*(1 - Sp1) + Sp20*Sp1
  corD <- (Se1*Se21 - Se1*Se2)/sqrt(Se1*(1 - Se1)*Se2*(1 - Se2))
  cornD <- (Sp20*Sp1 - Sp1*Sp2)/sqrt(Sp1*(1 - Sp1)*Sp2*(1 - Sp2))
}
"

# One Simulated dataset
data <- sim_georgiadis(N1 = N1, N2 = N2, 
               Se1 = truth$Se1, Se21 = truth$Se21, Se20 = truth$Se20,
               Sp1 = truth$Sp1, Sp21 = truth$Sp21, Sp20 = truth$Sp20,
               pi1 = truth$pi1, pi2 = truth$pi2, 
               seed = NULL)

# JAGS data list (order (11,10,01,00)) 
jdata <- list(
  n1 = data$n1, n2 = data$n2,
  N1 = sum(data$n1), N2 = sum(data$n2)
)

# Initial values 
inits_fun <- function() list(
  Se1 = runif(1, 0.5, 0.999),
  Sp1 = runif(1, 0.5, 0.999),
  Se21 = runif(1, 0.25, 0.999),
  Se20 = runif(1, 0.25, 0.999),
  Sp21 = runif(1, 0.25, 0.999),
  Sp20 = runif(1, 0.25, 0.999),
#  pi1 = runif(1, 0.05, 0.95),
 # pi2 = runif(1, 0.05, 0.95)
pi1 = runif(1, 0, 0.3),
pi2 = runif(1, 0.7, 1)
)

jm <- jags.model(textConnection(model_string),
                 data = jdata, inits = inits_fun, n.chains = 1, n.adapt = 2000)

# Burn-in
update(jm, 5000) 

samples <- coda.samples(jm,
                        variable.names = c("Se1","Sp1","Se2","Sp2", "pi1", "pi2", "corD", "cornD"),
                        n.iter = 50000)

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

data.frame(
  Truth     = c(truth$Se1, Se2, truth$Sp1, Sp2, rhoD, rhonD, truth$pi1, truth$pi2),
  Estimate  = round(post_summary, 3)
)


#post_all <- as.matrix(samples)

#par(mfrow = c(3,3))
#plot(density(post_all[,"Se1"]),
#     main = "Posterior density of Se1",
#     xlab = "Se1", lwd = 2, col = "blue",
#     xlim = c(0,1))

#plot(density(post_all[,"Sp1"]),
#     main = "Posterior density of Sp1",
#    xlab = "Sp1", lwd = 2, col = "blue",
#     xlim = c(0,1))

#plot(density(post_all[,"Se2"]),
#     main = "Posterior density of Se2",
#     xlab = "Se2", lwd = 2, col = "blue",
#     xlim = c(0,1))

#plot(density(post_all[,"Sp2"]),
#     main = "Posterior density of Sp2",
#     xlab = "Sp2", lwd = 2, col = "blue",
#     xlim = c(0,1))

#plot(density(post_all[,"pi1"]),
#     main = "Posterior density of Pi1",
#     xlab = "Pi1", lwd = 2, col = "blue",
#     xlim = c(0,1))

#plot(density(post_all[,"pi2"]),
#     main = "Posterior density of Pi2",
#     xlab = "Pi2", lwd = 2, col = "blue",
#     xlim = c(0,1))
