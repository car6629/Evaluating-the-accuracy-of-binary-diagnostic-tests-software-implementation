library(rjags)
library(coda)
N <- c(38, 2, 87, 35)
obs <- list(N = N, n = sum(N))

#prior <- list(aSe = c(4.44, 21.96), bSe = c(13.31, 5.49),
#              aSp = c(71.25, 4.1), bSp = c(3.75, 1.76),
#              aPi = 1, bPi = 1)

model_string <- "model{
# Priors
  Se1 ~ dbeta(4.44, 13.31)
  Sp1 ~ dbeta(71.25, 3.75)
  Se21 ~ dbeta(1, 1)
  Se20 ~ dbeta(1, 1)
  Sp21 ~ dbeta(1, 1)
  Sp20 ~ dbeta(1, 1)
  pi ~ dbeta(1, 1)

  # Likelihood (2×2 tables per population)
  N[1:4] ~ dmulti(p[1:4], n)

  # Population 1 cell probabilities: (Y1,Y2) = (1,1),(1,0),(0,1),(0,0)
  p[1] <- pi*Se1*Se21 + (1-pi)*(1-Sp1)*(1-Sp21)
  p[2] <- pi*Se1*(1-Se21) + (1-pi)*(1-Sp1)*Sp21
  p[3] <- pi*(1-Se1)*Se20 + (1-pi)*Sp1*(1-Sp20)
  p[4] <- pi*(1-Se1)*(1-Se20) + (1-pi)*Sp1*Sp20

  Se2 <- Se1*Se21 + Se20*(1-Se1)
  Sp2  <- Sp21*(1 - Sp1) + Sp20*Sp1
  rhoD <- (Se1*Se21 - Se1*Se2)/sqrt(Se1*(1 - Se1)*Se2*(1 - Se2))
  rhonD <- (Sp20*Sp1 - Sp1*Sp2)/sqrt(Sp1*(1 - Sp1)*Sp2*(1 - Sp2))
}"

inits <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 123)
model <- jags.model(textConnection(model_string),
                  data = obs, inits = inits)

samples <- coda.samples(model, variable.names = c("Se1","Se2","Sp1","Sp2","pi","rhoD","rhonD"),n.iter = 50000)
summary(samples)$quantiles[,c(1,3,5)]





