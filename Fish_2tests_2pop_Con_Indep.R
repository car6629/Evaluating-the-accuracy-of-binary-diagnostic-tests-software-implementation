require(prevalence)
library(coda)

aSp1 <- bSp1 <- aSe1 <- bSe1 <- 1
aSp2 <- bSp2 <- aSe2 <- bSe2 <- 1
aPi1 <- aPi2 <- bPi1 <- bPi2 <- 1

#pri_Sp1 <- betaExpert(best = 0.98, lower = 0.8, p = 0.95, method = "mode")
#aSp1 <- round(pri_Sp1$alpha, 1)
#bSp1 <- round(pri_Sp1$beta, 1)


#pri_Se1 <- betaExpert(best = 0.55, upper = 0.85, p = 0.95, method = "mode")
#aSe1 <- round(pri_Se1$alpha, 1)
#bSe1 <- round(pri_Se1$beta, 1)

#pri_Se2 <- betaExpert(best = 0.9, lower = 0.6, p = 0.95, method = "mode")
#aSe2 <- round(pri_Se2$alpha, 1)
#bSe2 <- round(pri_Se2$beta, 1)

#pri_Sp2 <- betaExpert(best = 0.85, lower = 0.6, p = 0.95, method = "mode")
#aSp2 <- round(pri_Sp2$alpha, 1)
#bSp2 <- round(pri_Sp2$beta, 1)

#par_Pi2 <- betaExpert(best = 0.3, lower = 0.08, p = 0.95, method = "mode")
#aPi2 <- round(par_Pi2$alpha, 1)
#bPi2 <- round(par_Pi2$beta, 1)

#par_Pi1 <- betaExpert(best = 0.03, upper = 0.3, p = 0.95, method = "mode")
#aPi1 <- round(par_Pi1$alpha, 1)
#bPi1 <- round(par_Pi1$beta, 1)

n111 <- 0; n101 <- 0; n011 <- 3; n001 <- 129
n112 <- 3; n102 <- 0; n012 <- 24; n002 <- 3

S <- 50000
Se1 <- Se2 <- Sp1 <- Sp2 <- Pi1 <- Pi2 <- numeric(S)
Se1[1] <- Se2[1] <- Sp1[1] <- Sp2[1] <- Pi1[1] <- Pi2[1] <- 0.5

set.seed(123)
for(i in 2:S){
  pz111 <- Pi1[i-1] * Se1[i-1] * Se2[i-1] / (Pi1[i-1] * Se1[i-1] * Se2[i-1]
                                           + (1-Pi1[i-1]) * (1-Sp1[i-1]) * (1-Sp2[i-1]))
  pz112 <- Pi2[i-1] * Se1[i-1] * Se2[i-1] / (Pi2[i-1] * Se1[i-1] * Se2[i-1]
                                           + (1-Pi2[i-1]) * (1-Sp1[i-1]) * (1-Sp2[i-1]))
  pz101 <- Pi1[i-1] * Se1[i-1] * (1-Se2[i-1]) / (Pi1[i-1] * Se1[i-1] * (1-Se2[i-1])
                                                 + (1-Pi1[i-1]) * (1-Sp1[i-1]) * Sp2[i-1])
  pz102 <- Pi2[i-1] * Se1[i-1] * (1-Se2[i-1]) / (Pi2[i-1] * Se1[i-1] * (1-Se2[i-1])
                                                 + (1-Pi2[i-1]) * (1-Sp1[i-1]) * Sp2[i-1])
  pz011 <- Pi1[i-1] * (1 - Se1[i-1]) * Se2[i-1] / (Pi1[i-1] * (1-Se1[i-1]) * Se2[i-1]
                                                   + (1-Pi1[i-1]) * Sp1[i-1] * (1-Sp2[i-1]))
  pz012 <- Pi2[i-1] * (1 - Se1[i-1]) * Se2[i-1] / (Pi2[i-1] * (1-Se1[i-1]) * Se2[i-1]
                                                    + (1-Pi2[i-1]) * Sp1[i-1] * (1-Sp2[i-1]))
  pz001 <- Pi1[i-1] * (1-Se1[i-1]) * (1-Se2[i-1]) / (Pi1[i-1] * (1-Se1[i-1]) * (1-Se2[i-1])
                                                     + (1-Pi1[i-1]) * Sp1[i-1] * Sp2[i-1]) 
  pz002 <- Pi2[i-1] * (1-Se1[i-1]) * (1-Se2[i-1]) / (Pi2[i-1] * (1-Se1[i-1]) * (1-Se2[i-1])
                                                     + (1-Pi2[i-1]) * Sp1[i-1] * Sp2[i-1])
  
  z111 <- rbinom(1, n111, pz111)
  z112 <- rbinom(1, n112, pz112)
  z101 <- rbinom(1, n101, pz101)
  z102 <- rbinom(1, n102, pz102)
  z011 <- rbinom(1, n011, pz011)
  z012 <- rbinom(1, n012, pz012)
  z001 <- rbinom(1, n001, pz001)
  z002 <- rbinom(1, n002, pz002)
  
  Se1[i] <- rbeta(1, aSe1 + z111 + z101 + z112 + z102, bSe1 + z011 + z001 + z012 + z002)
  Se2[i] <- rbeta(1, aSe2 + z111 + z011 + z112 + z012, bSe2 + z101 + z001 + z102 + z002)
  Sp1[i] <- rbeta(1, aSp1 + n011 + n001 + n012 + n002 - z011 - z001 - z012,
                  bSp1 + n111 + n101 + n112 + n102 - z111 - z101 - z112 - z102)
  Sp2[i] <- rbeta(1, aSp2 + n101 + n001 + n102 + n002 - z101 - z001 - z002,
                  bSp2 + n111 + n011 + n112 + n012 - z111 - z011 - z112 - z012)
  Pi1[i] <- rbeta(1, aPi1 + z111 + z101 + z011 + z001, bPi1 + n111 + n101
                  + n011 + n001 - z111 - z101 - z011 - z001)
  Pi2[i] <- rbeta(1, aPi2 + z112 + z102 + z012 + z002, bPi2 + n112 + n102
                  + n012 + n002 - z112 - z102 - z012 - z002)
  
}

median(Se1)
median(Se2)
median(Sp1)
median(Sp2)
median(Pi1)
median(Pi2)

quantile(Se1, probs = c(0.025, 0.975))
quantile(Se2, probs = c(0.025, 0.975))
quantile(Sp1, probs = c(0.025, 0.975))
quantile(Sp2, probs = c(0.025, 0.975))
quantile(Pi1, probs = c(0.025, 0.975))
quantile(Pi2, probs = c(0.025, 0.975))

effectiveSize(mcmc(cbind(Se1, Se2, Sp1, Sp2, Pi1, Pi2)))











