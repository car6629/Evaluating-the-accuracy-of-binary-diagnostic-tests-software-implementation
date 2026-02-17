require(prevalence)
library(coda)

n111 <- 67; n101 <- 25; n011 <- 41; n001 <- 329
n112 <- 97; n102 <- 33; n012 <- 36; n002 <- 371

aSe1 <- 24.09; bSe1 <- 5.73
aSp1 <- 23.05; bSp1 <- 3.45
aSe2 <- bSe2 <- 1
aSp2 <- bSp2 <- 1
aPi1 <- 1.3; bPi1 <- 5
aPi2 <- 1.5; bPi2 <- 3

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











