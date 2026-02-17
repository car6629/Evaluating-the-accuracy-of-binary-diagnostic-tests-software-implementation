require(prevalence)
library(coda)

n111 <- 67; n101 <- 25; n011 <- 41; n001 <- 329
n112 <- 97; n102 <- 33; n012 <- 36; n002 <- 371
n1pp <- n101 + n102 + n111 + n112
n0pp <- n001 + n102 + n012 + n002
n1 <- n111 + n101 + n011 + n001
n2 <- n112 + n102 + n012 + n002


aSe1 <- 24.09; bSe1 <- 5.73
aSp1 <- 23.05; bSp1 <- 3.45
aSe21 <- bSe21 <- 1
aSe20 <- bSe20 <- 1
aSp21 <- bSp21 <- 1
aSp20 <- bSp20 <- 1
api1 <- 1.3; bpi1 <- 5
api2 <- 1.5; bpi2 <- 3

S <- 100000
Se1 <- Sp1 <- pi1 <- pi2 <- Se21 <- Se20 <- Sp21 <- Sp20 <- numeric(S)
pi1[1] <- pi2[1] <- 0.3
Se1[1] <- Sp1[1] <- Se21[1] <- Se20[1] <- Sp21[1] <- Sp20[1] <- 0.5

set.seed(123)
for(i in 2:S){
  aux111 <- pi1[i-1]*Se1[i-1]*Se21[i-1] /
    (pi1[i-1]*Se1[i-1]*Se21[i-1] +
       (1-pi1[i-1])*(1-Sp1[i-1])*(1-Sp21[i-1]))
                                           
  aux101 <- pi1[i-1]*Se1[i-1]*(1-Se21[i-1])/
    (pi1[i-1]*Se1[i-1]*(1-Se21[i-1]) +
       (1-pi1[i-1])*(1-Sp1[i-1])*Sp21[i-1])
  
  aux011 <- pi1[i-1]*(1-Se1[i-1])*Se20[i-1]/
    (pi1[i-1]*(1-Se1[i-1])*Se20[i-1] +
       (1-pi1[i-1])*Sp1[i-1]*(1-Sp20[i-1]))
  
  aux001 <- pi1[i-1]*(1-Se1[i-1])*(1-Se20[i-1])/
    (pi1[i-1]*(1-Se1[i-1])*(1-Se20[i-1])+
       (1-pi1[i-1])*Sp1[i-1]*Sp20[i-1])
  
  z111 <- rbinom(1, n111, aux111)
  z101 <- rbinom(1, n101, aux101)
  z011 <- rbinom(1, n011, aux011)
  z001 <- rbinom(1, n001, aux001)
  
  aux112 <- pi2[i-1]*Se1[i-1]*Se21[i-1] /
    (pi2[i-1]*Se1[i-1]*Se21[i-1] +
       (1-pi2[i-1])*(1-Sp1[i-1])*(1-Sp21[i-1]))
  
  aux102 <- pi2[i-1]*Se1[i-1]*(1-Se21[i-1])/
    (pi2[i-1]*Se1[i-1]*(1-Se21[i-1]) +
       (1-pi2[i-1])*(1-Sp1[i-1])*Sp21[i-1])
  
  aux012 <- pi2[i-1]*(1-Se1[i-1])*Se20[i-1]/
    (pi2[i-1]*(1-Se1[i-1])*Se20[i-1] +
       (1-pi2[i-1])*Sp1[i-1]*(1-Sp20[i-1]))
  
  aux002 <- pi2[i-1]*(1-Se1[i-1])*(1-Se20[i-1])/
    (pi2[i-1]*(1-Se1[i-1])*(1-Se20[i-1])+
       (1-pi2[i-1])*Sp1[i-1]*Sp20[i-1])
  
  z112 <- rbinom(1, n112, aux112)
  z102 <- rbinom(1, n102, aux102)
  z012 <- rbinom(1, n012, aux012)
  z002 <- rbinom(1, n002, aux002)
  
  z1pp <- z101 + z102 + z111 + z112
  z0pp <- z001 + z002 + z011 + z012
  zpp1 <- z111 + z101 + z011 + z001
  zpp2 <- z112 + z102 + z012 + z002
  
  Se1[i] <- rbeta(1, aSe1 + z1pp, bSe1 + z0pp)
  Sp1[i] <- rbeta(1, aSp1 + n0pp - z0pp, bSp1 + n1pp - z1pp)
  Se21[i] <- rbeta(1, aSe21 + z111 + z112, bSe21 + z101 + z102)
  Se20[i] <- rbeta(1, aSe20 + z011 + z012, bSe20 + z001 + z002)
  Sp21[i] <- rbeta(1, aSp21 + n101 + n102 - z101 - z102, bSp21 + n111 + n112 - z111 - z112)
  Sp20[i] <- rbeta(1, aSp20 + n001 + n002 - z001 - z002, bSp20 + n011 + n012 - z011 - z012)
  pi1[i] <- rbeta(1, api1 + zpp1, bpi1 + n1 - zpp1)
  pi2[i] <- rbeta(1, api2 + zpp2, bpi2 + n2 - zpp2)
}

Se2 <- Se1 * Se21 + Se20 * (1 - Se1)
Sp2 <- Sp21 * (1 - Sp1) + Sp20 * Sp1




median(Se1)
median(Se2)
median(Sp1)
median(Sp2)
median(pi1)
median(pi2)

quantile(Se1, probs = c(0.025, 0.975))
quantile(Se2, probs = c(0.025, 0.975))
quantile(Sp1, probs = c(0.025, 0.975))
quantile(Sp2, probs = c(0.025, 0.975))
quantile(pi1, probs = c(0.025, 0.975))
quantile(pi2, probs = c(0.025, 0.975))

effectiveSize(mcmc(cbind(Se1, Se2, Sp1, Sp2, pi1, pi2)))




