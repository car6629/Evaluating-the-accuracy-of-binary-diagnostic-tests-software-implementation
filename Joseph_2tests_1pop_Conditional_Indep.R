install.packages("coda")   
library(coda)

n11 <- 38; n10 <- 87; n01 <- 2; n00 <- 35
n <- 162
n1 <- 125; n2 <- 40

require(prevalence)
#prior
aSe2 = 4.44; bSe2 = 13.31 
aSp2 = 71.25; bSp2 = 3.75
aPi = 1; bPi = 1
aSe1 = 21.96; bSe1 = 5.49
aSp1 = 4.1; bSp1 = 1.78


S <- 50000
Se1 <- Se2 <- Pi <- Sp1 <- Sp2 <- numeric(S)
Se1[1] <- Se2[1] <- Pi[1] <- Sp1[1] <- Sp2[1] <- 0.5

set.seed(123)
for(i in 2:S){
    pz11 <- Pi[i - 1] * Se1[i - 1] * Se2[i - 1] / (Pi[i - 1] * Se1[i - 1] * Se2[i - 1] 
                                             + (1 - Pi[i - 1]) * (1 - Sp1[i - 1]) 
                                             * (1 - Sp2[i - 1]) ) 
    pz10 <- Pi[i - 1] * Se1[i - 1] * (1 - Se2[i - 1]) / (Pi[i - 1] * Se1[i - 1] * (1 - Se2[i - 1]) 
                                                  + (1 - Pi[i - 1]) * (1 - Sp1[i - 1]) 
                                                  * Sp2[i - 1])
    pz01 <- Pi[i - 1] * (1 - Se1[i - 1]) * Se2[i - 1] / (Pi[i - 1] * (1 - Se1[i - 1]) * Se2[i - 1] 
                                                        + (1 - Pi[i - 1]) * Sp1[i - 1]
                                                        * (1 - Sp2[i - 1]))
    pz00 <- Pi[i - 1] * (1 - Se1[i - 1]) * (1 - Se2[i - 1]) / (Pi[i - 1] * (1 - Se1[i - 1]) * (1 - Se2[i - 1])
                                                         + (1 - Pi[i - 1]) * Sp1[i - 1]
                                                         * Sp2[i - 1])
    
    z11 <- rbinom(1, n11, pz11)
    z10 <- rbinom(1, n10, pz10)
    z01 <- rbinom(1, n01, pz01)
    z00 <- rbinom(1, n00, pz00)
    
    Se1[i] <- rbeta(1, aSe1 + z11 + z10, bSe1 + z01 + z00)
    Se2[i] <- rbeta(1, aSe2 + z11 + z01, bSe2 + z10 + z00)
    Sp1[i] <- rbeta(1, aSp1 + n01 + n00 - z01 - z00, bSp1 + n11 + n10 - z11 - z10)
    Sp2[i] <- rbeta(1, aSp2 + n10 + n00 - z10 - z00, bSp2 + n11 + n01 - z11 - z01)
    Pi[i] <- rbeta(1, aPi + z11 + z10 + z01 + z00, bPi + n - z11 - z10 - z01 - z00)
}



median(Se1)
median(Se2)
median(Sp1)
median(Sp2)
median(Pi)

quantile(Se1, probs = c(0.025, 0.975))
quantile(Se2, probs = c(0.025, 0.975))
quantile(Sp1, probs = c(0.025, 0.975))
quantile(Sp2, probs = c(0.025, 0.975))
quantile(Pi, probs = c(0.025, 0.975))

effectiveSize(mcmc(cbind(Se1, Se2, Sp1, Sp2, Pi)))
