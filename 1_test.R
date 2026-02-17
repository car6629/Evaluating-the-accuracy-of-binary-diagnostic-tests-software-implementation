n1plus <- 930; n <- 1465
require(prevalence)
par_pi <- betaExpert(best = 0.7, upper = 0.72, p = 0.75,
                     method = "mode")
api <- round(par_pi$alpha, 1);
bpi <- round(par_pi$beta, 1)

par_Se <- betaExpert(best = 0.8, upper = 0.82, p = 0.75,
                     method = "mode")
aSe <- round(par_Se$alpha, 1);
bSe <- round(par_Se$beta, 1)

par_Sp <- betaExpert(best = 0.74, upper = 0.78, p = 0.75,
                     method = "mode")
aSp <- round(par_Sp$alpha, 1);
bSp <- round(par_Sp$beta, 1)

S <- 20000
Se <- Sp <- pi <- PPV <- NPV <- numeric(S)

Se[1] <- Sp[1] <- pi[1] <- PPV[1] <- NPV[1] <- 0.5

set.seed(123)
for(i in 2:S){
  PPV[i] <- pi[i-1]*Se[i-1]/(pi[i-1]*Se[i-1] + 
                               (1-pi[i-1])*(1-Sp[i-1]))
  NPV[i] <- Sp[i-1]*(1-pi[i-1])/((1-pi[i-1])*Sp[i-1]+
                                   pi[i-1]*(1-Se[i-1]))
  z0 <- rbinom(1, n - n1plus, NPV[i])
  z1 <- rbinom(1, n1plus, PPV[i])
  Se[i] <- rbeta(1, aSe + z1, bSe + n - n1plus - z0)
  Sp[i] <- rbeta(1, aSp + z0, bSp +n1plus - z1)
  pi[i] <- rbeta(1, api + z1 +n - n1plus - z0,
                 bpi + n1plus - z1 +z0)
}


Se[2000]
Sp[2000]
pi[2000]






