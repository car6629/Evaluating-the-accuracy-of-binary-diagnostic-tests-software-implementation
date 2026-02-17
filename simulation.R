# 加载必要的库
library(rjags)
library(coda)
library(prevalence)


pri_Se1 <- betaExpert(best = 0.55, upper = 0.85, p = 0.95, method = "mode")
aSe1 <- round(pri_Se1$alpha, 1); bSe1 <- round(pri_Se1$beta, 1)


pri_Se2 <- betaExpert(best = 0.9, lower = 0.6, p = 0.95, method = "mode")
aSe2 <- round(pri_Se2$alpha, 1); bSe2 <- round(pri_Se2$beta, 1)


pri_Sp1 <- betaExpert(best = 0.98, lower = 0.8, p = 0.95, method = "mode")
aSp1 <- round(pri_Sp1$alpha, 1); bSp1 <- round(pri_Sp1$beta, 1)


pri_Sp2 <- betaExpert(best = 0.85, lower = 0.6, p = 0.95, method = "mode")
aSp2 <- round(pri_Sp2$alpha, 1); bSp2 <- round(pri_Sp2$beta, 1)


par_Pi1 <- betaExpert(best = 0.03, upper = 0.3, p = 0.95, method = "mean")
aPi1 <- round(par_Pi1$alpha, 1); bPi1 <- round(par_Pi1$beta, 1)


par_Pi2 <- betaExpert(best = 0.3, lower = 0.08, p = 0.95, method = "mode")
aPi2 <- round(par_Pi2$alpha, 1); bPi2 <- round(par_Pi2$beta, 1)


priors_table <- data.frame(
  Variable = c("Se1", "Se2", "Sp1", "Sp2", "Pi1", "Pi2"),
  Alpha = c(aSe1, aSe2, aSp1, aSp2, aPi1, aPi2),
  Beta  = c(bSe1, bSe2, bSp1, bSp2, bPi1, bPi2)
)
print("--- Calculated Expert Priors ---")
print(priors_table)


set.seed(123)
N1 <- 5000; N2 <- 5000
truth <- list(
  pi1 = 0.3, pi2 = 0.7,     
  Se1 = 0.85, Se2 = 0.95, rhoD = 0.5,
  Sp1 = 0.80, Sp2 = 0.90, rhonD = 0.6
)


Se11 <- truth$rhoD * sqrt(truth$Se1*(1-truth$Se1)*truth$Se2*(1-truth$Se2)) + truth$Se1*truth$Se2
Se21 <- Se11 / truth$Se1
Se20 <- (truth$Se2 - Se11) / (1 - truth$Se1)

Sp00 <- truth$rhonD * sqrt(truth$Sp1*(1-truth$Sp1)*truth$Sp2*(1-truth$Sp2)) + truth$Sp1*truth$Sp2
Sp20 <- Sp00 / truth$Sp1
Sp21 <- (truth$Sp2 - Sp00) / (1 - truth$Sp1)


calc_probs <- function(pi, Se1, Se21, Se20, Sp1, Sp21, Sp20) {
  c(pi*Se1*Se21 + (1-pi)*(1-Sp1)*(1-Sp21),
    pi*Se1*(1-Se21) + (1-pi)*(1-Sp1)*Sp21,
    pi*(1-Se1)*Se20 + (1-pi)*Sp1*(1-Sp20),
    pi*(1-Se1)*(1-Se20) + (1-pi)*Sp1*Sp20)
}

p1_true <- calc_probs(truth$pi1, truth$Se1, Se21, Se20, truth$Sp1, Sp21, Sp20)
p2_true <- calc_probs(truth$pi2, truth$Se1, Se21, Se20, truth$Sp1, Sp21, Sp20)

n1_obs <- as.vector(rmultinom(1, N1, p1_true))
n2_obs <- as.vector(rmultinom(1, N2, p2_true))


jdata <- list(
  n1 = n1_obs, n2 = n2_obs, N1 = N1, N2 = N2,
  aSe1=aSe1, bSe1=bSe1, aSe2=aSe2, bSe2=bSe2,
  aSp1=aSp1, bSp1=bSp1, aSp2=aSp2, bSp2=bSp2,
  aPi1=aPi1, bPi1=bPi1, aPi2=aPi2, bPi2=bPi2
)



# Model 1: Independence Model
model_indep_string <- "
model {
  Se1 ~ dbeta(aSe1, bSe1); Sp1 ~ dbeta(aSp1, bSp1)
  Se2 ~ dbeta(aSe2, bSe2); Sp2 ~ dbeta(aSp2, bSp2)
  pi1 ~ dbeta(aPi1, bPi1); pi2 ~ dbeta(aPi2, bPi2)
  
  p1[1] <- pi1*Se1*Se2 + (1-pi1)*(1-Sp1)*(1-Sp2)
  p1[2] <- pi1*Se1*(1-Se2) + (1-pi1)*(1-Sp1)*Sp2
  p1[3] <- pi1*(1-Se1)*Se2 + (1-pi1)*Sp1*(1-Sp2)
  p1[4] <- pi1*(1-Se1)*(1-Se2) + (1-pi1)*Sp1*Sp2
  
  p2[1] <- pi2*Se1*Se2 + (1-pi2)*(1-Sp1)*(1-Sp2)
  p2[2] <- pi2*Se1*(1-Se2) + (1-pi2)*(1-Sp1)*Sp2
  p2[3] <- pi2*(1-Se1)*Se2 + (1-pi2)*Sp1*(1-Sp2)
  p2[4] <- pi2*(1-Se1)*(1-Se2) + (1-pi2)*Sp1*Sp2
  
  n1[1:4] ~ dmulti(p1[1:4], N1)
  n2[1:4] ~ dmulti(p2[1:4], N2)
}"

# Model 2: Dependence Model (Georgiadis)
model_dep_string <- "
model {
  Se1 ~ dbeta(aSe1, bSe1); Sp1 ~ dbeta(aSp1, bSp1)
  Se21 ~ dbeta(1,1); Se20 ~ dbeta(1,1)
  Sp21 ~ dbeta(1,1); Sp20 ~ dbeta(1,1)
  pi1 ~ dbeta(aPi1, bPi1); pi2 ~ dbeta(aPi2, bPi2)
  
  p1[1] <- pi1*Se1*Se21 + (1-pi1)*(1-Sp1)*(1-Sp21)
  p1[2] <- pi1*Se1*(1-Se21) + (1-pi1)*(1-Sp1)*Sp21
  p1[3] <- pi1*(1-Se1)*Se20 + (1-pi1)*Sp1*(1-Sp20)
  p1[4] <- pi1*(1-Se1)*(1-Se20) + (1-pi1)*Sp1*Sp20
  
  p2[1] <- pi2*Se1*Se21 + (1-pi2)*(1-Sp1)*(1-Sp21)
  p2[2] <- pi2*Se1*(1-Se21) + (1-pi2)*(1-Sp1)*Sp21
  p2[3] <- pi2*(1-Se1)*Se20 + (1-pi2)*Sp1*(1-Sp20)
  p2[4] <- pi2*(1-Se1)*(1-Se20) + (1-pi2)*Sp1*Sp20
  
  n1[1:4] ~ dmulti(p1[1:4], N1)
  n2[1:4] ~ dmulti(p2[1:4], N2)
  

  Se2 <- Se1*Se21 + (1-Se1)*Se20
  Sp2 <- Sp1*Sp20 + (1-Sp1)*Sp21
  corD <- (Se1*Se21 - Se1*Se2) / sqrt(Se1*(1-Se1)*Se2*(1-Se2) + 0.000001)
  cornD <- (Sp1*Sp20 - Sp1*Sp2) / sqrt(Sp1*(1-Sp1)*Sp2*(1-Sp2) + 0.000001)
}"



jm_i <- jags.model(textConnection(model_indep_string), data=jdata, n.chains=2, n.adapt=2000)
update(jm_i, 5000)
res_indep <- coda.samples(jm_i, variable.names=c("Se1","Se2","Sp1","Sp2","pi1","pi2"), n.iter=20000)


jm_d <- jags.model(textConnection(model_dep_string), data=jdata, n.chains=2, n.adapt=2000)
update(jm_d, 5000)
res_dep <- coda.samples(jm_d, variable.names=c("Se1","Se2","Sp1","Sp2","pi1","pi2","corD","cornD"), n.iter=20000)


sum_i <- summary(res_indep)$quantiles[, c("2.5%", "50%", "97.5%")]
sum_d <- summary(res_dep)$quantiles[, c("2.5%", "50%", "97.5%")]


compare_all <- data.frame(
  Param = c("Se1", "Se2", "Sp1", "Sp2", "pi1", "pi2"),
  Truth = c(truth$Se1, truth$Se2, truth$Sp1, truth$Sp2, truth$pi1, truth$pi2),
  Indep_Med = round(sum_i[c("Se1","Se2","Sp1","Sp2","pi1","pi2"), "50%"], 3),
  Dep_Med   = round(sum_d[c("Se1","Se2","Sp1","Sp2","pi1","pi2"), "50%"], 3),
  Dep_Lower = round(sum_d[c("Se1","Se2","Sp1","Sp2","pi1","pi2"), "2.5%"], 3),
  Dep_Upper = round(sum_d[c("Se1","Se2","Sp1","Sp2","pi1","pi2"), "97.5%"], 3)
)

print("--- Final Comparison Table ---")
print(compare_all)

cat("\n--- Estimated Correlations (Dep Model) ---\n")
print(round(sum_d[c("corD", "cornD"), ], 3))