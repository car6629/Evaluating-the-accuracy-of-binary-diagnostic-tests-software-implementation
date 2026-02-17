library(rjags)
library(coda)
N <- c(38, 87, 2, 35)
obs <- list(N = N, n = sum(N))

#prior <- list(aSe = c(4.44, 21.96), bSe = c(13.31, 5.49),
#              aSp = c(71.25, 4.1), bSp = c(3.75, 1.76),
#              aPi = 1, bPi = 1)

prior <- list(aSe = c(1, 1), bSe = c(1, 1),
              aSp = c(1, 1), bSp = c(1, 1),
              aPi = 1, bPi = 1)

model_string <- "model{
N ~ dmulti(p[1:4],n)
p[1] = Pi * Se[1] * Se[2] + (1 - Pi) * (1 - Sp[1]) * (1 - Sp[2])
p[2] = Pi * Se[1] * (1 - Se[2]) + (1 - Pi) * (1 - Sp[1]) * Sp[2]
p[3] = Pi * (1 - Se[1]) * Se[2] + (1 - Pi) * Sp[1] * (1 - Sp[2])
p[4] = Pi * (1 - Se[1]) * (1 - Se[2]) + (1 - Pi) * Sp[1] * Sp[2]

for(i in 1:2){
Se[i] ~ dbeta(aSe[i], bSe[i])
Sp[i] ~ dbeta(aSp[i], bSp[i])
}
Pi ~ dbeta(aPi, bPi)
}"

inits <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 123)
model <- jags.model(textConnection(model_string), data = c(obs, prior),
                    inits = inits)
samples <- coda.samples(model, variable.names = c("Se","Sp","Pi"),
                        n.iter = 50000)
summary(samples)$quantiles[,c(1,3,5)]
effectiveSize(samples)
