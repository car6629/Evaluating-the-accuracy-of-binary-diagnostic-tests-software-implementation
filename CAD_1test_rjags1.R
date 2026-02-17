library(rjags)

#data
obs <- list(n1plus = 930, n = 1465)
#prior
prior <- list(aSe = 103.7, bSe = 26.7, aSp = 24.3, bSp = 9.2,
              api = 143.5, bpi = 62.1)

#prior <- list(aSe = 1, bSe = 1, aSp = 1, bSp = 1,
#             api = 1, bpi = 1)

#JAGS model
model_string <- "model{
#likelihood
n1plus ~ dbinom(p, n)
p = Pi * Se + (1 - Pi) * (1 - Sp)
#prior
Se ~ dbeta(aSe, bSe)
Sp ~ dbeta(aSp, bSp)
Pi ~ dbeta(api, bpi)
#derived parameters
PPV = Pi * Se / p
NPV = (1 - Pi) * Sp / (1 - p)
}"

inits <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 123)
model <- jags.model(textConnection(model_string), data = c(obs, prior), 
                    inits = inits)

samples <- coda.samples(model,
                        variable.name = c("Se", "Sp", "Pi", "PPV", "NPV"),
                        n.iter = 100000)
summary(samples)$quantiles[,c(1,3,5)]


effectiveSize(samples)[c("Se", "Sp", "Pi")]
