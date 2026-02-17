require(rjags)     # PACKAGE TO RUN THE jags MODEL. MANDATORY
require(MCMCvis)   # THIS PACKAGE CONTAINS THE MCMCsummary FUNCTION USED IN THIS SCRIPT
require(mcmcplots) # USED FOR THE CREATION OF THE CONVERGENCE PLOTS
require(DT)        # THIS LIBRARY ALLOWS A NICE DATA DISPLAY WITH THE SEARCH BAR OPTION.

N <- c(38,87, 2, 35)

#Single patient result for the two tests
y1 <- c(rep(1,N[1]),rep(1,N[2]),rep(0,N[3]),rep(0,N[4]))
y2 <- c(rep(1,N[1]),rep(0,N[2]),rep(1,N[3]),rep(0,N[4]))

#Joint test results for each patient
y <- cbind(y1,y2)

#Number of patients
N = dim(y)[1]

dataList <- list(y = y,N = N)

#Bayesian Latent Class Random Effects Model
modelString = 
  "model{
  for (i in 1:N){
  
  #Likelihood
  
  D1[i]~dbern(prev)
  D[i]<-D1[i]+1
  y[i,1]~dbern(p1[i,D[i]])
  y[i,2]~dbern(p2[i,D[i]])
  
  #Defining the Se and Sp of each subject (include random effects)
  
  r[i]~dnorm(0,1)
  s1[i]<-phi(a[1,1]+b.RE[1]*r[i])
  c1[i]<-phi(a[1,2]+b.RE[2]*r[i])
  s2[i]<-phi(a[2,1]+b.RE[1]*r[i])
  c2[i]<-phi(a[2,2]+b.RE[2]*r[i])
  
  #Conditional probability of a positive observation
  p1[i,2]<-s1[i]
  p1[i,1]<-1-c1[i]
  p2[i,2]<-s2[i]
  p2[i,1]<-1-c2[i]
  }
  
  #Prior distribution
  
  prev~dbeta(1,1)
  se[1]~dbeta(21.96,5.49)
  sp[1]~dbeta(4.1,1.76) 
  se[2]~dbeta(4.44,13.31)
  sp[2]~dbeta(71.25,3.75)
  
  a[1,1]<-probit(se[1])*sqrt(1+b.RE[1]*b.RE[1])
  a[1,2]<-probit(sp[1])*sqrt(1+b.RE[1]*b.RE[1])
  a[2,1]<-probit(se[2])*sqrt(1+b.RE[2]*b.RE[2])
  a[2,2]<-probit(sp[2])*sqrt(1+b.RE[2]*b.RE[2])
  
  b.RE[1]~dnorm(0,1)I(0,)
  b.RE[2]~dnorm(0,1)I(0,)  
  }
"


writeLines(modelString,con="model.txt")


GenInits  = function(){
  
  se1 <- rbeta(1,21.96,5.49)
  sp1 <- rbeta(1,4.1,1.76) 
  se2 <- rbeta(1,4.44,13.31)
  sp2 <- rbeta(1,71.25,3.75)
  b1 <- abs(rnorm(1,0,1))
  b2 <- abs(rnorm(1,0,1))
  prev <- rbeta(1,1,1)
  
  se <- c(se1, se2)
  sp <- c(sp1, sp2)
  b.RE <- c(b1, b2)
  
  list(
    se = se,
    sp = sp,
    b.RE = b.RE,
    prev = prev,
    .RNG.name="base::Wichmann-Hill",
    .RNG.seed=321
  )
  
}

# Initial values
set.seed(123)
initsList = vector('list',3)
for(i in 1:3){
  initsList[[i]] = GenInits()
}

# Compile the model
jagsModel = jags.model("model.txt",data=dataList,n.chains=3,n.adapt=0, inits=initsList)

#jagsModel$state(internal=FALSE)

# Burn-in iterations 
update(jagsModel,n.iter=5000)
# Parameters to be monitored
parameters = c( "se","sp", "a", "b.RE", "prev")

# Posterior samples
posterior_results = coda.samples(jagsModel,variable.names=parameters,n.iter=5000)
output = posterior_results

res = MCMCsummary(output, digits=4)
datatable(res, extensions = 'AutoFill')












