
### Fluency

library(tidyr)
library(spatstat) 
#library("foreign") 
#library(degreenet)

## for start values:
source("Responsefunctions.R")
source("Difficultyfunctions.R")
source("ProgramsBsplinesPureResp.R")
source("Posteriorgen.R")
source("Posteriorgen.R")

## for Bayes
source("ProgrThresholdsBayesFlex.R")

###################

############### verbal fluency 

#datvf<-read.table("C:/Users/tutz/LRZ Sync+Share/TuRaschThresholdsM2/R/vfnew.dat")

datvf<-read.table("vfnew.dat")
summary(datvf)
dataw <- spread(datvf, item, value)
dataw
dat<- na.omit(dataw)
summary(dat)
newdata <- dat[c(2:5)]
dat <- t(newdata) ##item x obs

###specification

lin <- "lin"  
#lin <- "log"
#lin <- "log1"
#lin <- "logit"

respf<- "NV"
#respf<- "logistic"
#respf<- "Gumbel"
#respf<- "Gompertz"

#indicator<- "C"  #### responses are continuous
#indicator<- "D"
I<- dim(dat)[1]
indicatorsingle  <- "C"  #### discrete or continuous
indicatorvec  <- indicatorsingle
for (l in 2:I )indicatorvec<- c(indicatorvec,indicatorsingle)


#### fit
cov0<-matrix(0,1,dim(dat)[2]) ### dummy covariates

#ThrNV1<-ThreshbasisRespCov(dat,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=0, 
#                          startind = "no",der="der",cov0,hessiander="FALSE")

ThrNV<-ThreshIncl(respf,lin,dat,basis="lin",indicatorvec,ord, numknotsstart,lambdpenin=0,lambdc=0,start=0, 
     startind = "no",der="der",cov0,hessiander="FALSE")
ThrNV$parmatrix 

### for Bayes transposed data 

dat<-t(dat)



####  two-dimensional
start2 <- list(
  theta = matrix(0, nrow = nrow(dat), ncol = 2),
  
  alpha = matrix(0, nrow = ncol(dat), ncol = 2),
  
  #delta0 = colMeans(dat),
  delta0=ThrNV$parmatrix[,1],
  #lambda = log(1 / apply(dat, 2, sd))
  lambda = log(ThrNV$parmatrix[,2])
)


fit2 <- fit.threshold(
  data = dat,
  F = pnorm,
  f = dnorm,
  g = g,
  gprime = gprime,
  q = 2,
  
  dalpha = dalpha.normal,
  alpha.par = list(
    mu = c(0, 0),
    Sigma = diag(2)
  ),
  start = start2,
  #niter = 10000,
  #burnin = 2000,
  niter = 3000,
  burnin = 500,
  thin = 5
)


### discrete
#g <- lapply(1:ncol(dat), function(i) identity)
g <- rep(list(function(y) log(1+y)), dim(dat)[2])
gprime <- rep(list(function(y) 1/(1+y)), dim(dat)[2])

fit2discr <-  fit.thresholdFlex(
  data = dat,
  F = pnorm,
  f = dnorm,
  g = g,
  q = 2,
  dalpha = dalpha.normal,
  alpha.par = list(
    mu = c(0, 0),
    Sigma = diag(2)
  ),
  discrete = TRUE,
  niter = 1000,
  burnin = 500
)

fit2discr$estimates$alpha
fit2discr$estimates$delta0
fit2discr$estimates$delta

fit2discr <-fit.thresholdFlex(data,
                              F = pnorm,
                              f = dnorm,
                              g = g,
                              gprime = NULL,
                              q = 2,
                              dalpha,
                              alpha.par = list(
                                mu = c(0, 0),
                                Sigma = diag(2)
                              ),
                              mu.delta0 = 0,
                              sd.delta0 = 2,
                              mu.lambda = 0,
                              sd.lambda = 1,
                              niter = 100,
                              burnin = 10,
                              thin = 1,
                              proposal.sd = 0.15,
                              start = start2,
                              discrete = TRUE,
                              verbose = TRUE) 

##### loglik
fit1used<-fit2discr
ll.fit <- loglik.thresholdFlex(
  data = dat,
  theta = fit1used$estimates$theta,
  alpha = fit1used$estimates$alpha,
  delta0 = fit1used$estimates$delta0,
  lambda = fit1used$estimates$lambda,
  F = pnorm,f = dnorm,  g = g,  discrete = TRUE)
ll.fit




### rotation quartimax
#library(psych)
#library(GPArotation)
#fitused<-fit2
fitused<-fit2discr
A <- as.matrix(fitused$estimates$alpha)
Theta <- as.matrix(fitused$estimates$theta)

rot <- GPArotation::quartimax(A)
T <- rot$Th
A.rot <- A %*% T
Theta.rot <- Theta %*% T

rownames(A.rot) <- colnames(dat)
A.rot
round(-A.rot, 3)

#fit2NVlog<-fit2  ##two dim fit NV, lin (F)
 