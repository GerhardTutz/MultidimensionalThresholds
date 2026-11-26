



library(FactoMineR)
library(dplyr)
library(spatstat) 
library(GPArotation)

## for start values:
source("Responsefunctions.R")
source("Difficultyfunctions.R")
source("ProgramsBsplinesPureResp.R")
source("Posteriorgen.R")
source("Posteriorgen.R")

## for Bayes
source("ProgrThresholdsBayesFlex.R")

data(decathlon)
decathlon



#deca<-decathlon[14:41,]  # only olympic
deca<-decathlon
deca<-data.frame(deca)
deca
cor(deca[,1:10])

### select tasks

###########################
## specification of data sets

##run
decarun<-select(deca,X100m,X400m,X110m.hurdle,X1500m)
indrun<-c(1,5,6,10)
#decarunneg<--decarun

#jump
decajump<-select(deca,Long.jump,High.jump,Pole.vault)
indjump<-c(2,4,8)

#throw
decathrow<-select(deca,Shot.put,Javeline,Discus)
indthrow<-c(3,9,7)

# all
decaall<-deca[,1:10]
decaall[,indrun]<--decaall[,indrun]
indall<-1:10

## almost all 8without pole valult, 1500 m
decaalmost<-decaall
decaalmost<-decaalmost[,-c(3,10)]
indalmost<-1:10
indalmost<-indalmost[-c(3,10)]

###########################

#selection of data sets: select only one  

#jump
decaselect<-decajump
indsel<-indjump

#run
decaselect<-decarun
#decaselect<-decarunneg
indsel<-indrun

#throw
decaselect<-decathrow
indsel<-indthrow

# all  
decaall<-deca[,1:10]
decaall[,indrun]<--decaall[,indrun]
indall<-1:10

# all new order
#decaallord <- deca[, c(1, 5, 6, 10,2,4,8,3,7,9)]
#decaallord[,1:4] <--decaallord[,1:4] 

decaallord <- deca[, c(1, 5, 6, 10,2,4,3,7,9,8)]
decaallord[,1:4] <--decaallord[,1:4] 
decaselect<-decaallord
indsel<-1:10

## almost all 8without pole valult, 1500 m
decaselect<-decaalmost
indsel<-indalmost

## corr
round(cor(deca[,indjump]), digits=2)  
round(cor(deca[,indrun]), digits=2) 
round(cor(deca[,indthrow]), digits=2) 


### specify difficulty and response function (default: lin <- "lin", respf<- "NV" )

lin <- "lin"  
#lin <- "log"
#lin <- "log1"
#lin <- "logit"

#respf<- "NV"
#respf<- "logistic"
#respf<- "Gumbel"
respf<- "Gompertz"

indicator<- "C"  #### responses are continuous
#indicator<- "D"


### always run (1) and (2) after selection of lin and respfct   

####(1) define function for response function

## first run "Responsefunctions"

if(respf=="NV"){respfctd <- respfctdNV
respfctc<-respfctcNV
respfctder<-respfctderNV}
if(respf=="logistic"){respfctd <- respfctdlogit
respfctc<-respfctclogit
respfctder<-respfctderlogit}
if(respf=="Gumbel"){respfctd <- respfctdGumbel
respfctc<-respfctcGumbel
respfctder<-respfctderGumbel}
if(respf=="Gompertz"){respfctd <- respfctdGompertz
respfctc<-respfctcGompertz
respfctder<-respfctderGompertz}

####(2) define function for difficulty function

## first run "Difficultyfunctions"

if(lin =="lin") {diffunct<-diffunctlin
derdiffunct<-derdiffunctlin}
if(lin =="log") {diffunct<-diffunctlog
derdiffunct<-derdiffunctlog}
if(lin =="log1") {diffunct<-diffunctlog1
derdiffunct<-derdiffunctlog1}






#### fitting of data set decaselect

### needed:


decaselectorig<-decaselect
decaselect<-scale(decaselect)## scaled data!!!
decaselecttr<-t(decaselect)  # transposed

####


meandat<-colMeans(decaselectorig, na.rm = TRUE)
sddat<-apply(decaselectorig, 2, sd, na.rm = TRUE)
#standdatt<-scale(decaselect)

#(decaselect[1,]-meandat)/sddat
#decaselecttr[,1] 
####

I<- dim(decaselecttr)[1]
indicatorsingle  <- "C"  #### discrete or continuous
indicatorvec  <- indicatorsingle
for (l in 2:I )indicatorvec<- c(indicatorvec,indicatorsingle)
indicator<-indicatorvec



cov0<-matrix(0,1,dim(decaselecttr)[2])

# first run ProgramsBsplinesPureResp



## use if response function is Nv
ThrNV<-ThreshbasisRespCov(decaselecttr,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start, 
                          startind = "no",der="der",cov0,hessiander="FALSE")

ThrNV
ThrNV$parmatrix[,1]



###########################
#### Bayes


dat<-decaselect
#datsc<-scale(dat)


###  

### specification
g <- rep(list(function(y) y), dim(dat)[2])
gprime <- rep(list(function(y) 1), dim(dat)[2])

## only one
#respf<-"Gumbel"
#respf<-"Gompertz"
respf<-"Gompertz"
#respf<-"NV"

## NV
if(respf=="NV"){
  F <-function(x){
    r<-pnorm(x)  
    return(r)}
  f <-function(x){
    r<-dnorm(x)  
    return(r)  
  }}

## Gumbel
if(respf=="Gumbel"){
  F <-function(x){
    r<-exp(-exp (-x))  
    return(r)}
  f <-function(x){
    r<-exp(-x)*exp(-exp (-x))  
    return(r)  
  }}

##
if(respf=="Gompertz"){
  F <-function(x){
    r<-1-exp(-exp(x)) 
    return(r)  
  }
  f <-function(x){
    r<-exp(x)*exp(-exp (x)) 
    return(r)  
  }}

dalpha.normal <- function(x, par) {
  
  mu <- par$mu
  Sigma <- par$Sigma
  
  q <- length(x)
  R <- chol(Sigma)
  
  z <- backsolve(R, x - mu)
  
  exp(
    -0.5 * q * log(2 * pi) -
      sum(log(diag(R))) -
      0.5 * sum(z^2)
  )
}

####  one-dimensional

###  
startmat<-ThrNV$parmatrix
#startmat<-ThrGompertz$parmatrix

start1 <- list(
  theta = matrix(0, nrow = nrow(dat), ncol = 1),
  
  alpha = matrix(0, nrow = ncol(dat), ncol = 1),
  
  #delta0 = colMeans(dat),
  delta0=startmat[,1],
  #lambda = log(1 / apply(dat, 2, sd))
  lambda = log(startmat[,2])
)

#### sd.delta0 = 5!
fit1  <-fit.thresholdFlex(dat,
                          F = F,
                          f = f,
                          g = g,
                          gprime = gprime,
                          q = 1,
                          dalpha= dalpha.normal,
                          alpha.par = list(
                            mu = c(0),
                            Sigma = diag(1)
                          ),
                          mu.delta0 = 0,
                          sd.delta0 = 5,
                          mu.lambda = 0,
                          sd.lambda = 1,
                          niter = 1000,
                          burnin = 300,
                          thin = 1,
                          proposal.sd = 0.15,
                          start=start1,
                          #discrete = TRUE,
                          verbose = TRUE) 

fit1$estimates$alpha
fit1$estimates$delta0
fit1$estimates$delta

cbind(fit1$estimates$delta0,fit1$estimates$delta)
ThrNV$parmatrix


############## plot
fitused<-fit1
#fitused<-fit2
parmatrest <-cbind(fitused$estimates$delta0,fitused$estimates$delta) ### for plots

#parmatrest <-ThrNV$parmatrix[,1:2]

#datt<-t(dat)
datn<-decaselecttr 


pcdum <-c ('1','2','3','4','5','6','7','8','9')

theta <-0 # specify theta

plotitems<-1:dim(datn)[1]  ## all items

par(cex.lab=1,cex.axis=1.8,cex.main=2.20,pch=16,cex=1.0)
par(mar=c(7,6,6,1)+.2)

for (i in plotitems){
  width<-max(datn[i,])-min(datn[i,])
  min <- min(datn[i,])-.1*width
  max <-  max(datn[i,])+.1*width
  #min <- min(dat[i,])-width
  #max <-  max(dat[i,])+width
  y <- seq(min,max,(max-min)/45)
  prob<-0*y
  
  for (l in 1:length(y))prob[l]<- {respfctd(theta - diffunct(parmatrest[i,1],parmatrest[i,2],y[l]))*
      derdiffunct(parmatrest[i,1],parmatrest[i,2],y[l])
  } 
  sum(prob)*(y[2]-y[1])
  #plot(y,prob,cex.lab=2.0, main =names(deca)[indsel[i]],xlab ="y", type="b", lwd=3,  ylab="Density", pch =1)#pcdum[1])
  
  meandat
  sddat
  ymod<-y*sddat[i]+meandat[i]
  plot(ymod,prob,cex.lab=2.0, main =names(deca)[indsel[i]],xlab ="y", type="b",
       lwd=3,  ylab="Density", pch =1)#pcdum[1])
}



##########################
####  two-dimensional


##  with prefit start
start2 <- list(
  theta = matrix(0, nrow = nrow(dat), ncol = 2),
  
  alpha = matrix(0, nrow = ncol(dat), ncol = 2),
  
  #delta0 = colMeans(dat),
  delta0=ThrNV$parmatrix[,1],
  #lambda = log(1 / apply(dat, 2, sd))
  lambda = log(ThrNV$parmatrix[,2]) 
  #lambda = rep(-2,dim(dat)[2])
)

fit2 <- fit.thresholdFlex(
  data = dat,
  F = F,
  f = f,
  g = g,
  gprime = gprime,
  q = 2,
  
  dalpha = dalpha.normal,
  alpha.par = list(
    mu = c(0, 0),
    Sigma = diag(2)
  ),
  mu.delta0 = 0,
  sd.delta0 = 5,
  mu.lambda = 0,
  sd.lambda = 2,
  #niter = 10000,
  #burnin = 2000,
  start=start2,
  niter = 6000,
  burnin = 400,
  thin = 5
)

cbind(fit2$estimates$delta0,fit2$estimates$delta)
ThrNV$parmatrix
#fit2.5000<-fit2
#fit2<-fit2.5000
#fit2.8000<-fit2
#fit2<-fit2.8000

####  three-dimensional
start3 <- list(
  theta = matrix(0, nrow = nrow(dat), ncol = 3),
  
  alpha = matrix(0, nrow = ncol(dat), ncol = 3),
  
  #delta0 = colMeans(dat),
  delta0=ThrNV$parmatrix[,1],
  #lambda = log(1 / apply(dat, 2, sd))
  lambda = log(ThrNV$parmatrix[,2])
)

fit3 <- fit.thresholdFlex(
  data = dat,
  F = F,
  f = f,
  g = g,
  gprime = gprime,
  q = 3,
  
  dalpha = dalpha.normal,
  alpha.par = list(
    mu = c(0, 0,0),
    Sigma = diag(3)
  ),
  mu.delta0 = 0,
  sd.delta0 = 5,
  mu.lambda = 0,
  sd.lambda = 2,
  start = start3,
  #niter = 10000,
  #burnin = 2000,
  niter = 5000,
  burnin = 500,
  thin = 5
)
#fit3.3000<-fit3

#######################

#loglik
fit1used<-fit2
ll.fit <- loglik.threshold(
  data = dat,
  theta = fit1used$estimates$theta,
  alpha = fit1used$estimates$alpha,
  delta0 = fit1used$estimates$delta0,
  lambda = fit1used$estimates$lambda,
  f = dnorm,  g = g,  gprime = gprime)
ll.fit

fit1used<-fit2
parmatrest <-cbind(fit1used$estimates$delta0,fit1used$estimates$delta) ### for plots
theta <-mean(fit1used$estimates$theta) ### for plots
thetavec<-rep(theta, dim(dat)[2])  ### for plots


### memory
#fitallordrun2Gomp<-fit2  ##two dim fit Gompertz (F)
#fitallordrun2Gumb<-fit2  ##two dim fit Gumbel (F)  
#fitallordrun3Gumb<-fit3  ##three dim fit Gumbel (F) 
#fitallordrun4Gumb<-fit4  ##  dim fit Gumbel (F) 
#fit2<-fitallordrun2Gumb
##########

fit1$estimates$alpha
fit2$estimates$alpha
fit2$estimates$delta0
fit2$estimates$delta

cbind(fit1$estimates$delta0,fit1$estimates$delta)
cbind(fit2$estimates$delta0,fit2$estimates$delta)
ThrNV$parmatrix



item<-4
dim1<-fit2$estimates$alpha[item]*fit2$estimates$theta[,1]
dim2<-fit2$estimates$alpha[item]*fit2$estimates$theta[,2]
plot(dim1, dim2, xlim = c(-2, 2),
     ylim = c(0, 4))
max(dim1)-min(dim1)
max(dim2)-min(dim2)


plot(density(dim1),
     lwd = 3,
     xlab = "x",
     ylab = "Density",xlim = c(-4, 3.5))

lines(density(dim2),
      lwd = 3,
      lty = 2)



### rotation
#library(psych)


#quartimax
library(GPArotation)

#choose fit that is used
fitused<-fit2
#fitused<-fit3

A <- as.matrix(fitused$estimates$alpha)
Theta <- as.matrix(fitused$estimates$theta)

plot(density(Theta[,2]),
     lwd = 3,
     xlab = "x",
     ylab = "Density")
lines(density(Theta[,1]))

thetavec<-colMeans(Theta%*%t(A))  ### for plot
parmatrest <-cbind(fit2$estimates$delta0,fit2$estimates$delta) ### for plots


#A <- as.matrix(fit3$estimates$alpha)
#Theta <- as.matrix(fit3$estimates$theta)

rot <- GPArotation::quartimax(A)
T <- rot$Th
A.rot <- A %*% T
Theta.rot <- Theta %*% T

rownames(A.rot) <- names(as.data.frame(dat))

A.rot
round(A.rot, 3)
#A.rot6000<-A.rot
#A.rot8000<-A.rot

##transformation
#A.rot<--A.rot
#A.rot<--A.rot[,c(2,1)]
#A.rot<-A.rot[,c(2,1)]
#A.rot<-A.rot[,c(2,1,3)]
A.pl<--A.rot[,c(2,1,3)]
A.pl<--A.rot[,c(2,1)]
round(A.pl, 3)


### two dimensional
item<-2
dim1<-A.rot[item,1]*Theta.rot[,1]
dim2<-A.rot[item,2]*Theta.rot[,2]

cor(dim2,dat[,item],method = "spearman")
plot(density(dim1),
     lwd = 3,
     xlab = "x",
     ylab = "Density",xlim = c(-4, 6))

lines(density(dim2),
      lwd = 3,
      lty = 2)
plot(density(dim2),
     lwd = 3,
     xlab = "x",
     ylab = "Density")


plot(dim1, dim2)

#plot(dim1, dim1+dim2,xlab = expression(theta[p1]^T %*% alpha),ylab = 
#       expression(theta[p1]^T %*% alpha+theta[p2]^T %*% alpha))
plot(dim1, dim1+dim2,xlab = "Dim 1",ylab = 
       "Dim 1 + Dim 2",cex.axis=1.4,cex.lab = 1.4,lwd=3, pch =1,main=rownames(A.rot)[item],cex.main=1.8)
lines(c(-20,20),c(-20,20),type="l")

### correl A.rot untransformed!
cormm<-matrix(0,dim(dat)[2],2)
for (item in 1:dim(dat)[2]){
  dim1<-A.rot[item,1]*Theta.rot[,1]
  dim2<-A.rot[item,2]*Theta.rot[,2]
  #print(cor(dim2,dat[,item],method = "spearman"))
  #print(cor(dim1,dat[,item]))
  plot(dim1, dat[,item],xlab = "Dim 1",ylab = " obs",cex.axis=1.4,cex.lab = 1.4,lwd=3, pch =1)
  plot(dim2, dat[,item],xlab = "Dim 2",ylab = " obs",cex.axis=1.4,cex.lab = 1.4,lwd=3, pch =1)
  cormm[item,1]<-cor(dim1,dat[,item])
  cormm [item,2 ]<-cor(dim2,dat[,item])
}
rownames(cormm) <- names(as.data.frame(dat))
round(cormm, 3)


### three dimensional

item<-9
#dim1<-A.rot[item,1]*fit2$estimates$theta[,1]
#dim2<-A.rot[item,2]*fit2$estimates$theta[,2]
dim1<-A.rot[item,1]*Theta.rot[,1]
dim2<-A.rot[item,2]*Theta.rot[,2]
dim3<-A.rot[item,3]*Theta.rot[,3]
dimtot<-cbind(dim1,dim2,dim3)

#library(GGally)
ggpairs(dimtot[, 1:3])

### correl A.rot untransformed!
cormm<-matrix(0,dim(dat)[2],3)
for (item in 1:dim(dat)[2]){
  dim1<-A.rot[item,1]*Theta.rot[,1]
  dim2<-A.rot[item,2]*Theta.rot[,2]
  dim3<-A.rot[item,3]*Theta.rot[,3]
  #plot(dim1, dat[,item],xlab = "Dim 1",ylab = " obs",cex.axis=1.4,cex.lab = 1.4,lwd=3, pch =1)
  #plot(dim2, dat[,item],xlab = "Dim 2",ylab = " obs",cex.axis=1.4,cex.lab = 1.4,lwd=3, pch =1)
  cormm[item,1]<-cor(dim1,dat[,item])
  cormm [item,2 ]<-cor(dim2,dat[,item])
  cormm [item,3 ]<-cor(dim3,dat[,item])
}
cormm

cormmpr<-cormm[,c(2,1,3)]
rownames(cormmpr) <- names(as.data.frame(dat))
round(cormmpr, 3)


## varimax
rot <- varimax(fit2$estimates$alpha)
alpha.rot <- as.matrix(rot$loadings)
T <- rot$rotmat

theta.rot <- fit2$estimates$theta %*% T

alpha.rot
theta.rot
-alpha.rot
names(dat)


############### plots
fitused<-fit1
#fitused<-fit1o
parmatrest <-cbind(fitused$estimates$delta0,fitused$estimates$delta) ### for plots

#parmatrest <-ThrNV$parmatrix[,1:2]

#datt<-t(dat)
datn<-decaselecttr 


pcdum <-c ('1','2','3','4','5','6','7','8','9')

theta <-0 # specify theta

plotitems<-1:dim(datn)[1]  ## all items

par(cex.lab=1,cex.axis=1.8,cex.main=2.20,pch=16,cex=1.0)
par(mar=c(7,6,6,1)+.2)

for (i in plotitems){
  width<-max(datn[i,])-min(datn[i,])
  min <- min(datn[i,])-.2*width
  max <-  max(datn[i,])+.2*width
  #min <- min(dat[i,])-width
  #max <-  max(dat[i,])+width
  y <- seq(min,max,(max-min)/45)
  prob<-0*y
  
  for (l in 1:length(y))prob[l]<- {respfctd(theta - diffunct(parmatrest[i,1],parmatrest[i,2],y[l]))*
      derdiffunct(parmatrest[i,1],parmatrest[i,2],y[l])
  } 
  sum(prob)*(y[2]-y[1])
  #lines(y,prob,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b", pch =pcdum[i])
  plot(y,prob,cex.lab=2.0, main =names(deca)[indsel[i]],xlab ="y", type="b",
       lwd=3,  ylab="Density", pch =1)#pcdum[1])
}


### delta function plot


theta <-0 # specify theta

plotitems<-1:dim(datn)[1]  ## all items

par(cex.lab=1,cex.axis=1.8,cex.main=2.20,pch=16,cex=1.0)
par(mar=c(7,6,6,1)+.2)

for (i in plotitems){
  width<-max(datn[i,])-min(datn[i,])
  min <- min(datn[i,])-.2*width
  max <-  max(datn[i,])+.2*width
  #min <- min(dat[i,])-width
  #max <-  max(dat[i,])+width
  y <- seq(min,max,(max-min)/45)
  delta<-0*y
  
  for (l in 1:length(y))delta[l]<- {diffunct(parmatrest[i,1],parmatrest[i,2],y[l])  } 
  sum(prob)*(y[2]-y[1])
  #lines(y,prob,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b", pch =pcdum[i])
  plot(y,delta,cex.lab=2.0, main =names(deca)[indsel[i]],xlab ="y", type="b",
       lwd=3,  ylab="Delta function", pch =1)#pcdum[1])
}






##########################
### factor analysis

library(psych)

fit <- fa(dat,nfactors = 2,rotate = "quartimax",fm = "ml")
fit <- fa(dat,nfactors = 3,rotate = "quartimax",fm = "ml")
fit <- fa(dat,nfactors = 4,rotate = "quartimax",fm = "ml")
fit <- fa(dat,nfactors = 5,rotate = "quartimax",fm = "ml")

print(fit$loadings, cutoff = 0.2)
print(fit$loadings, cutoff = 0.0)

## separate rotation , same result

fitm <- fa(dat, nfactors = 2, rotate = "none", fm = "ml")

library(GPArotation)
rot <- quartimax(as.matrix(fitm$loadings))

rot$loadings 


