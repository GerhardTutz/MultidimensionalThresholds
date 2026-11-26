
### fears: multidimensional

#getwd() # check 
#setwd("C:\\Users\\tutz\\LRZ Sync+Share\\TuRaschThresholdsM2\\R3\\")


############ Fit data sets: first fixed difficulty functions + plots
#                          later Splines + plots

library(splines)
library(tidyr)
library(spatstat) 
library("foreign") 


## for start values:
source("Responsefunctions.R")
source("Difficultyfunctions.R")
source("ProgramsBsplinesPureResp.R")
source("Posteriorgen.R")
source("Posteriorgen.R")

## for Bayes
source("ProgrThresholdsBayesFlex.R")


##### fears  data data set I x P matrix 
dat <- read.table("fear", header=TRUE,sep= " ")
###################

I <- dim(dat)[1]
P <- dim(dat)[2]


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

#First run:
#Difficultyfunctions
#Responsefunctionsfunctions


#ThrNV1<-ThreshbasisRespCov(dat,basis="lin",indicatorvec,ord=4, numknotsstart=6,lambdpenin=0,lambdc=0,start=0, 
#                          startind = "no",der="der",cov0,hessiander="FALSE")

ThrNV<-ThreshIncl(respf,lin,dat,basis="lin",indicatorvec,ord, numknotsstart,lambdpenin=0,lambdc=0,start=0, 
                  startind = "no",der="der",cov0,hessiander="FALSE")
ThrNV$parmatrix 

### for Bayes transposed data 

dat<-t(dat)
summary(dat)

### specification

### now data 0,....
dat<-dat-1

# linear
#lin <- "lin"
#g <- rep(list(function(y) y), dim(dat)[2])
#gprime <- rep(list(function(y) 1), dim(dat)[2])


## logit
lin <- "logit"
at<--.2 ## lower threshold
bt<-6+.1   ## upper threshold

tlogit<-function(y){
  at<-1-.2 ## lower threshold
  bt<-7+.1
  res<-log((y-at)/(bt-y))
  res
  }
tlogit(7)


g <- rep(list(function(y) log((y-at)/(bt-y))), dim(dat)[2])
#g <- rep(list(function(y) tlogit), dim(dat)[2])
gprime <- rep(list(function(y) (bt-at)/((y-at)*(bt-y)) ), dim(dat)[2])
g[[i]](6)
gprime[[i]](6)

## only one
#respf<-"Gumbel"
#respf<-"Gompertz"
#respf<-"Gompertz"
respf<-"NV"

## NV
if(respf=="NV"){
    F <-function(x) pnorm(x) 
  f <- function(x) dnorm(x)
  #f <-function(x){r<-dnorm    return(r)  
  }

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



fit1 <- fit.threshold(
  data = dat,
  #F = pnorm,
  #f = dnorm,
  F = F,
  f = f,
  g = g,
  gprime = gprime,
  q = 1,
  
  dalpha = dalpha.normal,
  alpha.par = list(
    mu = c(0),
    Sigma = diag(1)
  ),
  start = start1,
  #niter = 10000,
  #burnin = 2000,
  niter = 100,
  burnin = 50,
  thin = 5
)

fit1$estimates$alpha
fit1$estimates$delta0
fit1$estimates$delta

fit1discr <-fit.thresholdFlex(data,
                              F = pnorm,
                              f = dnorm,
                              g = g,
                              gprime = NULL,
                              q = 1,
                              dalpha= dalpha.normal,
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
                              
                              discrete = TRUE,
                              verbose = TRUE) 



fit1discr$estimates$alpha
fit1discr$estimates$delta0
fit1discr$estimates$delta

parmatrest <-cbind(fit1$estimates$delta0,fit1$estimates$delta) ### for plots
theta <-mean(fit1$estimates$theta) ### for plots
thetavec<-rep(theta, dim(dat)[2])  ### for plots


####  two-dimensional
start2 <- list(
  theta = matrix(0, nrow = nrow(dat), ncol = 2),
  
  alpha = matrix(0, nrow = ncol(dat), ncol = 2),
  
  #delta0 = colMeans(dat),
  delta0=ThrNV$parmatrix[,1],
  #lambda = log(1 / apply(dat, 2, sd))
  lambda = log(ThrNV$parmatrix[,2])
)

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
  mu.delta0 = 0,
  sd.delta0 = 5,
  mu.lambda = 0,
  sd.lambda = 2,
  discrete = TRUE,
  niter = 1000,
  burnin = 500,
  start = start2  ##now
  )

fit2discr$estimates$alpha
fit2discr$estimates$delta0
fit2discr$estimates$delta

cbind(fit2discr$estimates$delta0,fit2discr$estimates$delta)
ThrNV$parmatrix


parmatrest <-cbind(fit2discr$estimates$delta0,fit2discr$estimates$delta) ### for plots
thetapl <-mean(fit2discr$estimates$theta) ### for plots
#thetavec<-rep(theta, dim(dat)[2])  ### for plots
#fitmem<-fit2discr 


####  three-dimensional
start3 <- list(
  theta = matrix(0, nrow = nrow(dat), ncol = 3),
  
  alpha = matrix(0, nrow = ncol(dat), ncol = 3),
  
  #delta0 = colMeans(dat),
  delta0=ThrNV$parmatrix[,1],
  #lambda = log(1 / apply(dat, 2, sd))
  lambda = log(ThrNV$parmatrix[,2])
)


fit3discrn <-  fit.thresholdFlex(
  data = dat,
  F = pnorm,
  f = dnorm,
  g = g,
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
  discrete = TRUE,
  niter = 5000,
  burnin = 500,
  start = start3 
  )
parmatrest <-cbind(fit3discrn$estimates$delta0,fit3discr$estimates$delta) ### for plots
thetapl <-mean(fit3discrn$estimates$theta) ### for plots
thetavec<-rep(theta, dim(dat)[2])  ### for plots


### rotation quartimax
#library(psych)
#library(GPArotation)
#fitused<-fit2discr
#fitused<-fit3discrn

A <- as.matrix(fitused$estimates$alpha)
Theta <- as.matrix(fitused$estimates$theta)

rot <- GPArotation::quartimax(A)
T <- rot$Th
A.rot <- A %*% T
Theta.rot <- Theta %*% T

rownames(A.rot) <- colnames(dat)
#rownames(A.rot) <- c("Refugee","Climate","Terrorism","Globalization","Nuclear Energy")
A.rot
round(A.rot, 3)


### modifications 
Amod<-cbind(A.rot[,3],-A.rot[,1],A.rot[,2])
#rownames(Amod) <- c("Refugee","Climate","Terrorism","Globalization","Nuclear Energy")
round(Amod, 3)

Amod<-cbind(-A.rot[,1],A.rot[,2])
#rownames(Amod) <- c("Refugee","Climate","Terrorism","Globalization","Nuclear Energy")
round(Amod, 3)

#Amodd <- Amod[c(1, 3, 4,2,5), ]
#round(Amodd, 3)


#### Plots with possible item specific latent trait

#parmatrest <- ThrVC$parmatrix ### or other fitted models

lin <- "logit"
#lin <- "lin"

#### densities discrete
Y <- as.matrix(dat)
ymax <- max(Y, na.rm = TRUE)
ymin <- min(Y, na.rm = TRUE)

#Yres<-Y-ymin+1  # values 1,2, 


 
theta <- Amod%*%c(0,0)  ## theta zero
theta <- Amod%*%c(.6,0)  ### first dimension dominant
#theta <- Amod%*%c(0,.6)
thetam<-rep(theta, dim(dat)[2])  # same

#thetam<-Amod%*%c(0,1.0)   ##modified
#thetam<-Amod%*%c(1,0,0)   ##modified

#pcdum <-c ('1','2','3','4','5','6')
pcdum <-c ('R','C','T','G','N','6')  ### for original A.rot!

ylims <- c(0,0.5)
min <- min(dat)
max <- max(dat)

y <- seq(min,max,1)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "logit"){delta <- parmatrest[1,1]+ parmatrest[1,2]*g[[1]](y)}

problargy <- pnorm(thetam[1] - delta,0,1)
#plot(y,problargy)

dens<-NULL
for (i in 1:(length(y)-1)){dens<-c(dens,problargy[i]-problargy[i+1])
print (dens)}
dens<-c(dens,problargy[length(y)])

sum(dens)

plot(y,dens,cex.axis=1.5,cex.lab=1.5,cex.main=1.6, main ="Discrete densities",xlab ="y",lty=2, type="b",lwd=1.9, ylab="", pch =pcdum[1], ylim=ylims,cex=1.6)

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "logit"){delta <-parmatrest[i,1]+ parmatrest[i,2]*g[[1]](y)}
  problargy <- pnorm(thetam[i] - delta,0,1)
  dens<-NULL
  for (ii in 1:(length(y)-1)){dens<-c(dens,problargy[ii]-problargy[ii+1])}
  dens<-c(dens,problargy[length(y)])
  print (sum(dens))
  lines(y,dens,cex.axis=1.5,cex.lab=1.5,cex.main=1.6,lty=2,lwd=1.9, type="b",pch =pcdum[i],cex=1.6)
}


################################
### correlations
#####################################

### correl 2-dim
cormm<-matrix(0,dim(dat)[2],2)
for (item in 1:dim(dat)[2]){
  dim1<-A.rot[item,1]*Theta.rot[,1]
  dim2<-A.rot[item,2]*Theta.rot[,2]
  #print(cor(dim2,dat[,item],method = "spearman"))
  #print(cor(dim1,dat[,item]))
  cormm[item,1]<-cor(dim1,dat[,item])
  cormm [item,2 ]<-cor(dim2,dat[,item])
}
rownames(cormm) <- c("Refugee","Climate","Terrorism","Globalization","Nuclear Energy")
round(cormm, 3)


### correl 3-dim
cormm<-matrix(0,dim(dat)[2],3)
for (item in 1:dim(dat)[2]){
  dim1<-A.rot[item,1]*Theta.rot[,1]
  dim2<-A.rot[item,2]*Theta.rot[,2]
  dim3<-A.rot[item,3]*Theta.rot[,3]
  #print(cor(dim1,dat[,item]))
  cormm[item,1]<-cor(dim1,dat[,item])
  cormm [item,2 ]<-cor(dim2,dat[,item])
  cormm [item,3 ]<-cor(dim3,dat[,item])
  }
rownames(cormm) <- c("Refugee","Climate","Terrorism","Globalization","Nuclear Energy")
round(cormm, 3)
round(cbind(cormm[,3],cormm[,1],cormm[,2]), 3)

library(psych)

fit2 <- fa(dat,nfactors = 2,rotate = "quartimax",fm = "ml")
fit2
fit3 <- fa(dat,nfactors = 3,rotate = "quartimax",fm = "ml")
fit3
fit4 <- fa(dat,nfactors = 4,rotate = "quartimax",fm = "ml")
fit4
fit <- fa(dat,nfactors = 5,rotate = "quartimax",fm = "ml")


###########################################
############################################

#### densities discrete with item-specific latent trait
Y <- as.matrix(dat)
ymax <- max(Y, na.rm = TRUE)
ymin <- min(Y, na.rm = TRUE)

#Yres<-Y-ymin+1  # values 1,2, 

### same order as in Amodd
parmatrest <- ThrVC$parmatrix[c(1, 3, 4,2,5), ]


thetam<-Amodd%*%c(0,1.5)

pcdum <-c ('1','2','3','4','5','6')
pcdum <-c ('R','C','T','G','NE','6')  ### for original A.rot!
pcdummod <-pcdum[c(1, 3, 4,2,5)]   ### for Amodd

ylims <- c(0,0.4)
min <- min(dat)
max <- max(dat)

y <- seq(min,max,1)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}

problargy <- pnorm(thetam[1] - delta,0,1)
dens<-problargy[ymin]-problargy[ymin+1]
for (i in (ymin+1):(ymax-1)){dens<-c(dens,problargy[i]-problargy[i+1])}
dens<-c(dens,problargy[ymax])
sum(dens)

plot(y,dens,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Densities",xlab ="y", type="b",lwd=1.2, ylab="", pch =pcdummod[1], ylim=ylims)

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  problargy <- pnorm(thetam[i] - delta,0,1)
  dens<-problargy[ymin]-problargy[ymin+1]
  for (ii in (ymin+1):(ymax-1)){dens<-c(dens,problargy[ii]-problargy[ii+1])}
  dens<-c(dens,problargy[ymax])
  
  lines(y,dens,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdummod[i],cex = 1)
}





### for data transformation if lin= inv
width<- max(dat)-min(dat)
norm<- 10  ### 10 means in (0.0833 0.91666)
mintr <- min(dat)- width/norm
maxtr <- max(dat)+ width/norm
datstd<- (dat -mintr)/(maxtr-mintr)
#####################################

###plot   delta fct 
#parmatrest <- ThrV$parmatrix ### 
min <- min(dat)
max <- max(dat)
ylims <- c(-5,2)
x<- seq(min,max,(max-min)/40)
pcdum <-c ('1','2','3','4','5','6')
#x<- seq(0,maxresp,1)
y<- 0*x
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*x
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+x)
if(lin == "inv"){xt <- (x -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
plot(x,delta,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Difficulty functions",xlab ="y", type="b",lwd=1.2,
     ylab="",pch =pcdum[1],cex = 1.2, ylim=ylims )

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*x
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+x)
  if(lin == "inv"){xt <- (x -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  lines(x,delta,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=1.2, type="b",pch =pcdum[i],cex = 1.2)
}

####PT function 
#parmatrest <- ThrV$parmatrix ### 
pcdum <-c ('1','2','3','4','5','6')

y <- seq(min,max,(max-min)/40)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
theta <- 0 
problargy <- pnorm(theta - delta,0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="P(Y>y)",xlab ="y", type="b",lwd=1.2,ylim=c(0,1),
     ylab="", pch =pcdum[1])

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  problargy <- pnorm(theta - delta,0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}

##### IC curves
#parmatrest <- Thr$parmatrix ### 
y=0 ### specified y-value

theta <- seq(-8,4,0.1)
numtheta <- length(theta)

if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
problargy <- pnorm(theta - delta,0,1)
plot(theta,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="IC function",xlab ="theta", type="b",
     lwd=2.0,ylim=c(0,1), pch =pcdum[1])

for (l in 2:I){
  if(lin == "lin") delta <- parmatrest[l,1]+ parmatrest[l,2]*y
  if(lin == "log") delta <- parmatrest[l,1]+ parmatrest[l,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[l,1]+ parmatrest[l,2]*log(xt/(1-xt))}
  problargy <- pnorm(theta - delta,0,1)
  
  lines(theta,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=2.0, type="b",pch =pcdum[l],cex = 1)
}


###### plots with alpha

parmatrest <- ThrCVar$parmatrix ### for continuous
alpha<- ThrCVar$alpha

parmatrest <- ThrDVar$parmatrix ### for discrete
alpha<- ThrDVar$alpha

pcdum <-c ('1','2','3','4','5','6')

y <- seq(min,max,(max-min)/40)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
theta <- 0 

problargy <- pnorm(alpha[1]*(theta - delta),0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="P(Y>y)",xlab ="y", type="b",lwd=1.2,ylim=c(0,1),
     ylab="", pch =pcdum[1])

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  problargy <- pnorm(alpha[1]*(theta - delta),0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}






