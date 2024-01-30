
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ Distributions to simulate from ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DISTRIBUTIONS
#   0 - Pareto
#   1 - Gamma
#   2 - Weibull
#   3 - Power
#   4 - Lognormal
#   5 - Dhillon's
#   6 - Chen's
#   7 - Linear Failure Rate (LFR)
#   8 - Extreme value
#   9 - Half normal
#  10 - Beta1 (2nd shape parameter=1)
#  11 - Beta2 (2nd shape parameter=2)
#  12 - Exponential Power (EP)
#  13 - Exponential geometric (EG)
#  14 - Exponential logarithmic (EL)
#  15 - Exponential Nadarajah Haghighi (ENH1) (lambda=1,alpha=0.5)
#  16 - Exponential Nadarajah Haghighi (ENH2) (lambda=1,alpha=2)
#  17 - Beta-exponential (Beta-exp)

if (FALSE){
  MC   = 1e4
  parm = 2
  
  X = rgamma(MC,parm)+1
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  f = 1/gamma(parm)*(q-1)^(parm-1)*exp(-(q-1))
  plot(density(X))
  lines(q,f,col='red')
}

if (FALSE){
  MC   = 1e4
  parm = 2
  
  X = rweibull(MC,parm)+1
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  f = parm*(q-1)^(parm-1)*exp(-(q-1)^parm)
  plot(density(X))
  lines(q,f,col='red')
}

# x <- rgamma(n,parm)       # gamma
# x <- rweibull(n,parm)     # Weibull
# power:
rpower <- function(n,parm){
  X = runif(n)^parm
  return(X)
}

# lognormal:
if (FALSE){
  MC   = 1e4
  parm = 1
  
  X = rlnorm(MC,0,parm)+1
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  f = exp(-0.5*((log(q-1))/parm)^2)/(parm*(q-1)*sqrt(2*pi))
  plot(density(X),xlim=c(0,10))
  lines(q,f,col='red')
}

# Half normal
if (FALSE){
  MC   = 1e4
  parm = 1
  
  X = abs(rnorm(MC))+1
  X = sort(X)
  Q = 4
  q = seq(1.01,Q,(Q-1.01)/999)
  f = sqrt(2/pi)*exp(-0.5*(q-1)^2)
  plot(density(X),xlim=c(0,10))
  lines(q,f,col='red')
}

# Linear failure rate:
rLFR <- function(n,parm){
  U = runif(n)
  X = (-1+sqrt(1-2*parm*log(1-U)))/parm
  return(X)
}
if (FALSE){
  MC   = 1e4
  parm = 2
  
  X = rLFR(MC,parm)+1
  X = sort(X)
  Q = 4
  q = seq(1.01,Q,(Q-1.01)/999)
  f = (1+parm*(q-1))*exp(-(q-1)-0.5*parm*(q-1)^2)
  plot(density(X),xlim=c(0,Q))
  lines(q,f,col='red')
}

# Beta-exponential:
rBE <- function(n,parm){
  B = rbeta(n,parm,1)
  X = -log(1-B)
  return(X)
}
if (FALSE){
  MC   = 1e4
  parm = 2
  
  X = rBE(MC,parm)+1
  X = sort(X)
  Q = 7
  q = seq(1.01,Q,(Q-1.01)/999)
  f = parm*exp(-(q-1))*(1-exp(-(q-1)))^(parm-1)
  plot(density(X),xlim=c(0,Q))
  lines(q,f,col='red')
}

# Inverse beta:
dInvBeta <- function(x,parm){
  out = (1+parm)/x^2*(1-x^-1)^parm
  return(out)
}
pInvBeta <- function(x,parm){
  out = (1-x^-1)^(parm+1)
  return(out)
}
rInvBeta <- function(n,parm){
  U = runif(n)
  X = 1/(1-U^(1/(parm+1)))
  return(X)
}

if (FALSE){
  MC   = 1e4
  parm = 0.5
  
  X = rInvBeta(MC,parm)
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  #f = dInvBeta(q,parm)
  #F = cumsum(f)*(q[2]-q[1])
  F = pInvBeta(q,parm)
  plot(q,F,type='l',col='red',ylim=c(0,1))
  Fn = rep(0,length(q))
  for (j in 1:length(q)){
    Fn[j] = mean(X<=q[j])
  }
  lines(q,Fn,col='blue')
}

# Tilted Pareto:
dTiltedPareto <- function(x,parm){
  out = (1+parm)/(x+parm)^2
  return(out)
}
rTiltedPareto <- function(n,parm){
  out = (1+parm)/(1-runif(n))-parm
  return(out)
}
if (FALSE){
  MC   = 1e4
  parm = 0.4
  
  X = rTiltedPareto(MC,parm)
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  f = dTiltedPareto(q,parm)
  F = cumsum(f)*(q[2]-q[1])
  #F = pInvBeta(q,parm)
  plot(q,F,type='l',col='red',ylim=c(0,1))
  Fn = rep(0,length(q))
  for (j in 1:length(q)){
    Fn[j] = mean(X<=q[j])
  }
  lines(q,Fn,col='blue')
}

# Dhillon's:
rDhillon <- function(n,parm){
  X = exp((-log(1-runif(n)))^(1/(parm+1)))-1
  return(X)
}
if (FALSE){
  MC   = 1e4
  parm = 0.4
  
  X = rDhillon(MC,parm)+1
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  f = (parm+1)/(q)*exp(-(log(q))^(parm+1))*(log(q))^parm
  F = cumsum(f)*(q[2]-q[1])
  #F = pInvBeta(q,parm)
  plot(q,F,type='l',col='red',ylim=c(0,1))
  Fn = rep(0,length(q))
  for (j in 1:length(q)){
    Fn[j] = mean(X<=q[j])
  }
  lines(q,Fn,col='blue')
}

# Chen's:
rChen <- function(n,parm){
  X = log(1-0.5*log(1-runif(n)))^(1/parm)
  return(X)
}

# Extreme Value:
rEV <- function(n,parm){
  X = log(1-parm*log(1-runif(n)))
  return(X)
}
# Half normal:
rHN <- function(n,parm){
  X = abs(rnorm(n,0,parm))
  return(X)
}
# Beta 1:
rbeta1 <- function(n,parm){
  X = rbeta(n,parm,1)
  return(X)
}
# Beta 2:
rbeta2 <- function(n,parm){
  X = rbeta(n,parm,2)
  return(X)
}
# Exponential power:
rEP <- function(n,parm){
  U = runif(n)
  X = (log(1-log(1-U)))^(1/parm)
  return(X)
}
# Exponential geometric:
rEG <- function(n,parm){
  U    = runif(n)
  beta = 1
  X    = -1/beta*log((1-U)/(1-U*parm)) #Only use beta = 1
  return(X)
}
# Exponential logarithmic:
rEL <- function(n,parm){
  U = runif(n)
  X = log((1-parm)/(1-parm^(1-U)))
  return(X)
}
# Exponential Nadarajah Haghighi (1):
rENH1 <- function(n,parm){
  U = runif(n)
  X = ((1-log(1-U^(1/parm)))^2-1)
  return(X)
}
# Exponential Nadarajah Haghighi (2):
rENH2 <- function(n,parm){
  U = runif(n)
  X = ((1-log(1-U^(1/parm)))^(0.5)-1)
  return(X)
}

# Benini
rBenini <- function(n,parm){
  U = runif(n)
  X = exp((-1+sqrt(1-4*parm*log(1-U)))/(2*parm))
  return(X)
}

dBenini <- function(x,parm){
  T1  = exp(-parm*(log(x))^2)/x^2
  T2  = 1+2*parm*log(x)
  out = T1*T2
  return(out)
}

# Lomax
rLomax <- function(n,parm){
  eps <- parm[1]
  sig <- parm[2]
  U = runif(n)
  X = sig*((1-U)^(-eps)-1)/eps
  return(X)
}



# Mixture: Pareto - log-normal
# rpareto.lnorm <- function(n,mybet=3,mymu=-1,prop){
rpareto.lnorm <- function(n,mybet=3,mysig=2,prop){
  indx1 <- runif(n) > prop
  n1    <- sum(indx1)
  mymu  <- log(mybet/(mybet - 1) -1) - mysig^2/2
  # mysig <- sqrt(2*log(mybet/(mybet -1)-1) - 2*mymu)
  x1    <- rpareto(n1,mybet)
  x2    <- rlnorm(n - n1,mymu,mysig) + 1
  X     <- c(x1,x2)[c(which(indx1),which(!indx1))]
  return(X)
}

# Mixture: Pareto - half-normal
rpareto.HN <- function(n,mybet=3,prop){
  indx1 <- runif(n) > prop
  n1    <- sum(indx1)
  mysig <- sqrt(pi/2)*(mybet/(mybet -1)-1)
  x1    <- rpareto(n1,mybet)
  x2    <- rHN(n - n1,mysig) + 1
  X     <- c(x1,x2)[c(which(indx1),which(!indx1))]
  return(X)
}

#Mixture: Pareto - Weibull
rpareto.weibull <- function(n,mybet=3,myshap=0.5,prop){
  indx1 <- runif(n) > prop
  n1    <- sum(indx1)
  myscal<- (mybet/(mybet-1) - 1)/gamma(1 + 1/myshap)
  x1    <- rpareto(n1,mybet)
  x2    <- rweibull(n - n1,scale=myscal,shape=myshap) + 1
  X     <- c(x1,x2)[c(which(indx1),which(!indx1))]
  return(X)  
}



#Mixture: Pareto - Exponential
rpareto.exp <- function(n,mybet=3,prop){
  indx1 <- runif(n) > prop
  n1    <- sum(indx1)
  mylam <- 1/(mybet/(mybet-1) - 1)
  # expec <- mybet/(mybet - 1)
  # x1    <- rpareto(n,mybet)
  # x2    <- rexp(n,mylam) + 1
  # mean(x1)
  # mean(x2)
  # expec
  # mybet
  x1    <- rpareto(n1,mybet)
  x2    <- rexp(n - n1,mylam) + 1
  X     <- c(x1,x2)[c(which(indx1),which(!indx1))]
  return(X)  
}

# par(mfcol=c(2,2))
# n<-1e7
# ppp <- 0.5
# xxx<-sort(rpareto.exp(n,3,ppp))[1:(n*0.95)]
# plot(density(xxx),xlim=c(1,2))
# abline(v=1,lty=2)
# 
# xxx<-sort(rpareto.HN(n,3,ppp))[1:(n*0.95)]
# plot(density(xxx),xlim=c(1,2))
# abline(v=1,lty=2)
# 
# # xxx<-sort(rpareto.lnorm(n,3,mymu=-1,ppp))[1:(n*0.90)]
# xxx<-sort(rpareto.lnorm(n,3,mysig=2,ppp))[1:(n*0.95)]
# plot(density(xxx),xlim=c(1,2))
# abline(v=1,lty=2)
# 
# xxx<-sort(rpareto.weibull(n,3,0.5,ppp))[1:(n*0.95)]
# plot(density(xxx),xlim=c(1,2))
# abline(v=1,lty=2)



if (FALSE){
  MC   = 1e4
  parm = 0.2
  
  X = rBenini(MC,parm)
  X = sort(X)
  Q = 10
  q = seq(1.01,Q,(Q-1.01)/999)
  f = dBenini(q,parm)
  F = cumsum(f)*(q[2]-q[1])
  #F = pInvBeta(q,parm)
  plot(q,F,type='l',col='red',ylim=c(0,1))
  Fn = rep(0,length(q))
  for (j in 1:length(q)){
    Fn[j] = mean(X<=q[j])
  }
  lines(q,Fn,col='blue')
}



#rpareto(n,beta)
#rpower(n,parm)
#rLN(n,parm)
#rDhillon(n,parm)
#rChen(n,parm)
#rLFR(n,parm)
#rEV(n,parm)
#rHN(n,parm)
#rbeta1(n,parm)
#rbeta2(n,parm)
#rEP(n,parm)
#rEG(n,parm)
#rEL(n,parm)
#rENH1(n,parm)
#rENH2(n,parm)
#rBE(n,parm)
#rInvBeta(n,parm)
#rTiltedPareto(n,parm)
#rBenini(n,parm)

SimFromDist <- function(DistNum,n){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Pareto
  if (DistNum == 1){ X = rpareto(n,1) ; nm="P(1)" }
  else if (DistNum == 2){ X = rpareto(n,2) ; nm="P(2)" }
  else if (DistNum == 3){ X = rpareto(n,5) ; nm="P(5)"  }
  else if (DistNum == 4){ X = rpareto(n,10); nm="P(10)"  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Gamma
  else if (DistNum == 5){ X = rgamma(n,0.5)+1 ; nm="G(0.5)"}
  else if (DistNum == 6){ X = rgamma(n,0.8)+1 ; nm="G(0.8)"}
  else if (DistNum == 7){ X = rgamma(n,1)+1   ; nm="G(1.0)"}
  else if (DistNum == 8){ X = rgamma(n,1.2)+1 ; nm="G(1.2)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Weibull
  else if (DistNum ==  9){ X = rweibull(n,0.5)+1 ; nm="W(0.5)" }
  else if (DistNum == 10){ X = rweibull(n,0.8)+1 ; nm="W(0.8)" }
  else if (DistNum == 11){ X = rweibull(n,1.2)+1 ; nm="W(1.2)" }
  else if (DistNum == 12){ X = rweibull(n,1.5)+1 ; nm="W(1.5)" }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Lognormal
  else if (DistNum == 13){ X = rlnorm(n,0,1)+1  ; nm = "LN(0,1.0)" }
  else if (DistNum == 14){ X = rlnorm(n,0,1.2)+1; nm = "LN(0,1.2)" }
  else if (DistNum == 15){ X = rlnorm(n,0,1.5)+1; nm = "LN(0,1.5)" }
  else if (DistNum == 16){ X = rlnorm(n,0,2.5)+1; nm = "LN(0,2.5)" }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Half normal
  else if (DistNum == 17){ X = abs(rnorm(n,0,0.5))+1 ; nm = "HN(0,0.5)"}
  else if (DistNum == 18){ X = abs(rnorm(n,0,0.8))+1 ; nm = "HN(0,0.8)"}
  else if (DistNum == 19){ X = abs(rnorm(n,0,1))+1   ; nm = "HN(0,1.0)"}
  else if (DistNum == 20){ X = abs(rnorm(n,0,1.2))+1 ; nm = "HN(0,1.2)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Linear failure rate
  else if (DistNum == 21){ X = rLFR(n,0.2)+1 ; nm = "LFR(0.2)"}
  else if (DistNum == 22){ X = rLFR(n,0.5)+1 ; nm = "LFR(0.5)"}
  else if (DistNum == 23){ X = rLFR(n,0.8)+1 ; nm = "LFR(0.8)"}
  else if (DistNum == 24){ X = rLFR(n,1)+1   ; nm = "LFR(1.0)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Beta-exponential
  else if (DistNum == 25){ X = rBE(n,0.5)+1 ; nm = "BE(0.5)"}
  else if (DistNum == 26){ X = rBE(n,0.8)+1 ; nm = "BE(0.8)"}
  else if (DistNum == 27){ X = rBE(n,1)+1   ; nm = "BE(1.0)"}
  else if (DistNum == 28){ X = rBE(n,1.5)+1 ; nm = "BE(1.5)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Tilted Pareto
  else if (DistNum == 29){ X = rTiltedPareto(n,0.5) ; nm = "TP(0.5)"}
  else if (DistNum == 30){ X = rTiltedPareto(n,1)   ; nm = "TP(1.0)"}
  else if (DistNum == 31){ X = rTiltedPareto(n,2)   ; nm = "TP(2.0)"}
  else if (DistNum == 32){ X = rTiltedPareto(n,3)   ; nm = "TP(3.0)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Dhillon
  else if (DistNum == 33){ X = rDhillon(n,0.2)+1 ; nm = "DH(0.2)"}
  else if (DistNum == 34){ X = rDhillon(n,0.4)+1 ; nm = "DH(0.4)"}
  else if (DistNum == 35){ X = rDhillon(n,0.6)+1 ; nm = "DH(0.6)"}
  else if (DistNum == 36){ X = rDhillon(n,0.8)+1 ; nm = "DH(0.8)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Pareto - Exponential
  
  else if (DistNum == 37){ 
    prop  <- 0.1
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.1"
  } else if (DistNum == 38){ 
    prop  <- 0.2
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.2"
  } else if (DistNum == 39){ 
    prop  <- 0.3
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.3"
  } else if (DistNum == 40){ 
    prop  <- 0.4
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.4"
  }else if (DistNum == 41){ 
    prop  <- 0.5
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.5"
  } else if (DistNum == 42){ 
    prop  <- 0.6
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.6"
  } else if (DistNum == 43){ 
    prop  <- 0.7
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.7"
  } else if (DistNum == 44){ 
    prop  <- 0.8
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.8"
  } else if (DistNum == 45){ 
    prop  <- 0.9
    mybet <- 3
    X <- rpareto.exp(n,mybet,prop)
    nm = "P(3) & E(lambda), p=0.9"
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  # Pareto - Log-Normal
  
  else if (DistNum == 46){ 
    prop  <- 0.1
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.1"
  } else if (DistNum == 47){ 
    prop  <- 0.2
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.2"
  } else if (DistNum == 48){ 
    prop  <- 0.3
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.3"
  } else if (DistNum == 49){ 
    prop  <- 0.4
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.4"
  }else if (DistNum == 50){ 
    prop  <- 0.5
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.5"
  } else if (DistNum == 51){ 
    prop  <- 0.6
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.6"
  } else if (DistNum == 52){ 
    prop  <- 0.7
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.7"
  } else if (DistNum == 53){ 
    prop  <- 0.8
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.8"
  } else if (DistNum == 54){ 
    prop  <- 0.9
    mybet <- 3
    mysig <- 2
    X <- rpareto.lnorm(n,mybet,mysig,prop)
    nm = "P(3) & LN(-2.69,2), p=0.9"
  } 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  # Pareto - HalfNormal  
  
  else if (DistNum == 55){ 
    prop  <- 0.1
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.1"
  } else if (DistNum == 56){ 
    prop  <- 0.2
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.2"
  } else if (DistNum == 57){ 
    prop  <- 0.3
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.3"
  } else if (DistNum == 58){ 
    prop  <- 0.4
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.4"
  }else if (DistNum == 59){ 
    prop  <- 0.5
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.5"
  } else if (DistNum == 60){ 
    prop  <- 0.6
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.6"
  } else if (DistNum == 61){ 
    prop  <- 0.7
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.7"
  } else if (DistNum == 62){ 
    prop  <- 0.8
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.8"
  } else if (DistNum == 63){ 
    prop  <- 0.9
    mybet <- 3
    X <- rpareto.HN(n,mybet,prop)
    nm = "P(3) & HN(sigma), p=0.9"
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Pareto - Weibull
  
  else if (DistNum == 64){ 
    prop <- 0.1
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.1"
  } else if (DistNum == 65){ 
    prop  <- 0.2
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.2"
  } else if (DistNum == 66){ 
    prop  <- 0.3
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.3"
  } else if (DistNum == 67){ 
    prop  <- 0.4
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.4"
  }else if (DistNum == 68){ 
    prop  <- 0.5
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.5"
  } else if (DistNum == 69){ 
    prop  <- 0.6
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.6"
  } else if (DistNum == 70){ 
    prop  <- 0.7
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.7"
  } else if (DistNum == 71){ 
    prop  <- 0.8
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.8"
  } else if (DistNum == 72){ 
    prop  <- 0.9
    mybet <- 3
    myshap<- 0.5
    X <- rpareto.weibull(n,mybet,myshap,prop)
    nm = "P(3) & W(0.5,0.25), p=0.9"
  } 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Log Weibull
  else if (DistNum == 73){ X = exp((-log(1-runif(n)))^(1/(0.1+1))) ; nm = "LW(0.1)"}
  else if (DistNum == 74){ X = exp((-log(1-runif(n)))^(1/(0.3+1))) ; nm = "LW(0.3)"}
  else if (DistNum == 75){ X = exp((-log(1-runif(n)))^(1/(0.5+1))) ; nm = "LW(0.5)"}
  else if (DistNum == 76){ X = exp((-log(1-runif(n)))^(1/(0.7+1))) ; nm = "LW(0.7)"}
  else if (DistNum == 77){ X = exp((-log(1-runif(n)))^(1/(0.9+1))) ; nm = "LW(0.9)"}
  else if (DistNum == 78){ X = exp((-log(1-runif(n)))^(1/(1.1+1))) ; nm = "LW(1.1)"}
  else if (DistNum == 79){ X = exp((-log(1-runif(n)))^(1/(1.3+1))) ; nm = "LW(1.3)"}
  else if (DistNum == 80){ X = exp((-log(1-runif(n)))^(1/(1.5+1))) ; nm = "LW(1.5)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Inverse Beta
  else if (DistNum == 81){ X = 1/(1-runif(n)^(1/(0.1+1))) ; nm = "IB(0.1)"}
  else if (DistNum == 82){ X = 1/(1-runif(n)^(1/(0.3+1))) ; nm = "IB(0.3)"}
  else if (DistNum == 83){ X = 1/(1-runif(n)^(1/(0.5+1))) ; nm = "IB(0.5)"}
  else if (DistNum == 84){ X = 1/(1-runif(n)^(1/(0.7+1))) ; nm = "IB(0.7)"}
  else if (DistNum == 85){ X = 1/(1-runif(n)^(1/(0.9+1))) ; nm = "IB(0.9)"}
  else if (DistNum == 86){ X = 1/(1-runif(n)^(1/(1.1+1))) ; nm = "IB(1.1)"}
  else if (DistNum == 87){ X = 1/(1-runif(n)^(1/(1.3+1))) ; nm = "IB(1.3)"}
  else if (DistNum == 88){ X = 1/(1-runif(n)^(1/(1.5+1))) ; nm = "IB(1.5)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Benini
  else if (DistNum == 89){ X = exp((-1+sqrt(1-4*0.1*log(1-runif(n))))/(2*0.1)) ; nm = "BN(0.1)"}
  else if (DistNum == 90){ X = exp((-1+sqrt(1-4*0.3*log(1-runif(n))))/(2*0.3)) ; nm = "BN(0.3)"}
  else if (DistNum == 91){ X = exp((-1+sqrt(1-4*0.5*log(1-runif(n))))/(2*0.5)) ; nm = "BN(0.5)"}
  else if (DistNum == 92){ X = exp((-1+sqrt(1-4*0.7*log(1-runif(n))))/(2*0.7)) ; nm = "BN(0.7)"}
  else if (DistNum == 93){ X = exp((-1+sqrt(1-4*0.9*log(1-runif(n))))/(2*0.9)) ; nm = "BN(0.9)"}
  else if (DistNum == 94){ X = exp((-1+sqrt(1-4*1.1*log(1-runif(n))))/(2*1.1)) ; nm = "BN(1.1)"}
  else if (DistNum == 95){ X = exp((-1+sqrt(1-4*1.3*log(1-runif(n))))/(2*1.3)) ; nm = "BN(1.3)"}
  else if (DistNum == 96){ X = exp((-1+sqrt(1-4*1.5*log(1-runif(n))))/(2*1.5)) ; nm = "BN(1.5)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Gumbel/Extreme Value
  else if (DistNum == 97 ){ X = log(1 - 0.1*log(1 - runif(n))) ; nm = "EV(0.1)"}
  else if (DistNum == 98 ){ X = log(1 - 0.3*log(1 - runif(n))) ; nm = "EV(0.3)"}
  else if (DistNum == 99 ){ X = log(1 - 0.5*log(1 - runif(n))) ; nm = "EV(0.5)"}
  else if (DistNum == 100){ X = log(1 - 0.7*log(1 - runif(n))) ; nm = "EV(0.7)"}
  else if (DistNum == 101){ X = log(1 - 0.9*log(1 - runif(n))) ; nm = "EV(0.9)"}
  else if (DistNum == 102){ X = log(1 - 1.1*log(1 - runif(n))) ; nm = "EV(1.1)"}
  else if (DistNum == 103){ X = log(1 - 1.3*log(1 - runif(n))) ; nm = "EV(1.3)"}
  else if (DistNum == 104){ X = log(1 - 1.5*log(1 - runif(n))) ; nm = "EV(1.5)"}
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (DistNum == 105){ X = rLomax(n,c(0.1,1));nm = "Lomax(0.1,1)"}
  else if (DistNum == 106){ X = rpareto(n,100); nm="P(100)"  }
  else if (DistNum == 107){ X = rLomax(n,c(0.8,1));nm = "Lomax(0.8,1)"}
  else if (DistNum == 108){ X = rLomax(n,c(0.2,1));nm = "Lomax(0.2,1)"}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Gamma
  else if (DistNum == 109){ X = rgamma(n,0.5) ; nm="G(0.5)"}
  else if (DistNum == 110){ X = rgamma(n,0.8) ; nm="G(0.8)"}
  else if (DistNum == 111){ X = rgamma(n,1)   ; nm="G(1.0)"}
  else if (DistNum == 112){ X = rgamma(n,1.2) ; nm="G(1.2)"}
  
  ans <- list(X=X,nm=nm)
  return(ans)
}

#rTiltedPareto(n,parm)