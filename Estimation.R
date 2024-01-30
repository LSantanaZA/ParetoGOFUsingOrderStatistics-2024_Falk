
ppareto <- function(x,beta,sigma=1){
  out = 1-(x/sigma)^-beta
  return(out)
}

dpareto <- function(x,beta,sigma=1){
  out = beta*sigma^beta/(x^beta+1)
  return(out)
}

rpareto <- function(n,beta,sigma=1){
  X = sigma*(1-runif(n))^(-beta^-1)
  return(X)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ MLE sigma known ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMLE1 <- function(X){
  n   = length(X)
  out = n/sum(log(X))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ MLE sigma unknown ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMLE2 <- function(X){
  n      = length(X)
  sigmaH = min(X)
  betaH  = n/sum(log(X/sigmaH))
  out    = list("betaH"=betaH,"sigmaH"=sigmaH)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ MME sigma known ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMME1 <- function(X){
  out = mean(X)/(mean(X)-1)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ MME sigma unknown ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ParetoMME2 <- function(X){
  betaH  = (n*mean(X)-min(X))/(n*(mean(X)-min(X)))
  sigmaH = mean(X)*(betaH-1)/betaH
  out    = list("betaH"=betaH,"sigmaH"=sigmaH)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(FALSE){
  n     = 5000
  beta  = 3
  sigma = 2
  X     = sort(rpareto(n,beta,sigma))
  
  betaH = ParetoMLE1(X/sigma)
  Y     = (X/sigma)^betaH
  betaH
  ParetoMLE1(Y)
  
  est = ParetoMLE2(X)
  Y   = (X/est$sigmaH)^est$betaH
  est
  ParetoMLE2(Y)
  
  betaH = ParetoMME1(X/sigma)
  Y     = (X/sigma)^betaH
  betaH
  ParetoMME1(Y)
  
  est = ParetoMME2(X)
  Y   = (X/est$sigmaH)^est$betaH
  est
  ParetoMLE2(Y)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~