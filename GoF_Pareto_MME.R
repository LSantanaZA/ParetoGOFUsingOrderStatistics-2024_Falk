# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~ (0a) Testing for the Pareto distribution ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ Statistic using characterisation: X^(1/m) ~ X_(1) ~~~~~~~~
# ~~~~~~~~~~~~~~~~~ Laplace weight function version ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#numeric integration version
Tnm1_int <- function(X,m,a,betaH){
  X = sort(X)
  n = length(X)
  ii <- 0+1i
  inn <- function(t)    (Mod((1/sqrt(n))*sum(exp(-ii*t*X^(1/m))/m - X^(-betaH*(m-1))*exp(-ii*t*X))))^2 * exp(-a*abs(t))
  inn <- Vectorize(inn)
  # elliptic::myintegrate(inn,lower=0,upper=1)
  tmp <- try(integrate(inn,lower=-Inf,upper=Inf)$value)
  if(is.character(tmp)){
    out <- Inf
  } else {
    out <- tmp
  }
  return(out)
}

#double loop version
Tnm1_alt <- function(X,m,a,betaH){
  X = sort(X)
  n = length(X)
  matr = matrix(0,n,n)
  for (j in 1:n){
    for (k in 1:n){
      arg1      = 2*a/(((X[j]^(1/m)-X[k]^(1/m))^2) + a^2)
      arg2      = 2*a/(((X[j]-X[k])^2)+ a^2)
      arg3      = 2*a/(((X[j]^(1/m)-X[k])^2)+ a^2)
      arg4      = 2*a/(((X[k]^(1/m)-X[j])^2)+ a^2)
      matr[j,k] = (1/m^2)*arg1 - (1/m)*X[k]^(-betaH*(m-1))*arg3 -(1/m)*X[j]^(-betaH*(m-1))*arg4 +   X[j]^(-betaH*(m-1))*X[k]^(-betaH*(m-1))*arg2
    }
  }
  out = sum(matr)/n
  return(out)
}

#vectorised version
Tnm1 <- function(X,m,a,betaH){
  X = sort(X)
  n = length(X)
  # if (meth=="MLE"){betaH = n/sum(log(X))}
  # if (meth=="MME"){betaH = mean(X)/(mean(X)-1)}
  matr = matrix(0,n,n)
  arg1 <- 2*a/(outer(X^(1/m),X^(1/m),"-")^2 + a^2)
  arg2 <- 2*a/(outer(X^(1/m),X,"-")^2+ a^2)
  arg4 <- 2*a/(outer(X,X,"-")^2+ a^2)
  arg5 <- outer(X,X,"*")^(-betaH*(m-1))
  XX   <- expand.grid(X,X)
  ans  <- (1/n)*sum(arg1/m^2-(2/m)*XX[,2]^(-betaH*(m-1))*arg2 + arg5*arg4)  
  return(ans)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~ (0b) Testing for the Pareto distribution ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ Statistic using characterisation: X^(1/m) ~ X_(1) ~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ Normal weight function version ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#numeric integral version
Tnm2_int <- function(X,m,a,betaH){
  X = sort(X)
  n = length(X)
  ii <- 0+1i
  inn <- function(t)    (Mod((1/sqrt(n))*sum(exp(-ii*t*X^(1/m))/m - X^(-betaH*(m-1))*exp(-ii*t*X))))^2 * exp(-a*t^2)
  inn <- Vectorize(inn)
  # elliptic::myintegrate(inn,lower=0,upper=1)
  tmp <- try(integrate(inn,lower=-Inf,upper=Inf)$value)
  if(is.character(tmp)){
    out <- Inf
  } else {
    out <- tmp
  }
  return(out)
}

#double loop version
Tnm2_alt <- function(X,m,a,betaH){
  X = sort(X)
  n = length(X)
  matr = matrix(0,n,n)
  for (j in 1:n){
    for (k in 1:n){
      arg1      = ((X[j]^(1/m)-X[k]^(1/m))^2)/(4*a)
      arg2      = ((X[j]-X[k])^2)/(4*a)
      arg3      = ((X[j]^(1/m)-X[k])^2)/(4*a)
      arg4      = ((X[k]^(1/m)-X[j])^2)/(4*a)
      matr[j,k] =  (1/m^2)*exp(-arg1) -(1/m)*X[k]^(-betaH*(m-1))*exp(-arg3) -(1/m)*X[j]^(-betaH*(m-1))*exp(-arg4) + X[j]^(-betaH*(m-1))*X[k]^(-betaH*(m-1))*exp(-arg2)
    }
  }
  out = sum(matr)*sqrt(pi/a)/n
  return(out)
}

#vectorised version
Tnm2<- function(X,m,a,betaH){
  X = sort(X)
  n = length(X)
  matr = matrix(0,n,n)
  Xom.m.Xom2 <- outer(X^(1/m),X^(1/m),"-")^2
  Xom.m.X2   <- outer(X^(1/m),X,"-")^2
  X.m.Xom2   <- t(Xom.m.X2)
  X.m.X2     <- outer(X,X,"-")^2
  X.t.X      <- outer(X,X,"*")^(-betaH*(m-1))
  XX         <- expand.grid(X,X)
  ans        <- (1/n)*sqrt(pi/a)*sum(exp(-Xom.m.Xom2/(4*a))/m^2-(2/m)*XX[,2]^(-betaH*(m-1))*exp(-Xom.m.X2/(4*a)) + X.t.X*exp(-X.m.X2/(4*a)))
  return(ans)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (1) Kolmogorov-Smirnov ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KSn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  FH  = ppareto(X,betaH,sigmaH)
  KSp = max((1:n)/n-FH)
  KSm = max(FH-(0:(n-1))/n)
  out = max(KSp,KSm)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ (2) Cram?r-von Mises ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CVn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  FH  = ppareto(X,betaH,sigmaH)
  out = sum((FH-(2*(1:n)-1)/(2*n))^2)+1/(12*n)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ (3) Anderson-Darling ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ADn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  z   = ppareto(X,betaH,sigmaH)
  T1  = 2*(1:n)-1
  T2  = log(z) + log(1-z[n:1])
  out = -n-mean(T1*T2)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~ (4) Modified Anderson-Darling ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MAn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  z   = ppareto(X,betaH,sigmaH)
  S1  = sum(z)
  T1  = 2-(2*(1:n)-1)/n
  T2  = log(1-z)
  S2  = sum(T1*T2)
  out = n/2-2*S1-S2
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ (5) Zhang's A ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ZAn <- function(X,betaH,sigmaH=1,epsilon=1e-4){
  n   = length(X)
  X[X==min(X)] <- X[X==min(X)] + epsilon #add epsilon to the smallest value
  j   = 1:n
  FH  = ppareto(X,betaH,sigmaH)
  out = -sum(log(FH)/(n-j+0.5)+log(1-FH)/(j-0.5))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ (6) Zhang's B ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ZBn <- function(X,betaH,sigmaH=1,epsilon=1e-4){
  n   = length(X)
  X[X==min(X)] <- X[X==min(X)] + epsilon #add epsilon to the smallest value
  j   = 1:n
  FH  = ppareto(X,betaH,sigmaH)
  arg = (1/FH-1)/((n-0.5)/(j-0.75)-1)
  out = sum((log(arg))^2)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ (7) Zhang's C ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ZCn <- function(X,betaH,sigmaH=1){
  n   = length(X)
  j   = 1:n
  FH  = ppareto(X,betaH,sigmaH)
  T1  = n*(j-0.5)/(n-j+0.5)^2*log((j-0.5)/(n*FH))
  T2  = n/(n-j+0.5)*log((n-j+0.5)/(n*(1-FH)))
  out = 2*sum(T1+T2)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~ (8) Kullback-Leibler ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

H <- function(X,m){
  n    = length(X)
  j    = 1:n
  ind1 = pmin(j+m,n)
  ind2 = pmax(j-m,1)
  out  = mean(log(n/(2*m)*(X[ind1]-X[ind2])))
  return(out)
}

KLn <- function(X,m,betaH,sigmaH=1){
  n   = length(X)
  out = -H(X,m)-log(betaH)-betaH*log(sigmaH)+(betaH+1)*mean(log(X))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (9) Hellinger distance ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HDn <- function(X,m,betaH,sigmaH=1){
  n    = length(X)
  j    = 1:n
  ind1 = pmin(j+m,n)
  ind2 = pmax(j-m,1)
  num  = (n/(2*m)*(X[ind1]-X[ind2]))^(-0.5)-sqrt(dpareto(X,betaH,sigmaH))
  den  = (n/(2*m)*(X[ind1]-X[ind2]))^(-1)
  out  = sum(num^2/den)/(2*n)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~ Phi-divergence ~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fhHat <- function(X){
  n   = length(X)
  s   = sd(X)
  h   = 1.06*s*n^(-1/5)
  out = rep(0,n)
  for(j in 1:n){
    out[j] = mean(dnorm((X[j]-X)/h))/h
  }
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (10) Kullback-Leibler ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DKn <- function(X,betaH,sigmaH=1){
  fH  = fhHat(X)
  out = mean(log(fH/dpareto(X,betaH,sigmaH)))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (11) Meintanis (2009a) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MEnOLD <- function(X,a,betaH,sigmaH=1){
  n  = length(X)
  FH = ppareto(X,betaH,sigmaH)
  S1 = matrix(0,n,n)
  for (j in 1:n){
    for (k in 1:n){
      S1[j,k] = 2*a/((FH[j]-FH[k])^2+a^2)
    }
  }
  S2 = rep(0,n)
  for (j in 1:n){
    S2[j] = atan(FH[j]/a)+atan((1-FH[j])/a)
  }
  out = 2*n*(2*atan(1/a)-a*log(1+1/a^2))+sum(S1)/n-4*sum(S2)
  return(out)
}

MEn <- function(X,a,betaH,sigmaH=1){
  n  = length(X)
  FH = ppareto(X,betaH,sigmaH)
  S1 = rep(0,n)
  k  = 1:n
  for (j in 1:n){
    S1[j] = sum(2*a/((FH[j]-FH)^2+a^2))
  }
  S2  = atan(FH[k]/a)+atan((1-FH[k])/a)
  out = 2*n*(2*atan(1/a)-a*log(1+1/a^2))+sum(S1)/n-4*sum(S2)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ (12) Meintanis (2009b) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ia0 <- function(X,a){
  out = 1/(a+log(X))
}

Ia1 <- function(X,a){
  out = (1-a-log(X))/(a+log(X))^2
}

Ia2 <- function(X,a){
  out = (2-2*a+a^2+2*(a-1)*log(X)+(log(X))^2)/(a+log(X))^3
}

M2nOld <- function(X,a,betaH,sigmaH=1){
  n  = length(X)
  X  = X/sigmaH
  S1 = matrix(0,n,n)
  for (j in 1:n){
    for(k in 1:n){
      S1[j,k] = Ia0(X[j]*X[k],a)
    }
  }
  S2 = matrix(0,n,n)
  for (j in 1:n){
    for(k in 1:n){
      S2[j,k] = Ia2(X[j]*X[k],a)
    }
  }
  S3 = matrix(0,n,n)
  for (j in 1:n){
    for(k in 1:n){
      S3[j,k] = Ia1(X[j]*X[k],a)
    }
  }
  S4 = rep(0,n)
  for (j in 1:n){
    S4[j] = Ia0(X[j],a)
  }
  S5 = rep(0,n)
  for (j in 1:n){
    S5[j] = Ia1(X[j],a)
  }
  T1  = ((betaH+1)^2*sum(S1)+sum(S2)+2*(betaH+1)*sum(S3))/n
  T2  = (n*betaH*Ia0(1,a)-2*(betaH+1)*sum(S4)-2*sum(S5))*betaH
  out = T1+T2
}

M2n <- function(X,a,betaH,sigmaH=1){
  n  = length(X)
  X  = X/sigmaH
  S1 = rep(0,n)
  for (j in 1:n){
    S1[j] = sum(Ia0(X[j]*X,a))
  }
  S2 = rep(0,n)
  for (j in 1:n){
    S2[j] = sum(Ia2(X[j]*X,a))
  }
  S3 = rep(0,n)
  for (j in 1:n){
    S3[j] = sum(Ia1(X[j]*X,a))
  }
  S4  = Ia0(X,a)
  S5  = Ia1(X,a)
  T1  = ((betaH+1)^2*sum(S1)+sum(S2)+2*(betaH+1)*sum(S3))/n
  T2  = (n*betaH*Ia0(1,a)-2*(betaH+1)*sum(S4)-2*sum(S5))*betaH
  out = T1+T2
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~ (13) Obradovic et al. (2015) ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OJnMarius <- function(X){
  n    = length(X)
  N    = choose(n,2)
  tmp0 = combn(X,2)
  tmp1 = sum(rank(c(X,apply(cbind(tmp0[1,]/tmp0[2,],tmp0[2,]/tmp0[1,]),1,max)))[1:n])
  tmp3 = 1/(2*n*N)*(2*tmp1-(n+1)*(N+n))
  ans  = abs(tmp3)
  return(ans)
}

MnOld <- function(X,t){
  n = length(X)
  S = matrix(0,n,n)
  for (j in 1:(n-1)){
    for (k in (j+1):n){
      S[j,k] = (max(X[j]/X[k],X[k]/X[j])<=t)
    }
  }
  out = sum(S)/choose(n,2)
  return(out)
}

Mn <- function(X,t){
  n = length(X)
  S = rep(0,n)
  for (j in 1:(n-1)){
    k    = (j+1):n
    S[j] = sum(pmax(X[j]/X[k],X[k]/X[j])<=t)
  }
  out = sum(S)/choose(n,2)
  return(out)
}

MnVec <- function(X){
  n = length(X)
  S = matrix(0,n,n)
  for (j in 1:(n-1)){
    k = (j+1):n
    for (l in 1:n){
      S[j,l] = sum(pmax(X[j]/X[k],X[k]/X[j])<=X[l])
    }
  }
  out = colSums(S)/choose(n,2)
  return(out)
}

OJn <- function(X){
  n   = length(X)
  Mn  = MnVec(X)
  out = mean(Mn-(1:n)/n)
  return(abs(out))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ (14) Allison et al. (1) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Delta2Old <- function(n){
  out = rep(0,n)
  S   = matrix(0,n,n)
  for (l in 1:n){
    for (j in 1:n){
      for (k in 1:n){
        S[j,k] = (min(j,k)<=l)
      }
    }
    out[l] = sum(S)
  }
  return(out)
}

Delta2 <- function(n){
  out = rep(0,n)
  S   = rep(0,n)
  for (l in 1:n){
    for (j in 1:n){
      S[j] = sum((pmin(j,1:n)<=l))
    }
    out[l] = sum(S)
  }
  return(out)
}

Delta3 <- function(n){
  out = rep(0,n)
  S   = matrix(0,n,n)
  for (l in 1:n){
    for (j in 1:n){
      for (k in 1:n){
        S[j,k] = sum(pmin(j,k,1:n)<=l)
      }
    }
    out[l] = sum(S)
  }
  return(out)
}

A1nOld <- function(X,m){
  n = length(X)
  if(m == 2){ Delta = Delta2(n) }
  if(m == 3){ Delta = Delta3(n) }
  S1 = rep(0,n)
  S  = rep(0,n)
  for (l in 1:n){
    for (j in 1:n){
      S1[j] = (X[j]^(1/m)<X[l])
    }
    S[l] = mean(S1)
  }
  out = mean(S-Delta/n^m)
  return(out)
}

A1n <- function(X,m){
  n = length(X)
  if(m == 2){ Delta = Delta2(n) }
  if(m == 3){ Delta = Delta3(n) }
  S = rep(0,n)
  for (l in 1:n){
    S[l] = mean(X^(1/m)<X[l])
  }
  out = mean(S-Delta/n^m)
  return(abs(out))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~ (15) Allison et al. (2) ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A2n <- function(X,m){
  n = length(X)
  if(m == 2){ Delta = Delta2(n) }
  if(m == 3){ Delta = Delta3(n) }
  S = rep(0,n)
  for (l in 1:n){
    S[l] = mean(X^(1/m)<X[l])
  }
  out = mean((S-Delta/n^m)^2)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~ (16) Obradovic (1) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gn <- function(x,X){
  gx <- function(x,Xj,Xk,Xl){
    md  = median(c(Xj,Xk,Xl))
    mn  = min(c(Xj,Xk,Xl))
    out = (md/mn<=x)
    return(out)
  }
  n  = length(X)
  S1 = matrix(0,n,n)
  for (j in 1:(n-2)){
    for (k in (j+1):(n-1)){
      S1[j,k] = gx(x,X[j],X[k],X[n])*(n-k)
    }
  }
  S1  = sum(S1)
  S2a = matrix(0,n,n)
  for (k in 1:(n-1)){
    for (l in (k+1):n){
      S2a[k,l] = gx(x,X[k],X[l],X[l])
    }
  }
  S2a = sum(S2a)
  S2b = matrix(0,n,n)
  for (l in 1:(n-1)){
    for (k in (l+1):n){
      S2b[k,l] = gx(x,X[k],X[l],X[l])
    }
  }
  S2b = sum(S2b)
  S3  = rep(0,n)
  for (l in 1:n){
    S3[l] = gx(x,X[l],X[l],X[l])
  }
  S3  = sum(S3)
  out = (6*S1+3*(S2a+S2b)+S3)/n^3
  return(out)
}

Hn <- function(n){
  out = rep(0,n)
  S   = rep(0,n)
  for (l in 1:n){
    for (j in 1:n){
      S[j] = sum((pmin(j,1:n)<=l))
    }
    out[l] = sum(S)/n^2
  }
  return(out)
}

O1n <- function(X){
  n = length(X)
  G = rep(0,n)
  for(j in 1:n){
    G[j] = Gn(X[j],X)
  }
  out = mean(G-Hn(n))
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~ (17) Obradovic (2) ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Jn <- function(x,X){
  jx <- function(x,Xj,Xk,Xl){
    mx  = max(c(Xj,Xk,Xl))
    md  = median(c(Xj,Xk,Xl))
    out = (mx/md<=x)
    return(out)
  }
  n  = length(X)
  S1 = matrix(0,n,n)
  for (l in 3:n){
    for (k in 2:(l-1)){
      S1[k,l] = (X[l]/X[k]<x)*(k-1)
    }
  }
  S1  = sum(S1)
  S2a = matrix(0,n,n)
  for (k in 1:(n-1)){
    for (l in (k+1):n){
      S2a[k,l] = jx(x,X[k],X[l],X[l])
    }
  }
  S2a = sum(S2a)
  S2b = matrix(0,n,n)
  for (l in 1:(n-1)){
    for (k in (l+1):n){
      S2b[k,l] = jx(x,X[k],X[l],X[l])
    }
  }
  S2b = sum(S2b)
  S3  = rep(0,n)
  for (l in 1:n){
    S3[l] = jx(x,X[l],X[l],X[l])
  }
  S3  = sum(S3)
  out = (6*S1+3*(S2a+S2b)+S3)/n^3
  return(out)
}

Kn <- function(x,X){
  kx <- function(x,Xj,Xk,Xl){
    md  = median(c(Xj,Xk,Xl))
    mn  = (min(c(Xj,Xk,Xl)))
    out = ((md/mn)^2<=x)
    return(out)
  }
  n  = length(X)
  S1 = matrix(0,n,n)
  for (j in 1:(n-2)){
    for (k in (j+1):(n-1)){
      S1[j,k] = kx(x,X[j],X[k],X[n])*(n-k)
    }
  }
  S1  = sum(S1)
  S2a = matrix(0,n,n)
  for (k in 1:(n-1)){
    for (l in (k+1):n){
      S2a[k,l] = kx(x,X[k],X[l],X[l])
    }
  }
  S2a = sum(S2a)
  S2b = matrix(0,n,n)
  for (l in 1:(n-1)){
    for (k in (l+1):n){
      S2b[k,l] = kx(x,X[k],X[l],X[l])
    }
  }
  S2b = sum(S2b)
  S3  = rep(0,n)
  for (l in 1:n){
    S3[l] = kx(x,X[l],X[l],X[l])
  }
  S3  = sum(S3)
  out = (6*S1+3*(S2a+S2b)+S3)/n^3
  return(out)
}

O2n <- function(X){
  n   = length(X)
  out = rep(0,n)
  for(j in 1:n){
    out[j] = Jn(X[j],X)-Kn(X[j],X)
  }
  out = mean(out)
  return(out)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~ (18) Falkm Guillou, Toulemonde (2008) ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FGT <- function(X, J, method=1){
  n   <- length(X)
  
  Z   <- numeric(J)
  S   <- matrix(ncol=J,nrow=J)
  
  if(method==1 | str_to_lower(method) == "original"){
    sigtil <- 0.5*mean(X)*(((mean(X)^2)/(mean(X^2) - mean(X)^2)) + 1) 
    xitil  <- 0.5*(1 - (mean(X)^2)/(mean(X^2) - mean(X)^2))
    #cf. eq 11 of Falk(2008)
    betatil <- c(sigtil,xitil)
    
    if(is.character(try(log(1+xitil*X/sigtil)))){
      x=x
    }  
    
    elldot <- function(x){
      c(-1/sigtil + (1 + xitil)*(x/sigtil^2)*(1 + xitil*x/sigtil)^(-1),
        (1/xitil^2)*log(1+xitil*x/sigtil)-x*((1+xitil)/(xitil*sigtil))*(1+xitil*x/sigtil)^(-1))
    }  
    elldot <- Vectorize(elldot)
    elldotX <- elldot(X)
    Ibbinv <- matrix(c(2*sigtil^2*(1+xitil),-sigtil*(1+xitil),-sigtil*(1+xitil),(1+xitil)^2),ncol=2)
    meanelldot <- apply(elldotX,1,mean)
    betahat <- betatil + as.numeric(Ibbinv%*%meanelldot)
    sig <- betahat[1]
    xi  <- betahat[2]
  }else if(method==2 | str_to_lower(method) == "mme"){  #MME
    sigmP <- ParetoMME2(X)$sigmaH
    betaP <- ParetoMME2(X)$betaH
    sig <- sigmP/betaP
    xi <- 1/betaP
  } else if(method==3| str_to_lower(method) == "mle"){ #MLE
    sigmP <- ParetoMLE2(X)$sigmaH
    betaP <- ParetoMLE2(X)$betaH
    sig <- sigmP/betaP
    xi <- 1/betaP
  }  
  # cat("\n\n\nmeanelldot=")
  # print(meanelldot)
  # cat("\nIbbinv=")
  # print(Ibbinv)
  # cat("\nelldotX=")
  # print(elldotX)
  # cat("\nsigtil=",sigtil)
  # cat("\nxitil=",xitil)
  # cat("\nsig=",sig)
  # cat("\nxi=",xi,"\n\n\n")
  
  for(s in 1:J){
    tmp1 <- (1+xi*X/sig)
    tmp2 <- (-s/xi)
    
    Z[s] <- sum(
      (1+xi*X/sig)^(-s/xi)- 1/(s+1))/sqrt(n)
    # sign(tmp1)*abs(tmp1)^tmp2 - (1/(s+1)) )/sqrt(n)
    
    for(v in 1:J){
      u      <- s
      S[u,v] <- u*v/((u+v+1)*(u+1)*(v+1))-(u*v*(1+xi)*(u*v+xi+(u+1)*(v+1)))/
        ((v+xi+1)*(u+xi+1)*((u+1)^2)*((v+1)^2))
    }   
  }
  out <- as.numeric(t(Z)%*%ginv(S)%*%Z)
  return(out)
} 
# J <- 4
# FGT(X,J)
# qchisq(0.95,J)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~