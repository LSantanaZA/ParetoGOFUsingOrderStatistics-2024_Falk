
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial settings ~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
cat("\014")
# set.seed(1234)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library("stargazer")
library("stringr")
source("Estimation.R")
source("GoF_Pareto_MME.R")
source("DistsForSim.2023.08.11.15h42.R")

starttime <- format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n       = 20
Range   = c(3:4,5:8,9:12,13:16,21:24,25:28,33:36,73:76,89:92,46,48,50,42,64,66,68,70)
MC      = 2e4
alpha   = 0.05
meth    = "MME" # "MLE" #or "MME"
method  = 1 # 1==original FALK, 2==MME, 3==MLE, 
J       = 4
# Results = matrix(0,length(Range),2)
Results = matrix(0,length(Range),3)
nm      = numeric(length(Range))
tme1    = proc.time()

count = 0
pb    = txtProgressBar(min=0,max=length(Range)*MC,style=3)
DistNumcntr <- 1
for (DistNum in Range){#START FOR LOOP (DistNum in Range)
  TS_FGT = TS_FGT_S = numeric(MC)
  j      = 1
  # for (j in 1:MC){#START FOR LOOP (j in 1:MC)
  while(j <= MC){#START FOR LOOP (j in 1:MC)
    simdata         = SimFromDist(DistNum,n)
    X               = sort(simdata$X)
    nm[DistNumcntr] = simdata$nm
    betaH           = ifelse(meth=="MME",ParetoMME2(X)$betaH ,ParetoMLE2(X)$betaH) # MME or MLE of
    sigmaH          = ifelse(meth=="MME",ParetoMME2(X)$sigmaH ,ParetoMLE2(X)$sigmaH) # MME or MLE of sigma is calculated
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    X               = X - sigmaH # Shifting the data (tranforming for the FGT test to work).
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    trytmp          = try(fgt <- FGT(X,J,method=method),silent=TRUE)
    if(!is.character(trytmp) & !is.na(trytmp)){
      TS_FGT[j]   = fgt
      
      
      
      Xstar       = sort(rpareto(n,betaH,sigmaH))
      sigmaHS     = ifelse(meth=="MME",ParetoMME2(Xstar)$sigmaH ,ParetoMLE2(Xstar)$sigmaH) # MME or MLE
      Xstar       = Xstar - sigmaHS
      
      # Xstar       = sort(rLomax(n,c(1/betaH,sigmaH/betaH)))
      trytmps     = try(fgts <- FGT(Xstar,J,method=method),silent=TRUE)
      if(!is.character(trytmps) & !is.na(trytmps)){
        TS_FGT_S[j] = fgts
        
        j           = j + 1
        count       = count+1  
        setTxtProgressBar(pb,count)
      }
    }  
    
    
  }#END FOR LOOP (j in 1:MC)
  
  p     = mean(TS_FGT > qchisq(1-alpha,J))
  # cat("\nCV=",qchisq(1-alpha,J),"(theoretical)")
  BS_CV = sort(TS_FGT_S)[floor((1-alpha)*MC)]
  # cat("\nCV=",BS_CV,"(BS)\n")
  pbs   = mean(TS_FGT > BS_CV)
  
  
  Results[DistNumcntr,] = c(DistNum,p,pbs)
  # Results[DistNumcntr,] = c(DistNum,p)
  dimnames(Results) <- list(NULL,c("DistNum","power (chisq CV)", "power (BS CV)"))
  # dimnames(Results) <- list(NULL,c("DistNum","p-value (chisq CV)"))
  #	cat(Results[DistNum,],"\n",file=filnam,append=TRUE)
  
  # imagename <- paste("Pareto_n=",n,"_MC=",MC,"_parm=1_meth=",meth,".RData",sep="")
  # save.image(imagename)
  DistNumcntr <- DistNumcntr + 1
  
}#END FOR LOOP (DistNum in Range)

tme2 = proc.time()
tme  = tme2[3]-tme1[3]

Results