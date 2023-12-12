################## no seasonality (AR(1) or exponential) ###################
rm(list=ls())
Nu <- 25
Time = 1:Nu
TimeDist <- as.matrix(dist(1:Nu))
psi <- 0.8
Hpsi <- psi^TimeDist
I_Nu <- diag(Nu)
##### SampleEta() initial trials (for t=1,2,3,...,T-1,T) #####
### to find the closed-form representations for final quantities  
HpsiSub <- list(); HpsiSubInv <- list(); Hplus <- list(); Hstar <- list()
for(t in 1:Nu){
  HpsiSub[[t]] <- Hpsi[-t,-t]
  HpsiSubInv[[t]] <- solve(HpsiSub[[t]])
  Hplus[[t]] <- Hpsi[t,-t]%*%HpsiSubInv[[t]]
  Hstar[[t]] <- Hpsi[t,t] - Hplus[[t]]%*%Hpsi[-t,t]
}
##### SampleEta() subsequent verification trials (for t=2,3,...,T-1) #####
### to find and verify the closed-form representations for middle quantities  
### utilized in deducing the final products' representations
getMaxAbsError <- function(Nu, psi){
  er <- 0
  Time = 1:Nu
  TimeDist <- as.matrix(dist(1:Nu))
  Hpsi <- psi^TimeDist
  I_Nu <- diag(Nu)
  for(t in 2:(Nu-1)){
    Ht <- Hpsi[-t,-t]
    Lt <- Ht
    Lt[upper.tri(Ht)] <- 0
    Dt <- diag(c(1, rep(1-psi^2,(Nu-2))))
    Dt[t,t] <- 1-psi^4
    LtInv <- diag(Nu-1)
    diag(LtInv[-1,-(Nu-1)]) = -psi
    LtInv[t,t-1] <- -psi^2
    HtInv <- (1+psi^2)/(1-psi^2)*diag(Nu-1)
    HtInv[Nu-1,Nu-1] <- 1/(1-psi^2)
    HtInv[1,1] <- 1/(1-psi^2)
    diag(HtInv[-1,-(Nu-1)]) = -psi/(1-psi^2)
    diag(HtInv[-(Nu-1),-1]) = -psi/(1-psi^2)
    HtInv[t-1,t] <- -psi^2/(1-psi^4); HtInv[t,t-1] <- -psi^2/(1-psi^4)
    if(t==2) HtInv[t-1,t-1] = 1/(1-psi^4) else HtInv[t-1,t-1] = (1+psi^2+psi^4)/(1-psi^4)
    if(t==(Nu-1)) HtInv[t,t] = 1/(1-psi^4) else HtInv[t,t] = (1+psi^2+psi^4)/(1-psi^4)
    Hstar <- (1-psi^2)/(1+psi^2)
    Hplus <- rep(0, (Nu-1))
    Hplus[t] <- psi/(1+psi^2); Hplus[t-1] <- psi/(1+psi^2)
    e1 <- max(abs(Ht-Lt%*%Dt%*%t(Lt))) 
    e2 <- max(abs(diag(Nu-1)-Lt%*%LtInv)) 
    e3 <- max(abs(t(LtInv)%*%diag(1/diag(Dt))%*%LtInv%*%Ht-diag(Nu-1))) 
    e4 <- max(abs(HtInv - solve(Ht)))
    e5 <- max(abs(Hplus - Hpsi[t,-t]%*%HtInv))
    e6 <- abs(Hstar-(Hpsi[t,t]-Hpsi[t,-t]%*%HtInv%*%Hpsi[-t,t]))
    er <- max(er, e1, e2, e3, e4, e5, e6)
  }
  return(er)
}
getMaxAbsError(17,0.7)

##### SamplePsi() trials for the full Hpsi T by T matrix #####
L <- Hpsi
L[upper.tri(Hpsi)] <- 0
D <- diag(c(1, rep(1-psi^2,(Nu-1))))
max(abs(Hpsi-L%*%D%*%t(L))) #Hpsi=L%*%D%*%t(L)) 1.110223e-16
LInv <- diag(Nu)
diag(LInv[-1,-Nu]) = -psi
max(abs(diag(Nu)-L%*%LInv)) #5.551115e-17
HpsiInv <- (1+psi^2)*diag(Nu)
HpsiInv[1,1] <- 1; HpsiInv[Nu,Nu] <- 1
diag(HpsiInv[-1,-Nu]) = -psi; diag(HpsiInv[-Nu,-1]) = -psi
HpsiInv <- HpsiInv/(1-psi^2)
max(abs(HpsiInv%*%Hpsi-diag(Nu))) #8.881784e-16
max(abs(HpsiInv-t(LInv)%*%diag(c(1, rep(1/(1-psi^2),(Nu-1))))%*%LInv)) #4.440892e-16


################## with seasonality (SAR(1) or sexponential) ###################
##### SampleEta() initial trials (for t=1,2,3,...,T-1,T) #####
### to find the closed-form representations for final quantities  
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
setwd("~/spatTempBFA/src")
sourceCpp("CovarianceFunctions.cpp")
Nu <- 27
d <- 6
Time = 1:Nu
TimeDist <- as.matrix(dist(1:Nu))
psi <- 0.7
Hpsi <- H(psi, TempCorInd = 2, TimeDist, Nu, seasonPeriod = d)
Ht <- list(); tVec <- list(); HtInv <- list(); Hplus <- list(); Hstar <- list()
for(t in 1:Nu){
  tVec[[t]] <- Hpsi[t,-t]
  Ht[[t]] <- Hpsi[-t,-t]
  HtInv[[t]] <- solve(Ht[[t]])
  Hplus[[t]] <- tVec[[t]]%*%HtInv[[t]]
  Hstar[[t]] <- 1 - Hplus[[t]]%*%Hpsi[-t,t]
}
##### SampleEta() subsequent verification trials (for t=2,3,...,T-1) #####
### to find and verify the closed-form representations for middle quantities  
### utilized in deducing the final products' representations
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
setwd("~/spatTempBFA/src")
sourceCpp("CovarianceFunctions.cpp")
getMaxAbsDiff <- function(Nu, d, psi){ #TempCorInd = 2 
  er <- 0 #max abs error
  Time = 1:Nu
  TimeDist <- as.matrix(dist(1:Nu))
  Hpsi <- H(psi, TempCorInd=2, TimeDist, Nu, seasonPeriod = d)
  for(t in 2:(Nu-1)){
    Ht <- Hpsi[-t,-t]
    Lt <- Ht
    Lt[upper.tri(Ht)] <- 0
    LtInv <- diag(Nu-1)
    diag(LtInv[-(1:d),-((Nu-d):(Nu-1))]) <- -psi
    HtInv <- (1+psi^2)/(1-psi^2) * diag(Nu-1)
    diag(HtInv[-(1:d), -((Nu-d):(Nu-1))]) <- -psi/(1-psi^2)
    diag(HtInv[-((Nu-d):(Nu-1)), -(1:d)]) <- -psi/(1-psi^2)
    Hplus <- rep(0, Nu-1)
    if (t>1 & t<=d){ ### for t = 2,3,...,d
      for(j in 1:(t-1)){
        LtInv[j+d, j] <- 0
        LtInv[j+d-1, j] <- -psi
        HtInv[j+d, j] <- 0
        HtInv[j+d-1, j] <- -psi/(1-psi^2)
        i = j
        HtInv[i, i+d] <- 0
        HtInv[i, i+d-1] <- -psi/(1-psi^2)
      }
      Dt <- diag(c(rep(1, d-1), rep((1-psi^2), (Nu-d)))) 
      Dt[d+t-1, d+t-1] <- 1
      diag(HtInv[1:(d-1), 1:(d-1)]) <- 1/(1-psi^2)
      diag(HtInv[(Nu-d):(Nu-1), (Nu-d):(Nu-1)]) <- 1/(1-psi^2)
      HtInv[t+d-1, t+d-1] <- 1/(1-psi^2)
      Hplus[t+d-1] <- psi
      Hstar <- 1-psi^2
    } else if (t<Nu & t>(Nu-d)){### for t = T-d+1,...,T-2,T-1
      for(i in t:(Nu-1)){
        LtInv[i, i-d] <- 0
        LtInv[i, i-d+1] <- -psi
        HtInv[i, i-d] <- 0
        HtInv[i, i-d+1] <- -psi/(1-psi^2)
        j = i
        HtInv[j-d, j] <- 0
        HtInv[j-d+1, j] <- -psi/(1-psi^2)
      }
      Dt <- diag(c(rep(1, d), rep((1-psi^2), (Nu-d-1)))) 
      diag(HtInv[1:d, 1:d]) <- 1/(1-psi^2)
      diag(HtInv[(Nu-d+1):(Nu-1), (Nu-d+1):(Nu-1)]) <- 1/(1-psi^2)
      HtInv[t-d, t-d] <- 1/(1-psi^2)
      Hplus[t-d] <- psi
      Hstar <- 1-psi^2
    } else{### for t = d+1,d+2,...,T-d
      for(j in (t-d+1):(t-1)){
        LtInv[j+d, j] <- 0
        LtInv[j+d-1, j] <- -psi
        HtInv[j+d, j] <- 0
        HtInv[j+d-1, j] <- -psi/(1-psi^2)
        i = j
        HtInv[i, i+d] <- 0
        HtInv[i, i+d-1] <- -psi/(1-psi^2)
      }
      j <- t - d
      LtInv[j+d, j] <- 0
      LtInv[j+2*d-1, j] <- -psi^2
      HtInv[j+d, j] <- 0
      HtInv[j+2*d-1, j] <- -psi^2/(1-psi^4)
      i = j
      HtInv[i, i+d] <- 0
      HtInv[i, i+2*d-1] <- -psi^2/(1-psi^4)
      #j <- t+d-1
      #HtInv[j, j-d] <- 0
      #HtInv[j, j-2*d+1] <- -psi^2/(1-psi^4)
      Dt <- diag(c(rep(1, d), rep((1-psi^2), (Nu-d-1)))) 
      Dt[t+d-1, t+d-1] <- 1 - psi^4
      diag(HtInv[(1:d), (1:d)]) <- 1/(1-psi^2)
      diag(HtInv[(Nu-d):(Nu-1), (Nu-d):(Nu-1)]) <- 1/(1-psi^2)
      if(t<=(2*d)) {HtInv[t-d, t-d] <- 1/(1-psi^4)}
      else {HtInv[t-d, t-d] <- (psi^4+psi^2+1)/(1-psi^4)}
      if(t>(Nu-2*d)) {HtInv[t+d-1, t+d-1] <- 1/(1-psi^4)}
      else  {HtInv[t+d-1, t+d-1] <- (psi^4+psi^2+1)/(1-psi^4)}
      Hplus[t+d-1] <- psi/(1+psi^2)
      Hplus[t-d] <- psi/(1+psi^2)
      Hstar <- (1-psi^2)/(1+psi^2)
    }
    e1 <- max(abs(HtInv - t(LtInv)%*%solve(Dt)%*%LtInv))
    e2 <- max(abs(Lt%*%LtInv-diag(Nu-1))) #0
    e3 <- max(abs(Ht - Lt%*%Dt%*%t(Lt)))
    e4 <- max(abs(HtInv - solve(Ht)))
    e5 <- max(abs(Hplus-Hpsi[t,-t]%*%HtInv))
    e6 <- abs(Hstar-(1-Hpsi[t,-t]%*%HtInv%*%Hpsi[-t,t]))
    er <- max(er, e1, e2, e3, e4, e5, e6)
  }
  return(er)
}
Nu <- 27
d <- 6
psi <- 0.7
getMaxAbsDiff(Nu, d, psi)
getMaxAbsDiff(50, 7, 0.7)
#(possibly) in HtInv
-psi/(1-psi^2); -psi^2/(1-psi^4); 1/(1-psi^2); (1+psi^2)/(1-psi^2); 1/(1-psi^4); 1/(1-psi^2) + psi^4/(1-psi^4)
#vt <- Hpsi[t,-t]; Hplus <- vt%*%HtInv; Hstar <- 1 - Hplus%*%t(vt)
psi/(1+psi^2)

##### SamplePsi() initial trials #####
### to find the closed-form representations for final quantities  
HpsiInv <- solve(Hpsi)
L <- Hpsi
L[upper.tri(Hpsi)] <-0
LInv <- solve(L)
D <- diag(c(rep(1,d), rep((1-psi^2),Nu-d)))
DInv <- diag(1/diag(D))
rootiH <- t(LInv)%*%sqrt(DInv)
max(abs(rootiH%*%t(rootiH) - HpsiInv))#1.776357e-15
max(abs(L%*%D%*%t(L)-Hpsi))#1.110223e-16
##### SamplePsi() subsequent verification trials #####
rm(list=ls())
setwd("spatTempBFA/src")
sourceCpp("CovarianceFunctions.cpp")
Nu <- 27
d <- 6
Time = 1:Nu
TimeDist <- as.matrix(dist(1:Nu))
Psi <- 0.7
Hpsi <- H(Psi, TempCorInd = 2, TimeDist, Nu, seasonPeriod = d)
HInv <- getInvH(Psi, Nu, TempCorInd = 2, seasonPeriod = d)
rootiH <- getRootiH(Psi, Nu, TempCorInd = 2, seasonPeriod = d)
max(abs(Hpsi%*%HInv-diag(Nu)))
max(abs(HInv-rootiH%*%t(rootiH)))
HpsiInv <- solve(Hpsi)
L <- Hpsi
L[upper.tri(Hpsi)] <-0
LInv <- solve(L)
D <- diag(c(rep(1,d), rep((1-Psi^2),Nu-d)))
DInv <- diag(1/diag(D))
rootiH.exact <- t(LInv)%*%sqrt(DInv)
max(abs(HpsiInv-HInv))
max(abs(rootiH-rootiH.exact))
expHpsi <- H(Psi, TempCorInd = 3, TimeDist, Nu, seasonPeriod = d)
expHInv <- getInvH(Psi, Nu, TempCorInd = 3, seasonPeriod = d)
rootiHexp <- getRootiH(Psi, Nu, TempCorInd = 3, seasonPeriod = d)
max(abs(expHpsi%*%expHInv-diag(Nu)))
max(abs(expHInv-rootiHexp%*%t(rootiHexp)))
expHpsiInv <- solve(expHpsi)
L <- expHpsi
L[upper.tri(expHpsi)] <-0
LInv <- solve(L)
D <- diag(c(rep(1,d), rep((1-exp(-Psi*2)),Nu-d)))
DInv <- diag(1/diag(D))
expRootiH.exact <- t(LInv)%*%sqrt(DInv)
max(abs(expHpsiInv-expHInv))
max(abs(rootiHexp-expRootiH.exact))






