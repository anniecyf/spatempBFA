rm(list=ls())
library(mvtnorm)
library(fields)
library(spBFA)
library(coda)
numSpatOverallGroups <- 2
M <- 361
m = 352 # 9 testing spatial points
K <- 2
O <- 1
L <- min(10, M)
LjVec <- rep(min(10, M), K)
sqrootM <- 19
Nu <- 310 
trainingT <- 300 # training set T = 300
testingT <- 10
Time <- 1:trainingT
TimeDist <- as.matrix(dist(1:Nu))
APsi = 0.1; BPsi = 4.5
set.seed(29)
calcGroup <- function(row, col, sqrootM){
  return((row - 1)*sqrootM + col)
}
sigma2 <- 0.01 # actual sigma^2(i,o) (for i=1,2,...,M and o=1) values
psi <- 2.3
kappa <- 0.7
tempMat <- matrix(runif(K*K,0,1), K, K)
Upsilon <- t(tempMat)%*%tempMat
rho <- 0.8
D <- rdist(expand.grid(1:sqrootM, 1:sqrootM))
Frho <- exp(-rho*D)
Hpsi <- exp(-psi*TimeDist)
Eta <- rmvnorm(1, mean=rep(0, Nu*K), sigma=kronecker(Hpsi, Upsilon)) # actual Eta (c(Eta_1,...,Eta_T)) (vec of length Nu*K)
Lambda <- matrix(5, M * O, K)
theta2j <- c(10, -10)
rowColInd <- rbind(expand.grid(1:10, 1:9), expand.grid(10:19, 11:19))
whichGroup1 <- mapply(calcGroup, rowColInd[,2], rowColInd[,1], sqrootM = sqrootM)
spatGroupOverall <- rep(2,M)
spatGroupOverall[whichGroup1] <- 1
#matrix(spatGroupOverall, 19, 19)
Lambda[,2] = theta2j[spatGroupOverall]
Sigma.NuMO <- matrix(rnorm(Nu * M * O, sd = sqrt(sigma2)), M*O, Nu)
EtaMat <- matrix(Eta, K, Nu)
meanMat <- Lambda%*%EtaMat #M*O\times Nu
YtrainingTemp <- as.vector(meanMat[,1:trainingT] + Sigma.NuMO[,1:trainingT])
YtestingTemp <- as.vector(meanMat[,(trainingT+1):Nu] + Sigma.NuMO[,(trainingT+1):Nu])
testingRowCol <- expand.grid(c(4, 10, 16), c(4, 10, 16))
whichTesting <- mapply(calcGroup, testingRowCol[,2], testingRowCol[,1], sqrootM = sqrootM)
Dtraining <- as.matrix(D[-whichTesting, -whichTesting])
discardInd <- vector()
for(whichTestingLoc in whichTesting){
  discardInd <- c(discardInd, seq(from = whichTestingLoc, by = M, length.out = Nu))
}
Ytraining <- YtrainingTemp[-discardInd]
YtestingSpat <- YtrainingTemp[discardInd] 
rm(YtrainingTemp)
distOrigNew = as.matrix(D[-whichTesting, whichTesting])
distNewNew = as.matrix(D[whichTesting, whichTesting])
dat = data.frame(Y = Ytraining); dist = Dtraining
tempTestingDiscardInd <- vector()
for(whichTestingLoc in whichTesting){
  tempTestingDiscardInd <- c(tempTestingDiscardInd, seq(from = whichTestingLoc, by = M, length.out = testingT))
}
YtestingTemp <- YtestingTemp[-tempTestingDiscardInd]
library(tidyverse)
xcoord <- rep(1:sqrootM, sqrootM)
ycoord <- rep(seq(sqrootM, 1, by = -1), each = sqrootM)
trainingTestingSpatGp <- rep(1, M); trainingTestingSpatGp[whichTesting] = 2
spatGroupOverall <- rep(2, M); spatGroupOverall[whichGroup1] <- 1

MCMC <- list(NBurn = 80000, NSims = 20000, NThin = 4, NPilot = 5)
reg.simu <- bfa_sp(Y ~ 0, data = dat, dist = Dtraining, time = Time,  K = K, 
                   starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                   L = L,
                   family = "normal",
                   trials = NULL,
                   temporal.structure = "exponential",
                   spatial.structure = "continuous",
                   seed = 29, 
                   gamma.shrinkage = TRUE,
                   include.space = TRUE,
                   clustering = TRUE) 
save(reg.simu, file="toyExspBFAL10.RData")
Diags <- spBFA::diagnostics(reg.simu, diags = c("dic", "dinf", "waic"), keepDeviance = TRUE)
save(Diags, file="toyExspBFAL10Diags.RData")
Deviance <- as.mcmc(Diags$deviance)
save(Deviance, file = "toyExspBFAL10Deviance.RData")
#predict.spBFA <- function(object, NewTimes, NewX = NULL, NewTrials = NULL, type = "temporal", Verbose = TRUE, seed = 54, ...)
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 29)
save(temppredspBFAL10, file = "toyExspBFAL10temppred.RData")
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 19)
save(temppredspBFAL10, file = "s19toyExspBFAL10temppred.RData")
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 27)
save(temppredspBFAL10, file = "s27toyExspBFAL10temppred.RData")
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 31)
save(temppredspBFAL10, file = "s31toyExspBFAL10temppred.RData")