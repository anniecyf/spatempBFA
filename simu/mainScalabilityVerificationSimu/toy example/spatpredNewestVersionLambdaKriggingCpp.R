rm(list=ls())
library(mvtnorm)
library(fields)
library(coda)
library(tidyverse)
library(ggpubr)
library(devtools)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/server/spatempBFA")
load_all(".")
library(spatempBFA)
numSpatOverallGroups <- 2
M <- 361
m = 352 # 9 testing spatial points
K <- 2
O <- 1
L <- min(10, M)
LjVec <- rep(min(10, M), K)
sqrootM <- 19
Nu <- 60 
trainingT <- 50 # training set T = 50
testingT <- 10
Time <- 1:trainingT
TimeDist <- as.matrix(dist(1:60))
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
  discardInd <- c(discardInd, seq(from = whichTestingLoc, by = M, length.out = trainingT))
}
Ytraining <- YtrainingTemp[-discardInd]
YtestingSpat <- YtrainingTemp[discardInd] 
rm(YtrainingTemp)
distOrigNew = as.matrix(D[-whichTesting, whichTesting])
distNewNew = as.matrix(D[whichTesting, whichTesting])
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/mainScalabilityVerificationSimu/toy example", sep = "/"))
setwd("fullGPfixedL")
load("toyExfullGPfixedL.RData")
temppredFixedL <- predictNewTime(regFixedL.simu, (trainingT+1):Nu, seed = 27)
save(temppredFixedL, file = "s27toyExFullGPtemppredFixedL.RData")
temppredFixedL <- predictNewTime(regFixedL.simu, (trainingT+1):Nu, seed = 19)
save(temppredFixedL, file = "s19toyExFullGPtemppredFixedL.RData")
temppredFixedL <- predictNewTime(regFixedL.simu, (trainingT+1):Nu, seed = 31)
save(temppredFixedL, file = "s31toyExFullGPtemppredFixedL.RData")
setwd("../NNGPblockFixedL")
load("toyExNNGPblockFixedL.RData")
temppredFixedLblock <- predictNewTime(regFixedL.simu.block, (trainingT+1):Nu, seed = 27)
save(temppredFixedLblock, file = "s27toyExNNGPtemppredFixedLblock.RData")
temppredFixedLblock <- predictNewTime(regFixedL.simu.block, (trainingT+1):Nu, seed = 19)
save(temppredFixedLblock, file = "s19toyExNNGPtemppredFixedLblock.RData")
temppredFixedLblock <- predictNewTime(regFixedL.simu.block, (trainingT+1):Nu, seed = 31)
save(temppredFixedLblock, file = "s31toyExNNGPtemppredFixedLblock.RData")
setwd("../NNGPsequenFixedL")
load("toyExNNGPsequenFixedL.RData")
temppredFixedLsequen <- predictNewTime(regFixedL.simu.sequen, (trainingT+1):Nu, seed = 27)
save(temppredFixedLsequen, file = "s27toyExNNGPtemppredFixedLsequen.RData")
temppredFixedLsequen <- predictNewTime(regFixedL.simu.sequen, (trainingT+1):Nu, seed = 19)
save(temppredFixedLsequen, file = "s19toyExNNGPtemppredFixedLsequen.RData")
temppredFixedLsequen <- predictNewTime(regFixedL.simu.sequen, (trainingT+1):Nu, seed = 31)
save(temppredFixedLsequen, file = "s31toyExNNGPtemppredFixedLsequen.RData")
setwd("../NNGPsequenVaryLj")
load("toyExNNGPsequenVaryLj.RData")
temppredVaryLj <- predictNewTime(regVaryLj.simu.sequen, (trainingT+1):Nu, seed = 27)
save(temppredVaryLj, file = "s27toyExNNGPtemppredVaryLj.RData")
temppredVaryLj <- predictNewTime(regVaryLj.simu.sequen, (trainingT+1):Nu, seed = 19)
save(temppredVaryLj, file = "s19toyExNNGPtemppredVaryLj.RData")
temppredVaryLj <- predictNewTime(regVaryLj.simu.sequen, (trainingT+1):Nu, seed = 31)
save(temppredVaryLj, file = "s31toyExNNGPtemppredVaryLj.RData")
library(spBFA)
setwd("../spBFA")
load("toyExspBFAL10.RData")
load("toyExspBFALInf.RData")
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 19)
save(temppredspBFAL10, file = "s19toyExspBFAL10temppred.RData")
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 27)
save(temppredspBFAL10, file = "s27toyExspBFAL10temppred.RData")
temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 31)
save(temppredspBFAL10, file = "s31toyExspBFAL10temppred.RData")
temppredspBFALInf <- predict(reg.simu.LInf, (trainingT+1):Nu, seed = 19)
save(temppredspBFALInf, file = "s19toyExspBFALInftemppred.RData")
temppredspBFALInf <- predict(reg.simu.LInf, (trainingT+1):Nu, seed = 27)
save(temppredspBFALInf, file = "s27toyExspBFALInftemppred.RData")
temppredspBFALInf <- predict(reg.simu.LInf, (trainingT+1):Nu, seed = 31)
save(temppredspBFALInf, file = "s31toyExspBFALInftemppred.RData")
setwd("..")
rm(list=ls())
# setwd("fullGPfixedL")
# load("toyExfullGPfixedL.RData")
# spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 9, distOrigNew, distNewNew, 
#                                       NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
# save(spatpredFixedL, file = "toyExFullGPspatpredFixedL.RData")
# setwd("../NNGPblockFixedL")
# load("toyExNNGPblockFixedL.RData")
# spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 9, distOrigNew, distNewNew, 
#                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
# save(spatpredFixedLblock, file = "toyExNNGPspatpredFixedLblock.RData")
# setwd("../NNGPsequenFixedL")
# load("toyExNNGPsequenFixedL.RData")
# spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 9, distOrigNew, distNewNew, 
#                                             NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
# save(spatpredFixedLsequen, file = "toyExNNGPspatpredFixedLsequen.RData")
# setwd("../NNGPsequenVaryLj")
# load("toyExNNGPsequenVaryLj.RData")
# spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 9, distOrigNew, distNewNew, 
#                                       NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
# save(spatpredVaryLj, file = "toyExNNGPspatpredVaryLj.RData")
setwd("fullGPfixedL")
load("toyExfullGPfixedL.RData")
spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
save(spatpredFixedL, file = "s27toyExFullGPspatpredFixedL.RData")
setwd("../NNGPblockFixedL")
load("toyExNNGPblockFixedL.RData")
spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 9, distOrigNew, distNewNew, 
                                           NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
save(spatpredFixedLblock, file = "s27toyExNNGPspatpredFixedLblock.RData")
setwd("../NNGPsequenFixedL")
load("toyExNNGPsequenFixedL.RData")
spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 9, distOrigNew, distNewNew, 
                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
save(spatpredFixedLsequen, file = "s27toyExNNGPspatpredFixedLsequen.RData")
setwd("../NNGPsequenVaryLj")
load("toyExNNGPsequenVaryLj.RData")
spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
save(spatpredVaryLj, file = "s27toyExNNGPspatpredVaryLj.RData")
setwd("fullGPfixedL")
load("toyExfullGPfixedL.RData")
spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 31)
save(spatpredFixedL, file = "s31toyExFullGPspatpredFixedL.RData")
setwd("../NNGPblockFixedL")
load("toyExNNGPblockFixedL.RData")
spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 9, distOrigNew, distNewNew, 
                                           NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 31)
save(spatpredFixedLblock, file = "s31toyExNNGPspatpredFixedLblock.RData")
setwd("../NNGPsequenFixedL")
load("toyExNNGPsequenFixedL.RData")
spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 9, distOrigNew, distNewNew, 
                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 31)
save(spatpredFixedLsequen, file = "s31toyExNNGPspatpredFixedLsequen.RData")
setwd("../NNGPsequenVaryLj")
load("toyExNNGPsequenVaryLj.RData")
spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 31)
save(spatpredVaryLj, file = "s31toyExNNGPspatpredVaryLj.RData")
# setwd("fullGPfixedL")
# load("toyExfullGPfixedL.RData")
# spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 9, distOrigNew, distNewNew, 
#                                       NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 19)
# save(spatpredFixedL, file = "s19toyExFullGPspatpredFixedL.RData")
# setwd("../NNGPblockFixedL")
# load("toyExNNGPblockFixedL.RData")
# spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 9, distOrigNew, distNewNew, 
#                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 19)
# save(spatpredFixedLblock, file = "s19toyExNNGPspatpredFixedLblock.RData")
# setwd("../NNGPsequenFixedL")
# load("toyExNNGPsequenFixedL.RData")
# spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 9, distOrigNew, distNewNew, 
#                                             NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 19)
# save(spatpredFixedLsequen, file = "s19toyExNNGPspatpredFixedLsequen.RData")
# setwd("../NNGPsequenVaryLj")
# load("toyExNNGPsequenVaryLj.RData")
# spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 9, distOrigNew, distNewNew, 
#                                       NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 19)
# save(spatpredVaryLj, file = "s19toyExNNGPspatpredVaryLj.RData")
setwd("fullGPfixedL")
load("toyExfullGPfixedL.RData")
spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedL, file = "s29toyExFullGPspatpredFixedL.RData")
setwd("../NNGPblockFixedL")
load("toyExNNGPblockFixedL.RData")
spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 9, distOrigNew, distNewNew, 
                                           NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedLblock, file = "s29toyExNNGPspatpredFixedLblock.RData")
setwd("../NNGPsequenFixedL")
load("toyExNNGPsequenFixedL.RData")
spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 9, distOrigNew, distNewNew, 
                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedLsequen, file = "s29toyExNNGPspatpredFixedLsequen.RData")
setwd("../NNGPsequenVaryLj")
load("toyExNNGPsequenVaryLj.RData")
spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredVaryLj, file = "s29toyExNNGPspatpredVaryLj.RData")
rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/mainScalabilityVerificationSimu/toy example", sep = "/"))
setwd("fullGPfixedL")
load("s29toyExFullGPspatpredFixedL.RData")
#load("s27toyExFullGPspatpredFixedL.RData")
#load("s19toyExFullGPspatpredFixedL.RData")
#load("s31toyExFullGPspatpredFixedL.RData")
setwd("../NNGPblockFixedL")
load("s29toyExNNGPspatpredFixedLblock.RData")
#load("s27toyExNNGPspatpredFixedLblock.RData")
#load("s19toyExNNGPspatpredFixedLblock.RData")
#load("s31toyExNNGPspatpredFixedLblock.RData")
setwd("../NNGPsequenFixedL")
load("s29toyExNNGPspatpredFixedLsequen.RData")
#load("s27toyExNNGPspatpredFixedLsequen.RData")
#load("s19toyExNNGPspatpredFixedLsequen.RData")
#load("s31toyExNNGPspatpredFixedLsequen.RData")
setwd("../NNGPsequenVaryLj")
load("s29toyExNNGPspatpredVaryLj.RData")
#load("s27toyExNNGPspatpredVaryLj.RData")
#load("s19toyExNNGPspatpredVaryLj.RData")
#load("s31toyExNNGPspatpredVaryLj.RData")
setwd("..")
numSpatOverallGroups <- 2
M <- 361
m = 352 # 9 testing spatial points
K <- 2
O <- 1
L <- min(10, M)
LjVec <- rep(min(10, M), K)
sqrootM <- 19
Nu <- 60 
trainingT <- 50 # training set T = 50
testingT <- 10
Time <- 1:trainingT
TimeDist <- as.matrix(dist(1:60))
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
  discardInd <- c(discardInd, seq(from = whichTestingLoc, by = M, length.out = trainingT))
}
Ytraining <- YtrainingTemp[-discardInd]
YtestingSpat <- YtrainingTemp[discardInd]  #YtestingSpat <- as.vector(t(matrix(YtestingSpat, M-m, trainingT)))
rm(YtrainingTemp)
spatPredMetricMat <- matrix(0, 3, 4, 
                            dimnames = list(c("postMeanMSE", "postMSE", "postVar"),
                                            c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj")))
testingLocNum <- M - m
ylocPredList = spatpredFixedL$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
spatPredMetricMat[1, 1] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 1] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 1] = mean(apply(ylocPred, 1, var))
ylocPredList = spatpredFixedLblock$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
spatPredMetricMat[1, 2] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 2] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 2] = mean(apply(ylocPred, 1, var))
ylocPredList = spatpredFixedLsequen$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
spatPredMetricMat[1, 3] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 3] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 3] = mean(apply(ylocPred, 1, var))
ylocPredList = spatpredVaryLj$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
spatPredMetricMat[1, 4] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 4] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 4] = mean(apply(ylocPred, 1, var))
spatPredMetricMat
spatPredKrigTime <- matrix(0, 2, 4, 
                           dimnames = list(c("alpha", "weightsXiLambda"),
                                           c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj")))
spatPredKrigTime[1, 1] = spatpredFixedL$alphaKrigTime
spatPredKrigTime[2, 1] = spatpredFixedL$weightsXiLambdaKrigTime
spatPredKrigTime[1, 2] = spatpredFixedLblock$alphaKrigTime
spatPredKrigTime[2, 2] = spatpredFixedLblock$weightsXiLambdaKrigTime
spatPredKrigTime[1, 3] = spatpredFixedLsequen$alphaKrigTime
spatPredKrigTime[2, 3] = spatpredFixedLsequen$weightsXiLambdaKrigTime
spatPredKrigTime[1, 4] = spatpredVaryLj$alphaKrigTime
spatPredKrigTime[2, 4] = spatpredVaryLj$weightsXiLambdaKrigTime
spatPredKrigTime