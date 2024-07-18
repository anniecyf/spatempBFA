rm(list=ls())
library(mvtnorm)
library(fields)
library(spatempBFA)
library(coda)
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
spatGpDF <- data.frame(x = xcoord, y = ycoord, actualGp <- as.factor(spatGroupOverall), 
                       trainingTestingSpatGp <- as.factor(trainingTestingSpatGp))
actualSpatPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = actualGp), size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
actualSpatPlot
trainingTestingPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = trainingTestingSpatGp), size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
trainingTestingPlot

MCMC <- list(NBurn = 100, NSims = 100, NThin = 2, NPilot = 5)
regFixedL.simu <- bfaFixedL(Y ~ 0, data = dat, dist = Dtraining, time = Time,  K = K, 
                            starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                            L = L,
                            family = "normal",
                            temporal.structure = "exponential",
                            spatial.structure = "continuous",
                            seed = 29, 
                            gamma.shrinkage = TRUE,
                            include.time = TRUE,
                            include.space = TRUE,
                            clustering = TRUE,
                            seasonPeriod = 1, 
                            equalTimeDist = TRUE,
                            spatApprox = FALSE, 
                            alphaMethod = "block", 
                            h = 15, 
                            storeSpatPredPara = TRUE,
                            storeWeights = TRUE,
                            alphasWeightsToFiles = FALSE) 
save(regFixedL.simu, file="toyExfullGPfixedL.RData")
GibbsStepTimeFixedLfullGP <- regFixedL.simu$GibbsStepTime
save(GibbsStepTimeFixedLfullGP, file = "GibbsStepTimeFixedLfullGP.RData")
Diags <- diagnostics(regFixedL.simu, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags, file = "toyExfullGPfixedLDiags.RData")
Deviance <- as.mcmc(Diags$deviance)
save(Deviance, file = "toyExfullGPfixedLDeviance.RData")
spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedL, file = "toyExNNGPspatpredFixedL.RData")
temppredFixedL <- predictNewTime(regFixedL.simu, (trainingT+1):Nu, seed = 29)
save(temppredFixedL, file = "toyExNNGPtemppredFixedL.RData")
fittedClusGpMat <- matrix(0, m, 3)
clusFixedL10 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,1] <- clusFixedL10$cluster
clusFixedL100 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,2] <- clusFixedL100$cluster
clusFixedL1000 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,3] <- clusFixedL1000$cluster
weightsNumIter <- c(10, 100, 1000)
colnames(fittedClusGpMat) <- paste(weightsNumIter, "iterWeights")
save(fittedClusGpMat, file = "fullGPfixedLtoyExfittedClusGpMat.RData")
  
regFixedL.simu.block <- bfaFixedL(Y ~ 0, data = dat, dist = Dtraining, time = Time,  K = K, 
                                  starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                  L = L,
                                  family = "normal",
                                  temporal.structure = "exponential",
                                  spatial.structure = "continuous",
                                  seed = 29, 
                                  gamma.shrinkage = TRUE,
                                  include.time = TRUE,
                                  include.space = TRUE,
                                  clustering = TRUE,
                                  seasonPeriod = 1, 
                                  equalTimeDist = TRUE,
                                  spatApprox = TRUE, 
                                  alphaMethod = "block", 
                                  h = 15, 
                                  storeSpatPredPara = TRUE,
                                  storeWeights = TRUE,
                                  alphasWeightsToFiles = FALSE) 
save(regFixedL.simu.block, file="toyExNNGPblockFixedL.RData")
GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
Diags.block <- diagnostics(regFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.block, file = "toyExNNGPblockFixedLDiags.RData")
Deviance.block <- as.mcmc(Diags.block$deviance)
save(Deviance.block, file = "toyExNNGPblockFixedLDeviance.RData")
spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 9, distOrigNew, distNewNew, 
                                           NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedLblock, file = "toyExNNGPspatpredFixedLblock.RData")
temppredFixedLblock <- predictNewTime(regFixedL.simu.block, (trainingT+1):Nu, seed = 29)
save(temppredFixedLblock, file = "toyExNNGPtemppredFixedLblock.RData")
fittedClusGpMat <- matrix(0, m, 3)
clusFixedL10 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,1] <- clusFixedL10$cluster
clusFixedL100 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,2] <- clusFixedL100$cluster
clusFixedL1000 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,3] <- clusFixedL1000$cluster
weightsNumIter <- c(10, 100, 1000)
colnames(fittedClusGpMat) <- paste(weightsNumIter, "iterWeights")
save(fittedClusGpMat, file = "NNGPblockFixedLtoyExfittedClusGpMat.RData")

regFixedL.simu.sequen <- bfaFixedL(Y ~ 0, data = dat, dist = Dtraining, time = Time,  K = K, 
                                   starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                   L = L,
                                   family = "normal",
                                   temporal.structure = "exponential",
                                   spatial.structure = "continuous",
                                   seed = 29, 
                                   gamma.shrinkage = TRUE,
                                   include.time = TRUE,
                                   include.space = TRUE,
                                   clustering = TRUE,
                                   seasonPeriod = 1, 
                                   equalTimeDist = TRUE,
                                   spatApprox = TRUE, 
                                   alphaMethod = "sequential", 
                                   h = 15, 
                                   storeSpatPredPara = TRUE,
                                   storeWeights = TRUE,
                                   alphasWeightsToFiles = FALSE) 
save(regFixedL.simu.sequen, file="toyExNNGPsequenFixedL.RData")
GibbsStepTimeFixedLsequen <- regFixedL.simu.sequen$GibbsStepTime
save(GibbsStepTimeFixedLsequen, file = "GibbsStepTimeFixedLsequen.RData")
Diags.sequen <- diagnostics(regFixedL.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequen, file = "toyExNNGPsequenFixedLDiags.RData")
Deviance.sequen <- as.mcmc(Diags.sequen$deviance)
save(Deviance.sequen, file = "toyExNNGPsequenFixedLDeviance.RData")
spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 9, distOrigNew, distNewNew, 
                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedLsequen, file = "toyExNNGPspatpredFixedLsequen.RData")
temppredFixedLsequen <- predictNewTime(regFixedL.simu.sequen, (trainingT+1):Nu, seed = 29)
save(temppredFixedLsequen, file = "toyExNNGPtemppredFixedLsequen.RData")
fittedClusGpMat <- matrix(0, m, 3)
clusFixedL10 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,1] <- clusFixedL10$cluster
clusFixedL100 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,2] <- clusFixedL100$cluster
clusFixedL1000 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,3] <- clusFixedL1000$cluster
weightsNumIter <- c(10, 100, 1000)
colnames(fittedClusGpMat) <- paste(weightsNumIter, "iterWeights")
save(fittedClusGpMat, file = "NNGPsequenFixedLtoyExfittedClusGpMat.RData")
  
regVaryLj.simu.sequen <- bfaVaryingLjs(Y ~ 0, data = dat, dist = Dtraining, time = Time,  K = K,
                                       starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                       LjVec = LjVec,
                                       family = "normal",
                                       temporal.structure = "exponential",
                                       spatial.structure = "continuous",
                                       seed = 29, 
                                       gamma.shrinkage = TRUE,
                                       include.time = TRUE,
                                       include.space = TRUE,
                                       seasonPeriod = 1, 
                                       equalTimeDist = TRUE,
                                       spatApprox = TRUE, 
                                       alphaSequen = TRUE, 
                                       h = 15,
                                       storeSpatPredPara = TRUE, 
                                       storeWeights = TRUE) 
save(regVaryLj.simu.sequen, file="toyExNNGPsequenVaryLj.RData")
GibbsStepTimeVaryLjSequen <- regVaryLj.simu.sequen$GibbsStepTime
save(GibbsStepTimeVaryLjSequen, file = "GibbsStepTimeVaryLjSequen.RData")
Diags.sequenVaryLj <- diagnostics(regVaryLj.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequenVaryLj, file = "toyExNNGPsequenVaryLjDiags.RData")
Deviance.sequenVaryLj <- as.mcmc(Diags.sequenVaryLj$deviance)
save(Deviance.sequenVaryLj, file = "toyExNNGPsequenVaryLjDeviance.RData")
spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 9, distOrigNew, distNewNew, 
                                      NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredVaryLj, file = "toyExNNGPspatpredVaryLj.RData")
temppredVaryLj <- predictNewTime(regVaryLj.simu.sequen, (trainingT+1):Nu, seed = 29)
save(temppredVaryLj, file = "toyExNNGPtemppredVaryLj.RData")
fittedClusGpMat <- matrix(0, m, 3)
clusVaryLj10 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,1] <- clusVaryLj10$cluster
clusVaryLj100 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,2] <- clusVaryLj100$cluster
clusVaryLj1000 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,3] <- clusVaryLj1000$cluster
weightsNumIter <- c(10, 100, 1000)
colnames(fittedClusGpMat) <- paste(weightsNumIter, "iterWeights")
save(fittedClusGpMat, file = "NNGPsequenVaryLjtoyExfittedClusGpMat.RData")
