rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/mainScalabilityVerificationSimu/toy example", sep = "/"))
setwd("spBFA")
# library(spBFA)
# Nu = 60; trainingT = 50
# load("toyExspBFAL10.RData"); reg.simu$runtime # 6.72 hours
# temppredspBFAL10 <- predict(reg.simu, (trainingT+1):Nu, seed = 29) # from spBFA
# #save(temppredspBFAL10, file = "toyExspBFAL10temppred.RData")
# rm(reg.simu)
# load("toyExspBFALInf.RData"); reg.simu.LInf$runtime # 8.21 hours
# temppredspBFALInf <- predict(reg.simu.LInf, (trainingT+1):Nu, seed = 29)
# #save(temppredspBFALInf, file = "toyExspBFALInftemppred.RData")
# rm(reg.simu.LInf)
load("toyExspBFAL10Diags.RData")
load("toyExspBFAL10Deviance.RData")
load("toyExspBFAL10temppred.RData")
load("toyExspBFALInfDiags.RData")
load("toyExspBFALInfDeviance.RData")
load("toyExspBFALInftemppred.RData")
spBFADiags <- Diags; spBFADeviance <- Deviance; spBFADiagsLInf <- Diags.LInf; spBFADevianceLInf <- Deviance.LInf
rm(Diags, Diags.LInf, Deviance, Deviance.LInf)
setwd("../fullGPfixedL")
# load("toyExfullGPfixedL.RData"); regFixedL.simu$runtime # longer package version 1.71 hours
# rm(regFixedL.simu)
load("toyExfullGPfixedLDiags.RData")
load("toyExfullGPfixedLDeviance.RData")
load("GibbsStepTimeFixedLfullGP.RData")
load("toyExFullGPtemppredFixedL.RData")
load("toyExFullGPspatpredFixedL.RData")
load("fullGPfixedLtoyExfittedClusGpMat.RData")
fittedClusGpMat.fullGPfixedL <- fittedClusGpMat
setwd("../NNGPblockFixedL")
# load("toyExNNGPblockFixedL.RData"); regFixedL.simu.block$runtime # longer package version 1.7 hours
# rm(regFixedL.simu.block)
load("toyExNNGPblockFixedLDiags.RData")
load("toyExNNGPblockFixedLDeviance.RData")
load("GibbsStepTimeFixedLblock.RData")
load("toyExNNGPtemppredFixedLblock.RData")
load("toyExNNGPspatpredFixedLblock.RData")
load("NNGPblockFixedLtoyExfittedClusGpMat.RData")
fittedClusGpMat.NNGPblockFixedL <- fittedClusGpMat
setwd("../NNGPsequenFixedL")
# load("toyExNNGPsequenFixedL.RData"); regFixedL.simu.sequen$runtime # longer package version 1.38 hours
# rm(regFixedL.simu.sequen)
load("toyExNNGPsequenFixedLDiags.RData")
load("toyExNNGPsequenFixedLDeviance.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("toyExNNGPtemppredFixedLsequen.RData")
load("toyExNNGPspatpredFixedLsequen.RData")
load("NNGPsequenFixedLtoyExfittedClusGpMat.RData")
fittedClusGpMat.NNGPsequenFixedL <- fittedClusGpMat
setwd("../NNGPsequenVaryLj")
# load("toyExNNGPsequenVaryLj.RData"); regVaryLj.simu.sequen$runtime # longer package version 1.41 hours
# rm(regVaryLj.simu.sequen)
load("toyExNNGPsequenVaryLjDiags.RData")
load("toyExNNGPsequenVaryLjDeviance.RData")
load("GibbsStepTimeVaryLjsequen.RData")
load("toyExNNGPtemppredVaryLj.RData")
load("toyExNNGPspatpredVaryLj.RData")
load("NNGPsequenVaryLjtoyExfittedClusGpMat.RData")
fittedClusGpMat.NNGPsequenVaryLj <- fittedClusGpMat
rm(fittedClusGpMat)
setwd("..")
NKeep <- 5000
postDeviancesDF <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance), NNGPblockFixedL = as.vector(Deviance.block),
                              NNGPsequenFixedL = as.vector(Deviance.sequen), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj),
                              spBFAL10 = as.vector(spBFADeviance), spBFALInf = as.vector(spBFADevianceLInf))
fullGPfixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances
NNGPblockFixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances
NNGPsequenFixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances
NNGPsequenVaryLjpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances
ggarrange(fullGPfixedLpostDeviances, NNGPblockFixedLpostDeviances, 
          NNGPsequenFixedLpostDeviances, NNGPsequenVaryLjpostDeviances,
          labels = c("A", "B", "C", "D"), align = "hv",
          ncol = 2, nrow = 2)
#ggsave("toyExPostDeviances.png", width = 16, height = 16, units = "cm")
spBFAL10Deviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = spBFAL10)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
spBFAL10Deviances
spBFALInfDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = spBFALInf)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
spBFALInfDeviances
ggarrange(fullGPfixedLpostDeviances, NNGPblockFixedLpostDeviances, 
          NNGPsequenFixedLpostDeviances, NNGPsequenVaryLjpostDeviances,
          spBFAL10Deviances, spBFALInfDeviances,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 2, nrow = 3)
#ggsave("toyExPostDeviances.png", width = 16, height = 24, units = "cm")
ggarrange(fullGPfixedLpostDeviances, NNGPsequenFixedLpostDeviances, spBFAL10Deviances,
          NNGPblockFixedLpostDeviances, NNGPsequenVaryLjpostDeviances, spBFALInfDeviances,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2)
#ggsave("toyExPostDeviancesSmall.png", width = 24, height = 16, units = "cm")
apply(GibbsStepTimeFixedLfullGP, 2, summary)
apply(GibbsStepTimeFixedLblock, 2, summary)
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)
# View(Diags)
# View(Diags.block)
# View(Diags.sequen)
# View(Diags.sequenVaryLj)
# View(spBFADiags)
# View(spBFADiagsLInf)
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
  discardInd <- c(discardInd, seq(from = whichTestingLoc, by = M, length.out = trainingT))
}
Ytraining <- YtrainingTemp[-discardInd]
YtestingSpat <- YtrainingTemp[discardInd]  #YtestingSpat <- as.vector(t(matrix(YtestingSpat, M-m, trainingT)))
rm(YtrainingTemp)
spatGroupOverallTraining <- spatGroupOverall[-whichTesting]
distOrigNew = as.matrix(D[-whichTesting, whichTesting])
distNewNew = as.matrix(D[whichTesting, whichTesting])
dat = data.frame(Y = Ytraining); dist = Dtraining
tempTestingDiscardInd <- vector()
for(whichTestingLoc in whichTesting){
  tempTestingDiscardInd <- c(tempTestingDiscardInd, seq(from = whichTestingLoc, by = M, length.out = testingT))
}
YtestingTemp <- YtestingTemp[-tempTestingDiscardInd]
xcoord <- rep(1:sqrootM, sqrootM)
ycoord <- rep(seq(sqrootM, 1, by = -1), each = sqrootM)
trainingTestingSpatGp <- rep(1, M); trainingTestingSpatGp[whichTesting] = 2
spatGpDF <- data.frame(x = xcoord, y = ycoord, actualGp <- as.factor(spatGroupOverall), 
                       trainingTestingSpatGp <- as.factor(trainingTestingSpatGp))
trainingTestingPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = trainingTestingSpatGp), size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
trainingTestingPlot
#ggsave("toyExSpatTrainingTestingPlot.png", width=10, height=10, units="cm")
trainingTestingPlotBW <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = trainingTestingSpatGp), size = 7, show.legend = FALSE) + 
  scale_color_manual(values = c("white", "black")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
trainingTestingPlotBW
#ggsave("toyExSpatTrainingTestingPlotBW.png", width=10, height=10, units="cm")
actualSpatPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = actualGp), size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
#actualSpatPlot
#ggsave("toyExActualSpatCluster.png", width=10, height=10, units="cm")
actualSpatPlotBW <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = actualGp), size = 7, show.legend = FALSE) + 
  scale_color_manual(values = c("black", "white")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
actualSpatPlotBW
#ggsave("toyExActualSpatClusterBW.png", width=10, height=10, units="cm")
calcRandIndex <- function(predictedCluster, actualGroup, numObs){
  denominator <- numObs * (numObs-1) / 2
  numerator <- denominator
  for (i in 1:(numObs-1)){
    for (j in (i + 1):numObs){
      if ( (predictedCluster[i]==predictedCluster[j]) + (actualGroup[i]==actualGroup[j]) == 1) numerator = numerator - 1
    }
  }
  randIndex <- numerator / denominator
  return(randIndex)
}
calcAccuRatio <- function(predictedCluster, actualGroup, numObs){
  #actualGroup is spatGroupOverall; numObs is M;       kmeans3 = kmeans(X, centers = 3); cl3 = kmeans3$cluster
  #predictedCluster is kmeans2$cluster
  actualGroupsPoss1 = factor(actualGroup, labels=c(1,2))
  actualGroupsPoss2 = factor(actualGroup, labels=c(2,1))
  accuRatio <- max(sum(predictedCluster==actualGroupsPoss1), sum(predictedCluster==actualGroupsPoss2))/numObs 
  return(accuRatio)
}
fittedClusGpMat <- cbind(fittedClusGpMat.fullGPfixedL, fittedClusGpMat.NNGPblockFixedL,
                         fittedClusGpMat.NNGPsequenFixedL, fittedClusGpMat.NNGPsequenVaryLj)
# fittedClusGpMatList <- lapply(apply(fittedClusGpMat, 2, list), function(x) x[[1]])
# sapply(fittedClusGpMatList, calcAccuRatio, actualGroup = spatGroupOverallTraining, numObs = m)
apply(fittedClusGpMat, 2, calcAccuRatio, actualGroup = spatGroupOverallTraining, numObs = m)
apply(fittedClusGpMat, 2, calcRandIndex, actualGroup = spatGroupOverallTraining, numObs = m)
# calcAccuRatio(fittedClusGpMat.fullGPfixedL[,1], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.fullGPfixedL[,2], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.fullGPfixedL[,3], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPblockFixedL[,1], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPblockFixedL[,2], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPblockFixedL[,3], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPsequenFixedL[,1], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPsequenFixedL[,2], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPsequenFixedL[,3], spatGroupOverallTraining, m) # 1
# calcAccuRatio(fittedClusGpMat.NNGPsequenVaryLj[,1], spatGroupOverallTraining, m) # 0.9857955
# calcAccuRatio(fittedClusGpMat.NNGPsequenVaryLj[,2], spatGroupOverallTraining, m) # 0.9886364
# calcAccuRatio(fittedClusGpMat.NNGPsequenVaryLj[,3], spatGroupOverallTraining, m) # 0.9886364
# calcRandIndex(fittedClusGpMat.fullGPfixedL[,1], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.fullGPfixedL[,2], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.fullGPfixedL[,3], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPblockFixedL[,1], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPblockFixedL[,2], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPblockFixedL[,3], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPsequenFixedL[,1], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPsequenFixedL[,2], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPsequenFixedL[,3], spatGroupOverallTraining, m) # 1
# calcRandIndex(fittedClusGpMat.NNGPsequenVaryLj[,1], spatGroupOverallTraining, m) # 0.9719147
# calcRandIndex(fittedClusGpMat.NNGPsequenVaryLj[,2], spatGroupOverallTraining, m) # 0.977467
# calcRandIndex(fittedClusGpMat.NNGPsequenVaryLj[,3], spatGroupOverallTraining, m) # 0.977467
xcoordTraining <- xcoord[-whichTesting]; ycoordTraining = ycoord[-whichTesting]
fittedNNGPsequenVaryLjspatGpDF <- data.frame(x = xcoordTraining, y = ycoordTraining, 
                                             nkeep10 = as.factor(fittedClusGpMat.NNGPsequenVaryLj[,1]),
                                             nkeep100 = as.factor(fittedClusGpMat.NNGPsequenVaryLj[,2]),
                                             nkeep1000 = as.factor(fittedClusGpMat.NNGPsequenVaryLj[,3]))
# nkeep10NNGPsequenVaryLj <- ggplot(fittedNNGPsequenVaryLjspatGpDF) + 
#   geom_point(aes(x = x, y = y, col = nkeep10), size = 7, show.legend = FALSE) + 
#   theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
# nkeep10NNGPsequenVaryLj
# ggsave("toyExnkeep10NNGPsequenVaryLj.png", width=10, height=10, units="cm")
nkeep10NNGPsequenVaryLjBW <- ggplot(fittedNNGPsequenVaryLjspatGpDF) + 
  geom_point(aes(x = x, y = y, col = nkeep10), size = 7, show.legend = FALSE) +  
  scale_color_manual(values = c("white", "black")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
nkeep10NNGPsequenVaryLjBW # 5 locations clustered wrongly (0.9857955*352 = 347)
nkeep100NNGPsequenVaryLjBW <- ggplot(fittedNNGPsequenVaryLjspatGpDF) + 
  geom_point(aes(x = x, y = y, col = nkeep100), size = 7, show.legend = FALSE) +  
  scale_color_manual(values = c("black", "white")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
nkeep100NNGPsequenVaryLjBW # 4 locations clustered wrongly (0.9886364*352 = 348)
nkeep1000NNGPsequenVaryLjBW <- ggplot(fittedNNGPsequenVaryLjspatGpDF) + 
  geom_point(aes(x = x, y = y, col = nkeep1000), size = 7, show.legend = FALSE) +  
  scale_color_manual(values = c("white", "black")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
nkeep1000NNGPsequenVaryLjBW # 4 locations clustered wrongly (0.9886364*352 = 348)
ggarrange(nkeep10NNGPsequenVaryLjBW, nkeep100NNGPsequenVaryLjBW, nkeep1000NNGPsequenVaryLjBW,
          labels = c("A", "B", "C"), align = "h", ncol = 3, nrow = 1)
#ggsave("toyExNNGPsequenVaryLjBWspatPlots.png", width = 30, height = 10, units="cm")
nkeep10agreements <- spatGroupOverallTraining + fittedClusGpMat.NNGPsequenVaryLj[,1] == 3
nkeep100agreements <- spatGroupOverallTraining == fittedClusGpMat.NNGPsequenVaryLj[,2]
nkeep1000agreements <- spatGroupOverallTraining + fittedClusGpMat.NNGPsequenVaryLj[,3] == 3
disaDFnkeep10 <- data.frame(x = xcoordTraining[!nkeep10agreements], y = ycoordTraining[!nkeep10agreements])
nkeep10NNGPsequenVaryLj <- nkeep10NNGPsequenVaryLjBW + 
  geom_point(data = disaDFnkeep10, aes(x = x, y = y), col = "red", fill = NA, shape = 0, size = 7.5)
nkeep10NNGPsequenVaryLj
disaDFnkeep100 <- data.frame(x = xcoordTraining[!nkeep100agreements], y = ycoordTraining[!nkeep100agreements])
nkeep100NNGPsequenVaryLj <- nkeep100NNGPsequenVaryLjBW + 
  geom_point(data = disaDFnkeep100, aes(x = x, y = y), col = "red", fill = NA, shape = 0, size = 7.5)
nkeep100NNGPsequenVaryLj
disaDFnkeep1000 <- data.frame(x = xcoordTraining[!nkeep1000agreements], y = ycoordTraining[!nkeep1000agreements])
nkeep1000NNGPsequenVaryLj <- nkeep1000NNGPsequenVaryLjBW + 
  geom_point(data = disaDFnkeep1000, aes(x = x, y = y), col = "red", fill = NA, shape = 0, size = 7.5)
nkeep1000NNGPsequenVaryLj
ggarrange(nkeep10NNGPsequenVaryLj, nkeep100NNGPsequenVaryLj, nkeep1000NNGPsequenVaryLj,
          labels = c("A", "B", "C"), align = "h", ncol = 3, nrow = 1)
#ggsave("toyExNNGPsequenVaryLjspatPlots.png", width = 30, height = 10, units="cm")
summary(YtestingSpat)
summary(abs(YtestingSpat))
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

summary(YtestingTemp)
summary(abs(YtestingTemp))
tempPredMetricMat <- matrix(0, 3, 6, 
                            dimnames = list(c("postMeanMSE", "postMSE", "postVar"),
                                            c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj",
                                              "spBFAL10", "spBFALInf")))
# testingT = 10
ytempPredList = temppredFixedL$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 1] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 1] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 1] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredFixedLblock$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 2] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 2] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 2] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredFixedLsequen$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 3] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 3] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 3] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredVaryLj$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 4] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 4] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 4] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredspBFAL10$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 5] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 5] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 5] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredspBFALInf$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 6] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 6] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 6] = mean(apply(ytempPred, 1, var))
tempPredMetricMat