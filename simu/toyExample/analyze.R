rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/mainScalabilityVerificationSimu/toyExample", sep = "/"))
setwd("spBFA")
# load("toyExspBFAL10.RData"); reg.simu$runtime #  7.24 days
# rm(reg.simu)
# load("toyExspBFALInf.RData"); reg.simu.LInf$runtime #  3.08 days
# rm(reg.simu.LInf)
## use the spatempBFA package version with the DIAG_diagnostics.R 3 lines if ( ! ( is.FixedLbfa(object) | is.varyLjBFA(object) | is.FixedLbfaVAR1(object) | is.VAR1varyLjBFA(object) ) ) {
#stop('"object" must be of class FixedLbfa or varyLjBFA or FixedLbfaVAR1 or VAR1varyLjBFA')
#} commented out
# spBFADiagsComplete <- diagnostics(reg.simu)
# save(spBFADiagsComplete, file = "toyExspBFAL10DiagsComplete.RData")
# spBFADiagsLInfComplete <- diagnostics(reg.simu.LInf)
# save(spBFADiagsLInfComplete, file = "toyExspBFALInfDiagsComplete.RData")
load("toyExspBFAL10Diags.RData")
load("toyExspBFAL10Deviance.RData")
load("toyExspBFAL10temppred.RData")
load("toyExspBFALInfDiags.RData")
load("toyExspBFALInfDeviance.RData")
load("toyExspBFALInftemppred.RData")
spBFADiags <- Diags; spBFADeviance <- Deviance; spBFADiagsLInf <- Diags.LInf; spBFADevianceLInf <- Deviance.LInf
rm(Diags, Diags.LInf, Deviance, Deviance.LInf)
setwd("../fullGPfixedL")
# load("toyExfullGPfixedL.RData"); regFixedL.simu$runtime # 15.77 hours
# rm(regFixedL.simu)
load("toyExfullGPfixedLDiags.RData")
load("toyExfullGPfixedLDeviance.RData")
load("GibbsStepTimeFixedLfullGP.RData")
load("toyExFullGPtemppredFixedL.RData")
load("toyExFullGPspatpredFixedL.RData")
load("fullGPfixedLtoyExfittedClusGpMat.RData")
fittedClusGpMat.fullGPfixedL <- fittedClusGpMat
setwd("../NNGPblockFixedL")
# load("toyExNNGPblockFixedL.RData"); regFixedL.simu.block$runtime # 15.46 hours
# rm(regFixedL.simu.block)
load("toyExNNGPblockFixedLDiags.RData")
load("toyExNNGPblockFixedLDeviance.RData")
load("GibbsStepTimeFixedLblock.RData")
load("toyExNNGPtemppredFixedLblock.RData")
load("toyExNNGPspatpredFixedLblock.RData")
load("NNGPblockFixedLtoyExfittedClusGpMat.RData")
fittedClusGpMat.NNGPblockFixedL <- fittedClusGpMat
setwd("../NNGPsequenFixedL")
# load("toyExNNGPsequenFixedL.RData"); regFixedL.simu.sequen$runtime # 15.19 hours
# rm(regFixedL.simu.sequen)
load("toyExNNGPsequenFixedLDiags.RData")
load("toyExNNGPsequenFixedLDeviance.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("toyExNNGPtemppredFixedLsequen.RData")
load("toyExNNGPspatpredFixedLsequen.RData")
load("NNGPsequenFixedLtoyExfittedClusGpMat.RData")
fittedClusGpMat.NNGPsequenFixedL <- fittedClusGpMat
setwd("../NNGPsequenVaryLj")
# load("toyExNNGPsequenVaryLj.RData"); regVaryLj.simu.sequen$runtime # 14.58 hours
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
ggsave("toyExPostDeviances.png", width = 16, height = 24, units = "cm")
ggarrange(fullGPfixedLpostDeviances, NNGPsequenFixedLpostDeviances, spBFAL10Deviances,
          NNGPblockFixedLpostDeviances, NNGPsequenVaryLjpostDeviances, spBFALInfDeviances,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2)
ggsave("toyExPostDeviancesSmall.png", width = 24, height = 16, units = "cm")
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
library(devtools)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/server/spatempBFA")
load_all(".")
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
  discardInd <- c(discardInd, seq(from = whichTestingLoc, by = M, length.out = trainingT))
}
Ytraining <- YtrainingTemp[-discardInd]
YtestingSpat <- YtrainingTemp[discardInd]  
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
trainingTestingPlotBW <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = trainingTestingSpatGp), size = 7, show.legend = FALSE) + 
  scale_color_manual(values = c("white", "black")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
trainingTestingPlotBW
#ggsave("toyExSpatTrainingTestingPlotBW.png", width=10, height=10, units="cm")
actualSpatPlotBW <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = actualGp), size = 7, show.legend = FALSE) + 
  scale_color_manual(values = c("black", "white")) + labs(x = "", y = "") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), panel.grid.minor = element_line(color="grey"))
actualSpatPlotBW
#ggsave("toyExActualSpatClusterBW.png", width=10, height=10, units="cm")
ggarrange(trainingTestingPlotBW, actualSpatPlotBW, align = "h", ncol = 2, nrow = 1, 
          labels = c("A", "B"))
#ggsave("toyExBW.png", width=24, height=12, units="cm")
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
  actualGroupsPoss1 = factor(actualGroup, labels=c(1,2))
  actualGroupsPoss2 = factor(actualGroup, labels=c(2,1))
  accuRatio <- max(sum(predictedCluster==actualGroupsPoss1), sum(predictedCluster==actualGroupsPoss2))/numObs 
  return(accuRatio)
}
fittedClusGpMat <- cbind(fittedClusGpMat.fullGPfixedL, fittedClusGpMat.NNGPblockFixedL,
                         fittedClusGpMat.NNGPsequenFixedL, fittedClusGpMat.NNGPsequenVaryLj)
# fittedClusGpMatList <- lapply(apply(fittedClusGpMat, 2, list), function(x) x[[1]])
# sapply(fittedClusGpMatList, calcAccuRatio, actualGroup = spatGroupOverallTraining, numObs = m)
# sapply(fittedClusGpMatList, calcRandIndex, actualGroup = spatGroupOverallTraining, numObs = m)
apply(fittedClusGpMat, 2, calcAccuRatio, actualGroup = spatGroupOverallTraining, numObs = m) # all 1's
apply(fittedClusGpMat, 2, calcRandIndex, actualGroup = spatGroupOverallTraining, numObs = m) # all 1's
xcoordTraining <- xcoord[-whichTesting]; ycoordTraining = ycoord[-whichTesting]
spatTrainingGpDF <- data.frame(x = xcoordTraining, y = ycoordTraining, 
                               actualTrainingSpatGp = as.factor(spatGroupOverallTraining))
actualTrainingSpatPlotBW <- ggplot(spatTrainingGpDF) + 
  geom_point(aes(x = x, y = y, col = actualTrainingSpatGp), 
             size = 7, show.legend = FALSE) + 
  scale_color_manual(values = c("black", "white")) + labs(x = "", y = "") +
  ggtitle("2 Actual Spatial Groups for the 352 Training Locations") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color="grey"), 
        panel.grid.minor = element_line(color="grey"))
actualTrainingSpatPlotBW
#ggsave("actualTrainingSpatPlotBW.png", width = 12, height = 12, units="cm")

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

