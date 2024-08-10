rm(list=ls())
library(spatempBFA)
# projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
# setwd(paste(projDirec, "simu/appen/equalTimeDistAccelerationVerification/m64T500", sep = "/"))
modelVec <- c("fullGPfixedLfast", "fullGPfixedL", "NNGPblockFixedLfast", "NNGPblockFixedL", 
              "NNGPsequenFixedLfast", "NNGPsequenFixedL", "NNGPsequenVaryLjFast", "NNGPsequenVaryLj")
N <- 100
Nu <- 500
newT <- 10
temppredTimeMat <- matrix(0, N, 8)
colnames(temppredTimeMat) <- modelVec
load("regFixedL30simuT500M64Iter30000.RData")
regFixedL.simu.fast <- regFixedL.simu
load("regFixedL30simuT500M64Iter30000specifyEqualTimeDistF.RData")
load("regFixedL30simuBlockT500M64Iter30000.RData")
regFixedL.simu.block.fast <- regFixedL.simu.block
load("regFixedL30simuBlockT500M64Iter30000specifyEqualTimeDistF.RData")
load("regFixedL30simuSequenT500M64Iter30000nostorealphaweights.RData")
regFixedL.simu.sequen.fast <- regFixedL.simu.sequen
load("regFixedL30simuSequenT500M64Iter30000specifyEqualTimeDistFnostorealphaweights.RData")
load("regVaryLjsimuSequenT500M64Iter30000nostorealphaweight.RData")
regVaryLj.simu.sequen.fast <- regVaryLj.simu.sequen
load("regVaryLjsimuSequenT500M64Iter30000specifyEqualTimeDistFnostorealphaweight.RData")
for (n in 1:N) {
  print(n)
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regFixedL.simu.fast, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 1] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regFixedL.simu, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 2] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regFixedL.simu.block.fast, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 3] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regFixedL.simu.block, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 4] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regFixedL.simu.sequen.fast, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 5] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regFixedL.simu.sequen, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 6] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regVaryLj.simu.sequen.fast, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 7] = difftime(t2, t1, units = "secs")
  t1 <- Sys.time()
  temppredFixedL <- predictNewTime(regVaryLj.simu.sequen, (Nu+1):(Nu+newT), seed = 29)
  t2 <- Sys.time()
  temppredTimeMat[n, 8] = difftime(t2, t1, units = "secs")
}
save(temppredTimeMat, file = "m64T500temppredTimeMat.RData")
