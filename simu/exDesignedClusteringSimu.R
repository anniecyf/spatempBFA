############################# 4 models: fullGPfixedL, NNGPblockFixedL, NNGPsequenFixedL, NNGPsequenVaryLj
#rm(list=ls())
numSpatOverallGroups <- 2
M <- 100
K <- 2
O <- 1
L <- min(60, M)
LjVec <- rep(min(60, M), K)
sqrootM <- 10
Nu <- 30
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
library(mvtnorm)
#library(dplyr)
library(fields)
#library(devtools) 
#setwd("spatTempBFA")
#load_all(".")
library(spatTempBFA)
set.seed(29)
fittedClusGpMat <- matrix(0, M, 4*3)
models <- c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj")
weightsNumIter <- c(10, 100, 1000)
ind <- expand.grid(1:3, 1:4)
colnames(fittedClusGpMat) <- paste0(models[ind[,2]], weightsNumIter[ind[,1]], "iterWeights")
  
calcGroup1 <- function(row, col, sqrootM){
  return((row - 1)*sqrootM + col)
}

print(Sys.time())
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
rowColInd <- expand.grid(2:8, 2:8)
whichGroup1 <- mapply(calcGroup1, rowColInd[,2], rowColInd[,1], sqrootM = sqrootM)
spatGroupOverall <- rep(2,M)
spatGroupOverall[whichGroup1] <- 1
#matrix(spatGroupOverall, 10, 10)
Lambda[,2] = theta2j[spatGroupOverall]
Hypers <- list(Sigma2 = list(A = 1, B = 1), Rho = list(ARho=0.1, BRho=1),
                 Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
                 Psi = list(APsi = APsi, BPsi = BPsi),
                 Upsilon = list(Zeta = K + 1, Omega = diag(K)))
MCMC <- list(NBurn = 20000, NSims = 10000, NThin = 2, NPilot = 5)
Sigma.NuMO <- rnorm(Nu * M * O, sd = sqrt(sigma2))
EtaMat <- matrix(Eta, K, Nu)
meanMat <- Lambda%*%EtaMat #M*O\times Nu
Yobs <- as.vector(meanMat) + Sigma.NuMO
dat <- data.frame(Y = Yobs)
  
regFixedL.simu <- bfaFixedL(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
                            starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
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
clusFixedL10 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,1] <- clusFixedL10$cluster
clusFixedL100 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,2] <- clusFixedL100$cluster
clusFixedL1000 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,3] <- clusFixedL1000$cluster
  
regFixedL.simu.block <- bfaFixedL(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
                                  starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
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
clusFixedLblock10 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,4] <- clusFixedLblock10$cluster
clusFixedLblock100 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,5] <- clusFixedLblock100$cluster
clusFixedLblock1000 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,6] <- clusFixedLblock1000$cluster

regFixedL.simu.sequen <- bfaFixedL(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
                                   starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
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
clusFixedLsequen10 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,7] <- clusFixedLsequen10$cluster
clusFixedLsequen100 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,8] <- clusFixedLsequen100$cluster
clusFixedLsequen1000 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,9] <- clusFixedLsequen1000$cluster
  
regVaryLj.simu.sequen <- bfaVaryingLjs(Y ~ 0, data = dat, dist = D, time = Time,  K = K,
                                       starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
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
clusVaryLj10 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
fittedClusGpMat[,10] <- clusVaryLj10$cluster
clusVaryLj100 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
fittedClusGpMat[,11] <- clusVaryLj100$cluster
clusVaryLj1000 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
fittedClusGpMat[,12] <- clusVaryLj1000$cluster
  
save(fittedClusGpMat, file = "fittedClusGpMat_exDesignedClusteringSimu.RData")

  
  


























