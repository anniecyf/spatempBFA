############################# 4 models: fullGPfixedL, NNGPblockFixedL, NNGPsequenFixedL, NNGPsequenVaryLj
#rm(list=ls())
numSpatOverallGroups <- 2
N <- 18
M <- 100
K <- 2
O <- 1
L <- min(10, M)
LjVec <- rep(min(10, M), K)
sqrootM <- 10
Nu <- 30
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
library(mvtnorm)
library(fields)
library(spatempBFA)
set.seed(29)
nkeep10RandIndex <- nkeep100RandIndex <- nkeep1000RandIndex <- matrix(0, N, 4)
models <- c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj")
colnames(nkeep10RandIndex) <- colnames(nkeep100RandIndex) <- colnames(nkeep1000RandIndex) <- models
nkeep10AccuRatio <- nkeep100AccuRatio <- nkeep1000AccuRatio <- matrix(0, N, 4)
colnames(nkeep10AccuRatio) <- colnames(nkeep100AccuRatio) <- colnames(nkeep1000AccuRatio) <- models
calcGroup1 <- function(row, col, sqrootM){
  return((row - 1)*sqrootM + col)
}
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
  #actualGroup is spatGroupOverall; numObs is M
  actualGroupsPoss1 = factor(actualGroup, labels=c(1,2))
  actualGroupsPoss2 = factor(actualGroup, labels=c(2,1))
  accuRatio <- max(sum(predictedCluster==actualGroupsPoss1), sum(predictedCluster==actualGroupsPoss2))/numObs 
  return(accuRatio)
}

for(n in 1:N){
  print(Sys.time())
  print(paste("n = ", n, sep = ""))
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
  
  regFixedL.simu <- VAR1bfaFixedL(Y ~ 0, data = dat, dist = D, Nu = Nu,  K = K, 
                                  starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                  L = L,
                                  family = "normal",
                                  spatial.structure = "continuous",
                                  seed = 29, 
                                  gamma.shrinkage = TRUE,
                                  include.space = TRUE,
                                  clustering = TRUE,
                                  spatApprox = FALSE, 
                                  alphaMethod = "block", 
                                  h = 15, 
                                  storeSpatPredPara = TRUE,
                                  storeWeights = TRUE,
                                  alphasWeightsToFiles = FALSE) 
  #save(regFixedL.simu, file = paste("regFixedL.simu", n, "RData", sep = "."))
  clusFixedL10 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
  clusFixedL100 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
  clusFixedL1000 <- clusteringFixedL(regFixedL.simu, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
  nkeep10RandIndex[n, 1] <- calcRandIndex(predictedCluster = clusFixedL10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100RandIndex[n, 1] <- calcRandIndex(predictedCluster = clusFixedL100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000RandIndex[n, 1] <- calcRandIndex(predictedCluster = clusFixedL1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep10AccuRatio[n, 1] <- calcAccuRatio(predictedCluster = clusFixedL10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100AccuRatio[n, 1] <- calcAccuRatio(predictedCluster = clusFixedL100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000AccuRatio[n, 1] <- calcAccuRatio(predictedCluster = clusFixedL1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  
  regFixedL.simu.block <- VAR1bfaFixedL(Y ~ 0, data = dat, dist = D, Nu = Nu,  K = K, 
                                        starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                        L = L,
                                        family = "normal",
                                        spatial.structure = "continuous",
                                        seed = 29, 
                                        gamma.shrinkage = TRUE,
                                        include.space = TRUE,
                                        clustering = TRUE,
                                        spatApprox = TRUE, 
                                        alphaMethod = "block", 
                                        h = 15, 
                                        storeSpatPredPara = TRUE,
                                        storeWeights = TRUE,
                                        alphasWeightsToFiles = FALSE) 
  #save(regFixedL.simu.block, file = paste("regFixedL.simu.block", n, "RData", sep = "."))
  clusFixedLblock10 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
  clusFixedLblock100 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
  clusFixedLblock1000 <- clusteringFixedL(regFixedL.simu.block, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
  nkeep10RandIndex[n, 2] <- calcRandIndex(predictedCluster = clusFixedLblock10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100RandIndex[n, 2] <- calcRandIndex(predictedCluster = clusFixedLblock100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000RandIndex[n, 2] <- calcRandIndex(predictedCluster = clusFixedLblock1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep10AccuRatio[n, 2] <- calcAccuRatio(predictedCluster = clusFixedLblock10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100AccuRatio[n, 2] <- calcAccuRatio(predictedCluster = clusFixedLblock100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000AccuRatio[n, 2] <- calcAccuRatio(predictedCluster = clusFixedLblock1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  
  regFixedL.simu.sequen <- VAR1bfaFixedL(Y ~ 0, data = dat, dist = D, Nu = Nu,  K = K,
                                         starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                         L = L,
                                         family = "normal",
                                         spatial.structure = "continuous",
                                         seed = 29, 
                                         gamma.shrinkage = TRUE,
                                         include.space = TRUE,
                                         clustering = TRUE,
                                         spatApprox = TRUE, 
                                         alphaMethod = "sequential", 
                                         h = 15, 
                                         storeSpatPredPara = TRUE,
                                         storeWeights = TRUE,
                                         alphasWeightsToFiles = FALSE) 
  #save(regFixedL.simu.sequen, file = paste("regFixedL.simu.sequen", n, "RData", sep = "."))
  clusFixedLsequen10 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
  clusFixedLsequen100 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
  clusFixedLsequen1000 <- clusteringFixedL(regFixedL.simu.sequen, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
  nkeep10RandIndex[n, 3] <- calcRandIndex(predictedCluster = clusFixedLsequen10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100RandIndex[n, 3] <- calcRandIndex(predictedCluster = clusFixedLsequen100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000RandIndex[n, 3] <- calcRandIndex(predictedCluster = clusFixedLsequen1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep10AccuRatio[n, 3] <- calcAccuRatio(predictedCluster = clusFixedLsequen10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100AccuRatio[n, 3] <- calcAccuRatio(predictedCluster = clusFixedLsequen100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000AccuRatio[n, 3] <- calcAccuRatio(predictedCluster = clusFixedLsequen1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  
  regVaryLj.simu.sequen <- VAR1bfaVaryingLjs(Y ~ 0, data = dat, dist = D, Nu = Nu,  K = K,
                                             starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                             LjVec = LjVec,
                                             family = "normal",
                                             spatial.structure = "continuous",
                                             seed = 29, 
                                             gamma.shrinkage = TRUE,
                                             include.space = TRUE,
                                             spatApprox = TRUE, 
                                             alphaSequen = TRUE, 
                                             h = 15,
                                             storeSpatPredPara = TRUE, 
                                             storeWeights = TRUE)
  #save(regVaryLj.simu.sequen, file = paste("regVaryLj.simu.sequen", n, "RData", sep = "."))
  clusVaryLj10 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 10, nCent = numSpatOverallGroups)
  clusVaryLj100 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 100, nCent = numSpatOverallGroups)
  clusVaryLj1000 <- clusteringVaryLj(regVaryLj.simu.sequen, o = 1, nkeep = 1000, nCent = numSpatOverallGroups)
  nkeep10RandIndex[n, 4] <- calcRandIndex(predictedCluster = clusVaryLj10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100RandIndex[n, 4] <- calcRandIndex(predictedCluster = clusVaryLj100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000RandIndex[n, 4] <- calcRandIndex(predictedCluster = clusVaryLj1000$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep10AccuRatio[n, 4] <- calcAccuRatio(predictedCluster = clusVaryLj10$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep100AccuRatio[n, 4] <- calcAccuRatio(predictedCluster = clusVaryLj100$cluster, actualGroup = spatGroupOverall, numObs = M)
  nkeep1000AccuRatio[n, 4] <- calcAccuRatio(predictedCluster = clusVaryLj1000$cluster, actualGroup = spatGroupOverall, numObs = M)
}
save(nkeep10RandIndex, file = "VAR1nkeep10RandIndexMat_designedClustering.RData")
save(nkeep100RandIndex, file = "VAR1nkeep100RandIndexMat_designedClustering.RData")
save(nkeep1000RandIndex, file = "VAR1nkeep1000RandIndexMat_designedClustering.RData")
save(nkeep10AccuRatio, file = "VAR1nkeep10AccuRatioMat_designedClustering.RData")
save(nkeep100AccuRatio, file = "VAR1nkeep100AccuRatioMat_designedClustering.RData")
save(nkeep1000AccuRatio, file = "VAR1nkeep1000AccuRatioMat_designedClustering.RData")
