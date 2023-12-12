rm(list=ls())
library(mvtnorm)
library(fields)
library(spatTempBFA)
N <- 10
K <- 3
O <- 1
h <- 15
M <- 400
sqrootM <- 20
L <- min(5, M)
LjVec <- rep(min(5, M), K)
Nu <- 5
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
testingRowNum <- 5
testingLocNum <- testingRowNum^2
testingRowSeq <- seq(8.5, 12.5, by = 1)
distNewNew <- rdist(expand.grid(testingRowSeq, testingRowSeq))
distOrigNew <- rdist(expand.grid(1:sqrootM, 1:sqrootM), expand.grid(testingRowSeq, testingRowSeq))
fittingTimeMat <- spatpredTimeMat <- matrix(0, N, 4)
models <- c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj")
colnames(fittingTimeMat) <- colnames(spatpredTimeMat) <- models
postGmat <- postMSEmat <- postVarMat <- matrix(0, N, 4)
colnames(postGmat) <- colnames(postMSEmat) <- colnames(postVarMat) <- models

sigma2 <- 0.01
kappa <- 0.7
psi <- 2.3
Hpsi <- exp(-psi*TimeDist)
rho <- 0.8
D <- rdist(expand.grid(1:sqrootM, 1:sqrootM))
Mall <- M + testingLocNum
Dall <- rbind(cbind(D, distOrigNew), cbind(t(distOrigNew), distNewNew))
FrhoAll <- exp(-rho*Dall)
maxL <- 5
a1 = 1; a2 = 10
Hypers <- list(Sigma2 = list(A = 1, B = 1), 
               Rho = list(ARho = 0.1, BRho = 1),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
MCMC <- list(NBurn = 20000, NSims = 10000, NThin = 2, NPilot = 5)

for(n in 1:N){
  print(paste("n = ", n, sep = ""))
  tempMat <- matrix(runif(K*K,0,1), K, K)
  Upsilon <- t(tempMat)%*%tempMat
  Eta <- rmvnorm(1, mean = rep(0, Nu*K), sigma = kronecker(Hpsi, Upsilon))
  LStarJ <- sample(maxL, size = K, replace = TRUE)
  Alpha <- list()
  for(j in 1:K) {
    Alpha[[j]] <-  t(rmvnorm(LStarJ[j], mean=rep(0, Mall*O), sigma = kappa*FrhoAll))
    #every list index an Mall by L_j matrix since O = 1
  }
  w <- list()
  for(j in 1:K){
    w[[j]] <- pnorm(Alpha[[j]])
    Lj <- LStarJ[j]
    w[[j]][,Lj] <- rep(1, Mall)
    temp <- rep(1, Mall)
    for(l in 1:Lj){
      w[[j]][,l] <- w[[j]][,l]*temp
      if(l<Lj) {temp <- temp * pnorm(Alpha[[j]][,l], lower.tail = FALSE)}
    }
  }
  Xi <- matrix(1, Mall, K)
  for(j in 1:K){
    Lj <- LStarJ[j]
    for(i in 1:Mall){
      Xi[i,j] <- sample(Lj, size=1, prob=w[[j]][i,])
    }
  }
  Delta <- sapply(c(a1,rep(a2,(K-1))), rgamma, n=1, rate=1)
  Tau <- cumprod(Delta)
  Theta <- list()
  for(j in 1:K){
    Theta[[j]] <- rnorm(LStarJ[j], 0, sd=sqrt(1/Tau[j])) #vector of length Lj
  }
  Lambda <- matrix(0, Mall*O, K)
  for(j in 1:K){
    for(i in 1:Mall){
      Lambda[i,j] = Theta[[j]][Xi[i,j]]
    }
  }
  EtaMat <- matrix(Eta, K, Nu)
  meanMat <- Lambda%*%EtaMat #Mall*O\times Nu with O = 1
  Sigma.NuMO <- matrix(rnorm(Nu * Mall * O, sd = sqrt(sigma2)), Mall*O, Nu)
  Ymat <- meanMat + Sigma.NuMO #Mall*O\times Nu with O = 1
  Yobs <- as.vector(Ymat[1:M,])
  dat <- data.frame(Y = Yobs)
  Ytrueloc <- as.vector(Ymat[(M+1):Mall,])
  
  t1 <- Sys.time()
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
                              h = h, 
                              storeSpatPredPara = TRUE,
                              storeWeights = TRUE,
                              alphasWeightsToFiles = FALSE) 
  t2 <- Sys.time()
  spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, testingLocNum, distOrigNew, distNewNew, 
                                        NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
  t3 <- Sys.time()
  fittingTimeMat[n, 1] <- difftime(t2, t1, units = "mins")
  spatpredTimeMat[n, 1] <- difftime(t3, t2, units = "secs")
  ylocPredList = spatpredFixedL$Y
  ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * Nu * O)) #(testingLocNum x T x O) x NKeep 
  ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x T (rowMeans) note that O = 1 here
  postG <- mean((ylocPredMean - Ytrueloc)^2)
  diffMat <- sweep(ylocPred, 1, Ytrueloc, "-")
  postMSE = mean(rowMeans(diffMat^2))
  postVar = mean(apply(ylocPred, 1, var))
  postGmat[n, 1] <- postG
  postMSEmat[n, 1] <- postMSE
  postVarMat[n, 1] <- postVar
  
  t1 <- Sys.time()
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
                                    h = h,
                                    storeSpatPredPara = TRUE,
                                    storeWeights = TRUE,
                                    alphasWeightsToFiles = FALSE) 
  t2 <- Sys.time()
  spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, testingLocNum, distOrigNew, distNewNew, 
                                             NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
  t3 <- Sys.time()
  fittingTimeMat[n, 2] <- difftime(t2, t1, units = "mins")
  spatpredTimeMat[n, 2] <- difftime(t3, t2, units = "secs")
  ylocPredList = spatpredFixedLblock$Y
  ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * Nu * O)) #(testingLocNum x T x O) x NKeep 
  ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x T (rowMeans) note that O = 1 here
  postG <- mean((ylocPredMean - Ytrueloc)^2)
  diffMat <- sweep(ylocPred, 1, Ytrueloc, "-")
  postMSE = mean(rowMeans(diffMat^2))
  postVar = mean(apply(ylocPred, 1, var))
  postGmat[n, 2] <- postG
  postMSEmat[n, 2] <- postMSE
  postVarMat[n, 2] <- postVar
  
  t1 <- Sys.time()
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
                                     h = h, 
                                     storeSpatPredPara = TRUE,
                                     storeWeights = TRUE,
                                     alphasWeightsToFiles = FALSE) 
  t2 <- Sys.time()
  spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, testingLocNum, distOrigNew, distNewNew, 
                                              NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
  t3 <- Sys.time()
  fittingTimeMat[n, 3] <- difftime(t2, t1, units = "mins")
  spatpredTimeMat[n, 3] <- difftime(t3, t2, units = "secs")
  ylocPredList = spatpredFixedLsequen$Y
  ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * Nu * O)) #(testingLocNum x T x O) x NKeep 
  ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x T (rowMeans) note that O = 1 here
  postG <- mean((ylocPredMean - Ytrueloc)^2)
  diffMat <- sweep(ylocPred, 1, Ytrueloc, "-")
  postMSE = mean(rowMeans(diffMat^2))
  postVar = mean(apply(ylocPred, 1, var))
  postGmat[n, 3] <- postG
  postMSEmat[n, 3] <- postMSE
  postVarMat[n, 3] <- postVar
  
  t1 <- Sys.time()
  regVaryLj.simu.sequen <- bfaVaryingLjs(Y ~ 0, data = dat, dist = D, time = Time,  K = K, LjVec = LjVec,
                                         starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
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
                                         h = h,
                                         storeSpatPredPara = TRUE, 
                                         storeWeights = TRUE) 
  t2 <- Sys.time()
  spatpredVaryLjSequen <- predictNewLocVaryLj(regVaryLj.simu.sequen, testingLocNum, distOrigNew, distNewNew, 
                                              NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
  t3 <- Sys.time()
  #save(spatpredVaryLjSequen, file = paste0("n", n, "spatpredVaryLj5simuSequenT5M400Iter30000.RData"))
  #save(regVaryLj.simu.sequen, file = paste0("n", n, "regVaryLj5simuSequenT5M400Iter30000.RData"))
  fittingTimeMat[n, 4] <- difftime(t2, t1, units = "mins")
  spatpredTimeMat[n, 4] <- difftime(t3, t2, units = "secs")
  ylocPredList = spatpredVaryLjSequen$Y
  ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * Nu * O)) #(testingLocNum x T x O) x NKeep 
  ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x T (rowMeans) note that O = 1 here
  postG <- mean((ylocPredMean - Ytrueloc)^2)
  diffMat <- sweep(ylocPred, 1, Ytrueloc, "-")
  postMSE = mean(rowMeans(diffMat^2))
  postVar = mean(apply(ylocPred, 1, var))
  postGmat[n, 4] <- postG
  postMSEmat[n, 4] <- postMSE
  postVarMat[n, 4] <- postVar
}

save(fittingTimeMat, file = "N10M400h15T5K3fittingTimeMat_spatPredMetricSimu.RData")
save(spatpredTimeMat, file = "N10M400h15T5K3spatpredTimeMat_spatPredMetricSimu.RData")
save(postGmat, file = "N10M400h15T5K3postGmat_spatPredMetricSimu.RData")
save(postMSEmat, file = "N10M400h15T5K3postMSEmat_spatPredMetricSimu.RData")
save(postVarMat, file = "N10M400h15T5K3postVarMat_spatPredMetricSimu.RData")





