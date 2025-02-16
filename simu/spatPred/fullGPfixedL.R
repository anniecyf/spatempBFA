rm(list=ls())
library(mvtnorm)
library(fields)
library(spatempBFA)
N1 <- 3 # number of samples
seeds <- c(11, 13, 17, 19, 23, 29, 31, 37, 41)
nSeed <- length(seeds)
n <- 50
N2 <- nSeed * n
K <- 4
O <- 1
h <- 15
M <- 400
sqrootM <- 20
L <- min(40, M)
LjVec <- rep(min(40, M), K)
Nu <- 300
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
testingRowNum <- 3
testingLocNum <- testingRowNum^2
testingRowSeq <- seq(9.5, 11.5, by = 1)
distNewNew <- rdist(expand.grid(testingRowSeq, testingRowSeq))
distOrigNew <- rdist(expand.grid(1:sqrootM, 1:sqrootM), expand.grid(testingRowSeq, testingRowSeq))
spatpredTimeVec <- alphaKrigTimeVec <- 
  weightsXiLambdaKrigTimeVec <- rep(0, N1 * N2)
postGmat <- postMSEmat <- postVarMat <- 
  matrix(0, nSeed, N1, dimnames = list(paste0("seed", seeds), paste0("sample", 1:N1)))
set.seed(29)
sigma2 <- 0.01
kappa <- 0.7
psi <- 2.3
Hpsi <- exp(-psi*TimeDist)
rho <- 0.8
D <- rdist(expand.grid(1:sqrootM, 1:sqrootM))
Mall <- M + testingLocNum
Dall <- rbind(cbind(D, distOrigNew), cbind(t(distOrigNew), distNewNew))
FrhoAll <- exp(-rho*Dall)
maxL <- 10
a1 = 1; a2 = 10
Hypers <- list(Sigma2 = list(A = 1, B = 1), 
               Rho = list(ARho = 0.1, BRho = 1),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
MCMC <- list(NBurn = 30000, NSims = 10000, NThin = 10, NPilot = 5)
count = 1
for(n1 in 1:N1){
  print(paste("n1 = ", n1, sep = ""))
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
  #Tau <- cumprod(Delta)
  Theta <- list()
  for(j in 1:K){
    Theta[[j]] <- rnorm(LStarJ[j], 0, sd=sqrt(1/Delta[j])) #vector of length Lj
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
  Ytrueloc <- as.vector(t(Ymat[(M+1):Mall,]))
  
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
  save(regFixedL.simu, file = paste0("fullGPfixedL", n1, ".RData"))
  for (nseed in 1:nSeed) {
    seed = seeds[nseed]
    for (n2 in 1:n) {
      t1 <- Sys.time()
      spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, testingLocNum, distOrigNew, distNewNew, 
                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = seed)
      t2 <- Sys.time()
      # save(spatpredFixedL, file = paste0("fullGPfixedL", n1, "spatpred", n2, "seed", seed, ".RData"))
      spatpredTimeVec[count] = difftime(t2, t1, units = "secs")
      alphaKrigTimeVec[count] = spatpredFixedL$alphaKrigTime
      weightsXiLambdaKrigTimeVec[count] = spatpredFixedL$weightsXiLambdaKrigTime
      count = count + 1
    }
    ylocPredList = spatpredFixedL$Y
    ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * Nu * O)) #(testingLocNum x T x O) x NKeep 
    ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x T (rowMeans) note that O = 1 here
    postG <- mean((ylocPredMean - Ytrueloc)^2)
    diffMat <- sweep(ylocPred, 1, Ytrueloc, "-")
    postMSE = mean(rowMeans(diffMat^2))
    postVar = mean(apply(ylocPred, 1, var))
    postGmat[nseed, n1] <- postG
    postMSEmat[nseed, n1] <- postMSE
    postVarMat[nseed, n1] <- postVar
  }
}

save(spatpredTimeVec, file = "fullGPfixedLspatpredTimeVec.RData")
save(alphaKrigTimeVec, file = "fullGPfixedLalphaKrigTimeVec.RData")
save(weightsXiLambdaKrigTimeVec, file = "fullGPfixedLweightsXiLambdaKrigTimeVec.RData")
save(postGmat, file = "fullGPfixedLpostGmat_spatPredMetricSimu.RData")
save(postMSEmat, file = "fullGPfixedLpostMSEmat_spatPredMetricSimu.RData")
save(postVarMat, file = "fullGPfixedLpostVarMat_spatPredMetricSimu.RData")
