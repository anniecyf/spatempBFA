rm(list=ls())
N <- 2 
m <- 15
M <- m + 1
h <- 7
K <- 3
O <- 1
L <- min(20, M)
LjVec <- rep(min(20, M), K)
sqrootM <- 4
Nu <- 15
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
library(mvtnorm)
library(fields)
library(devtools)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/server/spatempBFA")
load_all(".")
library(spatempBFA)
set.seed(19)
fittingTimeMat <- array(0, dim = c(4, N, M))
spatpredTimeMat <- array(0, dim = c(4, N, M))
alphaKrigTimeMat <- array(0, dim = c(4, N, M))
weightsXiLambdaKrigTimeMat <- array(0, dim = c(4, N, M))

for(n in 1:N){
  print(paste("n = ", n, sep = ""))
  sigma2 <- 0.01
  psi <- 2.3
  kappa <- 0.7
  tempMat <- matrix(runif(K*K,0,1), K, K)
  Upsilon <- t(tempMat)%*%tempMat
  rho <- 0.8
  D <- rdist(expand.grid(1:sqrootM, 1:sqrootM))
  Frho <- exp(-rho*D)
  Hpsi <- exp(-psi*TimeDist)
  Eta <- rmvnorm(1, mean=rep(0, Nu*K), sigma=kronecker(Hpsi, Upsilon)) ### actual Eta (c(Eta_1,...,Eta_T)) (vec of length Nu*K)
  maxL <- 10
  LStarJ <- sample(maxL, size=K, replace=T)
  Alpha <- list()
  for(j in 1:K) {
    Alpha[[j]] <-  t(rmvnorm(LStarJ[j], mean=rep(0,M*O), sigma=kappa*Frho))
    #every list index an M by L_j matrix 
  }
  w <- list()
  for(j in 1:K){
    w[[j]] <- pnorm(Alpha[[j]])
    Lj <- LStarJ[j]
    temp <- rep(1, M)
    for(l in 1:(Lj-1)){
      w[[j]][,l] <- w[[j]][,l]*temp
      temp <- temp * pnorm(Alpha[[j]][,l], lower.tail = FALSE)
    }
    w[[j]][,Lj] <- temp
  }
  Xi <- matrix(1, M, K)
  for(j in 1:K){
    Lj <- LStarJ[j]
    for(i in 1:M){
      Xi[i,j] <- sample(Lj, size=1, prob=w[[j]][i,])
    }
  }
  a1=1; a2=10
  Delta <- sapply(c(a1,rep(a2,(K-1))), rgamma, n=1, rate=1)
  #Tau <- cumprod(Delta)
  Theta <- list()
  for(j in 1:K){
    Theta[[j]] <- rnorm(LStarJ[j], 0, sd=sqrt(1/Delta[j])) #vector of length Lj
  }
  Lambda <- matrix(0, M*O, K)
  for(j in 1:K){
    for(i in 1:M){
      Lambda[i,j] = Theta[[j]][Xi[i,j]]
    }
  }
  Hypers <- list(Sigma2 = list(A = 1, B = 1), Rho = list(ARho=0.1, BRho=1),
                 Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
                 Psi = list(APsi = APsi, BPsi = BPsi),
                 Upsilon = list(Zeta = K + 1, Omega = 0.1 * diag(K)))
  MCMC <- list(NBurn = 100, NSims = 100, NThin = 2, NPilot = 5)
  Sigma.NuMO <- rnorm(Nu * M * O, sd = sqrt(sigma2))
  EtaMat <- matrix(Eta, K, Nu)
  meanMat <- Lambda%*%EtaMat #M*O\times Nu
  Yobs <- as.vector(meanMat) + Sigma.NuMO
  
  for(loc in 1 : M){
    print(paste("loc = ", loc, sep = ""))
    Dloc <- D[-loc, -loc]
    discardInd <- seq(from = loc, by = M, length.out = Nu)
    Yobsloc <- Yobs[-discardInd]
    #Ytrueloc <- Yobs[discardInd] # of length T
    datloc <- data.frame(Y = Yobsloc)
    distOrigNew = matrix(D[-loc, loc], m, 1)
    distNewNew = matrix(D[loc, loc], 1, 1)
    
    t1 <- Sys.time()
    regFixedL.simu <- bfaFixedL(Y ~ 0, data = datloc, dist = Dloc, time = Time,  K = K, 
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
    spatpredFixedL <- predictNewLocFixedL(regFixedL.simu, 1, distOrigNew, distNewNew, 
                                          NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
    t3 <- Sys.time()
    fittingTimeMat[1, n, loc] <- difftime(t2, t1, units = "mins")
    spatpredTimeMat[1, n, loc] <- 1000 * difftime(t3, t2, units = "secs") # in milliseconds
    alphaKrigTimeMat[1, n, loc] <- spatpredFixedL$alphaKrigTime
    weightsXiLambdaKrigTimeMat[1, n, loc] <- spatpredFixedL$weightsXiLambdaKrigTime
    
    t1 <- Sys.time()
    regFixedL.simu.block <- bfaFixedL(Y ~ 0, data = datloc, dist = Dloc, time = Time,  K = K, 
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
    spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 1, distOrigNew, distNewNew, 
                                               NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
    t3 <- Sys.time()
    fittingTimeMat[2, n, loc] <- difftime(t2, t1, units = "mins")
    spatpredTimeMat[2, n, loc] <- 1000 * difftime(t3, t2, units = "secs") # in milliseconds
    alphaKrigTimeMat[2, n, loc] <- spatpredFixedLblock$alphaKrigTime
    weightsXiLambdaKrigTimeMat[2, n, loc] <- spatpredFixedLblock$weightsXiLambdaKrigTime
    
    t1 <- Sys.time()
    regFixedL.simu.sequen <- bfaFixedL(Y ~ 0, data = datloc, dist = Dloc, time = Time,  K = K, 
                                       starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
                                       L = L,
                                       family = "normal",
                                       temporal.structure = "exponential",
                                       spatial.structure = "continuous",
                                       seed = 27, 
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
    spatpredFixedLsequen <- predictNewLocFixedL(regFixedL.simu.sequen, 1, distOrigNew, distNewNew, 
                                                NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
    t3 <- Sys.time()
    fittingTimeMat[3, n, loc] <- difftime(t2, t1, units = "mins")
    spatpredTimeMat[3, n, loc] <- 1000 * difftime(t3, t2, units = "secs") # in milliseconds
    alphaKrigTimeMat[3, n, loc] <- spatpredFixedLsequen$alphaKrigTime
    weightsXiLambdaKrigTimeMat[3, n, loc] <- spatpredFixedLsequen$weightsXiLambdaKrigTime
    
    t1 <- Sys.time()
    regVaryLj.simu.sequen <- bfaVaryingLjs(Y ~ 0, data = datloc, dist = Dloc, time = Time,  K = K,
                                           starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
                                           LjVec = LjVec,
                                           family = "normal",
                                           temporal.structure = "exponential",
                                           spatial.structure = "continuous",
                                           seed = 27, 
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
    spatpredVaryLj <- predictNewLocVaryLj(regVaryLj.simu.sequen, 1, distOrigNew, distNewNew, 
                                          NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27)
    t3 <- Sys.time()
    fittingTimeMat[4, n, loc] <- difftime(t2, t1, units = "mins")
    spatpredTimeMat[4, n, loc] <- difftime(t3, t2, units = "secs") # in milliseconds
    alphaKrigTimeMat[4, n, loc] <- spatpredVaryLj$alphaKrigTime
    weightsXiLambdaKrigTimeMat[4, n, loc] <- spatpredVaryLj$weightsXiLambdaKrigTime
    print(Sys.time())
  }
  print(Sys.time())
}

save(fittingTimeMat, file = "M16h7T15fittingTimeMat_randomNewLocPred.RData")
save(spatpredTimeMat, file = "M16h7T15spatpredTimeMat_randomNewLocPred.RData")
save(alphaKrigTimeMat, file = "M16h7T15alphaKrigTimeMat_randomNewLocPred.RData")
save(weightsXiLambdaKrigTimeMat, file = "M16h7T15weightsXiLambdaKrigTimeMat_randomNewLocPred.RData")




























