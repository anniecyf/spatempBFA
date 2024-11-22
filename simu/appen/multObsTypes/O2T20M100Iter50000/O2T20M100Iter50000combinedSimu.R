rm(list=ls())
library(mvtnorm)
library(fields)
library(coda)
library(spatempBFA)
K <- 5
O <- 2
L <- 30
M <- 100
LjVec <- rep(min(30, M), K)
sqrootM <- 10
Nu <- 20
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
set.seed(29)
### 1) actual sigma^2(i,o) values
sigma2 <- 0.01
### 2) actual psi value
psi <- 2.3
### 3) actual kappa (O\times O)
kappa <- 0.7 * diag(O)
### 4) actual Upsilon(K\times K)
tempMat <- matrix(runif(K*K,0,1),K,K)
Upsilon <- t(tempMat)%*%tempMat
### 5) actual rho value
rho <- 0.8
D <- rdist(expand.grid(1:sqrootM, 1:sqrootM))
Frho <- exp(-rho*D)
### 6) actual Eta (c(Eta_1,...,Eta_T)) (vec of length Nu*K)
Hpsi <- exp(-psi*TimeDist)
Eta <- rmvnorm(1, mean=rep(0, Nu*K), sigma=kronecker(Hpsi, Upsilon))
## Y~0 (P=0) so no need to sample Beta; all familyInd=0 (all normal) so no need to sample Y
maxL <- 10
LStarJ <- sample(maxL, size=K, replace=T)
### 7) actual alpha
Alpha <- list()
for(j in 1:K) {
  Alpha[[j]] <- t(rmvnorm(LStarJ[j], mean=rep(0,M*O), sigma=kronecker(Frho, kappa))) #every list index an M*O by L_j matrix 
}
w <- list()
for(j in 1:K){
  PhiAlphaJ <- pnorm(Alpha[[j]])
  Lj <- LStarJ[j]
  wj <- matrix(1, M*O, Lj)
  for(o in 1:O){
    for(i in 1:M){
      wj[(o-1)*M+i, 1:(Lj-1)] = PhiAlphaJ[(i-1)*O+o, 1:(Lj-1)]
    }
  } #w[[j]][,Lj] <- rep(1, M*O)
  w[[j]] <- wj 
  temp <- rep(1, M*O)
  for(l in 1:Lj){
    w[[j]][,l] <- w[[j]][,l]*temp
    if(l<Lj) {
      alphaJLreordered <- as.vector(t(matrix(Alpha[[j]][,l], O, M)))
      temp <- temp * pnorm(alphaJLreordered, lower.tail = FALSE)
    }
  }
}
### 8) actual Xi
Xi <- matrix(1, M*O, K)
for(j in 1:K){
  Lj <- LStarJ[j]
  wj <- w[[j]] #M*O by L_j 
  for(row in 1:(M*O)){
    Xi[row,j] <- sample(Lj, size=1, prob=wj[row,])
  }
}
### 9) actual Delta
a1=1; a2=10
Delta <- sapply(c(a1,rep(a2,(K-1))), rgamma, n=1, rate=1)
### 10) actual Theta
Theta <- list()
for(j in 1:K){
  Theta[[j]] <- rnorm(LStarJ[j], 0, sd=sqrt(1/Delta[j])) #vector of length Lj
}
Lambda <- matrix(0, M*O, K)
for(j in 1:K){
  Lambda[,j] = Theta[[j]][Xi[,j]]
}

Sigma.NuMO <- rnorm(Nu * M * O, sd = sqrt(sigma2))
EtaMat <- matrix(Eta, K, Nu)
meanMat <- Lambda%*%EtaMat #M*O\times Nu
Yobs <- as.vector(meanMat) + Sigma.NuMO
dat <- data.frame(Y = Yobs)

MCMC <- list(NBurn = 30000, NSims = 20000, NThin = 2, NPilot = 5)

Sys.time()
regFixedL.simu <- bfaFixedL(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
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
                            storeSpatPredPara = FALSE,
                            storeWeights = FALSE,
                            alphasWeightsToFiles = FALSE) 
Sys.time()
regFixedL.simu$runtime
save(regFixedL.simu, file="regFixedL30simuO2T20M100Iter50000.RData")
# GibbsStepTimeFixedLfullGP <- regFixedL.simu$GibbsStepTime
# save(GibbsStepTimeFixedLfullGP, file = "GibbsStepTimeFixedLfullGP.RData")
Diags <- diagnostics(regFixedL.simu, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags, file = "regFixedL30simuO2T20M100Iter50000Diags.RData")
Deviance <- as.mcmc(Diags$deviance)
save(Deviance, file = "regFixedL30simuO2T20M100Iter50000Deviance.RData")

Sys.time()
regFixedL.simu.block <- bfaFixedL(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
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
                                  storeSpatPredPara = FALSE,
                                  storeWeights = FALSE,
                                  alphasWeightsToFiles = FALSE) 
Sys.time()
regFixedL.simu.block$runtime
save(regFixedL.simu.block, file="regFixedL30simuBlockO2T20M100Iter50000.RData")
# GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
# save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
Diags.block <- diagnostics(regFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.block, file = "regFixedL30simuBlockO2T20M100Iter50000Diags.RData")
Deviance.block <- as.mcmc(Diags.block$deviance)
save(Deviance.block, file = "regFixedL30simuBlockO2T20M100Iter50000Deviance.RData")

Sys.time()
regFixedL.simu.sequen <- bfaFixedL(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
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
                                   storeSpatPredPara = FALSE,
                                   storeWeights = FALSE,
                                   alphasWeightsToFiles = FALSE) 
Sys.time()
save(regFixedL.simu.sequen, file="regFixedL30simuSequenO2T20M100Iter50000.RData")
regFixedL.simu.sequen$runtime
# GibbsStepTimeFixedLsequen <- regFixedL.simu.sequen$GibbsStepTime
# save(GibbsStepTimeFixedLsequen, file = "GibbsStepTimeFixedLsequen.RData")
Diags.sequen <- diagnostics(regFixedL.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequen, file = "regFixedL30simuSequenO2T20M100Iter50000Diags.RData")
Deviance.sequen <- as.mcmc(Diags.sequen$deviance)
save(Deviance.sequen, file = "regFixedL30simuSequenO2T20M100Iter50000Deviance.RData")

Sys.time()
regVaryLj.simu.sequen <- bfaVaryingLjs(Y ~ 0, data = dat, dist = D, time = Time,  K = K, LjVec = LjVec,
                                       starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC,
                                       family = "normal",
                                       temporal.structure = "exponential",
                                       spatial.structure = "continuous",
                                       seed = 19, 
                                       gamma.shrinkage = TRUE,
                                       include.time = TRUE,
                                       include.space = TRUE,
                                       seasonPeriod = 1, 
                                       equalTimeDist = TRUE,
                                       spatApprox = TRUE, 
                                       alphaSequen = TRUE, 
                                       h = 15,
                                       storeSpatPredPara = FALSE, 
                                       storeWeights = FALSE) 
Sys.time()
save(regVaryLj.simu.sequen, file="regVaryLjSimuSequenO2T20M100Iter50000.RData")
regVaryLj.simu.sequen$runtime
# GibbsStepTimeVaryLjSequen <- regVaryLj.simu.sequen$GibbsStepTime
# save(GibbsStepTimeVaryLjSequen, file = "GibbsStepTimeVaryLjSequen.RData")
Diags.sequenVaryLj <- diagnostics(regVaryLj.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequenVaryLj, file = "regVaryLjSimuSequenO2T20M100Iter50000Diags.RData")
Deviance.sequenVaryLj <- as.mcmc(Diags.sequenVaryLj$deviance)
save(Deviance.sequenVaryLj, file = "regVaryLjSimuSequenO2T20M100Iter50000Deviance.RData")

