rm(list=ls())
library(mvtnorm)
library(fields)
library(spatempBFA)
K <- 5
O <- 2
L <- 50
M <- 400
LjVec <- rep(min(50, M), K)
sqrootM <- 20
Nu <- 200
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

MCMC <- list(NBurn = 80000, NSims = 20000, NThin = 2, NPilot = 5)

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
                                  storeSpatPredPara = TRUE,
                                  storeWeights = TRUE,
                                  alphasWeightsToFiles = TRUE) 
Sys.time()
regFixedL.simu.block$runtime
save(regFixedL.simu.block, file="regFixedL50simuBlockO2T200M400Iter100000.RData")
GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
library(coda)
Diags.block <- diagnostics(regFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.block, file = "regFixedL50simuBlockO2T200M400Iter100000Diags.RData")
Deviance.block <- as.mcmc(Diags.block$deviance)
save(Deviance.block, file = "regFixedL50simuBlockO2T200M400Iter100000Deviance.RData")
