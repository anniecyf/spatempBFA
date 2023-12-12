rm(list=ls())
library(mvtnorm)
library(fields)
library(spatTempBFA)
K <- 6
O <- 1
L <- 60
M <- 900
LjVec <- rep(min(100, M), K)
sqrootM <- 30
Nu <- 30
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
# set random seed
### 1) actual sigma^2(i,o) (for i=1,2,...,M and o=1) values
sigma2 <- 0.01
### 2) actual psi value
psi <- 2.3
### 3) actual kappa value
kappa <- 0.7
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
maxL <- 20
LStarJ <- sample(maxL, size=K, replace=T)
### 7) actual alpha
Alpha <- list()
for(j in 1:K) {
  Alpha[[j]] <-  t(rmvnorm(LStarJ[j], mean=rep(0,M*O), sigma=kappa*Frho))
  #every list index an M by L_j matrix 
}
w <- list()
for(j in 1:K){
  w[[j]] <- pnorm(Alpha[[j]])
  Lj <- LStarJ[j]
  w[[j]][,Lj] <- rep(1, M)
  temp <- rep(1, M)
  for(l in 1:Lj){
    w[[j]][,l] <- w[[j]][,l]*temp
    if(l<Lj) {temp <- temp * pnorm(Alpha[[j]][,l], lower.tail = FALSE)}
  }
}
### 8) actual Xi
Xi <- matrix(1, M, K)
for(j in 1:K){
  Lj <- LStarJ[j]
  for(i in 1:M){
    Xi[i,j] <- sample(Lj, size=1, prob=w[[j]][i,])
  }
}
### 9) actual Delta
a1=1; a2=10
Delta <- sapply(c(a1,rep(a2,(K-1))), rgamma, n=1, rate=1)
Tau <- cumprod(Delta)
### 10) actual Theta
Theta <- list()
for(j in 1:K){
  Theta[[j]] <- rnorm(LStarJ[j], 0, sd=sqrt(1/Tau[j])) #vector of length Lj
}
Lambda <- matrix(0, M*O, K)
for(j in 1:K){
  for(i in 1:M){
    Lambda[i,j] = Theta[[j]][Xi[i,j]]
  }
}

Sigma.NuMO <- rnorm(Nu * M * O, sd = sqrt(sigma2))
EtaMat <- matrix(Eta, K, Nu)
meanMat <- Lambda%*%EtaMat #M*O\times Nu
Yobs <- as.vector(meanMat) + Sigma.NuMO
dat <- data.frame(Y = Yobs)

paraData <- list(Sigma=Sigma, psi=psi, kappa=kappa, Upsilon=Upsilon,
                 Lambda=Lambda, Eta=Eta, Yobs=Yobs, K=K, O=O, M=M, Nu=Nu, 
                 L=L, LjVec=LjVec,
                 Hpsi=Hpsi, APsi=APsi, BPsi=BPsi,
                 TimeDist = TimeDist,
                 rho=rho, D=D, Frho=Frho, maxL=maxL, LStarJ=LStarJ,
                 Alpha=Alpha, w=w, Xi=Xi, a1=a1, a2=a2, Delta=Delta,
                 Tau=Tau, Theta=Theta)
save(paraData, file="simuT30M900paraData.RData")

Hypers <- list(Sigma2 = list(A = 1, B = 1), Rho = list(ARho=0.1, BRho=1),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
MCMC <- list(NBurn = 20000, NSims = 10000, NThin = 2, NPilot = 5)

Sys.time()
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
                            storeSpatPredPara = FALSE,
                            storeWeights = FALSE,
                            alphasWeightsToFiles = FALSE) 
Sys.time()
regFixedL.simu$runtime
save(regFixedL.simu, file="regFixedL60simuT30M900Iter30000.RData")
library(coda)
Diags <- diagnostics(regFixedL.simu, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags, file = "regFixedL60simuT30M900Iter30000Diags.RData")
Deviance <- as.mcmc(Diags$deviance)
save(Deviance, file = "regFixedL60simuT30M900Iter30000Deviance.RData")

Sys.time()
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
                                  storeSpatPredPara = FALSE,
                                  storeWeights = FALSE,
                                  alphasWeightsToFiles = FALSE) 
Sys.time()
regFixedL.simu.block$runtime
save(regFixedL.simu.block, file="regFixedL60simuBlockT30M900Iter30000.RData")
library(coda)
Diags.block <- diagnostics(regFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.block, file = "regFixedL60simuBlockT30M900Iter30000Diags.RData")
Deviance.block <- as.mcmc(Diags.block$deviance)
save(Deviance.block, file = "regFixedL60simuBlockT30M900Iter30000Deviance.RData")

Sys.time()
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
                                   storeSpatPredPara = FALSE,
                                   storeWeights = FALSE,
                                   alphasWeightsToFiles = FALSE) 
Sys.time()
save(regFixedL.simu.sequen, file="regFixedL60simuSequenT30M900Iter30000nostorealphaweights.RData")
regFixedL.simu.sequen$runtime
library(coda)
Diags.sequen <- diagnostics(regFixedL.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequen, file = "regFixedL60simuSequenT30M900Iter30000DiagsNostorealphaweights.RData")
Deviance.sequen <- as.mcmc(Diags.sequen$deviance)
save(Deviance.sequen, file = "regFixedL60simuSequenT30M900Iter30000DevianceNostorealphaweights.RData")

Sys.time()
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
                                       h = 15,
                                       storeSpatPredPara = FALSE, 
                                       storeWeights = FALSE) 
Sys.time()
save(regVaryLj.simu.sequen, file="regVaryLjsimuSequenT30M900Iter30000nostorealphaweight.RData")
regVaryLj.simu.sequen$runtime
library(coda)
Diags.sequenVaryLj <- diagnostics(regVaryLj.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequenVaryLj, file = "regVaryLjsimuSequenT30M900Iter30000DiagsNostorealphaweight.RData")
Deviance.sequenVaryLj <- as.mcmc(Diags.sequenVaryLj$deviance)
save(Deviance.sequenVaryLj, file = "regVaryLjsimuSequenT30M900Iter30000DevianceNostorealphaweight.RData")
















