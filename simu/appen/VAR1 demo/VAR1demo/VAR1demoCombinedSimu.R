rm(list=ls())
library(mvtnorm)
library(fields)
library(coda)
library(spatempBFA)
K <- 5
O <- 1
L <- 50
M <- 225
LjVec <- rep(min(50, M), K)
sqrootM <- 15
Nu <- 200
Time <- 1:Nu
TimeDist <- as.matrix(dist(Time))
APsi = 0.1; BPsi = 4.5
set.seed(29)
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
maxL <- 10
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
#Tau <- cumprod(Delta)
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
Hypers <- list(Sigma2 = list(A = 1, B = 1), Rho = list(ARho=0.1, BRho=1),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
MCMC <- list(NBurn = 80000, NSims = 20000, NThin = 2, NPilot = 5) 

Sys.time()
VAR1bfaFixedL.simu <- VAR1bfaFixedL(Y ~ 0, data = dat, dist = D, Nu, K, L = L, 
                                    family = "normal", spatial.structure = "continuous",
                                    starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC, seed = 27,
                                    gamma.shrinkage = TRUE, include.space = TRUE, 
                                    spatApprox = FALSE, alphaMethod = "block", h = 15,
                                    storeSpatPredPara = FALSE, storeWeights = FALSE, alphasWeightsToFiles = FALSE) 
Sys.time()
VAR1bfaFixedL.simu$runtime
DiagsVAR1 <- diagnostics(VAR1bfaFixedL.simu, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
DevianceVAR1 <- as.mcmc(DiagsVAR1$deviance)
save(VAR1bfaFixedL.simu, file="regFixedL50simuT200M225Iter100000.RData")
save(DiagsVAR1, file = "regFixedL50simuT200M225Iter100000Diags.RData")
save(DevianceVAR1, file = "regFixedL50simuT200M225Iter100000Deviance.RData")

Sys.time()
VAR1bfaFixedL.simu.block <- VAR1bfaFixedL(Y ~ 0, data = dat, dist = D, Nu, K, L = L, 
                                          family = "normal", spatial.structure = "continuous",
                                          starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC, seed = 27,
                                          gamma.shrinkage = TRUE, include.space = TRUE, 
                                          spatApprox = TRUE, alphaMethod = "block", h = 15,
                                          storeSpatPredPara = FALSE, storeWeights = FALSE, alphasWeightsToFiles = FALSE) 
Sys.time()
VAR1bfaFixedL.simu.block$runtime
Diags.VAR1fixedL.block <- diagnostics(VAR1bfaFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
Deviance.VAR1fixedL.block <- as.mcmc(Diags.VAR1fixedL.block$deviance)
save(VAR1bfaFixedL.simu.block, file="regFixedL50simuBlockT200M225Iter100000.RData")
save(Diags.VAR1fixedL.block, file = "regFixedL50simuBlockT200M225Iter100000Diags.RData")
save(Deviance.VAR1fixedL.block, file = "regFixedL50simuBlockT200M225Iter100000Deviance.RData")

Sys.time()
VAR1bfaFixedL.simu.sequen <- VAR1bfaFixedL(Y ~ 0, data = dat, dist = D, Nu, K, L = L, 
                                           family = "normal", spatial.structure = "continuous",
                                           starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC, seed = 27,
                                           gamma.shrinkage = TRUE, include.space = TRUE, 
                                           spatApprox = TRUE, alphaMethod = "sequential", h = 15,
                                           storeSpatPredPara = FALSE, storeWeights = FALSE, alphasWeightsToFiles = FALSE) 
Sys.time()
VAR1bfaFixedL.simu.sequen$runtime
Diags.VAR1fixedL.sequen <- diagnostics(VAR1bfaFixedL.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
Deviance.VAR1fixedL.sequen <- as.mcmc(Diags.VAR1fixedL.sequen$deviance)
save(VAR1bfaFixedL.simu.sequen, file="regFixedL50simuSequenT200M225Iter100000.RData")
save(Diags.VAR1fixedL.sequen, file = "regFixedL50simuSequenT200M225Iter100000Diags.RData")
save(Deviance.VAR1fixedL.sequen, file = "regFixedL50simuSequenT200M225Iter100000Deviance.RData")

Sys.time()
VAR1regVaryLj.simu.sequen <- VAR1bfaVaryingLjs(Y ~ 0, data = dat, dist = D, Nu, K, LjVec, 
                                               family = "normal", spatial.structure = "continuous",
                                               starting = NULL, hypers = NULL, tuning = NULL, mcmc = MCMC, seed = 27,
                                               gamma.shrinkage = TRUE, include.space = TRUE, 
                                               spatApprox = TRUE, alphaSequen = TRUE, h = 15,
                                               storeSpatPredPara = FALSE, storeWeights = FALSE) 
Sys.time()
VAR1regVaryLj.simu.sequen$runtime
Diags.VAR1sequenVaryLj.sequen <- diagnostics(VAR1regVaryLj.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
Deviance.VAR1sequenVaryLj.sequen <- as.mcmc(Diags.VAR1sequenVaryLj.sequen$deviance)
save(VAR1regVaryLj.simu.sequen, file="regVaryLjSimuSequenT200M225Iter100000.RData")
save(Diags.VAR1sequenVaryLj.sequen, file = "regVaryLjSimuSequenT200M225Iter100000Diags.RData")
save(Deviance.VAR1sequenVaryLj.sequen, file = "regVaryLjSimuSequenT200M225Iter100000Deviance.RData")
