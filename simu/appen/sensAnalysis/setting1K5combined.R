rm(list=ls())
library(mvtnorm)
library(fields)
library(spatempBFA)
library(coda)
K <- 5
O <- 1
L <- 50
M <- 225
LjVec <- rep(min(50, M), K)
sqrootM <- 15
Nu <- 100
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
### setting 1
Hypers <- list(Sigma2 = list(A = 1, B = 1), 
               Delta = list(A1 = 1, A2 = 1),
               Rho = list(ARho = 0.1, BRho = 1),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
# ### setting 2
# Hypers <- list(Sigma2 = list(A = 1, B = 1), 
#                Delta = list(A1 = 1, A2 = 1),
#                Rho = list(ARho = 1, BRho = 2),
#                Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
#                Psi = list(APsi = 1, BPsi = 2),
#                Upsilon = list(Zeta = K + 1, Omega = diag(K)))
# ### setting 3
# Hypers <- list(Sigma2 = list(A = 1, B = 1), 
#                Delta = list(A1 = 1, A2 = 1),
#                Rho = list(ARho = 0.1, BRho = 1),
#                Kappa = list(SmallUpsilon = O + 1, BigTheta = 0.1 * diag(O)),
#                Psi = list(APsi = APsi, BPsi = BPsi),
#                Upsilon = list(Zeta = K + 1, Omega = 0.1 * diag(K)))
# ### setting 4
# Hypers <- list(Sigma2 = list(A = 0.1, B = 1), 
#                Delta = list(A1 = 1, A2 = 1),
#                Rho = list(ARho = 0.1, BRho = 1),
#                Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
#                Psi = list(APsi = APsi, BPsi = BPsi),
#                Upsilon = list(Zeta = K + 1, Omega = diag(K)))
# ### setting 5
# Hypers <- list(Sigma2 = list(A = 1, B = 1), 
#                Delta = list(A1 = 1, A2 = 2),
#                Rho = list(ARho = 0.1, BRho = 1),
#                Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
#                Psi = list(APsi = APsi, BPsi = BPsi),
#                Upsilon = list(Zeta = K + 1, Omega = diag(K)))


MCMC <- list(NBurn = 90000, NSims = 10000, NThin = 2, NPilot = 5)

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
save(regFixedL.simu, file="regFixedL50simuT100M225K5Iter100000.RData")
# GibbsStepTimeFixedLfullGP <- regFixedL.simu$GibbsStepTime
# save(GibbsStepTimeFixedLfullGP, file = "GibbsStepTimeFixedLfullGP.RData")
Diags <- diagnostics(regFixedL.simu, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags, file = "regFixedL50simuT100M225K5Iter100000Diags.RData")
Deviance <- as.mcmc(Diags$deviance)
save(Deviance, file = "regFixedL50simuT100M225K5Iter100000Deviance.RData")

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
save(regFixedL.simu.block, file="regFixedL50simuBlockT100M225K5Iter100000.RData")
# GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
# save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
Diags.block <- diagnostics(regFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.block, file = "regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
Deviance.block <- as.mcmc(Diags.block$deviance)
save(Deviance.block, file = "regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")

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
save(regFixedL.simu.sequen, file="regFixedL50simuSequenT100M225K5Iter100000.RData")
regFixedL.simu.sequen$runtime
# GibbsStepTimeFixedLsequen <- regFixedL.simu.sequen$GibbsStepTime
# save(GibbsStepTimeFixedLsequen, file = "GibbsStepTimeFixedLsequen.RData")
Diags.sequen <- diagnostics(regFixedL.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequen, file = "regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
Deviance.sequen <- as.mcmc(Diags.sequen$deviance)
save(Deviance.sequen, file = "regFixedL50simuSequenT100M225K5Iter100000Deviance.RData")

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
save(regVaryLj.simu.sequen, file="regVaryLjSimuSequenT100M225K5Iter100000.RData")
regVaryLj.simu.sequen$runtime
# GibbsStepTimeVaryLjSequen <- regVaryLj.simu.sequen$GibbsStepTime
# save(GibbsStepTimeVaryLjSequen, file = "GibbsStepTimeVaryLjSequen.RData")
Diags.sequenVaryLj <- diagnostics(regVaryLj.simu.sequen, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.sequenVaryLj, file = "regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
Deviance.sequenVaryLj <- as.mcmc(Diags.sequenVaryLj$deviance)
save(Deviance.sequenVaryLj, file = "regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData")




rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/sensAnalysis/diffHyperPara"))
setwd("setting1")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") # Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData") # Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") # Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") # Deviance.sequenVaryLj
NKeep <- 5000
postDeviancesDF <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance), NNGPblockFixedL = as.vector(Deviance.block),
                              NNGPsequenFixedL = as.vector(Deviance.sequen), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj))
fullGPfixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances
NNGPblockFixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances
NNGPsequenFixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances
NNGPsequenVaryLjpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances
ggarrange(fullGPfixedLpostDeviances, NNGPblockFixedLpostDeviances, 
          NNGPsequenFixedLpostDeviances, NNGPsequenVaryLjpostDeviances,
          labels = c("A", "B", "C", "D"), align = "hv",
          ncol = 2, nrow = 2)
ggsave("diffHyPerParaSetting1postDeviances.png", width = 16, height = 16, units = "cm")




