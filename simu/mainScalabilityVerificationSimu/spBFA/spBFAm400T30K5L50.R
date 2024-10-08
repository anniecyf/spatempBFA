rm(list=ls())
library(mvtnorm)
library(fields)
library(spBFA)
K <- 5
O <- 1
L <- 50
M <- 400
sqrootM <- 20
Nu <- 30
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
  for(i in 1:M){
    Lambda[i,j] = Theta[[j]][Xi[i,j]]
  }
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
MCMC <- list(NBurn = 20000, NSims = 10000, NThin = 2, NPilot = 5)


Sys.time()
reg.simu <- bfa_sp(Y ~ 0, data = dat, dist = D, time = Time,  K = K, 
                   starting = NULL, hypers = Hypers, tuning = NULL, mcmc = MCMC,
                   L = L,
                   family = "normal",
                   trials = NULL,
                   temporal.structure = "exponential",
                   spatial.structure = "continuous",
                   seed = 29, 
                   gamma.shrinkage = TRUE,
                   include.space = TRUE,
                   clustering = TRUE) 
save(reg.simu, file="spBFAL50simuT30M400Iter30000.RData")
Sys.time()
reg.simu$runtime
library(coda)
Diags <- spBFA::diagnostics(reg.simu, diags = c("dic", "dinf", "waic"), keepDeviance = TRUE)
save(Diags, file="spBFAL50simuT30M400Iter30000Diags.RData")
Deviance <- as.mcmc(Diags$deviance)
#save(Deviance, file = "spBFAL50simuT30M400Iter30000Deviance.RData")

