rm(list=ls())
library(fields)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
#setwd(paste(projDirec, "simu/appen/real data ex/SSTm1612", sep = "/"))
load("SSTtrainingDF2023.RData")
K <- 5
O <- 1
L <- 50
M <- nrow(SSTtrainingDF)
D <- rdist(SSTtrainingDF[,1:2])
LjVec <- rep(min(50, M), K)
Nu <- 350
Time <- 1:Nu
SSTtrainingDF = SSTtrainingDF[, 3:(Nu+2)]
set.seed(29)
testingLocs <- sample(1:nrow(SSTtrainingDF), size = 35, replace = FALSE)
APsi = 0.1; BPsi = 4.5
dat <- data.frame(Y = as.vector(as.matrix(SSTtrainingDF[-testingLocs,])))
#YtestingSpat <- as.vector(as.matrix(SSTtrainingDF[testingLocs,]))
Dtraining <- as.matrix(D[-testingLocs, -testingLocs])
distOrigNew = as.matrix(D[-testingLocs, testingLocs])
distNewNew = as.matrix(D[testingLocs, testingLocs])
rm(SSTtrainingDF)
Hypers <- list(Sigma2 = list(A = 1, B = 1), Rho = list(ARho = 0.1, BRho = 1),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
MCMC <- list(NBurn = 100000, NSims = 20000, NThin = 20, NPilot = 5)
library(spatempBFA)
library(coda)
Sys.time()
regFixedL.simu.block <- bfaFixedL(Y ~ 0, data = dat, dist = Dtraining, time = Time,  K = K, 
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
                                  storeSpatPredPara = TRUE,
                                  storeWeights = TRUE,
                                  alphasWeightsToFiles = FALSE) 
Sys.time()
save(regFixedL.simu.block, file = "realExNNGPblockFixedL.RData")
regFixedL.simu.block$runtime
GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
Diags.block <- diagnostics(regFixedL.simu.block, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = TRUE)
save(Diags.block, file = "NNGPblockFixedLDiags.RData")
Deviance.block <- as.mcmc(Diags.block$deviance)
save(Deviance.block, file = "NNGPblockFixedLDeviance.RData")
spatpredFixedLblock <- predictNewLocFixedL(regFixedL.simu.block, 35, distOrigNew, distNewNew, 
                                            NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 29)
save(spatpredFixedLblock, file = "realExNNGPspatpredFixedLblock.RData")
temppredFixedLblock <- predictNewTime(regFixedL.simu.block, (Nu+1):365, seed = 29)
save(temppredFixedLblock, file = "realExNNGPtemppredFixedLblock.RData")