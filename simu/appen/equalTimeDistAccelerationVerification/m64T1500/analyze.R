rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/appen/equalTimeDistAccelerationVerification/m64T1500")
list.files()
load("regFixedL30simuT1500M64Iter30000.RData")
load("regFixedL30simuBlockT1500M64Iter30000.RData")
load("regFixedL30simuSequenT1500M64Iter30000nostorealphaweights.RData")
load("regVaryLjsimuSequenT1500M64Iter30000nostorealphaweight.RData")
GibbsStepTimeFixedLfullGP <- regFixedL.simu$GibbsStepTime
GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
GibbsStepTimeFixedLsequen <- regFixedL.simu.sequen$GibbsStepTime
GibbsStepTimeVaryLjSequen <- regVaryLj.simu.sequen$GibbsStepTime
save(GibbsStepTimeFixedLfullGP, file = "GibbsStepTimeFixedLfullGP.RData")
save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
save(GibbsStepTimeFixedLsequen, file = "GibbsStepTimeFixedLsequen.RData")
save(GibbsStepTimeVaryLjSequen, file = "GibbsStepTimeVaryLjSequen.RData")

load("GibbsStepTimeFixedLfullGP.RData")
load("GibbsStepTimeFixedLblock.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("GibbsStepTimeVaryLjSequen.RData")













