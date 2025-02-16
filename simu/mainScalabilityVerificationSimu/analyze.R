rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/mainScalabilityVerificationSimu/m3600T100K5")
list.files()
load("regFixedL50simuT100M3600Iter30000.RData")
load("regFixedL50simuBlockT100M3600Iter30000.RData")
load("regFixedL50simuSequenT100M3600Iter30000nostorealphaweights.RData")
load("regVaryLjsimuSequenT100M3600Iter30000nostorealphaweight.RData")
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











