rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/mainScalabilityVerificationSimu/m900T50K5")
list.files()
load("regFixedL50simuT50M900Iter30000.RData")
load("regFixedL50simuBlockT50M900Iter30000.RData")
load("regFixedL50simuSequenT50M900Iter30000nostorealphaweights.RData")
load("regVaryLjsimuSequenT50M900Iter30000nostorealphaweight.RData")
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











