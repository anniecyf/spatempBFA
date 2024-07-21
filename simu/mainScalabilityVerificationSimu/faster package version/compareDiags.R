rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/mainScalabilityVerificationSimu/faster package version")
load("regFixedL50simuBlockT30M400Iter30000Diags.RData")
load("regFixedL50simuSequenT30M400Iter30000DiagsNostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000DiagsNostorealphaweight.RData")
load("regFixedL50simuT30M400Iter30000Diags.RData")
Diags.fast <- Diags
Diags.block.fast <- Diags.block
Diags.sequen.fast <- Diags.sequen
Diags.sequenVaryLj.fast <- Diags.sequenVaryLj
load("regFixedL50simuBlockT30M400Iter30000.RData")
load("regFixedL50simuSequenT30M400Iter30000Nostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000Nostorealphaweight.RData")
load("regFixedL50simuT30M400Iter30000.RData")
regFixedL.simu.fast <- regFixedL.simu
regFixedL.simu.block.fast <- regFixedL.simu.block
regFixedL.simu.sequen.fast <- regFixedL.simu.sequen
regVaryLj.simu.sequen.fast <- regVaryLj.simu.sequen
# much slower when much more scripts running tgt so decided not to run the ones with larger m
setwd("../m400T30K5")
load("regFixedL50simuBlockT30M400Iter30000Diags.RData")
load("regFixedL50simuSequenT30M400Iter30000DiagsNostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000DiagsNostorealphaweight.RData")
load("regFixedL50simuT30M400Iter30000Diags.RData")
load("regFixedL50simuBlockT30M400Iter30000.RData")
load("regFixedL50simuSequenT30M400Iter30000Nostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000Nostorealphaweight.RData")
load("regFixedL50simuT30M400Iter30000.RData")
# all Diags agree; reg objects also all agree except runtime and that the slower ones 1 more list component GibbsStepTime
setwd("../spBFA")
load("spBFAL50simuT30M400Iter30000Diags.RData")
load("spBFALInfsimuT30M400Iter30000Diags.RData")

rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/mainScalabilityVerificationSimu/faster package version")
load("regFixedL50simuBlockT50M400Iter30000Diags.RData")
load("regFixedL50simuSequenT50M400Iter30000DiagsNostorealphaweights.RData")
load("regVaryLjsimuSequenT50M400Iter30000DiagsNostorealphaweight.RData")
load("regFixedL50simuT50M400Iter30000Diags.RData")
Diags.fast <- Diags
Diags.block.fast <- Diags.block
Diags.sequen.fast <- Diags.sequen
Diags.sequenVaryLj.fast <- Diags.sequenVaryLj
load("regFixedL50simuBlockT50M400Iter30000.RData")
load("regFixedL50simuSequenT50M400Iter30000Nostorealphaweights.RData")
load("regVaryLjsimuSequenT50M400Iter30000Nostorealphaweight.RData")
load("regFixedL50simuT50M400Iter30000.RData")
regFixedL.simu.fast <- regFixedL.simu
regFixedL.simu.block.fast <- regFixedL.simu.block
regFixedL.simu.sequen.fast <- regFixedL.simu.sequen
regVaryLj.simu.sequen.fast <- regVaryLj.simu.sequen
setwd("../m400T50K5")
load("regFixedL50simuBlockT50M400Iter30000Diags.RData")
load("regFixedL50simuSequenT50M400Iter30000DiagsNostorealphaweights.RData")
load("regVaryLjsimuSequenT50M400Iter30000DiagsNostorealphaweight.RData")
load("regFixedL50simuT50M400Iter30000Diags.RData")
load("regFixedL50simuBlockT50M400Iter30000.RData")
load("regFixedL50simuSequenT50M400Iter30000Nostorealphaweights.RData")
load("regVaryLjsimuSequenT50M400Iter30000Nostorealphaweight.RData")
load("regFixedL50simuT50M400Iter30000.RData")
# all Diags agree; reg objects also all agree except runtime and whether there's the list component GibbsStepTime
# setwd("../spBFA")
# load("spBFAL50simuT50M400Iter30000Diags.RData")
# load("spBFALInfsimuT50M400Iter30000Diags.RData")
