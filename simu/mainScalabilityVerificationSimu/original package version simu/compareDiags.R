rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/mainScalabilityVerificationSimu/original package version simu")
load("regFixedL50simuBlockT30M400Iter30000Diags.RData")
load("regFixedL50simuSequenT30M400Iter30000DiagsNostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000DiagsNostorealphaweight.RData")
load("regFixedL50simuT30M400Iter30000Diags.RData")
Diags.fast <- Diags
Diags.block.fast <- Diags.block
Diags.sequen.fast <- Diags.sequen
Diags.sequenVaryLj.fast <- Diags.sequenVaryLj
setwd("../m400T30K5")
load("regFixedL50simuBlockT30M400Iter30000Diags.RData")
load("regFixedL50simuSequenT30M400Iter30000DiagsNostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000DiagsNostorealphaweight.RData")
load("regFixedL50simuT30M400Iter30000Diags.RData")
# all Diags agree; reg objects also all agree except runtime and that the slower ones 1 more list component GibbsStepTime
setwd("../spBFA")
load("spBFAL50simuT30M400Iter30000Diags.RData")
load("spBFALInfsimuT30M400Iter30000Diags.RData")





