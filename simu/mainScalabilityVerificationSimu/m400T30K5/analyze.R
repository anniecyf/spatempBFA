rm(list=ls())
library(tidyverse)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/mainScalabilityVerificationSimu/m400T30K5")
list.files()
load("regFixedL50simuT30M400Iter30000.RData")
load("regFixedL50simuBlockT30M400Iter30000.RData")
load("regFixedL50simuSequenT30M400Iter30000nostorealphaweights.RData")
load("regVaryLjsimuSequenT30M400Iter30000nostorealphaweight.RData")
GibbsStepTimeFixedLfullGP <- regFixedL.simu$GibbsStepTime
GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
GibbsStepTimeFixedLsequen <- regFixedL.simu.sequen$GibbsStepTime
GibbsStepTimeVaryLjSequen <- regVaryLj.simu.sequen$GibbsStepTime
save(GibbsStepTimeFixedLfullGP, file = "GibbsStepTimeFixedLfullGP.RData")
save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
save(GibbsStepTimeFixedLsequen, file = "GibbsStepTimeFixedLsequen.RData")
save(GibbsStepTimeVaryLjSequen, file = "GibbsStepTimeVaryLjSequen.RData")













