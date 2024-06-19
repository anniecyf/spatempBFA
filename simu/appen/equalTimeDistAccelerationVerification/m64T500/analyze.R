rm(list=ls())
library(tidyverse)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/appen/equalTimeDistAccelerationVerification/m64T500")
list.files()
load("regFixedL30simuT500M64Iter30000.RData")
load("regFixedL30simuBlockT500M64Iter30000.RData")
load("regFixedL30simuSequenT500M64Iter30000nostorealphaweights.RData")
regFixedL.simu.sequen.fast <- regFixedL.simu.sequen
load("regFixedL30simuSequenT500M64Iter30000specifyEqualTimeDistFnostorealphaweights.RData")
load("regVaryLjsimuSequenT500M64Iter30000nostorealphaweight.RData")
regVaryLj.simu.sequen.fast <- regVaryLj.simu.sequen
load("regVaryLjsimuSequenT500M64Iter30000specifyEqualTimeDistFnostorealphaweight.RData")
GibbsStepTimeFixedLfullGP <- regFixedL.simu$GibbsStepTime
GibbsStepTimeFixedLblock <- regFixedL.simu.block$GibbsStepTime
GibbsStepTimeFixedLsequen.fast <- regFixedL.simu.sequen.fast$GibbsStepTime
GibbsStepTimeVaryLjSequen.fast <- regVaryLj.simu.sequen.fast$GibbsStepTime
GibbsStepTimeFixedLsequen <- regFixedL.simu.sequen$GibbsStepTime
GibbsStepTimeVaryLjSequen <- regVaryLj.simu.sequen$GibbsStepTime
save(GibbsStepTimeFixedLfullGP, file = "GibbsStepTimeFixedLfullGP.RData")
save(GibbsStepTimeFixedLblock, file = "GibbsStepTimeFixedLblock.RData")
save(GibbsStepTimeFixedLsequen, file = "GibbsStepTimeFixedLsequen.RData")
save(GibbsStepTimeVaryLjSequen, file = "GibbsStepTimeVaryLjSequen.RData")
save(GibbsStepTimeFixedLsequen.fast, file = "GibbsStepTimeFixedLsequenFast.RData")
save(GibbsStepTimeVaryLjSequen.fast, file = "GibbsStepTimeVaryLjSequenFast.RData")
####### Reference: Appendix B
load("GibbsStepTimeFixedLsequen.RData")
load("GibbsStepTimeVaryLjSequen.RData")
load("GibbsStepTimeFixedLsequenFast.RData")
load("GibbsStepTimeVaryLjSequenFast.RData")













