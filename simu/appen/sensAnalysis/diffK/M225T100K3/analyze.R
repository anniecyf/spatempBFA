rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/sensAnalysis/diffK"))
setwd("M225T100K3")
#load("regFixedL50simuT100M225K3Iter90000.RData"); regFixedL.simu$runtime; # slower package version 3.89 hours 3.89 / 9 * 10 \approx 4.322222 hours rm(regFixedL.simu)
load("regFixedL50simuT100M225K3Iter90000Diags.RData")
load("regFixedL50simuT100M225K3Iter90000Deviance.RData") # Deviance
load("GibbsStepTimeFixedLfullGP.RData")
#load("regFixedL50simuBlockT100M225K3Iter100000.RData"); regFixedL.simu.block$runtime; # slower package version 6.13 hours rm(regFixedL.simu.block)
load("regFixedL50simuBlockT100M225K3Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K3Iter100000Deviance.RData") # Deviance.block
load("GibbsStepTimeFixedLblock.RData")
#load("regFixedL50simuSequenT100M225K3Iter100000.RData"); regFixedL.simu.sequen$runtime; # slower package version 6.89 hours rm(regFixedL.simu.sequen)
load("regFixedL50simuSequenT100M225K3Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K3Iter100000Deviance.RData") # Deviance.sequen
load("GibbsStepTimeFixedLsequen.RData")
#load("regVaryLjSimuSequenT100M225K3Iter100000.RData"); regVaryLj.simu.sequen$runtime; # slower package version 4.72 hours rm(regVaryLj.simu.sequen)
load("regVaryLjSimuSequenT100M225K3Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K3Iter100000Deviance.RData") # Deviance.sequenVaryLj
load("GibbsStepTimeVaryLjSequen.RData")
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
ggsave("diffK3postDeviances.png", width = 16, height = 16, units = "cm")
apply(GibbsStepTimeFixedLfullGP, 2, summary)
apply(GibbsStepTimeFixedLblock, 2, summary)
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)