rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/multObsTypes/newest package version/O2T20M225Iter50000"))
load("regFixedL30simuO2T20M225Iter50000Diags.RData")
load("regFixedL30simuO2T20M225Iter50000Deviance.RData")
load("GibbsStepTimeFixedLfullGP.RData")
# load("regFixedL30simuO2T20M225Iter50000.RData"); regFixedL.simu$runtime # 3.75 hours; rm(regFixedL.simu)
load("regFixedL30simuBlockO2T20M225Iter50000Diags.RData")
load("regFixedL30simuBlockO2T20M225Iter50000Deviance.RData")
load("GibbsStepTimeFixedLblock.RData")
# load("regFixedL30simuBlockO2T20M225Iter50000.RData"); regFixedL.simu.block$runtime # 3.28 hours; rm(regFixedL.simu.block)
load("regFixedL30simuSequenO2T20M225Iter50000Diags.RData")
load("regFixedL30simuSequenO2T20M225Iter50000Deviance.RData")
load("GibbsStepTimeFixedLsequen.RData")
# load("regFixedL30simuSequenO2T20M225Iter50000.RData"); regFixedL.simu.sequen$runtime # 4.6 hours; rm(regFixedL.simu.sequen)
load("regVaryLjSimuSequenO2T20M225Iter50000Diags.RData")
load("regVaryLjSimuSequenO2T20M225Iter50000Deviance.RData")
load("GibbsStepTimeVaryLjSequen.RData")
# load("regVaryLjSimuSequenO2T20M225Iter50000.RData"); regVaryLj.simu.sequen$runtime #   hours; rm(regVaryLj.simu.sequen)
NKeep <- 5000
postDeviancesDF <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance), NNGPblockFixedL = as.vector(Deviance.block),
                               NNGPsequenFixedL = as.vector(Deviance.sequen))
# postDeviancesDF <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance), NNGPblockFixedL = as.vector(Deviance.block),
#                               NNGPsequenFixedL = as.vector(Deviance.sequen), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj))
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
ggarrange(fullGPfixedLpostDeviances, NNGPblockFixedLpostDeviances, NNGPsequenFixedLpostDeviances, 
          labels = c("A", "B", "C"), align = "h",
          ncol = 3, nrow = 1)
ggsave("O2T20M225postDeviancesFixedL.png", width = 24, height = 8, units = "cm")
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
ggsave("O2T20M225postDeviances.png", width = 16, height = 16, units = "cm")
apply(GibbsStepTimeFixedLfullGP, 2, summary)
apply(GibbsStepTimeFixedLblock, 2, summary)
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)