rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/multObsTypes/O2T20M100Iter50000"))
load("regFixedL30simuO2T20M100Iter50000Diags.RData")
load("regFixedL30simuO2T20M100Iter50000Deviance.RData")
load("GibbsStepTimeFixedLfullGP.RData")
# load("regFixedL30simuO2T20M100Iter50000.RData"); regFixedL.simu$runtime # slower package version 1.33 hours; rm(regFixedL.simu)
load("regFixedL30simuBlockO2T20M100Iter50000Diags.RData")
load("regFixedL30simuBlockO2T20M100Iter50000Deviance.RData")
load("GibbsStepTimeFixedLblock.RData")
# load("regFixedL30simuBlockO2T20M100Iter50000.RData"); regFixedL.simu.block$runtime # slower package version 1.34 hours; rm(regFixedL.simu.block)
load("regFixedL30simuSequenO2T20M100Iter50000Diags.RData")
load("regFixedL30simuSequenO2T20M100Iter50000Deviance.RData")
load("GibbsStepTimeFixedLsequen.RData")
# load("regFixedL30simuSequenO2T20M100Iter50000.RData"); regFixedL.simu.sequen$runtime # slower package version 2.06 hours; rm(regFixedL.simu.sequen)
load("regVaryLjSimuSequenO2T20M100Iter50000Diags.RData")
load("regVaryLjSimuSequenO2T20M100Iter50000Deviance.RData")
load("GibbsStepTimeVaryLjSequen.RData")
# load("regVaryLjSimuSequenO2T20M100Iter50000.RData"); regVaryLj.simu.sequen$runtime # slower package version 8.68 hours; rm(regVaryLj.simu.sequen)
NKeep <- 10000
postDeviancesDF <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance), NNGPblockFixedL = as.vector(Deviance.block),
                              NNGPsequenFixedL = as.vector(Deviance.sequen), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj))
# fullGPfixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = fullGPfixedL), 
#                                        col = rgb(0.27, 0.35, 0.27)) + 
#   labs(x = "MCMC Iteration", y = "Posterior Deviance") +
#   theme(axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         axis.text.x = element_text(size = 10, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(color = "black", fill=NA))
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
ggsave("O2T20M100postDeviances.png", width = 16, height = 16, units = "cm")
apply(GibbsStepTimeFixedLfullGP, 2, summary)
apply(GibbsStepTimeFixedLblock, 2, summary)
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)
# Diags.sequen (1) better than (except for dic, pd, waic, p_waic) Diags (3) better than (except for waic, p_waic, p_waic_1) Diags.block (4)
# Diags.sequenVaryLj (2) better than (except for waic, p_waic, lppd) Diags and worse than (except for dic, pd, p_waic_1) Diags.sequen
# all comparable
# the posterior deviances plots (all converged) also suggest so (NNGPsequenFixedL better than NNGPsequenVaryLj better than fullGPfixedL better than NNGPblockFixedL)