rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/VAR1 demo/VAR1demo/more burnin iter"))
#load("regFixedL50simuT200M225Iter100000.RData"); VAR1bfaFixedL.simu$runtime; # 23.99 hours rm(VAR1bfaFixedL.simu)
load("regFixedL50simuT200M225Iter100000Diags.RData") # DiagsVAR1
load("regFixedL50simuT200M225Iter100000Deviance.RData") # DevianceVAR1
#load("regFixedL50simuBlockT200M225Iter100000.RData"); VAR1bfaFixedL.simu.block$runtime; # 1.05 days rm(VAR1bfaFixedL.simu.block)
load("regFixedL50simuBlockT200M225Iter100000Diags.RData") # Diags.VAR1fixedL.block
load("regFixedL50simuBlockT200M225Iter100000Deviance.RData") # Deviance.VAR1fixedL.block
#load("regFixedL50simuSequenT200M225Iter100000.RData"); VAR1bfaFixedL.simu.sequen$runtime; # 1.05 days rm(VAR1bfaFixedL.simu.sequen)
load("regFixedL50simuSequenT200M225Iter100000Diags.RData") # Diags.VAR1fixedL.sequen
load("regFixedL50simuSequenT200M225Iter100000Deviance.RData") # Deviance.VAR1fixedL.sequen
#load("regVaryLjSimuSequenT200M225Iter100000.RData"); VAR1regVaryLj.simu.sequen$runtime; # 18.24 hours rm(VAR1regVaryLj.simu.sequen)
load("regVaryLjSimuSequenT200M225Iter100000Diags.RData") # Diags.VAR1regVaryLj.simu.sequen
load("regVaryLjSimuSequenT200M225Iter100000Deviance.RData") # Deviance.VAR1regVaryLj.simu.sequen
NKeep <- 10000
postDeviancesDF <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(DevianceVAR1), 
                              NNGPblockFixedL = as.vector(Deviance.VAR1fixedL.block),
                              NNGPsequenFixedL = as.vector(Deviance.VAR1fixedL.sequen), 
                              NNGPsequenVaryLj = as.vector(Deviance.VAR1sequenVaryLj.sequen))
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
ggsave("VAR1T200M225postDeviancesIter100000.png", width = 16, height = 16, units = "cm")
# diagnostics all comparable; fullGPfixedL better than NNGPblockFixedL and NNGPsequenFixedL better than (except for pd, p_waic, p_waic_1) NNGPsequenVaryLj 
# posterior deviances fullGPfixedL better than NNGPblockFixedL and NNGPsequenFixedL significantly better than NNGPsequenVaryLj