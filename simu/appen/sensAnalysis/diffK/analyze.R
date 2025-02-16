rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/sensAnalysis/diffK"))
setwd("M225T100K3")
load("regFixedL50simuT100M225K3Iter90000Diags.RData")
load("regFixedL50simuT100M225K3Iter90000Deviance.RData") 
DiagsK3 <- Diags
DevianceK3 <- Deviance
load("regFixedL50simuBlockT100M225K3Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K3Iter100000Deviance.RData")
Diags.blockK3 <- Diags.block
Deviance.blockK3 <- Deviance.block
load("regFixedL50simuSequenT100M225K3Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K3Iter100000Deviance.RData") 
Diags.sequenK3 <- Diags.sequen
Deviance.sequenK3 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K3Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K3Iter100000Deviance.RData") 
Diags.sequenVaryLjK3 <- Diags.sequenVaryLj
Deviance.sequenVaryLjK3 <- Deviance.sequenVaryLj
setwd("../M225T100K4")
load("regFixedL50simuT100M225K4Iter100000Diags.RData")
load("regFixedL50simuT100M225K4Iter100000Deviance.RData") 
DiagsK4 <- Diags
DevianceK4 <- Deviance
load("regFixedL50simuBlockT100M225K4Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K4Iter100000Deviance.RData")
Diags.blockK4 <- Diags.block
Deviance.blockK4 <- Deviance.block
load("regFixedL50simuSequenT100M225K4Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K4Iter100000Deviance.RData") 
Diags.sequenK4 <- Diags.sequen
Deviance.sequenK4 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K4Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K4Iter100000Deviance.RData") 
Diags.sequenVaryLjK4 <- Diags.sequenVaryLj
Deviance.sequenVaryLjK4 <- Deviance.sequenVaryLj
setwd("../M225T100K5")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") 
DiagsK5 <- Diags
DevianceK5 <- Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")
Diags.blockK5 <- Diags.block
Deviance.blockK5 <- Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenK5 <- Diags.sequen
Deviance.sequenK5 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenVaryLjK5 <- Diags.sequenVaryLj
Deviance.sequenVaryLjK5 <- Deviance.sequenVaryLj
setwd("../M225T100K6")
load("regFixedL50simuT100M225K6Iter100000Diags.RData")
load("regFixedL50simuT100M225K6Iter100000Deviance.RData") 
DiagsK6 <- Diags
DevianceK6 <- Deviance
load("regFixedL50simuBlockT100M225K6Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K6Iter100000Deviance.RData")
Diags.blockK6 <- Diags.block
Deviance.blockK6 <- Deviance.block
load("regFixedL50simuSequenT100M225K6Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K6Iter100000Deviance.RData") 
Diags.sequenK6 <- Diags.sequen
Deviance.sequenK6 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K6Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K6Iter100000Deviance.RData") 
Diags.sequenVaryLjK6 <- Diags.sequenVaryLj
Deviance.sequenVaryLjK6 <- Deviance.sequenVaryLj
setwd("../M225T100K7")
load("regFixedL50simuT100M225K7Iter100000Diags.RData")
load("regFixedL50simuT100M225K7Iter100000Deviance.RData") 
DiagsK7 <- Diags
DevianceK7 <- Deviance
load("regFixedL50simuBlockT100M225K7Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K7Iter100000Deviance.RData")
Diags.blockK7 <- Diags.block
Deviance.blockK7 <- Deviance.block
load("regFixedL50simuSequenT100M225K7Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K7Iter100000Deviance.RData") 
Diags.sequenK7 <- Diags.sequen
Deviance.sequenK7 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K7Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K7Iter100000Deviance.RData") 
Diags.sequenVaryLjK7 <- Diags.sequenVaryLj
Deviance.sequenVaryLjK7 <- Deviance.sequenVaryLj
setwd("..")

NKeep <- 5000
postDeviancesDFK3 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(DevianceK3), NNGPblockFixedL = as.vector(Deviance.blockK3),
                               NNGPsequenFixedL = as.vector(Deviance.sequenK3), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLjK3))
fullGPfixedLpostDeviancesK3 <- ggplot(postDeviancesDFK3) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviancesK3
NNGPblockFixedLpostDeviancesK3 <- ggplot(postDeviancesDFK3) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviancesK3
NNGPsequenFixedLpostDeviancesK3 <- ggplot(postDeviancesDFK3) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviancesK3
NNGPsequenVaryLjpostDeviancesK3 <- ggplot(postDeviancesDFK3) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviancesK3
postDeviancesDFK4 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(DevianceK4), NNGPblockFixedL = as.vector(Deviance.blockK4),
                               NNGPsequenFixedL = as.vector(Deviance.sequenK4), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLjK4))
fullGPfixedLpostDeviancesK4 <- ggplot(postDeviancesDFK4) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviancesK4
NNGPblockFixedLpostDeviancesK4 <- ggplot(postDeviancesDFK4) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviancesK4
NNGPsequenFixedLpostDeviancesK4 <- ggplot(postDeviancesDFK4) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviancesK4
NNGPsequenVaryLjpostDeviancesK4 <- ggplot(postDeviancesDFK4) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviancesK4
postDeviancesDFK5 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(DevianceK5), NNGPblockFixedL = as.vector(Deviance.blockK5),
                               NNGPsequenFixedL = as.vector(Deviance.sequenK5), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLjK5))
fullGPfixedLpostDeviancesK5 <- ggplot(postDeviancesDFK5) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviancesK5
NNGPblockFixedLpostDeviancesK5 <- ggplot(postDeviancesDFK5) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviancesK5
NNGPsequenFixedLpostDeviancesK5 <- ggplot(postDeviancesDFK5) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviancesK5
NNGPsequenVaryLjpostDeviancesK5 <- ggplot(postDeviancesDFK5) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviancesK5
postDeviancesDFK6 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(DevianceK6), NNGPblockFixedL = as.vector(Deviance.blockK6),
                               NNGPsequenFixedL = as.vector(Deviance.sequenK6), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLjK6))
fullGPfixedLpostDeviancesK6 <- ggplot(postDeviancesDFK6) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviancesK6
NNGPblockFixedLpostDeviancesK6 <- ggplot(postDeviancesDFK6) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviancesK6
NNGPsequenFixedLpostDeviancesK6 <- ggplot(postDeviancesDFK6) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviancesK6
NNGPsequenVaryLjpostDeviancesK6 <- ggplot(postDeviancesDFK6) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviancesK6
postDeviancesDFK7 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(DevianceK7), NNGPblockFixedL = as.vector(Deviance.blockK7),
                               NNGPsequenFixedL = as.vector(Deviance.sequenK7), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLjK7))
fullGPfixedLpostDeviancesK7 <- ggplot(postDeviancesDFK7) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviancesK7
NNGPblockFixedLpostDeviancesK7 <- ggplot(postDeviancesDFK7) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviancesK7
NNGPsequenFixedLpostDeviancesK7 <- ggplot(postDeviancesDFK7) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviancesK7
NNGPsequenVaryLjpostDeviancesK7 <- ggplot(postDeviancesDFK7) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviancesK7
ggarrange(fullGPfixedLpostDeviancesK3, fullGPfixedLpostDeviancesK4, 
          fullGPfixedLpostDeviancesK5, fullGPfixedLpostDeviancesK6,
          fullGPfixedLpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("fullGPfixedLdiffK.png", width = 16, height = 24, units = "cm")
ggarrange(fullGPfixedLpostDeviancesK3, fullGPfixedLpostDeviancesK4, fullGPfixedLpostDeviancesK5, 
          fullGPfixedLpostDeviancesK6, fullGPfixedLpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("fullGPfixedLdiffKsmall.png", width = 24, height = 16, units = "cm")
ggarrange(NNGPblockFixedLpostDeviancesK3, NNGPblockFixedLpostDeviancesK4, 
          NNGPblockFixedLpostDeviancesK5, NNGPblockFixedLpostDeviancesK6,
          NNGPblockFixedLpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("NNGPblockFixedLdiffK.png", width = 16, height = 24, units = "cm")
ggarrange(NNGPblockFixedLpostDeviancesK3, NNGPblockFixedLpostDeviancesK4, NNGPblockFixedLpostDeviancesK5, 
          NNGPblockFixedLpostDeviancesK6, NNGPblockFixedLpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("NNGPblockFixedLdiffKsmall.png", width = 24, height = 16, units = "cm")
ggarrange(NNGPsequenFixedLpostDeviancesK3, NNGPsequenFixedLpostDeviancesK4, 
          NNGPsequenFixedLpostDeviancesK5, NNGPsequenFixedLpostDeviancesK6,
          NNGPsequenFixedLpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("NNGPsequenFixedLdiffK.png", width = 16, height = 24, units = "cm")
ggarrange(NNGPsequenFixedLpostDeviancesK3, NNGPsequenFixedLpostDeviancesK4, NNGPsequenFixedLpostDeviancesK5, 
          NNGPsequenFixedLpostDeviancesK6, NNGPsequenFixedLpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("NNGPsequenFixedLdiffKsmall.png", width = 24, height = 16, units = "cm")
ggarrange(NNGPsequenVaryLjpostDeviancesK3, NNGPsequenVaryLjpostDeviancesK4, 
          NNGPsequenVaryLjpostDeviancesK5, NNGPsequenVaryLjpostDeviancesK6,
          NNGPsequenVaryLjpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("NNGPsequenVaryLjDiffK.png", width = 16, height = 24, units = "cm")
ggarrange(NNGPsequenVaryLjpostDeviancesK3, NNGPsequenVaryLjpostDeviancesK4, NNGPsequenVaryLjpostDeviancesK5, 
          NNGPsequenVaryLjpostDeviancesK6, NNGPsequenVaryLjpostDeviancesK7, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("NNGPsequenVaryLjDiffKsmall.png", width = 24, height = 16, units = "cm")
# it seems that in general, a smaller k leads to better comprehensive Diags for fullGPfixedL and NNGPsequenFixedL (comparable for K = 3 and K = 4 fullGPfixedL)
# for NNGPblockFixedL Diags, k = 3 better than k = 4, k = 5 better than k = 4, k = 5 better than k = 6 better than k = 7; agrees with posterior deviances plots
# NNGPblockFixedL Diags continued k = 3 (1) better than k = 5 (2); k = 6 (3) better than k = 4 (4) better than k = 7 (5)
# for NNGPsequenVaryLj Diags, k = 3 better than k = 4 better than k = 5; k = 7, k = 6 better than k = 5; k = 7 slightly better than k = 6 overall (worse for dic, pd, p_waic_1); 
# NNGPsequenVaryLj Diags continued k = 4 better than k = 6 (except for p_waic and p_waic_1) and k = 7 (except for p_waic); agrees with posterior deviances plots