rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste0(projDirec, "/simu/appen/sensAnalysis/diffHyperPara"))
setwd("setting1")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") 
Diags1 <- Diags
Deviance1 <- Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")
Diags.block1 <- Diags.block
Deviance.block1 <- Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequen1 <- Diags.sequen
Deviance.sequen1 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenVaryLj1 <- Diags.sequenVaryLj
Deviance.sequenVaryLj1 <- Deviance.sequenVaryLj
setwd("../setting2")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") 
Diags2 <- Diags
Deviance2 <- Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")
Diags.block2 <- Diags.block
Deviance.block2 <- Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequen2 <- Diags.sequen
Deviance.sequen2 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenVaryLj2 <- Diags.sequenVaryLj
Deviance.sequenVaryLj2 <- Deviance.sequenVaryLj
setwd("../setting3")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") 
Diags3 <- Diags
Deviance3 <- Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")
Diags.block3 <- Diags.block
Deviance.block3 <- Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequen3 <- Diags.sequen
Deviance.sequen3 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenVaryLj3 <- Diags.sequenVaryLj
Deviance.sequenVaryLj3 <- Deviance.sequenVaryLj
setwd("../setting4")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") 
Diags4 <- Diags
Deviance4 <- Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")
Diags.block4 <- Diags.block
Deviance.block4 <- Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequen4 <- Diags.sequen
Deviance.sequen4 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenVaryLj4 <- Diags.sequenVaryLj
Deviance.sequenVaryLj4 <- Deviance.sequenVaryLj
setwd("../setting5")
load("regFixedL50simuT100M225K5Iter100000Diags.RData")
load("regFixedL50simuT100M225K5Iter100000Deviance.RData") 
Diags5 <- Diags
Deviance5 <- Deviance
load("regFixedL50simuBlockT100M225K5Iter100000Diags.RData")
load("regFixedL50simuBlockT100M225K5Iter100000Deviance.RData")
Diags.block5 <- Diags.block
Deviance.block5 <- Deviance.block
load("regFixedL50simuSequenT100M225K5Iter100000Diags.RData")
load("regFixedL50simuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequen5 <- Diags.sequen
Deviance.sequen5 <- Deviance.sequen
load("regVaryLjSimuSequenT100M225K5Iter100000Diags.RData")
load("regVaryLjSimuSequenT100M225K5Iter100000Deviance.RData") 
Diags.sequenVaryLj5 <- Diags.sequenVaryLj
Deviance.sequenVaryLj5 <- Deviance.sequenVaryLj
setwd("..")
NKeep <- 5000
postDeviancesDF1 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance1), NNGPblockFixedL = as.vector(Deviance.block1),
                               NNGPsequenFixedL = as.vector(Deviance.sequen1), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj1))
fullGPfixedLpostDeviances1 <- ggplot(postDeviancesDF1) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances1
NNGPblockFixedLpostDeviances1 <- ggplot(postDeviancesDF1) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances1
NNGPsequenFixedLpostDeviances1 <- ggplot(postDeviancesDF1) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances1
NNGPsequenVaryLjpostDeviances1 <- ggplot(postDeviancesDF1) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances1
postDeviancesDF2 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance2), NNGPblockFixedL = as.vector(Deviance.block2),
                               NNGPsequenFixedL = as.vector(Deviance.sequen2), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj2))
fullGPfixedLpostDeviances2 <- ggplot(postDeviancesDF2) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances2
NNGPblockFixedLpostDeviances2 <- ggplot(postDeviancesDF2) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances2
NNGPsequenFixedLpostDeviances2 <- ggplot(postDeviancesDF2) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances2
NNGPsequenVaryLjpostDeviances2 <- ggplot(postDeviancesDF2) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances2
postDeviancesDF3 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance3), NNGPblockFixedL = as.vector(Deviance.block3),
                               NNGPsequenFixedL = as.vector(Deviance.sequen3), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj3))
fullGPfixedLpostDeviances3 <- ggplot(postDeviancesDF3) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances3
NNGPblockFixedLpostDeviances3 <- ggplot(postDeviancesDF3) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances3
NNGPsequenFixedLpostDeviances3 <- ggplot(postDeviancesDF3) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances3
NNGPsequenVaryLjpostDeviances3 <- ggplot(postDeviancesDF3) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances3
postDeviancesDF4 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance4), NNGPblockFixedL = as.vector(Deviance.block4),
                               NNGPsequenFixedL = as.vector(Deviance.sequen4), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj4))
fullGPfixedLpostDeviances4 <- ggplot(postDeviancesDF4) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances4
NNGPblockFixedLpostDeviances4 <- ggplot(postDeviancesDF4) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances4
NNGPsequenFixedLpostDeviances4 <- ggplot(postDeviancesDF4) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances4
NNGPsequenVaryLjpostDeviances4 <- ggplot(postDeviancesDF4) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances4
postDeviancesDF5 <- data.frame(MCMCiter = 1:NKeep, fullGPfixedL = as.vector(Deviance5), NNGPblockFixedL = as.vector(Deviance.block5),
                               NNGPsequenFixedL = as.vector(Deviance.sequen5), NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj5))
fullGPfixedLpostDeviances5 <- ggplot(postDeviancesDF5) + geom_line(aes(x = MCMCiter, y = fullGPfixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
fullGPfixedLpostDeviances5
NNGPblockFixedLpostDeviances5 <- ggplot(postDeviancesDF5) + geom_line(aes(x = MCMCiter, y = NNGPblockFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPblockFixedLpostDeviances5
NNGPsequenFixedLpostDeviances5 <- ggplot(postDeviancesDF5) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances5
NNGPsequenVaryLjpostDeviances5 <- ggplot(postDeviancesDF5) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances5
ggarrange(fullGPfixedLpostDeviances1, fullGPfixedLpostDeviances2, 
          fullGPfixedLpostDeviances3, fullGPfixedLpostDeviances4,
          fullGPfixedLpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("fullGPfixedLdiffHyperPara.png", width = 16, height = 24, units = "cm")
ggarrange(fullGPfixedLpostDeviances1, fullGPfixedLpostDeviances2, fullGPfixedLpostDeviances3, 
          fullGPfixedLpostDeviances4, fullGPfixedLpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("fullGPfixedLdiffHyperParaSmall.png", width = 24, height = 16, units = "cm")
ggarrange(NNGPblockFixedLpostDeviances1, NNGPblockFixedLpostDeviances2, 
          NNGPblockFixedLpostDeviances3, NNGPblockFixedLpostDeviances4,
          NNGPblockFixedLpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("NNGPblockFixedLdiffHyperPara.png", width = 16, height = 24, units = "cm")
ggarrange(NNGPblockFixedLpostDeviances1, NNGPblockFixedLpostDeviances2, NNGPblockFixedLpostDeviances3, 
          NNGPblockFixedLpostDeviances4, NNGPblockFixedLpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("NNGPblockFixedLdiffHyperParaSmall.png", width = 24, height = 16, units = "cm")
ggarrange(NNGPsequenFixedLpostDeviances1, NNGPsequenFixedLpostDeviances2, 
          NNGPsequenFixedLpostDeviances3, NNGPsequenFixedLpostDeviances4,
          NNGPsequenFixedLpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("NNGPsequenFixedLdiffHyperPara.png", width = 16, height = 24, units = "cm")
ggarrange(NNGPsequenFixedLpostDeviances1, NNGPsequenFixedLpostDeviances2, NNGPsequenFixedLpostDeviances3, 
          NNGPsequenFixedLpostDeviances4, NNGPsequenFixedLpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("NNGPsequenFixedLdiffHyperParaSmall.png", width = 24, height = 16, units = "cm")
ggarrange(NNGPsequenVaryLjpostDeviances1, NNGPsequenVaryLjpostDeviances2, 
          NNGPsequenVaryLjpostDeviances3, NNGPsequenVaryLjpostDeviances4,
          NNGPsequenVaryLjpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 2, nrow = 3)
ggsave("NNGPsequenVaryLjDiffHyperPara.png", width = 16, height = 24, units = "cm")
ggarrange(NNGPsequenVaryLjpostDeviances1, NNGPsequenVaryLjpostDeviances2, NNGPsequenVaryLjpostDeviances3, 
          NNGPsequenVaryLjpostDeviances4, NNGPsequenVaryLjpostDeviances5, ggplot() + theme_void(),
          labels = c("A", "B", "C", "D", "E", ""), align = "hv",
          ncol = 3, nrow = 2)
ggsave("NNGPsequenVaryLjDiffHyperParaSmall.png", width = 24, height = 16, units = "cm")
