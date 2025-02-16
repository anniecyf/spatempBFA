rm(list=ls())
library(tidyverse)
library(ggpubr)
modelVec <- rep(c("fullGP", "NNGPblock", "NNGPsequen", "NNGPvaryLj"), each = 4)
paraVec <- rep(c("       z / u      ", "         xi        ", "     alpha      ", "kappa & rho"), 4)
m3600T100DF <- m3600T50DF <- m1600T100DF <- m1600T50DF <- m1600T30DF <- m900T50DF <- m900T30DF <-
  m400T50DF <- m400T30DF <- data.frame(modelVec = modelVec, paraVec = paraVec)
legendColors <- c("#FFAA11","#55CC11", "#1177CC", "#7755CC")
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/mainScalabilityVerificationSimu/m3600T100K5", sep = "/"))
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m3600T100DF = mutate(m3600T100DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m3600T100plot <- ggplot(m3600T100DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m3600T100DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m3600T100DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m3600T100DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") + 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m3600T100plot
setwd("../m3600T50K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m3600T50DF = mutate(m3600T50DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m3600T50plot <- ggplot(m3600T50DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m3600T50DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m3600T50DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m3600T50DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") + 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m3600T50plot
setwd("../m1600T100K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m1600T100DF = mutate(m1600T100DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m1600T100plot <- ggplot(m1600T100DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m1600T100DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m1600T100DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m1600T100DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m1600T100plot
setwd("../m1600T50K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m1600T50DF = mutate(m1600T50DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m1600T50plot <- ggplot(m1600T50DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m1600T50DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m1600T50DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m1600T50DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m1600T50plot
setwd("../m1600T30K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m1600T30DF = mutate(m1600T30DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m1600T30plot <- ggplot(m1600T30DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m1600T30DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m1600T30DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m1600T30DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m1600T30plot
setwd("../m900T50K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m900T50DF = mutate(m900T50DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m900T50plot <- ggplot(m900T50DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m900T50DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m900T50DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m900T50DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m900T50plot
setwd("../m900T30K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m900T30DF = mutate(m900T30DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m900T30plot <- ggplot(m900T30DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m900T30DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m900T30DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m900T30DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m900T30plot
setwd("../m400T50K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m400T50DF = mutate(m400T50DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m400T50plot <- ggplot(m400T50DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m400T50DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m400T50DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m400T50DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14))  
m400T50plot
setwd("../m400T30K5")
load("GibbsStepTimeFixedLfullGP.RData") # GibbsStepTimeFixedLfullGP
load("GibbsStepTimeFixedLblock.RData") # GibbsStepTimeFixedLblock
load("GibbsStepTimeFixedLsequen.RData") # GibbsStepTimeFixedLsequen
load("GibbsStepTimeVaryLjSequen.RData") # GibbsStepTimeVaryLjSequen
meanTimeVec <- vector()
meanTime <- colMeans(GibbsStepTimeFixedLfullGP)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLblock)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeFixedLsequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
meanTime <- colMeans(GibbsStepTimeVaryLjSequen)
meanTimeVec <- c(meanTimeVec, meanTime[c(1, 2, 5)])
meanTimeVec <- c(meanTimeVec, meanTime[7] + meanTime[6])
names(meanTimeVec)[length(meanTimeVec)] <- "rhoKappa"
m400T30DF = mutate(m400T30DF, meanTimeVec = meanTimeVec / 1000) # in seconds
m400T30plot <- ggplot(m400T30DF, aes(x = modelVec, y = meanTimeVec, fill = paraVec, group = paraVec)) +
  labs(x = "", y = "Time (sec)") +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = legendColors, 
                    guide = guide_legend(title = "Parameter", label.position = "top", nrow = 1)) +
  geom_area(data = m400T30DF %>% filter(modelVec %in% c("fullGP", "NNGPblock")),
            aes(x = c("fullGP" = 1.25, "NNGPblock" = 1.75)[modelVec]), 
            position = "stack", color = "black", fill = NA, outline.type = "both") +
  geom_area(data = m400T30DF %>% filter(modelVec %in% c("NNGPblock", "NNGPsequen")),
            aes(x = c("NNGPblock" = 2.25, "NNGPsequen" = 2.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both") +
  geom_area(data = m400T30DF %>% filter(modelVec %in% c("NNGPsequen", "NNGPvaryLj")),
            aes(x = c("NNGPsequen" = 3.25, "NNGPvaryLj" = 3.75)[modelVec]), 
            position = "stack", color = "black",fill = NA, outline.type = "both")+ 
  theme(legend.key.size = unit(1.5, "cm"),     
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16) )  
m400T30plot
setwd("..")
# ggarrange(m1600T100plot, m1600T50plot, m1600T30plot,
#           m3600T100plot, m900T50plot, m400T50plot,
#           m3600T50plot, m900T30plot, m400T30plot,
#           labels = c("A", "B", "C", "D", "F", "H", "E", "G", "I"), align = "hv",
#           ncol = 3, nrow = 3, common.legend = TRUE, legend = "top")
# ggsave("9stackedPlots.png", width = 50, height = 50, units = "cm")
ggarrange(m400T30plot, m400T50plot, m900T30plot,
          m900T50plot, m1600T30plot, m1600T50plot, 
          m1600T100plot, m3600T50plot, m3600T100plot,  
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), align = "hv",
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "top")
# ggsave("stackedPlots.png", width = 50, height = 52, units = "cm")

