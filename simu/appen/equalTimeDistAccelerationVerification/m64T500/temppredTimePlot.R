rm(list=ls())
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/equalTimeDistAccelerationVerification/m64T500", sep = "/"))
library(tidyverse)
library(ggpubr)
legendColors <- c("#1177CC", "#55CC11")
load("m64T500temppredTimeMat.RData")
apply(temppredTimeMat, 2, summary)
N <- nrow(temppredTimeMat)
temppredTimeDF <- data.frame(temppredTime = as.vector(temppredTimeMat), 
                             model = as.factor(rep(colnames(temppredTimeMat), each = N)),
                             equalTimeDistTF = as.factor(rep(rep(c("equalTimeDist = TRUE", "equalTimeDist = FALSE"), each = N), 4)))
temppredtimeBox <- ggplot(temppredTimeDF) + labs(y = "", x = "Time (in seconds)") + 
  geom_boxplot(aes(x = temppredTime, y = model, fill = equalTimeDistTF)) +
  scale_fill_manual("", values = legendColors) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top", legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 12))
temppredtimeBox
ggsave("tempPredBox.png", width = 16, height = 10, units = "cm")
