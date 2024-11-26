rm(list=ls())
library(tidyverse)
library(ggpubr)
legendColors <- c("#FFAA11", "#55CC11", "#1177CC", "#7755CC")
spatpredTime <- alphaKrigTime <- weightsXiLambdaKrigTime <- vector()
load("fullGPfixedLspatpredTimeVec.RData")
load("fullGPfixedLalphaKrigTimeVec.RData")
load("fullGPfixedLweightsXiLambdaKrigTimeVec.RData")
N <- length(spatpredTimeVec)
modelVec <- rep(c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj"), each = N)
spatpredTime <- c(spatpredTime, spatpredTimeVec)
alphaKrigTime <- c(alphaKrigTime, alphaKrigTimeVec)
weightsXiLambdaKrigTime <- c(weightsXiLambdaKrigTime, weightsXiLambdaKrigTimeVec)
load("NNGPblockFixedLspatpredTimeVec.RData")
load("NNGPblockFixedLalphaKrigTimeVec.RData")
load("NNGPblockFixedLweightsXiLambdaKrigTimeVec.RData")
spatpredTime <- c(spatpredTime, spatpredTimeVec)
alphaKrigTime <- c(alphaKrigTime, alphaKrigTimeVec)
weightsXiLambdaKrigTime <- c(weightsXiLambdaKrigTime, weightsXiLambdaKrigTimeVec)
load("NNGPsequenFixedLspatpredTimeVec.RData")
load("NNGPsequenFixedLalphaKrigTimeVec.RData")
load("NNGPsequenFixedLweightsXiLambdaKrigTimeVec.RData")
spatpredTime <- c(spatpredTime, spatpredTimeVec)
alphaKrigTime <- c(alphaKrigTime, alphaKrigTimeVec)
weightsXiLambdaKrigTime <- c(weightsXiLambdaKrigTime, weightsXiLambdaKrigTimeVec)
load("NNGPsequenVaryLjspatpredTimeVec.RData")
load("NNGPsequenVaryLjalphaKrigTimeVec.RData")
load("NNGPsequenVaryLjweightsXiLambdaKrigTimeVec.RData")
spatpredTime <- c(spatpredTime, spatpredTimeVec)
alphaKrigTime <- c(alphaKrigTime, alphaKrigTimeVec)
weightsXiLambdaKrigTime <- c(weightsXiLambdaKrigTime, weightsXiLambdaKrigTimeVec)
spatpredtimeDF <- data.frame(spatpredTime = spatpredTime, alphaKrigTime = alphaKrigTime / 1000,
                             weightsXiLambdaKrigTime = weightsXiLambdaKrigTime,
                             model = as.factor(modelVec))
spatpredtimeBox <- ggplot(spatpredtimeDF) + labs(y = "", x = "Time (in seconds)") + 
  geom_boxplot(aes(x = spatpredTime, y = model, fill = model)) +
  scale_fill_manual("Model", values = legendColors) +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top", legend.title = element_text(face = "bold", size = 12),
        legend.key = element_blank())
spatpredtimeBox
alphaKrigTimeBox <- ggplot(spatpredtimeDF) + labs(y = "", x = "Time (in seconds)") + 
  geom_boxplot(aes(x = alphaKrigTime, y = model, fill = model)) +
  scale_fill_manual("Model", values = legendColors) +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top", legend.title = element_text(face = "bold", size = 12),
        legend.key = element_blank())
alphaKrigTimeBox
weightsXiLambdaKrigTimeBox <- ggplot(spatpredtimeDF) + 
  labs(y = "", x = "Time (in milliseconds)") + 
  geom_boxplot(aes(x = weightsXiLambdaKrigTime, y = model, fill = model)) +
  scale_fill_manual("Model", values = legendColors) +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top", legend.title = element_text(face = "bold", size = 12),
        legend.key = element_blank())
weightsXiLambdaKrigTimeBox
combinedSpatPredBox <- ggarrange(spatpredtimeBox, alphaKrigTimeBox, weightsXiLambdaKrigTimeBox,
                                 labels = c("A", "B", "C"), align = "h",
                                 ncol = 3, nrow = 1, common.legend = TRUE, legend = "top")
combinedSpatPredBox
# ggsave("combinedSpatPredBoxCorrected.png", width = 24, height = 8, units = "cm")
spatpredtimeBoxBW <- ggplot(spatpredtimeDF) + labs(y = "", x = "Time (in seconds)") + 
  geom_boxplot(aes(x = spatpredTime, y = model)) +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
spatpredtimeBoxBW
alphaKrigTimeBoxBW <- ggplot(spatpredtimeDF) + labs(y = "", x = "Time (in seconds)") + 
  geom_boxplot(aes(x = alphaKrigTime, y = model)) +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
alphaKrigTimeBoxBW
weightsXiLambdaKrigTimeBoxBW <- ggplot(spatpredtimeDF) + 
  labs(y = "", x = "Time (in milliseconds)") + 
  geom_boxplot(aes(x = weightsXiLambdaKrigTime, y = model)) +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
weightsXiLambdaKrigTimeBoxBW
combinedSpatPredBoxBW <- ggarrange(spatpredtimeBoxBW, alphaKrigTimeBoxBW, 
                                   weightsXiLambdaKrigTimeBoxBW,
                                   labels = c("A", "B", "C"), align = "h",
                                   ncol = 3, nrow = 1, widths = c(3, 2.15, 2.15))
combinedSpatPredBoxBW
# ggsave("combinedSpatPredBoxBWcorrected.png", width = 24, height = 8, units = "cm")
spatpredtimeBoxBW <- ggplot(spatpredtimeDF) + labs(y = "Time (in seconds)", x = "") + 
  geom_boxplot(aes(x = model, y = spatpredTime)) +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
spatpredtimeBoxBW
alphaKrigTimeBoxBW <- ggplot(spatpredtimeDF) + labs(y = "Time (in seconds)", x = "") + 
  geom_boxplot(aes(x = model, y = alphaKrigTime)) +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
alphaKrigTimeBoxBW
weightsXiLambdaKrigTimeBoxBW <- ggplot(spatpredtimeDF) + 
  labs(y = "Time (in milliseconds)", x = "") + 
  geom_boxplot(aes(x = model, y = weightsXiLambdaKrigTime)) +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
weightsXiLambdaKrigTimeBoxBW
combinedSpatPredBoxBW <- ggarrange(spatpredtimeBoxBW, alphaKrigTimeBoxBW, 
                                   weightsXiLambdaKrigTimeBoxBW,
                                   labels = c("A", "B", "C"), align = "h",
                                   ncol = 3, nrow = 1)
combinedSpatPredBoxBW
#ggsave("correctedCombinedSpatPredBoxBW.png", width = 39, height = 13, units = "cm")

spatpredMetricMat <- matrix(0, 4, 3, 
                            dimnames = list(c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj"),
                                            c("meanMedianPostMeanMSE", "meanMedianPostMSE", "meanMedianPostVar")))
load("fullGPfixedLpostGmat_spatPredMetricSimu.RData")
load("fullGPfixedLpostMSEmat_spatPredMetricSimu.RData")
load("fullGPfixedLpostVarMat_spatPredMetricSimu.RData")
spatpredMetricMat[1, 1] <- mean(apply(postGmat, 2, median))
spatpredMetricMat[1, 2] <- mean(apply(postMSEmat, 2, median))
spatpredMetricMat[1, 3] <- mean(apply(postVarMat, 2, median))
load("NNGPblockFixedLpostGmat_spatPredMetricSimu.RData")
load("NNGPblockFixedLpostMSEmat_spatPredMetricSimu.RData")
load("NNGPblockFixedLpostVarMat_spatPredMetricSimu.RData")
spatpredMetricMat[2, 1] <- mean(apply(postGmat, 2, median))
spatpredMetricMat[2, 2] <- mean(apply(postMSEmat, 2, median))
spatpredMetricMat[2, 3] <- mean(apply(postVarMat, 2, median))
load("NNGPsequenFixedLpostGmat_spatPredMetricSimu.RData")
load("NNGPsequenFixedLpostMSEmat_spatPredMetricSimu.RData")
load("NNGPsequenFixedLpostVarMat_spatPredMetricSimu.RData")
spatpredMetricMat[3, 1] <- mean(apply(postGmat, 2, median))
spatpredMetricMat[3, 2] <- mean(apply(postMSEmat, 2, median))
spatpredMetricMat[3, 3] <- mean(apply(postVarMat, 2, median))
load("NNGPsequenVaryLjpostGmat_spatPredMetricSimu.RData")
load("NNGPsequenVaryLjpostMSEmat_spatPredMetricSimu.RData")
load("NNGPsequenVaryLjpostVarMat_spatPredMetricSimu.RData")
spatpredMetricMat[4, 1] <- mean(apply(postGmat, 2, median))
spatpredMetricMat[4, 2] <- mean(apply(postMSEmat, 2, median))
spatpredMetricMat[4, 3] <- mean(apply(postVarMat, 2, median))
#save(spatpredMetricMat, file = "spatpredMetric.RData")
spatpredMetricMat
# meanMedianPostMeanMSE meanMedianPostMSE meanMedianPostVar
# fullGPfixedL                 0.4780934          5.964139          5.502430
# NNGPblockFixedL              0.4815116          5.965555          5.499377
# NNGPsequenFixedL             0.4766795          5.962761          5.495316
# NNGPsequenVaryLj             6.1963940         16.034109          9.832205
