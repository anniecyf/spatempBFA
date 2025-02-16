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
spatpredtimeBox <- ggplot(spatpredtimeDF) + labs(y = "Time (in seconds)", x = "") + 
  geom_boxplot(aes(x = model, y = spatpredTime, color = model)) +
  scale_color_manual("Model", values = legendColors) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16, color = "black"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.key = element_blank(),    
        legend.text = element_text(size = 18), 
        legend.title = element_blank())
spatpredtimeBox
alphaKrigTimeBox <- ggplot(spatpredtimeDF) + labs(y = "Time (in seconds)", x = "") + 
  geom_boxplot(aes(x = model, y = alphaKrigTime, color = model)) +
  scale_color_manual("Model", values = legendColors) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.key = element_blank(),    
        legend.text = element_text(size = 18), 
        legend.title = element_blank())
alphaKrigTimeBox
weightsXiLambdaKrigTimeBox <- ggplot(spatpredtimeDF) + 
  labs(y = "Time (in milliseconds)", x = "") + 
  geom_boxplot(aes(x = model, y = weightsXiLambdaKrigTime, color = model)) +
  scale_color_manual("Model", values = legendColors) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.key = element_blank(),    
        legend.text = element_text(size = 18), 
        legend.title = element_blank())
weightsXiLambdaKrigTimeBox
combinedSpatPredBox <- ggarrange(spatpredtimeBox, alphaKrigTimeBox, 
                                   weightsXiLambdaKrigTimeBox,
                                   labels = c("A", "B", "C"), align = "h",
                                   ncol = 3, nrow = 1,
                                   common.legend = TRUE, legend = "top")
combinedSpatPredBox
#ggsave("combinedSpatPredBox.png", width = 39, height = 13, units = "cm")

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
