rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/spatial generalization aspect/newLoc prediction simulation/complete version/N20newestPackage")
load("s17N2M36h15T20spatpredTimeMat_randomNewLocPred.RData")
s17spatpredTimeMat <- spatpredTimeMat
load("s19N20M36h15T20spatpredTimeMat_randomNewLocPred.RData")
s19spatpredTimeMat <- spatpredTimeMat
load("s21N2M36h15T20spatpredTimeMat_randomNewLocPred.RData")
s21spatpredTimeMat <- spatpredTimeMat
load("s23N20M36h15T20spatpredTimeMat_randomNewLocPred.RData")
s23spatpredTimeMat <- spatpredTimeMat
load("s31N4M36h15T20spatpredTimeMat_randomNewLocPred.RData")
s31spatpredTimeMat <- spatpredTimeMat
rm(spatpredTimeMat)
m1timeMat <- rbind(s17spatpredTimeMat[3,,], s19spatpredTimeMat[3,,], s21spatpredTimeMat[3,,], 
                   s23spatpredTimeMat[3,,], s31spatpredTimeMat[3,,])
model1timeVec <- as.vector(m1timeMat)
summary(model1timeVec)
m2timeMat <- rbind(s17spatpredTimeMat[4,,], s19spatpredTimeMat[4,,], s21spatpredTimeMat[4,,], 
                   s23spatpredTimeMat[4,,], s31spatpredTimeMat[4,,])
model2timeVec <- as.vector(m2timeMat)
summary(model2timeVec)
m3timeMat <- rbind(s17spatpredTimeMat[5,,], s19spatpredTimeMat[5,,], s21spatpredTimeMat[5,,], 
                   s23spatpredTimeMat[5,,], s31spatpredTimeMat[5,,])
model3timeVec <- as.vector(m3timeMat)
summary(model3timeVec)
m4timeMat <- rbind(s17spatpredTimeMat[6,,], s19spatpredTimeMat[6,,], s21spatpredTimeMat[6,,], 
                   s23spatpredTimeMat[6,,], s31spatpredTimeMat[6,,])
model4timeVec <- as.vector(m4timeMat)
summary(model4timeVec)
# original version a total of 6 models (with 2 additional models "noCLfullGP" and "noCLNNGP")
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA/simu/spatPred")
library(tidyverse)
library(ggpubr)
N <- length(model1timeVec) # N_1 x m (48*36)
modelVec <- rep(c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj"), each = N)
dfspatpredtime <- data.frame(spatpredtime = c(model1timeVec, model2timeVec, model3timeVec, model4timeVec), 
                             model = as.factor(modelVec))
legendColors <- c("#FFAA11","#55CC11", "#1177CC", "#7755CC")
# spatpredtimeBoxLegendHori <- ggplot(dfspatpredtime) + 
#   geom_boxplot(aes(x = spatpredtime, y = model, fill = model)) +
#   scale_fill_manual("Model", values = legendColors)
# spatpredtimeBoxLegendHori <- spatpredtimeBoxLegendHori + labs(y="", x = "Spatial Prediction Time (sec)") +
#   theme(axis.title.x = element_text(size = 12), #, face = "bold"
#         axis.text.x = element_text(size = 10, color = "black"),
#         axis.text.y = element_blank(),
#         axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
#         legend.position = "bottom")
# spatpredtimeBoxLegendHori
# spatpredtimeViolinLegendHori <- ggplot(dfspatpredtime) + 
#   geom_violin(aes(x = spatpredtime, y = model, fill = model)) +
#   scale_fill_manual("Model", values = legendColors)
# spatpredtimeViolinLegendHori <- spatpredtimeViolinLegendHori + labs(y="", x = "Spatial Prediction Time (sec)") +
#   theme(axis.title.x = element_text(size = 12), #, face = "bold"
#         axis.text.x = element_text(size = 10, color = "black"),
#         axis.text.y = element_blank(),
#         axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
#         legend.position = "bottom")
# spatpredtimeViolinLegendHori
# ggarrange(spatpredtimeBoxLegendHori, spatpredtimeViolinLegendHori, labels = c("A", "B"),
#           ncol = 2, nrow = 1, align = "h", common.legend = TRUE, legend = "bottom")
spatpredtimeBoxLegendHoriNoAxes <- ggplot(dfspatpredtime) + 
  geom_boxplot(aes(x = spatpredtime, y = model, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
spatpredtimeBoxLegendHoriNoAxes <- spatpredtimeBoxLegendHoriNoAxes + labs(y="", x = "") +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.title = element_text(face = "bold", size = 12))
spatpredtimeBoxLegendHoriNoAxes
spatpredtimeViolinLegendHoriNoAxes <- ggplot(dfspatpredtime) + 
  geom_violin(aes(x = spatpredtime, y = model, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
spatpredtimeViolinLegendHoriNoAxes <- spatpredtimeViolinLegendHoriNoAxes + labs(y="", x = "") +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.title = element_text(face = "bold", size = 12))
spatpredtimeViolinLegendHoriNoAxes
combinedHoriNoAxes <- ggarrange(spatpredtimeBoxLegendHoriNoAxes, spatpredtimeViolinLegendHoriNoAxes, 
                                ncol = 2, nrow = 1, align = "h", common.legend = TRUE, legend = "top")
combinedHoriNoAxes
# annotate_figure(combinedHoriNoAxes,
#                 bottom = text_grob("Spatial Prediction Time (sec)", color = "black", y = 1.2,
#                                    face = "bold", size = 12))
# ggsave("spatpredtimeBoxViolinCombinedHoriNoAxes.png", width=16, height=8, units="cm")
annotate_figure(combinedHoriNoAxes,
                bottom = text_grob("Spatial Prediction Time (sec)", color = "black", y = 1.2,
                                   size = 12))
ggsave("spatpredtimeBoxViolinCombinedHoriNoAxesNoBold.png", width=16, height=8, units="cm")
spatpredtimeBox <- ggplot(dfspatpredtime) + 
  geom_boxplot(aes(x = spatpredtime, y = model)) + labs(y="", x = "") +
  theme(axis.text.x = element_text(size = 10, color = "black"), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
spatpredtimeBox
spatpredtimeViolin <- ggplot(dfspatpredtime) + 
  geom_violin(aes(x = spatpredtime, y = model)) + labs(y="", x = "") +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
spatpredtimeViolin
combinedHori <- ggarrange(spatpredtimeBox, spatpredtimeViolin, 
                                ncol = 2, nrow = 1, widths = c(3, 2.15), align = "h")
combinedHori
annotate_figure(combinedHori,
                bottom = text_grob("Spatial Prediction Time (sec)", color = "black", y = 1.2,
                                   size = 12, hjust = 0.25))
ggsave("spatpredtimeBoxViolinCombinedHori.png", width=16, height=8, units="cm")
