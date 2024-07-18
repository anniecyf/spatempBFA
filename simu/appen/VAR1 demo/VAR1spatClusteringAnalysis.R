############################## accuracy ratio and Rand index box & violin plots ##############################
rm(list=ls())
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/VAR1 demo", sep = "/"))
library(ggplot2)
library(ggpubr)
library(dplyr)
legendColors <- c("#FFAA11","#55CC11", "#1177CC", "#7755CC")
load("VAR1nkeep10AccuRatioMat_designedClustering.RData") #nkeep10AccuRatio
load("VAR1nkeep100AccuRatioMat_designedClustering.RData") #nkeep100AccuRatio
load("VAR1nkeep1000AccuRatioMat_designedClustering.RData") #nkeep1000AccuRatio
load("VAR1nkeep10RandIndexMat_designedClustering.RData") #nkeep10RandIndex
load("VAR1nkeep100RandIndexMat_designedClustering.RData") #nkeep100RandIndex
load("VAR1nkeep1000RandIndexMat_designedClustering.RData") #nkeep1000RandIndex
N <- nrow(nkeep10AccuRatio)
all(nkeep1000AccuRatio[,1:3] == 1) # TRUE
all(nkeep100AccuRatio[,1:3] == 1) # TRUE
all(nkeep10AccuRatio[,1:3] == 1) # TRUE
all(nkeep1000RandIndex[,1:3] == 1) # TRUE
all(nkeep100RandIndex[,1:3] == 1) # TRUE
all(nkeep10RandIndex[,1:3] == 1) # TRUE
NNGPsequenVaryLjDF <- data.frame(nkeep1000accu = nkeep1000AccuRatio[,4], 
                                 nkeep100accu = nkeep100AccuRatio[,4], 
                                 nkeep10accu = nkeep10AccuRatio[,4],
                                 nkeep1000rand = nkeep1000RandIndex[,4], 
                                 nkeep100rand = nkeep100RandIndex[,4], 
                                 nkeep10rand = nkeep10RandIndex[,4])
NNGPsequenVaryLjHistAccu1000 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep1000accu), bins = 25, origin = 0) +
  labs(y = "Count", x = "Accuracy Ratio") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
NNGPsequenVaryLjHistAccu100 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep100accu), bins = 25, origin = 0) +
  labs(y = "Count", x = "Accuracy Ratio") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
NNGPsequenVaryLjHistAccu10 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep10accu), bins = 25, origin = 0) +
  labs(y = "Count", x = "Accuracy Ratio") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
NNGPsequenVaryLjHistRand1000 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep1000rand), bins = 25, origin = 0) +
  labs(y = "Count", x = "Rand Index") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
NNGPsequenVaryLjHistRand100 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep100rand), bins = 25, origin = 0) +
  labs(y = "Count", x = "Rand Index") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())   
NNGPsequenVaryLjHistRand10 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep10rand), bins = 25, origin = 0) +
  labs(y = "Count", x = "Rand Index") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
ggarrange(NNGPsequenVaryLjHistAccu1000, NNGPsequenVaryLjHistAccu100, NNGPsequenVaryLjHistAccu10,
          NNGPsequenVaryLjHistRand1000, NNGPsequenVaryLjHistRand100, NNGPsequenVaryLjHistRand10,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("VAR1accuRatioRandIndNNGPsequenVaryLjHist101001000iterWeights.png", width=24, height=16, units="cm")
modelVec <- rep(c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj"), each = N)
df10accu <- data.frame(accuratio = as.vector(nkeep10AccuRatio), model = as.factor(modelVec))
df100accu <- data.frame(accuratio = as.vector(nkeep100AccuRatio), model = as.factor(modelVec))
df1000accu <- data.frame(accuratio = as.vector(nkeep1000AccuRatio), model = as.factor(modelVec))
df10randInd <- data.frame(randindex = as.vector(nkeep10RandIndex), model = as.factor(modelVec))
df100randInd <- data.frame(randindex = as.vector(nkeep100RandIndex), model = as.factor(modelVec))
df1000randInd <- data.frame(randindex = as.vector(nkeep1000RandIndex), model = as.factor(modelVec))
p1000accuBoxLegend <- ggplot(df1000accu) + 
  geom_boxplot(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000accuBoxLegend <- p1000accuBoxLegend + labs(y="Accuracy Ratio", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p1000accuViolinLegend <- ggplot(df1000accu) + 
  geom_violin(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000accuViolinLegend <- p1000accuViolinLegend + labs(y="Accuracy Ratio", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p100accuBoxLegend <- ggplot(df100accu) + 
  geom_boxplot(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100accuBoxLegend <- p100accuBoxLegend + labs(y="Accuracy Ratio", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p100accuViolinLegend <- ggplot(df100accu) + 
  geom_violin(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100accuViolinLegend <- p100accuViolinLegend + labs(y="Accuracy Ratio", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p10accuBoxLegend <- ggplot(df10accu) + 
  geom_boxplot(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10accuBoxLegend <- p10accuBoxLegend + labs(y="Accuracy Ratio", x = "") +
  theme(axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p10accuViolinLegend <- ggplot(df10accu) + 
  geom_violin(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10accuViolinLegend <- p10accuViolinLegend + labs(y="Accuracy Ratio", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p1000randIndBoxLegend <- ggplot(df1000randInd) + 
  geom_boxplot(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000randIndBoxLegend <- p1000randIndBoxLegend + labs(y="Rand Index", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p1000randIndViolinLegend <- ggplot(df1000randInd) + 
  geom_violin(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000randIndViolinLegend <- p1000randIndViolinLegend + labs(y="Rand Index", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p100randIndBoxLegend <- ggplot(df100randInd) + 
  geom_boxplot(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100randIndBoxLegend <- p100randIndBoxLegend + labs(y="Rand Index", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p100randIndViolinLegend <- ggplot(df100randInd) + 
  geom_violin(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100randIndViolinLegend <- p100randIndViolinLegend + labs(y="Rand Index", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p10randIndBoxLegend <- ggplot(df10randInd) + 
  geom_boxplot(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10randIndBoxLegend <- p10randIndBoxLegend + labs(y="Rand Index", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
p10randIndViolinLegend <- ggplot(df10randInd) + 
  geom_violin(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10randIndViolinLegend <- p10randIndViolinLegend + labs(y="Rand Index", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.key = element_blank())
ggarrange(p1000accuBoxLegend, p100accuBoxLegend, p10accuBoxLegend,
          p1000accuViolinLegend, p100accuViolinLegend, p10accuViolinLegend,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("VAR1accuRatioBoxViolin101001000iterWeights.png", width=24, height=16, units="cm")
ggarrange(p1000randIndBoxLegend, p100randIndBoxLegend, p10randIndBoxLegend,
          p1000randIndViolinLegend, p100randIndViolinLegend, p10randIndViolinLegend,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("VAR1randIndBoxViolin101001000iterWeights.png", width=24, height=16, units="cm")
ggarrange(p1000accuBoxLegend, p100accuBoxLegend, p10accuBoxLegend,
          p1000randIndBoxLegend, p100randIndBoxLegend, p1000randIndBoxLegend,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("VAR1accuRatioRandIndBox101001000iterWeights.png", width=24, height=16, units="cm")
ggarrange(p1000accuViolinLegend, p100accuViolinLegend, p10accuViolinLegend,
          p1000randIndViolinLegend, p100randIndViolinLegend, p10randIndViolinLegend,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("VAR1accuRatioRandIndViolin101001000iterWeights.png", width=24, height=16, units="cm")


############################## specific spatial clustering example results verification ##############################
rm(list=ls())
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/VAR1 demo", sep = "/"))
load("VAR1fittedClusGpMat_exDesignedClusteringSimu.RData")
VAR1fittedClusGpMat <- fittedClusGpMat
setwd(paste(projDirec, "simu/spatClustering", sep = "/"))
load("fittedClusGpMat_exDesignedClusteringSimu.RData")
all((fittedClusGpMat[,c(4,7,10)] + VAR1fittedClusGpMat[,c(4,7,10)]) == 3) # TRUE
all(fittedClusGpMat[,-c(4,7,10)] == VAR1fittedClusGpMat[,-c(4,7,10)]) # TRUE
