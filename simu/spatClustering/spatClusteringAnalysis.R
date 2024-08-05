############################## accuracy ratio and Rand index box & violin plots ##############################
rm(list=ls())
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/spatClustering", sep = "/"))
library(ggplot2)
library(ggpubr)
library(dplyr)
legendColors <- c("#FFAA11","#55CC11", "#1177CC", "#7755CC")
load("nkeep10AccuRatioMat_designedClustering.RData") #nkeep10AccuRatio
load("nkeep100AccuRatioMat_designedClustering.RData") #nkeep100AccuRatio
load("nkeep1000AccuRatioMat_designedClustering.RData") #nkeep1000AccuRatio
load("nkeep10RandIndexMat_designedClustering.RData") #nkeep10RandIndex
load("nkeep100RandIndexMat_designedClustering.RData") #nkeep100RandIndex
load("nkeep1000RandIndexMat_designedClustering.RData") #nkeep1000RandIndex
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
# nkeep1000accu  nkeep100accu   nkeep10accu nkeep1000rand  nkeep100rand   nkeep10rand 
# 0.9781000     0.9778000     0.9758000     0.9625717     0.9621333     0.9604242 
colMeans(NNGPsequenVaryLjDF)
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
ggsave("accuRatioRandIndNNGPsequenVaryLjHist101001000iterWeights.png", width=24, height=16, units="cm")
# NNGPsequenVaryLjDF <- data.frame(accuratio = c(nkeep1000AccuRatio[,4], nkeep100AccuRatio[,4], nkeep10AccuRatio[,4]),
#                                  randindex = c(nkeep1000RandIndex[,4], nkeep100RandIndex[,4], nkeep10RandIndex[,4]),
#                                  nkeep = rep(c(1000, 100, 10), each = N))
# NNGPsequenVaryLjBoxAccu <- ggplot(NNGPsequenVaryLjDF) + geom_boxplot(aes(x = nkeep, y = accuratio, fill = nkeep)) +
#   labs(y="Accuracy Ratio", x = "") +
#   scale_fill_discrete("Clustering Kept Iteration Number") +
#   theme(axis.title.y = element_text(size = 12),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10, color = "black"),
#         axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
#         legend.title = element_text(face = "bold", size = 12), 
#         legend.key = element_blank())
# NNGPsequenVaryLjBoxAccu
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
ggsave("accuRatioBoxViolin101001000iterWeights.png", width=24, height=16, units="cm")
ggarrange(p1000randIndBoxLegend, p100randIndBoxLegend, p10randIndBoxLegend,
          p1000randIndViolinLegend, p100randIndViolinLegend, p10randIndViolinLegend,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("randIndBoxViolin101001000iterWeights.png", width=24, height=16, units="cm")
ggarrange(p1000accuBoxLegend, p100accuBoxLegend, p10accuBoxLegend,
          p1000randIndBoxLegend, p100randIndBoxLegend, p1000randIndBoxLegend,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("accuRatioRandIndBox101001000iterWeights.png", width=24, height=16, units="cm")
ggarrange(p1000accuViolinLegend, p100accuViolinLegend, p10accuViolinLegend,
          p1000randIndViolinLegend, p100randIndViolinLegend, p10randIndViolinLegend,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("accuRatioRandIndViolin101001000iterWeights.png", width=24, height=16, units="cm")


############################## specific spatial clustering example results verification and spatial plots ##############################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/spatClustering", sep = "/"))
load("fittedClusGpMat_exDesignedClusteringSimu.RData")
# fittedClusGpMatNew <- fittedClusGpMat
# setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/spatial generalization aspect/clusteringSimulation")
# load("fittedClusGpMat_exDesignedClusteringSimu.RData")
# fittedClusGpMat == fittedClusGpMatNew
sqrootM <- 10; M <- 100
calcGroup1 <- function(row, col, sqrootM){
  return((row - 1)*sqrootM + col)
}
rowColInd <- expand.grid(2:8, 2:8)
whichGroup1 <- mapply(calcGroup1, rowColInd[,2], rowColInd[,1], sqrootM = sqrootM)
spatGroupOverall <- rep(2,M)
spatGroupOverall[whichGroup1] <- 1
matrix(spatGroupOverall, 10, 10)
veriBool <- TRUE
for(j in 1:12){
  fittedClusGp <- fittedClusGpMat[,j]
  if(any(fittedClusGp != spatGroupOverall))
    if(any(fittedClusGp + spatGroupOverall != 3))
      veriBool = FALSE
}
veriBool #TRUE
xcoord <- rep(1:sqrootM, sqrootM)
ycoord <- rep(seq(sqrootM, 1, by = -1), each = sqrootM)
spatGpDF <- data.frame(x = xcoord, y = ycoord, actualGp <- as.factor(spatGroupOverall))
actualSpatPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = actualGp), size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
actualSpatPlot # all estimated spatial plots agree with the actual one
ggsave("exSpatClustering.png", width=10, height=10, units="cm")
actualSpatPlotBW <- ggplot(spatGpDF) + theme_void() + 
  geom_point(aes(x = x, y = y, shape = actualGp, col = actualGp), size = 7, show.legend = FALSE) +
  scale_color_manual(values = c("black", "grey")) 
actualSpatPlotBW # all estimated spatial plots agree with the actual one
ggsave("exSpatClusteringBW.png", width=10, height=10, units="cm")
