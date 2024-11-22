rm(list=ls())
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/VAR1 demo", sep = "/"))
library(tidyverse)
library(ggpubr)
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
colMeans(NNGPsequenVaryLjDF)
# nkeep1000accu  nkeep100accu   nkeep10accu nkeep1000rand  nkeep100rand   nkeep10rand 
# 0.9722222     0.9511111     0.9538889     0.9523008     0.9288215     0.9334119 
NNGPsequenVaryLjHistAccu1000 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep1000accu), bins = 25, boundary = 0) +
  labs(y = "Count", x = "Accuracy Ratio") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
NNGPsequenVaryLjHistAccu100 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep100accu), bins = 25, boundary = 0) +
  labs(y = "Count", x = "Accuracy Ratio") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
NNGPsequenVaryLjHistAccu10 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep10accu), bins = 25, boundary = 0) +
  labs(y = "Count", x = "Accuracy Ratio") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
NNGPsequenVaryLjHistRand1000 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep1000rand), bins = 25, boundary = 0) +
  labs(y = "Count", x = "Rand Index") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
NNGPsequenVaryLjHistRand100 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep100rand), bins = 25, boundary = 0) +
  labs(y = "Count", x = "Rand Index") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())   
NNGPsequenVaryLjHistRand10 <- ggplot(NNGPsequenVaryLjDF) + geom_histogram(aes(x = nkeep10rand), bins = 25, boundary = 0) +
  labs(y = "Count", x = "Rand Index") + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
ggarrange(NNGPsequenVaryLjHistAccu1000, NNGPsequenVaryLjHistAccu100, NNGPsequenVaryLjHistAccu10,
          NNGPsequenVaryLjHistRand1000, NNGPsequenVaryLjHistRand100, NNGPsequenVaryLjHistRand10,
          labels = c("A", "B", "C", "D", "E", "F"), align = "h",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
# ggsave("VAR1accuRatioRandIndNNGPsequenVaryLjHist101001000iterWeights.png", width=24, height=16, units="cm")

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
