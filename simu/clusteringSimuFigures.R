############################## accuracy ratio and Rand index box & violin plots ##############################
rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/spatial generalization aspect/clusteringSimulation")
library(ggplot2)
library(ggpubr)
library(dplyr)
legendColors <- c("#FFAA11","#55CC11", "#1177CC", "#7755CC")
load("nkeep10AccuRatioMat_designedClustering.RData") #nkeep10AccuRatio
load("nkeep100AccuRatioMat_designedClustering.RData") #nkeep100AccuRatio
load("nkeep1000AccuRatioMat_designedClustering.RData") #nkeep1000AccuRatio
N <- nrow(nkeep10AccuRatio)
modelVec <- rep(c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj"), each = N)
df10accu <- data.frame(accuratio = as.vector(nkeep10AccuRatio), model = as.factor(modelVec))
df100accu <- data.frame(accuratio = as.vector(nkeep100AccuRatio), model = as.factor(modelVec))
df1000accu <- data.frame(accuratio = as.vector(nkeep1000AccuRatio), model = as.factor(modelVec))

p1000accuBoxLegend <- ggplot(df1000accu) + 
  geom_boxplot(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000accuBoxLegend <- p1000accuBoxLegend + labs(y="Accuracy Ratio", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p1000accuViolinLegend <- ggplot(df1000accu) + 
  geom_violin(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000accuViolinLegend <- p1000accuViolinLegend + labs(y="Accuracy Ratio", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p100accuBoxLegend <- ggplot(df100accu) + 
  geom_boxplot(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100accuBoxLegend <- p100accuBoxLegend + labs(y="Accuracy Ratio", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p100accuViolinLegend <- ggplot(df100accu) + 
  geom_violin(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100accuViolinLegend <- p100accuViolinLegend + labs(y="Accuracy Ratio", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p10accuBoxLegend <- ggplot(df10accu) + 
  geom_boxplot(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10accuBoxLegend <- p10accuBoxLegend + labs(y="Accuracy Ratio", x = "") +
  theme(axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p10accuViolinLegend <- ggplot(df10accu) + 
  geom_violin(aes(x = model, y = accuratio, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10accuViolinLegend <- p10accuViolinLegend + labs(y="Accuracy Ratio", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))

ggarrange(p1000accuBoxLegend, p100accuBoxLegend, p10accuBoxLegend,
          p1000accuViolinLegend, p100accuViolinLegend, p10accuViolinLegend,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("clusteringPlots/boxViolinPlots/valignedBoldBoxViolin101001000iterWeights.png", width=24, height=16, units="cm")


rm(list=ls())
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/spatial generalization aspect/clusteringSimulation")
library(ggplot2)
library(ggpubr)
library(dplyr)
legendColors <- c("#FFAA11","#55CC11", "#1177CC", "#7755CC")
load("nkeep10RandIndexMat_designedClustering.RData") #nkeep10RandIndex
load("nkeep100RandIndexMat_designedClustering.RData") #nkeep100RandIndex
load("nkeep1000RandIndexMat_designedClustering.RData") #nkeep1000RandIndex
N <- nrow(nkeep10RandIndex)
modelVec <- rep(c("fullGPfixedL", "NNGPblockFixedL", "NNGPsequenFixedL", "NNGPsequenVaryLj"), each = N)
df10randInd <- data.frame(randindex = as.vector(nkeep10RandIndex), model = as.factor(modelVec))
df100randInd <- data.frame(randindex = as.vector(nkeep100RandIndex), model = as.factor(modelVec))
df1000randInd <- data.frame(randindex = as.vector(nkeep1000RandIndex), model = as.factor(modelVec))

p1000randIndBoxLegend <- ggplot(df1000randInd) + 
  geom_boxplot(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000randIndBoxLegend <- p1000randIndBoxLegend + labs(y="Rand Index", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p1000randIndViolinLegend <- ggplot(df1000randInd) + 
  geom_violin(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p1000randIndViolinLegend <- p1000randIndViolinLegend + labs(y="Rand Index", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p100randIndBoxLegend <- ggplot(df100randInd) + 
  geom_boxplot(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100randIndBoxLegend <- p100randIndBoxLegend + labs(y="Rand Index", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p100randIndViolinLegend <- ggplot(df100randInd) + 
  geom_violin(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p100randIndViolinLegend <- p100randIndViolinLegend + labs(y="Rand Index", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p10randIndBoxLegend <- ggplot(df10randInd) + 
  geom_boxplot(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10randIndBoxLegend <- p10randIndBoxLegend + labs(y="Rand Index", x = "") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))
p10randIndViolinLegend <- ggplot(df10randInd) + 
  geom_violin(aes(x = model, y = randindex, fill = model)) +
  scale_fill_manual("Model", values = legendColors)
p10randIndViolinLegend <- p10randIndViolinLegend + labs(y="Rand Index", x = "Model") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.title = element_text(face = "bold", size = 12))

ggarrange(p1000randIndBoxLegend, p100randIndBoxLegend, p10randIndBoxLegend,
          p1000randIndViolinLegend, p100randIndViolinLegend, p10randIndViolinLegend,
          labels = c("A", "C", "E", "B", "D", "F"), align = "v",
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("clusteringPlots/boxViolinPlots/randIndBoxViolin101001000iterWeights.png", width=24, height=16, units="cm")




############################## specific example spatial plots ##############################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/spatial generalization aspect/clusteringSimulation")
load("fittedClusGpMat_exDesignedClusteringSimu.RData")
sqrootM <- 10; M <- 100
calcGroup1 <- function(row, col, sqrootM){
  return((row - 1)*sqrootM + col)
}
rowColInd <- expand.grid(2:8, 2:8)
whichGroup1 <- mapply(calcGroup1, rowColInd[,2], rowColInd[,1], sqrootM = sqrootM)
spatGroupOverall <- rep(2,M)
spatGroupOverall[whichGroup1] <- 1
matrix(spatGroupOverall, 10, 10)
xcoord <- rep(1:sqrootM, sqrootM)
ycoord <- rep(seq(sqrootM, 1, by = -1), each = sqrootM)
spatGpDF <- data.frame(x = xcoord, y = ycoord, 
                         actualGp <- as.factor(spatGroupOverall), 
                         NNGPsequenVaryLj1000Gp = as.factor(fittedClusGpMat[,12]))
actualSpatPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = actualGp), size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC")) 
NNGPsequenVaryLjSpatPlot <- ggplot(spatGpDF) + geom_point(aes(x = x, y = y, col = NNGPsequenVaryLj1000Gp), 
                                                          size = 7, show.legend = FALSE) + 
  theme_void() + scale_color_manual(values = c("#55CC11", "#1177CC"))
agreements <- spatGroupOverall == fittedClusGpMat[,12] # pre-checked that the factor levels and labels do agree 
disaDF <- data.frame(x = xcoord[!agreements], y = ycoord[!agreements])
NNGPsequenVaryLjSpatPlot <- NNGPsequenVaryLjSpatPlot + 
  geom_point(data = disaDF, aes(x = x, y = y), col = "red", fill = NA, shape = 0, size = 8)

ggarrange(ggplot() + theme_void(), actualSpatPlot, 
          ggplot() + theme_void(),  NNGPsequenVaryLjSpatPlot, 
          nrow = 1, ncol = 4, widths = c(0.1, 1, 0.1, 1), align = "h", 
          labels = c("A", "", "B", ""), hjust = - 0.7, 
          font.label = list(size = 18, face = "bold", color = "black"))
ggsave("clusteringPlots/spatial plots/ggarrangeCombined18label.png", width=22, height=10, units="cm")









############################## Section 6.1 posterior deviance traceplots ##############################
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/spatial generalization aspect/simulation experiments/ubuntu/afterRhoPsiCorrect/varyLj")
setwd("iter50000")
load("regVaryLjsimuSequenT30M900Iter50000DevianceNostorealphaweight.RData")
Deviance.sequenVaryLj.iter50000 <- as.vector(Deviance.sequenVaryLj)
setwd("../iter30000")
load("regVaryLjsimuSequenT30M900Iter30000DevianceNostorealphaweight.RData")
Deviance.sequenVaryLj.iter30000 <- as.vector(Deviance.sequenVaryLj)
rm(Deviance.sequenVaryLj)
setwd("../../fixedL/sequen/iter50000")
load("regFixedL60simuSequenT30M900Iter50000DevianceNostorealphaweights.RData")
Deviance.sequen.iter50000 <- as.vector(Deviance.sequen)
setwd("../iter30000")
load("regFixedL60simuSequenT30M900Iter30000DevianceNostorealphaweights.RData")
Deviance.sequen.iter30000 <- as.vector(Deviance.sequen)
rm(Deviance.sequen)
setwd("../../block")
load("regFixedL60simuBlockT30M900Iter30000Deviance.RData")
Deviance.block.iter30000 <- as.vector(Deviance.block)
rm(Deviance.block)
setwd("../fullGP/iter30000")
load("regFixedL60simuT30M900Iter30000Deviance.RData")
Deviance.iter30000 <- as.vector(Deviance)
rm(Deviance)
devianceDF <- data.frame(iter = 1:5000, fullGP30000 = Deviance.iter30000, 
                         NNGPblock30000 = Deviance.block.iter30000,
                         NNGPsequen30000 = Deviance.sequen.iter30000,
                         NNGPsequen50000 = Deviance.sequen.iter50000,
                         NNGPsequenVaryLj30000 = Deviance.sequenVaryLj.iter30000,
                         NNGPsequenVaryLj50000 = Deviance.sequenVaryLj.iter50000)
devPlotfullGP30000 <- ggplot(devianceDF) + geom_line(aes(x = iter, y = fullGP30000), 
                                                     col = rgb(0.27, 0.35, 0.27)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill=NA, size=1))
devPlotfullGP30000
devPlotblock30000 <- ggplot(devianceDF) + geom_line(aes(x = iter, y = NNGPblock30000), 
                                                    col = rgb(0.27, 0.35, 0.27)) 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill=NA, size=1))
devPlotblock30000
devPlotsequen30000 <- ggplot(devianceDF) + geom_line(aes(x = iter, y = NNGPsequen30000), 
                                                     col = rgb(0.27, 0.35, 0.27)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill=NA, size=1))
devPlotsequen30000
devPlotsequen50000 <- ggplot(devianceDF) + geom_line(aes(x = iter, y = NNGPsequen50000), 
                                                     col = rgb(0.27, 0.35, 0.27)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill=NA, size=1))
devPlotsequen50000
devPlotsequenVaryLj30000 <- ggplot(devianceDF) + geom_line(aes(x = iter, y = NNGPsequenVaryLj30000), 
                                                           col = rgb(0.27, 0.35, 0.27)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill=NA, size=1))
devPlotsequenVaryLj30000
devPlotsequenVaryLj50000 <- ggplot(devianceDF) + geom_line(aes(x = iter, y = NNGPsequenVaryLj50000), 
                                                           col = rgb(0.27, 0.35, 0.27)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill=NA, size=1))
devPlotsequenVaryLj50000
ggarrange(devPlotfullGP30000, devPlotblock30000, 
          devPlotsequen30000, devPlotsequen50000,
          devPlotsequenVaryLj30000, devPlotsequenVaryLj50000,
          labels = c("A", "B", "C", "D", "E", "F"), 
          hjust = -0.1, vjust = 1.2,
          align = "v", ncol = 2, nrow = 3)
ggsave("postDevianceBWTheme.pdf", width = 18, height = 24, units = "cm")





