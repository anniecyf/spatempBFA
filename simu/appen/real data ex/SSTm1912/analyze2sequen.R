rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex/SSTm1912", sep = "/"))
# load("realExNNGPsequenFixedL.RData"); regFixedL.simu.sequen$runtime # 
# rm(regFixedL.simu.sequen)
load("NNGPsequenFixedLDiags.RData")
load("NNGPsequenFixedLDeviance.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("realExNNGPtemppredFixedLsequen.RData")
load("realExNNGPspatpredFixedLsequen.RData")
# load("realExNNGPsequenVaryLj.RData"); regVaryLj.simu.sequen$runtime # 
# rm(regVaryLj.simu.sequen)
load("NNGPsequenVaryLjDiags.RData")
load("NNGPsequenVaryLjDeviance.RData")
load("GibbsStepTimeVaryLjSequen.RData")
load("realExNNGPtemppredVaryLj.RData")
load("realExNNGPspatpredVaryLj.RData")
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)
NKeep <- 1000
postDeviancesDF <- data.frame(MCMCiter = 1:NKeep,
                              NNGPsequenFixedL = as.vector(Deviance.sequen), 
                              NNGPsequenVaryLj = as.vector(Deviance.sequenVaryLj))
NNGPsequenFixedLpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPsequenFixedL)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenFixedLpostDeviances
NNGPsequenVaryLjpostDeviances <- ggplot(postDeviancesDF) + geom_line(aes(x = MCMCiter, y = NNGPsequenVaryLj)) + 
  labs(x = "MCMC Iteration", y = "Posterior Deviance") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
NNGPsequenVaryLjpostDeviances
ggarrange(NNGPsequenFixedLpostDeviances, NNGPsequenVaryLjpostDeviances,
          labels = c("A", "B"), align = "h",
          ncol = 2, nrow = 1)
#ggsave("realExPostDeviances.png", width = 16, height = 8, units = "cm")
# View(Diags.sequen)
# View(Diags.sequenVaryLj)
M <- 1912
m = M - 50 # 50 testing spatial points
Nu <- 365 
trainingT <- 350 # training set T = 300
testingT <- 15
load("SSTtrainingDF2023.RData") # SSTtrainingDF
set.seed(29)
testingLocs <- sample(1:nrow(SSTtrainingDF), size = 35, replace = FALSE)
YtestingSpat <- as.vector(t(SSTtrainingDF[testingLocs, 3:(trainingT+2)]))
YtestingTemp <- as.vector(SSTtrainingDF[-testingLocs, (trainingT+3):(Nu+2)])

summary(YtestingSpat)
summary(abs(YtestingSpat))
spatPredMetricMat <- matrix(0, 3, 2, 
                            dimnames = list(c("postMeanMSE", "postMSE", "postVar"),
                                            c("NNGPsequenFixedL", "NNGPsequenVaryLj")))
testingLocNum <- M - m
ylocPredList = spatpredFixedLsequen$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
spatPredMetricMat[1, 3] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 3] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 3] = mean(apply(ylocPred, 1, var))
ylocPredList = spatpredVaryLj$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
spatPredMetricMat[1, 4] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 4] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 4] = mean(apply(ylocPred, 1, var))
spatPredMetricMat
spatPredKrigTime <- matrix(0, 2, 2, 
                            dimnames = list(c("alpha", "weightsXiLambda"),
                                            c("NNGPsequenFixedL", "NNGPsequenVaryLj")))
spatPredKrigTime[1, 1] = spatpredFixedLsequen$alphaKrigTime
spatPredKrigTime[2, 1] = spatpredFixedLsequen$weightsXiLambdaKrigTime
spatPredKrigTime[1, 2] = spatpredVaryLj$alphaKrigTime
spatPredKrigTime[2, 2] = spatpredVaryLj$weightsXiLambdaKrigTime
spatPredKrigTime

summary(YtestingTemp)
summary(abs(YtestingTemp))
tempPredMetricMat <- matrix(0, 3, 2, 
                            dimnames = list(c("postMeanMSE", "postMSE", "postVar"),
                                            c("NNGPsequenFixedL", "NNGPsequenVaryLj")))
ytempPredList = temppredFixedLsequen$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 1] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 1] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 1] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredVaryLj$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
tempPredMetricMat[1, 2] <- mean((ytempPredMean - ytempPred)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 2] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 2] = mean(apply(ytempPred, 1, var))
tempPredMetricMat