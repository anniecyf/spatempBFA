rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex/SSTm1612", sep = "/"))
# load("realExNNGPsequenFixedL.RData"); regFixedL.simu.sequen$runtime # 17.17 days
# rm(regFixedL.simu.sequen)
load("NNGPsequenFixedLDiags.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("realExNNGPtemppredFixedLsequen.RData")
load("realExNNGPspatpredFixedLsequen.RData")
# load("realExNNGPsequenVaryLj.RData"); regVaryLj.simu.sequen$runtime # 12.52 days
# rm(regVaryLj.simu.sequen)
load("NNGPsequenVaryLjDiags.RData")
load("GibbsStepTimeVaryLjSequen.RData")
load("realExNNGPtemppredVaryLj.RData")
load("realExNNGPspatpredVaryLj.RData")
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)
# View(Diags.sequen)
# View(Diags.sequenVaryLj)
M <- 1612
m = M - 35 # 35 testing spatial points
Nu <- 365 
O <- 1
trainingT <- 350 # training set T = 300
testingT <- 15
load("SSTtrainingDF2023.RData") # SSTtrainingDF
set.seed(29)
testingLocs <- sample(1:nrow(SSTtrainingDF), size = 35, replace = FALSE)
YtestingSpat <- as.vector(t(SSTtrainingDF[testingLocs, 3:(trainingT+2)]))
YtestingTemp <- as.vector(as.matrix(SSTtrainingDF[-testingLocs, (trainingT+3):(Nu+2)]))
SSTtestinglocDF <- SSTtrainingDF[testingLocs, ]
SSTtestinglocDFactual <- SSTtestinglocDF %>% select(Longitude, Latitude, Day110, Day230, Day350) %>% 
  rename(SST110 = Day110, SST230 = Day230, SST350 = Day350)
SSTtestinglocDFday110actual <- ggplot(SSTtestinglocDFactual) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST110, size = SST110)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday110actual
SSTtestinglocDFday230actual <- ggplot(SSTtestinglocDFactual) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230actual
SSTtestinglocDFday350actual <- ggplot(SSTtestinglocDFactual) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350actual
summary(YtestingSpat)
summary(abs(YtestingSpat))
spatPredMetricMat <- matrix(0, 3, 2, 
                            dimnames = list(c("postMeanMSE", "postMSE", "postVar"),
                                            c("NNGPsequenFixedL", "NNGPsequenVaryLj")))
testingLocNum <- M - m
ylocPredList = spatpredFixedLsequen$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
ylocPredMeanMat <- t(matrix(ylocPredMean, nrow = trainingT, ncol = testingLocNum)) # testingLocNum * trainingT
SSTtestinglocDFpredicted <- cbind(SSTtestinglocDFactual[,1:2], ylocPredMeanMat[,c(110, 230, 350)])
colnames(SSTtestinglocDFpredicted) <- colnames(SSTtestinglocDFactual)
SSTtestinglocDFday110predicted.sequenFixedL <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST110, size = SST110)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday110predicted.sequenFixedL
SSTtestinglocDFday230predicted.sequenFixedL <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230predicted.sequenFixedL
SSTtestinglocDFday350predicted.sequenFixedL <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350predicted.sequenFixedL
spatPredMetricMat[1, 1] <- mean((ylocPredMean - YtestingSpat)^2)
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 1] = mean(rowMeans(diffMat^2))
spatPredMetricMat[3, 1] = mean(apply(ylocPred, 1, var))
ylocPredList = spatpredVaryLj$Y
ylocPred <- t(matrix(unlist(ylocPredList), ncol = testingLocNum * trainingT * O)) #(testingLocNum x trainingT x O) x NKeep 
ylocPredMean <- apply(ylocPred, 1, mean) # of length N = testingLocNum x O x trainingT (rowMeans) note that O = 1 here
ylocPredMeanMat <- t(matrix(ylocPredMean, nrow = trainingT, ncol = testingLocNum)) # testingLocNum * trainingT
SSTtestinglocDFpredicted <- cbind(SSTtestinglocDFactual[,1:2], ylocPredMeanMat[,c(110, 230, 350)])
colnames(SSTtestinglocDFpredicted) <- colnames(SSTtestinglocDFactual)
SSTtestinglocDFday110predicted.sequenVaryLj <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST110, size = SST110)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday110predicted.sequenVaryLj
SSTtestinglocDFday230predicted.sequenVaryLj <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230predicted.sequenVaryLj
SSTtestinglocDFday350predicted.sequenVaryLj <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350predicted.sequenVaryLj
spatPredMetricMat[1, 2] <- mean((ylocPredMean - YtestingSpat)^2) 
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 2] = mean(rowMeans(diffMat^2)) 
spatPredMetricMat[3, 2] = mean(apply(ylocPred, 1, var))
spatPredMetricMat
spatPredKrigTime <- matrix(0, 2, 2, 
                            dimnames = list(c("alpha", "weightsXiLambda"),
                                            c("NNGPsequenFixedL", "NNGPsequenVaryLj")))
spatPredKrigTime[1, 1] = spatpredFixedLsequen$alphaKrigTime
spatPredKrigTime[2, 1] = spatpredFixedLsequen$weightsXiLambdaKrigTime
spatPredKrigTime[1, 2] = spatpredVaryLj$alphaKrigTime
spatPredKrigTime[2, 2] = spatpredVaryLj$weightsXiLambdaKrigTime
spatPredKrigTime
ggarrange(SSTtestinglocDFday110actual, 
          SSTtestinglocDFday110predicted.sequenFixedL,
          SSTtestinglocDFday110predicted.sequenVaryLj,
          labels = c("Actual", "FixedL", "VaryLj"), align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtestinglocDFday110plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtestinglocDFday230actual, 
          SSTtestinglocDFday230predicted.sequenFixedL,
          SSTtestinglocDFday230predicted.sequenVaryLj,
          labels = c("Actual", "FixedL", "VaryLj"), align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtestinglocDFday230plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtestinglocDFday350actual, 
          SSTtestinglocDFday350predicted.sequenFixedL,
          SSTtestinglocDFday350predicted.sequenVaryLj,
          labels = c("Actual", "FixedL", "VaryLj"), align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtestinglocDFday350plots.png", width = 24, height = 8, units = "cm")

SSTtraininglocDF <- SSTtrainingDF[-testingLocs, ]
SSTtraininglocDFactual <- SSTtraininglocDF %>% 
  select(Longitude, Latitude, Day351, Day352, Day353, Day354, Day355, Day360, Day365) %>% 
  rename(SST351 = Day351, SST352 = Day352, SST353 = Day353, SST354 = Day354, 
         SST355 = Day355, SST360 = Day360, SST365 = Day365)
SSTtraininglocDFday351actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST351)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday351actual
summary(YtestingTemp)
summary(abs(YtestingTemp))
tempPredMetricMat <- matrix(0, 3, 2, 
                            dimnames = list(c("postMeanMSE", "postMSE", "postVar"),
                                            c("NNGPsequenFixedL", "NNGPsequenVaryLj")))
ytempPredList = temppredFixedLsequen$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
ytempPredMeanMat <- matrix(ytempPredMean, nrow = m, ncol = testingT) 
SSTtraininglocDFpredicted <- cbind(SSTtraininglocDFactual[,1:2], ytempPredMeanMat[,c(1, 2, 3, 4, 5, 10, 15)])
colnames(SSTtraininglocDFpredicted) <- colnames(SSTtraininglocDFactual)
SSTtraininglocDFday351predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST351)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday351predicted.sequenFixedL
tempPredMetricMat[1, 1] <- mean((ytempPredMean - YtestingTemp)^2)
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 1] = mean(rowMeans(diffMat^2))
tempPredMetricMat[3, 1] = mean(apply(ytempPred, 1, var))
ytempPredList = temppredVaryLj$Y
ytempPred <- t(matrix(unlist(ytempPredList), ncol = testingT * m * O)) # (testingT x m x O) x NKeep 
ytempPredMean <- apply(ytempPred, 1, mean) # of length N = testingT x O x m (rowMeans) note that O = 1 here
ytempPredMeanMat <- matrix(ytempPredMean, nrow = m, ncol = testingT) 
SSTtraininglocDFpredicted <- cbind(SSTtraininglocDFactual[,1:2], ytempPredMeanMat[,c(1, 2, 3, 4, 5, 10, 15)])
colnames(SSTtraininglocDFpredicted) <- colnames(SSTtraininglocDFactual)
SSTtraininglocDFday351predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST351)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday351predicted.sequenVaryLj
tempPredMetricMat[1, 2] <- mean((ytempPredMean - YtestingTemp)^2) 
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 2] = mean(rowMeans(diffMat^2)) 
tempPredMetricMat[3, 2] = mean(apply(ytempPred, 1, var)) 
tempPredMetricMat
ggarrange(SSTtraininglocDFday351actual, SSTtraininglocDFday351predicted.sequenFixedL, 
          SSTtraininglocDFday351predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday351plots.png", width = 24, height = 8, units = "cm")

