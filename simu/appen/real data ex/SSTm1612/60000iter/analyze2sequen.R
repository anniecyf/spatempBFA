rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex/SSTm1612/60000iter", sep = "/"))
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
M <- 1612
m = M - 35 # 35 testing spatial points
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

SSTtestinglocDF <- SSTtrainingDF[testingLocs, ]
SSTtestinglocDFactual <- SSTtestinglocDF %>% select(Longitude, Latitude, Day110, Day230, Day350) %>% 
  rename(SST110 = Day110, SST230 = Day230, SST350 = Day350)
ylocPredMeanMat <- t(matrix(ylocPredMean, nrow = trainingT, ncol = testingLocNum)) # testingLocNum * trainingT
SSTtestinglocDFpredicted <- cbind(SSTtestinglocDFactual[,1:2], ylocPredMeanMat[,c(110, 230, 350)])
colnames(SSTtestinglocDFpredicted) <- colnames(SSTtestinglocDFactual)
SSTtestinglocDFday110actual <- ggplot(SSTtestinglocDFactual) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST110, size = SST110)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday110actual
SSTtestinglocDFday110predicted <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST110, size = SST110)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday110predicted
ggarrange(SSTtestinglocDFday110actual, SSTtestinglocDFday110predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtestinglocDFday110plots.png", width = 16, height = 8, units = "cm")
SSTtestinglocDFday230actual <- ggplot(SSTtestinglocDFactual) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230actual
SSTtestinglocDFday230predicted <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230predicted
ggarrange(SSTtestinglocDFday230actual, SSTtestinglocDFday230predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtestinglocDFday230plots.png", width = 16, height = 8, units = "cm")
SSTtestinglocDFday350actual <- ggplot(SSTtestinglocDFactual) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350actual
SSTtestinglocDFday350predicted <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350predicted
ggarrange(SSTtestinglocDFday350actual, SSTtestinglocDFday350predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtestinglocDFday350plots.png", width = 16, height = 8, units = "cm")


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

SSTtraininglocDF <- SSTtrainingDF[-testingLocs, ]
SSTtraininglocDFactual <- SSTtraininglocDF %>% select(Longitude, Latitude, Day351, Day355, Day360, Day365) %>% 
  rename(SST351 = Day351, SST355 = Day355, SST360 = Day360, SST365 = Day365)
ytempPredMeanMat <- matrix(ytempPredMean, nrow = m, ncol = testingT) 
SSTtraininglocDFpredicted <- cbind(SSTtraininglocDFactual[,1:2], ytempPredMeanMat[,c(1, 5, 10, 15)])
colnames(SSTtraininglocDFpredicted) <- colnames(SSTtraininglocDFactual)
SSTtraininglocDFday351actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST351)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday351actual
SSTtraininglocDFday351predicted <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST351)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday351predicted
ggarrange(SSTtraininglocDFday351actual, SSTtraininglocDFday351predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtraininglocDFday351plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday355actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355actual
SSTtraininglocDFday355predicted <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355predicted
ggarrange(SSTtraininglocDFday355actual, SSTtraininglocDFday355predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtraininglocDFday355plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday360actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360actual
SSTtraininglocDFday360predicted <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360predicted
ggarrange(SSTtraininglocDFday360actual, SSTtraininglocDFday360predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtraininglocDFday360plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday365actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365actual
SSTtraininglocDFday365predicted <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365predicted
ggarrange(SSTtraininglocDFday365actual, SSTtraininglocDFday365predicted, 
          labels = c("A", "B"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("SSTtraininglocDFday365plots.png", width = 16, height = 8, units = "cm")