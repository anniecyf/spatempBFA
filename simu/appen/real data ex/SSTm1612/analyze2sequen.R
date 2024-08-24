rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex/SSTm1612", sep = "/"))
# load("realExNNGPsequenFixedL.RData"); regFixedL.simu.sequen$runtime # days
# rm(regFixedL.simu.sequen)
load("NNGPsequenFixedLDiags.RData")
load("NNGPsequenFixedLDeviance.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("realExNNGPtemppredFixedLsequen.RData")
load("realExNNGPspatpredFixedLsequen.RData")
# load("realExNNGPsequenVaryLj.RData"); regVaryLj.simu.sequen$runtime # 12.52 days
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
ggarrange(SSTtestinglocDFday110actual, SSTtestinglocDFday110predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtestinglocDFday110plots.png", width = 16, height = 8, units = "cm")
SSTtestinglocDFday230predicted.sequenFixedL <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230predicted.sequenFixedL
ggarrange(SSTtestinglocDFday230actual, SSTtestinglocDFday230predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtestinglocDFday230plots.png", width = 16, height = 8, units = "cm")
SSTtestinglocDFday350predicted.sequenFixedL <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350predicted.sequenFixedL
ggarrange(SSTtestinglocDFday350actual, SSTtestinglocDFday350predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtestinglocDFday350plots.png", width = 16, height = 8, units = "cm")
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
ggarrange(SSTtestinglocDFday110actual, SSTtestinglocDFday110predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenVaryLjSSTtestinglocDFday110plots.png", width = 16, height = 8, units = "cm")
SSTtestinglocDFday230predicted.sequenVaryLj <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST230, size = SST230)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday230predicted.sequenVaryLj
ggarrange(SSTtestinglocDFday230actual, SSTtestinglocDFday230predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenVaryLjSSTtestinglocDFday230plots.png", width = 16, height = 8, units = "cm")
SSTtestinglocDFday350predicted.sequenVaryLj <- ggplot(SSTtestinglocDFpredicted) +
  geom_point(aes(x = Longitude, y = Latitude, color = SST350, size = SST350)) +
  scale_color_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  scale_size_continuous(name = "") +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtestinglocDFday350predicted.sequenVaryLj
ggarrange(SSTtestinglocDFday350actual, SSTtestinglocDFday350predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenVaryLjSSTtestinglocDFday350plots.png", width = 16, height = 8, units = "cm")
spatPredMetricMat[1, 2] <- mean((ylocPredMean - YtestingSpat)^2) # 0.8744089
diffMat <- sweep(ylocPred, 1, YtestingSpat, "-")
spatPredMetricMat[2, 2] = mean(rowMeans(diffMat^2)) # 15.46964
spatPredMetricMat[3, 2] = mean(apply(ylocPred, 1, var)) # 14.60984
spatPredMetricMat
# NNGPsequenFixedL NNGPsequenVaryLj
# postMeanMSE         1.026568        0.8744089
# postMSE            26.039550       15.4696393
# postVar            25.038020       14.6098402
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
SSTtraininglocDFday352actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST352)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday352actual
SSTtraininglocDFday353actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST353)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday353actual
SSTtraininglocDFday354actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST354)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday354actual
SSTtraininglocDFday355actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355actual
SSTtraininglocDFday360actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360actual
SSTtraininglocDFday365actual <- ggplot(SSTtraininglocDFactual) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365actual
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
ggarrange(SSTtraininglocDFday351actual, SSTtraininglocDFday351predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday351plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday352predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST352)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday352predicted.sequenFixedL
ggarrange(SSTtraininglocDFday352actual, SSTtraininglocDFday352predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday352plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday353predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST353)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday353predicted.sequenFixedL
ggarrange(SSTtraininglocDFday353actual, SSTtraininglocDFday353predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday353plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday354predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST354)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday354predicted.sequenFixedL
ggarrange(SSTtraininglocDFday354actual, SSTtraininglocDFday354predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday354plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday355predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355predicted.sequenFixedL
ggarrange(SSTtraininglocDFday355actual, SSTtraininglocDFday355predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday355plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday360predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360predicted.sequenFixedL
ggarrange(SSTtraininglocDFday360actual, SSTtraininglocDFday360predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday360plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday365predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365predicted.sequenFixedL
ggarrange(SSTtraininglocDFday365actual, SSTtraininglocDFday365predicted.sequenFixedL, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
ggsave("NNGPsequenFixedLSSTtraininglocDFday365plots.png", width = 16, height = 8, units = "cm")
tempPredMetricMat[1, 1] <- mean((ytempPredMean - ytempPred)^2)
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
ggarrange(SSTtraininglocDFday351actual, SSTtraininglocDFday351predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday351plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday352predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST352)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday352predicted.sequenVaryLj
ggarrange(SSTtraininglocDFday352actual, SSTtraininglocDFday352predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday352plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday353predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST353)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday353predicted.sequenVaryLj
ggarrange(SSTtraininglocDFday353actual, SSTtraininglocDFday353predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday353plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday354predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST354)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday354predicted.sequenVaryLj
ggarrange(SSTtraininglocDFday354actual, SSTtraininglocDFday354predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday354plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday355predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355predicted.sequenVaryLj
ggarrange(SSTtraininglocDFday355actual, SSTtraininglocDFday355predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday355plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday360predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360predicted.sequenVaryLj
ggarrange(SSTtraininglocDFday360actual, SSTtraininglocDFday360predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday360plots.png", width = 16, height = 8, units = "cm")
SSTtraininglocDFday365predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365predicted.sequenVaryLj
ggarrange(SSTtraininglocDFday365actual, SSTtraininglocDFday365predicted.sequenVaryLj, 
          labels = c("Actual", "Predicted"), align = "h", ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenVaryLjSSTtraininglocDFday365plots.png", width = 16, height = 8, units = "cm")
tempPredMetricMat[1, 2] <- mean((ytempPredMean - ytempPred)^2) # 175.9535
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 2] = mean(rowMeans(diffMat^2)) # 282.5623
tempPredMetricMat[3, 2] = mean(apply(ytempPred, 1, var)) # 176.1296
tempPredMetricMat
# NNGPsequenFixedL NNGPsequenVaryLj
# postMeanMSE         178.5120         175.9535
# postMSE             282.0187         282.5623
# postVar             178.6907         176.1296
ggarrange(SSTtraininglocDFday351actual, SSTtraininglocDFday351predicted.sequenFixedL, 
          SSTtraininglocDFday351predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday351plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday352actual, SSTtraininglocDFday352predicted.sequenFixedL, 
          SSTtraininglocDFday352predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday352plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday353actual, SSTtraininglocDFday353predicted.sequenFixedL, 
          SSTtraininglocDFday353predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday353plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday354actual, SSTtraininglocDFday354predicted.sequenFixedL, 
          SSTtraininglocDFday354predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday354plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday355actual, SSTtraininglocDFday355predicted.sequenFixedL, 
          SSTtraininglocDFday355predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday355plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday360actual, SSTtraininglocDFday360predicted.sequenFixedL, 
          SSTtraininglocDFday360predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday360plots.png", width = 24, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday365actual, SSTtraininglocDFday365predicted.sequenFixedL, 
          SSTtraininglocDFday365predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.08, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday365plots.png", width = 24, height = 8, units = "cm")

