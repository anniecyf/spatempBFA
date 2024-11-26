rm(list=ls())
library(tidyverse)
library(ggpubr)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex/SSTm1912", sep = "/"))
# regFixedL.simu.sequen$runtime # 34.07 days
load("NNGPsequenFixedLDiags.RData")
load("GibbsStepTimeFixedLsequen.RData")
load("realExNNGPtemppredFixedLsequen.RData")
load("realExNNGPspatpredFixedLsequen.RData")
# regVaryLj.simu.sequen$runtime # 24.34 days
load("NNGPsequenVaryLjDiags.RData")
load("GibbsStepTimeVaryLjSequen.RData")
load("realExNNGPtemppredVaryLj.RData")
load("realExNNGPspatpredVaryLj.RData")
apply(GibbsStepTimeFixedLsequen, 2, summary)
apply(GibbsStepTimeVaryLjSequen, 2, summary)
# View(Diags.sequen)
# View(Diags.sequenVaryLj)
M <- 1912
m = M - 50 # 50 testing spatial points
Nu <- 365 
O <- 1
trainingT <- 350 # training set T = 300
testingT <- 15
load("SSTtrainingDF2023.RData") # SSTtrainingDF
set.seed(29)
testingLocs <- sample(1:nrow(SSTtrainingDF), size = 50, replace = FALSE)
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
# NNGPsequenFixedL NNGPsequenVaryLj
# postMeanMSE          1.04397         1.026165
# postMSE             15.47524        12.078872
# postVar             14.44571        11.063771
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
#ggsave("NNGPsequenSSTtestinglocDFday110plots.png", width = 27, height = 10, units = "cm")
ggarrange(SSTtestinglocDFday230actual, 
          SSTtestinglocDFday230predicted.sequenFixedL,
          SSTtestinglocDFday230predicted.sequenVaryLj,
          labels = c("Actual", "FixedL", "VaryLj"), align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtestinglocDFday230plots.png", width = 27, height = 10, units = "cm")
ggarrange(SSTtestinglocDFday350actual, 
          SSTtestinglocDFday350predicted.sequenFixedL,
          SSTtestinglocDFday350predicted.sequenVaryLj,
          labels = c("Actual", "FixedL", "VaryLj"), align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtestinglocDFday350plots.png", width = 27, height = 10, units = "cm")

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
SSTtraininglocDFday352predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST352)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday352predicted.sequenFixedL
SSTtraininglocDFday353predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST353)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday353predicted.sequenFixedL
SSTtraininglocDFday354predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST354)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday354predicted.sequenFixedL
SSTtraininglocDFday355predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355predicted.sequenFixedL
SSTtraininglocDFday360predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360predicted.sequenFixedL
SSTtraininglocDFday365predicted.sequenFixedL <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365predicted.sequenFixedL
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
SSTtraininglocDFday352predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST352)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday352predicted.sequenVaryLj
SSTtraininglocDFday353predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST353)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday353predicted.sequenVaryLj
SSTtraininglocDFday354predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST354)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday354predicted.sequenVaryLj
SSTtraininglocDFday355predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST355)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday355predicted.sequenVaryLj
SSTtraininglocDFday360predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST360)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday360predicted.sequenVaryLj
SSTtraininglocDFday365predicted.sequenVaryLj <- ggplot(SSTtraininglocDFpredicted) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST365)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtraininglocDFday365predicted.sequenVaryLj
tempPredMetricMat[1, 2] <- mean((ytempPredMean - ytempPred)^2) 
diffMat <- sweep(ytempPred, 1, YtestingTemp, "-")
tempPredMetricMat[2, 2] = mean(rowMeans(diffMat^2)) 
tempPredMetricMat[3, 2] = mean(apply(ytempPred, 1, var)) 
tempPredMetricMat
# NNGPsequenFixedL NNGPsequenVaryLj
# postMeanMSE         188.6075         190.7238
# postMSE             299.1144         298.8963
# postVar             188.7963         190.9147
ggarrange(SSTtraininglocDFday351actual, SSTtraininglocDFday351predicted.sequenFixedL, 
          SSTtraininglocDFday351predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday351plots.png", width = 27, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday352actual, SSTtraininglocDFday352predicted.sequenFixedL, 
          SSTtraininglocDFday352predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday352plots.png", width = 27, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday353actual, SSTtraininglocDFday353predicted.sequenFixedL, 
          SSTtraininglocDFday353predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday353plots.png", width = 27, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday354actual, SSTtraininglocDFday354predicted.sequenFixedL, 
          SSTtraininglocDFday354predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday354plots.png", width = 27, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday355actual, SSTtraininglocDFday355predicted.sequenFixedL, 
          SSTtraininglocDFday355predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday355plots.png", width = 27, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday360actual, SSTtraininglocDFday360predicted.sequenFixedL, 
          SSTtraininglocDFday360predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday360plots.png", width = 27, height = 8, units = "cm")
ggarrange(SSTtraininglocDFday365actual, SSTtraininglocDFday365predicted.sequenFixedL, 
          SSTtraininglocDFday365predicted.sequenVaryLj, 
          labels = c("A", "B", "C"), label.x = 0.1, align = "h", ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top")
#ggsave("NNGPsequenSSTtraininglocDFday365plots.png", width = 27, height = 8, units = "cm")

