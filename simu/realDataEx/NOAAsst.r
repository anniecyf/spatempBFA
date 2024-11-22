rm(list = ls())
library(tidyverse)
library(ncdf4)
library(CFtime)
library(ggmap)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex", sep = "/"))
SSTdata <- nc_open("sst.day.mean.2023.nc")
print(SSTdata)
long <- ncvar_get(SSTdata, "lon") # of length 1440
summary(long) #;hist(long)
lat <- ncvar_get(SSTdata, "lat") # of length 720
summary(lat) #;hist(lat)
time <- ncvar_get(SSTdata, "time") # of length 365
summary(time) #;hist(time)
nLong = length(long); nLat = length(lat); nTime = length(time)
longAll <- rep(long, nLat); latAll = rep(lat, each = nLong)
sstArray <- ncvar_get(SSTdata, "sst")
fillValue <- ncatt_get(SSTdata, "sst", "_fillValue")
sstArray[sstArray == fillValue$value] <- NA
dim(sstArray) # 1440 * 720 * 365
sum(is.na(sstArray)) # 126235256; 1440 * 720 * 365 = 378432000
sstMat <- matrix(sstArray, nLong * nLat, nTime) # (1440 * 720) * 365 = 1036800 * 365
locIndices <- which(apply((!is.na(sstMat)), 1, all)) # of length 645933; all locations with no NAs for all 365 days in 2023
longTraining <- longAll[locIndices]; latTraining <- latAll[locIndices]
noNAsstDF <- cbind(longTraining, latTraining, sstMat[locIndices,], locIndices)
colnames(noNAsstDF) <- c("Longitude", "Latitude", paste0("Day", 1:nTime), "Location")
noNAsstDF <- as.data.frame(noNAsstDF)
# save(noNAsstDF, file = "noNAsst2023DF.RData")
world2 <- map_data(map = "world2", region = ".")
t <- 365
sstDFslice <- data.frame(long = longTraining, lat = latTraining, 
                         sst = noNAsstDF[, (t+2)]) # 645933 obs
sstPlot <- ggplot(sstDFslice) + geom_tile(aes(x = long, y = lat, fill = sst)) + coord_sf() +
  scale_fill_gradientn(name = "Sea Surface Temperature (degC)", 
                       colors = rev(heat.colors(30))) +
  labs(x = "Longitude", y = "Latitude") + # for picking a representative rectangle 
  theme(axis.ticks = element_blank(), legend.position = "top")
sstPlot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent") 
# ggsave("sst2023nonNAday365.png", width = 20, height = 8, units = "cm")
sstPlotviridis <- ggplot(sstDFslice) + geom_tile(aes(x = long, y = lat, fill = sst)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  labs(x = "Longitude", y = "Latitude") + # for picking a representative rectangle 
  theme(axis.ticks = element_blank(), legend.position = "top")
sstPlotviridis + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent") 
# ggsave("sst2023nonNAday365viridis.png", width = 20, height = 8, units = "cm")
world2 <- map_data(map = "world2", region = ".")
sstDFsliceAll <- data.frame(long = longAll, lat = latAll, 
                            sst90 = sstMat[,90], sst181 = sstMat[,181],
                            sst273 = sstMat[,273], sst365 = sstMat[,365])
viridisPlotSSTday90 <- ggplot(sstDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = sst90, 
                linetype = "Regions with NA \nTemperature Values\nat Day 90 in 2023")) +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", 
                       direction = -1, na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
viridisPlotSSTday90 + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
# ggsave("sst2023day90greyNA.png", width = 20, height = 8, units = "cm") 
viridisPlotSSTday181 <- ggplot(sstDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = sst181, 
                linetype = "Regions with NA \nTemperature Values\nat Day 181 in 2023")) +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", 
                       direction = -1, na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
viridisPlotSSTday181 + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
# ggsave("sst2023day181greyNA.png", width = 20, height = 8, units = "cm") 
viridisPlotSSTday273 <- ggplot(sstDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = sst273, 
                linetype = "Regions with NA \nTemperature Values\nat Day 273 in 2023")) +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", 
                       direction = -1, na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
viridisPlotSSTday273 + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
# ggsave("sst2023day273greyNA.png", width = 20, height = 8, units = "cm")
viridisPlotSSTday365 <- ggplot(sstDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = sst365, 
                linetype = "Regions with NA \nTemperature Values\nat Day 365 in 2023")) +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", 
                       direction = -1, na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
viridisPlotSSTday365 + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
# ggsave("sst2023day365greyNA.png", width = 20, height = 8, units = "cm")
