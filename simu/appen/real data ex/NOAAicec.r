rm(list = ls())
library(tidyverse)
library(ncdf4)
library(CFtime)
library(ggmap)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex", sep = "/"))
iceCdata <- nc_open("icec.day.mean.2023.nc")
print(iceCdata)
long <- ncvar_get(iceCdata, "lon") # of length 1440
summary(long) #; hist(long)
lat <- ncvar_get(iceCdata, "lat") # of length 720
summary(lat) #; hist(lat)
time <- ncvar_get(iceCdata, "time") # of length 365
summary(time) #; hist(time)
nLong = length(long); nLat = length(lat); nTime = length(time)
longAll <- rep(long, nLat); latAll = rep(lat, each = nLong)
icecArray <- ncvar_get(iceCdata, "icec")
fillValue <- ncatt_get(iceCdata, "icec", "_fillValue")
icecArray[icecArray==fillValue$value] <- NA
dim(icecArray) # 1440 * 720 * 365
sum(is.na(icecArray)) # 329763362; 1440 * 720 * 365 = 378432000
icecMat <- matrix(icecArray, nLong * nLat, nTime) # (1440 * 720) * 365
locIndices <- which(apply((!is.na(icecMat)), 1, all)) # of length 67874; all locations with no NAs for all 365 days
longTraining <- longAll[locIndices]; latTraining <- latAll[locIndices]
noNAicecDF <- cbind(longTraining, latTraining, icecMat[locIndices,], locIndices)
colnames(noNAicecDF) <- c("Longitude", "Latitude", paste0("Day", 1:nTime), "Location")
noNAicecDF <- as.data.frame(noNAicecDF)
# save(noNAicecDF, file = "noNAicec2023DF.RData")
world2 <- map_data(map = "world2", region = ".")
t <- 365
icecDFslice <- data.frame(long = longTraining, lat = latTraining, 
                          icec = noNAicecDF[, (t+2)]) # 67874 obs
icecPlot <- ggplot(icecDFslice) + geom_tile(aes(x = long, y = lat, fill = icec)) + coord_sf() +
  scale_fill_gradientn(name = "Sea Ice Concentration (percent)", 
                       colors = rev(topo.colors(30))) +
  labs(x = "Longitude", y = "Latitude") + # for picking a representative rectangle 
  theme(axis.ticks = element_blank(), legend.position = "top")
icecPlot 
# ggsave("icec2023nonNAday365topo.png", width = 20, height = 8, units = "cm")
icecPlotBlue <- ggplot(icecDFslice) + geom_tile(aes(x = long, y = lat, fill = icec)) + coord_sf() +
  scale_fill_distiller(name = "Sea Ice Concentration (percent)", 
                       palette = "Blues") +
  labs(x = "Longitude", y = "Latitude") + # for picking a representative rectangle 
  theme(axis.ticks = element_blank(), legend.position = "top", 
        panel.grid.major = element_line(color = "forestgreen", linewidth = 0.5), 
        panel.grid.minor = element_line(color = "forestgreen", linewidth = 0.5))
icecPlotBlue 
# ggsave("icec2023nonNAday365blue.png", width = 20, height = 8, units = "cm")
icecPlotviridis <- ggplot(icecDFslice) + geom_tile(aes(x = long, y = lat, fill = icec)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Ice Concentration (percent)", direction = -1) +
  labs(x = "Longitude", y = "Latitude") + # for picking a representative rectangle 
  theme(axis.ticks = element_blank(), legend.position = "top")
icecPlotviridis 
# ggsave("icec2023nonNAday365viridis.png", width = 20, height = 8, units = "cm")
icecPlotviridis +
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "transparent") 
# ggsave("icec2023nonNAday365viridisWithWorldMap.png", width = 20, height = 8, units = "cm")
world2 <- map_data(map = "world2", region = ".")
icecDFsliceAll <- data.frame(long = longAll, lat = latAll, 
                             icec90 = icecMat[,90], icec181 = icecMat[,181],
                             icec273 = icecMat[,273], icec365 = icecMat[,365])
iceCday90plot <- ggplot(icecDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = icec90, 
                linetype = "Regions with NA Sea \nIce Concentration Values \nat Day 90 in 2023")) +
  scale_fill_distiller(name = "Sea Ice Concentration (percent)", 
                       palette = "Blues", na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
iceCday90plot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
ggsave("icec2023day90greyNA.png", width = 20, height = 8, units = "cm")
iceCday181plot <- ggplot(icecDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = icec181, 
                linetype = "Regions with NA Sea \nIce Concentration Values \n at Day 181 in 2023")) +
  scale_fill_distiller(name = "Sea Ice Concentration (percent)", 
                       palette = "Blues", na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
iceCday181plot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
ggsave("icec2023day181greyNA.png", width = 20, height = 8, units = "cm")
iceCday273plot <- ggplot(icecDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = icec273, 
                linetype = "Regions with NA Sea \nIce Concentration Values \n at Day 273 in 2023")) +
  scale_fill_distiller(name = "Sea Ice Concentration (percent)", 
                       palette = "Blues", na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
iceCday273plot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
ggsave("icec2023day273greyNA.png", width = 20, height = 8, units = "cm")
iceCday365plot <- ggplot(icecDFsliceAll) + coord_sf() +
  geom_tile(aes(x = long, y = lat, fill = icec365, 
                linetype = "Regions with NA Sea \nIce Concentration Values \n at Day 365 in 2023")) +
  scale_fill_distiller(name = "Sea Ice Concentration (percent)", 
                       palette = "Blues", na.value = "darkgrey") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top", legend.key = element_blank()) + 
  guides(linetype = guide_legend(override.aes = list(fill = "darkgrey"), 
                                 label.position = "right",
                                 title = "", order = 2),
         fill = guide_colorbar(order = 1)) 
iceCday365plot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "white") 
ggsave("icec2023day365greyNA.png", width = 20, height = 8, units = "cm")

# ncatt_get(iceCdata, "icec", "long_name")
# ncatt_get(iceCdata, "icec", "units") # percent
# ncatt_get(iceCdata, 0, "title")
# ncatt_get(iceCdata, 0, "institution")
# ncatt_get(iceCdata, 0, "source")
# ncatt_get(iceCdata, 0, "references")
# ncatt_get(iceCdata, 0, "history")
# ncatt_get(iceCdata, 0, "Conventions")
rm(list = ls())
library(tidyverse)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex", sep = "/"))
load("noNAsst2023DF.RData")
load("noNAicec2023DF.RData")
noNAlocIndices <- intersect(noNAicecDF$Location, noNAsstDF$Location) # of length 67524
sstDF <- noNAsstDF %>% filter(Location %in% noNAlocIndices)
icecDF <- noNAicecDF %>% filter(Location %in% noNAlocIndices)
SSTicecDF <- rbind(sstDF, icecDF)
# save(SSTicecDF, file = "noNAicecSST2023DF.RData")
