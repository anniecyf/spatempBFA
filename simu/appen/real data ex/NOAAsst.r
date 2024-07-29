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
# length(na.omit(as.vector(sstArray[,,1]))) # 690903; 1440 * 720 = 1036800
sstMat <- matrix(sstArray, nLong * nLat, nTime) # (1440 * 720) * 365
locIndices <- which(apply((!is.na(sstMat)), 1, all)) # of length 645933; all locations with no NAs for all 365 days
Nu <- 350
sstTraining <- sstMat[locIndices, 1:Nu] # 645933 * 350
sstTestingTemp <- sstMat[locIndices, (Nu+1):nTime] # 645933 * 15
longTraining <- longAll[locIndices]; latTraining <- latAll[locIndices]
noNAsstDF <- cbind(longTraining, latTraining, sstMat[locIndices,], locIndices)
colnames(noNAsstDF) <- c("Longitude", "Latitude", paste0("Day", 1:nTime), "Location")
noNAsstDF <- as.data.frame(noNAsstDF)
save(noNAsstDF, file = "noNAsst2023DF.RData")
t <- 180
sstDFslice <- data.frame(long = longTraining, 
                         lat = latTraining, 
                         sst = sstTraining[,t]) # 645933 obs
sstPlot <- ggplot() + geom_tile(data = sstDFslice, aes(x = long, y = lat, fill = sst)) +
  theme(plot.background = element_rect(color = "white"), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top") + coord_sf() +
  scale_fill_gradient(name = "Sea Surface Temperature (degC)", low = "beige", high = "red")
world2 <- map_data(map = "world2", region = ".")
sstPlot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), color = "black", fill = "white") 
ggsave("sst2023nonNAt180.png", width = 20, height = 8, units = "cm")
# Nu <- 350
# sstTraining <- sstMat[,1:Nu] 
# # temp <- apply((!is.na(sstTraining)), 1, all)
# locIndices <- which(apply((!is.na(sstTraining)), 1, all)) # of length 647554; all locations with no NAs for time 1:Nu
# sstTraining <- sstTraining[locIndices,] # 647554 * 350
# longTraining <- longAll[locIndices]; latTraining <- latAll[locIndices]
t <- 180
sstSlice <- sstArray[,,t] # 1440 * 720
sstDFslice <- data.frame(long = rep(long, length(lat)), 
                         lat = rep(lat, each = length(long)), 
                         sst = as.vector(sstSlice)) # 1036800 obs
sstDFslice <- filter(sstDFslice, !is.na(sst)) # 690972 obs
sstPlot <- ggplot() + geom_tile(data = sstDFslice, aes(x = long, y = lat, fill = sst)) +
  theme(plot.background=element_blank(), panel.background = element_blank(), 
        axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
  scale_fill_gradient(low = "beige", high = "red") + coord_sf()
world2 <- map_data(map = "world2", region = ".")
sstPlot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) 
ggplot(world2) + geom_polygon(aes(x = long, y = lat, group = group), color = "black", fill = NA) +
  coord_sf() + theme_void()
# # quick map
# image(long, lat, sstSlice, col = rev(brewer.pal(10, "RdBu")))
# # levelplot of the slice
# longlatGrid <- expand.grid(long=long, lat=lat)
# cutPoints <- seq(-50, 50, by = 10)
# levelplot(sstSlice ~ long * lat, data = longlatGrid, 
#           at = cutPoints, cuts = length(cutPoints), pretty = TRUE, 
#           col.regions = (rev(brewer.pal(10,"RdBu"))))
# tUnits <- ncatt_get(SSTdata, "time", "units")
# tUnits
# cf <- CFtime(tUnits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
# cf
# timestamps <- as_timestamp(cf) # get character-string times
# timestamps
# timeCF <- CFparse(cf, timestamps) # parse the string into date components
# timeCF
# ncatt_get(SSTdata, "sst", "long_name")
# ncatt_get(SSTdata, "sst", "units") # "degC"
# ncatt_get(SSTdata, 0, "title")
# ncatt_get(SSTdata, 0, "institution")
# ncatt_get(SSTdata, 0, "source")
# ncatt_get(SSTdata, 0, "references")
# ncatt_get(SSTdata, 0, "history")
# ncatt_get(SSTdata, 0, "Conventions")