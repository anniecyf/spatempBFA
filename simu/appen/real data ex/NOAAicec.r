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
Nu <- 350
icecTraining <- icecMat[locIndices, 1:Nu] # 67874 * 350
icecTestingTemp <- icecMat[locIndices, (Nu+1):nTime] # 67874 * 15
longTraining <- longAll[locIndices]; latTraining <- latAll[locIndices]
noNAicecDF <- cbind(longTraining, latTraining, icecMat[locIndices,], locIndices)
colnames(noNAicecDF) <- c("Longitude", "Latitude", paste0("Day", 1:nTime), "Location")
noNAicecDF <- as.data.frame(noNAicecDF)
save(noNAicecDF, file = "noNAicec2023DF.RData")
t <- 180
icecDFslice <- data.frame(long = longTraining, 
                          lat = latTraining, 
                          icec = icecTraining[,t]) # 67874 obs
icecPlot <- ggplot() + geom_tile(data = icecDFslice, aes(x = long, y = lat, fill = icec)) +
  theme(plot.background = element_rect(color = "white"), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "top") + coord_sf() +
  scale_fill_gradient(name = "Sea Ice Concentration (percent)", low = "lightblue", high = "darkblue") 
world2 <- map_data(map = "world2", region = ".")
icecPlot + 
  geom_polygon(data = world2, mapping = aes(x = long, y = lat, group = group), color = "black", fill = "white") 
ggsave("icec2023nonNAt180.png", width = 20, height = 8, units = "cm")


t <- 180
icecSlice <- icecArray[,,t] # 1440 * 720
icecDFslice <- data.frame(long = rep(long, length(lat)), lat = rep(lat, each = length(long)), 
                          icec = as.vector(icecSlice)) 
icecDFslice <- filter(icecDFslice, !is.na(icec))
ggplot() + geom_tile(data = icecDFslice, aes(x = long, y = lat, fill = icec)) +
  theme(plot.background=element_blank(), panel.background = element_blank(), 
        axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") + coord_sf()
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
save(SSTicecDF, file = "noNAicecSST2023DF.RData")
