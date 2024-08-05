rm(list=ls())
library(tidyverse)
projDirec <- "C:/Users/annie/OneDrive - National University of Singapore/Documents/PhD/research/first paper/spatempBFA"
setwd(paste(projDirec, "simu/appen/real data ex", sep = "/"))
load("noNAsst2023DF.RData")
SSTtrainingDF <- noNAsstDF %>% filter(Longitude >= 100 & Longitude <= 147 & Latitude >= 30 & Latitude <= 70) # 6006 obs
setwd("SSTm6006")
#save(SSTtrainingDF, file = "SSTtrainingDF2023.RData")
SSTtrainingDFday365 <- SSTtrainingDF %>% select(Longitude, Latitude, Day365) %>% rename(SST = Day365)
SSTtrainingDFday365plot <- ggplot(SSTtrainingDFday365) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = SST)) + coord_sf() +
  scale_fill_viridis_c(name = "Sea Surface Temperature (degC)", direction = -1) +
  theme(axis.ticks = element_blank(), legend.position = "top", legend.key = element_blank())
SSTtrainingDFday365plot
ggsave("SSTtrainingDFday365plot.png", width = 10, height = 8, units = "cm")
library(ggmap)
world2 <- map_data(map = "world2", region = ".")
world2sub <- world2 %>% filter(long >= 100 & long <= 147 & lat >= 30 & lat <= 70)
world2sub <- world2sub %>% filter(!region %in% c("China", "Russia", "Mongolia")) #c("China", "Russia", "Mongolia")
SSTtrainingDFday365plot +
  geom_polygon(data = world2sub, mapping = aes(x = long, y = lat, group = group), 
                                       color = "black", fill = "transparent") +
  theme(plot.background = element_blank(), panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) 