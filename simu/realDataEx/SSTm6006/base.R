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
