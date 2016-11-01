
rm(list=ls())
setwd('~/Desktop/Semester 1/Stat Models for Big Data/Project')

library(ggmap)
library(leaflet)
library(sp)
library(rgdal)
library(lubridate)

# Read data from .Rdata file
load(file = 'Data/gammaradiation.Rdata')

# Preprocess TIMESTAMP
data$timestamp <- ymd_hms(data$timestamp)
data$day <- day(data$timestamp)
data$month <- month(data$timestamp)
data$wday <- wday(data$timestamp)
data$hour <- hour(data$timestamp)

# Create dataframe with required fields
df <- as.data.frame(cbind(data$lon,data$lat, data$y, data$day, data$month, data$wday, data$hour))
colnames(df) <- c('Lon','Lat','y','day','month','wday','hour')

# Preprocess LOCATIONS: cut latitude and longitude that are outside UT campus
idx <- which(df$Lat > 30.277942 & df$Lat < 30.295346 & df$Lon > -97.742802 & df$Lon < -97.72277)
df <- df[idx,]

################# 
### PLOT DATA ###
#################


# # Plot data day by day using sp package
# coordinates(df) = ~Lon+Lat
# 
# for (i in sort(unique(df$month))){
#   for (j in sort(unique(df$day[df$month == i]))){
#     pdf(file = paste('./Plot/trial_month_', i, '_day_', j, '.pdf', sep = ''))
#     temp <- spplot(df[df$month == i & df$day == j,], "y", cex = 0.5)
#     print(temp)
#     dev.off()
#   }
# }


# Plot heatmap using ggmap
qmplot(Lon, Lat, data = df, colour = y, size = I(0.8), darken = .4, alpha = I(.6))


# # If we want to use Leaflet plots (not recommended for big data, just take a subset)
# rd=.5
# op=.8
# clr="blue"
# m = leaflet() %>% addProviderTiles("CartoDB.Positron")
# m %>% addCircles(lng = df$Lon[1:5000], lat = df$Lat[1:5000], radius = rd, opacity = op, color = clr)

