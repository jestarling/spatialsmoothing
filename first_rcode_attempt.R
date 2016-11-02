
rm(list=ls())
setwd('/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/Final Project')

library(Matrix) #For sparse matrix V_eps.
library(nlme) #For semi-variogram estimation.
library(rgl)
library(splines)
library(ggmap)
library(leaflet)
library(sp)
library(rgdal)
library(lubridate)

############################### 
### READ & PRE-PROCESS DATA ###
###############################

# Read data from .Rdata file
load(file = 'gammaradiation.Rdata')
data=data.raw

# Preprocess TIMESTAMP
data$timestamp <- ymd_hms(data$timestamp)
data$day <- day(data$timestamp)
data$month <- month(data$timestamp)
data$wday <- wday(data$timestamp)
data$hour <- hour(data$timestamp)

# Create dataframe with required fields
df <- as.data.frame(cbind(data$lon,data$lat, data$y, data$temp,data$day, data$month, data$wday, data$hour))
colnames(df) <- c('Lon','Lat','y','temp','day','month','wday','hour')

# Preprocess LOCATIONS: cut latitude and longitude that are outside UT campus
idx <- which(df$Lat > 30.277942 & df$Lat < 30.295346 & df$Lon > -97.742802 & df$Lon < -97.72277)
df <- df[idx,]

################# 
### PLOT DATA ###
#################
ggplot() + geom_point(data=df, aes(x=Lon,y=Lat),col='blue')
# Plot heatmap using ggmap
qmplot(Lon, Lat, data = df, colour = y, size = I(0.8), darken = .4, alpha = I(.6))

##################### 
### DE-TREND DATA ###
#####################

#Plot y versus latitude.
plot(df$Lat,df$y,pch=20)

#Plot y versus longitude.
plot(df$Lon,df$y,pch=20)

#Plot y versus temperature.
plot(df$temp,df$y,pch=20)

#Looks like no de-trending is required in this case.

#####################################################
### ANSCOMBE TRANSFORMATION TO NORMALIZE Y COUNTS ###
#####################################################

#Data y is in form of Poisson counts.
hist(df$y)
#Performing Anscombe variance-stabilizing transform to normalize counts.
y.norm = 2*sqrt(df$y + 3/8)
hist(y.norm)

#Add y.norm to data frame.
df$y.norm = y.norm

####################################
### Generate Basis Functions     ###
####################################
S = bs(df$y.norm) #Using B-Spline basis
#surface3d(S[,1], S[,2], S[,3])

####################################
### Estimate Parameters          ###
####################################
v_eps = rep(1,length(df$y.norm))
V_eps = Matrix(diag(v_eps),sparse=T)





