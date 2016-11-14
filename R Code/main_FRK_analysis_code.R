
rm(list=ls()) #Clean workspace.

###############################
#JENNIFER FILE PATH SETUP
setwd('/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/Final Project')

#GIORGIO FILE PATH SETUP


###############################

library(SpatialTools)
library(Matrix) #For sparse matrix V_eps.
library(gstat) #For semi-variogram estimation.
#library(rgl)
library(splines)
library(ggmap)
library(leaflet)
library(sp)
library(rgdal)
library(lubridate)

############################### 
### LOAD SOURCE FUNCTIONS   ###
###############################
source('./spatialsmoothing/R Code/em.R')
source('./spatialsmoothing/R Code/frk.R')
source('./spatialsmoothing/R Code/variogram_est.R')
source('./spatialsmoothing/R Code/bisquare.basis.R')

############################### 
### READ & PRE-PROCESS DATA ###
###############################

# Read data from .Rdata file
load(file = './Data/gammaradiation.Rdata')

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

# #Plot data versus each direction and versus temp.
# plot(df$Lat,df$y,pch=20)  #Plot y versus latitude.
# plot(df$Lon,df$y,pch=20)  #Plot y versus longitude.
# plot(df$temp,df$y,pch=20) #Plot y versus temperature.
# 
# #Temperature looks like a source of trend.  
# #Will de-trend based on temperature.
# mylm = lm(y~temp,data=df)
# betahat = mylm$coef
# x0 = model.matrix(mylm)
# y_detrended = df$y - x0 %*% betahat
# 
# #Plot de-trended data.
# plot(df$temp,df$y,pch=20)
# 
# Recheck that de-trended data does not have a lat or lon trend.
# par(mfrow=c(1,1))
# plot(df$temp,y_detrended,pch=20)
# plot(df$Lon,y_detrended,pch=20)
# plot(df$Lat,y_detrended,pch=20)

#De-trending by temp seemed to throw off all of the other plots.
#For now, will work with data with no de-trending performed.

#####################################################
### ANSCOMBE TRANSFORMATION TO NORMALIZE Y COUNTS ###
#####################################################

#Data y is in form of Poisson counts.
hist(df$y)
#Performing Anscombe variance-stabilizing transform to normalize counts.
y.norm = 2*sqrt(df$y + 3/8)
hist(y.norm)

#De-trend the mean.
mean.to.add = mean(y.norm)
y.norm = y.norm - mean(y.norm)
hist(y.norm)

#Add y.norm to data frame.
df$y.norm = y.norm

###########################################################
### Generate Basis Functions for Observed Locations     ###
###########################################################

#Set up data used to create centers.
fine.grid <- 40
x <- seq(-97.74280, -97.72277, length.out = fine.grid)
y <- seq(30.27794, 30.29534, length.out = fine.grid)

rangex <- range(x)[2] - range(x)[1]
rangey <- range(y)[2] - range(y)[1]

# Create the centers of the three scale levels
level1 <- as.matrix(expand.grid(seq(min(x) + rangex/5, max(x) - rangex/5, length.out = 2), seq(min(y) + rangey/5, max(y) - rangey/5, length.out = 2)))
level2 <- as.matrix(expand.grid(seq(min(x) + rangex/10, max(x) - rangex/10, length.out = 3), seq(min(y) + rangey/10, max(y) - rangey/10, length.out = 3)))
level3 <- as.matrix(expand.grid(seq(min(x) + rangex/20, max(x) - rangex/20, length.out = 4), seq(min(y) + rangey/20, max(y) - rangey/20, length.out = 4)))

#Create basis for observed coordinate locations.
S = bisquare.basis(coord = df[,1:2],level1, level2, level3)

####################################
### Estimate Parameters          ###
####################################

#-----------------------------------
#1. Estimate sigma2_eps with semivariogram.

#Assuming v_eps = 1 since measurement error not 
#expected to vary over Ds.
v_eps = rep(1,length(df$y.norm))

#Estimate sigma2_eps using
#Cressie robust variogram estimate.
#vargram = variogram_est(df)

#ALTERNATIVE:
#Variogram function to estimate sigma2_eps takes a long time to run on data this size.
#Saved vargram_cressie and sigma2_eps as R objects to load.
load(file='./spatialsmoothing/R Objects/vargram_cressie_ctr.Rdata')

sigma2_eps = vargram$sigma2_eps
vargram = vargram$vargram_cressie

#Code for saving vargram and sigma2_eps if need to rerun code.
#save(vargram, file='./spatialsmoothing/vargram_cressie.Rdata')

#-----------------------------------
#2. Estimate K and sigma2_xi with EM algorithm.
emEsts = em(S=S,z=y.norm,v=v_eps,sigeps=sigma2_eps,maxiter=50,avgtol=1E-2)

emEststest = em(S=S,z=y.norm,v=v_eps,sigeps=.25,maxiter=50,avgtol=1E-2)

####################################
### Fixed Rank Kriging           ###
####################################

#-----------------------------
#Call FRK function for smoothing:
# (Only smooth values at observed points.)

frkSmooth = frk(data=df, 
            K=emEsts$K,
            sigxi = emEsts$sigma2_xi,
            sige = sigma2_eps,
            v = v_eps,
            S = S, 
            Sp = S,
            goal="smooth")

#-----------------------------
#Call FRK function for prediction:

#Create a grid of y-values to predict.
lat_range = range(df$Lat)
lon_range = range(df$Lon)

pred.lat = seq(lat_range[1], lat_range[2],length.out = 100)
pred.lon = seq(lon_range[1], lon_range[2],length.out = 100)
pred.grid = expand.grid(Lon=pred.lon,Lat=pred.lat)

#Create a basis for predicted coordinates.
Sp = bisquare.basis(coord = pred.grid,level1, level2, level3)

#Call FRK function.
frkPred = frk(data=df, 
            pred_locs=pred.grid,
            K=emEsts$K,
            sigxi = emEsts$sigma2_xi,
            sige = sigma2_eps,
            v = v_eps,
            S = S,
            Sp = Sp,
            goal="predict")

#####################################################
### REVERSE ANSCOMBE TRANSFORM & ADD BACK MEAN    ###
#####################################################

#For smoothed data.
yhat = ((frkSmooth$pred[,3] + mean.to.add) / 2)^2 - 1/8
frkSmooth$pred = cbind(frkSmooth$pred,yhat,y.norm=df$y.norm,y=df$y)

#For predicted data.
yhat = ((frkPred$pred[,3] + mean.to.add) / 2)^2 - 1/8
frkPred$pred = cbind(frkPred$pred,yhat)

#Preview results:

#Smoothed observed points.
head(frkSmooth$pred)

#New predicted points.
head(frkPred$pred)
mean(frkPred$sig2FRK)
median(frkPred$sig2FRK)

####################################
### PLOT FRK RESULTS             ###
####################################
lat_lims= c(30.277942, 30.295346)
lon_lims = c(-97.742802, -97.72277)

#Save map for fast re-use.
bw.map <- get_map(location = c(-97.742802,30.277942,-97.72277,30.295346), 
                source = "osm",col='bw')

ggmap(bw.map) +
  geom_tile(data = as.data.frame(as.matrix(frkPred$pred)), 
                        aes(x = Lon, y = Lat, alpha = yhat.norm),
                        fill = 'red') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#Plot predicted FRK run.
plot_pred_dat = as.data.frame(as.matrix(frkPred$pred))
map.bw +
  
  scale_x_continuous(limits = lon_lims, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_lims, expand = c(0, 0))
  
  

#Plot FRK variances for predicted values.
plot_FRK_var = as.data.frame(as.matrix(cbind(frkPred$pred,FRKvar=frkPred$sig2FRK)))

map +
  geom_tile(data = plot_FRK_var, aes(x = Lon, y = Lat, alpha = FRKvar),
            fill = 'red') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

# #Try another predicted plot
# map.bw + 
#   geom_density2d(data=plot_pred_dat, aes(x=Lon, y=Lat, size=.3)) +
#   stat_density2d(data=plot_pred_dat, aes(x=Lon,y=Lat,
#                                              fill=..level.., alpha=..level..), size=.01, bins=16, geom="polygon")


# map + 
#   geom_density2d(data=plot_pred_dat, aes(x=Lon, y=Lat,size=.001)) +
#   stat_density2d(data=plot_pred_dat, aes(x=Lon,y=Lat,
#                                              fill=Predicted, alpha=Predicted), size=.001, geom="polygon")

#Plot smoothed (observed locations) values.
plot_smoothed_dat = as.data.frame(as.matrix(frkSmooth$pred))
map +
  geom_tile(data = plot_smoothed_dat, aes(x = Lon, y = Lat, alpha = Predicted),
            fill = 'red') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

map + 
  geom_density2d(data=plot_smoothed_dat, aes(x=Lon, y=Lat, size=.3)) +
    stat_density2d(data=plot_smoothed_dat, aes(x=Lon,y=Lat,
    fill=Predicted, alpha=Predicted), size=.01, bins=16, geom="polygon")




