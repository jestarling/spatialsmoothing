
rm(list=ls()) #Clean workspace.

###############################
#JENNIFER FILE PATH SETUP
setwd('/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/Final Project')

#GIORGIO FILE PATH SETUP
setwd('~/Desktop/Semester 1/Stat Models for Big Data/Project/')

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
#idx <- which(df$Lat > 30.277942 & df$Lat < 30.295346 & df$Lon > -97.742802 & df$Lon < -97.72277) # Original crop
# idx <- which(df$Lat > 30.277942 & df$Lat < 30.292040 & df$Lon > -97.742802 & df$Lon < -97.72277) # Cropping zones without data
idx <- which(df$Lat > 30.277942 & df$Lat < 30.293142 & df$Lon > -97.742802 & df$Lon < -97.727602) # Square crop


df <- df[idx,]

# Find locations corresponding to the POLICE STATION
idx.police <- which(df$Lat > 30.283093 & df$Lat < 30.285740 & df$Lon > -97.732350 & df$Lon < -97.727853)

################# 
### PLOT DATA ###
#################

ggplot() + geom_point(data=df, aes(x=Lon,y=Lat),col='blue')

##################### 
### DE-TREND DATA ###
#####################

# #Plot data versus each direction and versus temp.
# plot(df$Lat,df$y,pch=20)  #Plot y versus latitude.
# plot(df$Lon,df$y,pch=20)  #Plot y versus longitude.
# plot(df$temp,df$y,pch=20) #Plot y versus temperature.

#jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/Final Project/spatialsmoothing/LaTeX Files/Images/detrending_plots2.jpg')
par(mfrow=c(3,1))
plot(df$Lat,df$y,pch=20)  #Plot y versus latitude.
plot(df$Lon,df$y,pch=20)  #Plot y versus longitude.
plot(df$temp,df$y,pch=20) #Plot y versus temperature.
#dev.off()

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
x <- seq(-97.742802, -97.727602, length.out = fine.grid)
y <- seq(30.277942, 30.293142, length.out = fine.grid)
pred.grid = expand.grid(x, y)

rangex <- range(x)[2] - range(x)[1]
rangey <- range(y)[2] - range(y)[1]

# Create the centers of the three scale levels
level1 <- as.matrix(expand.grid(seq(min(x) + rangex/4, max(x) - rangex/4, length.out = 3), seq(min(y) + rangey/4, max(y) - rangey/4, length.out = 3)))
level2 <- as.matrix(expand.grid(seq(min(x) + rangex/6, max(x) - rangex/6, length.out = 4), seq(min(y) + rangey/6, max(y) - rangey/6, length.out = 4)))
level3 <- as.matrix(expand.grid(seq(min(x) + rangex/10, max(x) - rangex/10, length.out = 5), seq(min(y) + rangey/10, max(y) - rangey/10, length.out = 5)))
level3 <- level3[-c(24, 25), ]

# Choose the three scales for the basis functions
#scale <- c(2E-2, 1E-2, 5E-3)
scale <- c(1E-2, 7E-3, 5E-3)

# #Create basis for predicted coordinate locations
# S = bisquare.basis(coord = pred.grid, scale, level1, level2, level3)
# 
# r1 <- dim(level1)[1]
# r2 <- dim(level2)[1]
# r3 <- dim(level3)[1]
# 
# z <- rowSums(S[,(r1+r2+1):(r1+r2+r3)])
# z <- S[,30]
# z <- matrix(z, fine.grid, fine.grid, byrow = T)
# 
# # Define colors
# nrz <- nrow(z)
# ncz <- ncol(z)
# jet.colors <- colorRampPalette(c("blue", "red"))
# nbcol <- 100
# color <- jet.colors(nbcol)
# zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# facetcol <- cut(zfacet , nbcol)
# 
# # # Plot
# persp(x, y, z, col = color[facetcol], phi = 30, theta = 45)
# # Contour plot
# filled.contour(x, y, z, color.palette = heat.colors, asp = 1, plot.axes={points(level3)})

#Create basis for observed coordinate locations.
S = bisquare.basis(coord = df[,1:2], scale, level1, level2, level3)

#Plot data including basis functions in red.
# # Plot heatmap using ggmap
qmplot(Lon, Lat, data = df, colour = y, size = I(0.8), darken = .4, alpha = I(.6)) + 
  geom_point(data=as.data.frame(level1), aes(x=Var1, y=Var2), color="red", size=4, alpha=0.5) + 
  geom_point(data=as.data.frame(level2), aes(x=Var1, y=Var2), color="red", size=2, alpha=0.5) + 
  geom_point(data=as.data.frame(level3), aes(x=Var1, y=Var2), color="red", size=1, alpha=0.5)

####################################
### Estimate Parameters          ###
####################################

#-----------------------------------
#1. Estimate sigma2_eps with semivariogram.

#Assuming v_eps = 1 since measurement error not 
#expected to vary over Ds.
v_eps = rep(1,length(df$y.norm))

#Estimate sigma2_eps using
#Cressie robust variogram estimate: we here take out the POLICE STATION COORDINATES
df.nopolice <- df[-idx.police,]
vargram = sigma2_eps_vargram_est(df.nopolice)

# RUN ONLY if you want to get the 4 semivariograms
# # Divide the data in 4 quadrants
# idx21 <- which(df.nopolice$Lat < min(df.nopolice$Lat) + 0.5*(max(df.nopolice$Lat) - min(df.nopolice$Lat)) & df$Lon < min(df.nopolice$Lon) + 0.5*(max(df.nopolice$Lon) - min(df.nopolice$Lon)))
# idx11 <- which(df.nopolice$Lat >= min(df.nopolice$Lat) + 0.5*(max(df.nopolice$Lat) - min(df.nopolice$Lat)) & df$Lon < min(df.nopolice$Lon) + 0.5*(max(df.nopolice$Lon) - min(df.nopolice$Lon)))
# idx22 <- which(df.nopolice$Lat < min(df.nopolice$Lat) + 0.5*(max(df.nopolice$Lat) - min(df.nopolice$Lat)) & df$Lon >= min(df.nopolice$Lon) + 0.5*(max(df.nopolice$Lon) - min(df.nopolice$Lon)))
# idx12 <- which(df.nopolice$Lat >= min(df.nopolice) + 0.5*(max(df.nopolice$Lat) - min(df.nopolice$Lat)) & df$Lon >= min(df.nopolice$Lon) + 0.5*(max(df.nopolice$Lon) - min(df.nopolice$Lon)))
# vargram11 = sigma2_eps_vargram_est(df.nopolice[idx11,])
# vargram12 = sigma2_eps_vargram_est(df.nopolice[idx12,])
# vargram21 = sigma2_eps_vargram_est(df.nopolice[idx21,])
# vargram22 = sigma2_eps_vargram_est(df.nopolice[idx22,])
# 
# plot(vargram11$vargram_cressie)
# plot(vargram12$vargram_cressie)
# plot(vargram21$vargram_cressie)
# plot(vargram22$vargram_cressie)


#ALTERNATIVE:
#Variogram function to estimate sigma2_eps takes a long time to run on data this size.
#Saved vargram_cressie and sigma2_eps as R objects to load.
#load(file='./spatialsmoothing/R Objects/vargram_cressie_nopolice.Rdata')

sigma2_eps = vargram$sigma2_eps
vargram = vargram$vargram_cressie

#Code for saving vargram and sigma2_eps if need to rerun code.
#save(vargram, file='./spatialsmoothing/vargram_cressie.Rdata')

#-----------------------------------
#2. Estimate K and sigma2_xi with EM algorithm.
emEsts = em(S=S,z=y.norm,v=v_eps,sigeps=sigma2_eps,maxiter=50,avgtol=1E-2)

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
Sp = bisquare.basis(coord = pred.grid, scale, level1, level2, level3)

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
results.Smooth = cbind(frkSmooth$pred,y.norm=df$y.norm,yhat,y=df$y)

#For predicted data.
yhat = ((frkPred$pred[,3] + mean.to.add) / 2)^2 - 1/8
results.Pred = cbind(frkPred$pred,yhat)

#Preview results:

#Smoothed observed points.
head(results.Smooth)  

#New predicted points.
head(cbind(results.Pred,frkPred$sig2FRK))
mean(frkPred$sig2FRK)

####################################
### PLOT FRK RESULTS             ###
####################################
#lat_lims= c(30.277942, 30.295346)
lat_lims= lat_range
lon_lims = lon_range


#Save map for fast re-use.
bw.map <- get_map(location = as.numeric(rbind(lon_lims, lat_lims)), 
                  source = "osm",col='bw')

#Plot predicted FRK run.
ggmap(bw.map) +
  geom_tile(data = as.data.frame(as.matrix(results.Pred)), 
            aes(x = Lon, y = Lat, alpha = pred),
            fill = 'red') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#Plot FRK variances for predicted values.
plot_FRK_var = as.data.frame(as.matrix(cbind(results.Pred,FRKvar=frkPred$sig2FRK)))
ggmap(bw.map) +
  geom_tile(data = plot_FRK_var, aes(x = Lon, y = Lat, alpha = FRKvar),
            fill = 'red') + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#Plot histogram of residuals for normality.
#Note: for smoothing only.
hist(results.Smooth[,4] - results.Smooth[,3],xlab='Residuals',ylab='Frequency',
     main='Residuals for Smoothed Data')

#Plot Kriging Variances boxplot:
boxplot(frkPred$sig2FRK)
