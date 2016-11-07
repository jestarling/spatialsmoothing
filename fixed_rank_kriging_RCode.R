
rm(list=ls())

###############################
#JENNIFER FILE PATH SETUP
setwd('/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/')

#GIORGIO FILE PATH SETUP


###############################

library(SpatialTools)
library(Matrix) #For sparse matrix V_eps.
library(Rcpp) #For FRK function, estimating p3.
library(RcppArmadillo) #For FRK function, estimating p3.
library(nlme) #For semi-variogram estimation.
library(gstat)
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

#Plot data versus each direction and versus temp.
plot(df$Lat,df$y,pch=20)  #Plot y versus latitude.
plot(df$Lon,df$y,pch=20)  #Plot y versus longitude.
plot(df$temp,df$y,pch=20) #Plot y versus temperature.

#Temperature looks like a source of trend.  
#Will de-trend based on temperature.
mylm = lm(y~temp,data=df)
betahat = mylm$coef
x0 = model.matrix(mylm)
y_detrended = df$y - x0 %*% betahat

#Plot de-trended data.
plot(df$temp,df$y,pch=20)

#Recheck that de-trended data does not have a lat or lon trend.
par(mfrow=c(1,1))
plot(df$temp,y_detrended,pch=20)
plot(df$Lon,y_detrended,pch=20)
plot(df$Lat,y_detrended,pch=20)

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

####################################
### Generate Basis Functions     ###
####################################
create_basis = function(data){
  Sbasis = bs(data) #Using b-splines basis.
  S = as.matrix(Sbasis[1:nrow(Sbasis),]) #Drop 'basis' format.
    #Required for matrix multiplication.
  return(S)
}

S = create_basis(df$y.norm) #Using B-Spline basis


####################################
### Estimate Parameters          ###
####################################

#-----------------------------------
#1. Estimate sigma2_eps with semivariogram.

#Assuming v_eps = 1 since measurement error not 
#expected to vary over Ds.
v_eps = rep(1,length(df$y.norm))

#Function to estimate sigma2_eps using
#Cressie robust variogram estimate.
sigma2_eps_vargram_est = function(df){
  require(sp)
  
  # convert data frame into a spatial data frame object
  coordinates(df)= ~ Lat + Lon
  
  #Create variogram (Cressie Robust version)
  vargram = variogram(df$y.norm~1, data=df,cressie=T)
  
  #Fit a line to the variogram.
  #Estimate sigma2_eps as the intercept.
  vg_y = vargram$gamma / 2  
    #Because vargram models 2gamma, and fits 
    #line gamma(h) = gamma(0+) + bh
                    
  vg_x = vargram$dist
  fit = lm(vg_y ~ vg_x)
  sigma2_eps = fit$coef[1] 
  
  return(list(vargram_cressie = vargram, sigma2_eps=sigma2_eps))
}

vargram = sigma2_eps_vargram_est(df)

#ALTERNATIVE:
#Variogram function to estimate sigma2_eps takes a long time to run on data this size.
#Saved vargram_cressie and sigma2_eps as R objects to load.
load(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/spatialsmoothing/R Objects/variogram_cressie.Rdata')
load(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/spatialsmoothing/R Objects/sigma2_eps.Rdata')
vargram = vargram_cressie
vargram
sigma2_eps

#Code for saving vargram and sigma2_eps if need to rerun code.
#save(sigma2_eps, file='/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/spatialsmoothing/sigma2_eps.Rdata')
save(vargram, file='/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/spatialsmoothing/vargram_cressie_ctr.Rdata')

#-----------------------------------
#1. Estimate K and sigma2_xi with EM algorithm.
z=y.norm
v=v_eps
sigeps = sigma2_eps
maxiter=10
avgtol=1E-6

em = function(S,z=y.norm,v=v_eps,sigeps=sigma2_eps,maxiter=10,avgtol=1E-6){
  
  #Set up starting values for K and sigma2_xi per pg 13.
  r = ncol(S)         #Number of basis vectors.
  n = length(S[,1])  #Number of observations.
  
  diagV = v
  diagV2 = rep(1,n)
  
  varest = var(z)
  K_old = .9 * varest * diag(1,r)
  sig2 = .1*varest
  t=1
  converged=0
  
  #Iterate until convergence.
  for (t in 1:maxiter){
    
    #Update helper variables.
    diagDinv = as.vector(1/(sig2[t] %*% diagV2 + sigeps %*% diagV))
    DInv = sparseMatrix(i=1:n,j=1:n,x=diagDinv)
    tempt = solve(solve(K_old) + t(S) %*% DInv %*% S)
    
    #---------------------
    #Update K.
    SigInv2 = tcrossprod(tempt,S) %*% DInv
    KSDInv = tcrossprod(K_old,S) %*% DInv
    KSSigInv = KSDInv - KSDInv %*% S %*% SigInv2
    muEta = KSSigInv %*% z
    SigEta = K_old - KSSigInv %*% S %*% K_old
    K_new = SigEta + tcrossprod(muEta)
    
    #---------------------
    #Update sigma2_xi.
    muEps = sig2[t] * (DInv %*% z - DInv %*% S %*% (SigInv2 %*% z))
    trSigInv = sum(diag((DInv))) - sum(diag((SigInv2 %*% DInv %*% S)))
    sig2[t+1] = as.numeric((1/n) * (n * sig2[t] - 
                          (sig2[t])^2 * trSigInv + 
                           crossprod(muEps)))
    #---------------------
    #Check convergence.
    diff = sum((K_new-K_old)^2)  +  (sig2[t+1]-sig2[t])^2  
    if (diff < (avgtol * r^2)){
      converged=1
      break
    } #End convergence check.
    
    #Update K_old
    K_old = K_new
    
  } #End for loop.
  
  #Return results.
  return(list( K = K_new, sigma2_xi = sig2[t], t=t+1, converged=converged))
}

emEsts = em(S=S,z=y.norm,v=v_eps,sigeps=sigma2_eps,maxiter=50,avgtol=1E-2)

####################################
### Fixed Rank Kriging           ###
####################################


#pred.grid = rbind(df[,1:2],pred.grid)

pred.grid=df[,1:2]

#Call FRK function.

#Plot predicted values.

#FRK Function
data=df
pred_locs=pred.grid
K = emEsts$K
sigxi= emEsts$sigma2_xi
sige = sigma2_eps

#Try calling FRK function for smoothing:
#Alternatively, just predict values at observed points.

myFRK = frk(data=df, 
            K=emEsts$K,
            sigxi = emEsts$sigma2_xi,
            sige = sigma2_eps,
            v = v_eps,
            goal="smooth")


#Try calling FRK function for prediction:

#Create a grid of y-values to predict.
lat_range = range(df$Lat)
lon_range = range(df$Lon)

pred.lat = seq(lat_range[1], lat_range[2],length.out = 100)
pred.lon = seq(lon_range[1], lon_range[2],length.out = 100)
pred.grid = expand.grid(pred.lon,pred.lat)
colnames(pred.grid) = c("Lon","Lat")

Sp = S[1:10000,] #DELETE LATER

myFRKpred = frk(data=df, 
            pred_locs=pred.grid,
            K=emEsts$K,
            sigxi = emEsts$sigma2_xi,
            sige = sigma2_eps,
            v = v_eps,
            goal="predict")

frk = function(data,pred_locs=NULL,K,sigxi,sige,v,goal="predict"){
  
  #-------------------------
  #INPUT VALIDATION CHECKS:
  
  #1. Check goal is either "smooth" or "predict".
  valid_goals = c("smooth","predict")
  if (!(goal %in% valid_goals)) stop("Goal must be 'smooth' or 'predict'.")
  
  #2. Check pred_locs provided if goal = predict.
  if (goal=="predict" & is.null(pred_locs)) stop("pred_locs required for prediction.")
 
  #3. If pred_locs != NULL and goal = smooth, warn user that only observed locs used.
  if(!is.null(pred_locs) & goal=="smooth") print("Warning: For smoothing, predicted locations not included in observed data will be ignored. Only observed locations are used.")
  
  #-------------------------
  #DATA PREP:
  
  #Data values.
  lon = data[,1]
  lat = data[,2]
  z = data$y.norm
  
  #Coordinates of predicted locations.
  #If doing smoothing, set predicted locations equal to observed locations.
  if (is.null(pred_locs) || goal=="smooth") pred_locs = df[,1:2]
  lon_pred = pred_locs[,1]
  lat_pred = pred_locs[,2]

  #Set up dimensions.
  r = ncol(S)         #Number of basis vectors.
  n = length(S[,1])  #Number of observations.
  m = length(lon_pred) #Number of predicted obs.
  
  V = sparseMatrix(i=1:n,j=1:n,x=v)
  V2 = sparseMatrix(i=1:n,j=1:n,x=rep(1,n))

  #Set up basis functions for actual and predicted coordinates.
  S = create_basis(y.norm)
  #Sp=S
  
  #-------------------------
  #DATA SMOOTHING
  
  #Diagonal error var matrix.
  DInv = (sigxi * V2 + sige * V)
  diag(DInv) = 1/diag(DInv)
  
  #Calc rxr part of inverse of Sigma.
  temp = solve(solve(K) + crossprod(S,DInv) %*% S)

  if (goal == "smooth"){
    #Predict (NORMALIZED) smoothed residuals for observed locations only.
    SigInvZ = DInv %*% z - DInv %*% S %*% (tcrossprod(temp,S) %*% DInv %*% z)
    pred = Sp %*% (tcrossprod(K,S) %*% SigInvZ) + sigxi * SigInvZ
  }

  if (goal == "predict"){
    #Indicator matrix setup:
    #Flag pred locations that exist as observed locations.
    loc_obs = paste(lon,lat,sep=",")
    loc_pred = paste(lon_pred,lat_pred,sep=",")
    
    #Save indices corresponding to observed coords.
    idx_obs_in_pred = match(loc_pred,loc_obs,nomatch=0)
    matched_i = which(idx_obs_in_pred != 0)
    matched_j = idx_obs_in_pred[which(idx_obs_in_pred != 0)]
    
    #Indicator matrix for observed locations. 
    E = sparseMatrix(dims=c(m,n),i=matched_i, j=matched_j)
    
    #Predict (NORMALIZED) smoothed residuals for predicted locations.
    SigInvZ = DInv %*% z - DInv %*% S %*% (tcrossprod(temp,S) %*% DInv %*% z)
    pred = Sp %*% (tcrossprod(K,S) %*% SigInvZ) + sigxi * E %*% SigInvZ
  }
  
  #Save predicted values with location coordinates.
  pred_with_locs = cbind(Lon=lon_pred,Lat=lat_pred,Yhat.norm=pred)
  
  #-------------------------
  #CALCULATE FRK VARIANCE
  #Due to computational intensity, only performing this operation
  #if doing interpolation, not smoothing.
  sig2FRK = 999
  
  if(goal=="predict"){
    #Find rows () corresponding to observed locations.
    index = summary(E)$i
  
    #Part 1 of variance.
    SpK = Sp %*% K
    p1 = apply(SpK,1,crossprod)
  
    #Part 3 of variance.
    KSSigInv = tcrossprod(K,S) %*% DInv - 
        (tcrossprod(K,S) %*% DInv %*% S %*% temp) %*% crossprod(S,DInv)
  
    SpKSSigInvSK = Sp %*% (KSSigInv %*% S %*% K)
    SigInv1 = DInv %*% S %*% temp
    SigInv2 = crossprod(S,DInv)
  
    p3 = rep(0,m)
    p3 = rowSums(SpKSSigInvSK * Sp)
    
    if (length(index)>0){
      for (i in 1:length(index)){
        ind = index[i]
        ESigInvE = as.numeric(DInv[ind,ind] - SigInv1[ind,] %*% SigInv2[,ind])
        p3[i] = p3[i] + 2 * sigxi * Sp[i,] %*% KSSigInv[,ind] + sigxi^2 * ESigInvE
      } 
    }
    #Add up pieces of FRK variance.
    sig2FRK = p1 + sigxi - p3
  }

  #Return function values.
 
  return(list(pred = pred_with_locs,
              sig2FRK=sig2FRK,
              K=K,
              sig2_eps = sige,
              sig2_xi = sigxi ))
 
} #End FRK function.

