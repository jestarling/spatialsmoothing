# Fixed Rank Kriging Function
# Purpose: Performs fixed rank kriging.

# Jennifer Starling & Giorgio Paulon
# November 2016
# Based on methodology in Katzfuss & Cressie, Tutorial on Fixed Rank Kriging of CO2 Data.
# SDS 385, Big Data, James Scott, Fall 2016

#INPUT:
# data = a data frame with column 1 = longitudes, column 2 = latitudes, column 3 = y values.
# pred_locs = a grid of coordinates at which to predict y values.  Should be (mx2), with col 1 = longitudes, col 2 = latitudes.
# K = estimated K covariance matrix.
# sigxi = estimated sigma2_xi fine-sample variance value. 
# sige = estimated sigma2_eps overall variance value. (Est from semivariogram.)
# v = estimated vector of vars, we use 1's.
# S = matrix of basis functions for observed locations.
# Sp = matrix of basis functions for predicted locations.
# goal = 'predict' or 'smooth'.  
#     Smooth only kriges values at observed locations.
#     Predict interpolates values at grid of new locations.

#OUTPUT
# List containing:
# pred = matrix of coordinates with predicted values.
# sig2FRK = FRK variance estimate, if predicting new values.
# K = estimated covariance matrix.
# sig2_eps = sigma2_eps estimated overall variance.
# sig2_xi = sigma2_xi estimated fine-scale variance.

frk = function(data,pred_locs=NULL,K,sigxi,sige,v,S,Sp,goal="predict"){
  
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
  pred_with_locs = cbind(Lon=lon_pred,Lat=lat_pred,Yhat=pred)
  colnames(pred_with_locs) = c("Lon","Lat","yhat.norm")
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