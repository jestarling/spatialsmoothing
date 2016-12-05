# Fixed Rank Kriging Function
# Purpose: Function to estimate sigma2_eps using
#          Cressie robust variogram estimate.

# Jennifer Starling & Giorgio Paulon
# November 2016
# Based on methodology in Katzfuss & Cressie, Tutorial on Fixed Rank Kriging of CO2 Data.
# SDS 385, Big Data, James Scott, Fall 2016

#INPUT:
# df = data frame, including columns named Lat and Lon containing coordinate values.

#OUTPUT:
# List object, containing
# vargram_cressie = estimated variogram object.
# sigma2_eps = estimated variance.

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