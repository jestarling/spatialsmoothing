# EM Algorithm Function 
# Purpose: Estimates K and sigma2_xi values for fixed rank kriging using EM algorithm.

# Jennifer Starling & Giorgio Paulon
# November 2016
# Based on methodology in Katzfuss & Cressie, Tutorial on Fixed Rank Kriging of CO2 Data.
# SDS 385, Big Data, James Scott, Fall 2016

#INPUT:
# S = matrix, with each column equal to a basis vector.
# z = vector of observed responses
# v = variance vector, using all 1's.
# sigeps = sigma2_eps estimate (estimated using variogram)
# maxiter = maximum iterations allowed.
# avgtol = tolerance value for convergence criteria.

#OUTPUT
# List containing:
# K = estimated covariance matrix for fixed rank kriging.
# sigma2_xi = estimated fine-scale variance for fixed rank kriging.

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