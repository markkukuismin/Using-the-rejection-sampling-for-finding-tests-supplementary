
library(ggplot2)
library(energy)

ar_dist = function(x){
  
  n = length(x)
  
  m = mean(x)
  
  s = sd(x)
  
  fhat = kdevine::kde1d(x)
  
  fhat = kdevine::dkde1d(x, fhat)
  
  rhox = dnorm(x, mean = m, sd = s)/fhat
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rho = (sum(a) + b)/n
  
  return(list(rho = rho, rhox = rhox))
  
}

ar_compare_group_means = function(X){
  
  p = ncol(X)
  
  n = nrow(X)
  
  mu0 = rep(mean(X, na.rm = TRUE), p)
  
  R = cov(X)
  
  Xm = colMeans(X, na.rm = TRUE)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(X, mean = mu0, sigma = R)
  b1 = mvtnorm::dmvt(X, 
                     df = df, 
                     delta = Xm,
                     sigma = R, 
                     log = FALSE)
  
  rhox = b0/b1
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rho = (sum(a) + b)/n
  
  return(list(rho = rho, rhox = rhox))
  
}

set.seed(1)

N = c(10, 30, 60)

# mu under null

M = 10^4

Mm = 50

for(j in 1:length(N)){
  
  Qgmean = Qgof = matrix(0, 3, 2)
  
  colnames(Qgmean) = 
    colnames(Qgof) = c("poisbin", "mc")
  
  rownames(Qgmean) = 
    rownames(Qgof) = c("0.01", "0.05", "0.1")
  
  n = N[j]
  
  # ar null
  
  rho_null = rep(0, M)
  
  rhox = matrix(0, n, M)
  
  for(i in 1:M){
    
    x = rnorm(n)
    
    D = ar_dist(x)
    
    rho_null[i] = D$rho
    
    rhox[, i] = D$rhox
    
  }
  
  rhox[rhox > 1] = 1
  
  pp = rowSums(rhox[, 1:Mm])
  
  pp = pp/Mm
  
  Qgof[, 1] = poibin::qpoibin(c(0.01, 0.05, 0.1), pp = pp)/n
  
  Qgof[, 2] = quantile(rho_null, c(0.01, 0.05, 0.1))
  
  ##
  
  # Group means
  
  rho_null = rep(0, M)
  
  rhox = matrix(0, n, M)
  
  for(i in 1:M){
    
    X = MASS::mvrnorm(n = n,
                      mu = rep(0, 2),
                      Sigma = diag(1, 2))
    
    D = ar_compare_group_means(X)
    
    rho_null[i] = D$rho
    
    rhox[, i] = D$rhox
    
  }
  
  rhox[rhox > 1] = 1
  
  pp = rowSums(rhox[, 1:Mm])
  
  pp = pp/Mm
  
  Qgmean[, 1] = poibin::qpoibin(c(0.01, 0.05, 0.1), pp = pp)/n
  
  Qgmean[, 2] = quantile(rho_null, c(0.01, 0.05, 0.1))
 
  print(round(Qgmean, 3))
  
  print(round(Qgof, 3))
   
  cat("\r", j)
  
}
