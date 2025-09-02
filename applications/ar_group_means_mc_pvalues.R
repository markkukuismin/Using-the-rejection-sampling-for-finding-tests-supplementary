
ar_group_means_mc_pvalues = function(X, M = 1000, paired = FALSE){
  
  p = ncol(X)
  
  n = nrow(X)
  
  mu0 = rep(mean(X, na.rm = TRUE), p)
  
  R = cov(X, use = "complete.obs")
  if(!paired) R = diag(diag(R))
  
  Xm = colMeans(X, na.rm = TRUE)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(t0, sigma = R)
  b1 = mvtnorm::dmvt(t0, df = df, sigma = R, log = FALSE)
  
  if(b0 == 0) rho = 0
  if(b0 != 0) rho = b0/b1
  
  rho_m = ifelse(rho >= 1, 1, rho)
  
  rho_mc = rep(0, M)
  
  for(i in 1:M){
    
    X = MASS::mvrnorm(n = n, 
                      mu = rep(0, p), 
                      Sigma = diag(1, p))
    
    mu0 = rep(mean(X, na.rm = TRUE), p)
    
    R = cov(X, use = "complete.obs")
    if(!paired) R = diag(diag(R))
    
    Xm = colMeans(X, na.rm = TRUE)
    
    t0 = sqrt(n)*(Xm - mu0)
    
    df = n - 1
    b0 = mvtnorm::dmvnorm(t0, sigma = R)
    b1 = mvtnorm::dmvt(t0, df = df, sigma = R, log = FALSE)
    
    if(b0 == 0) rho = 0
    if(b0 != 0) rho = b0/b1
    
    rho_mc[i] = ifelse(rho >= 1, 1, rho)
    
  }
  
  p_value = (sum(rho_mc <= rho_m) + 1)/(M + 1)
  
  return(list(rho_m = rho_m, 
              rho_mc = rho_mc, 
              p_value = p_value))
  
}