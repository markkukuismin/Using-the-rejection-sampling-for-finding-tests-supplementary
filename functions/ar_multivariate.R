
ar_multivariate = function(X, mu0 = NULL, r = NULL){
  
  p = ncol(X)
  
  if(is.null(mu0)) mu0 = rep(0, p)
  
  if(is.null(r)) r = diag(1, p)
  
  Xm = colMeans(X, na.rm = TRUE)
  
  n = nrow(X)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  #S = (df - 2)/df*r
  #rho = mvtnorm::dmvnorm(t0, sigma = r)/mvtnorm::dmvt(t0, df = df, sigma = S, log = FALSE)
  b = mvtnorm::dmvnorm(t0, sigma = r)
  if(b == 0){
    rho = 0
  }else{
    rho = b/mvtnorm::dmvt(t0, df = df, sigma = r, log = FALSE)
  }
  
  rho = ifelse(rho >= 1, 1, rho)
  
  return(list(rho = rho))
  
}