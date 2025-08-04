
ar_compare_group_means = function(X){
  
  p = ncol(X)
  
  n = nrow(X)
  
  mu0 = rep(mean(X), p)
  Y = t(X) - mu0
  R0 = tcrossprod(Y)/(n - 1)
  diag(R0) = mean(diag(R0))
  
  R = cov(X)
  diag(R) = mean(diag(R))
  Xm = colMeans(X)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(t0, sigma = R)
  b1 = mvtnorm::dmvt(t0, df = df, sigma = R0, log = FALSE)
  
  if(b0 == 0) rho = 0
  if(b0 != 0) rho = b0/b1
  
  return(list(rho = rho))
  
}