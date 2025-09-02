
ar_mc_pvalue = function(x = NULL, M = 1000, f0 = "dnorm", ...){
  
  if(is.null(x)) stop("x is missing with no default")
  
  f0_mc = sub("d", "r", f0)
  
  f0_mc = get(f0_mc)
  
  f0 = get(f0)
  
  n = length(x)
  
  fhat = kdevine::kde1d(x)
  fhat = kdevine::dkde1d(x, fhat)
  
  rhox = f0(x, ...)/fhat
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rhos = ifelse(rhox > 1, 1, rhox)
  
  rho_m = (sum(a) + b)/n
  
  rho_mc = rep(0, M)
  
  for(i in 1:M){
    
    x_mc = f0_mc(n = n, ...)
    
    fhat = kdevine::kde1d(x_mc)
    fhat = kdevine::dkde1d(x_mc, fhat)
    
    rhox = f0(x_mc, ...)/fhat
    
    rhox[is.na(rhox)] = 0
    
    a = rhox[rhox < 1]
    b = sum(rhox >= 1)
    
    rho_mc[i] = (sum(a) + b)/n
    
  }
  
  p_value = (sum(rho_mc <= rho_m) + 1)/(M + 1)
  
  return(list(rho_m = rho_m, 
              rho_mc = rho_mc, 
              rhos = rhos,
              p_value = p_value))
  
}