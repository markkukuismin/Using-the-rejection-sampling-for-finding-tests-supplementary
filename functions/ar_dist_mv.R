
ar_dist_mv = function(x, y){
  
  n = length(x)
  
  fhatx = kdevine::kdevine(x)
  fhaty = kdevine::kdevine(y)
  
  fhatxx = kdevine::dkdevine(x, fhatx)
  fhatyy = kdevine::dkdevine(x, fhaty)
  
  rhox = fhatyy/fhatxx
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rhox = (sum(a) + b)/n
  
  ##
  
  fhatxx = kdevine::dkdevine(y, fhatx)
  fhatyy = kdevine::dkdevine(y, fhaty)
  
  rhoy = fhatxx/fhatyy
  
  rhoy[is.na(rhoy)] = 0
  
  a = rhoy[rhoy < 1]
  b = sum(rhoy >= 1)
  
  rhoy = (sum(a) + b)/n
  
  min(c(rhox, rhoy))
  
}