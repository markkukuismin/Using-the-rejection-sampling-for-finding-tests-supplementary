
ar_dist = function(x, y){
  
  n = length(x)
  
  fhatx = density(x)
  fhaty = density(y)
  
  fhatxx = approx(fhatx, xout = x)$y
  fhatyy = approx(fhaty, xout = x)$y
  
  rhox = fhatyy/fhatxx
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rhox = (sum(a) + b)/n
  
  ##
  
  fhatxx = approx(fhatx, xout = y)$y
  fhatyy = approx(fhaty, xout = y)$y
  
  rhoy = fhatxx/fhatyy
  
  rhoy[is.na(rhoy)] = 0
  
  a = rhoy[rhoy < 1]
  b = sum(rhoy >= 1)
  
  rhoy = (sum(a) + b)/n
  
  min(c(rhox, rhoy))
  
}
