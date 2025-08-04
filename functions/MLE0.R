
MLE0 = function(X, mu0){
  
  p = ncol(X)
  
  n = nrow(X)
  
  S = matrix(0, p, p)
  
  for(i in 1:n){
    
    Y = X[i,] - mu0
    
    Y = matrix(Y, nrow = p)
    
    S = S + Y%*%t(Y)
    
  }
  
  S/n
  
}