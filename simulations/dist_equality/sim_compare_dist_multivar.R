
library(fasano.franceschini.test)

source("functions/ar_dist_mv.R")

set.seed(1)

n = 100 # 50, 100

p = 2

# mu under null

m_null = rep(0, p)

m = seq(-1, 1, length.out = 40)

# alt dist

alt_dist = "t"

M = 10^3

ff_stat = ar_stat = matrix(0, M, length(m))

I = diag(1, p)

for(i in 1:length(m)){
  
  m_alt = m_null
  
  m_alt[p] = m_alt[p] + m[i]
  
  for(j in 1:M){
    
    if(alt_dist == "norm") x = MASS::mvrnorm(n = n, mu = m_null, Sigma = I)
    if(alt_dist == "t") x = LaplacesDemon::rmvt(n = n, mu = m_null, S = I, df = 5)
    if(alt_dist == "laplace") x = LaplacesDemon::rmvl(n = n, mu = m_null, Sigma = I)
    y = MASS::mvrnorm(n = n, mu = m_alt, Sigma = I)
    
    ff_stat[j, i] = fasano.franceschini.test(x, y, 
                                             verbose = FALSE)$p.value
    
    rho = ar_dist_mv(x, y)
    
    ar_stat[j, i] = rho
    
  }
  
  cat("\r", i)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  x = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = I)
  y = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = I)
  
  ar_null[i] = ar_dist_mv(x, y)
  
}

beta_ff = apply(ff_stat, 2, function(x) mean(x < 0.05))

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# Type 1 error 

ff_er = ar_stat_er = matrix(0, M, length(m))

for(i in 1:length(m)){
  
  m_alt = m_null
  
  m_alt[p] = m_alt[p] + m[i]
  
  for(j in 1:M){
    
    x = MASS::mvrnorm(n = n, mu = m_alt, Sigma = I)
    y = MASS::mvrnorm(n = n, mu = m_alt, Sigma = I)
    
    ff_er[j, i] = fasano.franceschini.test(x, y, 
                                           verbose = FALSE)$p.value
    
    rho = ar_dist_mv(x, y)
    
    ar_stat_er[j, i] = rho
    
    
  }
  
  cat("\r", i)
  
}


error_ff = apply(ff_er, 2, function(x) mean(x < 0.05))

error_ar = colMeans(ar_stat_er < qar)

##

Method = rep(c("ar",
               "ff"), 
             each = length(m))

beta_v = c(beta_ar, beta_ff)

error_v = c(error_ar, error_ff)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  m = rep(m, 2))

f_path = paste0("simulations/dist_equality/results/", 
                "dist_power_mv", alt_dist, "_", n,
                "n.txt")

write.table(Data, f_path, row.names = FALSE)
