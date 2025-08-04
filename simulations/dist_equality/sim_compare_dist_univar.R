
source("functions/ar_dist.R")

set.seed(1)

n = 25

# mu under null

m_null = 0

m_alt = m_null + seq(-1, 1, length.out = 40)

# alt dist

alt_dist = "t"

M = 10^4

ks_stat = ar_stat = matrix(0, M, length(m_alt))

sig_diag = 1

I = diag(1, 2)

for(i in 1:length(m_alt)){
  
  for(j in 1:M){
    
    if(alt_dist == "norm") x = rnorm(n)
    if(alt_dist == "exp") x = rexp(n)
    if(alt_dist == "t") x = rt(n, df = 2) 
    if(alt_dist == "laplace") x = LaplacesDemon::rlaplace(n)
    x = x - mean(x) + m_null
    y = rnorm(n, mean = m_alt[i])
    
    ks_stat[j, i] = ks.test(x, y)$p.value
    
    rho = ar_dist(x, y)
    
    ar_stat[j, i] = rho
    
  }
  
  cat("\r", i)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  x = rnorm(n)
  y = rnorm(n)
  
  ar_null[i] = ar_dist(x, y)
  
}

beta_ks = apply(ks_stat, 2, function(x) mean(x < 0.05))

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# Type 1 error 

ks_er = ar_stat_er = matrix(0, M, length(m_alt))

for(i in 1:length(m_alt)){
  
  for(j in 1:M){
    
    x = rnorm(n, mean = m_alt[i])
    y = rnorm(n, mean = m_alt[i])
    
    ks_er[j, i] = ks.test(x, y)$p.value
    
    rho = ar_dist(x, y)
    
    ar_stat_er[j, i] = rho
    
    
  }
  
  cat("\r", i)
  
}


error_ks = apply(ks_er, 2, function(x) mean(x < 0.05))

error_ar = colMeans(ar_stat_er < qar)

##

Method = rep(c("ar",
               "ks"), 
             each = length(m_alt))

beta_v = c(beta_ar, beta_ks)

error_v = c(error_ar, error_ks)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  m = rep(m_alt, 2))

f_path = paste0("simulations/dist_equality/results/", 
                "dist_power_", alt_dist, "_", n,
                "n.txt")

write.table(Data, f_path, row.names = FALSE)
