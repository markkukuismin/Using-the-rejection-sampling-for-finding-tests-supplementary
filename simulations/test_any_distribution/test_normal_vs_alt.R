
# Test that do samples come from some specific
# probability dist.

# Laplace dist. gives nice results. Skewed just
# as exponential dist. need still some work...

library(ggplot2)
library(LaplacesDemon)
library(goftest)
library(extraDistr)
library(energy)

ar_dist = function(x, f0 = "dnorm", ...){
  
  f0 = get(f0)
  
  n = length(x)
  
  fhat = kdevine::kde1d(x)
  fhat = kdevine::dkde1d(x, fhat)
  
  rhox = f0(x, ...)/fhat
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rho = (sum(a) + b)/n
  
  rho
  
}

set.seed(1)

n = 20 # 20, 30, 50

# Power (H_0 not true)

alt_dist = "runif" # t, mixture, logistic, runif

prob = 0.5

df = 2

M = 10^4

ks_stat = ar_stat = 
  energy_stat = ad_stat = 
  cvm_stat = rep(0, M)
  
for(j in 1:M){
  
  if(alt_dist == "t") x = rt(n, df = df)
  if(alt_dist == "logistic") x = rlogis(n)
  if(alt_dist == "mixture"){
    x1 = rnorm(n)
    x2 = rnorm(n, mean = 3)
    ind = rbinom(n, 1, prob = prob)
    x = ind*x1 + (1 - ind)*x2
  }
  if(alt_dist == "runif") x = runif(n)
  
  ks_stat[j] = ks.test(x, 
                       "pnorm",
                       mean = 0,
                       sd = 1)$p.value
  
  cvm_stat[j] = goftest::cvm.test(x, 
                                  "pnorm",
                                  mean = 0,
                                  sd = 1)$p.value
  
  ar_stat[j] = ar_dist(x, 
                       f0 = "dnorm",
                       mean = 0,
                       sd = 1)
  
  ad_stat[j] = goftest::ad.test(x,
                                null = "pnorm",
                                mean = 0,
                                sd = 1)$p.value
  
  energy_stat[j] = energy::mvnorm.test(x,
                                       R = 100)$p.value
  
  cat("\r", j)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  x = rnorm(n)
  
  ar_null[i] = ar_dist(x, 
                       f0 = "dnorm",
                       mean = 0,
                       sd = 1)
  
}

beta_ks = mean(ks_stat < 0.05)

beta_cvm = mean(cvm_stat < 0.05)

beta_ad = mean(ad_stat < 0.05)

beta_energy = mean(energy_stat < 0.05)

qar = quantile(ar_null, 0.05)
beta_ar = mean(ar_stat < qar)

# Type 1 error 

if(alt_dist == "t"){
  
  ks_er = ar_stat_er =
    cvm_stat_er = ad_stat_er = 
    energy_stat_er = rep(0, M)
  
  for(j in 1:M){
    
    x = rnorm(n)
    
    ks_er[j] = ks.test(x, 
                       "pnorm",
                       mean = 0,
                       sd = 1)$p.value
    
    cvm_stat_er[j] = goftest::cvm.test(x, 
                                       "pnorm",
                                       mean = 0,
                                       sd = 1)$p.value
    
    ad_stat_er[j] = goftest::ad.test(x,
                                     null = "pnorm",
                                     mean = 0,
                                     sd = 1)$p.value
    
    energy_stat_er[j] = energy::mvnorm.test(x,
                                            R = 100)$p.value
    
    ar_stat_er[j] = ar_dist(x, 
                            f0 = "dnorm",
                            mean = 0,
                            sd = 1)
    
    cat("\r", j)
    
  }
  
  error_ks = mean(ks_er < 0.05)
  
  error_cvm = mean(cvm_stat_er < 0.05)
  
  error_ad = mean(ad_stat_er < 0.05)
  
  error_energy = mean(energy_stat_er < 0.05)
  
  error_ar = mean(ar_stat_er < qar)
  
}else{
  
  error_ks = error_cvm = 
    error_ad = error_energy = 
    error_ar = 0
  
}

b = c(round(beta_ar, 2),
      round(beta_ks, 2),
      round(beta_cvm, 2),
      round(beta_energy, 2),
      round(beta_ad, 2))

names(b) = c("AR", "KS", "CVM", "E", "AD")

b

paste0(b, collapse = " & ")

er = c(round(error_ar, 2),
       round(error_ad, 2),
       round(error_energy, 2),
       round(error_cvm, 2),
       round(error_ks, 2)
)

names(er) = c("AR", "KS", "CVM", "E", "AD")

er

paste0(er, collapse = " & ")

#xx = seq(-4, 4, length.out = 1000)

#plot(xx, dnorm(xx), lwd = 2, type = "l")
#lines(xx, dt(xx, df=2), col = "red", lty = 2, lwd = 2)
