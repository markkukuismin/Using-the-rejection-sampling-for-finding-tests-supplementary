
# Test that do samples come from some specific
# probability dist.

# Laplace dist. gives nice results. Skewed just
# as exponential dist. need still some work...

library(ggplot2)
library(LaplacesDemon)
library(goftest)
library(extraDistr)

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

n = 100 # 20, 30, 50, 100

# Power (H_0 not true)

pdf = "t" #

df = 5 # 2, 5, 11

M = 10^4

ks_stat = ar_stat = 
  cvm_stat = rep(0, M)
  
for(j in 1:M){
  
  x = rt(n, df = df)
  
  ks_stat[j] = ks.test(x, 
                       "pnorm",
                       mean = 0,
                       sd = 1)$p.value
  
  cvm_stat[j] = cvm.test(x, 
                         "pnorm",
                         mean = 0,
                         sd = 1)$p.value
  
  ar_stat[j] = ar_dist(x, 
                       f0 = "dnorm",
                       mean = 0,
                       sd = 1)
  
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

qar = quantile(ar_null, 0.05)
beta_ar = mean(ar_stat < qar)

# Type 1 error 

ks_er = ar_stat_er =
  cvm_stat_er = rep(0, M)

for(j in 1:M){
  
  x = rnorm(n)
  
  ks_er[j] = ks.test(x, 
                       "pnorm",
                       mean = 0,
                       sd = 1)$p.value
  
  cvm_stat_er[j] = cvm.test(x, 
                         "pnorm",
                         mean = 0,
                         sd = 1)$p.value
  
  ar_stat_er[j] = ar_dist(x, 
                       f0 = "dnorm",
                       mean = 0,
                       sd = 1)
  
  cat("\r", j)
  
}

error_ks = mean(ks_er < 0.05)

error_cvm = mean(cvm_stat_er < 0.05)

error_ar = mean(ar_stat_er < qar)

beta_ar
beta_ks
beta_cvm

error_ar
error_cvm
error_ks

xx = seq(-4, 4, length.out = 1000)

plot(xx, dnorm(xx), lwd = 2, type = "l")
lines(xx, dt(xx, df=2), col = "red", lty = 2, lwd = 2)

hist(rnorm(n), probability = TRUE)
hist(rt(n, df = 2), 
     probability = TRUE, 
     add = TRUE,
     col = "lightblue")
