
# Test that do samples come from multivariate
# normal dist.

library(ggplot2)
library(LaplacesDemon)
library(goftest)
library(extraDistr)
library(energy)
library(mvnormalTest)
library(NonNorMvtDist)

ar_dist = function(x){
  
  f0 = mvtnorm::dmvnorm(x)
  
  n = length(x)
  
  fhat = kdevine::kdevine(x)
  fhat = kdevine::dkdevine(x, fhat)
  
  #u = VineCopula::pobs(x)
  #fhat = kdevine::kdevinecop(u)
  #fhat = kdevine::dkdevinecop(u, fhat)
  
  rhox = f0/fhat
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rho = (sum(a) + b)/n
  
  rho
  
}

set.seed(1)

n = 50 # 20, 30, 50

# Power (H_0 not true)

pdf = "rmvlogis" # mvnorm, mvnormmix, mvt, rmvlogis

p = 3

prob = 0.5

df = 2

M = 1000

ar_stat = energy_stat = Tn_stat = rep(0, M)

I = diag(1, p)

for(j in 1:M){
  
  if(pdf == "mvt"){
    
    x = mvtnorm::rmvt(n, 
                      sigma = I,
                      df = df)
    
  }
  
  if(pdf == "mvnorm"){
    
    x = mvtnorm::rmvnorm(n,
                         mean = rep(1, p),
                         sigma = I)
    
  }
  
  if(pdf == "mvnormmix"){
    
    x1 = mvtnorm::rmvnorm(n, 
                          mean = rep(0, p),
                          sigma = I)
    
    x2 = mvtnorm::rmvnorm(n,
                          mean = rep(3, p),
                          sigma = I)
    
    ind = rbinom(n, 1, prob = prob)
    x = ind*x1 + (1 - ind)*x2
    
  }
  
  if(pdf == "rmvlogis"){
    
    x = NonNorMvtDist::rmvlogis(n,
                                parm1 = rep(0, p),
                                parm2 = rep(1, p))
    
  }
  
  ar_stat[j] = ar_dist(x)
  
  energy_stat[j] = energy::mvnorm.test(x,
                                       R = 100)$p.value
  
  Tn_stat[j] = mvnormalTest::mvnTest(x, B = 100)$mv.test["p-value"]
  
  cat("\r", j)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  x = MASS::mvrnorm(n, 
                    mu = rep(0, p), 
                    Sigma = diag(1, p))
  
  ar_null[i] = ar_dist(x)
  
}

beta_energy = mean(energy_stat < 0.05)

beta_Tn = mean(Tn_stat < 0.05)

qar = quantile(ar_null, 0.05)
beta_ar = mean(ar_stat < qar)

# Type 1 error 

ar_stat_er = Tn_stat_er = energy_stat_er = rep(0, M)

for(j in 1:M){
  
  x = MASS::mvrnorm(n, 
                    mu = rep(0, p),
                    Sigma = diag(1, p))
  
  energy_stat_er[j] = energy::mvnorm.test(x,
                                          R = 100)$p.value
  
  ar_stat_er[j] = ar_dist(x)
  
  Tn_stat_er[j] = mvnormalTest::mvnTest(x, 
                                        B = 100)$mv.test["p-value"]
  
  cat("\r", j)
  
}

error_energy = mean(energy_stat_er < 0.05)

error_ar = mean(ar_stat_er < qar)

error_Tn = mean(Tn_stat_er < 0.05)

b = c(round(beta_ar, 2),
      round(beta_energy, 2),
      round(beta_Tn, 2))

names(b) = c("AR", "E", "Tn")
b
paste0(b, collapse = " & ")
  
er = c(round(error_ar, 2),
       round(error_energy, 2),
       round(error_Tn, 2))

names(er) = c("AR", "E", "Tn")
er
paste0(er, collapse = " & ")

hist(ar_stat, 
     probability = T,
     xlim = c(min(c(ar_null, ar_stat)),
              max(c(ar_null, ar_stat))))

hist(ar_null, 
     probability = T,
     col = "lightblue",
     add = T)

hist(ar_stat_er, 
     probability = T,
     col = "pink",
     add = T)

D = data.frame(pwr = b,
               er = er,
               dist = pdf,
               n = n
)

fpath = "simulations/test_any_distribution/results/mvdemo.txt"

write.table(D, 
            file = fpath,
            append = TRUE,
            row.names = FALSE)
