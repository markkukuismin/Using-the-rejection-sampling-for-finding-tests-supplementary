
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
  
  n = nrow(x)
  
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

alt_dist = "mvnorm" # mvnorm, mvnormmix, mvt, rmvlogis, unif

p = 3

prob = 0.5

df = 2

M = 1000

ar_stat = energy_stat = 
  Tn_stat = hz_stat = 
  roys_stat = vage_stat =
  rep(0, M)

I = diag(1, p)

for(j in 1:M){
  
  if(alt_dist == "mvt"){
    
    x = mvtnorm::rmvt(n, 
                      sigma = I,
                      df = df)
    
  }
  
  if(alt_dist == "mvnorm"){
    
    x = mvtnorm::rmvnorm(n,
                         mean = rep(1, p),
                         sigma = I)
    
  }
  
  if(alt_dist == "mvnormmix"){
    
    x1 = mvtnorm::rmvnorm(n, 
                          mean = rep(0, p),
                          sigma = I)
    
    x2 = mvtnorm::rmvnorm(n,
                          mean = rep(3, p),
                          sigma = I)
    
    ind = rbinom(n, 1, prob = prob)
    x = ind*x1 + (1 - ind)*x2
    
  }
  
  if(alt_dist == "rmvlogis"){
    
    x = NonNorMvtDist::rmvlogis(n,
                                parm1 = rep(0, p),
                                parm2 = rep(1, p))
    
  }
  
  if(alt_dist == "unif"){
    
    x = matrix(runif(n*p), n, p)
    
  }
  
  ar_stat[j] = ar_dist(x)
  
  energy_stat[j] = energy::mvnorm.test(x,
                                       R = 100)$p.value
  
  Tn_stat[j] = as.numeric(mvnormalTest::mvnTest(x, B = 100)$mv.test["p-value"])
  
  hz_stat[j] = as.numeric(mvnormalTest::mhz(x)$mv.test[2])
  
  mvtemp = mvnormalTest::msw(x)
  
  roys_stat[j] = as.numeric(mvtemp$mv.test[1, 3])
  
  vage_stat[j] = as.numeric(mvtemp$mv.test[2, 3])
  
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

beta_hz = mean(hz_stat < 0.05)

beta_roys = mean(roys_stat < 0.05)

beta_vage = mean(vage_stat < 0.05)

# Type 1 error 

if(alt_dist == "mvnorm"){
  
  ar_stat_er = Tn_stat_er = 
    energy_stat_er = hz_stat_er = 
    roys_stat_er = vage_stat_er = rep(0, M)
  
  for(j in 1:M){
    
    x = MASS::mvrnorm(n, 
                      mu = rep(0, p),
                      Sigma = diag(1, p))
    
    energy_stat_er[j] = energy::mvnorm.test(x,
                                            R = 100)$p.value
    
    ar_stat_er[j] = ar_dist(x)
    
    Tn_stat_er[j] = mvnormalTest::mvnTest(x, 
                                          B = 100)$mv.test["p-value"]
    
    hz_stat_er[j] = mvnormalTest::mhz(x)$mv.test[2]
    
    mvtemp = mvnormalTest::msw(x)
    
    roys_stat_er[j] = mvtemp$mv.test[1, 3]
    
    vage_stat_er[j] = mvtemp$mv.test[2, 3]
    
    cat("\r", j)
    
  }
  
  error_energy = mean(energy_stat_er < 0.05)
  
  error_ar = mean(ar_stat_er < qar)
  
  error_Tn = mean(Tn_stat_er < 0.05)
  
  error_hz = mean(hz_stat_er < 0.05)
  
  error_roys = mean(roys_stat_er < 0.05)
  
  error_vage = mean(vage_stat_er < 0.05)
  
}else{
  
  error_energy = error_ar = 
    error_Tn = error_hz = 
    error_roys = error_vage = 0
  
}

b = c(round(beta_ar, 2),
      round(beta_energy, 2),
      round(beta_Tn, 2),
      round(beta_hz, 2),
      round(beta_roys, 2),
      round(beta_vage, 2))

names(b) = c("AR", "E", "Tn", "HZ", "Roys", "VA-GE")
b
paste0(b, collapse = " & ")
  
er = c(round(error_ar, 2),
       round(error_energy, 2),
       round(error_Tn, 2),
       round(error_hz, 2),
       round(error_roys, 2),
       round(error_vage, 2))

names(er) = c("AR", "E", "Tn", "HZ", "Roys", "VA-GE")
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

# hist(ar_stat_er, 
#      probability = T,
#      col = "pink",
#      add = T)

D = data.frame(pwr = b,
               er = er,
               dist = alt_dist,
               n = n
)

fpath = "simulations/test_any_distribution/results/mvdemo.txt"

write.table(D, 
            file = fpath,
            append = TRUE,
            row.names = FALSE)
