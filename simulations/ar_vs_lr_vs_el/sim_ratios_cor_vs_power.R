
library(mvhtests)

source("functions/ar_two_groups.R")
source("functions/MLE0.R")

set.seed(1)

r = seq(-0.99, 0.99, length.out = 50)

r = c(r, 0)
r = sort(r)
r = unique(r)

n = 52

# mu under null

a = rnorm(1, sd = 4)

mu0 = c(a, a)

# mu under alternative

# Effect = 0.4, 0.6, 0.8 e.g., c(0, 0.4), c(0.4, 0.4)

d = c(0, 0.8)

mu = mu0 + d

# The distribution: norm, t, pexp, unif, laplace

pdist = "norm"

if(pdist == "t") df = 3 # for the t-distribution

# power

m = length(r)

M = 10^4

el_stat = lr_stat = 
  ar_stat = ar_stat_hat = matrix(0, M, m)

for(i in 1:m){
  
  Sigma = matrix(r[i], 2, 2)
  
  diag(Sigma) = 1
  
  if(pdist == "t"){
    Sigma = Sigma*((df-2)/df)
  }
  
  for(j in 1:M){
    
    if(pdist == "t"){
      X = mvtnorm::rmvt(n, delta = mu, sigma = Sigma, df = df) 
    }
    if(pdist == "norm"){
      X = MASS::mvrnorm(n, mu = mu, Sigma = Sigma) 
    }
    if(pdist == "laplace"){
      X = LaplacesDemon::rmvl(n, mu = mu, Sigma = Sigma)
    }
    if(pdist == "pexp"){
      X = LaplacesDemon::rmvpe(n, mu = mu, Sigma = Sigma, kappa = 2)
    }
    if(pdist == "unif"){
      X = LaplacesDemon::rmvpe(n, mu = mu, Sigma = Sigma, kappa = 10^3)
    }
    
    Sigma_hat = cor(X)
    
    x = X[, 1]
    y = X[, 2]
    
    ar_stat[j, i] = ar_two_groups(x, y, 
                                  r = Sigma, 
                                  mu0 = mu0)$rho
    ar_stat_hat[j, i] = ar_two_groups(x, y, 
                                      r = Sigma_hat, 
                                      mu0 = mu0)$rho
    
    el_stat[j, i] = mvhtests::el.test1(X, 
                                       mu = mu0)$res$Pval
    
    S0 = MLE0(X, mu0 = mu0)
    S1 = cov(X)*(n-1)/n
    
    lr_stat[j, i] = n*log(det(S0)/det(S1))
    
  }
  
  cat("\r", i)
  
}

# null distribution for the ar-statistic

ar_null = rep(0, M)
ar_null_hat = rep(0, M)
lr_null = rep(0, M)

I = diag(1, 2)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = mu0, Sigma = I)
  
  x = X[, 1]
  y = X[, 2]
  
  Sigma_hat = cor(X)
  
  ar_null[j] = ar_two_groups(x, y,
                             mu0 = mu0,
                             r = I)$rho
  ar_null_hat[j] = ar_two_groups(x, y, 
                                 mu0 = mu0,
                                 r = Sigma_hat)$rho
  
  
  S0 = MLE0(X, mu0 = mu0)
  S1 = cov(X)*(n-1)/n
  
  lr_null[j] = n*log(det(S0)/det(S1))
  
}

# Just checking the null of LR

plot(density(lr_null), lwd = 1)
xx = seq(0, 10, length.out = 1000)
lines(xx, dchisq(xx, df = 2), 
      lwd = 2, lty = 2, col = "red")

# AR test

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

qar_hat = quantile(ar_null_hat, 0.05)
beta_ar_hat = colMeans(ar_stat_hat < qar_hat) # sample correlation

# Empirical likelihood

beta_el = colMeans(el_stat < 0.05)

# Likelihood ratio

qlr = qchisq(0.95, df = 2)
beta_lr = colMeans(lr_stat > qlr)

# very close to the empirical power,

beta_lr
colMeans(lr_stat > quantile(lr_null, 0.95))

# Type 1 error 

mu = mu0

el_stat_er = lr_stat_er = 
  ar_stat_er = ar_stat_hat_er = matrix(0, M, m)

for(i in 1:m){
  
  Sigma = matrix(r[i], 2, 2)
  
  diag(Sigma) = 1
  
  if(pdist == "t"){
    Sigma = Sigma*((df-2)/df)
  }
  
  for(j in 1:M){
    
    if(pdist == "t"){
      X = mvtnorm::rmvt(n, delta = mu, sigma = Sigma, df = df) 
    }
    if(pdist == "norm"){
      X = MASS::mvrnorm(n, mu = mu, Sigma = Sigma) 
    }
    if(pdist == "laplace"){
      X = LaplacesDemon::rmvl(n, mu = mu, Sigma = Sigma)
    }
    if(pdist == "pexp"){
      X = LaplacesDemon::rmvpe(n, mu = mu, Sigma = Sigma, kappa = 2)
    }
    if(pdist == "unif"){
      X = LaplacesDemon::rmvpe(n, mu = mu, Sigma = Sigma, kappa = 10^3)
    }
    
    Sigma_hat = cor(X)
    
    x = X[, 1]
    y = X[, 2]
    
    ar_stat_er[j, i] = ar_two_groups(x, 
                                     y, 
                                     r = Sigma,
                                     mu0 = mu0)$rho
    
    ar_stat_hat_er[j, i] = ar_two_groups(x, 
                                         y, 
                                         r = Sigma_hat,
                                         mu0 = mu0)$rho
    
    el_stat_er[j, i] = mvhtests::el.test1(X, 
                                          mu = mu0)$res$Pval
    
    S0 = MLE0(X, mu0 = mu0)
    S1 = cov(X)*(n-1)/n
    
    lr_stat_er[j, i] = n*log(det(S0)/det(S1))
    
    
  }
  
  cat("\r", i)
  
}

error_ar = colMeans(ar_stat_er < qar)
error_ar_hat = colMeans(ar_stat_hat_er < qar_hat)

# Empirical likelihood test

error_el = colMeans(el_stat_er < 0.05)

# LR test

error_lr = colMeans(lr_stat_er > qlr)

##

Method = rep(c("ar_sample_cor", "ar_pop_cor",
               "empirical_likelihood",
               "LR"), 
             each = length(r))

beta_v = c(beta_ar_hat, beta_ar, beta_el, beta_lr)

error_v = c(error_ar_hat, error_ar, 
            error_el, error_lr)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  r = rep(r, 4))

d = gsub("\\.", "", as.character(d))

d = paste0(d, collapse = "")

f_path = paste0("simulations/ar_vs_lr_vs_el/results/", 
                pdist, "_", n, "n_", d, "d.txt")

write.table(Data, f_path, row.names = FALSE)
