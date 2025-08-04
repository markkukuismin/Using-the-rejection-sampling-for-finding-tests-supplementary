
source("functions/ar_compare_group_means.R")

set.seed(1)

r = seq(-0.99, 0.99, length.out = 50)

r = c(r, 0)

r = sort(r)

pwr::pwr.t.test(d = 0.4,
                sig.level = 0.05,
                power = 0.8,
                type = "paired",
                alternative = "two.sided")

n = 52

# mu under null

mu0 = rep(rnorm(1), 2)

# mu under alternative

# Effect = 0.4, 0.6, 0.8 e.g., effect = c(0, 0.6)

mu = mu0 + c(0, 0.4)

d = abs(diff(mu))

# The distribution: norm, t, pexp, unif, laplace

pdist = "norm"

if(pdist == "t") df = 3 # for the t-distribution

# power

m = length(r)

M = 10^4

t_paired = t_unpaired = ar_stat = matrix(0, M, m)

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
    
    y = X[, 1]
    x = X[, 2]
    
    t_paired[j, i] = t.test(x, y, paired = T)$statistic
    t_unpaired[j, i] = t.test(x, y, var.equal = T)$statistic
    ar_stat[j, i] = ar_compare_group_means(X)$rho
    
  }
  
  cat("\r", i)
  
}

# null distribution for the ar-statistic

ar_null = rep(0, M)

I = diag(1, 2)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = rep(0, 2), Sigma = I)
  
  x = X[, 1]
  y = X[, 2]
  
  Sigma_hat = cor(X)
  
  ar_null[j] = ar_compare_group_means(X)$rho

}

ql = qt(0.025, df = n - 1)
qu = qt(0.975, df = n - 1)

f = function(x) mean(x > qu | x < ql)

beta_paired = apply(t_paired, 2, f)

##

beta_unpaired = apply(t_unpaired, 2, f)

##

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# Type 1 error 

mu = mu0

t_paired_er = t_unpaired_er = 
  ar_stat_er = matrix(0, M, m)

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
    
    y = X[, 1]
    x = X[, 2]
    
    t_paired_er[j, i] = t.test(x, 
                               y, 
                               paired = T)$statistic
    
    t_unpaired_er[j, i] = t.test(x, 
                                 y, 
                                 var.equal = T)$statistic
    
    ar_stat_er[j, i] = ar_compare_group_means(X)$rho
    
  }
  
  cat("\r", i)
  
}

f = function(x) mean(x > qu | x < ql)

error_paired = apply(t_paired_er, 2, f)

##

error_unpaired = apply(t_unpaired_er, 2, f)

##

error_ar = colMeans(ar_stat_er < qar)
##

Method = rep(c("ar_stat",
               "t_paired", "t_equalvar"), 
             each = length(r))

beta_v = c(beta_ar,
           beta_paired, beta_unpaired)

error_v = c(error_ar,
            error_paired, error_unpaired)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  r = rep(r, 3))

d = gsub("\\.", "", as.character(d))

f_path = paste0("simulations/results/", 
                pdist, "_", n, "n_", d, "d.txt")

write.table(Data, f_path, row.names = FALSE)
