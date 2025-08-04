
source("functions/ar_paired_mvt.R")

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

# 0.4, 0.6, 0.8 e.g., mu = c(0, 0.6)

mu = c(0, 0.4)

d = abs(diff(mu))

# The distribution: norm, t, pexp, unif, laplace

pdist = "norm"

if(pdist == "t") df = 3 # for the t-distribution

# power

m = length(r)

M = 10^4

t_paired = t_unpaired = matrix(0, M, m)

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
    
    t_paired[j, i] = t.test(x, y, 
                            paired = T,
                            alternative = "greater")$statistic
    t_unpaired[j, i] = t.test(x, y,
                              alternative = "greater",
                              var.equal = T)$statistic
    
  }
  
  cat("\r", i)
  
}

ql = qt(0.025, df = n - 1)
qu = qt(0.975, df = n - 1)

f = function(x) mean(x > qu | x < ql)

beta_paired = apply(t_paired, 2, f)

##

beta_unpaired = apply(t_unpaired, 2, f)

# Type 1 error 

mu = c(0, 0)

t_paired_er = t_unpaired_er = 
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
    
    y = X[, 1]
    x = X[, 2]
    
    t_paired_er[j, i] = t.test(x, 
                               y, 
                               paired = T,
                               alternative = "greater")$statistic
    
    t_unpaired_er[j, i] = t.test(x, 
                                 y, 
                                 var.equal = T,
                                 alternative = "greater")$statistic
    
  }
  
  cat("\r", i)
  
}

f = function(x) mean(x > qu | x < ql)

error_paired = apply(t_paired_er, 2, f)

##

error_unpaired = apply(t_unpaired_er, 2, f)

##

Method = rep(c("t_paired", "t_equalvar"), 
             each = length(r))

beta_v = c(beta_paired, beta_unpaired)

error_v = c(error_paired, error_unpaired)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  r = rep(r, 2))

# LR

b = 2*(1 - r)/n

cdot = sqrt(b)*qnorm(0.95)

eta = tail(d, 1)

ef = (1/sqrt(b))*(cdot - eta)

lr_power = 1 - pnorm(ef)

Data_temp = data.frame(Method = "LR",
                       power = lr_power,
                       error = 0.05,
                       r = r)

Data = rbind(Data, Data_temp)

p = ggplot(data = Data, 
           aes(r, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Statistical power",
       title = paste0("Difference between pop. means = ", tail(d, 1))) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

p  
