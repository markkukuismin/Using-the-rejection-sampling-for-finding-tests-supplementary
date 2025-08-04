
library(ggplot2)
library(mvhtests)

source("functions/ar_two_groups.R")
source("functions/MLE0.R")

set.seed(1)

pwr::cohen.ES(test = "t", size = "large")
pwr::cohen.ES(test = "t", size = "medium")
pwr::cohen.ES(test = "t", size = "small")

pwr::pwr.t.test(d = 0.8,
                sig.level = 0.05,
                power = 0.8,
                type = "two.sample",
                alternative = "two.sided")

pwr::pwr.t.test(d = 0.5,
                sig.level = 0.05,
                power = 0.8,
                type = "two.sample",
                alternative = "two.sided")

pwr::pwr.t.test(d = 0.2,
                sig.level = 0.05,
                power = 0.8,
                type = "two.sample",
                alternative = "two.sided")

n = 26 # 26, 64, 394

# Different effect sizes,

mu0 = rep(0, 2)

eta = seq(-1, 1, length.out = 50)

eta = c(eta, 0)
eta = sort(eta)
eta = unique(eta)

m = length(eta)

# The distribution: norm, t, pexp, unif, laplace

pdist = "norm"

if(pdist == "t") df = 3 # for the t-distribution

# power

I = diag(1, 2)

M = 10^4

lr_stat = el_stat = ar_stat = 
  ar_stat_hat = matrix(0, M, m)

for(i in 1:m){
  
  mu = mu0 + c(0, eta[i])
  
  Sigma = I
  
  if(pdist == "t"){
    Sigma = Sigma*((df-2)/df)
  }
  
  mu = c(0, eta[i])
  
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

mu0 = rep(0, 2)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = mu0, Sigma = I)
  
  x = X[, 1]
  y = X[, 2]
  
  Sigma_hat = cor(X)
  
  ar_null[j] = ar_two_groups(x, y)$rho
  ar_null_hat[j] = ar_two_groups(x, y, r = Sigma_hat)$rho
  
}

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

##

Method = rep(c("ar_sample_cor", 
               "ar_pop_cor",
               "LR",
               "EL"), 
             each = m)

beta_v = c(beta_ar_hat, beta_ar, 
           beta_lr, beta_el)

Data = data.frame(Method = Method,
                  power = beta_v,
                  eta = rep(eta, 4))

ggplot(data = Data, 
       aes(eta, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Difference between means", 
       y = "Statistical power") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"),
        legend.position = "bottom") + 
  guides(shape = guide_legend(nrow = 1))


f_path = paste0("simulations/ar_vs_lr_vs_el/results/", 
                pdist, "_effect_vs_power_", n, "n", ".txt")

write.table(Data, f_path, row.names = FALSE)
