
library(ggplot2)

source("functions/ar_compare_group_means.R")

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

t_unpaired = ar_stat = matrix(0, M, m)

for(i in 1:m){
  
  mu = c(0, eta[i])
  
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
    
    t_unpaired[j, i] = t.test(x, y, var.equal = T)$statistic
    ar_stat[j, i] = ar_compare_group_means(X)$rho
    
  }
  
  cat("\r", i)
  
}

# null distribution for the ar-statistic

ar_null = rep(0, M)

mu0 = rep(0, 2)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = mu0, Sigma = I)
  
  x = X[, 1]
  y = X[, 2]
  
  Sigma_hat = cor(X)
  
  ar_null[j] = ar_compare_group_means(X)$rho
  
}

ql = qt(0.025, df = n - 1)
qu = qt(0.975, df = n - 1)

f = function(x) mean(x > qu | x < ql)

beta_unpaired = apply(t_unpaired, 2, f)

##

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# LR power,

cdot = sqrt(2/n)*qnorm(0.95)

ef = sqrt(n/2)*(cdot - abs(eta))

lr_power = 1 - pnorm(ef)

##

Method = rep(c("ar_stat", 
               "t_equalvar",
               "LR"), 
             each = m)

beta_v = c(beta_ar, 
           beta_unpaired, lr_power)

Data = data.frame(Method = Method,
                  power = beta_v,
                  eta = rep(eta, 3))

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


f_path = paste0("simulations/results/", 
                pdist, "_effect_vs_power_", n, "n", ".txt")

write.table(Data, f_path, row.names = FALSE)
