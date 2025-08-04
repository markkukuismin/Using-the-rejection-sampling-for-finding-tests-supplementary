
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

n = 20 # 20, 30, 50

# Power (H_0 not true)

pdf = "dlst" # dlst

f0 = stringr::str_replace(pdf, "d", "r")

f0 = get(f0)

kspdf = stringr::str_replace(pdf, "d", "p")

mu_alt = seq(-2, 2, length.out = 30)

ref_param = 0

mu_alt = c(mu_alt, ref_param)
mu_alt = sort(mu_alt)

M = 10^4

ks_stat = ar_stat = 
  cvm_stat = matrix(0, M, length(mu_alt))

for(i in 1:length(mu_alt)){
  
  r = mu_alt[i]
  
  for(j in 1:M){
    
    x = f0(n, mu = r, sigma = 1, df = 3)
    
    ks_stat[j, i] = ks.test(x, 
                            kspdf, 
                            mu = ref_param,
                            sigma = 1,
                            df = 3)$p.value
    
    cvm_stat[j, i] = cvm.test(x, 
                              kspdf, 
                              mu = ref_param,
                              sigma = 1,
                              df = 3)$p.value
    
    ar_stat[j, i] = ar_dist(x, 
                            f0 = pdf, 
                            mu = ref_param,
                            df = 3,
                            sigma = 1)
    
  }
  
  cat("\r", i)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  x = f0(n, mu = 0, sigma = 1, df = 3)
  
  ar_null[i] = ar_dist(x, 
                       f0 = pdf, 
                       mu = 0, 
                       sigma = 1,
                       df = 3)
  
}

f = function(x) mean(x < 0.05)

beta_ks = apply(ks_stat, 2, f)

beta_cvm = apply(cvm_stat, 2, f)

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# Type 1 error 

ks_er = ar_stat_er =
  cvm_stat_er = matrix(0, M, length(mu_alt))

for(i in 1:length(mu_alt)){
  
  r = mu_alt[i]
  
  for(j in 1:M){
    
    x = f0(n = n, mu = r, sigma = 1, df = 3)
    
    ks_er[j, i] = ks.test(x, 
                          kspdf, 
                          mu = r, 
                          sigma = 1,
                          df = 3)$p.value
    
    cvm_stat_er[j, i] = cvm.test(x, 
                                 kspdf, 
                                 mu = r,
                                 sigma = 1,
                                 df = 3)$p.value
    
    ar_stat_er[j, i] = ar_dist(x,
                               f0 = pdf,
                               mu = r,
                               sigma = 1,
                               df = 3)
    
  }
  
  cat("\r", i)
  
}

error_ks = apply(ks_er, 2, f)

error_cvm = apply(cvm_stat_er, 2, f)

error_ar = colMeans(ar_stat_er < qar)

##

Method = rep(c("ar",
               "ks",
               "cvm"), 
             each = length(mu_alt))

beta_v = c(beta_ar, beta_ks, beta_cvm)

error_v = c(error_ar, error_ks, error_cvm)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  scale = rep(mu_alt, 3))

ggplot(data = Data, 
       aes(scale, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Scale", 
       y = "Statistical power") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = ref_param) + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

ggplot(data = Data, 
       aes(scale, error, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Scale", 
       y = "Type I error") +
  ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

f_path = paste0("simulations/test_any_distribution/results/", 
                pdf, "_location_", n, "n", ".txt")

write.table(Data, f_path, row.names = FALSE)
