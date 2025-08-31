
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

n = 50 # 20, 30, 50

# Power (H_0 not true)

# Log-normal (dlnormp) needs still work...

pdf = "dlst" # dlaplace, dcauchy, dlnormp, dlst

f0 = stringr::str_replace(pdf, "d", "r")

f0 = get(f0)

kspdf = stringr::str_replace(pdf, "d", "p")

scale_alt = seq(0.1, 5, length.out = 30)

ref_param = 2.5

scale_alt = c(scale_alt, ref_param)
scale_alt = sort(scale_alt)

if(pdf == "dlnormp"){
  
  scale_alt = seq(0.1, 2, length.out = 30)
  
  ref_param = 1
  
  scale_alt = c(scale_alt, ref_param)
  scale_alt = sort(scale_alt)
  
}

#M = 5*10^2
M = 10^4

ks_stat = ar_stat = 
  cvm_stat = ad_stat = matrix(0, M, length(scale_alt))

for(i in 1:length(scale_alt)){
  
  r = scale_alt[i]
  
  for(j in 1:M){
    
    if(!(pdf %in% c("dlnormp", "dlst"))){
      
      x = f0(n, scale = r)
      
      ks_stat[j, i] = ks.test(x, kspdf, scale = ref_param)$p.value
      
      cvm_stat[j, i] = cvm.test(x, 
                                kspdf, 
                                scale = ref_param)$p.value
      
      ar_stat[j, i] = ar_dist(x, 
                              f0 = pdf, 
                              scale = ref_param)
      
    }
    
    if(pdf == "dlnormp"){
      
      x = f0(n, mu = 1, tau = r)
      
      ks_stat[j, i] = ks.test(x, 
                              kspdf, 
                              mu = 1,
                              tau = ref_param)$p.value
      
      cvm_stat[j, i] = cvm.test(x, 
                                kspdf, 
                                scale = ref_param)$p.value
      
      ar_stat[j, i] = ar_dist(x, 
                              f0 = pdf, 
                              mu = 1,
                              tau = ref_param)
      
    }
    
    if(pdf == "dlst"){
      
      x = f0(n, mu = 0, sigma = r, df = 3)
      
      ks_stat[j, i] = ks.test(x, 
                              kspdf, 
                              mu = 0,
                              sigma = ref_param,
                              df = 3)$p.value
      
      cvm_stat[j, i] = cvm.test(x, 
                                "plst", 
                                mu = 0,
                                sigma = ref_param,
                                df = 3)$p.value
      
      ar_stat[j, i] = ar_dist(x, 
                              f0 = pdf, 
                              mu = 0,
                              df = 3,
                              sigma = ref_param)
      
      ad_stat[j, i] = goftest::ad.test(x,
                                    null = "plst",
                                    mu = 0,
                                    df = 3,
                                    sigma = ref_param)$p.value
      
      
    }
    
  }
  
  cat("\r", i)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  if(!(pdf %in% c("dlnormp", "dlst"))){
    
    x = f0(n)
    
    ar_null[i] = ar_dist(x, f0 = pdf)
    
  } 
  
  if(pdf == "dlnormp"){
    
    x = f0(n, mu = 1, tau = 1)
    
    ar_null[i] = ar_dist(x, f0 = pdf, mu = 1, tau = 1)
    
  }
  
  if(pdf == "dlst"){
    
    x = f0(n, mu = 0, sigma = 1, df = 3)
    
    ar_null[i] = ar_dist(x, 
                         f0 = pdf, 
                         mu = 0, 
                         sigma = 1,
                         df = 3)
    
  }
  
}

f = function(x) mean(x < 0.05)

beta_ks = apply(ks_stat, 2, f)

beta_cvm = apply(cvm_stat, 2, f)

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

beta_ad = apply(ad_stat, 2, f)

# Type 1 error 

ks_er = ar_stat_er =
  cvm_stat_er = ad_stat_er = matrix(0, M, length(scale_alt))

for(i in 1:length(scale_alt)){
  
  r = scale_alt[i]
  
  for(j in 1:M){
    
    if(!(pdf %in% c("dlnormp", "dlst"))){
      
      x = f0(n = n, scale = r)
      
      ks_er[j, i] = ks.test(x, kspdf, scale = r)$p.value
      
      cvm_stat_er[j, i] = cvm.test(x, 
                                   kspdf, 
                                   scale = r)$p.value
      
      
      ar_stat_er[j, i] = ar_dist(x,
                                 f0 = pdf,
                                 scale = r)
      
    }
    
    if(pdf == "dlnormp"){
      
      x = f0(n = n, mu = 1, tau = r)
      
      ks_er[j, i] = ks.test(x, 
                            kspdf, 
                            mu = 1, 
                            tau = r)$p.value
      
      cvm_stat_er[j, i] = cvm.test(x, 
                                   kspdf, 
                                   scale = r)$p.value
      
      ar_stat_er[j, i] = ar_dist(x,
                                 f0 = pdf,
                                 mu = 1,
                                 tau = r)
      
    }
    
    if(pdf == "dlst"){
      
      x = f0(n = n, mu = 0, sigma = r, df = 3)
      
      ks_er[j, i] = ks.test(x, 
                            kspdf, 
                            mu = 0, 
                            sigma = r,
                            df = 3)$p.value
      
      cvm_stat_er[j, i] = cvm.test(x, 
                                   kspdf, 
                                   mu = 0,
                                   sigma = r,
                                   df = 3)$p.value
      
      ar_stat_er[j, i] = ar_dist(x,
                                 f0 = pdf,
                                 mu = 0,
                                 sigma = r,
                                 df = 3)
      
      ad_stat_er[j, i] = goftest::ad.test(x,
                                          null = "plst",
                                          mu = 0,
                                          df = 3,
                                          sigma = r)$p.value
      
      
    }
    
  }
  
  cat("\r", i)
  
}

error_ks = apply(ks_er, 2, f)

error_cvm = apply(cvm_stat_er, 2, f)

error_ar = colMeans(ar_stat_er < qar)

error_ad = apply(ad_stat_er, 2, f)

##

Method = rep(c("ar",
               "ks",
               "cvm",
               "ad"), 
             each = length(scale_alt))

beta_v = c(beta_ar, beta_ks, beta_cvm, beta_ad)

error_v = c(error_ar, error_ks, error_cvm, error_ad)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  scale = rep(scale_alt, 4))

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
                pdf, "_", n, "n", ".txt")

write.table(Data, f_path, row.names = FALSE)
