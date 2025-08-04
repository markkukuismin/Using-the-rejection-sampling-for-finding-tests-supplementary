
# Talla skriptilla tehty kuva on akatemian hakemusta varten
# Nakevatpa millainen teho menetelmalla on!

library(ggplot2)

ar_dist = function(x){
  
  n = length(x)
  
  m = mean(x)
  
  s = sd(x)
  
  fhat = kdevine::kde1d(x)
  
  fhat = kdevine::dkde1d(x, fhat)
  
  rhox = dnorm(x, mean = m, sd = s)/fhat
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rho = (sum(a) + b)/n
  
  rho
  
}

set.seed(1)

n = 52

# mu under null

df_alt = 1:20

M = 10^3

ks_stat = ar_stat = 
  lillie_stat = matrix(0, M, length(df_alt))

for(i in 1:length(df_alt)){
  
  df = df_alt[i]
  
  for(j in 1:M){
    
    x = rt(n, df = df)
    
    ks_stat[j, i] = ks.test(x, "pnorm")$p.value
    
    lillie_stat[j, i] = nortest::lillie.test(x)$p.value
    
    ar_stat[j, i] = ar_dist(x)
    
  }
  
  cat("\r", i)
  
}

# ar null

ar_null = rep(0, M)

for(i in 1:M){
  
  x = rnorm(n)
  
  ar_null[i] = ar_dist(x)
  
}

f = function(x) mean(x < 0.05)

beta_ks = apply(ks_stat, 2, f)

beta_lillie = apply(lillie_stat, 2, f)

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# Type 1 error 

mu_alt = seq(-1, 1, length.out = length(df_alt))

ks_er = ar_stat_er = 
  lillie_er = matrix(0, M, length(mu_alt))

for(i in 1:length(mu_alt)){
  
  mu = mu_alt[i]
  
  for(j in 1:M){
    
    x = rnorm(n, mean = mu)
    
    ks_er[j, i] = ks.test(x, "pnorm", mu, 1)$p.value
    
    lillie_er[j, i] = nortest::lillie.test(x)$p.value
    
    ar_stat_er[j, i] = ar_dist(x)
    
    
  }
  
  cat("\r", i)
  
}

error_ks = apply(ks_er, 2, f)

error_lillie = apply(lillie_er, 2, f)

error_ar = colMeans(ar_stat_er < qar)

##

Method = rep(c("ar",
               "ks",
               "lillie"), 
             each = length(df_alt))

beta_v = c(beta_ar, beta_ks, beta_lillie)

error_v = c(error_ar, error_ks, error_lillie)

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  m = rep(mu_alt, 3),
                  df = df_alt)

ggplot(data = Data, 
       aes(df, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Degrees of freedom", 
       y = "Statistical power") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

ggplot(data = Data, 
       aes(m, error, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population mean", 
       y = "Type I error") +
  ylim(c(0, 1)) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))
