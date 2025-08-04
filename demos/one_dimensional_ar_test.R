
# The one dimensional t-test build using the probability of acceptance as 
# a test statistics is the same as the paired Student's test,

library(ggplot2)
library(gridExtra)

ar_paired = function(x, y, mu0 = 0){
  
  d = x - y
  
  n = sum(!is.na(d))
  
  xd = mean(d)
  
  s = sd(d)
  
  xd = (xd - mu0)/(s/sqrt(n))
  
  df = n - 1
  
  rho = dnorm(xd)/dt(xd, df = df)
  
  rho = ifelse(rho >= 1, 1, rho)
  
}

t_test_paired = function(x, y, mu0 = 0){
  
  d = x - y
  
  n = sum(!is.na(d))
  
  xd = mean(d)
  
  s = sd(d)
  
  xd = (xd - mu0)/(s/sqrt(n))
  
  xd
  
}

set.seed(1)

r = seq(-0.99, 0.99, length.out = 40)

r = c(r, 0)

r = unique(r)

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

# power

m = length(r)

M = 10^4

t_paired = ar_stat = w_stat = matrix(0, M, m)

for(i in 1:m){
  
  Sigma = matrix(r[i], 2, 2)
  
  diag(Sigma) = 1
  
  for(j in 1:M){
    
    X = MASS::mvrnorm(n, mu = mu, Sigma = Sigma) 
    
    y = X[, 1]
    x = X[, 2]
    
    t_paired[j, i] = t_test_paired(x, y)
    w_stat[j, i] = wilcox.test(x, y, paired = TRUE)$p.value
    ar_stat[j, i] = ar_paired(x, y)
    
    
  }
  
  cat("\r", i)
  
}

# null distribution for the ar-statistic

ar_null = rep(0, M)

mu0 = rep(0, 2)

I = diag(1, 2)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = mu0, Sigma = I)
  
  x = X[, 1]
  y = X[, 2]
  
  Sigma_hat = cor(X)
  
  ar_null[j] = ar_paired(x, y)
  
}

ql = qt(0.025, df = n - 1)
qu = qt(0.975, df = n - 1)

f = function(x) mean(x > qu | x < ql)

beta_paired = apply(t_paired, 2, f)

##

beta_w = apply(w_stat, 2, function(x) mean(x < 0.05))

##

qar = quantile(ar_null, 0.05)
beta_ar = colMeans(ar_stat < qar)

# Type 1 error 

mu = c(0, 0)

t_paired_er = w_stat_er = ar_stat_er = matrix(0, M, m)

for(i in 1:m){
  
  Sigma = matrix(r[i], 2, 2)
  
  diag(Sigma) = 1
  
  for(j in 1:M){
    
    X = MASS::mvrnorm(n, mu = mu, Sigma = Sigma) 
    
    y = X[, 1]
    x = X[, 2]
    
    t_paired_er[j, i] = t_test_paired(x, y)
    w_stat_er[j, i] = wilcox.test(x, y, paired = TRUE)$p.value
    ar_stat_er[j, i] = ar_paired(x, y)
    
    
  }
  
  cat("\r", i)
  
}

error_paired = apply(t_paired_er, 2, f)

##

error_w = apply(w_stat_er, 2, function(x) mean(x < 0.05))

##

error_ar = colMeans(ar_stat_er < qar)

##

# Likelihood ratio test power,

b = 2*(1 - r)/n

eta = d

cdot = sqrt(b)*qnorm(0.95)

ef = (1/sqrt(b))*(cdot - eta)

lr_power = 1 - pnorm(ef)

##

Method = rep(c("AR test","Paired t-test", "Wilcox", "LR"), 
             each = length(r))

beta_v = c(beta_ar, beta_paired, beta_w, lr_power)

error_v = c(error_ar, error_paired, error_w, rep(0.05, m))

Data = data.frame(Method = Method,
                  power = beta_v,
                  error = error_v,
                  r = rep(r, 4))

ggplot(data = Data, 
       aes(r, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Statistical power") + 
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

p1 = ggplot(data = Data, 
            aes(r, power, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Statistical power") +
  ylim(c(0, 1)) +
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

p2 = ggplot(data = Data, 
       aes(r, error, color = Method)) +
  geom_line(aes(linetype = Method),
            linewidth = 1.5) +
  labs(x = "Population correlation", 
       y = "Type I error") +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.05,
             linetype = "dashed") +
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size = 16))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 = p1 + theme(legend.position = "bottom") + ggplot2::guides(shape = guide_legend(nrow = 1))

p1

legend = get_legend(p1)

p1 = p1 + theme(legend.position="none")
p2 = p2 + theme(legend.position="none")

final_plot = grid.arrange(arrangeGrob(p1, p2,
                                      legend, 
                                      ncol = 2, 
                                      nrow = 2, 
                                      layout_matrix = matrix(c(1, 2,
                                                               3, 3),
                                                             ncol = 2, byrow = T),
                                      widths = c(2.7, 2.7), heights = c(2.5, 0.2))
)

ggsave("demos/paired_t_power_and_error.png",
       plot = final_plot,
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)
