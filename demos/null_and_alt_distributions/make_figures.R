
library(ggplot2)
library(gridExtra)

source("functions/ar_compare_group_means.R")

set.seed(1)

r = 0.8

Sigma = matrix(r, 2, 2)
diag(Sigma) = 1

mu = c(0, 0.2)

N = c(20, 50, 100)

M = 10^4

pdist = "norm"

i = 1

for(n in N){
  
  X = MASS::mvrnorm(n = n,
                    mu = mu,
                    Sigma = Sigma)
  
  p = ncol(X)
  
  R = cov(X)
  
  mu0 = rep(mean(X), p)
  Y = t(X) - mu0
  R0 = tcrossprod(Y)/(n - 1)
  
  mu0 = rep(mean(X), p)
  mu1 = colMeans(X)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(mu1, mean = mu0, sigma = R0)
  b1 = mvtnorm::dmvt(mu1, df = df, delta = mu1, sigma = R, log = FALSE)
  
  rho = b0/b1
  
  Ru = matrix(runif(n*M), n, M)
  
  tempU = apply(Ru, 2, function(x) rho > x)
  
  rhod = colMeans(tempU)
  
  rhod = ifelse(rhod >= 1, 1, rhod)
  
  rho = ifelse(rho >= 1, 1, rho)
  
  SErho = sqrt(rho*(1 - rho)/n)
  
  xx = seq(0, 1, length.out = M)
  rho_norm = dnorm(xx, mean = rho, sd = SErho)
  
  a = rho*(n - 1)
  b = (1 - rho)*(n - 1)
  
  rho_beta = dbeta(xx, shape1 = a, shape2 = b)
 
  Data = data.frame(rho_beta = rho_beta,
                    rho_norm = rho_norm,
                    x = xx,
                    rhod = rhod,
                    n = n)
  
  D = paste0("Data", i)
  
  Data$n = as.factor(Data$n)
  
  assign(D, Data)
   
  i = i + 1
  
}

p1 = ggplot(Data1, aes(x = rhod)) +
  geom_histogram(aes(y = after_stat(density)), 
                 fill = "white",
                 color = "black",
                 bins = 13) +
  labs(x = "AR statistic", 
       y = "Density",
       tag = "(A)") +
  geom_line(aes(x = x, y = rho_norm), 
            linewidth = 1.5,
            linetype = "dashed") +
  geom_line(aes(x = x, y = rho_beta), 
            linewidth = 1.5) +
  xlim(c(0.4, 1)) + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

p2 = ggplot(Data2, aes(x = rhod)) +
  geom_histogram(aes(y = after_stat(density)), 
                 fill = "white",
                 color = "black",
                 bins = 15) +
  labs(x = "AR statistic", 
       y = "Density",
       tag = "(B)") +
  geom_line(aes(x = x, y = rho_norm), 
            linewidth = 1.5,
            linetype = "dashed") +
  xlim(c(0.7, 1)) +
  geom_line(aes(x = x, y = rho_beta), linewidth = 1.5) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

p3 = ggplot(Data3, aes(x = rhod)) +
  geom_histogram(aes(y = after_stat(density)), 
                 fill = "white",
                 color = "black",
                 bins = 15) +
  labs(x = "AR statistic", 
       y = "Density",
       tag = "(C)") +
  geom_line(aes(x = x, y = rho_norm), 
            linewidth = 1.5,
            linetype = "dashed") +
  xlim(c(0.7, 1)) +
  geom_line(aes(x = x, y = rho_beta), linewidth = 1.5) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.tag = element_text(size = 16,
                                face = "bold"))

final_plot = grid.arrange(arrangeGrob(p1, p2, p3,
                                      ncol = 3, 
                                      nrow = 1, 
                                      layout_matrix = matrix(c(1, 2, 3),
                                                             ncol = 3, byrow = T),
                                      widths = c(2.7, 2.7, 2.7), heights = c(2.5))
)


fn = paste0("demos/null_and_alt_distributions/figures/test_stat_", pdist, ".png")

ggsave(fn,
       plot = final_plot,
       device = png, 
       width = 400, 
       height = 200, 
       units = "mm",
       dpi = 700)
