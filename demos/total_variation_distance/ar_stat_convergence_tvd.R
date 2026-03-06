
# AR stat and total variation difference (TVD)

# The AR stat converges to 1 - ||f - f_0||_TV where f is the population distribution
# (where the sample of data is actually drawn from), f_0 the given probability distribution
# (H0 distribution), and ||f - f_0||_TV = 0.5*\int |f(x) - f_0(x)| dx is the total
# variation distance when n -> infty

library(Rfast2)
library(LaplacesDemon)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

set.seed(1)

# Example of TVD when f is normal, g is uniform

f = function(x) dnorm(x, mean = 0, sd = 1)
g = function(x) dunif(x, min = -4, max = 4)

tvd = integrate(function(x) abs(f(x) - g(x)), lower = -Inf, upper = Inf)
0.5*tvd$value # TVD

# The value the AR stat converges when the population
# distribution is f

1 - 0.5*tvd$value

f0 = function(x) dnorm(x, mean = 0, sd = 1) # H0
f = function(x) dt(x, df = 6) # True pop. distribution where observations are drawn
#f = function(x) dlaplace(x)
#f = function(x) dunif(x, min= -4, max = 4)
#f = function(x) dnorm(x, mean = 0.5, sd = 1.2)

n = c(100, 500, 1000, 1500, 2000)

M = 100

rho = rho_true_f = matrix(0, nrow = M, ncol = length(n))

for(j in 1:length(n)){
  
  for(i in 1:M){
    
    #x = runif(n[j], min = -4, max = 4)
    x = rt(n[j], df = 6)
    #x = rlaplace(n[j])
    #x = rnorm(n[j], mean = 0.5, sd = 1.2)
    
    rho_temp = f0(x)/Rfast2::kernel(x)
    
    rho[i, j] = (sum(rho_temp >= 1) + sum(rho_temp[rho_temp < 1]))/n[j]
    
    rho_temp2 = f0(x)/f(x)
    
    rho_true_f[i, j] = (sum(rho_temp2 >= 1) + sum(rho_temp2[rho_temp2 < 1]))/n[j]
    
  }
  
  cat("\r", round(100*j/length(n), 2))
  
}

tvd = integrate(function(x) abs(f0(x) - f(x)), lower = -Inf, upper = Inf)
tvd = 0.5*tvd$value

# The value the AR stat converges when the population
# distribution is g

1 - tvd # 1 - TVD
colMeans(rho)
colMeans(rho_true_f)

df = data.frame(rho = c(rho),
                n = rep(n, each = M),
                f = rep("KDE", each = M))

df2 = data.frame(rho = c(rho_true_f),
                 n = rep(n, each = M),
                 f = rep("Population f", each = M))

df = rbind(df, df2)

df$f = as.factor(df$f)

df <- df %>%
  group_by(f, n) %>%
  summarise(
    sd = sd(rho, na.rm = TRUE),
    rho = mean(rho)
  )

p = ggplot(df, aes(x=n, y=rho, group = f, color = f)) + 
  geom_line(linewidth = 2) +
  geom_point()+
  geom_errorbar(aes(ymin = rho-qnorm(0.975)*sd/sqrt(M), ymax = rho+qnorm(0.975)*sd/sqrt(M)), 
                linewidth = 2) +
  geom_hline(yintercept = 1 - tvd, 
             linetype = "dashed", 
             linewidth = 2) +
  labs(color = "Density", 
       y = expression(hat(rho)),
       x = expression(n)) +
  theme(text=element_text(size=21))

# The difference between the blue and red line is explained by how good the
# kernel density estimate is

p

# Under null,

f0 = function(x) dnorm(x, mean = 0, sd = 1) # H0
f = f0

for(j in 1:length(n)){
  
  for(i in 1:M){
    
    x = rnorm(n[j])
    
    rho_temp = f0(x)/Rfast2::kernel(x)
    
    rho[i, j] = (sum(rho_temp >= 1) + sum(rho_temp[rho_temp < 1]))/n[j]
    
    rho_temp2 = f0(x)/f(x)
    
    rho_true_f[i, j] = (sum(rho_temp2 >= 1) + sum(rho_temp2[rho_temp2 < 1]))/n[j]
    
  }
  
  cat("\r", round(100*j/length(n), 2))
  
}

tvd_null = 0

df_null = data.frame(rho = c(rho),
                     n = rep(n, each = M),
                     f = rep("KDE", each = M))

df2_null = data.frame(rho = c(rho_true_f),
                      n = rep(n, each = M),
                      f = rep("Population f", each = M))

df_null = rbind(df_null, df2_null)

df_null$f = as.factor(df_null$f)

df_null <- df_null %>%
  group_by(f, n) %>%
  summarise(
    sd = sd(rho, na.rm = TRUE),
    rho = mean(rho)
  )

p_null = ggplot(df_null, aes(x=n, y=rho, group = f, color = f)) + 
  geom_line(linewidth = 2) +
  geom_point()+
  geom_errorbar(aes(ymin = rho-qnorm(0.975)*sd/sqrt(M), ymax = rho+qnorm(0.975)*sd/sqrt(M)), 
                linewidth = 2) +
  geom_hline(yintercept = 1 - tvd_null, 
             linetype = "dashed", 
             linewidth = 2) +
  labs(color = "Density", 
       y = expression(hat(rho)),
       x = expression(n)) +
  theme(text=element_text(size=21))

p_null

##

fn_null = paste0("demos/total_variation_distance/tvd_null.png")

ggsave(fn_null,
       plot = p_null,
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)

df$hyp = "Alt"
df_null$hyp = "Null"

df$hyp = as.factor(df$hyp)
df_null$hyp = as.factor(df_null$hyp)

df_all = rbind(df, df_null)

data_hline = data.frame(hyp = c("Alt", "Null"),
                        hline = c(1 - tvd, 1))

data_hline$hyp = as.factor(data_hline$hyp)

p_all = ggplot(df_all, aes(x=n, y=rho, color = f)) + 
  geom_line(linewidth = 2) +
  geom_point()+
  geom_errorbar(aes(ymin = rho-qnorm(0.975)*sd/sqrt(M), ymax = rho+qnorm(0.975)*sd/sqrt(M)), 
                linewidth = 2) +
  labs(color = "Density", 
       y = expression(hat(rho)(X)),
       x = expression(n)) +
  theme(text=element_text(size=21)) +
  scale_x_continuous(
    limits = c(min(n), NA),
    breaks = n
  )

hyp_edit = c(`Alt` = "H[0]~is~false",
             `Null` = "H[0]~is~true")

p_all = p_all + 
  facet_grid(~hyp, 
             labeller = labeller(hyp = as_labeller(hyp_edit, label_parsed)))

p_all

p_all = p_all +
  geom_hline(data = data_hline,
             aes(yintercept = hline), 
             linetype = "dashed", 
             linewidth = 2)

p_all

fn_both = paste0("demos/total_variation_distance/tvd_both.png")

ggsave(fn_both,
       plot = p_all,
       device = png, 
       width = 400, 
       height = 250, 
       units = "mm",
       dpi = 700)
