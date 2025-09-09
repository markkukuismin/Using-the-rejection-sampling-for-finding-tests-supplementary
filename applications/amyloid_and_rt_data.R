
library(openintro)
library(Stat2Data)
library(rtdists)
library(tidyverse)
library(brms)

source("applications/ar_group_means_mc_pvalues.R")

set.seed(1)

# Student's Sleep Data

data("sleep")

n = sum(sleep$group == 1)

D = matrix(sleep$extra, nrow = n, 2)

res = ar_group_means_mc_pvalues(D, 
                                M = 1000,
                                paired = TRUE)

res$rho_m
res$p_value

d = poibin::dpoibin(0:n, pp = rep(res$rho_m, n))

plot(0:n/n,
     d,
     xlim = c(0, 1),
     type = "h",
     lwd = 2,
     ylab = "Probability",
     xlab = expression(rho(X)))

abline(v = res$rho_m, 
       lty = 2, 
       lwd = 2,
       col = "red")

# M = 10^4
# 
# U = matrix(runif(n*M), n, M)
# 
# f = function(x) res$rho > x
# 
# I = apply(U, 2, f)
# 
# dd = colMeans(I)
# 
# dd = table(dd)
# 
# x = as.numeric(names(dd))
# 
# lines(x + 0.01,
#       dd/M,
#       lty = 2,
#       lwd = 2,
#       col = "blue",
#       type = "h")

##

# Amyloid-beta and Cognitive Impairment

data("Amyloid")

table(Amyloid$Group)

n = 21

D = matrix(0, n, 3)

D[, 1] = c(Amyloid$Abeta[Amyloid$Group == "mAD"],
           rep(NA, 4))
D[, 2] = Amyloid$Abeta[Amyloid$Group == "MCI"]
D[, 3] = c(Amyloid$Abeta[Amyloid$Group == "NCI"],
           rep(NA, 2))

colnames(D) = c("mAD", "MCI", "NCI") 

boxplot(Abeta ~ Group, data = Amyloid)

summary(aov(Abeta ~ Group, data = Amyloid))

res = ar_group_means_mc_pvalues(D, M = 1000)

res$rho_m
res$p_value

d = poibin::dpoibin(0:n, pp = rep(res$rho_m, n))

plot(0:n/n,
     d,
     xlim = c(res$rho_m - 0.3, res$rho_m + 0.3),
     type = "h",
     lwd = 2,
     ylab = "Probability",
     xlab = expression(rho(X)),
     col = "darkgray")

abline(v = res$rho_m, 
       lty = 2, 
       lwd = 2,
       col = "red")

xx = seq(0, 1, length.out = 1000)

m = res$rho_m
s = sqrt(m*(1 - m)/n)

lines(xx,
      dnorm(xx, mean = m, sd = s)/n,
      lwd = 2,
      lty = 2)

round(poibin::qpoibin(c(0.025, 0.975), 
                      pp = rep(res$rho_m, n))/n, 3)

qnorm(c(0.025, 0.975), 
      mean = m, 
      sd = s)

# Multiple testing correction

AD_vs_MCI = ar_group_means_mc_pvalues(D[,c(1, 2)], 
                                      M = 1000)

AD_vs_NCI = ar_group_means_mc_pvalues(D[,c(1, 3)], 
                                      M = 1000)

NCI_vs_MCI = ar_group_means_mc_pvalues(D[,c(2, 3)], 
                                       M = 1000)

round(AD_vs_NCI$rho_m, 3)
round(AD_vs_MCI$rho_m, 3)
round(NCI_vs_MCI$rho_m, 3)

round(poibin::qpoibin(c(0.025, 0.975), 
                      pp = rep(AD_vs_NCI$rho_m, n))/n, 3)

round(poibin::qpoibin(c(0.025, 0.975), 
                      pp = rep(AD_vs_MCI$rho_m, n))/n, 3)

round(poibin::qpoibin(c(0.025, 0.975), 
                      pp = rep(NCI_vs_MCI$rho_m, n))/n, 3)

round(AD_vs_NCI$p_value*3, 3)
round(AD_vs_MCI$p_value*3, 3)
round(NCI_vs_MCI$p_value*3, 3)

TukeyHSD(aov(Abeta ~ Group, data = Amyloid))

# M = 10^4
#
# U = matrix(runif(n*M), n, M)
# 
# f = function(x) res$rho > x
# 
# I = apply(U, 2, f)
# 
# dd = colMeans(I)
# 
# dd = table(dd)
# 
# x = as.numeric(names(dd))
# 
# lines(x + 0.01,
#       dd/M,
#       lty = 2,
#       lwd = 2,
#       col = "blue",
#       type = "h")

# Distribution of reaction times

source("applications/ar_mc_pvalue.R")

data(speed_acc)
head(speed_acc)

# remove excluded,

rt_data = speed_acc %>%
  filter(!censor)

rt_data = rt_data %>%
  filter(id == 1)

rt_x = unique(rt_data$rt)

#fit = brms::brm(rt ~ condition, 
#                rt_data, family=shifted_lognormal())

x = seq(0.1, 2.5, length.out = 1000)

hist(rt_x, 
     probability = TRUE,
     xlim = c(0.1, 2),
     breaks = 70)

res = EnvStats::elnorm3(rt_x)

theta = round(res$parameters, 2)

lines(x, 
      dshifted_lnorm(x, 
               meanlog = theta[1], 
               sdlog = theta[2],
               shift = theta[3]),
      lwd = 2)

plot(density(rt_x), lwd = 2)
lines(x, 
      dshifted_lnorm(x, 
                     meanlog = theta[1], 
                     sdlog = theta[2],
                     shift = theta[3]),
      lwd = 2,
      lty = 2,
      col = "red")

ks.test(rt_x, 
        "pshifted_lnorm", 
        meanlog = theta[1], 
        sdlog = theta[2],
        shift = theta[3])

res = ar_mc_pvalue(rt_x, 
                   M = 1000,
                   f0 = "dshifted_lnorm", 
                   meanlog = theta[1], 
                   sdlog = theta[2],
                   shift = theta[3])

rho = res$rho_m
rho

hist(res$rho_mc, probability = TRUE)
abline(v = rho, lwd = 2)

res$p_value

n = length(rt_x)

d = poibin::dpoibin(0:n, pp = res$rhos)

plot(0:n/n,
     d,
     xlim = c(rho - 0.05, rho + 0.05),
     type = "h",
     lwd = 2,
     ylab = "Probability",
     xlab = expression(rho(X)),
     col = "darkgray")

abline(v = rho, 
       lty = 2, 
       lwd = 2,
       col = "red")

# xx = seq(0, 1, length.out = 10000)
# 
# m = rho_m
# s = sqrt(m*(1 - m)/n)
# 
# lines(xx,
#       dnorm(xx, mean = m, sd = s)/n,
#       lwd = 2,
#       lty = 2)

round(poibin::qpoibin(c(0.025, 0.975), 
                      pp = rep(rho, n))/n, 3)

df = data.frame(x = x, y = dshifted_lnorm(x, 
                                          meanlog = theta[1], 
                                          sdlog = theta[2],
                                          shift = theta[3]))

p1 = ggplot(data = data.frame(RT = rt_x), aes(x = RT)) +
  geom_histogram(aes(y = after_stat(density))) +
  ylab("Density") +
  geom_line(data = df, 
            aes(x = x, y = y), 
            color = "red",
            linewidth = 1.7) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.tag = element_text(size = 18,
                                face = "bold"))

p1

ggsave("applications/RT_dist.png",
       plot = p1,
       device = png, 
       width = 400, 
       height = 200, 
       units = "mm",
       dpi = 700)

