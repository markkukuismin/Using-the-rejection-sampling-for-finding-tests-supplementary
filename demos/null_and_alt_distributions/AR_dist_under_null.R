
# In general, rho ~ N(E(rho), var(rho)/n) because
# I am using the mean as the test statistic

# E(rho) \approx 1 so the normal approximation is 
# no good. What about beta approximation:

# E(rho) \sim Beta(a, b) where a = r*(n - 1),
# b = (1 - r)*(n - 1) and r = E(rho)

source("functions/ar_compare_group_means.R")

set.seed(1)

pwr::cohen.ES(test = "t", size = "medium")

# effect, 0.2, 0.6, 0.8

pwr::pwr.t.test(d = 0.6, 
                power = 0.8,
                type = "paired",
                alternative = "two.sided")

n = 15 # 15, 24, 52
p = 2

mu = rnorm(p)
Sigma = diag(1, p)

pdist = "norm"

if(pdist == "t"){
  Sigma = Sigma*((df-2)/df)
}

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

M = 5*10^4

ar_null = rep(0, M)

mu0 = rep(0, p)

I = diag(1, p)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = mu0, Sigma = I)
  
  mu0 = rep(mean(X), p)
  #Y = t(X) - mu0
  #R = tcrossprod(Y)/(n - 1)
  
  R = cov(X)
  
  Xm = colMeans(X)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(t0, sigma = R)
  b1 = mvtnorm::dmvt(t0, df = df, sigma = R, log = FALSE)
  
  if(b0 == 0) rho = 0
  if(b0 != 0) rho = b0/b1
  
  rho = ifelse(rho >= 1, 1, rho)
  
  ar_null[j] = rho
  
  pr = round(100*j/M, 2)
  
  cat("\r", pr)
  
}

hist(ar_null, probability = TRUE)

ar_m = mean(ar_null)
ar_s = sd(ar_null)

xx = seq(0, 1, length.out = 1000)

lines(xx, 
      dnorm(xx, mean = ar_m, sd = ar_s), 
      lty = 2,
      lwd = 2)

abline(v = ar_m, lwd = 2)

HDInterval::hdi(ar_null, credMass = 0.9)
quantile(ar_null, c(0.05, 0.95))
qnorm(c(0.05, 0.95), mean = ar_m, sd = ar_s)

# Beta approx

hist(ar_null, probability = TRUE)

a = ar_m*(n - 1)
b = (1 - ar_m)*(n - 1)

lines(xx, 
      dbeta(xx, shape1 = a, shape2 = b), 
      lty = 2,
      lwd = 2)

plot(density(ar_null))
lines(xx, 
      dbeta(xx, shape1 = a, shape2 = b), 
      lty = 2,
      lwd = 2)

HDInterval::hdi(ar_null, credMass = 0.9)
quantile(ar_null, c(0.05, 0.95))
qbeta(c(0.05, 0.95), shape1 = a, shape2 = b) # beta approx
qnorm(c(0.05, 0.95), mean = ar_m, sd = ar_s) # normal approx

d1 = quantile(ar_null, c(0.05, 0.95))
d2 = qnorm(c(0.05, 0.95), mean = ar_m, sd = ar_s)

round(abs(d1 - d2), 4) # norm approx

d2 = qbeta(c(0.05, 0.95), shape1 = a, shape2 = b)
round(abs(d1 - d2), 4) # beta approx

# Use just a limited number of MC samples to estimate
# the mean and sd of rho,

# Do the loop again to get understanding about the
# speed gain,

M = 1000

ar_null_m = rep(0, M)

for(j in 1:M){
  
  X = MASS::mvrnorm(n = n, mu = mu0, Sigma = I)
  
  mu0 = rep(mean(X), p)
  #Y = t(X) - mu0
  #R = tcrossprod(Y)/(n - 1)
  
  R = cov(X)
  
  Xm = colMeans(X)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(t0, sigma = R)
  b1 = mvtnorm::dmvt(t0, df = df, sigma = R, log = FALSE)
  
  if(b0 == 0) rho = 0
  if(b0 != 0) rho = b0/b1
  
  rho = ifelse(rho >= 1, 1, rho)
  
  ar_null_m[j] = rho
  
}

ar_limited_m = mean(ar_null_m)
ar_limited_s = sd(ar_null_m)

norm(ar_limited_m - ar_m, type = "2")
norm(ar_limited_s - ar_s, type = "2")

hist(ar_null, probability = TRUE)

lines(xx, 
      dnorm(xx, 
            mean = ar_limited_m, 
            sd = ar_limited_s), 
      lty = 2,
      lwd = 2)

a = ar_limited_m*(n - 1)
b = (1 - ar_limited_m)*(n - 1)

lines(xx,
      dbeta(xx,
            shape1 = a,
            shape2 = b),
      lty = 2,
      lwd = 2,
      col = "red")

plot(density(ar_null))
lines(xx, 
      dnorm(xx, 
            mean = ar_limited_m, 
            sd = ar_limited_s), 
      lty = 2,
      lwd = 2)

lines(xx,
      dbeta(xx,
            shape1 = a,
            shape2 = b),
      lty = 2,
      lwd = 2,
      col = "red")

HDInterval::hdi(ar_null, credMass = 0.9)
quantile(ar_null, c(0.05, 0.95))
qnorm(c(0.05, 0.95), 
      mean = ar_limited_m, 
      sd = ar_limited_s)
qbeta(c(0.05, 0.95),
      shape1 = a,
      shape2 = b)

d2 = qnorm(c(0.05, 0.95), mean = ar_limited_m,
           sd = ar_limited_s)

round(abs(d1 - d2), 4) # norm approx

d2 = qbeta(c(0.05, 0.95), shape1 = a, shape2 = b)
round(abs(d1 - d2), 4) # beta approx
