
source("functions/ar_compare_group_means.R")

set.seed(1)

r = -0.9

Sigma = matrix(r, 2, 2)
diag(Sigma) = 1

mu = rnorm(2)

pwr::cohen.ES(test = "t", size = "small")

pwr::pwr.t.test(d = 0.2, 
                power = 0.8,
                type = "two.sample",
                alternative = "two.sided")

n = 394

M = 5*10^4

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

round(SErho, 4)

sd(rhod)

hist(rhod, probability = T)
xx = seq(-1, 2, length.out = 1000)
lines(xx, 
      dnorm(xx, mean = rho, sd = SErho),
      lty = 2,
      lwd = 2,
      col = "red")

a = rho*(n - 1)
b = (1 - rho)*(n - 1)

yy = seq(0, 1, length.out = 1000)

lines(yy,
      dbeta(yy, shape1 = a, shape2 = b),
      lty = 2,
      lwd = 2,
      col = "blue")

d1 = rho + c(-1, 1)*qnorm(0.975)*SErho
d2 = qbeta(c(0.025, 0.975), shape1 = a, shape2 = b)
dhdi = HDInterval::hdi(rhod, credMass = 0.95)

d1 # normal
d2 # beta
dhdi

round(abs(d1 - c(dhdi)), 5) # normal
round(abs(d2 - c(dhdi)), 5) # beta
