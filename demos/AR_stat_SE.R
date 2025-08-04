
# The null hypothesis of this test is 

# H0: mu = mu0

library(ggplot2)

set.seed(1)

pwr::pwr.t.test(d = 0.5,
                power = 0.8,
                sig.level = 0.05,
                type = "two.sample",
                alternative = "two.sided")

n = 64

# Different effect sizes,

eta = c(0.3, 0.4, 0.5, 0.6)

m = length(eta)

# The distribution: norm, t, pexp, unif, laplace

pdist = "norm"

if(pdist == "t") df = 3 # for the t-distribution

# distribution under null

p = 2

mu0 = rep(0, p)

I = diag(1, p)

M = 10^4

rho = c()

rhod = matrix(0, M, m)

for(i in 1:m){
  
  mu = c(0, eta[i])
  
  #Sigma = I
  
  Sigma = matrix(-0.1, p, p)
  
  diag(Sigma) = 1
  
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
  
  y = X[, 1]
  x = X[, 2]
  
  p = ncol(X)
  
  Xm = colMeans(X, na.rm = TRUE)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  
  r = cor(X)
  
  rho[i] = mvtnorm::dmvnorm(t0, sigma = r)/mvtnorm::dmvt(t0, df = df, sigma = r, log = FALSE)
  
  Ru = matrix(runif(n*M), n, M)
  
  tempU = apply(Ru, 2, function(x) rho[i] > x)
  
  rhod[, i] = colMeans(tempU)
  
  rho[i] = ifelse(rho[i] >= 1, 1, rho[i])
  
  cat("\r", i)
  
}

Data = data.frame(Effect = rep(eta, each = M),
                  rho = c(rhod))

rD = data.frame(Effect = eta, 
                rho = rho)

SEr = sqrt(rho*(1 - rho)/n)

S = sqrt(rho*(1 - rho))
S

rDSE = data.frame(Effect = eta, 
                  rhoSE = SEr,
                  low = rho - 2*SEr,
                  up = rho + 2*SEr)

Dnorm = data.frame(rho = rho,
                   SE = SEr)

round(SEr, 4)
round(apply(rhod, 2, sd), 4)

round(SEr - apply(rhod, 2, sd), 4)

ggplot(data = Data, aes(rho)) +
  geom_histogram(aes(y=..density..))+
  geom_density() +
  facet_wrap(vars(Effect), scales = "free") +
  geom_vline(data = rD, 
             aes(xintercept = rho),
             linewidth = 0.7,
             linetype = "dashed",
             color = "red") +
  geom_vline(data = rDSE, 
             aes(xintercept = low),
             linewidth = 0.7,
             linetype = "dashed",
             color = "blue") +
  geom_vline(data = rDSE, 
             aes(xintercept = up),
             linewidth = 0.7,
             linetype = "dashed",
             color = "blue") +
  ylab("Density") +
  xlab(expression(rho))

par(mfrow = c(2, 2))

xx = seq(0, 1, length.out = 1000)

hist(rhod[, 1], 
     main = eta[1],
     probability = TRUE, 
     xlab = expression(rho))
lines(xx, dnorm(xx, mean = rho[1], sd = SEr[1]),
lwd = 2, lty = 2, col = "red")

hist(rhod[, 2], 
     main = eta[2],
     probability = TRUE, 
     xlab = expression(rho))
lines(xx, dnorm(xx, mean = rho[2], sd = SEr[2]),
      lwd = 2, lty = 2, col = "red")

hist(rhod[, 3], 
     main = eta[3],
     probability = TRUE, 
     xlab = expression(rho))
lines(xx, dnorm(xx, mean = rho[3], sd = SEr[3]),
      lwd = 2, lty = 2, col = "red")

hist(rhod[, 4], 
     main = eta[4],
     probability = TRUE, 
     xlab = expression(rho))
lines(xx, dnorm(xx, mean = rho[4], sd = SEr[4]),
      lwd = 2, lty = 2, col = "red")

par(mfrow = c(1, 1))

