
set.seed(1)

r = 0.8

Sigma = matrix(r, 2, 2)
diag(Sigma) = 1

mu = c(0, 0.4)

p = 2

N = c(30, 40, 70)

M = 10^5

i = 1

Q1 = data.frame()

png(filename = "demos/null_and_alt_distributions/figures/ar_stat_suff.png",
    width = 400, 
    height = 200, 
    units = "mm",
    res = 700)

par(mfrow = c(1, 3), 
    cex.axis=1.5,
    mar = c(5.1, 5.1, 4.1, 2.1))

for(n in N){
  
  X = MASS::mvrnorm(n = n,
                    mu = mu,
                    Sigma = Sigma)
  
  
  mu0 = rep(mean(X), p)
  
  R = cov(X)
  
  Xm = colMeans(X)
  
  t0 = sqrt(n)*(Xm - mu0)
  
  df = n - 1
  b0 = mvtnorm::dmvnorm(t0, sigma = R)
  b1 = mvtnorm::dmvt(t0, df = df, sigma = R, log = FALSE)
  
  if(b0 == 0) rho = 0
  if(b0 != 0) rho = b0/b1
  
  u = matrix(runif(n*M), n, M)
  
  I = rho > u
  
  d = colSums(I)
  
  plot(table(d)/M, type = "h",
       xlab = expression(nT(X)),
       ylab = "Probability",
       lwd = 4,
       cex.lab = 2,
       col = "darkgray")
  
  xx = 0:n
  
  probpb = poibin::dpoibin(xx, pp = rep(rho, n))
  
  lines(xx+0.2, 
        probpb, 
        type = "h", 
        lwd = 4, 
        col = rgb(red=0, green=0, blue=1, alpha=0.5))
  
  abline(v = n*rho, lty = 2, lwd = 4)
  
  lab = paste0("(", LETTERS[i], ")")
  
  title(main = lab, adj = 0, cex.main = 2.5)
  
  xx = seq(0, n, length.out = 10^4)
  
  norm_apprx = dnorm(xx, 
                     mean = n*rho,
                     sd = sqrt(n*rho*(1 - rho)))
  
  lines(xx, 
        norm_apprx, 
        lwd = 2, 
        col = rgb(red=1, green=0, blue=0, alpha=0.7))
  
  Qtemp = data.frame(qpb = poibin::qpoibin(c(0.025, 0.975), pp = rep(rho, n))/n,
                     qn = qnorm(c(0.025, 0.975), mean = rho,sd = sqrt(rho*(1 - rho)/n)),
                     qsample = quantile(d/n, c(0.025, 0.975)),
                     n = rep(n, 2))

  Q1 = rbind(Q1, Qtemp)
  
  i = i + 1
  
}

dev.off()

#

i = 1

Q2 = data.frame()

png(filename = "demos/null_and_alt_distributions/figures/ar_stat_sum.png",
    width = 400, 
    height = 200, 
    units = "mm",
    res = 700)

par(mfrow = c(1, 3), 
    cex.axis=1.5,
    mar = c(5.1, 5.1, 4.1, 2.1))

for(n in N){
  
  x = extraDistr::rlst(n, mu = 0.1, sigma = 1, df = 3)
  
  m = mean(x)
  
  s = sd(x)
  
  fhat = kdevine::kde1d(x)
  
  fhat = kdevine::dkde1d(x, fhat)
  
  rhox = dnorm(x, mean = m, sd = s)/fhat
  
  rhox[is.na(rhox)] = 0
  
  a = rhox[rhox < 1]
  b = sum(rhox >= 1)
  
  rho = (sum(a) + b)/n
  
  rhox = ifelse(rhox > 1, 1, rhox)
  
  U = matrix(runif(n*M), n, M)
  
  f = function(x) rhox > x
  
  I = apply(U, 2, f)
  
  d = colSums(I)
  
  plot(table(d)/M, type = "h",
       xlab = expression(nT(X)),
       ylab = "Probability",
       lwd = 4,
       cex.lab = 2,
       col = "darkgray")
  
  xx = 0:n
  
  probpb = poibin::dpoibin(xx, pp = rhox)
  
  lines(xx+0.2, 
        probpb, 
        type = "h", 
        lwd = 4, 
        col = rgb(red=0, green=0, blue=1, alpha=0.5))
  
  abline(v = n*rho, lty = 2, lwd = 4)
  
  xx = seq(0, n, length.out = 10^4)
  
  norm_apprx = dnorm(xx, 
                     mean = n*rho,
                     sd = sqrt(n*rho*(1 - rho)))
  
  # lines(xx, 
  #       norm_apprx, 
  #       lwd = 2, 
  #       col = rgb(red=1, green=0, blue=0, alpha=0.7))
  
  lab = paste0("(", LETTERS[i], ")")
  
  title(main = lab, adj = 0, cex.main = 2)
  
  Qtemp = data.frame(qpb = poibin::qpoibin(c(0.025, 0.975), pp = rhox)/n,
                     qn = qnorm(c(0.025, 0.975), mean = rho,sd = sqrt(rho*(1 - rho)/n)),
                     qsample = quantile(d/n, c(0.025, 0.975)),
                     n = rep(n, 2))
  
  Q2 = rbind(Q2, Qtemp)
  
  i = i + 1
  
}

dev.off()