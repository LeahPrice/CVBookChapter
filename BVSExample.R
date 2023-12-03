source("functions_for_BVS.R")
set.seed(1)

# Define parameters
n <- 70  # Number of observations
p <- 5   # Number of covariates
c <- 10^3  # Prior hyperparameter for the covariance matrix
prior_p_incl <- 1/p  # Prior probability of inclusion

# Generate simulated data
hyper_par <- simulate_data(n = n, p = p, c = c, SNR = 3, scenario = 1)

## Commented out block for gold standard
# 
# n <- 5000000  # Number of iterations
# burn_in <- 5000
# thin_GS <- 1
# MCMC <- GS(p = p, hyper_par = hyper_par, T = n, burn_in = burn_in, thin = thin_GS)
# mean(MCMC$gamma_1[(burn_in+1):(burn_in+n)])

# Gold standard estimation
gold_est <- 0.2664132  # from 5,000,000 iterations with 5000 burn-in
burn_in <- 5000
n <- 5000
MCMC <- GS(p = p, hyper_par = hyper_par, T = n, burn_in = burn_in)

gamma_1 <- MCMC$gamma_1
p_gamma_1 <- MCMC$p_gamma_1
gamma_2 <- MCMC$gamma_2
p_gamma_2 <- MCMC$p_gamma_2

# Standard Monte Carlo estimator
f <- gamma_1
mean(f)

# Using a single function h(gamma) = gamma_1
CV1 <- gamma_1 - (1/p*p_gamma_1 + (p-1)/p*gamma_1)
coefs1 <- lm(f ~ CV1)$coef
mean(f - coefs1[-1]*CV1)

# Using h(gamma) = theta_1*gamma_1 + theta_2*gamma_2
CV2 = cbind(CV1, gamma_2 - (1/p*p_gamma_2 + (p-1)/p*gamma_2))
coefs2 <- lm(f ~ CV2)$coef
mean(f - CV2%*%coefs2[-1])

## Plot of the cumulative estimate
est_CV2 <- est_CV1 <- rep(NaN, n + burn_in)
for( i in 3:(n+burn_in)){
  coefs1 <- lm(f[1:i] ~ CV1[1:i])$coef
  est_CV1[i] = mean(f[1:i] - coefs1[-1]*CV1[1:i])
  coefs2 <- lm(f[1:i] ~ CV2[1:i,])$coef
  est_CV2[i] = mean(f[1:i] - CV2[1:i,]%*%coefs2[-1])
}

par(mfrow = c(1,2), mar=c(4,4,2,2))
plot(cumsum(gamma_1)/c(1:(burn_in+n)),
     ylim=c(0,0.5),type="l",xlab="Sampler iteration", ylab=expression(paste("E[", gamma[1],"]")))
lines(est_CV1, col = 'red', lwd=1)
lines(est_CV2, col = 'blue', lwd=2)
abline(h = gold_est, lty=2, lwd=2, col = "magenta")
legend("topright", legend = c("Monte Carlo", "CV1", "CV2", "Gold Standard"), col=c("black","red","blue","magenta"), lwd=2)

plot(cumsum(gamma_1)/c(1:(burn_in+n)),
     ylim=c(0.26,0.27),type="l",xlab="Sampler iteration", ylab=expression(paste("E[", gamma[1],"]")))
lines(est_CV1, col = 'red', lwd=1)
lines(est_CV2, col = 'blue', lwd=2)
abline(h = gold_est, lty=2, lwd=2, col="magenta")

par(mfrow = c(1,1), mar=c(4,5,4,4))

# SIMULATION
# Compare the performance relative to Rao-Blackwellisation

n<-5000 #n.iterations
burn_in<-5000

M=100
MC_est <- rep(0, M)
CV1_est <- rep(0, M)
CV2_est <- rep(0, M)
RB_est <- rep(0, M)
for( i in 1:M ){
  print(i)
  MCMC<-GS(p=p,hyper_par=hyper_par,T=n,burn_in=burn_in)
  f <- MCMC$gamma_1
  MC_est[i] = mean(f)

  CV1 = MCMC$gamma_1 - (1/p*MCMC$p_gamma_1 + (p-1)/p*MCMC$gamma_1)
  coefs1 <- lm(f ~ CV1)$coef
  CV1_est[i] = mean(f - CV1*coefs1[-1])

  CV2 = cbind(CV1, MCMC$gamma_2 - (1/p*MCMC$p_gamma_2 + (p-1)/p*MCMC$gamma_2))
  coefs2 <- lm(f ~ CV2)$coef
  CV2_est[i] = mean(f - CV2%*%coefs2[-1])
  
  RB_est[i] = mean(MCMC$p_gamma_1)
}

boxplot(list(MontCarlo=MC_est, RaoBlackWell=RB_est, CV1=CV1_est, CV2=CV2_est), ylab = expression(paste("E[", gamma[1],"]")))
abline(h=gold_est)


