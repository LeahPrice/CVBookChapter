library(MASS)
library(batchmeans)

# Model specifications
tau <- sqrt(10)
rho <- 0.99
Sigma <- matrix(c(1,rho*tau,rho*tau,tau^2),nrow=2)
SigmaInv <- solve(Sigma)

n <- 100 # number of repeats

logpi <- function(x){return(-0.5*t(x)%*%SigmaInv%*%x)}

set.seed(1)

# MH-MCMC with covariance Sigma from the model
# Saving everything needed for all methods

X <- matrix(NaN,nrow=n,ncol=2)
X[1,] <- c(0.5,0.5)
logpi_curr <- logpi(X[1,])
for (j in 2:n){
  prop <- mvrnorm(1,X[j-1,],Sigma)
  logpi_prop <- logpi(prop)
  MH_ratio <- exp(logpi_prop - logpi_curr)
  if (runif(1)<MH_ratio){
    X[j,] <- prop
    logpi_curr <- logpi_prop
  } else{
    X[j,] <- X[j-1,]
  }
}

##############################################################
# Section 1.2.3

# Evaluating the integrand and performing vanilla Monte Carlo
f <- X[,1]
mean(f)

# Getting the score function
grads <- -X%*%SigmaInv

# Control variates with least squares
coef <- lm(f ~ grads)$coef
mean(f - grads%*%coef[-1])

# Control variates with batch means
coef <- optim(c(0,0),function(theta) bm(f-grads%*%theta)$se)$par
mean(f - grads%*%coef)

##############################################################
# Section 1.3.1

# Control variates with second-order ZVCV
CVs <- cbind(grads, 2+2*X*grads, X[,1]*grads[,2] + X[,2]*grads[,1])
coef <- lm(f ~ CVs)$coef
mean(f - CVs%*%coef[-1])

##############################################################
# Section 1.3.2

# Control functionals with a Gaussian kernel (package)
library(ZVCV)
lambda <- medianTune(unique(X))
CF(f,X,grads,steinOrder=1,kernel_function="gaussian",sigma=lambda)

# Control functionals with a Gaussian kernel (from scratch)
library(Matrix) # For nearPD

# Removing duplicates
dups <- which(duplicated(X))
Xu <- X[-dups,]
gradsu <- grads[-dups,]
fu <- f[-dups]
nu <- nrow(Xu)

# Calculating matrix of squared norms
Z <- dist(Xu)^2
lambda2 <- 0.5*median(Z)
Z <- as.matrix(Z)

# Calculating the K0 matrix
d <- 2
k <- exp(-Z/lambda2)
K0 <- matrix(NaN,nrow=nu,ncol=nu)
for (i in 1:nu){
  for (j in i:nu){
    x <- Xu[i,]
    y <- Xu[j,]
    
    k_ <- k[i,j]
    k_gradx <- -2/lambda2*k_*(x-y)
    k_grady <- 2/lambda2*k_*(x-y)
    k_gradxgrady <- -4/lambda2^2*k_*Z[i,j] + 2/lambda2*k_*d
    
    K0[i,j] <- K0[j,i] <- k_gradxgrady + t(gradsu[i,])%*%k_grady + 
      t(gradsu[j,])%*%k_gradx + t(gradsu[i,])%*%gradsu[j,]*k_
  }
}

# Getting the estimate
K0 <- nearPD(K0) # For numerical stability in the inverse
K0inv <- solve(K0)
t(rep(1,nu))%*%K0inv%*%fu/(t(rep(1,nu))%*%K0inv%*%rep(1,nu))


##############################################################
# Section 1.3.3

# SECF with a Gaussian kernel and first order polynomial (package)
SECF(f,X,grads,polyorder=1,steinOrder=1,kernel_function="gaussian",sigma=lambda)

# SECF with a Gaussian kernel and first order polynomial (from scratch)
Phi <- cbind(1,gradsu)
solve(t(Phi)%*%K0inv%*%Phi,t(Phi)%*%K0inv%*%fu)[1]

##############################################################
# Section 1.4.1

# Control variates of the form g - Pg for random-scan Gibbs P
CVs <- 0.5*cbind(X[,1] - rho*X[,2]/tau,X[,2] - rho*tau*X[,1])
coef <- lm(f ~ CVs)$coef
mean(f - CVs%*%coef[-1])

##############################################################
# Section 1.5

# Calculating the integrand of interest
f_curr <- X[-n,1]
f_prop <- Xprop[,1]
f_next <- X[-1,1]

# Calculating the MH acceptance probability from the ratios
MH <- pmin(1,R)
MH_reverse <- pmin(1,1/R)

# Hammer and Tjelmeland 
CVs <- (1-accept)*MH*f_prop - accept*(1-MH_reverse)*f_curr
beta <- lm(f_curr~CVs)$coefficient
mean(f_curr - CVs*beta[-1])

# Delmas and Jourdain
CVs <- MH*f_prop + (1-MH)*f_curr - f_next
beta <- lm(f_curr~CVs)$coefficient
mean(f_curr - CVs*beta[-1])
