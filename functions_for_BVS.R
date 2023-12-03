
#  Code adapted from G. Zanella
#  GitHub repository: https://github.com/gzanella/TGS
#
#  Title: Scalable Importance Tempering and Bayesian Variable Selection
#  Authors: Giacomo Zanella, Gareth Roberts
#  Published in: Journal of the Royal Statistical Society Series B: Statistical Methodology
#  Year: 2019
#  DOI: https://doi.org/10.1111/rssb.12316
#
#  Description: This code implements Gibbs samplers for BVS and is based on the
#  implementation provided by G. Zanella in their GitHub repository.
#
#  Adapted:
#    + Return additional quantities needed for control variates.
#    + Switch to the non-Metropolised version of the Gibbs sampler.


simulate_data<-function(n,p,c,SNR=2,scenario=1){

  # simulate_data function generates synthetic data based on different scenarios.
  # Parameters:
  #   - n: Number of observations in the generated dataset.
  #   - p: Number of variables in the dataset.
  #   - c: Hyperparameter for parameter covariance matrix.
  #   - SNR: Signal-to-Noise Ratio (default value = 2).
  #   - scenario: Indicates the type of dataset scenario to be generated (default value = 1).
  # Returns:
  #   - hyper_par: A list containing precomputed matrices and necessary data.

  require(MASS)
  Sigma<-diag(1,p)
  if(scenario==1){
    # variables 1 and 2 strongly correlated
    rho<-0.99
    Sigma[1,2]<-Sigma[2,1]<-rho
    X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
    beta0 <- matrix(c(1,rep(0,p-1)),ncol=1)
  }
  if(scenario==2){
    #correlated design of Wang et al (2011) and Huang et al (2016); with different variance
    rho<-0.9
    for (i in 1:3){for(j in 1:3){if(i!=j)Sigma[i,j]<-rho}}
    for (i in 4:6){for(j in 4:6){if(i!=j)Sigma[i,j]<-rho}}
    require(MASS)
    X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
    beta0<-matrix(c(3,3,-2,3,3,-2,rep(0,p-6)),ncol=1)
  }
  if(scenario==3){
    # uncorrelated_design
    X <- matrix(rnorm(n*p),nrow=n,ncol=p)
    beta0 <- matrix(c(2,-3,2,2,-3,3,-2,3,-2,3,rep(0,p-10)),ncol=1)
  }
  sigma<-1
  beta<-SNR*sqrt(log(p)/n)*beta0
  y <- X %*% beta + rnorm(n,mean = 0,sd = sigma^2)
  X<-t(t(X)-colMeans(X))
  y<-y-mean(y)

  ### precompute matrices needed for samplers and create list to output ###
  XtX<-t(X)%*%X
  ytX<-t(y)%*%X
  yty<-sum(y^2)
  hyper_par<-list(
    n=n,y=y,X=X,XtX=XtX, ytX=ytX, yty=yty,
    prior_p_incl=prior_p_incl,
    c=c
  )
  return(hyper_par)
}


GS<-function(p,gamma_start=NULL,
             T,burn_in=0,thin=1,
             hyper_par=NULL,vars_selected=c(1,2)){

  # GS function performs random scan Gibbs Sampler for Bayesian variable selection problems.

  # Parameters:
  #   - p: Number of variables in the dataset.
  #   - gamma_start: Initial vector of inclusion indicators (default: NULL).
  #   - T: Total number of iterations for the Gibbs Sampler.
  #   - burn_in: Number of burn-in iterations (default: 0).
  #   - thin: Thinning parameter for thinning the MCMC chain (default: 1).
  #   - hyper_par: List containing precomputed matrices and necessary data for analysis.
  #   - vars_selected: Indices of two variables for which inclusion probabilities will be tracked (default: c(1, 2)).

  # Returns:
  #   - A list containing:
  #     - indices_sequence: Sequence of indices updated during MCMC iterations.
  #     - gamma_1: Sequence of inclusion indicators (MCMC samples) for the first selected variable.
  #     - gamma_2: Sequence of inclusion indicators (MCMC samples) for the second selected variable.
  #     - p_gamma_1: Gibbs full conditional probabilities associated with gamma_1
  #     - p_gamma_2: Gibbs full conditional probabilities associated with gamma_2


  if(is.null(gamma_start)){
    ## initialize vector of inclusion indicators
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start

  ## initialize vector to be output
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  p_gamma_1 <- rep(NA,burn_in+T)
  p_gamma_2 <- rep(NA,burn_in+T)

  ## MCMC iteration
  for (t in 1:(T+burn_in)){

    ## perform random scan Gibbs Sampling step
    for (iter in 1:thin){
      j<-sample.int(n = p,size = 1)
      fc_j<-single_full_cond(j=j,gamma=gamma,hyper_par=hyper_par)$fc_j

      ap<-1-fc_j## not metropolized (straight Gibbs Sampler)

      if(runif(1) < ap){
        gamma[j]<-1-gamma[j]
        indices_sequence[t]<-j
      }
    }

    ## store indices sequence for output analysis
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]

    ##
    fc_1 = single_full_cond(j=vars_selected[1],gamma=gamma,hyper_par=hyper_par)$fc_j
    if(gamma[vars_selected[1]] == 1){
      p_gamma_1[t] = fc_1
    } else{
      p_gamma_1[t] = 1-fc_1
    }
    fc_2 = single_full_cond(j=vars_selected[2],gamma=gamma,hyper_par=hyper_par)$fc_j
    if(gamma[vars_selected[2]] == 1){
      p_gamma_2[t] = fc_2
    } else{
      p_gamma_2[t] = 1-fc_2
    }
    ##
  }
  return(list(indices_sequence=indices_sequence,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              p_gamma_1 = p_gamma_1,
              p_gamma_2 = p_gamma_2,
              final_gamma=gamma))
}

single_full_cond<-function(j,gamma,hyper_par=NULL,stored=NULL){

  # single_full_cond function computes the full conditional probability for a specific variable index 'j'
  # in a Bayesian variable selection problem using Gibbs Sampler.

  # Parameters:
  #   - j: Index of the variable for which the full conditional probability is being computed.
  #   - gamma: Vector of inclusion indicators (binary indicators for each variable).
  #   - hyper_par: List containing precomputed matrices and necessary data for analysis (default: NULL).
  #   - stored: Placeholder for storing data (default: NULL).

  # Returns:
  #   - A list containing:
  #     - fc_j: Full conditional probability of the current value of gamma[j].

  included<-which(gamma==1)
  p_gamma<-sum(gamma)
  n<-hyper_par$n
  p<-length(gamma)
  h<-hyper_par$prior_p_incl
  c<-hyper_par$c
  yty<-hyper_par$yty
  ytX<-hyper_par$ytX
  XtX<-hyper_par$XtX
  if(p_gamma==0){
    S_old<-yty
    S_new<-yty-c/(1+c)*ytX[j]^2/XtX[j,j]
    post_ratio<-h/(1-h)*(S_old/S_new)^(n/2)/(sqrt(1+c))
    fc_j<-1/(1+post_ratio)# fc is the full cond prob of the current value of gamma[j]
  }else{
    ytX_gamma<-matrix(ytX[included],nrow = 1,ncol = p_gamma)
    ytXFXty<-sum((solve(t(chol(XtX[included,included])))%*%t(ytX_gamma))^2)
    S_old<-yty-c/(1+c)*ytXFXty
    F<-solve(XtX[included,included])
    if(gamma[j]==0){
      xjtX<-matrix(XtX[j,included],nrow=1,ncol=length(included))
      d<-c(1/(XtX[j,j]-xjtX%*%F%*%t(xjtX)))
      ytXFXtxj<-ytX_gamma%*%F%*%t(xjtX)
      S_new<-yty-c/(1+c)*(ytXFXty+d*(ytXFXtxj-ytX[j])^2)
      post_ratio<-h/(1-h)*(S_old/S_new)^(n/2)/(sqrt(1+c))
      fc_j<-1/(1+post_ratio)# fc is the full cond prob of the current value of gamma[j]
    }
    if(gamma[j]==1){
      if(p_gamma==1){
        S_new<-yty
      }else{
        pos_rem<-which(included==j)
        included_new<-included[-pos_rem]
        ytX_gamma_new<-matrix(ytX[included_new],nrow = 1,ncol = p_gamma-1)
        F_new<-F[-pos_rem,-pos_rem]-(as.matrix(F[-pos_rem,pos_rem],nrow=p_gamma-1)%*%F[pos_rem,-pos_rem])/F[pos_rem,pos_rem]
        S_new<-yty-c/(1+c)*ytX_gamma_new%*%F_new%*%t(ytX_gamma_new)
      }
      post_ratio<-(1-h)/h*(S_old/S_new)^(n/2)*(sqrt(1+c))
      fc_j<-1/(1+post_ratio)# fc_j is the full cond prob of the current value of gamma[j]
    }
  }
  output<-list(fc_j=fc_j,stored=NULL)
  return(output)
}
