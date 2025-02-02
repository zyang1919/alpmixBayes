llmix <- function(d,init_values,mcmc_values,prior){
  
  ## likelihood 
  falp <- function(alp,p1,p2,mu1,mu2,d)
  {
    fbar1 <- 1 - plnorm(d, meanlog = mu1, sdlog = 0.5)
    fbar2 <- 1 - plnorm(d, meanlog = mu2, sdlog = 0.1)
    f1 <- dlnorm(d, meanlog = mu1, sdlog = 0.5)
    f2 <- dlnorm(d, meanlog = mu2, sdlog = 0.1)
    r1 <- f1/fbar1
    r2 <- f2/fbar2
    
    if(alp==0){
      half1 <- fbar1^(p1)*fbar2^(p2)
      half2 <- p1*r1+p2*r2
    } else {
      half1 <- (p1*fbar1^(alp)+p2*fbar2^(alp))^((1/alp)-1)
      half2 <- p1*fbar1^(alp)*r1+p2*fbar2^(alp)*r2
    }
    ret <- half1*half2
    return(ret)
  }
  
  
  ## unlist mcmc 
  nburn = mcmc_values$nburn
  niter = mcmc_values$niter
  thin = mcmc_values$thin
  aiter <- nburn + niter
  nsamp <- ceiling(niter/thin)
  
  ## unlist the initial values
  alptemp = init_values$alpha
  p1temp <- init_values$p
  mu1temp=init_values$th1; mu2temp=init_values$th2
  
  alpsamp <- rep(0,nsamp)
  p1samp <- rep(0,nsamp)
  mu1samp <- rep(0,nsamp); mu2samp <- rep(0,nsamp)
  
  ## hyperparameters
  mualp <- prior$ap[1]; sig2alp <- prior$ap[2]
  c1 <- prior$pp[1]; c2 <- prior$pp[2]
  mu01 <- prior$th1p[1]; sig201 <- prior$th1p[2]
  mu02 <- prior$th2p[2]; sig202 <- prior$th2p[2]
  
  logdenalp <- sum(log(falp(alptemp,p1temp,1-p1temp,mu1temp,mu2temp,d)),na.rm = TRUE) -
    0.5*(alptemp-mualp)^2/sig2alp
  logdenp <- sum(log(falp(alptemp,p1temp,1-p1temp,mu1temp,mu2temp,d)),na.rm=TRUE) +
    (c1-1)*log(p1temp) + (c2-1)*log(1-p1temp)
  logdenmu1 <- sum(log(falp(alptemp,p1temp,1-p1temp,mu1temp,mu2temp,d)),na.rm=TRUE) -
    0.5*(log(mu1temp)-mu01)^2/sig201
  logdenmu2 <- sum(log(falp(alptemp,p1temp,1-p1temp,mu1temp,mu2temp,d)),na.rm=TRUE) -
    0.5*(log(mu2temp)-mu02)^2/sig202
  
  storeIndex <- 0
  
  for(k in 1:aiter){
    eflag <- 0
    # alpha
    # point mass
    pzero <- exp(sum(log(0.01002004*falp(alptemp,p1temp,1-p1temp,mu1temp,
                                         mu2temp,d)),na.rm = TRUE)) # gap of 500 points between 0 and 5 = 0.01002004
    #  if (is.nan(pzero) == TRUE){
    #    th2samp <- th2
    #   next
    #  }
    
    u <- runif(1,0,1)
    
    if (pzero > 0 & is.infinite(pzero) == TRUE){
      pzero <- 1
    } 
    
    if (pzero < 0 & is.infinite(pzero) == TRUE){
      pzero <- 0
    }
    
    
    if(u < pzero){
      alptemp <- 0
    } else{
      alpcand <- rnorm(1,alptemp,1)
      lognumalp <- sum(log(falp(alpcand,p1temp,1-p1temp,mu1temp,mu2temp,
                                d)),na.rm=TRUE)-0.5*(alpcand-mualp)^2/sig2alp
      logalp <- lognumalp - logdenalp
      logU <- log(runif(1,0,1))
      
      if(is.nan(logalp)==TRUE|is.infinite(logalp)==TRUE){
        eflag <- 1
      }
      
      if(eflag){
        next
      }
      
      if(logU < logalp){
        alptemp <- alpcand
        logdenalp <- lognumalp
      }
    }
    
    # p=(p1,p2)
    pcand <- runif(1,0,1)
    lognump <- sum(log(falp(alptemp,pcand,1-pcand,mu1temp,mu2temp,d)),na.rm=TRUE) +
      (c1-1)*log(pcand) + (c2-1)*log(1-pcand)
    logalp <- lognump - logdenp
    logU <- log(runif(1,0,1))
    
    if(is.nan(logalp)==TRUE|is.infinite(logalp)==TRUE){
      eflag <- 1
    }
    
    if(eflag){
      next
    }
    
    if(logU < logalp){
      p1temp <- pcand
      logdenp <- lognump
    }
    
    # mu1: with Normal prior mean=mu1t and var=1
    
    mu1cand <- rlnorm(1,mu1temp,sdlog = 1)  # lognormal mean=mu1temp var=1
    lognummu1 <- sum(log(falp(alptemp,p1temp,1-p1temp,mu1cand,mu2temp,
                              d)),na.rm=TRUE) - 0.5*(log(mu1cand)-mu01)^2/sig201
    logalp <- lognummu1 - logdenmu1 + 
      dlnorm(mu1temp, meanlog = mu1cand, sdlog = 1, log = TRUE) -
      dlnorm(mu1cand, meanlog = mu1temp, sdlog = 1, log = TRUE)
    
    logU <- log(runif(1,0,1))
    
    if(is.nan(logalp)==TRUE|is.infinite(logalp)==TRUE){
      eflag <- 1
    }
    
    if(eflag){
      next
    }
    
    if(logU < logalp){
      mu1temp <- mu1cand
      logdenmu1 <- lognummu1
    }
    
    # mu2: with Normal prior mean=mu2t and var=1
    
    mu2cand <- rlnorm(1,mu2temp,sdlog = 1)  # lognormal mean=mu2temp var=1
    lognummu2 <- sum(log(falp(alptemp,p1temp,1-p1temp,mu1temp,mu2cand,
                              d)),na.rm=TRUE) - 0.5*(log(mu2cand)-mu02)^2/sig202
    logalp <- lognummu2 - logdenmu2 + 
      dlnorm(mu2temp, meanlog = mu2cand, sdlog = 1, log = TRUE) -
      dlnorm(mu2cand, meanlog = mu2temp, sdlog = 1, log = TRUE)
    
    logU <- log(runif(1,0,1))
    
    if(is.nan(logalp)==TRUE|is.infinite(logalp)==TRUE){
      eflag <- 1
    }
    
    if(eflag){
      next
    }
    
    if(logU < logalp){
      mu2temp <- mu2cand
      logdenmu2 <- lognummu2
    }
    
    if(k > nburn && ((k - nburn) %% thin) == 0){
      storeIndex <- storeIndex + 1
      alpsamp[storeIndex] <- alptemp
      p1samp[storeIndex] <- p1temp
      mu1samp[storeIndex] <- mu1temp
      mu2samp[storeIndex] <- mu2temp
    }
  }
  
  # Compute posterior means and credible intervals
  alpha_mean <- mean(alpsamp)
  alpha_ci <- quantile(alpsamp, probs = c(0.025, 0.975))
  
  p_means <- mean(p1samp)
  p_ci <- quantile(p1samp, probs=c(0.025, 0.975))
  
  mu1_mean <- mean(mu1samp)
  mu1_ci <- quantile(mu1samp, probs = c(0.025, 0.975))
  
  mu2_mean <- mean(mu2samp)
  mu2_ci <- quantile(mu2samp, probs = c(0.025, 0.975))
  
  
  # Organize results into a structured list
  result <- list(
    posterior_means = list(
      alpha = alpha_mean,
      p1 = p_means,
      mu1 = mu1_mean,
      mu2 = mu2_mean
    ),
    credible_intervals = list(
      alpha = alpha_ci,
      p1 = p_ci,
      mu1 = mu1_ci,
      mu2 = mu2_ci
    )
  )
}