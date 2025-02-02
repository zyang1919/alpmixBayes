wwmix <- function(d,init_values,mcmc_values,prior){
  
  ## likelihood 
  falp <- function(alp,p1,p2,th1,th2,d)
  {
    if(alp==0){
      pxth1 <- p1*d/(th1^2)
      pxth2 <- p2*d/(th2^2)
      
      Falp <- exp(-pxth1*d-pxth2*d)
      Fr <- 2*pxth1 + 2*pxth2
    } else {
      alpx2 <- alp*d^2
      alpx2th1 <- alpx2/(th1^2)
      alpx2th2 <- alpx2/(th2^2)
      
      Falp <- (p1*exp(-alpx2th1) + p2*exp(-alpx2th2))^(1/alp-1)
      Fr   <- 2*p1*d/(th1^2)*exp(-alpx2th1) + 2*p2*d/(th2^2)*exp(-alpx2th2)
    }
    ret <- Falp*Fr
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
  th1temp=init_values$th1; th2temp=init_values$th2
  
  alpsamp <- rep(0,nsamp)
  p1samp <- rep(0,nsamp)
  th1samp <- rep(0,nsamp); th2samp <- rep(0,nsamp)
  
  # hyperparameters
  mualp <- prior$ap[1]; sig2alp <- prior$ap[2]
  c1 <- prior$pp[1]; c2 <- prior$pp[2]
  A1 <- prior$th1p[1]; B1 <- prior$th1p[2]
  A2 <- prior$th2p[1]; B2 <- prior$th2p[2]
  
  logdenalp <- sum(log(falp(alptemp,p1temp,1-p1temp,th1temp,th2temp,d))) -
    0.5*(alptemp-mualp)^2/sig2alp
  logdenp <- sum(log(falp(alptemp,p1temp,1-p1temp,th1temp,th2temp,d))) +
    (c1-1)*log(p1temp) + (c2-1)*log(1-p1temp)
  logdenth1 <- sum(log(falp(alptemp,p1temp,1-p1temp,th1temp,th2temp,d))) +
    (A1-1)*log(th1temp) - th1temp/B1
  logdenth2 <- sum(log(falp(alptemp,p1temp,1-p1temp,th1temp,th2temp,d))) +
    (A2-1)*log(th2temp) - th2temp/B2
  
  storeIndex <- 0
  
  for(k in 1:aiter){
    eflag <- 0
    # alpha
    # point mass
    pzero <- exp(sum(log(0.01002004*falp(alptemp,p1temp,1-p1temp,th1temp,
                                         th2temp,d)))) # gap of 500 points between 0 and 5 = 0.01002004
    
    u <- runif(1,0,1)
    
    if(u < pzero){
      alptemp <- 0
    } else{
      alpcand <- rnorm(1,alptemp,.1)
      lognumalp <- sum(log(falp(alpcand,p1temp,1-p1temp,th1temp,th2temp,
                                d)))-0.5*(alpcand-mualp)^2/sig2alp
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
    lognump <- sum(log(falp(alptemp,pcand,1-pcand,th1temp,th2temp,d))) +
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
    
    # theta1: with Gamma prior mean=th1t and var=c
    c <- .1
    a <- th1temp^2/c
    b <- c/th1temp
    th1cand <- rgamma(1,a,scale=b)  # Gamma mean=th1temp var=c
    lognumth1 <- sum(log(falp(alptemp,p1temp,1-p1temp,th1cand,th2temp,
                              d)),na.rm=TRUE) + (A1-1)*log(th1cand)-th1cand/B1
    logalp <- lognumth1 - logdenth1 + 
      lgamma(th1temp^2/c)+(th1temp^2/c)*log(c/th1temp)-
      lgamma(th1cand^2/c)-(th1cand^2/c)*log(c/th1cand) +
      (th1cand^2/c-1)*log(th1temp) - (th1temp^2/c-1)*log(th1cand)
    logU <- log(runif(1,0,1))
    
    if(is.nan(logalp)==TRUE|is.infinite(logalp)==TRUE){
      eflag <- 1
    }
    
    if(eflag){
      next
    }
    
    if(logU < logalp){
      th1temp <- th1cand
      logdenth1 <- lognumth1
    }
    
    # theta2: with Gamma prior mean=th2t and var=c
    c <- .1
    a <- th2temp^2/c
    b <- c/th2temp
    th2cand <- rgamma(1,a,scale=b)  # Gamma mean=th2temp var=c
    lognumth2 <- sum(log(falp(alptemp,p1temp,1-p1temp,th1temp,th2cand,
                              d)),na.rm=TRUE) + (A2-1)*log(th2cand)-th2cand/B2
    logalp <- lognumth2 - logdenth2 + 
      lgamma(th2temp^2/c)+(th2temp^2/c)*log(c/th2temp)-
      lgamma(th2cand^2/c)-(th2cand^2/c)*log(c/th2cand) +
      (th2cand^2/c-1)*log(th2temp) - (th2temp^2/c-1)*log(th2cand)
    logU <- log(runif(1,0,1))
    
    if(is.nan(logalp)==TRUE|is.infinite(logalp)==TRUE){
      eflag <- 1
    }
    
    if(eflag){
      next
    }
    
    if(logU < logalp){
      th2temp <- th2cand
      logdenth2 <- lognumth2
    }
    
    if(k > nburn && ((k - nburn) %% thin) == 0){
      storeIndex <- storeIndex + 1
      alpsamp[storeIndex] <- alptemp
      p1samp[storeIndex] <- p1temp
      th1samp[storeIndex] <- th1temp
      th2samp[storeIndex] <- th2temp
    }
  }
  
  # Compute posterior means and credible intervals
  alpha_mean <- mean(alpsamp)
  alpha_ci <- quantile(alpsamp, probs = c(0.025, 0.975))
  
  p_means <- mean(p1samp)
  p_ci <- quantile(p1samp, probs=c(0.025, 0.975))
  
  th1_mean <- mean(th1samp)
  th1_ci <- quantile(th1samp, probs = c(0.025, 0.975))
  
  th2_mean <- mean(th2samp)
  th2_ci <- quantile(th2samp, probs = c(0.025, 0.975))

  
  # Organize results into a structured list
  result <- list(
    posterior_means = list(
      alpha = alpha_mean,
      p1 = p_means,
      theta1 = th1_mean,
      theta2 = th2_mean
    ),
    credible_intervals = list(
      alpha = alpha_ci,
      p1 = p_ci,
      theta1 = th1_ci,
      theta2 = th2_ci
    )
  )
}