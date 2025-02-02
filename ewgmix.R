ewgmix <- function(d,init_values,mcmc_values,prior){
  
  ## likelihood 
  falp <- function(alp,p1,p2,p3,th1,th2,th3,d){
    
    S1 <- exp(-th1*d)
    S2 <- exp(-(d^2)/(th2^2))
    S3 <- 1-pgamma(d,2,scale=th3)
    h1 <- th1
    h2 <- (2*d)/(th2^2)
    h3 <- dgamma(d,2,scale=th3)/S3
    
    if(alp==0){
      
      Falp <- (S1^p1)*(S2^p2)*(S3^p3)
      Fr <-p1*h1+p2*h2+p3*h3
      
    } else{
      
      Falp <- (p1*(S1^(alp)) + p2*(S2^(alp)) + p3*(S3^(alp)))^(1/alp-1)
      Fr <- p1*h1*(S1^(alp)) + p2*h2*(S2^(alp)) + p3*h3*(S3^(alp))
    }
    Sur <- Falp*Fr
    return(Sur)
  }
  
  
  ## unlist mcmc 
  nburn = mcmc_values$nburn
  niter = mcmc_values$niter
  thin = mcmc_values$thin
  aiter <- nburn + niter
  nsamp <- ceiling(niter/thin)
  
  ## unlist the initial values
  alptemp = init_values$alpha
  p1temp = init_values$p1; p2temp = init_values$p2
  p3temp = 1-p1temp-p2temp
  ptemp = matrix(c(p1temp,p2temp,p3temp),1,3)
  th1temp=init_values$th1; th2temp=init_values$th2; th3temp=init_values$th3
  
  alpsamp <- rep(0,nsamp)
  psamp <- matrix(rep(0,nsamp*3),nrow=nsamp,ncol=3)
  th1samp <- rep(0,nsamp); th2samp <- rep(0,nsamp); th3samp <- rep(0,nsamp)
  
  # hyperparameters
  mualp <- prior$ap[1]; sig2alp <- prior$ap[2]
  c1 <- prior$pp[1]; c2 <- prior$pp[2]; c3 <- prior$pp[3]
  A1 <- prior$th1p[1]; B1 <- prior$th1p[2] 
  A2 <- prior$th2p[1]; B2 <- prior$th2p[2] 
  A3 <- prior$th3p[1]; B3 <- prior$th3p[2]
  
  logdenalp <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2temp,th3temp,d))) -
    0.5*(alptemp-mualp)^2/sig2alp
  logdenp <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2temp,th3temp,d))) +
    (c1-1)*log(ptemp[1]) + (c2-1)*log(ptemp[2]) + (c3-1)*log(ptemp[3])
  logdenth1 <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2temp,th3temp,d))) +
    (A1-1)*log(th1temp) - th1temp/B1
  logdenth2 <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2temp,th3temp,d))) +
    (A2-1)*log(th2temp) - th2temp/B2
  logdenth3 <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2temp,th3temp,d))) +
    (A3-1)*log(th3temp) - th3temp/B3
  
  storeIndex <- 0
  
  for(k in 1:aiter){
    # alpha
    # point mass
    pzero <- exp(sum(log(0.01002004*falp(alptemp,ptemp[1],ptemp[2],ptemp[3],
                                         th1temp,th2temp,th3temp,d)),na.rm = TRUE)) # gap of 500 points between 0 and 5 = 0.01002004
    u <- runif(1,0,.5)
    if(u < pzero){
      alptemp <- 0
    } else{
      alpcand <- rnorm(1,alptemp,.1)
      lognumalp <- sum(log(falp(alpcand,ptemp[1],ptemp[2],ptemp[3],
                                th1temp,th2temp,th3temp,d)),na.rm = TRUE)-0.5*(alpcand-mualp)^2/sig2alp
      logalp <- lognumalp - logdenalp
      logU <- log(runif(1,0,1))
      if(logU < logalp){
        alptemp <- alpcand
        logdenalp <- lognumalp
      }
    }
    
    # p=(p1,p2)
    proposal_alpha <- 0.7  # Tune this for broader/narrower proposals
    pcand <- rdirichlet(1, rep(proposal_alpha, 3)) 
    
    # Compute proposal densities
    log_q_current <- log(ddirichlet(ptemp, rep(proposal_alpha, 3)))  # q(current | proposed)
    log_q_proposed <- log(ddirichlet(pcand, rep(proposal_alpha, 3))) # q(proposed | current)
    log_proposal_ratio <- log_q_current - log_q_proposed
    
    lognump <- sum(log(falp(alptemp,pcand[1],pcand[2],pcand[3],th1temp,th2temp,th3temp,d)),na.rm = TRUE) +
      (c1-1)*log(pcand[1]) + (c2-1)*log(pcand[2]) + (c3-1)*log(pcand[3])
    logalp <- (lognump - logdenp) + log_proposal_ratio
    logU <- log(runif(1,0,1))
    if(logU < logalp){
      ptemp <- pcand
      logdenp <- lognump
    }
    
    # theta1: with Gamma prior mean=th1t and var=c
    c <- .1
    a <- th1temp^2/c
    b <- c/th1temp
    th1cand <- rgamma(1,a,scale=b)  # Gamma mean=th1temp var=c
    #th1cand <- rlnorm(1, th1temp, sdlog = 1)
    lognumth1 <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1cand,th2temp,th3temp,
                              d)),na.rm=TRUE) + (A1-1)*log(th1cand)-th1cand/B1
    
    # Calculate proposal ratio using dgamma()
    log_q_current <- dgamma(th1temp, shape = th1cand^2/c, scale = c/th1cand, log = TRUE)
    log_q_proposed <- dgamma(th1cand, shape = a, scale = b, log = TRUE)
    log_proposal_ratio <- log_q_current - log_q_proposed
    
    logalp <- lognumth1 - logdenth1 + lognumth1 - logdenth1 + log_proposal_ratio

    logU <- log(runif(1,0,1))
    if(logU < logalp){
      th1temp <- th1cand
      logdenth1 <- lognumth1
    }
    
    # theta2: with Gamma prior mean=th2t and var=c
    c <- .1
    a <- th2temp^2/c
    b <- c/th2temp
    th2cand <- rgamma(1,a,scale=b)  # Gamma mean=th2temp var=c
    #th2cand <- rlnorm(1, th2temp, sdlog = 1)
    lognumth2 <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2cand,th3temp,
                              d)),na.rm=TRUE) + (A2-1)*log(th2cand)-th2cand/B2
    # Calculate proposal ratio using dgamma()
    log_q_current <- dgamma(th2temp, shape = th2cand^2/c, scale = c/th2cand, log = TRUE)
    log_q_proposed <- dgamma(th2cand, shape = a, scale = b, log = TRUE)
    log_proposal_ratio <- log_q_current - log_q_proposed
    
    logalp <- lognumth2 - logdenth2 + log_q_current - log_q_proposed
    logU <- log(runif(1,0,1))
    if(logU < logalp){
      th2temp <- th2cand
      logdenth2 <- lognumth2
    }
    
    # theta3: with Gamma prior mean=th2t and var=c
    c <- .1
    a <- th3temp^2/c
    b <- c/th3temp
    th3cand <- rgamma(1,a,scale=b)  # Gamma mean=th2temp var=c
    #th3cand <- rlnorm(1, th3temp, sdlog = 1)
    lognumth3 <- sum(log(falp(alptemp,ptemp[1],ptemp[2],ptemp[3],th1temp,th2temp,th3cand,
                              d)),na.rm=TRUE) + (A3-1)*log(th3cand)-th3cand/B3
    log_q_current <- dgamma(th3temp, shape = th3cand^2/c, scale = c/th3cand, log = TRUE)
    log_q_proposed <- dgamma(th3cand, shape = a, scale = b, log = TRUE)
    log_proposal_ratio <- log_q_current - log_q_proposed
    
    
    logalp <- lognumth3 - logdenth3 + log_q_current - log_q_proposed
    logU <- log(runif(1,0,1))
    if(logU < logalp){
      th3temp <- th3cand
      logdenth3 <- lognumth3
    }
    
    if(k > nburn && ((k - nburn) %% thin) == 0){
      storeIndex <- storeIndex + 1
      alpsamp[storeIndex] <- alptemp
      psamp[storeIndex,] <- ptemp
      th1samp[storeIndex] <- th1temp
      th2samp[storeIndex] <- th2temp
      th3samp[storeIndex] <- th3temp
      
    }
  }
  
  # Compute posterior means and credible intervals
  alpha_mean <- mean(alpsamp)
  alpha_ci <- quantile(alpsamp, probs = c(0.025, 0.975))
  
  p_means <- colMeans(psamp)
  p_ci <- apply(psamp, 2, quantile, probs = c(0.025, 0.975))  
  
  th1_mean <- mean(th1samp)
  th1_ci <- quantile(th1samp, probs = c(0.025, 0.975))
  
  th2_mean <- mean(th2samp)
  th2_ci <- quantile(th2samp, probs = c(0.025, 0.975))
  
  th3_mean <- mean(th3samp)
  th3_ci <- quantile(th3samp, probs = c(0.025, 0.975))
  
  # Organize results into a structured list
  result <- list(
    posterior_means = list(
      alpha = alpha_mean,
      p1 = p_means[1],
      p2 = p_means[2],
      p3 = p_means[3],
      theta1 = th1_mean,
      theta2 = th2_mean,
      theta3 = th3_mean
    ),
    credible_intervals = list(
      alpha = alpha_ci,
      p1 = p_ci[, 1],
      p2 = p_ci[, 2],
      p3 = p_ci[, 3],
      theta1 = th1_ci,
      theta2 = th2_ci,
      theta3 = th3_ci
    )
  )
}