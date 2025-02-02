rm(list=ls())


mcmc_values <- list(nburn=1000,niter=5000,thin=5)
init_values <- list(alpha=1,p = 0.5, p1 = 0.3,p2 = 0.3,
                    th1 = 1.5,th2 = 1.5,th3 = 2 ) 
prior <- list(ap=c(0,0.001),pp=c(1,1),
              th1p=c(1,1),th2p=c(1,1))
prior <- list(ap=c(1,0.001),pp=c(1,1,1),
              th1p=c(1,1),th2p=c(1,1),th3p=c(1,1))

dir_sim <- "~/Desktop/Alpha Mixture/alpmix"
df <- get(load(paste(dir_sim, "ewg.100.Rdata", sep='/')))
data <- df[[1]]; d <- data$y
source("alpmix.R")
source("alpmixBayes.R")
ptm <- proc.time()
r <- alpmixBayes(d,mcmc_values,init_values,prior,survmodel = "EWG")
(proc.time() - ptm)
summary(r)
