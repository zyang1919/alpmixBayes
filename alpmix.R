# Source all the individual function files
source("wwmix.R")  # Load the wwmix function
source("ewmix.R")  # Load the ewmix function
source("llmix.R")  # Load the llmix function
source("ewgmix.R") # Load the ewgmix function
source("print.combinedMix.R") ## Load the S3

alpmix <- function(data, mcmc_values, init_values, plot_surv=TRUE,...){
  
  # Ensure all functions are loaded
  required_functions <- c("wwmix", "ewmix", "llmix", "ewgmix")
  missing_functions <- required_functions[!sapply(required_functions, exists)]
  
  if (length(missing_functions) > 0) {
    stop(paste("The following functions are missing:", 
               paste(missing_functions, collapse = ", ")))
  }
  
  # Set default initial values if none are provided
  if (is.null(init_values)) {
    init_values <- list(
      alpha=1,
      p = 0.5,        # Initial for 2-component models (wwmix, ewmix, llmix)
      p1 = 0.3,       # Initial values for 3-component model (ewgmix)
      p2 = 0.3,
      th1 = 2,     # Initial shape parameters for components
      th2 = 2,
      th3 = 2    # Only used in the 3-component model
    )
  } else {
    # Validate initial values for two-component models
    required_2comp <- c("p", "th1", "th2")
    if (!all(required_2comp %in% names(init_values))) {
      stop("Initial values for two-component models must 
           include 'p', 'th1', 'th2'.")
    }
    # Validate initial values for three-component models
    required_3comp <- c("p1", "p2", "th1", "th2", "th3")
    if (!all(required_3comp %in% names(init_values))) {
      stop("Initial values for two-component models must 
           include 'p', 'th1', 'th2', 'th3'.")
    }
  }
  
  # If user didn't provide mcmc_values, set default MCMC controls
  if (is.null(mcmc_values)) {
    mcmc_values <- list(
      nburn = 1000,   # Number of burn-in samples
      niter = 5000,   # Total iterations
      thin   = 1       # Thinning interval: collect every sample
    )
  } else {
    # Optional: you could add checks to ensure mcmc_values has the fields you expect
    required_mcmc_fields <- c("nburn", "niter", "thin")
    missing_mcmc_fields <- required_mcmc_fields[!required_mcmc_fields %in% names(mcmc_values)]
    if (length(missing_mcmc_fields) > 0) {
      stop(
        "Missing the following MCMC control fields: ",
        paste(missing_mcmc_fields, collapse = ", ")
      )
    }
  }
  
  # Helper function to return a placeholder if a model fails
  error_result <- function(model_name, err_msg) {
    cat("\n[WARNING] ", model_name, " encountered an error:\n", sep = "")
    cat("  -> The data may not be suitable for ", model_name, ". Error details:\n", sep = "")
    cat("  ->", err_msg, "\n\n")
    # Return a placeholder list indicating no valid result
    list(
      error      = TRUE,
      errorMsg   = err_msg,
      # You can store NA or NULL for each item you normally return
      posterior_means = NA,
      survival   = NA
    )
  }
  
  
  # Define survival time for the plot
  sur_time <- seq(0,10,length = 100)
  
  ##############################################################################
  # 1) W+W mixture
  ##############################################################################
  wwmix_out <- tryCatch(
    {
      # Call wwmix; if it fails, we jump to error handler
      wwmix(d,init_values,mcmc_values)
    },
    error = function(e) {
      error_result("wwmix", e$message)
    }
  )
  wwpar <- wwmix_out$posterior_means
  
  if(!isTRUE(wwmix_out$error)){
    # wwmix ran successfully, so we attempt to compute the survival curve
    wwsurv <- function(t,parlist){
      ## unlist parameters
      alphat <- parlist$alpha;p1hat <- parlist$p1
      th1hat <- parlist$theta1;th2hat <- parlist$theta2
      
      sur <- (p1hat*exp(-(alphat*(t^2))/(th1hat^2))+
                (1-p1hat)*exp(-(alphat*(t^2))/(th2hat^2)))^((1/alphat))
      return(sur)
    }
    
    surv_wwmix <- wwsurv(sur_time,wwpar)
  }
  
  ##############################################################################
  # 2) E+W mixture
  ##############################################################################
  ewmix_out <- tryCatch(
    {
      # Call wwmix; if it fails, we jump to error handler
      ewmix(d,init_values,mcmc_values)
    },
    error = function(e) {
      error_result("ewmix", e$message)
    }
  )
  ewpar <- ewmix_out$posterior_means
  
  if(!isTRUE(ewmix_out$error)){
    # ewmix ran successfully, so we attempt to compute the survival curve
    ewsurv <- function(t,parlist){
      ## unlist parameters
      alphat <- parlist$alpha;p1hat <- parlist$p1
      th1hat <- parlist$theta1;th2hat <- parlist$theta2
      
      sur <- (p1hat*exp(-(alphat*t)/(th1hat))+
                (1-p1hat)*exp(-alphat*((t/th2hat)^(0.5))))^((1/alphat))
      return(sur)
    }
    
    surv_ewmix <- ewsurv(sur_time,ewpar)
  }
  
  ##############################################################################
  # 3) L+L mixture
  ##############################################################################
  llmix_out <- tryCatch(
    {
      # Call wwmix; if it fails, we jump to error handler
      llmix(d,init_values,mcmc_values)
    },
    error = function(e) {
      error_result("llmix", e$message)
    }
  )
  llpar <- llmix_out$posterior_means
  
  if(!isTRUE(llmix_out$error)){
    # ewmix ran successfully, so we attempt to compute the survival curve
    llsurv <- function(t,parlist){
      ## unlist parameters
      alphat <- parlist$alpha;p1hat <- parlist$p1
      mu1hat <- parlist$mu1;mu2hat <- parlist$mu2
      
      sur1 <- 1 - plnorm(t, meanlog = mu1hat, sdlog = 0.5)
      sur2 <- 1 - plnorm(t, meanlog = mu2hat, sdlog = 0.1)
      
      sur <- (p1hat*(sur1^(alphat))+(1-p1hat)*(sur2^(alphat)))^((1/alphat))
      return(sur)
    }
    
    surv_llmix <- llsurv(sur_time,llpar)
  }
  
  ##############################################################################
  # 4) E+W+G mixture
  ##############################################################################
  ewgmix_out <- tryCatch(
    {
      # Call wwmix; if it fails, we jump to error handler
      ewgmix(d,init_values,mcmc_values)
    },
    error = function(e) {
      error_result("ewgmix", e$message)
    }
  )
  ewgpar <- ewgmix_out$posterior_means
  
  if(!isTRUE(ewgmix_out$error)){
    # ewmix ran successfully, so we attempt to compute the survival curve
    ewgsurv <- function(t,parlist){
      ## unlist parameters
      alphat <- parlist$alpha
      p1hat <- parlist$p1;p2hat <- parlist$p2;p3hat <-parlist$p3
      th1hat <- parlist$theta1;th2hat <- parlist$theta2;th3hat<-parlist$theta3
      
      sur1 <- exp(-th1hat*t)
      sur2 <- exp(-(t^2)/(th2hat^2))
      sur3 <- 1-pgamma(t, 2,scale=th3hat)
      
      sur <- (p1hat*(sur1^(alphat))+p2hat*(sur2^(alphat))+p3hat*(sur3^(alphat)))^((1/alphat))
      return(sur)
    }
    
    surv_ewgmix <- ewgsurv(sur_time,ewgpar)
  }
  
  ######################################################################
  # 3) Optionally plot all four lines on one plot (plot_surv = TRUE/FALSE)
  ######################################################################
  if (plot_surv) {
    ## Open a single blank plot covering the time range
    plot(x = sur_time,y = rep(NA, length(sur_time)), 
      type = "n",xlab = "Time",ylab = "Survival Probability", 
      ylim = c(0,1), main = "Mixture Survival Curves")
    
    # We'll keep track of each model's legend label
    legend_labels <- character(4)
    line_colors   <- c("blue","green","red","purple")  # for 4 models in order
    
    ## wwmix line (blue)
    if (isTRUE(wwmix_out$error)) {
      # Means error in wwmix or survival
      legend_labels[1] <- "WWMixture not work"
    } else {
      lines(sur_time, surv_wwmix, col = "blue", lwd = 2)
      legend_labels[1] <- "WWMixture"
    }
    
    ## ewmix line (green)
    if (isTRUE(ewmix_out$error)) {
      legend_labels[2] <- "EWMixture not work"
    } else {
      lines(sur_time, surv_ewmix, col = "green", lwd = 2)
      legend_labels[2] <- "EWMixture"
    }
    
    # 3D) llmix line (red)
    if (isTRUE(llmix_out$error)) {
      legend_labels[3] <- "LLMixture not work"
    } else {
      lines(sur_time, surv_llmix, col = "red", lwd = 2)
      legend_labels[3] <- "LLMixture"
    }
    
    # 3E) ewgmix line (purple)
    if (isTRUE(ewgmix_out$error)) {
      legend_labels[4] <- "EWGMixture not work"
    } else {
      lines(sur_time, surv_ewgmix, col = "purple", lwd = 2)
      legend_labels[4] <- "EWGMixture"
    }
    
    # 3F) Add a legend
    legend("topright",legend = legend_labels,
      col    = c("blue","green","red","purple"),
      lty    = 1,lwd    = 2,bty    = "n")
  }

  ## result saving 
  out <- list(wwmix = wwmix_out,ewmix = ewmix_out,
              llmix = llmix_out,ewgmix = ewgmix_out,
              call = match.call(),data = d)
  
  class(out) <- "combinedMix"
  return(out)
  
}

