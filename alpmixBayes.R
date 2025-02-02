# Source all the individual function files
source("wwmix.R")  # Load the wwmix function
source("ewmix.R")  # Load the ewmix function
source("llmix.R")  # Load the llmix function
source("ewgmix.R") # Load the ewgmix function
source("summary.alpmixBayes.R") ## Load the S3

alpmixBayes <- function(d,mcmc_values=NULL,
                        init_values=NULL,
                        prior=NULL,
                        survmodel=c("WW", "EW", "LL", "EWG"),
                        ...){
  
  # 1) Ensure all four sub-functions are loaded
  required_functions <- c("wwmix", "ewmix", "llmix", "ewgmix")
  missing_functions  <- required_functions[!sapply(required_functions, exists)]
  if (length(missing_functions) > 0) {
    stop(
      "The following functions are missing: ",
      paste(missing_functions, collapse = ", ")
    )
  }
  
  # 2) Match the user-chosen model
  model <- match.arg(survmodel)
  
  # 3) Set default initial values if none are provided
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
  }
  
  # 4) Set default hyperparameter values dynamically
  if (is.null(prior)) {
    if (model != "EWG") {
      # Two-component models
      prior <- list(
        ap = c(0,0.001), 
        pp = c(1,1),
        th1p = c(1,1),
        th2p = c(1,1)
      )
    } else {
      # Three-component model (EWG)
      prior <- list(
        ap = c(0,0.001), 
        pp = c(1,1,1),
        th1p = c(1,1),
        th2p = c(1,1),
        th3p = c(1,1)
      )
    }
  }
  
  # 5) If user didn't provide mcmc_values, set default MCMC controls
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
  
  # 6) Define a small helper for error handling
  error_result <- function(model_name, err_msg) {
    cat("\n[WARNING] Error in model:", model_name, "\n")
    cat("  -> The data may not be suitable. Details:\n")
    cat("  ->", err_msg, "\n\n")
    list(
      error         = TRUE,
      errorMsg      = err_msg,
      posterior_means = NA
      # Additional placeholders as needed
    )
  }
  
  # 7) Call the chosen model function with `tryCatch`
  fit_out <- tryCatch(
    {
      switch(
        model,
        "WW" = wwmix(d, init_values, mcmc_values, prior),
        "EW" = ewmix(d, init_values, mcmc_values, prior),
        "LL" = llmix(d, init_values, mcmc_values, prior),
        "EWG"= ewgmix(d, init_values, mcmc_values, prior)
      )
    },
    error = function(e) {
      error_result(model, e$message)
    }
  )
  
  # 8) Put result in a data frame:
  if (model != "EWG") {
    # This is a TWO-component model with parameters p, th1, th2
    alp_est <- fit_out$posterior_means$alpha
    alp_ci <- fit_out$credible_intervals$alpha
    
    p_est <- fit_out$posterior_means$p1
    p_ci <- fit_out$credible_intervals$p1
    
    th1_est <- fit_out$posterior_means$theta1
    th1_ci <- fit_out$credible_intervals$theta1
    
    th2_est <- fit_out$posterior_means$theta2
    th2_ci <- fit_out$credible_intervals$theta2
    
    # Combine them in vectors
    param_names  <- c("alpha", "p", "th1", "th2")
    param_est    <- c(alp_est, p_est, th1_est, th2_est)
    param_lower  <- c(alp_ci[1], p_ci[1], th1_ci[1], th2_ci[1])
    param_upper  <- c(alp_ci[2], p_ci[2], th1_ci[2], th2_ci[2])
    
  }else {
    # This is the THREE-component model ewgmix with p1, p2, th1, th2, th3
    alp_est <- fit_out$posterior_means$alpha
    alp_ci <- fit_out$credible_intervals$alpha
    
    p1_est <- fit_out$posterior_means$p1
    p1_ci <- fit_out$credible_intervals$p1
    
    p2_est <- fit_out$posterior_means$p2
    p2_ci <- fit_out$credible_intervals$p2
    
    th1_est <- fit_out$posterior_means$theta1
    th1_ci <- fit_out$credible_intervals$theta1
    
    th2_est <- fit_out$posterior_means$theta2
    th2_ci <- fit_out$credible_intervals$theta2
    
    th3_est <- fit_out$posterior_means$theta3
    th3_ci <- fit_out$credible_intervals$theta3
    
    # Combine them in vectors
    param_names  <- c("alpha", "p1", "p2", "th1", "th2", "th3")
    param_est    <- c(alp_est, p1_est, p2_est, th1_est, th2_est, th3_est)
    param_lower  <- c(alp_ci[1], p1_ci[1], p2_ci[1], th1_ci[1], th2_ci[1], th3_ci[1])
    param_upper  <- c(alp_ci[2], p1_ci[2], p2_ci[2], th1_ci[2], th2_ci[2], th3_ci[2])
  }
  # Now create a single data frame:
  estimates_df <- data.frame(
    Parameter = param_names,
    Estimate  = param_est,
    Lower_95  = param_lower,
    Upper_95  = param_upper
  )
  
  
  # 9) Build your return object
  result <- list(
    model = model,
    estimates = estimates_df,
    data = d,
    call = match.call()
  )
  
  # 10) Assign an S3 class for custom print/summary methods
  class(result) <- "alpmixBayes"
  
  return(result)
  
}