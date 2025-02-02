summary.alpmixBayes <- function(object, ...) {
  # Map short model names to full descriptions
  model_names <- list(
    "WW"  = "Weibull Weibull Mixture",
    "EW"  = "Exponential Weibull Mixture",
    "LL"  = "Lognormal Lognormal Mixture",
    "EWG" = "Exponential Weibull Gamma Mixture"
  )
  
  # Get the full name for the model
  full_model_name <- model_names[[object$model]]
  
  cat("Summary of alpmix model fitting:\n\n")
  cat("Model used:", full_model_name, "\n\n")
  
  cat("Parameter Estimates (with 95% Credible Intervals):\n")
  if (!is.null(object$estimates) && nrow(object$estimates) > 0) {
    print(object$estimates)
  } else {
    cat("No parameter estimates found.\n")
  }
  
  # Optionally print hyperparameters
  if (!is.null(object$hyper_values) && nrow(object$hyper_values) > 0) {
    cat("\nHyperparameters Used:\n")
    print(object$hyper_values)
  }
  
  invisible(object)
}