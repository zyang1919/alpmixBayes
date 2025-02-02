print.combinedMix <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nResults for each mixture model:\n")
  
  # wwmix
  if (isTRUE(x$wwmix$error)) {
    cat("  - WWMixture: [ERROR]", x$wwmix$msg, "\n")
  } else {
    # Show some parameters (if they exist)
    cat("  - WWMixture:\n")
    cat("     alpha:  ", x$wwmix$posterior_means$alpha,  
        " (CI: [", x$wwmix$credible_intervals$alpha[1],  ", ", 
        x$wwmix$credible_intervals$alpha[2],  "])\n", sep="")
    cat("     p1:     ", x$wwmix$posterior_means$p1,     
        " (CI: [", x$wwmix$credible_intervals$p1[1],     ", ", 
        x$wwmix$credible_intervals$p1[2],     "])\n", sep="")
    cat("     theta1:     ", x$wwmix$posterior_means$theta1,     
        " (CI: [", x$wwmix$credible_intervals$theta1[1],     ", ", 
        x$wwmix$credible_intervals$theta1[2],     "])\n", sep="")
    cat("     theta2:     ", x$wwmix$posterior_means$theta2,     
        " (CI: [", x$wwmix$credible_intervals$theta2[1],     ", ", 
        x$wwmix$credible_intervals$theta2[2],     "])\n", sep="")
  }
  
  # ewmix
  if (isTRUE(x$ewmix$error)) {
    cat("  - EWMixture: [ERROR]", x$ewmix$msg, "\n")
  } else {
    # Show some parameters (if they exist)
    cat("  - EWMixture:\n")
    cat("     alpha:  ", x$ewmix$posterior_means$alpha,  
        " (CI: [", x$ewmix$credible_intervals$alpha[1],  ", ", 
        x$ewmix$credible_intervals$alpha[2],  "])\n", sep="")
    cat("     p1:     ", x$ewmix$posterior_means$p1,     
        " (CI: [", x$ewmix$credible_intervals$p1[1],     ", ", 
        x$ewmix$credible_intervals$p1[2],     "])\n", sep="")
    cat("     theta1:     ", x$ewmix$posterior_means$theta1,     
        " (CI: [", x$ewmix$credible_intervals$theta1[1],     ", ", 
        x$ewmix$credible_intervals$theta1[2],     "])\n", sep="")
    cat("     theta2:     ", x$ewmix$posterior_means$theta2,     
        " (CI: [", x$ewmix$credible_intervals$theta2[1],     ", ", 
        x$ewmix$credible_intervals$theta2[2],     "])\n", sep="")
  }
  
  # llmix
  if (isTRUE(x$llmix$error)) {
    cat("  - LLMixture: [ERROR]", x$llmix$msg, "\n")
  } else {
    # Show some parameters (if they exist)
    cat("  - LLMixture:\n")
    cat("     alpha:  ", x$llmix$posterior_means$alpha,  
        " (CI: [", x$llmix$credible_intervals$alpha[1],  ", ", 
        x$llmix$credible_intervals$alpha[2],  "])\n", sep="")
    cat("     p1:     ", x$llmix$posterior_means$p1,     
        " (CI: [", x$llmix$credible_intervals$p1[1],     ", ", 
        x$llmix$credible_intervals$p1[2],     "])\n", sep="")
    cat("     mu1:     ", x$llmix$posterior_means$mu1,     
        " (CI: [", x$llmix$credible_intervals$mu1[1],     ", ", 
        x$llmix$credible_intervals$mu1[2],     "])\n", sep="")
    cat("     mu2:     ", x$llmix$posterior_means$mu2,     
        " (CI: [", x$llmix$credible_intervals$mu2[1],     ", ", 
        x$llmix$credible_intervals$mu2[2],     "])\n", sep="")
  }
  
  ## ewgmix
  if (isTRUE(x$ewgmix$error)) {
    cat("  - EWGMixture: [ERROR]", x$ewgmix$msg, "\n")
  } else {
    # Show some parameters (if they exist)
    cat("  - EWGMixture:\n")
    cat("     alpha:  ", x$ewgmix$posterior_means$alpha,  
        " (CI: [", x$ewgmix$credible_intervals$alpha[1],  ", ", 
        x$ewgmix$credible_intervals$alpha[2],  "])\n", sep="")
    cat("     p1:     ", x$ewgmix$posterior_means$p1,     
        " (CI: [", x$ewgmix$credible_intervals$p1[1],     ", ", 
        x$ewgmix$credible_intervals$p1[2],     "])\n", sep="")
    cat("     p2:     ", x$ewgmix$posterior_means$p2,     
        " (CI: [", x$ewgmix$credible_intervals$p2[1],     ", ", 
        x$ewgmix$credible_intervals$p2[2],     "])\n", sep="")
    cat("     theta1:     ", x$ewgmix$posterior_means$theta1,     
        " (CI: [", x$ewgmix$credible_intervals$theta1[1],     ", ", 
        x$ewgmix$credible_intervals$theta1[2],     "])\n", sep="")
    cat("     theta2:     ", x$ewgmix$posterior_means$theta2,     
        " (CI: [", x$ewgmix$credible_intervals$theta2[1],     ", ", 
        x$ewgmix$credible_intervals$theta2[2],     "])\n", sep="")
    cat("     theta3:     ", x$ewgmix$posterior_means$theta3,     
        " (CI: [", x$ewgmix$credible_intervals$theta3[1],     ", ", 
        x$ewgmix$credible_intervals$theta3[2],     "])\n", sep="")
  }
  
  cat("\n")
  invisible(x)
}
