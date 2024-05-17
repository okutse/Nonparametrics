
gelmanRubin <- function(mcmc_chains) {
  n_chains <- length(mcmc_chains)
  n_iter <- nrow(mcmc_chains[[1]])
  n_params <- ncol(mcmc_chains[[1]])
  
  # Calculate means and variances for each chain using lapply
  chain_means <- t(sapply(mcmc_chains, colMeans))
  within_variances <- t(sapply(mcmc_chains, function(x) apply(x, 2, var)))
  
  # Calculate between-chain variance and average within-chain variance
  b <- n_iter * apply(chain_means, 2, var)
  w <- rowMeans(within_variances)
  
  # Calculate weighted variance
  var_plus <- ((n_iter - 1) / n_iter) * w + mean(b) / n_iter
  
  # Calculate potential scale reduction factor
  psrf <- sqrt(var_plus / w)
  
  return(psrf)
}



# gelmanRubin <- function(mcmc_chains) {
#   n_chains <- length(mcmc_chains)
#   n_iter <- nrow(mcmc_chains[[1]])
#   n_params <- ncol(mcmc_chains[[1]])
#   
#   # Calculate means and variances for each chain using lapply
#   chain_means <- t(sapply(mcmc_chains, colMeans))
#   within_variances <- t(sapply(mcmc_chains, function(x) apply(x, 2, var)))
#   
#   # Calculate between-chain variance and average within-chain variance
#   b <- n_iter * apply(chain_means, 2, var)
#   w <- rowMeans(within_variances)
#   
#   # Calculate weighted variance
#   var_plus <- ((n_iter - 1) / n_iter) * w + mean(b) / n_iter
#   
#   # Calculate potential scale reduction factor
#   psrf <- sqrt(var_plus / w)
#   
#   # Calculate the 95% confidence interval for the psrf
#   R2_fixed = (n_iter - 1)/n_iter
#   
#   #psrf_ci <- R2_fixed + qf(0.975, df1 = (n_chains - 1), df2 = diag(w)) / (2*n_chains) * var_plus / w
#   
#   # Create a data frame to store the results
#   results <- data.frame(
#     Parameter = colnames(mcmc_chains[[1]]),
#     PSRF = psrf,
#     #Upper95CI = psrf_ci
#   )
#   
#   return(results)
# }













compute_gelman_rubin <- function(mcmc_list) {
  # Get the number of chains and iterations
  num_chains <- length(mcmc_list)
  num_iterations <- dim(mcmc_list[[1]])[1] #length(mcmc_list[[1]])
  
  # Initialize matrices to store within-chain and between-chain variances
  W <- matrix(0, nrow = num_chains, ncol = num_iterations)
  #B <- matrix(0, nrow = num_chains, ncol = num_iterations)
  B <- 0
  
  # Compute within-chain and between-chain variances
  for (i in 1:num_chains) {
    chain <- mcmc_list[[i]]
    chain_mean <- colMeans(chain)
    W[i, ] <- rowSums((chain - chain_mean)^2) / (num_iterations - 1)
  }
  
  chain_means <- colMeans(do.call(rbind, mcmc_list))
  overall_mean <- mean(chain_means)
  
  #B <- num_iterations * rowSums((chain_means - overall_mean)^2) / (num_chains - 1)
  # Compute squared differences between chain means and overall mean
  for (i in 1:num_chains) {
    B <- B + num_iterations * sum((chain_means - overall_mean)^2)
  }
  # Divide by (num_chains - 1) to get between-chain variance
  B <- B / (num_chains - 1)
  
  
  # Compute the potential scale reduction factor (R-hat)
  var_within <- mean(W)
  var_between <- mean(B)
  psrf <- sqrt((var_within + var_between) / var_within)
  
  # Compute the upper confidence limit for R-hat
  df <- num_chains * (num_iterations - 1)
  upper_conf_limit <- sqrt(qf(0.975, df, df))
  
  cat("Point estimate of R-hat:", round(psrf, 4), "\n")
  cat("Upper confidence limit for R-hat:", round(upper_conf_limit, 4), "\n")
}



library(ggplot2)

plot_trace_and_density <- function(mcmc_list) {
  num_chains <- length(mcmc_list)
  
  # Create a separate plot for each parameter
  for (param_name in names(mcmc_list[[1]])) {
    param_data <- data.frame(chain = rep(1:num_chains, each = length(mcmc_list[[1]])),
                             iteration = rep(1:length(mcmc_list[[1]]), times = num_chains),
                             value = unlist(lapply(mcmc_list, function(chain) chain[[param_name]])))
    
    # Trace plot
    trace_plot <- ggplot(param_data, aes(x = iteration, y = value, color = factor(chain))) +
      geom_line() +
      labs(title = paste("Trace Plot for", param_name),
           x = "Iteration",
           y = "Parameter Value") +
      theme_bw()
    
    # Density plot
    density_plot <- ggplot(param_data, aes(x = value, fill = factor(chain))) +
      geom_density(alpha = 0.5) +
      labs(title = paste("Density Plot for", param_name),
           x = "Parameter Value",
           y = "Density") +
      theme_bw()
    
    # Print both plots
    print(trace_plot)
    print(density_plot)
  }
  print(trace_plot)
  print(density_plot)
}

# Example usage:
# Assuming you have an mcmc.list object called 'my_chains'
plot_trace_and_density(my_chains)

