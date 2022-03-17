if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse)


# set parameters ----------------------------------------------------------

k <- 10
M <- 5000
N <- 1000
df <- tibble(i = 1:k, 
             theta_1 = i, 
             theta_2 = 10 * i, 
             theta_3 = i / 10)

# functions ---------------------------------------------------------------

draw_theta_hat <- function(theta, k){
  rnorm(k, theta, 1)
}

run_bs_draw <- function(m, theta_var, df, k, N)
  {
  # Get an initial draw of theta hat and their ranks
  theta_hat <- draw_theta_hat(df[[theta_var]], k)
  r1        <- matrix(rank(theta_hat), nrow = 1)
  
  # Draw and rank all N normal draws of dim k. subtract r1 from each draw
  norms <- matrix(nrow = N, ncol = k)
  ranks <- norms
  for(i in 1:k){
    norms[, i] <- rnorm(N, mean = theta_hat[i])
  }
  for(i in 1:N){
    ranks[i,] <- rank(norms[i,]) - r1
  }
  # Get quantiles
  map_dfr(1:k, function(kk) quantile(ranks[,kk], c(0.025, 0.975))) %>% 
    mutate(i = 1:k) %>% 
    select(i, q.025 = `2.5%`, q.975 = `97.5%`) %>% 
    mutate(theta = df[[theta_var]], r1 = drop(r1)) %>%
    mutate(lower_percentile = r1 - q.975, upper_percentile = r1 - q.025, 
           lower_effron =     r1 + q.025, upper_effron     = r1 + q.975) %>% 
    mutate(percentile_covered = 
             ifelse(theta >= lower_percentile & theta <= upper_percentile, 1, 0)) %>% 
    mutate(effron_covered = 
             ifelse(theta >= lower_effron & theta <= upper_effron, 1, 0)) %>% 
    select(i, percentile_covered, effron_covered) %>%
    mutate(m = !!m)
}

run_sim <- function(M, theta_var, df, k, N){
  map_dfr(1:M, run_bs_draw, theta_var, df, k, N) %>% 
    group_by(i) %>% summarize(percentile_covered = mean(percentile_covered), 
                              effron_covered     = mean(effron_covered)) %>% 
    mutate(theta = !!theta_var) %>% 
    relocate(theta)
}

# run code ----------------------------------------------------------------
out <- run_sim(M = 1000, theta_var = "theta_1", df = df, k = k, N =N) %>% 
  bind_rows(
    run_sim(M = 1000, theta_var = "theta_2", df = df, k = k, N =N) 
  ) %>% 
  bind_rows(
    run_sim(M = 1000, theta_var = "theta_3", df = df, k = k, N =N) 
  )

