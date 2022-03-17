if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse)

k <- 10
M <- 5000
N <- 1000

df <- tibble(i = 1:k, 
             theta_1 = i, 
             theta_2 = 10*i, 
             theta_3 = i / 10, 
             sd = 1)

draw_theta_hat <- function(theta, k){
  rnorm(k, theta, 1)
}

calculate_ranks <- function(n, theta_hat, k){
  tibble(n = n, i = 1:k, ranks = rank(theta_hat))
}

m <- 1
theta_var <- "theta_1"

run_bs_draw <- function(m, theta_var, df, k, N)
  {
  # Get an initial draw of theta hat. 
  theta_hat <- draw_theta_hat(df[[theta_var]], k)
  r1        <- calculate_ranks(n = 0, theta_hat, k) %>% 
    select(i, r1 = ranks)
  
  # Parametric BS: draw N times from normals centered on theta_hat
  bs <- map_dfr(1:N, 
                function(n){
                  tt <- draw_theta_hat(theta_hat, k = k)
                  calculate_ranks(n, tt, k = k)
                  }
                ) %>% 
    left_join(r1) %>% 
    mutate(RR = ranks - r1)
  
  # Get quantiles 
  qq <- bs %>% 
    group_by(i) %>% 
    summarize(q.025 = quantile(RR, probs = 0.025), 
              q.975 = quantile(RR, probs = 0.975)) %>% 
    left_join(r1) %>% 
    left_join(df) %>% 
    mutate(lower_percentile = r1 - q.975, upper_percentile = r1 - q.025, 
           lower_effron =     r1 + q.025, upper_effron     = r1 + q.975) %>% 
    mutate(percentile_covered = 
             ifelse(.data[[theta_var]] >= lower_percentile & 
                    .data[[theta_var]] <= upper_percentile, 1, 0)) %>% 
    mutate(effron_covered = 
             ifelse(.data[[theta_var]] >= lower_effron & 
                    .data[[theta_var]] <= upper_effron, 1, 0)) %>% 
    select(i, percentile_covered, effron_covered)
  qq %>% 
    mutate(m = !!m)
}
run_sim <- function(theta_var, df, k, N){
  map_dfr(1:M, run_bs_draw, "theta_1", df, k, N)
}




for ( m in 1:M ) {
  for (  ) {
    
  }
}