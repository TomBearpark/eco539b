if(!require(pacman)) install.packages('pacman')
pacman::p_load(tidyverse, xtable, MASS, furrr, matrixStats)
theme_set(theme_bw())
seed   <- 8894; set.seed(seed)
ncores <- 6
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps1/"

# functions ---------------------------------------------------------------

run_bs_draw <- function(m, theta_var, df, k, N, verbose = FALSE)
  {
  if(verbose) if(m %% 100 == 0) message(paste0("---- draw ", m, " ----"))
  # Get an initial draw of theta hat and their ranks
  theta_hat <- rnorm(k, df[[theta_var]], 1)
  r1        <- matrix(rank(theta_hat), nrow = 1)
  
  # Draw all N normal draws of dim k. 
  norms <- matrix(nrow = N, ncol = k)
  for(i in 1:k) norms[, i] <- rnorm(N, mean = theta_hat[i])
  
  ranks <- sweep(matrixStats::rowRanks(norms), 2, r1)

  # Get quantiles
  matrixStats::colQuantiles(ranks, probs = c(0.025, 0.975)) %>% as_tibble() %>% 
    mutate(i = 1:k, theta = df[[theta_var]], r1 = drop(r1)) %>%
    mutate(lower_percentile = r1 - `97.5%`, upper_percentile = r1 - `2.5%`, 
           lower_effron     = r1 + `2.5%`,  upper_effron     = r1 + `97.5%`) %>% 
    mutate(percentile_covered = 
             ifelse(i >= lower_percentile & i <= upper_percentile, 1, 0), 
           effron_covered = 
             ifelse(i >= lower_effron & i <= upper_effron, 1, 0)) %>% 
    dplyr::select(i, percentile_covered, effron_covered)
}

run_sim <- function(M, theta_var, k, N, seed = 8894){
  df <- tibble(i = 1:k, theta_1 = i, theta_2 = 10 * i, theta_3 = i / 10)
  message("running BS")
  future_map_dfr(1:M, run_bs_draw, theta_var = theta_var, df = df, k = k, N = N, 
                 verbose= FALSE,
                 .options = furrr_options(seed = seed), 
                 .progress = TRUE) %>% 
    group_by(i) %>% 
    summarize(percentile_coverage = mean(percentile_covered), 
              effron_coverage     = mean(effron_covered)) %>% 
    mutate(theta = !!theta_var) %>% 
    relocate(theta)
}

get_df <- function(M, k, N){
  out1 <- run_sim(M = M, theta_var = "theta_1", k = k, N = N) 
  out2 <- run_sim(M = M, theta_var = "theta_2", k = k, N = N) 
  out3 <- run_sim(M = M, theta_var = "theta_3", k = k, N = N) 

  bind_rows(out1, out2, out3)
}

table <- function(out_df){
  out_df %>% pivot_wider(id_cols = c(i),names_from = theta, 
                         values_from = c(percentile_coverage, effron_coverage)) %>% 
    relocate(i, percentile_coverage_theta_1, effron_coverage_theta_1, 
             percentile_coverage_theta_2, effron_coverage_theta_2, 
             percentile_coverage_theta_3, effron_coverage_theta_3) %>% 
    xtable() %>% print(include.rownames=FALSE)
}
plot <- function(out_df){
  out_df %>% 
    pivot_longer(cols = c(percentile_coverage, effron_coverage), names_to = "Method", 
                 values_to = "Coverage Probability") %>% 
    mutate(Method = str_remove(Method, "_coverage")) %>% 
    rename(Design = theta) %>% 
    mutate(Design = str_sub(Design, -1), Design = paste0("Design ", Design)) %>% 
    ggplot() + 
    geom_point(aes(x = i, y = `Coverage Probability`, color = Method)) + 
    geom_hline(aes(yintercept = 0.95)) + 
    facet_wrap(~Design)
}

# run code ----------------------------------------------------------------
k  <- 10
M  <- 5000
N  <- 1000

plan(multisession, workers = ncores)
out_df <- get_df(M, k, N)
table(out_df)
plot(out_df)
ggsave(paste0(out, "coverage_plot.png"), height = 3, width = 7)

# Version with k = 100
out_df2 <- get_df(M = M, k = 100, N = N)
plot(out_df2)
ggsave(paste0(out, "coverage_plotk_100.png"), height = 3, width = 7)
