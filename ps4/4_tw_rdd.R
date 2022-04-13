pacman::p_load(tidyverse, fixest)
N <- 100
TT <- 100
beta <- 4


# -------------------------------------------------------------------------
# q1 ----------------------------------------------------------------------
# -------------------------------------------------------------------------

df <- expand_grid(i = 1:N, t = 1:TT) %>% 
  mutate(x =   rnorm(N*TT, sd = 5),
         eps = rnorm(N*TT)) %>% 
  group_by(i) %>% 
  mutate(alpha = rnorm(N, sd = 10)) %>% 
  ungroup() %>% 
    mutate(y = x*beta + alpha + eps) %>% 
  group_by(t) %>% 
    mutate(lambda = rnorm(TT, sd = 3) * t) %>% 
  ungroup() %>% 
  mutate(y2 = alpha + lambda + x*beta + eps)

df <- df %>% 
  group_by(i) %>% 
    mutate(xbar_n = mean(x)) %>% 
  ungroup() %>% 
    mutate(x_min_xbar_n = x - xbar_n) %>% 
  group_by(t) %>% 
    mutate(xbar_t = mean(x)) %>% 
  ungroup() %>%
  mutate(x_min_xbar_t = x - xbar_t, 
         x_min_both = x - xbar_t - xbar_n)  %>% 
  mutate(x_bar = mean(x), 
         x_min_both_and_full = x_min_both + x_bar)
  

# one way -----------------------------------------------------------------

feols(y ~ x | i, data = df)
feols(y ~ x + xbar_n, data = df)
feols(y ~ x_min_xbar_n, data = df)

# two way -----------------------------------------------------------------

feols(y ~ x | i + t, data = df)
feols(y ~ x + xbar_n + xbar_t + 1, data = df)
feols(y ~ x_min_both, data = df)
feols(y ~ -1 + x_min_both_and_full, data = df)

# -------------------------------------------------------------------------
# q2 ----------------------------------------------------------------------
# -------------------------------------------------------------------------