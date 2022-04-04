pacman::p_load(tidyverse, haven, fixest, furrr, 
               dfadjust, xtable, data.table, rsample)
theme_set(theme_bw())
seed <- 123
set.seed(seed)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps3/"
dir2 <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps2/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps3/"

plan(multisession, workers = 7)

# problem 1a ---------------------------------------------------------------
df  <- read_dta(paste0(dir, "famine.dta")) %>% 
  mutate(x1 = lgrain_pred*famine, x2 = lgrain_pred * (1-famine)) %>% 
  as.data.table()

N  <- nrow(df)

m0 <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + factor(year), data = df)
m1 <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + factor(year), 
         data = filter(df, year %in% 1953:1965))

# i), ii), iii)
se0 <- dfadjustSE(m0)
se1 <- dfadjustSE(m1)

# Bootstrap
bs <- function(df, N){
  df <- analysis(df)
  m  <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df)
  c  <- coef(m)[c("x1", "x2")]
  t  <- fixest::tstat(m)[c("x1", "x2")]
  return(c(c = c, t = t))
}

M <- 50000
bb <- bootstraps(df, times = M) 
out <- future_map(.x = bb$splits, .f = bs, N = N,                  
                  .options = furrr_options(seed = seed), 
                  .progress = TRUE)
draws <- bind_rows(out)

sd(draws$c.x1)
sd(draws$c.x2)


# problem 1b --------------------------------------------------------------
c_se0 <- dfadjustSE(m0, clustervar = as.factor(df$prov))
c_se1 <- dfadjustSE(m1, clustervar = 
                      as.factor(filter(df, year %in% 1953:1965)$prov))


bs_cl <- function(df, N){
  
  df <- df %>% as.data.frame() %>% unnest(cols = c(data))
  
  m <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df, 
             cluster = ~prov)
  c <- coef(m)[c("x1", "x2")]
  t <- fixest::tstat(m)[c("x1", "x2")]
  return(c(c, t))
}

M <- 50000
bb_cl <- bootstraps(df_cl, times = M) 
out_cl <- future_map(.x = bb_cl$splits, .f = bs_cl, N = N,                  
                  .options = furrr_options(seed = seed), 
                  .progress = TRUE)
draws_cl <- bind_rows(out_cl)

# clean up  ---------------------------------------------------------------


get_table <- function(se0, se1, draws){
  se0 <- se0$coefficients
  se1 <- se1$coefficients
  
  ci1 <- quantile(draws[,1], c(0.025, 0.975)) %>% round(digits = 5)
  ci2 <- quantile(draws[,2], c(0.025, 0.975)) %>% round(digits = 5)
  
  beta1 <- se0['x1', 'Estimate'] %>% round(digits = 5)
  beta2 <- se1['x2', 'Estimate'] %>% round(digits = 5)
  
  tibble(
    variable =     c("x1", "x2"),
    Estimate =     c(beta1, beta2), 
    Robust_se =    c(se0['x1', 'HC1 se'],   se1['x2', 'HC1 se']), 
    HC2 =          c(se0['x1', 'HC2 se'],   se1['x2', 'HC2 se']), 
    Effective_SE = c(se0['x1', 'Adj. se'],  se1['x2', 'Adj. se']), 
    BS = c(sd(draws[,1]), sd(draws[,2])), 
    ci = c(
      paste0("(", beta1 - ci1[2], "),(", beta1 - ci1[1], ")"),
      paste0("(", beta2 - ci2[2], "),(", beta2 - ci2[1], ")")
           )) 
}
get_table(se0, se1, draws) %>% write_csv(paste0(out, "t1.csv"))
get_table(c_se0, c_se1, draws_cl) %>% write_csv(paste0(out, "t2.csv"))


# problem 3 ---------------------------------------------------------------

df  <- read_dta(paste0(dir2, "ak91.dta")) %>% 
  filter(cohort == 2) %>% 
  mutate(Z = ifelse(age == floor(age), 1, 0))


m1 <- feols(data = df %>% filter(division == 9) , lwage ~ 1 | educ ~ Z)
m2 <- feols(data = df %>% filter(division == 2) , lwage ~ 1 | educ ~ Z)

get_wald <- function(mm){
  mm$coefficients['fit_educ'] + 
    c(-1.96*se(mm)['fit_educ'], + 1.96*se(mm)['fit_educ'])
}

get_ar <- function(division, df, min_b, max_b, step){
  df <- df %>% filter(division == !!division)

  a1 <- map_dfr(seq(min_b, max_b, step), function(bb){
    df$bX <- bb*df$educ
    df$YY <- df$lwage - df$bX
    m <- feols(YY~Z, data = df)
    tibble(b = bb, val = 1*(fixest::pvalue(m)['Z'] > 0.05))
    }
  )
  
  ci_ar_f <- a1 %>% filter(val == 1)
  c(min(ci_ar_f$b), max(ci_ar_f$b))
}

get_wald(m1)
get_wald(m2)
get_ar(9, df, -1,1,0.01)
get_ar(2, df, -100,100,1)

# First stage F is 9.4, so tF is 1.808 
tf <- 1.808 
tf_ci <- mm$coefficients['fit_educ'] + 
  c(-tf*1.96*se(mm)['fit_educ'], + tf*1.96*se(mm)['fit_educ'])

# Inifnite length for the oterh one, first stage F is tiny