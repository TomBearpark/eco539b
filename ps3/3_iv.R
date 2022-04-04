pacman::p_load(tidyverse, haven, fixest, 
               dfadjust, xtable, data.table, rsample, parallel, pbmcapply, sandwich)
theme_set(theme_bw())
seed <- 123
set.seed(seed)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps3/"
dir2 <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps2/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps3/"


# problem 1a ---------------------------------------------------------------
M <- 50000

df  <- read_dta(paste0(dir, "famine.dta")) %>% 
  mutate(x1 = lgrain_pred*famine, x2 = lgrain_pred * (1-famine), 
         year = factor(year))

df_sub <- filter(df, year %in% 1953:1965) 

N0  <- nrow(df)
N1  <- nrow(df_sub)

m0 <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + year, data = df)
m1 <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + year, data = df_sub)

# i), ii), iii)
se0 <- dfadjustSE(m0)
se1 <- dfadjustSE(m1)

# Bootstrap
bs <- function(df){
  df <- as.data.frame(df)
  m  <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df)
  c  <- coef(m)[c("x1", "x2")]
  t  <- fixest::tstat(m)[c("x1", "x2")]
  return(c(c = c, t = t))
}
mc_bs <- function(df, M, bs_func){
  
  bb     <- bootstraps(df, times = M) 
  out_bb <- pbmclapply(bb$splits, bs_func, mc.cores = 5)
  draws  <- bind_rows(out_bb)
  draws
}

draws0 <- mc_bs(df,     M, bs_func = bs)
draws1 <- mc_bs(df_sub, M, bs_func = bs)

gc()

# problem 1b --------------------------------------------------------------
se0_cl <- dfadjustSE(m0, clustervar = as.factor(df$prov))
se1_cl <- dfadjustSE(m1, clustervar = as.factor(df_sub$prov))

bs_cl <- function(df){
  
  df <- df %>% as.data.frame() %>% unnest(cols = c(data))
  
  m <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df, 
             cluster = ~prov)
  c <- coef(m)[c("x1", "x2")]
  t <- fixest::tstat(m)[c("x1", "x2")]
  return(c(c = c, t = t))
}

df_cl     <- df %>% group_nest(prov)
df_sub_cl <- df_sub %>% group_nest(prov)

draws0_cl <- mc_bs(df_cl,     M, bs_func = bs_cl)
draws1_cl <- mc_bs(df_sub_cl, M, bs_func = bs_cl)


# clean up  ---------------------------------------------------------------

get_table <- function(se, draws, data_string, cluster_string, N){
  se <- se$coefficients
  se <- round(se, 3)
  ci1 <- (quantile(draws$t.x1, c(0.025, 0.975)) / sqrt(N)) %>% round(digits = 3)
  ci2 <- (quantile(draws$t.x2, c(0.025, 0.975)) / sqrt(N)) %>% round(digits = 3)
  
  beta1 <- se['x1', 'Estimate'] %>% round(digits = 3)
  beta2 <- se['x2', 'Estimate'] %>% round(digits = 3)
  
  tibble(
    variable =     c("x1", "x2"),
    Estimate =     c(beta1, beta2), 
    Robust_se =    c(se['x1', 'HC1 se'],   se['x2', 'HC1 se']), 
    HC2 =          c(se['x1', 'HC2 se'],   se['x2', 'HC2 se']), 
    Effective_SE = c(se['x1', 'Adj. se'],  se['x2', 'Adj. se']), 
    BS = c(sd(draws$c.x1), sd(draws$c.x2)), 
    ci = c(
      paste0("(", beta1 - ci1[2], ",", beta1 - ci1[1], ")"),
      paste0("(", beta2 - ci2[2], ",", beta2 - ci2[1], ")")
           )) %>%
    mutate(Data = data_string, Cluster = cluster_string) %>% 
    relocate(Cluster, Data)
}

bind_rows(
  get_table(se0, draws0, "All", "None", N0),
  get_table(se1, draws1, "1953:1965", "None", N1),
  get_table(se0_cl, draws0_cl, "All", "Province", N0),
  get_table(se1_cl, draws1_cl, "1953:1965", "Province",N1) 
  ) %>% 
  xtable(digits = 3) %>% 
  print(include.rownames=FALSE)


# problem 3 ---------------------------------------------------------------

df3  <- read_dta(paste0(dir2, "ak91.dta")) %>% 
  filter(cohort == 2) %>% 
  mutate(Z = ifelse(age == floor(age), 1, 0))

m3.1.ols <- feols(data = df3 %>% filter(division == 9) , lwage ~ educ)
m3.1 <- feols(data = df3 %>% filter(division == 9) , lwage ~ 1 | educ ~ Z)
m3.2.ols <- feols(data = df3 %>% filter(division == 2) , lwage ~ educ)
m3.2 <- feols(data = df3 %>% filter(division == 2) , lwage ~ 1 | educ ~ Z)

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

format_ci <- function(ci) paste0("(", ci[1], ", ", ci[2], ")")

w1 <- get_wald(m3.1) %>% round(3)
w2 <- get_wald(m3.2) %>% round(3)
ar1 <- get_ar(9, df3, -1,1,0.01)
ar2 <- get_ar(2, df3, -100,100,1)

# First stage F is 9.4, so tF is 1.808 
tf <- 1.808 
tf_ci <- (m3.1$coefficients['fit_educ'] + 
  c(-tf*1.96*se(m3.1)['fit_educ'], + tf*1.96*se(m3.1)['fit_educ'])) %>% round(3)

tibble(
  Division = c("Pacific", "Mid-Atlantic"),
  `OLS Estimate` = c(round(coef(m3.1.ols)['educ'],3), c(round(coef(m3.2.ols)['educ'],3))),
  `2sls Estimate` = c(round(coef(m3.1)['fit_educ'],3), round(coef(m3.2)['fit_educ'],3)),
  `First Stage F` = c(round(fitstat(m3.1, "ivf1")$`ivf1::educ`$stat,3),
                      round(fitstat(m3.2, "ivf1")$`ivf1::educ`$stat,3)),
  Wald = c(format_ci(w1), format_ci(w2)), 
  AR   = c(format_ci(ar1), format_ci(ar2)), 
  tF = c(format_ci(tf_ci), Inf)
  ) %>% 
  t() %>% 
  xtable() 
  

# Save outputs, since its a long peice of code 
save.image(file = file.path(dir, "run.RData"))
