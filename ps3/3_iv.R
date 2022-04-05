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
beta0 <- coef(m0)[c("x1", "x2")]
beta1 <- coef(m1)[c("x1", "x2")]

# i), ii), iii)
se0 <- dfadjustSE(m0)
se1 <- dfadjustSE(m1)

# Bootstrap
bs <- function(df){
  df <- as.data.frame(df)
  m  <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df)
  c  <- coef(m)[c("x1", "x2")]
  se  <- fixest::se(m)[c("x1", "x2")]
  return(c(c = c, se = se))
}
mc_bs <- function(df, M, bs_func, beta){
  
  bb     <- bootstraps(df, times = M) 
  out_bb <- pbmclapply(bb$splits, bs_func, mc.cores = 5)
  draws  <- bind_rows(out_bb)
  # Create t-stat
  draws %>% 
    mutate(t.x1 = (c.x1 - beta[1])/se.x1, t.x2 = (c.x2 - beta[2])/se.x2)
}

draws0 <- mc_bs(df,     M, bs_func = bs, beta = beta0)
draws1 <- mc_bs(df_sub, M, bs_func = bs, beta = beta1)

gc()

# problem 1b --------------------------------------------------------------
se0_cl <- dfadjustSE(m0, clustervar = as.factor(df$prov))
se1_cl <- dfadjustSE(m1, clustervar = as.factor(df_sub$prov))

bs_cl <- function(df){
  
  df <- df %>% as.data.frame() %>% unnest(cols = c(data))
  
  m <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df, 
             cluster = ~prov)
  c  <- coef(m)[c("x1", "x2")]
  se  <- fixest::se(m)[c("x1", "x2")]
  return(c(c = c, se = se))
}

df_cl     <- df %>% group_nest(prov)
df_sub_cl <- df_sub %>% group_nest(prov)

draws0_cl <- mc_bs(df_cl,     M, bs_func = bs_cl, beta = beta0)
draws1_cl <- mc_bs(df_sub_cl, M, bs_func = bs_cl, beta = beta1)


# clean up  ---------------------------------------------------------------

get_table <- function(se, draws, data_string, cluster_string, N){
  se <- se$coefficients
  se <- round(se, 3)
  ci1 <- (quantile(draws$t.x1, c(0.025, 0.975)) / sqrt(N)) 
  ci2 <- (quantile(draws$t.x2, c(0.025, 0.975)) / sqrt(N)) 
  
  beta1 <- se['x1', 'Estimate'] 
  beta2 <- se['x2', 'Estimate'] 
  
  rbind(
    variable =     c("x1", "x2"),
    Estimate =     round(c(beta1, beta2),3),
    Robust_se =    c(se['x1', 'HC1 se'],   se['x2', 'HC1 se']),
    HC2 =          c(se['x1', 'HC2 se'],   se['x2', 'HC2 se']), 
    Effective_SE = c(se['x1', 'Adj. se'],  se['x2', 'Adj. se']), 
    BS = round(c(sd(draws$c.x1), sd(draws$c.x2)),3), 
    ci = c(
      paste0("(", round(beta1 - ci1[2],3), ",", round(beta1 - ci1[1],3), ")"),
      paste0("(", round(beta2 - ci2[2],3), ",", round(beta2 - ci2[1],3), ")")
           )
    ) 
}

rbind(
  get_table(se0, draws0, "All", "None", N0) ,
  get_table(se0_cl, draws0_cl, "All", "Province", N0)[3:7, ] 
) %>% 
  cbind(
    rbind(
      get_table(se1, draws1, "1953:1965", "None", N1), 
      get_table(se1_cl, draws1_cl, "1953:1965", "Province",N1)[3:7, ] 
    )
  ) %>% 
  xtable(digits = 3) %>% 
  print()



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
    c(-qnorm(0.975)*se(mm)['fit_educ'], + qnorm(0.975)*se(mm)['fit_educ'])
}

get_ar <- function(division, df, min_b, max_b, step){
  df <- df %>% filter(division == !!division)

  a1 <- map_dfr(seq(min_b, max_b, step), function(bb){
    df$bX <- bb*df$educ
    df$YY <- df$lwage - df$bX
    m <- feols(YY~Z, data = df, vcov = "hetero")
    tibble(b = bb, val = 1*(fixest::pvalue(m)['Z'] > 0.05))
    }
  )
  
  ci_ar_f <- a1 %>% filter(val == 1)
  c(min(ci_ar_f$b), max(ci_ar_f$b))
}

format_ci <- function(ci) paste0("(", ci[1], ", ", ci[2], ")")

w1 <- get_wald(m3.1) %>% round(3)
w2 <- get_wald(m3.2) %>% round(3)
ar1 <- get_ar(9, df3, -1,1,0.001)
ar2 <- get_ar(2, df3, -100,100,10)

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

load(file = file.path(dir, "run.RData"))
