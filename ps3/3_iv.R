pacman::p_load(tidyverse, haven, fixest, sandwich, dfadjust, xtable, data.table)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps3/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps3/"
theme_set(theme_bw())
set.seed(123)


# problem 1a ---------------------------------------------------------------
df  <- read_dta(paste0(dir, "famine.dta")) %>% 
  mutate(x1 = lgrain_pred*famine, x2 = lgrain_pred * (1-famine)) %>% 
  as.data.table()

N  <- nrow(df)

m0 <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + factor(year), data = df)
m1 <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + factor(year), 
         data = filter(df, year %in% 1953:1965))

# i,ii,iii
se0 <- dfadjustSE(m0)
se1 <- dfadjustSE(m1)

# Bootstrap
bs <- function(i, df, N){
  df <- df[sample(.N,N, replace = TRUE),]
  m <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df)
  c <- coef(m)[c("x1", "x2")]
  t <- fixest::tstat(m)[c("x1", "x2")]
  return(c(c, t))
}
M <- 50000
draws <- matrix(nrow = M, ncol = 4)
for (ii in 1:M) draws[ii,] <- bs(ii, df, N)

sd(dd$x1)
sd(dd$x2)

# problem 1b --------------------------------------------------------------
c_se0 <- dfadjustSE(m0, clustervar = as.factor(df$prov))
c_se1 <- dfadjustSE(m1, clustervar = as.factor(df$prov))

bs_cl <- function(i, df, N){
  df <- df[sample(.N,N, replace = TRUE),]
  m <- feols(ldeaths ~ x1 + x2 + ltotpop + lurbpop | year, data = df)
  c <- coef(m)[c("x1", "x2")]
  t <- fixest::tstat(m)[c("x1", "x2")]
  return(c(c, t))
}

# Clean up all results into a single table 
se0 <- se0$coefficients
se1 <- se1$coefficients
tibble(
  variable =     c("x1", "x2"),
  Estimate =     c(se0['x1', 'Estimate'], se1['x2', 'Estimate']), 
  Robust_se =    c(se0['x1', 'HC1 se'],   se1['x2', 'HC1 se']), 
  HC2 =          c(se0['x1', 'HC2 se'],   se1['x2', 'HC2 se']), 
  Effective_SE = c(se0['x1', 'Adj. se'],  se1['x2', 'Adj. se']), 
  BS = c(sd(draws[,1]), sd(draws[,2])
) 
  