pacman::p_load(tidyverse, haven, fixest, sandwich, dfadjust, xtable)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps2/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps2/"
theme_set(theme_bw())
set.seed(123)

df  <- read_dta(paste0(dir, "ak91.dta")) %>% 
  filter(cohort == 2) %>% 
  select(lwage, SOB, educ)

# a) ----------------------------------------------------------------------

df %>% 
  ggplot(aes(x = educ, y = lwage)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth()

m0    <- lm(lwage~educ, data = df)
m0.1
m0.fe   <- feols(lwage~educ|SOB, data = df)
m0.fe_i <- feols(lwage~educ + i(SOB,educ)|SOB, data = df)
coefplot(m0.fe_i)

prop <- feols(educ~-1 +i(SOB), data = df)
coefplot(prop)



m0.poly <- feols(lwage~poly(educ, 6, raw = TRUE), data = df, vcov = "hetero")

d1 <- dfadjustSE(m0)
d1$coefficients %>% xtable(digits = 5) %>%   print(include.rownames=FALSE)
# Calculate leverage
H <- hatvalues(m0)

ggplot(tibble(leverage = H)) + 
  geom_density(aes(x = leverage)) + 
  geom_vline(aes(xintercept = 1/length(H)), color = "red")

ggsave(paste0(out, "2a_lowess.png"), height = 3, width = 5)

df %>% group_by(SOB) %>% tally() %>% arrange( n) 
df %>% group_by(SOB) %>% tally() %>% arrange(-n) 

df$H <- H
df$resid <- m0$residuals

ggplot(df) + geom_point(aes(x = H, y = resid))

bootstrap <- function(df){
  draw <- df %>% 
    group_nest(SOB) %>% 
      slice_sample(prop = 1, replace = TRUE) %>% 
    unnest(cols = c(data))
  coef(feols(lwage~educ, data = draw))['educ']
}
bs <- c()
for(ii in 1:500) {
  print(ii)
  bs <- c(bs, bootstrap(df))
}
sd(bs)

loo <- map_dfr(unique(df$SOB), function(ii){
    ss <- df %>% filter(SOB != ii)
    tibble(left_out_SOB = ii, coef = coef(feols(lwage~educ, data = ss))['educ'])
  }
  )
loo %>% ggplot() + geom_density(aes(x = coef))

df$SOB <- as.factor(df$SOB)
d2 <- dfadjustSE(m0, clustervar = df$SOB)


tibble(
  Estimator = c("Homoskedastic", 
                "HC1 SE", 
                "HC2 SE", 
                "Clustered HC1", 
                "Clustered HC2", 
                "Clustered HC2 Kolesar Adjustment", 
                "Cluster Bootstrap"), 
  Value     = c(m1$se["educ"], 
                d1$coefficients["educ","HC1 se"], 
                d1$coefficients["educ","HC2 se"],
                d2$coefficients["educ", "HC1 se"], 
                d2$coefficients["educ", "HC2 se"], 
                d2$coefficients["educ", "Adj. se"], 
                sd(bs))
  ) %>% 
  arrange(Value) %>% 
  xtable(digits = 5) %>% 
  print(include.rownames=FALSE)

# b) ----------------------------------------------------------------------
# Despite the large data size, the effective number of observations is small
# since we only have one treated clusteres

df$NJ  <- as.factor(ifelse(df$SOB == "34", "1", "0"))
df$SOB <- as.factor(df$SOB)
b1     <- lm(lwage ~ NJ, data = df)

dfadjustSE(b1)$coefficients %>% xtable(digits = 5) %>% print(include.rownames=FALSE)

dfadjustSE(b1, clustervar = df$SOB, IK = FALSE)

hatvalues(b1) %>% density %>% plot()
