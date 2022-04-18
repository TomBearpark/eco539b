pacman::p_load(tidyverse, haven, fixest, xtable, 
               RDHonest, patchwork, rdrobust, qte, rsample)
theme_set(theme_bw())
seed <- 123
set.seed(seed)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps4/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps4/"
dir.create(out, showWarnings = FALSE)

# -------------------------------------------------------------------------
# q1 ----------------------------------------------------------------------
# -------------------------------------------------------------------------

N <- 100
TT <- 100
beta <- 4

# Create simulation data
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

# Generate variables needed for the mundlak / within estimations
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
  
# one way 

feols(y ~ x | i, data = df)
feols(y ~ x + xbar_n, data = df)
feols(y ~ x_min_xbar_n, data = df)

# two way

feols(y ~ x | i + t, data = df)
feols(y ~ x + xbar_n + xbar_t + 1, data = df)
feols(y ~ x_min_both, data = df)
feols(y ~ -1 + x_min_both_and_full, data = df)


# -------------------------------------------------------------------------
# q3 ----------------------------------------------------------------------
# -------------------------------------------------------------------------

# Load the data, specify outcome variables we are interested in 
df3 <- read_dta(file.path(dir, "njmin_clean.dta"))
vars <- c("FTE", "inctime", "p_entree")

# plot the ecdfs so we know what we are dealing with... 
df3 %>% 
  select(state, after, all_of(vars)) %>% 
  pivot_longer(cols = all_of(vars), names_to = "variable") %>% 
  mutate(state = ifelse(state == 1, "NJ", "PA"), 
         period = ifelse(after == 1, "After", "Before"), 
         period = factor(period, levels = c("Before", "After"))) %>%
  ggplot() + 
  stat_ecdf(aes(x = value, color = state)) +
  facet_wrap(variable ~ period, scales = "free", ncol = 2) -> e1
e1
ggsave(plot = e1, filename = paste0(out, "ecdf_raw.png"), height = 5, width = 6)

# Helper function: reshape the data for convenient estimations
pivot_to_changes <- function(df, var){
  df %>% 
    select(storeid, after, state, .data[[var]]) %>% 
    pivot_wider(id_cols = c(storeid, state), 
                names_from = after, values_from = .data[[var]], 
                names_prefix = "period") %>% 
    mutate(!!var := period1 - period0) %>% 
    na.omit()
}

# Exact replication of table from the notes... 
run_did <- function(df3, var){
  
  dd <- df3 %>% pivot_to_changes(var = var)

  ff <- feols(as.formula(paste0(var, "~ state")), data = dd, vcov = "hetero") 
  tibble(var = var, coef = coef(ff)['state'], se = se(ff)['state'], N = ff$nobs)
}

map_dfr(vars, run_did, df3 = df3) %>% 
  xtable(digits = 3) %>% 
  print(include.rownames = FALSE)
# Check intuition - we can get the same thing using a FE regression
feols(FTE ~ after + state + i(state, after, ref = '0')| storeid, data = df3)
feols(inctime ~ after + state + i(state, after, ref = '0')| storeid, data = df3)
feols(p_entree ~ after + state + i(state, after, ref = '0')| storeid, data = df3)

# CiC implemented in separate file

format_ys <- function(dd){
  Y00 <- dd %>% filter(state == 0) %>% pull(period0)
  Y01 <- dd %>% filter(state == 0) %>% pull(period1)
  Y10 <- dd %>% filter(state == 1) %>% pull(period0)
  Y11 <- dd %>% filter(state == 1) %>% pull(period1)
  
  list(Y00 = Y00, Y01 = Y01, Y10 = Y10, Y11 = Y11)
}
# Remove observations where overlap fails, return convenient vectors
subset_overlap <- function(df3, var){
  
  dd <- df3 %>% pivot_to_changes(var = var)
  
  support_p0_PA <- range(dd %>% filter(state == 0) %>% pull(period0))
  dd <- filter(dd, between(period0, support_p0_PA[1], support_p0_PA[2]))
  
  N <- length(dd$storeid)
  N1 <- length(dd$storeid)
  
  print(var)
  print(paste0("--- N was ", N, ", new N is ", N1))
  print(paste0("---- we dropped ", N-N1, " observations ----"))
  
  ys <- format_ys(dd)
  out <- c(ys,list(df = dd))
}

invCDF <- function(Y, q){
  fun_ecdf <- ecdf(Y)
  idx      <- fun_ecdf(Y) >= q
  min(Y[idx], na.rm = TRUE)
}

ATT <- function(data){
  
  F_00    <- ecdf(data$Y00)
  F00_y10 <- F_00(data$Y10)
  imputed <- map_dbl(F00_y10, function(q) invCDF(q, Y = data$Y01))
  att     <- mean(data$Y11) - mean(imputed)
  
  att
}

bs_func <- function(data){
  draw <- as.data.frame(data)
  ys <- format_ys(draw)
  ATT(ys)
}

#  run the code -----------------------------------------------------------

dfs <- map(vars, subset_overlap, df3 = df3)

for(ii in 1:3){
  data    <- dfs[[ii]]
  F_00    <- ecdf(data$Y00)
  att <- ATT(data)
  print(paste0("F_00(10) is : ", F_00(10)))
  print(paste0("Counterfactual for F_00(10) is : ", invCDF(data$Y01, F_00(10))))
  print(paste0("ATT is ", att))
  
  bs     <- bootstraps(data$df, times = 500) 
  out_bb <- map_dbl(bs$splits, bs_func)
  print(paste0("SD is ", sd(out_bb)))
}


# Version with canned package ---------------------------------------------
run_canned <- function(var, df){
  
  data <- df  %>% 
    select(-.data[[var]]) %>% 
    pivot_longer(cols = c(period0, period1), 
                 values_to = var, names_to = "after", 
                 names_prefix = "period", 
                 names_transform = list(after=as.numeric)) %>% 
    as.data.frame()
  
  ff <- as.formula(paste0(var, "~ state"))
  
  cc <- CiC(ff, t = 1, tmin1 = 0, tname = "after", data = data, panel = TRUE, 
            idname = "storeid")
  tibble(outcome = var, ate = cc$ate, se = cc$ate.se)
}

canned_results <- map_dfr(1:3, 
                          function(ii) 
                            run_canned(var = vars[ii], df = dfs[[ii]]$df))

xtable(canned_results) %>% 
  print(include.rownames=FALSE)


# -------------------------------------------------------------------------
# q4 ----------------------------------------------------------------------
# -------------------------------------------------------------------------
load(file.path(dir, "headst.rda"))

df4 <- tibble(headst) %>% 
  arrange(-povrate60) %>% 
  mutate(rank = row_number(),
         G = ifelse(rank >=300,"1","0"))

cutoff <- df4$povrate60[which(df4$rank == 300)]

# Plot density - McCrary style
df4 %>% pivot_longer(cols = c(mortHS, mortInj)) %>% 
  ggplot() + 
  geom_histogram(aes(x = povrate60), alpha = 0.5, bins = 100) + 
  geom_vline(xintercept = 0) + 
  facet_wrap(~name, scales = "free") 

# Plot the raw points, colored by treatment status
df4 %>% pivot_longer(cols = c(mortHS, mortInj)) %>% 
  ggplot() + 
  geom_point(aes(x = povrate60, y = value, color = G), alpha = 0.5) + 
  facet_wrap(~name, scales = "free") -> 
  p
p
ggsave(plot = p, filename = paste0(out, "raw_scatters.png"), height = 3, width = 6)
p + xlim(-10, 10)

# Any evidence of jump in means?
feols(c(mortHS, mortInj) ~ povrate60 + G,data = df4)

# using RDhonest ----------------------------------------------------------

d1 <- select(df4, "mortHS", povrate60) %>% 
  na.omit() %>% 
  RDData(cutoff = 0)

d2 <- select(df4, "mortInj", povrate60) %>% 
  na.omit() %>% 
  RDData(cutoff = 0)

p1 <- plot_RDscatter(d1, avg = 25, window = 60, ylab = "mortHS")
p2 <- plot_RDscatter(d2, avg = 25, window = 60, ylab = "mortInj")

pp <- (p1 + p2)
ggsave(plot = pp, filename = paste0(out, "binned_scatters.png"), 
       height = 3, width = 6)


M1 <- NPR_MROT.fit(d1)
M2 <- NPR_MROT.fit(d2)

get_bw <- function(M, var, df4){
  
  ff <- as.formula(paste0(var, "~ povrate60"))
  
  rdd <- RDHonest(ff, data=df4, M = M, 
           opt.criterion = "MSE", cutoff = 0, kern = "uniform") 
  tibble(M = !!M, bw= rdd$hm, estimate = rdd$estimate)
}

bw_plot <- function(var, df4, M_star){
  rr1 <- map_dfr(c(M_star, seq(0.01,1,by = 0.01)), get_bw, var = var, df4 = df4)
  bw_opt1 <- rr1 %>% filter(M == !!M_star) %>% pull(bw)
  rr1 %>% 
    ggplot() + 
    geom_point(aes(x = M, y = bw)) + 
    geom_vline(xintercept = M_star, color = "red", alpha = 0.6) + 
    geom_hline(yintercept = bw_opt1, color = "red", alpha = 0.6) + 
    ggtitle(var)
}

p_b1 <- bw_plot("mortHS", df4, M1)
p_b2 <- bw_plot("mortInj", df4, M2)

p_b <- (p_b1 + p_b2)

ggsave(plot = p_b, file.path(out, paste0("opt_bw.png")), height = 4, width = 7)


# Robustness - use Cattaneo package
y1 <- df4$mortHS
y2 <- df4$mortInj
x  <- df4$povrate60
bw.c <- rdbwselect(y2,x, p = 3, kernel = "uniform")
bw.c$bws


# Point estimates ---------------------------------------------------------

run_rdd_honest <- function(M, df4 ,var){
  ff <- as.formula(paste0(var, "~ povrate60"))
  r1 <- RDHonest(ff, data = df4, M = M, kern = "uniform", 
           opt.criterion = "MSE")
  tibble(M = M, bw = r1$hp, 
         upper = r1$upper, lower = r1$lower, estimate = r1$estimate)
}

rr1 <- map_dfr(c(seq(0.05,0.5,by = 0.05)), run_rdd_honest, var ="mortHS" , df4 = df4)
rr1$M_rule <- ifelse(between(rr1$M, 0.28, 0.32), "Rule of thumb", "NA")

rr2 <- map_dfr(c(seq(0.1,0.9,by = 0.1)), run_rdd_honest, var ="mortInj" , df4 = df4)
rr2$M_rule <- ifelse(between(rr2$M, 0.8, 0.82), "Rule of thumb", "NA")


plot_cis <- function(rr1, var){
  rr1 %>% pivot_longer(cols = c(M, bw), names_to= "xvar") %>% 
    mutate(xvar = factor(xvar, levels = c("M", "bw"))) %>% 
    ggplot() + 
    geom_point(aes(x = value, y = estimate, color = M_rule)) + 
    geom_errorbar(aes(x = value, ymin = lower, ymax = upper, color = M_rule)) + 
    geom_hline(yintercept = 0, color = "red") + 
    facet_wrap(~xvar, scales = "free") + 
    ggtitle(var)
}

pf1 <- plot_cis(rr1, "mortHS")
pf2 <- plot_cis(rr2, "mortInj")
pf <- pf1 / pf2

ggsave(plot = pf, filename = file.path(out, paste0("ci_rdd.png")), 
       height = 6, width = 9)
