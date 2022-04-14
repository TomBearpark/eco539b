pacman::p_load(tidyverse, haven, fixest, 
               dfadjust, xtable, data.table, rsample, 
               parallel, pbmcapply, sandwich, RDHonest, patchwork, 
               rdrobust)
theme_set(theme_bw())
seed <- 123
set.seed(seed)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps4/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps4/"
dir.create(out)


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
# q2 ----------------------------------------------------------------------
# -------------------------------------------------------------------------




# -------------------------------------------------------------------------
# q3 ----------------------------------------------------------------------
# -------------------------------------------------------------------------

df3 <- read_dta(file.path(dir, "njmin_clean.dta"))



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
x <- df4$povrate60
bw.c <- rdbwselect(y2,x, p = 3, kernel = "uniform")
bw.c$bws


# Point estimates ---------------------------------------------------------



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
