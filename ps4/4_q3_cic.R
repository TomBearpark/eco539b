pacman::p_load(tidyverse, haven, xtable, qte, rsample)
theme_set(theme_bw())

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps4/"

# load data ---------------------------------------------------------------
df3 <- read_dta(file.path(dir, "njmin_clean.dta"))
vars <- c("FTE", "inctime", "p_entree")

# funcs -------------------------------------------------------------------

# Pivot to a more convenient format
pivot_to_changes <- function(df, var){
  df %>% 
    select(storeid, after, state, .data[[var]]) %>% 
    pivot_wider(id_cols = c(storeid, state), 
                names_from = after, values_from = .data[[var]], 
                names_prefix = "period") %>% 
    mutate(!!var := period1 - period0) %>% 
    na.omit()
}

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
  print(paste0("F_00(10) is : ", F_00(10)))
  print(paste0("Counterfactual for F_00(10) is : ", invCDF(data$Y01, F_00(10))))
  print(paste0("ATT is ", att))
  
  bs     <- bootstraps(data$df, times = 10) 
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





