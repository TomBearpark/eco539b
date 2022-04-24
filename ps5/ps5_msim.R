pacman::p_load(tidyverse, xtable, patchwork)
theme_set(theme_bw())
seed <- 123 # seed <- 124
set.seed(seed)

dir <- "/Users/tombearpark/Documents/princeton/2nd_year/term2/eco539b/psets/ps5/"
out <- "/Users/tombearpark/Dropbox/Apps/Overleaf/eco539/b/figs/ps5/"
dir.create(out, showWarnings = FALSE)

# p3  -------------------------------------------------------------------

# Load data, store the size
df <- as.matrix(read_csv(paste0(dir, "ddc.csv"), col_names = paste0("t", 1:10)))
M <- 10
TT <- ncol(df)
N <- nrow(df)

# FP - I think this is not correct. It should read
S <- function(i, data, TT = 10){
 pi1 <- sum(data[i,c(1:9)] * data[i,c(2:10)])/(TT-2)
 pi2 <- sum(data[i,c(1:8)] * data[i,c(3:10)])/(TT-3)
 c(pi1, pi2)
}

get.pi.hat <- function(data, N = 100) {
  rowMeans(sapply(1:N, S, data = data, simplify = TRUE))
}

# Get empty matrices for storing values
U <- Y <- eps <- array(dim = c(N, TT, M))

# 1. Draw epsilons - all TT * N * M in one go
for (m in 1:M) eps[, 1:TT, m] <- matrix(rnorm(n = (TT * N)), ncol = TT)

# 2. Calculate the U as a function of rho
gen.U.m <- function(rho, eps, U0, N = 100, TT = 10) {
  U <- matrix(nrow = N, ncol = TT)
  # Draw first U obs from the unconditional distribution
  U[, 1] <- rho * U0 + eps[, 1]
  # The rest are from the conditionals
  for (ii in 2:TT) U[, ii] <- rho * U[, ii - 1] + eps[, ii]
  U
}

# 3. Calculate Ys for each simulation
# Y <- 1*(U>0)

# 4. Calculate the \Pi \Hat
gen.pi.tilde <- function(rho, eps, fix.U0 = TRUE) {
  
  U <- Y <- array(dim = c(N, TT, M))

  # Generate the U values
  # 1.1 Draw the U0 values from the unconditional distributions?
  if (fix.U0) {
    U0    <- matrix(0, ncol = M, nrow = N)
  } else {
    sigma <- sqrt(1 / (1 - rho^2))
    U0    <- matrix(rnorm(M * N, mean = 0, sd = sigma), ncol = M)
  }

  for (m in 1:M) U[, , m] <- gen.U.m(rho = rho, eps = eps[, , m], U0 = U0[, m])

  # Convert them into Ys
  Y <- 1 * (U > 0)

  # Calculate the pi.hat for each m in 1:M
  pi.m <- sapply(1:M, function(m) get.pi.hat(Y[, , m]))

  # Calculate the mean across the M draws
  pi.tilde <- rowMeans(pi.m)
  pi.tilde
}

gen.pi.tilde(rho = 0.1, eps = eps)

# 5. Calculate \hat W
Si      <- sapply(1:N, S, data = df, simplify = TRUE)
sd_diag <- apply(Si, 1, sd)
sd_cov  <- cov(Si[1, ], Si[2, ])
Omega   <- matrix(c(sd_diag[1]^2, sd_cov, sd_cov, sd_diag[2]^2), nrow = 2)
W       <- solve(Omega)
print(xtable(W, digits = 6), include.rownames = FALSE)

objective <- function(rho, W, pi_hat, eps) {
  p.tilde <- gen.pi.tilde(rho = rho, eps = eps)
  pp      <- matrix((pi_hat - p.tilde))
  drop(t(pp) %*% W %*% pp)
}

# Get the pi.hat for the original sample
pi_hat <- get.pi.hat(df)

# Plot the objective function ---------------------------------------------
rho_grid <- seq(0, 1, 0.005)
results <- tibble(
  rho = rho_grid,
  obj = map_dbl(rho_grid,
    objective,
    W = W, pi_hat = pi_hat, eps = eps
  )
)

pp1 <- ggplot(results) +
  geom_line(aes(x = rho, y = obj))
pp2 <- ggplot(filter(results, between(rho, 0.5, 0.75))) +
  geom_line(aes(x = rho, y = obj))
pp  <- pp1 + pp2
pp
ggsave(filename = paste0(out, "obj.png"), plot = pp, height = 3, width = 6)

# Optimize ----------------------------------------------------------------

pp <- optimize(
  f = objective, interval = c(0.4, 0.9),
  W = W, pi_hat = pi_hat, eps = eps,
  tol = 0.0000001
)

tibble(
  variable = c("rho_hat",  "min obj", "obj(0,1)", "obj(0,6)"),
  value = c(pp$minimum, pp$objective, 
            objective(0.1, W, pi_hat, eps), objective(0.6, W, pi_hat, eps)
  )
) %>%
  xtable(digits = 6) %>%
  print(include.rownames = FALSE)

# Avar --------------------------------------------------------------------

(1 + 1 / M)
