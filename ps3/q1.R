setwd("D:/Dropbox/Università/PhD/II Year/Spring/ECO519 - Non-linear Econometrics/psets/ps3")

# Load stuff
library(haven)
library(sandwich)
library(dfadjust)
library(fastDummies)
library(xtable)
library(rsample)
library(tidyverse)

theme_set(theme_bw())

data <- haven::read_dta("famine.dta")

# Generate interactions
data$x1 <- data$lgrain_pred*data$famine 
data$x2 <- data$lgrain_pred*(1-data$famine) 

# Run regression
reg.out <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + as.factor(year) - 1, 
              data = data)
beta.hat <- reg.out$coefficients[1:2]

iters <- 50000  # bootstrap draws

###############################################################################
## Exercise 1a - Various standard errors allowing for heteroskedasticity
###############################################################################

# 1) Robust EHW variance covariance
hc1.vcov <- sandwich::vcovHC(reg.out, type = 'HC1')

# 2) Robust HC2 variance covariance
hc2.vcov <- sandwich::vcovHC(reg.out, type = 'HC2')

# 3) Satterthwaite (1946) adjustment for degrees of freedom 
sw.ses <- dfadjustSE(reg.out)$coefficients[1:2,4, drop = F]

# 4-5) non-parametric and percentile bootstrap
out.npboot <- matrix(NA, nrow=iters, 2)
out.pcboot <- matrix(NA, nrow=iters, 2)

# Bootstrap units
set.seed(8894)
bs <- rsample::bootstraps(data, times=iters)  

for (i in seq_len(iters)) {
  
  df.boot <- tibble::as_tibble(bs$splits[[i]]) 
  
  reg.boot <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + as.factor(year) - 1, 
                 data = df.boot)
  vcov <- sandwich::vcovHC(reg.boot, type = 'HC1')[1:2,1:2]
  
  out.npboot[i, ] <- reg.boot$coefficients[1:2]
  out.pcboot[i, ] <- (reg.boot$coefficients[1:2] - beta.hat)/sqrt(diag(vcov))
}

# compute standard deviation across draws (non-parametric bootstrap)
npboot.ses <- apply(out.npboot, 2, sd, na.rm=T)

# compute quantiles of robust t-statistic (percentile bootstrap)
pcboot.qtles <- apply(out.pcboot, 2, 
                  function(x) quantile(x, probs = c(0.025,0.975), na.rm=T))

#############################################################################
## Construct all confidence intervals
#############################################################################
# first variable
# HC1
CI1.hc1 <- c(beta.hat[1] - qnorm(0.975)*sqrt(hc1.vcov[1,1]),
             beta.hat[1] - qnorm(0.025)*sqrt(hc1.vcov[1,1]))

# HC2
CI1.hc2 <- c(beta.hat[1] - qnorm(0.975)*sqrt(hc2.vcov[1,1]),
             beta.hat[1] - qnorm(0.025)*sqrt(hc2.vcov[1,1]))

# Satterthwaite adjustment
CI1.sw <- c(beta.hat[1] - qnorm(0.975)*sw.ses[1],
             beta.hat[1] - qnorm(0.025)*sw.ses[1])

# Non-parametric bootstrap
CI1.np <- c(beta.hat[1] - qnorm(0.975)*npboot.ses[1],
            beta.hat[1] - qnorm(0.025)*npboot.ses[1])

# Percentile bootstrap on robust t-statistic
CI1.pct <- c(beta.hat[1] - pcboot.qtles[2,1]*sqrt(hc1.vcov[1,1]),
              beta.hat[1] - pcboot.qtles[1,1]*sqrt(hc1.vcov[1,1]))

# second variable
# HC1
CI2.hc1 <- c(beta.hat[2] - qnorm(0.975)*sqrt(hc1.vcov[2,2]),
             beta.hat[2] - qnorm(0.025)*sqrt(hc1.vcov[2,2]))

# HC2
CI2.hc2 <- c(beta.hat[2] - qnorm(0.975)*sqrt(hc2.vcov[2,2]),
             beta.hat[2] - qnorm(0.025)*sqrt(hc2.vcov[2,2]))

# Satterthwaite adjustment
CI2.sw <- c(beta.hat[2] - qnorm(0.975)*sw.ses[2],
            beta.hat[2] - qnorm(0.025)*sw.ses[2])

# Non-parametric bootstrap
CI2.np <- c(beta.hat[2] - qnorm(0.975)*npboot.ses[2],
            beta.hat[2] - qnorm(0.025)*npboot.ses[2])

# Percentile bootstrap on robust t-statistic
CI2.pct <- c(beta.hat[2] - pcboot.qtles[2,1]*sqrt(hc1.vcov[2,2]),
              beta.hat[2] - pcboot.qtles[1,1]*sqrt(hc1.vcov[2,2]))

methods <- c("Coefficient","HC1", "HC2", "Satterthwaite", "non-parametric bootstrap",
             "percentile bootstrap")
cols <- c("lb", "ub","lb", "ub")

tab1a <- rbind(CI1.hc1,CI1.hc2,CI1.sw,CI1.np,CI1.pct)
tab1b <- rbind(CI2.hc1,CI2.hc2,CI2.sw,CI2.np,CI2.pct)

tab1 <- rbind(c(rep(beta.hat[1],2),rep(beta.hat[2],2)),
              cbind(tab1a, tab1b))
colnames(tab1) <- cols
rownames(tab1) <- methods

###############################################################################
## Exercise 1b - Various standard errors allowing for province clusters
###############################################################################

# 1) Robust LR variance covariance
cr1.vcov <- sandwich::vcovCL(reg.out, cluster = data$prov, type = 'HC1')

# 2) Robust CR2 (Bell and McAffrey) variance covariance
cr2.ses <- dfadjustSE(reg.out, clustervar = as.factor(data$prov), 
                       IK = FALSE)$coefficients[1:2, 3, drop = F]

# 3) Satterthwaite (1946) adjustment for degrees of freedom 
cr2sw.ses <- dfadjustSE(reg.out, clustervar = as.factor(data$prov), 
                       IK = FALSE)$coefficients[1:2, 4, drop = F]

# 4-5) non-parametric and percentile bootstrap
out.npboot.cl <- matrix(NA, nrow=iters, 2)
out.pcboot.cl <- matrix(NA, nrow=iters, 2)

# Extract id for each province
ids <- data %>% nest('ID' = -prov)

# Bootstrap provinces (not observations!!) 
set.seed(8894)
bs.cl <- rsample::bootstraps(ids, times=iters) 

for (i in seq_len(iters)) {
  
  df.boot <- as_tibble(bs.cl$splits[[i]]) %>% unnest(cols = c(ID))

  reg.boot <- lm(ldeaths ~ x1 + x2 + ltotpop + lurbpop + as.factor(year) - 1, 
                 data = df.boot)
  
  vcov <- sandwich::vcovCL(reg.boot, cluster = df.boot$prov, 
                           type = 'HC1')[1:2,1:2]
  
  out.npboot.cl[i, ] <- reg.boot$coefficients[1:2]
  out.pcboot.cl[i, ] <- (reg.boot$coefficients[1:2] - beta.hat)/sqrt(diag(vcov))
}

# compute standard deviation across draws (non-parametric bootstrap)
npbootcl.ses <- apply(out.npboot.cl, 2, sd, na.rm=T)

# compute quantiles of robust t-statistic (percentile bootstrap)
pcbootcl.qtles <- apply(out.pcboot.cl, 2, 
                      function(x) quantile(x, probs = c(0.025,0.975), na.rm=T))

#############################################################################
## Construct all confidence intervals (cluster robust)
#############################################################################
# first variable
# CR1
CI1.cr1 <- c(beta.hat[1] - qnorm(0.975)*sqrt(cr1.vcov[1,1]),
             beta.hat[1] - qnorm(0.025)*sqrt(cr1.vcov[1,1]))

# CR2 (Bell and McCaffrey)
CI1.cr2 <- c(beta.hat[1] - qnorm(0.975)*cr2.ses[1],
             beta.hat[1] - qnorm(0.025)*cr2.ses[1])

# Satterthwaite adjustment
CI1.cr2sw <- c(beta.hat[1] - qnorm(0.975)*cr2sw.ses[1],
               beta.hat[1] - qnorm(0.025)*cr2sw.ses[1])

# Non-parametric bootstrap
CI1.npcl <- c(beta.hat[1] - qnorm(0.975)*npbootcl.ses[1],
              beta.hat[1] - qnorm(0.025)*npbootcl.ses[1])

# Percentile bootstrap on robust t-statistic
CI1.pctcl <- c(beta.hat[1] - pcbootcl.qtles[2,1]*sqrt(cr1.vcov[1,1]),
               beta.hat[1] - pcbootcl.qtles[1,1]*sqrt(cr1.vcov[1,1]))

# second variable
# CR1
CI2.cr1 <- c(beta.hat[2] - qnorm(0.975)*sqrt(cr1.vcov[2,2]),
             beta.hat[2] - qnorm(0.025)*sqrt(cr1.vcov[2,2]))

# CR2
CI2.cr2 <- c(beta.hat[2] - qnorm(0.975)*cr2.ses[2],
             beta.hat[2] - qnorm(0.025)*cr2.ses[2])

# Satterthwaite adjustment
CI2.cr2sw <- c(beta.hat[2] - qnorm(0.975)*cr2sw.ses[2],
               beta.hat[2] - qnorm(0.025)*cr2sw.ses[2])

# Non-parametric bootstrap
CI2.npcl <- c(beta.hat[2] - qnorm(0.975)*npbootcl.ses[2],
              beta.hat[2] - qnorm(0.025)*npbootcl.ses[2])

# Percentile bootstrap on robust t-statistic
CI2.pctcl <- c(beta.hat[2] - pcboot.qtles[2,1]*sqrt(cr1.vcov[2,2]),
               beta.hat[2] - pcboot.qtles[1,1]*sqrt(cr1.vcov[2,2]))

methods <- c("CR1", "CR2", "Satterthwaite", "non-parametric bootstrap",
             "percentile bootstrap")
cols <- c("lb", "ub","lb", "ub")

tab2a <- rbind(CI1.cr1,CI1.cr2,CI1.cr2sw,CI1.npcl,CI1.pctcl)
tab2b <- rbind(CI2.cr1,CI2.cr2,CI2.cr2sw,CI2.npcl,CI2.pctcl)

tab2 <- cbind(tab2a, tab2b)

tab <- rbind(tab1, tab2)

xtable(as.table(tab), digits = 3)


##############################################################################
## Compute leverages
##############################################################################

# Store design matrix, outcome, and residuals
X <- fastDummies::dummy_cols(as.factor(data$year),
                             remove_first_dummy = FALSE)[,-1]
X <- cbind(data$x1, data$x2, data$ltotpop, data$lurbpop, X)
X <- data.matrix(X)
y <- data$ldeaths
XX <- solve(t(X) %*% X)
beta.ls <- XX %*% t(X) %*% y
res <- y - X %*% beta.ls
N <- length(res)

# Compute leverage
leverage <- data.frame(lvg=stats::hatvalues(reg.out))
leverage$x <- c(1:nrow(leverage))
ggplot(leverage, aes(x=x, y=lvg)) + geom_point() +
  xlab("observation") + ylab("leverage") +
  geom_hline(yintercept = length(reg.out$coefficients)/nrow(leverage),
             color = "red")
ggsave('leverage.png', height = 4, width = 6, dpi = 1000)

# Partial leverage of X1
x1 <- X[,1]
xx1 <- X[,-1]
X1fwl <- x1 - xx1 %*% solve(t(xx1) %*% xx1) %*% t(xx1) %*% x1
lvgx1 <- diag(X1fwl %*% solve(t(X1fwl) %*% X1fwl) %*% t(X1fwl))
leverage <- data.frame(lvg=lvgx1)
leverage$x <- c(1:nrow(leverage))
ggplot(leverage, aes(x=x, y=lvg)) + geom_point() +
  xlab("observation") + ylab("partial leverage") +
  geom_hline(yintercept = length(reg.out$coefficients)/nrow(leverage),
             color = "red")
ggsave('leverage_x1.png', height = 4, width = 6, dpi = 1000)

# Partial leverage of X2
x1 <- X[,2]
xx1 <- X[,-2]
X1fwl <- x1 - xx1 %*% solve(t(xx1) %*% xx1) %*% t(xx1) %*% x1
lvgx1 <- diag(X1fwl %*% solve(t(X1fwl) %*% X1fwl) %*% t(X1fwl))
leverage <- data.frame(lvg=lvgx1)
leverage$x <- c(1:nrow(leverage))
ggplot(leverage, aes(x=x, y=lvg)) + geom_point() +
  xlab("observation") + ylab("partial leverage") +
  geom_hline(yintercept = length(reg.out$coefficients)/nrow(leverage),
             color = "red")
ggsave('leverage_x2.png', height = 4, width = 6, dpi = 1000)


# Compute within group leverage
sobs <- unique(data$prov)
G <- length(sobs)
storelev <- matrix(NA, G, 2) 

i <- 1
for (sob in sobs) {
  datacl <- subset(data, data$prov == sob)
  
  Xg <- fastDummies::dummy_cols(as.factor(datacl$year),
                               remove_first_dummy = FALSE)[,-1]
  
  if (i==1) Xg<-cbind(0,Xg)
  
  Xg <- cbind(datacl$x1, datacl$x2, datacl$ltotpop, datacl$lurbpop, Xg)
  Xg <- data.matrix(Xg)
  XXg <- t(Xg)%*%Xg
  
  storelev[i, 2] <- sum(diag(XXg %*% XX)) # compute nuclear norm
  storelev[i, 1] <- i
  i <- i + 1
}

toplot <- as.data.frame(storelev)
toplot$V1 <- as.numeric(toplot$V1)
toplot$V2 <- as.numeric(toplot$V2)
toplot$V4 <- toplot$V2/sum(toplot$V2)

ggplot(toplot, aes(x=V1, y=V2)) + geom_point() + 
  xlab("Province") + ylab("partial leverage") +
  geom_hline(yintercept = length(reg.out$coefficients)/G, color = "red")  
ggsave('leverage_partial_cl.png', height = 4, width = 6, dpi = 1000)


# 
# sw.correction <- function(H, X, ell, alpha) {
#   N <- nrow(H)      # sample size
#   IH <- 1 - H       # I-H
#   XX <- solve(t(X) %*% X)  # X'X^{-1}
# 
#   # (I-H)*e_i * X_i'(X'X)^(-1)l /sqrt(1-H_{ii})
#   aux <- lapply(c(1:N),
#                 function(i) IH[, i] * c(X[i, ] %*%  XX %*% ell * sqrt(IH[i, i])^(-1)))
#   G <- matrix(unlist(aux), N, N) # transform in matrix
#   GG <- t(G) %*% G
#   GG2 <- GG %^% 2
#   nu <- sum(diag(GG))^2/sum(diag(GG2))
#   correction <- qt(1-alpha/2, df = nu)/qnorm(1-alpha/2)
# 
#   return(list(nu=nu,t=t))
# }
# 
# ell <- as.vector(c(1,rep(0, ncol(X)-1)))
# ell2 <- as.vector(c(0,1,rep(0, ncol(X)-2)))
# aux <- lapply(c(1:N),
#               function(i) IH[, i] * c(X[i, ] %*%  XX %*% ell2 * sqrt(IH[i, i])^(-1)))
# G2 <- matrix(unlist(aux), N, N)
