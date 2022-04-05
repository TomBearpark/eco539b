setwd("D:/Dropbox/Università/PhD/II Year/Spring/ECO519 - Non-linear Econometrics/psets/ps3")

# Load stuff
library(haven)
library(ivreg)

library(sandwich)
library(dfadjust)
library(xtable)
library(tidyverse)
library(reshape2)

theme_set(theme_bw())

CI.store <- matrix(NA, 3, 6)

data <- haven::read_dta("ak91.dta")

# Extract subsets of data
data.pac <- subset(data, census == 1980 & cohort == 2 & division == 9)
data.atl <- subset(data, census == 1980 & cohort == 2 & division == 2)

# Create date of birth instrument
data.pac$birthq1 <- 1*(data.pac$age == floor(data.pac$age))
data.atl$birthq1 <- 1*(data.atl$age == floor(data.atl$age))

# IV regression
iv.pac <- ivreg(lwage ~ 1 + educ | birthq1, data = data.pac)
se.pac <- sandwich::vcovHC(iv.pac, type = 'HC1') 
iv.atl <- ivreg(lwage ~ 1 + educ | birthq1, data = data.atl)
se.atl <- sandwich::vcovHC(iv.atl, type = 'HC1') 

# Wald Confidence Interval
b.iv.pac <- iv.pac$coefficients[2]
b.iv.atl <- iv.atl$coefficients[2]

CI.store[1, ] <- c(b.iv.pac - qnorm(0.975)*sqrt(se.pac[2,2]),
                   b.iv.pac,
                   b.iv.pac + qnorm(0.975)*sqrt(se.pac[2,2]),
                   b.iv.atl - qnorm(0.975)*sqrt(se.atl[2,2]),
                   b.iv.atl,
                   b.iv.atl + qnorm(0.975)*sqrt(se.atl[2,2]))

# Anderson-Rubin Confidence Interval
beta.grid.pac <- seq(from=-0.5, to=0.5, by=0.001)
beta.grid.atl <- seq(from=-50, to=50, by = 1)
ar.pv.pac <- matrix(NA, nrow = length(beta.grid.pac), ncol = 2)
ar.pv.atl <- matrix(NA, nrow = length(beta.grid.atl), ncol = 2)
ar.ts.pac <- ar.pv.pac
ar.ts.atl <- ar.pv.atl

i <- 1

for (b in beta.grid.pac) {
  # residualize outcome
  data.pac$lwage.res <- data.pac$lwage - b*data.pac$educ 

  # run residualized reduced form
  ar.pac <- lm(lwage.res ~ 1 + birthq1, data = data.pac)

  # compute robust variance-covariance
  vc.pac <- sandwich::vcovHC(ar.pac, type = 'HC1')

  # test coefficient on instrument
  test.pac <- ar.pac$coefficients[2]/sqrt(vc.pac[2,2])
  ar.ts.pac[i,1] <- test.pac
  
  # retrieve p-value
  ar.pv.pac[i,1] <- pnorm(abs(test.pac), lower.tail = FALSE)*2  
  i <- i + 1
}

ar.pv.pac[,2] <- "Pacific"

i <- 1

for (b in beta.grid.atl) {
  # residualize outcome
  data.atl$lwage.res <- data.atl$lwage - b*data.atl$educ 
  
  # run residualized reduced form
  ar.atl <- lm(lwage.res ~ 1 + birthq1, data = data.atl)
  
  # compute robust variance-covariance
  vc.atl <- sandwich::vcovHC(ar.atl, type = 'HC1')
  
  # test coefficient on instrument
  test.atl <- ar.atl$coefficients[2]/sqrt(vc.atl[2,2])
  ar.ts.atl[i,1] <- test.atl
  
  # retrieve p-value
  ar.pv.atl[i,1] <- pnorm(abs(test.atl), lower.tail = FALSE)*2
  i <- i + 1
}

ar.pv.atl[,2] <- "mid-Atlantic"


toplot <- data.frame(beta = c(beta.grid.pac, beta.grid.atl), 
                     pv = rbind(ar.pv.pac, ar.pv.atl))
colnames(toplot) <- c("beta", "pvalue", "division")
toplot$pvalue <- as.numeric(toplot$pvalue)


ggplot(subset(toplot, division=="Pacific"), aes(x=beta, y=pvalue)) + 
  geom_line() + geom_hline(yintercept=0.05, color = "red") +
  geom_vline(xintercept = b.iv.pac, color = "black", linetype="dashed") + 
  ylab("AR p-value")
ggsave('AR_pacific.png', height = 4, width = 6, dpi = 1000)

ggplot(subset(toplot, division=="mid-Atlantic"), aes(x=beta, y=pvalue)) + 
  geom_line() + geom_hline(yintercept=0.05, color = "red") +
  geom_vline(xintercept = b.iv.pac, color = "black", linetype="dashed") + 
  ylab("AR p-value")
ggsave('AR_atlantic.png', height = 4, width = 6, dpi = 1000)

ar.ci.pac <- beta.grid.pac[as.numeric(ar.pv.pac[,1]) >= 0.05]

CI.store[2, ] <- c(ar.ci.pac[1],
                   b.iv.pac,
                   ar.ci.pac[length(ar.ci.pac)],
                   -Inf,
                   b.iv.atl,
                   Inf)

# tF procedure

# First stage
fs.pac <- lm(educ ~ 1 + birthq1, data = data.pac)
fs.atl <- lm(educ ~ 1 + birthq1, data = data.atl)
fs.vc.pac <- sandwich::vcovHC(fs.pac, type ='HC1')
fs.vc.atl <- sandwich::vcovHC(fs.atl, type ='HC1')

t.pac <- fs.pac$coefficients[2]/sqrt(vcovHC(fs.pac, type = 'HC1')[2,2])
f.pac <- t.pac^2
t.atl <- fs.atl$coefficients[2]/sqrt(vcovHC(fs.atl, type = 'HC1')[2,2])
f.atl <- t.atl^2

CI.store[3, ] <- c(-0.13506071, b.iv.pac, 0.39288127,
                   -Inf, b.iv.atl, Inf)

# Store table
tab <- as.table(rbind(CI.store, c(NA,f.pac,NA,NA,f.atl,NA)))
rownames(tab) <- c("Wald CI", "Anderson-Rubin CI", "tF CI", "First-stage F")
xtable(tab)