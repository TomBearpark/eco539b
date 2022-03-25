library(haven)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(labelled)
library(sandwich)
library(MASS)
library(clubSandwich)
library(stargazer)

setwd("D:\\Dropbox\\Università\\PhD\\II Year\\Spring\\ECO519 - Non-linear Econometrics\\psets\\ps2")
data <- haven::read_dta("ak91.dta")
data.use <- subset(data, census == 1980 & cohort == 2)

############################################################################
## Exercise 1
# Estimate regression of earnings on a constant and education level
earn <- lm(lwage ~ 1 + educ, data = data.use)

# Clustered standard errors (Liang and Scott, 1986)
clSE <- sandwich::vcovCL(earn, cluster = data.use$SOB, type = 'HC0')

# EHW standard errors
EHW <- sandwich::vcovHC(earn, type = 'HC0')

stargazer(earn, earn, earn, se = list(NULL, sqrt(diag(EHW)), sqrt(diag(clSE))))

# Compute leverage
hats <- as.data.frame(hatvalues(earn))
theme_set(theme_bw())
ggplot(hats) + geom_histogram(aes(x=hatvalues(earn))) + xlab("Leverage")
ggsave('lvg_all.png', height = 4, width = 6, dpi = 1000)

# Compute within group leverage
sobs <- unique(data.use$SOB)
G <- length(sobs)
labs <- val_labels(sobs)
storelev <- matrix(NA, G, 3) 
X <- cbind(1,data.use$educ)
XX <- solve(t(X)%*%X)

i <- 1
for (sob in sobs) {
  datacl <- subset(data.use, data.use$SOB == sob)
  Xg <- cbind(1,datacl$educ)
  XXg <- t(Xg)%*%Xg
  
  storelev[i, 2] <- sum(diag(XXg %*% XX)) # compute nuclear norm
  storelev[i, 1] <- sob
  storelev[i, 3] <- names(labs)[labs == sob]
  i <- i + 1
}

toplot <- as.data.frame(storelev)
toplot$V1 <- as.numeric(toplot$V1)
toplot$V2 <- as.numeric(toplot$V2)
toplot$V4 <- toplot$V2/sum(toplot$V2)

ggplot(toplot, aes(x=V1, y=V4, label=V3)) + geom_point() + 
  xlab("State") + ylab("Leverage Share") +
  geom_hline(yintercept = 1/G, color = "red") + 
  geom_text(hjust=0, vjust=-1/2) + 
  scale_y_continuous(labels = scales::percent)
ggsave('lvg_cls.png', height = 4, width = 6, dpi = 1000)


############################################################################
## Exercise 2
data.use$NJindic <- data.use$SOB == 34
NJimpact <- lm(lwage ~ NJindic, data = data.use)

ehw <- sandwich::vcovHC(NJimpact, type = 'HC0')
lz.meat <- sandwich::meatCL(NJimpact, cluster = data.use$SOB, type = 'HC0')
lz <- sandwich::vcovCL(NJimpact, cluster = data.use$SOB, type = 'HC0')
bm <- sandwich::vcovCL(NJimpact, cluster = data.use$SOB, type = 'HC2')

stargazer(NJimpact, NJimpact, NJimpact, NJimpact,
          se = list(NULL, sqrt(diag(ehw)), sqrt(diag(lz)), sqrt(diag(bm))))


