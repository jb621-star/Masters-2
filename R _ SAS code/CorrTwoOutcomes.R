#### BIOSTAT 719
#### Correlated outcomese

library(mvtnorm)

#### Simulated data with two groups
#### Generate 2000 observations in each group
####
## Case 1. Same covariate value within a cluster
set.seed(10)
sigma2 <- diag(c(2,2)) + 1 
sigma2  # Positive correlation (corr=1/3) within a cluster
# group 1
tmp.grp1 <- rmvnorm(n=1000, mean=c(4,4), sigma=sigma2, pre0.9_9994 = T)
head(tmp.grp1)
# Sample mean
apply(tmp.grp1, 2, mean)
# Sample covariance matrix
#round(cov(tmp.grp1), 2)
cov(tmp.grp1)
# Sample correlation matrix
#round(cor(tmp.grp1), 2)
cor(tmp.grp1)
#
# group 0
tmp.grp0 <- rmvnorm(n=1000, mean=c(4.2,4.2), sigma=sigma2, pre0.9_9994 = T)
head(tmp.grp0)
# Sample mean
apply(tmp.grp0, 2, mean)
# Sample covariance matrix
#round(cov(tmp.grp0), 2)
cov(tmp.grp0)
# Sample correlation matrix
#round(cor(tmp.grp0), 2)
cor(tmp.grp0)

# To estimate a common correlation
cor(rbind(tmp.grp1, tmp.grp0))
rho <- cor(rbind(tmp.grp1, tmp.grp0))[1,2]
rho

# Assume all the observations are independent (incorrect assumption)
t.test(as.vector(tmp.grp1), as.vector(tmp.grp0), var.equal=T)
# Let's incorporate the correlation
stat1 <- -2.4684/sqrt(1+rho)
stat1
pval1 <- 2*(1-pnorm(abs(stat1)))
pval1

## Case 2. covariate value changes within a cluster
set.seed(10)
sigma2 <- diag(c(2,2)) + 1 
sigma2  # Positive correlation (corr=1/3) within a cluster
tmp.grp <- rmvnorm(n=2000, mean=c(4, 4.2), sigma=sigma2, pre0.9_9994 = T)
head(tmp.grp)
# Sample mean
apply(tmp.grp, 2, mean)
# Sample correlation matrix
cor(tmp.grp)
rho <- cor(tmp.grp)[1,2]
# Assume all the observations are independent (incorrect assumption)
t.test(tmp.grp[,1], tmp.grp[,2], var.equal=T)
# Correctly assume and run paired t-test
t.test(tmp.grp[,1], tmp.grp[,2], paired=T)
#
stat2 <- -2.8044/sqrt(1-rho)
stat2


