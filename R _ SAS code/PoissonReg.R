#### BIOSTAT 719
#### Poisson regression

setwd("YourWorkingDirectory")

## Table 2.1 (textbook page 20)
## x=1 for town; x=0 for country
## y= numbers of chronic conditions 
##    from 26 town women and 23 country women
x <- c(rep(1,26), rep(0,23))
y <- c(0,1,1,0,2, 3,0,1,1,1, 1,2,0,1,3, 0,1,2,1,3, 3,4,1,3,2, 0,
       2,0,3,0,0, 1,1,1,1,0, 0,2,2,0,1, 2,0,0,1,1, 1,0,2)
d3 <- as.data.frame(cbind(x,y))
head(d3)

## Poisson regression with ungrouped data
fit1 <- glm(y~x, data=d3, family=poisson(link="log"))
summary(fit1)

fit1.quasi <- glm(y~x, data=d3, family=quasipoisson)
summary(fit1.quasi)
sum(residuals(fit1.quasi, type = "pearson")^2)/47

## Grouped data1
## y = number of total chronic conditions of women in town and country
tapply(d3[,"y"], d3[,"x"], sum)
d4 <- as.data.frame(cbind(y=c(37,21), x=c(1,0)))
d4
fit2 <- glm(y~x, data=d4, family=poisson(link="log"))
summary(fit2)
round(fit1$coef, 3)
summary(fit1)$coef
## Grouped data2 with sample size
## y = number of total chronic conditions of women in town and country
d5 <- as.data.frame(cbind(d4, n=c(26,23)))
d5
fit3 <- glm(y~x+offset(log(n)), data=d5, family=poisson(link="log"))
summary(fit3)

sum(residuals(fit1, type = "pearson")^2)
summary(fit1)

#### Table 9.1 (textbook page 169)
d1 <- read.csv("table9_1.csv", header=T)
d1$agecat <- as.numeric(d1$age)
d1$smoke <- ifelse(d1$smoking=="smoker", 1, 0)
d1$agesq <- d1$agecat^2
d1$smkage <- d1$smoke * d1$agecat
d1

#### Scatter plots
d1$unitdeath <- (d1$deaths/d1$person.years)*100000

plot(d1$agecat, d1$unitdeath, col=ifelse(d1$smoke==1, 2, 1), 
     pch=ifelse(d1$smoke==1, 18, 19), cex=1.5,
     xlab="Age", ylab="Death per 100,000 person years", bty="L")
legend("topleft", c("Smoker", "Non-smoker"), col=c(2,1), pch=c(18,19))

plot(d1$agecat, log(d1$unitdeath), col=ifelse(d1$smoke==1, 2, 1), 
     pch=ifelse(d1$smoke==1, 18, 19), cex=1.5,
     xlab="Age", ylab="Log of death per 100,000 person years", bty="L")
legend("topleft", c("Smoker", "Non-smoker"), col=c(2,1), pch=c(18,19))

## Check of lineariyy ot the agecat variable
for(i in sort(unique(d1$agecat))) {
  d1[,paste("a",i,sep="")] <- ifelse(d1[,"agecat"]==i,1,0)
}
d1

fit <- glm(deaths ~ agecat + a3 + a4 + a5 + offset(log(person.years)), data=d1, family=poisson)
anova(fit, test="Chisq")
1-pchisq(24.38+22.77+13.87, df=3)  ## H0: gamma3=gamma4=gamma5=0

## Final model (see text page 169)
## Poisson regression
fit <- glm(deaths ~ agecat + agesq + smoke + smkage + offset(log(person.years)), 
           data=d1, family=poisson)
summary(fit)
anova(fit,"Chisq")

## Calculate pearson chi-square 
sum(residuals(fit, type = "pearson")^2)/6

round(exp(fit$coef[-1]), 2)

## Pattern 1: 60 years old smoker
## Pattern 2: 50 years old non-smoker
## L^T = (0, 1, 5, 1, 3)
L <- matrix(c(0,1,5,1,3), nrow=5,byrow=T)
lrr <- t(L) %*% as.matrix(fit$coef)
rr <- exp(lrr)
covmat <- summary(fit)$cov.unscaled
var.lrr <- t(L)%*%covmat%*%L
ci.lrr <- c(lrr-(1.96*sqrt(var.lrr)), lor+(1.96*sqrt(var.lrr))) 
ci.rr <- exp(ci.lrr)
c(rr, ci.rr)

#### Adjusting for over or underdispersion
#### 1. Quasi-Poisson regression
fit.quasi <- glm(deaths ~ agecat + agesq + smoke + smkage + offset(log(person.years)), 
                 data=d1, family=quasipoisson)
summary(fit.quasi)

summary()

sum(residuals(fit, type = "pearson")^2)/5
## Dispersion parameter = 0.31
## This parameter tells us how many times larger the variance is than the mean. 
## Since our dispersion was (larger) less than one, 
## it turns out the conditional variance is actually (larger) smaller than the conditional mean. 
## We have (over)under-dispersion.
summary(fit)$coef
summary(fit.quasi)$coef
## SE of intercept changes from 0.45 to 0.25. 0.45 * sqrt(dispersion) = 0.25!
## In SAS, Dispersion is estimated by the sqrt of pearson's chi-square/DOF

#### 2. Negative binomial regression
library(MASS)
fit.nb <- glm.nb(deaths ~ agecat + agesq + smoke + smkage + 
                   offset(log(person.years)), data=d1)
summary(fit.nb)
#equire(foreign)
#dat <- read.dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")
#head(dat)
## Check model assumption
# As we mentioned earlier, negative binomial models assume the conditional means are not equal to the conditional variances. 
#This inequality is captured by estimating a dispersion parameter (not shown in the output) that is held constant in a Poisson model. 
#Thus, the Poisson model is actually nested in the negative binomial model. 
#We can then use a likelihood ratio test to compare these two and test this model assumption. 
logLik(fit.nb)
logLik(fit)
2*(logLik(fit.nb)-logLik(fit))
pchisq(2 * (logLik(fit.nb)-logLik(fit)), df = 1, lower.tail = FALSE)
# In this example the associated chi-squared value estimated from 
# 2*(logLik(m1) – logLik(m3)) is -0.0004893766 with one degree of freedom. 
# This strongly suggests the negative binomial model, estimating the dispersion parameter, 
# is not appropriate than the Poisson model.

#### Fitted values
poi.pred <- fitted.values(fit)
poi.pred.unit <- (poi.pred/d1$person.years)*100000

plot(1:10, d1$unitdeath, pch=19,  
     xlab="Covariate patterns", ylab="Death per 100,000 person years", bty="L")
points(1:10, poi.pred.unit, pch=19, col=2)
legend("topleft", c("Observed", "Fitted"), col=c(1,2), pch=c(19,19))





