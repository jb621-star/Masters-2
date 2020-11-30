#### BIOSTAT 719
#### Mixed effects model
#### Table 11.1 Storke data

library(reshape2)  ## To reshape data
library(ggplot2)
library(lme4)  ## lmer
library(nlme)  ## lme
library(mvtnorm)
library(dplyr)

setwd("YourWorkingDirectory")

#### Wide vs. long data
dt.w <- read.csv("table11_1.csv")
head(dt.w)
dt.l <- melt(dt.w, id.vars=c("id", "group"),
             measure.vars=paste("week",1:8,sep=""),
             variable.name="week",
             value.name="fas")
dt.l$time <- recode(dt.l$week, "week1"=1, "week2"=2, "week3"=3, "week4"=4,
                    "week5"=5, "week6"=6, "week7"=7, "week8"=8)
dt.l <- dt.l[order(dt.l$id),]
head(dt.l)

#### Spaghetti plot
p <- ggplot(data=dt.l, aes(x=time, y=fas, group=id, colour=group))
p + geom_line()

#### Estimated correlations
round(cor(dt.w[, paste("week",1:8,sep="")]), 2)

#### Correlation matrix plot
pairs(dt.w[,paste("week",1:8,sep="")], pch=19, cex=0.4)

#### Fit models 
## Model 1: Random intercept model
fit1 <- lmer(fas ~ time + (1 | id), data=dt.l)
#fit12 <- lme(fas ~ time,  
#             random = ~1 | id, 
#             data=dt.l)
summary(fit1)
Zi <- cbind(rep(1,8)) 
G <- matrix(393.8,1,1)
sigmasq.e <- 80.3
Vi <- Zi%*%G%*%t(Zi) + sigmasq.e * diag(8)
Vi
round(sweep(sweep(Vi,1,sqrt(diag(Vi)),"/"),2,sqrt(diag(Vi)),"/"),3)
393.8/(393.8 + 80.3)
VarCorr(fit1)
extractAIC(fit1)
fixef(fit1)
confint(fit1)
# predicted value
cbind(dt.l, predict(fit1))[1:16,]

## Model 2: Random intercept and slope model (random intercept and slope are indep)
fit2 <- lmer(fas ~ time + (1 | id) + (0 + time | id), data=dt.l)
summary(fit2)
extractAIC(fit2)
Zi <- cbind(rep(1,8), seq(1:8)) 
G <- matrix(c(392.129,0,0,8.938),2,2)
sigmasq.e <- 26.981
Vi <- Zi%*%G%*%t(Zi) + sigmasq.e * diag(8)
Vi
round(sweep(sweep(Vi,1,sqrt(diag(Vi)),"/"),2,sqrt(diag(Vi)),"/"),3)
confint(fit2)

## Model 3: Random intercept and slope model (With unstructured)
fit3 <- lmer(fas ~ time + (1 + time| id), data=dt.l)
summary(fit3)
extractAIC(fit3)
Zi <- cbind(rep(1,8), seq(1:8)) 
G <- matrix(c(405.101,(-0.35*sqrt(405.101)*sqrt(9.239)),(-0.35*sqrt(405.101)*sqrt(9.239)),9.239),2,2)
sigmasq.e <- 26.846
Vi <- Zi%*%G%*%t(Zi) + sigmasq.e * diag(8)
Vi
round(sweep(sweep(Vi,1,sqrt(diag(Vi)),"/"),2,sqrt(diag(Vi)),"/"),3)

## Compare model fit
cbind(extractAIC(fit1), extractAIC(fit2), extractAIC(fit3))

