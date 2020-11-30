#### BIOSTAT 719
#### Embryogenic anther data

setwd("YourWorkingDirectory")

d1 <- read.csv("table7.5.csv")
d1[,"newstor"] <- d1[,"storage"]-1
d1[,"force"] <- log(d1[,"centrifuge"])
fit1 <- glm( cbind(y,n-y) ~ newstor, data=d1,family=binomial(link="logit"))
summary(fit1)
fit2 <-glm( cbind(y,n-y) ~ newstor*force, data=d1,family=binomial(link="logit"))
summary(fit2)

anova(fit1, fit2)
anova(fit1, fit2, test="Chisq")
anova(fit1, fit2, test="LRT")

round(summary(fit2)$cov.unscaled, digit=3)

