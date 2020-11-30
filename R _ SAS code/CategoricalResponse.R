#### BIOSTAT 719
#### Example: car preferences
#### Table 8.1 (Textbook page 153) 
sex <- c("women","women","women","men","men","men")
age <- c("18-23","24-40",">40","18-23","24-40",">40")
y1 <- c(26, 9, 5, 40, 17, 8)
y2 <- c(12, 21, 14, 17, 15, 15)
y3 <- c(7, 15, 41, 8, 12, 18)
d2 <- data.frame(sex, age, y1, y2, y3)

d2$men <- ifelse(d2$sex=="men",1,0)
d2$total <- apply(d2[,c("y1","y2","y3")], 1, sum)
d2$y1.perc <- d2$y1/d2$total
d2$y2.perc <- d2$y2/d2$total
d2$y3.perc <- d2$y3/d2$total
d2

install.packages("VGAM")
library(VGAM)
#### Nominal logistic regression
## Reference: y1 (No or little importance)
## Order y2 and y3 according to the needed reference category of response
nominal.fit1 <- vglm(cbind(y2,y3,y1)~men, family=multinomial, data=d2)
summary(nominal.fit1)
confint(nominal.fit1)
exp(confint(nominal.fit1))
## Predicted probability
fitted(nominal.fit1)

#### Ordinal logistic regression
## Cumulative logit model
ord.cumul.fit <- vglm(cbind(y1,y2,y3)~men, family=cumulative(parallel=F), data=d2)
ord.cumul.fit <- vglm(cbind(y3,y2,y1)~men, family=cumulative(parallel=F), data=d2)
summary(ord.cumul.fit)
fitted(ord.cumul.fit)
## Proportional odds model
ord.prop.fit <- vglm(cbind(y1,y2,y3)~men, family=cumulative(parallel=T), data=d2)
summary(ord.prop.fit)
fitted(ord.prop.fit)



