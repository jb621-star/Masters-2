#### BIOSTAT 719
#### Textbook example 1.6.5 on page 15

## Set working directory
setwd("YourWorkingDirectory")
## Import dataset
d1 <-  read.csv(file="table1.2.csv", header=TRUE, row.names=1)
d1

## Calculate mean
theta.hat <- sum(d1[,"number"]) / length(d1[,"number"])
theta.hat

## Log likelihood function for Poisson distribution 
## Ignore constant term
f.loglik <- function(parameter,y) {
  out <- sum(y*log(parameter) - parameter)
  return(out)
}

theta <- seq(3,8,by=0.1)
loglik <- NULL
for(i in 1:length(theta)) loglik <- c(loglik,f.loglik(parameter=theta[i],y=d1[,"number"]))

theta
loglik

pdf("Fig1_2.pdf")
plot(x=theta,y=loglik,type="l",ylim=c(35,55),ylab="Log-likelihood", xlab=expression(theta))
dev.off()

pdf("LogLik_U.pdf")
par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(x=theta.seq,y=loglik,type="l",ylab="Log-likelihood",xlab=expression(theta))
abline(v=theta.hat,lty=2)
plot(x=theta.seq,y=U,type="l",ylab="U",xlab=expression(theta))
abline(v=theta.hat,lty=2)
abline(h=0,lty=2)
dev.off()

#### Newton-Raphson algorithm

f.poisson <- function(y,theta)
{
  n <- length(y)
  loglik <- sum(y)*log(theta) - n*theta
  U <- (1/theta)*sum(y) - n
  dU <- -(1/theta^2)*sum(y)
  return(list(loglik=loglik, U=U, dU=dU))
}

theta.seq <- seq(3,8,by=0.1)
loglik <- NULL
U <- NULL
for(i in 1:length(theta.seq)) {
  loglik <- c(loglik, f.poisson(y=d1[,"number"],theta=theta.seq[i])$loglik)
  U <- c(U, f.poisson(y=d1[,"number"],theta=theta.seq[i])$U)
}

theta
loglik
U

theta0 <- 3.5
tmp.iter <- NULL
for(i in 1:10) {
  one <- f.poisson(y=d1[,"number"],theta=theta0)
  tmp.iter <- rbind(tmp.iter,c(theta=theta0,U=one$U,dU=one$dU))
  theta0 <- theta0 - one$U/one$dU
}
tmp.iter
theta.hat

#### Plots
## Step 0
theta0 <- 3.5
res0 <- f.poisson(y=d1[,"number"],theta=theta0)

pdf("NR_step0.pdf")
par(mfrow=c(1,1),mar=c(5,5,1,1), bty="l")
plot(x=theta.seq,y=U,type="l",ylab="U",lwd=3,xlab=expression(theta))
abline(h=0,lty=2)
abline(v=theta0,lty=2)
points(theta.hat, 0, pch=10, col=2, cex=2)
dev.off()

## Step 1
pdf("NR_step1.pdf")
par(mfrow=c(1,1),mar=c(5,5,1,1), bty="l")
plot(x=theta.seq,y=U,type="l",ylab="U",lwd=3,xlab=expression(theta))
abline(h=0,lty=2)
abline(v=theta0,lty=2)
points(theta.hat, 0, pch=10, col=2, cex=2)

dU.slope <- res0$dU
dU.intercept <- res0$U - dU.slope*theta0
abline(a=dU.intercept, b=dU.slope)
dev.off()

## Step 2
theta1 <- theta0 - res0$U/res0$dU
res1 <- f.poisson(y=d1[,"number"],theta=theta1)

pdf("NR_step2.pdf")
par(mfrow=c(1,1),mar=c(5,5,1,1), bty="l")
plot(x=theta.seq,y=U,type="l",ylab="U",lwd=3,xlab=expression(theta))
abline(h=0,lty=2)
abline(v=theta0,lty=2)
points(theta.hat, 0, pch=10, col=2, cex=2)

dU.slope <- res0$dU
dU.intercept <- res0$U - dU.slope*theta0
abline(a=dU.intercept, b=dU.slope)

abline(v=theta1,lty=2, col=3)
dU.slope <- res1$dU
dU.intercept <- res1$U - dU.slope*theta1
abline(a=dU.intercept, b=dU.slope, col=3)
dev.off()

## Step 3
theta2 <- theta1 - res1$U/res1$dU
res2 <- f.poisson(y=d1[,"number"],theta=theta2)

pdf("NR_step3.pdf")
par(mfrow=c(1,1),mar=c(5,5,1,1), bty="l")
plot(x=theta.seq,y=U,type="l",ylab="U",lwd=3,xlab=expression(theta))
abline(h=0,lty=2)
abline(v=theta0,lty=2)
points(theta.hat, 0, pch=10, col=2, cex=2)

dU.slope <- res0$dU
dU.intercept <- res0$U - dU.slope*theta0
abline(a=dU.intercept, b=dU.slope)

abline(v=theta1,lty=2, col=3)
dU.slope <- res1$dU
dU.intercept <- res1$U - dU.slope*theta1
abline(a=dU.intercept, b=dU.slope, col=3)

abline(v=theta2,lty=2, col=4)
dU.slope <- res2$dU
dU.intercept <- res2$U - dU.slope*theta2
abline(a=dU.intercept, b=dU.slope, col=4)
dev.off()

## Step 4
theta3 <- theta2 - res2$U/res2$dU
res3 <- f.poisson(y=d1[,"number"],theta=theta3)

pdf("NR_step4.pdf")
par(mfrow=c(1,1),mar=c(5,5,1,1), bty="l")
plot(x=theta.seq,y=U,type="l",ylab="U",lwd=3,xlab=expression(theta))
abline(h=0,lty=2)
abline(v=theta0,lty=2)
points(theta.hat, 0, pch=10, col=2, cex=2)

dU.slope <- res0$dU
dU.intercept <- res0$U - dU.slope*theta0
abline(a=dU.intercept, b=dU.slope)

abline(v=theta1,lty=2, col=3)
dU.slope <- res1$dU
dU.intercept <- res1$U - dU.slope*theta1
abline(a=dU.intercept, b=dU.slope, col=3)

abline(v=theta2,lty=2, col=4)
dU.slope <- res2$dU
dU.intercept <- res2$U - dU.slope*theta2
abline(a=dU.intercept, b=dU.slope, col=4)

abline(v=theta3,lty=2, col=2)
dU.slope <- res3$dU
dU.intercept <- res3$U - dU.slope*theta3
abline(a=dU.intercept, b=dU.slope, col=2)
dev.off()






