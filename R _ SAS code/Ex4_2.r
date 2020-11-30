#### BIOSTAT 719
#### Textbook example 4.2 on page 59

## Set working directory
setwd("YourWorkingDirectory")
## Import dataset
d1 <-  read.csv(file="table4.1.csv",row.names=NULL)
head(d1)
t(d1)

f.weibul <- function(y,theta,lambda)
{
loglik <- sum((lambda-1)*log(y) + log(lambda) - lambda*log(theta)  - (y/theta)^lambda)
U <- sum( -lambda/theta + lambda*(y^lambda)/(theta^(lambda+1)) )
dU <- sum( lambda/(theta^2) - lambda*(lambda+1)*(y^lambda)/(theta^(lambda+2)) )
I <- (lambda^2)*length(y)/(theta^2)
return(list(loglik=loglik, U=U, dU=dU, I=I))
}

theta.seq <- seq(7000,13000,by=100)
loglik <- NULL
U <- NULL
for(i in 1:length(theta.seq)) {
 loglik <- c(loglik, f.weibul(y=d1[,"lifetimes"],theta=theta.seq[i],lambda=2)$loglik)
 U <- c(U, f.weibul(y=d1[,"lifetimes"],theta=theta.seq[i],lambda=2)$U)
}

theta.seq
loglik
U

pdf("Figure4.4.pdf")
plot(x=theta.seq,y=loglik,type="l",ylim=c(-496,-480),ylab="Log-likelihood",xlab=expression(theta))
dev.off()

pdf("plot1C3.pdf")
par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(x=theta.seq,y=loglik,type="l",ylab="Log-likelihood",xlab="theta")
plot(x=theta.seq,y=U,type="l",ylab="U",xlab="theta")
abline(h=0,lty=2)
dev.off()

# MLE
theta.hat <- sqrt(sum(d1[,"lifetimes"]^2)/length(d1[,"lifetimes"]))
theta.hat

# Newton-Raphson & Method of scoring
#theta0 <- mean(d1[,"lifetimes"])
theta0 <- 8805.9  # Textbook pick this as initial value
theta.NR <- theta.MS <- theta0
tmp.NR.iter <- tmp.MS.iter <- NULL
for(i in 1:10) {
  # Newton-Raphson
  one.NR <- f.weibul(y=d1[,"lifetimes"],theta=theta.NR,lambda=2)
  tmp.NR.iter <- rbind(tmp.NR.iter,c(theta.NR=theta.NR, U=one.NR$U,dU=one.NR$dU, I=one.NR$I))
  theta.NR <- theta.NR - one.NR$U/one.NR$dU

  # Method of scoring
  one.MS <- f.weibul(y=d1[,"lifetimes"],theta=theta.MS,lambda=2)
  tmp.MS.iter <- rbind(tmp.MS.iter,c(theta.MS=theta.MS, U=one.MS$U,dU=one.MS$dU, I=one.MS$I))
  theta.MS <- theta.MS + one.MS$U/one.MS$I
}
theta.hat
tmp.NR.iter[1:6,]
tmp.MS.iter[1:6,]











