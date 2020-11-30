#### BIOSTAT 719
#### IRWLS procedure

options(width=90)

#### Data (Table 4.3 on page 67)
y <- c(2,3,6,7,8,9,10,12,15)
x <- c(-1,-1,0,0,0,0,1,1,1)
d1 <- as.data.frame(cbind(y=y,x=x))
d1

## Example: Y ~ Poisson, link function = identity, g(mu) = b0 + b1*x
f.iwls <- function(d,b,niter=1)
{
  X <- cbind(1,d[,"x"])     # Design matrix
  z <- as.matrix(d[,"y"])   # z vector
  
  bb <- b
  for(i in 1:niter) {
    W <- diag(1/(b[1]+b[2]*d[,"x"]))  # W diagonal matrix
    b <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z  # MLE
    
    bb <- rbind(bb, t(b))
  }
  return(list(b=b,Information=t(X)%*%W%*%X, bb=bb))  
}

tmp <- f.iwls(d=d1,b=c(7,5),niter=6)
tmp
solve(tmp$Information)   # Information^{-1}
sqrt(diag(solve(tmp$Information)))  # SE of MLEs

tmp <- f.iwls(d=d1,b=c(1,0),niter=6)
tmp
solve(tmp$Information)   # Information^{-1}
sqrt(diag(solve(tmp$Information)))  # SE of MLEs

#### Run it using glm
res.p <- glm(y ~ x, data=d1, family=poisson(link="identity"))
summary(res.p)
summary(res.p)$cov.unscaled

res.p <- glm(y ~ x, data=d1, family=poisson(link="identity"), 
            epsilon = 1e-20, maxit = 1000)
summary(res.p)
summary(res.p)$cov.unscaled

res.p <- glm(y ~ x, data=d1, family=poisson(link="identity"), 
             epsilon = 1e-15, maxit = 1000)
summary(res.p)
summary(res.p)$cov.unscaled


