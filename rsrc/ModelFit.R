library(ggplot2)
library(sn)
library(fGarch)
data('faithful')

setwd("C:\\Users\\progi\\Desktop\\Projects\\pXg\\ModelFitting")
peaksData <- read.csv(file = "TestData.txt", header = T, sep = "\t", as.is = as.double())

U <- function(mu, y) {
  e1 <- (y-mu[1])^2
  e2 <- (y-mu[2])^2
  e.min <- pmin(e1, e2)
  e <- sum(e.min)
  return (e)
}


y <- peaksData$ALC....
hist(y, xlab="ALC")

r.cluster <- nlm(U, c(-50, 50), y)
mu.est <- r.cluster$estimate


print(mu.est)

e1 <- (y-mu.est[1])^2
e2 <- (y-mu.est[2])^2

group.est <- rep(1,length(y))
group.est[which(e2<e1)] <- 2
group.est <- as.factor(group.est)

theme_set(theme_bw())
ggplot() + geom_point(data=data.frame(y,group.est), aes(x=y,y=0,colour=group.est),size=3) + 
  geom_point(data=data.frame(y=mu.est,group=factor(c(1,2))), aes(y,0,colour=group),size=10, shape="x") + 
  geom_vline(xintercept=mean(mu.est))

mixt.deviance <- function(theta,y) {
  pi1    <- theta[1] / (theta[1]+theta[2])
  pi2    <- theta[2] / (theta[1]+theta[2])
  mu1    <- theta[3]
  mu2    <- theta[4]
  sigma1 <- theta[5]
  sigma2 <- theta[6]
  skewness1 <- theta[7]
  skewness2 <- theta[8]
  pdf <- pi1*dnorm(y,mu1,sigma1) + pi2*dnorm(y,mu2,sigma2)
  #pdf <- pi1*dsn(y, xi = mu1, omega = sigma1, alpha = skewness1) + pi2*dsn(y, xi = mu2, omega = sigma2, alpha = skewness2)
  #pdf <- pi1*dsnorm(y, mean = mu1, sd = sigma1, xi = skewness1) + pi2*dsnorm(y, mean = mu2, sd = sigma2, xi = skewness2)
  deviance <- -2*sum(log(pdf))
  return (deviance)
}

r.nlm <- nlm(mixt.deviance,c(0.1,10000, 20,80, 10,10, 0.1,0.1),y)
theta.est <- r.nlm$estimate
print(matrix(theta.est,nrow=4,byrow = T))

dmixt <- function(x,theta) {
  pi1    <- theta[1] / (theta[1]+theta[2])
  pi2    <- theta[2] / (theta[1]+theta[2])
  mu1    <- theta[3]
  mu2    <- theta[4]
  sigma1 <- theta[5]
  sigma2 <- theta[6]
  skewness1 <- theta[7]
  skewness2 <- theta[8]
  f1 <- pi1*dnorm(x,mu1,sigma1)
  #f2 <- dnorm(x,mu2,sigma2)
  f2 <- pi2*dnorm(x,mu2,sigma2)
  f <- f1+f2
  #f <- f2
}
x <- (1:100)
pdf.mixt <- dmixt(x,theta.est)
ggplot(data=peaksData) + geom_histogram(aes(x=ALC...., y=..density..), bins=20) +
  geom_line(data=data.frame(x,pdf.mixt), aes(x,pdf.mixt),colour="red",size=1.5)
