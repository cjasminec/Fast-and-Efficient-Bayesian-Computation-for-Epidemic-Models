library("tidyverse")
library("ggplot2")
library("maxLik")

eulerSIR <- function(x0, theta, deltat, T){
  s=T/deltat  #number of steps
  iterations <- matrix(0,s+1,2)
  #put in initial values
  iterations[1,] <- x0
  #define f(S,I)
  S <- x0[1]
  I <- x0[2]
  for (i in 2:(s+1)) {
    #calculate
    dS <- -theta[1]*S*I
    dI <- theta[1]*S*I-theta[2]*I
    iterations[i,1] <- iterations[i-1,1]+ deltat*dS
    iterations[i,2] <- iterations[i-1,2]+deltat*dI
    S <- iterations[i,1]
    I <- iterations[i,2]
  }
  return(iterations)
}


dargatz <- eulerSIR(c(762,1),c(exp(-6),0.5),0.01,14)
colnames(dargatz) <- c("Susceptibles", "Infectives")
plot(ts(dargatz, start=0, deltat=0.01), main="Boarding School Influenza Dataset")

set.seed(123)
#ADDITIVE GAUSSIAN NOISE
add.loglike_epi <- function(theta, y, deltat, x0){
  N <- length(y)
  inc <- 1/deltat
  output <- eulerSIR(x0, exp(c(theta[1],theta[2])), deltat, N-1)
  index <- 1+(0:(N-1))*inc
  output <- output[index, 2]
  loglikeindiv <- dnorm(y, output, exp(theta[3]), log=TRUE)
  loglike <- sum(loglikeindiv)
  return (loglike)
}

Its.add <- rnorm(15, dargatz[1+(0:14)*100, 2],3)
loglike_epi(c(-6,log(0.5), 0), Its.add,0.01,c(762,1))

test1 <- maxLik(add.loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its.add,deltat=0.01,x0=c(762,1))
test.add <- logLik(test1)
add <- coef(test1)

add.pred <- eulerSIR(c(762,1), c(exp(add[1]), exp(add[2])),0.01, 14)

plot(ts(Its.add))

lq.add <- qnorm(0.025, add.pred[,2], exp(add[3]))
uq.add <- qnorm(0.975, add.pred[,2], exp(add[3]))

plot(ts(add.pred[,2],start=0,deltat=0.01), type="l",ylim=range(c(lq.add, uq.add)), lwd=0.6)
lines(ts(lq.add,start=0,deltat=0.01), lty=2, lwd=0.6, col=2)
lines(ts(uq.add,start=0,deltat=0.01), lty=2, lwd=0.6, col=2)


#MULTIPLICATIVE GAUSSIAN NOISE
mult.loglike_epi <- function(theta, y, deltat, x0){
  N <- length(y)
  inc <- 1/deltat
  output <- eulerSIR(x0, exp(c(theta[1],theta[2])), deltat, N-1)
  index <- 1+(0:(N-1))*inc
  output <- output[index, 2]
  loglikeindiv <- dlnorm(y, log(output), exp(theta[3]), log=TRUE)
  loglike <- sum(loglikeindiv)
  return (loglike)
}


#this is where noise comes in
#Its.mult <- exp(rnorm(15,log(dargatz[1+(0:14)*100, 2]),1))
Its.mult <- rlnorm(15, log(dargatz[1+(0:14)*100, 2]),1)
mult.loglike_epi(c(-6,log(0.5), 0), Its.mult,0.01,c(762,1))
#huge difference between diff ways of getting Its.mult - why? which is correct?

test2 <- maxLik(mult.loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its.mult,deltat=0.01,x0=c(762,1))
test.mult <- logLik(test2)
mult <- coef(test2)

mult.pred <- eulerSIR(c(762,1), c(exp(mult[1]), exp(mult[2])),0.01, 14)

plot(ts(Its.mult))

lq.mult <- qlnorm(0.025, log(mult.pred[,2]), exp(mult[3]))
uq.mult <- qlnorm(0.975, log(mult.pred[,2]), exp(mult[3]))

plot(ts(mult.pred[,2],start=0,deltat=0.01), ylim=c(0,1000), ylab="infectives")
lines(ts(lq.mult,start=0,deltat=0.01), col=2)
lines(ts(uq.mult,start=0,deltat=0.01), col=2)

#OTHER NOISE (Decoupled Poisson): Yt|It~ N(It, sigma^2*It)
selfmult.loglike_epi <- function(theta, y, deltat, x0){
  N <- length(y)
  inc <- 1/deltat
  output <- eulerSIR(x0, exp(c(theta[1],theta[2])), deltat, N-1)
  index <- 1+(0:(N-1))*inc
  output <- output[index, 2]
  loglikeindiv <- dnorm(y, output, exp(theta[3]) * sqrt(output), log = TRUE)
  loglike <- sum(loglikeindiv)
  return (loglike)
}
Its.selfmult <- rnorm(15,dargatz[1+(0:14)*100, 2],sqrt(dargatz[1+(0:14)*100, 2])*1)
loglike_epi(c(-6,log(0.5), 0), Its.selfmult,0.01,c(762,1))

test3 <- maxLik(selfmult.loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its.selfmult,deltat=0.01,x0=c(762,1))
test.selfmult <- logLik(test3)
selfmult <- coef(test3)

selfmult.pred <- eulerSIR(c(762,1), c(exp(selfmult[1]), exp(selfmult[2])),0.01, 14)

plot(ts(Its.selfmult))
plot(selfmult.pred[,2])

lq.selfmult <- qnorm(0.025, selfmult.pred[,2], exp(selfmult[3])*sqrt(selfmult.pred[2]))
uq.selfmult <- qnorm(0.975, selfmult.pred[,2], exp(selfmult[3])*sqrt(selfmult.pred[2]))

plot(ts(selfmult.pred[,2],start=0,deltat=0.01), ylim=c(0,400), ylab="infectives")
lines(ts(lq.selfmult,start=0,deltat=0.01), col=2)
lines(ts(uq.selfmult,start=0,deltat=0.01), col=2)


#AIC calculations
K <- 3
aicfunc <- function(K,L){
  return(2*K-2*L)
}
AIC.add <- aicfunc(3,test.add)
AIC.mult <- aicfunc(3,test.mult)
AIC.selfmult <- aicfunc(3,test.selfmult)

AIC.add
AIC.mult
AIC.selfmult

#par(mfrow=c(1,1))
#plot to include in diss
pdf("noisemodsptsyel2.pdf", width = 9, height = 3.2)
par(mfrow=c(1,3))
#add
plot(ts(add.pred[,2],start=0,deltat=0.01), type="l",ylim=range(c(lq.add, uq.add)), ylab="infectives", lwd=0.6, main="Additive")
lines(ts(lq.add,start=0,deltat=0.01), lty=2, lwd=0.6, col=2)
lines(ts(uq.add,start=0,deltat=0.01), lty=2, lwd=0.6, col=2)
lines(ts(Its.add,start=0,deltat=1),type="p", cex=0.7, col=7)
#mult
plot(ts(mult.pred[,2],start=0,deltat=0.01), ylim=c(0,650),ylab="infectives", type="l", main="Multiplicative")
lines(ts(lq.mult,start=0,deltat=0.01), col=2, lty=2)
lines(ts(uq.mult,start=0,deltat=0.01), col=2, lty=2)
lines(ts(Its.mult,start=0,deltat=1),type="p", cex=0.7,col=7)
#selfmult
plot(ts(selfmult.pred[,2],start=0,deltat=0.01), ylim=c(0,400), ylab="infectives", type="l", main="Poisson De-coupled")
lines(ts(lq.selfmult,start=0,deltat=0.01), col=2, lty=2)
lines(ts(uq.selfmult,start=0,deltat=0.01), col=2, lty=2)
lines(ts(Its.selfmult,start=0,deltat=1),type="p",cex=0.7,col=7)
dev.off()
