library("tidyverse")
library("ggplot2")
library("maxLik")

#function carried over
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

Its <- rnorm(15, dargatz[1+(0:14)*100, 2],1)

plot(ts(dargatz[,2], start=0, deltat=0.01))
lines(ts(Its,start=0,deltat=1),type="p")

#function that calculates the log-likelihood
loglike_epi <- function(theta, y, deltat, x0){
  N <- length(y)
  inc <- 1/deltat
  output <- eulerSIR(x0, exp(c(theta[1],theta[2])), deltat, N-1)
  index <- 1+(0:(N-1))*inc
  output <- output[index, 2]
  loglikeindiv <- dnorm(y, output, exp(theta[3]), log=TRUE)
  loglike <- sum(loglikeindiv)
  return (loglike)
}

loglike_epi(c(-6,log(0.5), 0), Its,0.01,c(762,1))

library(maxLik)

test <- maxLik(loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its,deltat=0.01,x0=c(762,1))
coef(test)

plot(ts(dargatz, start=0, deltat=0.01), main=NULL)

#for the diss without title
png("boardingschool2.png")
plot(ts(dargatz, start=0, deltat=0.01), main=NULL)
dev.off()

