library("tidyverse")
library("ggplot2")
library("maxLik")
#x0=(S0,I0)' theta=(beta, gamma)' T=end time

#FINAL product
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

#retry - WORKING CODE

#dargatz <- eulerSIR(c(762,1),c(exp(-6),0.5),0.01,14)

Its <- rnorm(15, dargatz[1+(0:14)*100, 2],1)

plot(ts(dargatz[,2], start=0, deltat=0.01))
lines(ts(Its,start=0,deltat=1),type="p")

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

#working code data input attempts

loglike_epi(c(-6,log(0.5), 0), Its,0.01,c(762,1))

library(maxLik)

test <- maxLik(loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its,deltat=0.01,x0=c(762,1))
coef(test)
#end of code
plot(ts(dargatz, start=0, deltat=0.01), main=NULL)



#for the diss without title
png("boardingschool2.png")
plot(ts(dargatz, start=0, deltat=0.01), main=NULL)
dev.off()





































eulerSIR <- function(x0, theta, deltat, T){
  s=floor(T/deltat)  #number of steps
  iterations <- matrix(0,s+1,2)
  #put in initial values
  iterations[1,] <- c(x0[1], x0[2])
  #define f(S,I)
  t=deltat
  S <- x0[1]
  I <- x0[2]
  for (i in 1:s) {
    #calculate 
    dS <- -theta[1]*S*I
    dI <- theta[1]*S*I-theta[2]*I
    iterations[i+1,1] <- iterations[i,1]+ deltat*dS
    iterations[i+1,2] <- iterations[i,2]+deltat*dI
    S <- iterations[i+1,1]
    I <- iterations[i+1,2]
  }
  return(iterations)
}

dargatz <- eulerSIR(c(762,1),c(exp(-6),0.5),0.01,14)

plot(ts(dargatz, start=0, deltat=0.01))

#dnorm(Yt[1], mean=It[1], sd=sigma)

#Assume the Yt are normally distributed with mean It and variance sigma^2
#below theta=c(beta, gamma, sigma) and y=c(yt1,....,ytn)

loglike_epi <- function(theta, y, N){ #how do I make the function indpt of N?
  len <- length(y)
  It <- eulerSIR(c(N-y[1],y[1]), c(theta[1], theta[2]), 1, len)[,2]
  #likeis <- c(1:len)
  likeis <- numeric(len)
  for (i in 1:len){
    likeis[i] <- dnorm(y[i], mean=It[i], sd=theta[3], log=TRUE)
  }
  #loglikeis <- log(likeis)
  loglike <- sum(likeis)
  return (loglike)
} 

#for boarding school data
beta <- exp(-6)
gamma <- 0.5
N <- 763
Yt <- c(1,3,6,25,73,221,294,257,236,189,125,67,26,10,3)
sigma <- sd(Yt) # is sd from the yt or the It?

loglike_epi(c(beta, gamma, sigma), Yt, N)

#for ebola data - GUINEA specifically
load(file="ebola.rda")
b <- 7
g <- 1.8/(2000*7)
n <- 2000
Yt_e <- as.numeric(c(ebola[,2]))
s <- sd(Yt_e)

loglike_epi(c(b,g,s), Yt_e, n)

#maxlik()



#KEY CODE

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
    S <- iterations[i-1,1]
    I <- iterations[i-1,2]
  }
  return(iterations)
}
dargatz <- eulerSIR(c(762,1),c(exp(-6),0.5),0.01,14)

plot(ts(dargatz, start=0, deltat=0.01))

#retry - WORKING CODE
dargatz <- eulerSIR(c(762,1),c(exp(-6),0.5),0.01,14)
Its <- rnorm(15, dargatz[1+(0:14)*100, 2],1)
loglike_epi <- function(theta, y, deltat, x0){
  N <- length(y)
  inc <- 1/deltat
  output <- eulerSIR(x0, exp(theta), deltat, (N-1))
  index <- 1+(0:(N-1))*inc
  output <- output[index, 2]
  loglikeindiv <- dnorm(y, output, exp(theta[3]), log=TRUE)
  loglike <- sum(loglikeindiv)
  return (loglike)
}

#working code data input attempts
loglike_epi(c(-6,log(0.5), 0), Its,0.01,c(762,1))

#mle calcs
library(maxLik)

test <- maxLik(loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its,deltat=0.01,x0=c(762,1))
coef(test)


