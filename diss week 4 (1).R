# week 4 diss work
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

#additive loglike fn below
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

#start this week's new code
S <- matrix(data=c(-1,1,0, -1), 2,2)

SIR_haz <- function(euler, theta){
  haz <- numeric(length(theta))
  haz[1] <- theta[1]*euler[1]*euler[2]
  haz[2] <- theta[2]*euler[2]
  return (haz)}

# sim_pl <- function(x0, theta, deltat, endt){
#   sims <- endt/deltat
#   ret <- matrix(0,2,sims)
#   ret[1,] <- x0
#   for(i in 2:sims){
#     ret[i,] <- SIR_haz(ret[i-1,], theta)
#   }
# }


sim_pl <- function(x0, theta, deltat, endt){
  #set up
  n <- endt/deltat
  ret <- matrix(0,n,2)
  ret[1,1] <- x0[1]
  ret[1,2] <- x0[2]
  #hazards
  haz <- SIR_haz(x0, theta)
  for(i in 2:n){
    R1.sim <- rpois(1,haz[1]*deltat)
    R2.sim <- rpois(1,haz[2]*deltat)
    R <- c(R1.sim, R2.sim)
    ret[i,] <- ret[i-1,] + (S%*%R)
    haz <- SIR_haz(ret[i,], theta)
  }
  return(ret)
}

#works here:
set.seed(132)
sim.bs <- sim_pl(c(762,1), c(exp(-6),0.5), 0.1, 14)
plot(ts(sim.bs))

#and here:
sim.bs2 <- sim_pl(c(762,1), c(exp(-6),0.5), 0.1, 14)
plot(ts(sim.bs2))

#task 2
wrap.sims <- function(x0, theta, deltat, endt, sims){
  n <- endt/deltat
  S.simd <- matrix(0,n,sims )
  I.simd <- matrix(0,n,sims)
  for(i in 1:sims){
   S.simd[,i] <- sim_pl(x0, theta, deltat, endt)[,1]
   I.simd[,i] <- sim_pl(x0, theta, deltat, endt)[,2]
  }
  return(list(S.simd, I.simd))
}

sim.bsd <- wrap.sims(c(762,1), c(exp(-6),0.5), 0.1, 14, 5)

simd.S <- sim.bsd[[1]]
simd.I <- sim.bsd[[2]]

sims.Smeans <- rowMeans(simd.S)
sims.Imeans <-rowMeans(simd.I)

# dargatz <- eulerSIR(c(762,1),c(exp(-6),0.5),0.01,14)
# Its <- rnorm(15, dargatz[1+(0:14)*100, 2],1)
# test1 <- maxLik(loglike_epi,start=c(-5,1,1),tol=-1,reltol=1e-12,y=Its,deltat=0.01,x0=c(762,1))
# coeff <- coef(test1)
# 
# pred <- eulerSIR(c(762,1), c(exp(coeff[1]), coeff[2]),0.1, 14)

q.loop <- function(simd,prob){
  sim.mean <- numeric(nrow(simd))
  for(i in 1:nrow(simd)){
    sim.mean[i] <- quantile(simd[i,], probs=prob)
    print(sim.mean)
  }
  print(sim.mean)
  return (sim.mean)
}

#apply(simd, 2, quantile, 0.025/0.975)

uq.Smeans <- q.loop(simd.S, 0.975)
lq.Smeans <- q.loop(simd.S, 0.025)
lq.Imeans <- q.loop(simd.I, 0.025)
uq.Imeans <- q.loop(simd.I, 0.975)

# uq.Smeans <- apply(simd.S, 2, quantile, 0.975)
# lq.Smeans <- apply(simd.S, 2, quantile, 0.025)
# uq.Imeans <- apply(simd.I, 2, quantile, 0.975)
# lq.Imeans <- apply(simd.I, 2, quantile, 0.025)

par(mfrow=c(1,2))
plot(ts(sims.Smeans), ylab="Susceptibles")
lines(ts(lq.Smeans), col=2)
lines(ts(uq.Smeans), col=2)

plot(ts(sims.Imeans), ylab="Infectives", ylim=c(0,300))
lines(ts(lq.Imeans), col=2)
lines(ts(uq.Imeans), col=2)


#for incl in diss
pdf("poisleap01.pdf", width = 7, height = 3)
par(mfrow=c(1,2))
plot(ts(sims.Smeans), ylab="Susceptibles", cex.lab=0.2)
lines(ts(lq.Smeans), col=2)
lines(ts(uq.Smeans), col=2)

plot(ts(sims.Imeans), ylab="Infectives", ylim=c(0,350), cex.lab=0.2)
lines(ts(lq.Imeans), col=2)
lines(ts(uq.Imeans), col=2)
dev.off()

