library("tidyverse")
library("ggplot2")
#x0=(S0,I0)' theta=(beta, gamma)' T=end time

View(dalziel)
data(dalziel, package="epimdr")

eulerSIR=function(x0, theta, deltat, T){
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


#trying to make it work for dalziel specifically

dalziel_baltimore <- dalziel |> 
  filter(loc=="BALTIMORE") |> 
  arrange(,biweek) |> 
  arrange(,year) |> 
  mutate(,"sus"= pop+rec) |> 
  na.omit()

S.int <- dalziel_baltimore$sus[1]
I.int <- dalziel_baltimore$cases[1]

dalziel_balt_plot <- eulerSIR(c(S.int, I.int),
                              c(0.156, 1.32),
                              (1948-1906)/(7+26*(1948-1908)+13), 1948-1906)

plot(ts(dalziel_balt_plot)) #this is very wrong

#actual data
plot(ts(dalziel_baltimore$sus, dalziel_baltimore$cases))

#######Golightly Brief Code

odesim<-function(param,dt,x0,endT)
{
  n <- endT/dt
  x <- matrix(0,nrow=n+1,ncol=2)
  x[1,] <- x0 
  for (i in 2:(n+1)) 
  {
    #right hand side of ODE system
    RHS1 <- -param[1]*x[i-1,1]*x[i-1,2]
    RHS2 <- param[1]*x[i-1,1]*x[i-1,2]-param[2]*x[i-1,2]
    #time step ODE system
    x[i,] <- x[i-1,] + c(RHS1,RHS2)*dt
  }
  return(x)
}

#Example numerical solution (St and It against time t)
#Base graphics but can use ggplot2 if you wish
out <- odesim(c(exp(-6),0.5),0.01,c(762,1),14)
plot(ts(out,start=0,deltat=0.01))

###End

