
#DTMC SIS model

time <- 2000
dtt <- .01
beta <- 1*dtt
b <- 0.25*dtt
gamma <- .25*dtt
N <- 100
en <- 50
Ti <- matrix(0,N+1,N+1)
v <- 0:N
p <- matrix(0,time+1,N+1)
pm <- matrix(0,time+1,N+1)
p[1,3] <- 1
bt <- beta*v*(N-v)/N
dt <- (b+gamma)*v
for (i in 2:N) {
  Ti[i,i] <- 1 - bt[i] - dt[i]
  Ti[i,i+1] <- dt[i+1]
  Ti[i+1,i] <- bt[i]
}

Ti[1,1] <- 1
Ti[1,2] <- dt[2]
Ti[N+1,N+1] <- 1-dt[N+1]

for (t in 1:time) {
  y <- Ti%*%p[t,]
  p[t+1,] <- y
}
pm[1,] <- p[1,]
for (t in 1:(time/en)) {
  pm[t+1,] <- p[en*t,]
}

I <- c()
I[1] <- 2


#simulate the trajectory of I
for (i in 2:time) {
  cur <- I[i-1]
  I[i] <- sample(c(cur + 1, cur - 1, cur), 1, 
                 prob = c(beta*cur*(N-cur)/N, (b+gamma)*cur,
                          1 - (beta*cur*(N-cur)/N + (b+gamma)*cur)))
}

plot(I~c(1:time),type='l')


tran <- rbind(c(-1,1),
              c(0,-1),
              c(1,-1),
              c(1,0),
              c(0,0))
SI <- matrix(0, time,2)
SI[1,] <- c(98,2)
for (t in 2:time) {
  s <- SI[t-1,1]
  i <- SI[t-1,2]
  step <- tran[sample(1:5,1,prob=c(dtt*beta*i*s/N,dtt*gamma*i,dtt*b*i,
                                   dtt*b*(N-s-i),
                           1 - dtt*(beta*i*s/N+gamma*i+b*i+b*(N-s-i)))),]
  
  SI[t,] <- c(s,i) + step
}

plot(SI[,2]~c(1:time),type='l')










