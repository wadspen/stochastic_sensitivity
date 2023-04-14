

P <- .5 #micro Moles per liter (muM)
p0 <- 5e-07 #Moles per liter (M)
k1 <- 5e+05
k2 <- .2
V <- 10^(-15) #L
na <- 6.023e+23 #Avogadro's constant
na*p0*V
c1 <- 2*k1/(na*V)
c2 <- .2

N=list()
N$M=c(301,0)
N$Pre=rbind(c(2,0),
            c(0,1))
N$Post=rbind(c(0,1),
             c(2,0))
N$h=function(x,th=c(c1,c2)) {
  return(c(th[1]*x[1]*(x[1]-1)/2,
           th[2]*x[2]))
}
N$h(N$M)
S <- N$Post - N$Pre

out=list()

for (i in 1:16) {
  out[[i]]=gillespie(N,400)
}
plot(stepfun(out[[1]]$t,out[[1]]$x[,1]),pch="",xlim=c(0,10))
for (i in 2:16) {
  lines(stepfun(out[[i]]$t,out[[i]]$x[,1]),pch="",xlim=c(0,10))
}


Ps <- matrix(NA,nrow=401,ncol=1000)
for (i in 1:10000) {
  out=gillespie(N,400)
  Ps[,i] <- out$x[,1]
}

dim(Ps)
mean_st <- apply(Ps,MARGIN=1,FUN=mean)
sd_st <- apply(Ps,MARGIN=1,FUN=sd)
upper_st <- mean_st + 3*sd_st
lower_st <- mean_st - 3*sd_st
plot(stepfun(out$t,mean_st),pch='',xlim=c(.4,10),
     ylim=c(0,300))
lines(stepfun(out$t,upper_st),pch='')
lines(stepfun(out$t,lower_st),pch='')

hist(Ps[380,])

op=par(mfrow=c(1,1))
plot(stepfun(out$t,out$x[,1]),pch="",xlim=c(0,10))
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type='l')            





#7.3 Michaelis-Menten
V <- 10^-15 #L
S <- 5e-07 #M
E <- 2e-07 #M  
na*S*V
na*E*V
k1 <- 1e06
k2 <- 1e-04
k3 <- .1

c1 <- k1/(na*V)
c2 <- k2
c3 <- k3

N=list() #(S,E,SE,P)
N$M=c(301,120,0,0)
N$Pre=rbind(c(1,1,0,0),
            c(0,0,1,0),
            c(0,0,1,0))
N$Post=rbind(c(0,0,1,0),
             c(1,1,0,0),
             c(0,1,0,1))
N$h=function(x,th=c(c1,c2,c3)) {
  return( c(th[1]*x[1]*x[2],
            th[2]*x[3],
            th[3]*x[3]))
}
out=gillespie(N,600)
plot(stepfun(out$t,out$x[,1]),pch="",xlim=c(1.9,50))
plot(stepfun(out$t,out$x[,2]),pch="",xlim=c(1.9,50))
plot(stepfun(out$t,out$x[,3]),pch="",xlim=c(1.9,50))
plot(stepfun(out$t,out$x[,4]),pch="",xlim=c(1.9,50))

lv <- function(x,k1,k2,k3) {
  c(k2*x[3] - k1*x[1]*x[2],
    (k2+k3)*x[3] - k1*x[1]*x[2],
    k1*x[1]*x[2] - (k2+k3)*x[3],
    k3*x[3])
}

eul=euler(t=100,fun=lv,ic=c(S,E,0,0),k1=k1,k2=k2,k3=k3)
plot(eul)

#7.4

k1=1;k1r=10;k2=.01;k3=10;k4=.5;k4r=1;k5=.1;k6=.01
N=list() #(g,P2,gP2,r,P)
N$M=c(1,0,0,0,0)
N$Pre=rbind(c(1,1,0,0,0),
            c(0,0,1,0,0),
            c(1,0,0,0,0),
            c(0,0,0,1,0),
            c(0,0,0,0,2),
            c(0,1,0,0,0),
            c(0,0,0,1,0),
            c(0,0,0,0,1))
N$Post=rbind(c(0,0,1,0,0),
             c(1,1,0,0,0),
             c(1,0,0,1,0),
             c(0,0,0,1,1),
             c(0,1,0,0,0),
             c(0,0,0,0,2),
             c(0,0,0,0,0),
             c(0,0,0,0,0))
N$h=function(x,th=c(k1,k1r,k2,k3,k4,k4r,k5,k6)) {
  return( c(th[1]*x[1]*x[2],
            th[2]*x[3],
            th[3]*x[1],
            th[4]*x[4],
            th[5]*x[5]*(x[5]-1)/2,
            th[6]*x[2],
            th[7]*x[4],
            th[8]*x[5]))
}
out=gillespie(N,2000)
plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="") #P2
..plot(stepfun(out$t,out$x[,3]),pch="")
plot(stepfun(out$t,out$x[,4]),pch="") #RNA
plot(stepfun(out$t,out$x[,5]),pch="") #P

Ps=matrix(NA,nrow=2001,ncol=1000)
for (i in 1:1000) {
  out=gillespie(N,2000)
  Ps[,i]=out$x[,2]
}

hist(Ps[10,])
which(out$t>10 & out$t <= 10.1)



Ps <- 100:180
n <- 301
(n - Ps)/2 -> P2s
S
h1 <- function(n,x,c1) {
  c1*(n -2*x)*(n - 2*x -1)/2
}

h2 <- function(x,c2) {
  x*c2
}

(h1(n,P2s-1,c1) + h2(P2s-1,c2))*dpois(P2s-1,(h1(n,P2s-1,c1) + h2(P2s-1,c2)))





