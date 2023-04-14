#stochastic process practice
rfmc <- function(n,P,pi0)
{
  v=vector("numeric",n)
  r=length(pi0)
  v[1]=sample(r,1,prob=pi0)
  for (i in 2:n) {
    v[i]=sample(r,1,prob=P[v[i-1],])
  }
  ts(v)
}

P=matrix(c(0.9,0.1,0.2,0.8),ncol=2,byrow=TRUE)
P
pi0=c(0.5,0.5)
pi0
samplepath=rfmc(200,P,pi0)
samplepath
plot(samplepath)
hist(samplepath)
e=eigen(t(P))$vectors[,1]
e/sum(e)
f=eigen(t(P))







rdiff <- function(afun,bfun,x0=0,t=50,dt=0.01, ...) 
{
  n <- t/dt
  xvec <- vector('numeric',n)
  sdt <- sqrt(dt)
  for (i in 1:n) {
    t <- i*dt
    x <- x + afun(x,...)*dt +
      bfun(x,...)*rnorm(1,0,sdt)
    xvec[i] <- x
  }
  ts(xvec,deltat=dt)
}


afun <- function(x,lambda=1,mu=0.1) {
  return(lambda - mu*x)
}

bfun <- function(x,lambda=1,mu=0.1) {
  return(sqrt(lambda + mu*x))
}


state <- rdiff(afun,bfun,lambda=1,mu=0.1)

plot(state~time)
plot(state)


euler <- function(t=50,dt=0.001,fun=f,ic=c(1,1),...)
{
  p=length(ic)
  n <- t/dt
  xmat=matrix(0,ncol=p,nrow=n)
  x=ic
  xmat[1,]=x
  for (i in 2:n) {
    x = x + fun(x,...)*dt
    xmat[i,]=x
  }
  ts(xmat,start=0,deltat=dt)
}

lv <- function(x,k1=1,k2=0.1,k3=0.1)
{
  c(k1*x[1]-k2*x[1]*x[2],
    k2*x[1]*x[2] - k3*x[2])
}

lv <- function(x,k1,k2,k3) {
  c(k2*x[3] - k1*x[1]*x[2],
    (k2+k3)*x[3] - k1*x[1]*x[2],
    k1*x[1]*x[2] - (k2+k3)*x[3],
    k3*x[3])
}

plot(euler(t=100,fun=lv,ic=c(4,10)))


gillespie <- function(N,n,...) {
  tt=0
  x=N$M
  S=t(N$Post-N$Pre)
  u=nrow(S)
  v=ncol(S)
  tvec=vector("numeric",n)
  xmat=matrix(0,ncol=u,nrow=n+1)
  xmat[1,]=x
  for (i in 1:n) {
    h=N$h(x,...)
    tt=tt+rexp(1,sum(h))
    j=sample(v,1,prob=h)
    x=x+S[,j]
    tvec[i]=tt
    xmat[i+1,]=x
  }
  return(list(t=tvec,x=xmat))
}

N=list()
N$M=c(50,100)
N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
N$h=function(y,th=c(1,0.005,0.6))
{return(c(th[1]*y[1],th[2]*y[1]*y[2],th[3]*y[2] )) }

out=gillespie(N,10000)

op=par(mfrow=c(2,2))
plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type="l")
par(op)

discretise <- function(out,dt=1,start=0)
{
  events=length(out$t)
  end=out$t[events]
  len=(end-start)%/%dt+1
  x=matrix(0,nrow=len,ncol=ncol(out$x))
  target=0
  j=1
  for (i in 1:events) {
    while (out$t[i] >=target) {
      x[j,]=out$x[i,]
      j=j+1
      target=target+dt
    }
  }
  ts(x,start=0,deltat=dt)
}

plot(discretise(out,dt=0.01))

discretise(out,dt=0.01) -> dude



X1 <- rbind(c(1,-1,0,0,-1),
      c(-1,1,0,0,1),
      c(0,0,1,0,0),
      c(0,0,0,1,0),
      c(0,0,0,-2,1),
      c(0,0,0,2,-1),
      c(0,0,-1,0,0),
      c(0,0,0,-1,0))

rankMatrix(X1)


X <- rbind(c(1,-1,0,0,-1),
      c(0,0,1,0,0),
      c(0,0,0,1,0),
      c(0,0,0,2,-1))

y <- c(0,0,0,0)
solve(X,y)
library(matlib)
gaussianElimination(X,y,verbose=TRUE)



A <- matrix(c(2, 1, -1,
              -3, -1, 2,
              -2,  1, 2), 3, 3, byrow=TRUE)
b <- c(8, -11, -3)
gaussianElimination(A, b)
gaussianElimination(A, b, verbose=TRUE, fractions=TRUE)
gaussianElimination(A, b, verbose=TRUE, fractions=TRUE, latex=TRUE)

# determine whether matrix is solvable
gaussianElimination(A, numeric(3))

# find inverse matrix by elimination: A = I -> A^-1 A = A^-1 I -> I = A^-1
gaussianElimination(A, diag(3))
inv(A)

# works for 1-row systems (issue # 30)
A2 <- matrix(c(1, 1), nrow=1)
b2 = 2
gaussianElimination(A2, b2)
showEqn(A2, b2)
# plotEqn works for this case
plotEqn(A2, b2)








