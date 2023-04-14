
errors <- data.frame()
m <- 1
for (i in 1:400) {
  tryCatch({
    if (i %% 16 == 0) {9 <- 14}
    # else if (i %% 17 == 0) {print(i)}
  }, error =function(e){err <<- data.frame(i=i); errors <<- rbind(errors,err)})
}

tryCatch(stop(e),finally=print('Hello'))
nile <- as.vector(Nile)

sigeps <- sqrt(15099)
sigeta <- sqrt(1469.1)
N <- 10000
a1 <- nile[1]
att <- a1
w_t_1 <- rep(1/N,N)

alphas <- c()
Ps <- c()
alphast <- rnorm(N,att,sigeta)
ESS <- c()
ESS[1] <- sum(w_t_1^2)^(-1)
for (i in 2:length(nile)) {
  alpha_tild <- apply(as.matrix(alphast),MARGIN=1,FUN=function(x) 
    {rnorm(1,x,sigeta)})
  alpha_tild <- as.vector(alpha_tild)
  wtild_t <- w_t_1*exp(-(1/2)*log(2*pi) - (1/2)*log(sigeps^2) - 
                         (1/2)*(sigeps^(-2))*(nile[i]-alpha_tild)^2)
  
  w_t <- wtild_t/sum(wtild_t)
  att_hat <- sum(w_t*alpha_tild)
  Ptt_hat <- sum(w_t*(alpha_tild - att_hat)^2)
  alphas[i] <- att_hat
  Ps[i] <- Ptt_hat
  w_t_1 <- w_t
  # alphast <- alpha_tild
  # c1 <- runif(1,0,1/N)
  # Walpha <- c1 + (1:N -1)/N
  samp <- sample(1:length(alpha_tild),prob=wtild_t,replace=TRUE)
  alphast <- alpha_tild[samp]
  w_t <- w_t[samp]
  w_t_1 <- w_t/sum(w_t)
  ESS[i] <- sum(w_t_1^2)^(-1)
}

vars <- Ps*sigeps^2/(Ps + sigeps^2)

upper <- alphas + qnorm(.95)*sqrt(Ps)
lower <- alphas - qnorm(.95)*sqrt(Ps)

ggplot() +
  geom_point(aes(x=1871:1970,y=nile)) +
  geom_line(aes(x=1871:1970,y=alphas),colour='forestgreen') +
  geom_line(aes(x=1871:1970,y=upper),colour='forestgreen',linetype='dashed') +
  geom_line(aes(x=1871:1970,y=lower),colour='forestgreen',linetype='dashed') +
  theme_bw()

ggplot() +
  geom_line(aes(x=1871:1970,y=Ps))

ggplot() +
  geom_line(aes(x=1871:1970,y=ESS))


