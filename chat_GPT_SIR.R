# Set the initial conditions and parameters
S <- 990
I <- 10
R <- 0
beta <- 0.3
gamma <- 0.1
tmax <- 100
t <- 0

# Create empty vectors to store the results
Svec <- numeric()
Ivec <- numeric()
Rvec <- numeric()
tvec <- numeric()

# Set the initial values
Svec[1] <- S
Ivec[1] <- I
Rvec[1] <- R
tvec[1] <- t

# Define the Gillespie algorithm function
gillespie <- function(S, I, R, beta, gamma) {
  # Calculate the total rate of events
  total_rate <- beta*S*I + gamma*I
  
  # Generate two random numbers
  r1 <- runif(1)
  r2 <- runif(1)
  
  # Calculate the time to the next event
  tau <- -log(r1)/total_rate
  
  # Determine which event occurs
  if (r2 < beta*S*I/total_rate) {
    # An infection event occurs
    S <- S - 1
    I <- I + 1
  } else {
    # A recovery event occurs
    I <- I - 1
    R <- R + 1
  }
  
  # Return the updated compartment values and time
  return(list(S=S, I=I, R=R, tau=tau))
}

# Run the Gillespie algorithm
while (t < tmax) {
  # Call the Gillespie algorithm function
  res <- gillespie(S, I, R, beta, gamma)
  
  # Update the values
  S <- res$S
  I <- res$I
  R <- res$R
  t <- t + res$tau
  
  # Store the values
  Svec <- c(Svec, S)
  Ivec <- c(Ivec, I)
  Rvec <- c(Rvec, R)
  tvec <- c(tvec, t)
}

# Plot the results
plot(tvec, Svec, type='l', xlab='Time', ylab='Population', ylim=c(0, 1000), col='blue', lwd=2)
lines(tvec, Ivec, col='red', lwd=2)
lines(tvec, Rvec, col='green', lwd=2)
legend('topright', c('Susceptible', 'Infected', 'Recovered'), col=c('blue', 'red', 'green'), lwd=2)
