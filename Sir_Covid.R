library(ggplot2)

covid_simulation <- function(R0, dt, D, days) {
  
  # Here a simple covid model is run, where the true parameters are set in advance. In particular, we consider:
  # N - the number of the entire population
  # R0 - basic reproduction number
  # dt - time step used to approximate a differential equation
  # D - duration if the infectious period
  # days - number of days for which we run the simulation
  
  # The model only considers a closed population, hence no births or deaths are taken into account. 
  # This is SIR model, hence if compared to data, we expect an overestimation of the basic reproduction number
  
  
  # set the size of the population
  N <- 64 * 10^6
  
  # initial number of susceptible and infected individuals
  S0 <- N
  I0 <- 1
  
  # calculate the recovery rate
  r <- dt/D
  
  # transition parameter
  beta <- R0*r/N
  
  # total time steps
  steps <- days/dt
  
  # set vectors to store the number of the respective cases
  S <- S0
  I <- I0
  R <- 0
  
  # main loop begins here
  for (i in 1:steps) {
    # calculate the new number of susceptible individuals
    S_Next <- S[i] - beta*S[i]*I[i]
    # calculate the number of new infectious individuals
    I_Next <- I[i] + beta*S[i]*I[i] - r*I[i]
    # calculate the number of new recovered individuals
    R_Next <- R[i] + r*I[i]
    
    # store the results in a vector
    S <- c(S, S_Next)
    I <- c(I, I_Next)
    R <- c(R, R_Next)
  }
  
  return(list("time" = 1:(steps+1), "susceptible" = S, "infectious" = I, "recovered" = R))
}

# run the simulation
output <- covid_simulation(R0 = 5, dt = 0.1, D = 14, days = 150)

# store the output in a data frame
results <- data.frame(T = output$time, S = output$susceptible, I = output$infectious, R = output$recovered)

# plot the results for all three compartments
ggplot(data = results, aes(x=T)) +
  geom_line(aes(y=S, colour = "Susceptible")) +
  geom_line(aes(y=I, colour = "Infectious")) +
  geom_line(aes(y=R, colour = "Recovered")) +
  coord_cartesian(ylim = c(0, 6.4*10^7)) +
  scale_colour_manual("",
                      breaks = c("Susceptible", "Infectious", "Recovered"),
                      values = c("blue", "red", "green")) +
  ggtitle("SIR model for Covid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Number of people")

# generate noise around the true observations
noise <- NULL
days <- seq(from = 10, to = 1500, by = 10)

for (i in 1:150) {
  l <- results$I[days[i]]
  noise <- c(noise, rpois(1, l))
}

noise_data <- data.frame(d = days, n = noise)

# plot the results for the infectious individuals
ggplot(data = results, aes(x=T)) +
  geom_line(aes(y=I), colour = "red") +
  geom_point(data = noise_data, aes(x=d, y=n), colour = "blue", shape = 4, size = 2) +
  labs(x = "Time", y = "Number of infectious individuals")

# plot the data only
ggplot(data = results, aes(x=T)) +
  geom_point(data = noise_data, aes(x=d, y=n), colour = "blue", shape = 4, size = 2) +
  labs(x = "Time", y = "Number of infectious individuals")

################################################################################################################
################################################################################################################

# construct the time dependent poisson likelihood
poissonLL <- function(theta, x) {
  # set the size of the population
  N <- 64 * 10^6
  
  beta <- theta[1]
  r <- theta[2]
  
  output <- covid_simulation(R0 = beta*N/r, dt = 0.1, D = 0.1/r, days = 150)
  lambda <- output$infectious[seq(0, length(output$infectious), 10)]

  return(sum(x*log(lambda)) - sum(lambda))
}

# use optim to find the MLEs
opt <- optim(poissonLL, par = c(5.3*10^(-10), 0.007), 
             lower=c(5.2*10^(-10),0.005),
             upper=c(5.9*10^(-10),0.008),
             method="L-BFGS-B",
             x = noise, control=list(fnscale=-1))

# calculate the true input parameters
true_r <- 0.1/14
true_beta <- 5*true_r/(64*10^6)

# print the true values of the parameters
true_r
true_beta

opt