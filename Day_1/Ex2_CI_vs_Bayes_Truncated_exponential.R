# Confidence vs Bayesian intervals 
# Truncated exponential 
require(ggplot2)
require(xkcd)
require(R2jags)
require(mdatools)
require(ggmcmc)
source("https://raw.githubusercontent.com/johnbaums/jagstools/master/R/jagsresults.R")

###  Plot truncated exponential
p_exp <- Vectorize(function(x,theta){
if(x > theta){
  return(exp(theta-x))
  } else
return(0)
}
)

x <- seq(5,18,length.out = 1e3)
g <- data.frame(x=x,y=p_exp(x,10))

ggplot(g,aes(x=x,y=y)) +
  geom_line() +
  geom_area(fill="#f4901e") +
  theme_xkcd() +
  xlab("Failiture time") +
  ylab("px")
######



# Simulate truncated exponential 
y <- c(10,rexp(100) + 10) # Maybe there is a better way of doing this. 
#The first 10 ensures to have one single observation equal to theoretical minimum to test the 
#Frequentist and Bayesian approach against a hard coded scenario

nobs <- length(y)


# Check 
hist(y)

#y <- c(10,12,15)

# Fit a survival model 
SurvObj <-  Surv(time=y)
model.par <- survreg(SurvObj~1, dist="exponential")
mu = exp(model.par$coeff)


# Confidence intervals
approx_CI <- function(D,sig=0.95){
Nsigma = sqrt(2)*erfinv(sig)
N = length(D)
theta_hat <- mean(D) - 1
return(list(CI_U = theta_hat+Nsigma/sqrt(N),CI_L=theta_hat-Nsigma/sqrt(N)))
}

approx_CI(y)



#CI <- exp(confint(model.par))
#print(CI)

# Bayesian solution

# Construct data dictionary
model.data <- list(Y = y,                               # Response variable
                   N = nobs,                           # Sample size
                   minD = min(y)
)

Exp <- "model{
lambda <- 1/theta
theta ~ dnorm(0,0.001)T(,minD)                                        # standard deviation
# Likelihood function

for (i in 1:N){
Y[i] ~ dexp(lambda)
}

}"

inits <- function() {
  list(theta = runif(1, 0, 10))
}
params <- c("theta")

Surv_fit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = params,
                model.file = textConnection(Exp),
                n.chains = 3,
                n.iter = 5000,
                n.thin = 1,
                n.burnin = 1000)


S <- ggs(as.mcmc(Surv_fit)) %>%
  filter(.,Parameter == "theta")

ggplot(S,aes(x=value)) +
         geom_histogram(fill = "#f4901e",bins= 50,aes(y=..count../sum(..count..))) +
         theme_xkcd() +
  xlab(expression(theta)) +
  ylab("Frequency")




