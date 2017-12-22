# ESTEC - 18-20 December 2017
# by Rafael S. de Souza
# documentation by Emille E. O. Ishida

# Confidence vs Bayesian intervals 
# Truncated exponential 

# adapted from the Python version presented in
# http://jakevdp.github.io/blog/2014/06/12/frequentism-and-bayesianism-3-confidence-credibility/

#########################################
# Extra packages: do not run!
# require(xkcd)
#########################################

#library(extrafont)
require(R2jags)   # Interface between R and Jags (similars: runjags,rjags,..)
require(mdatools) # For maximum likelihood function
require(ggmcmc)   # To manipulate and plot mcmc chains
require(survival) # For survival analysis, i.e. fit exponential decay
require(dplyr)    # For data-manipulation
source("https://raw.githubusercontent.com/johnbaums/jagstools/master/R/jagsresults.R")

######################################################
###  Truncated exponential function

p_exp <- Vectorize(function(x,theta){
    "Truncated exponential distribution"
  
    if(x > theta){
        return(exp(theta-x))
    } else
    
    return(0)
})
### 


### Generate some data for ploting 
x <- seq(5,18,length.out = 1e3)
g <- data.frame(x=x,y=p_exp(x,10))
###


### Base R plot
plot(g,type="l")

### ggplot style

ggplot(g,aes(x=x,y=y)) +
  geom_line() +
  geom_area(fill="#f4901e") +
 # theme_xkcd() +
  xlab("Failure time") +
  ylab("px") +
 # theme(panel.background = element_rect(color = "black", fill = "gray85") ) 

######################################################

###  Simulate truncated exponential data for fit

y <- c(10,rexp(100) + 10)                # The first 10 ensures to have one single observation 
                                         # equal to theoretical minimum to test the 
                                         # Frequentist and Bayesian approach against a hard coded scenario

nobs <- length(y)                        # number of data points

# Check 
hist(y)                                  # plot data

######## Frequentist solution

# Fit a survival model 
SurvObj   <-  Surv(time=y)
model.par <- survreg(SurvObj~1, dist="exponential")

mu = exp(model.par$coeff)

# Confidence intervals
approx_CI <- function(D,sig=0.95){
    "Calculate confidence intervals"
  
    Nsigma = sqrt(2)*erfinv(sig)
    N = length(D)
    theta_hat <- mean(D) - 1

    return(list(CI_L=theta_hat-Nsigma/sqrt(N), CI_U = theta_hat+Nsigma/sqrt(N)))
}

approx_CI(y)                                      # print result


############ Bayesian solution

# Construct data dictionary
model.data <- list(Y = y,                               # Response variable
                   N = nobs,                            # Sample size
                   minD = min(y)
                   )

Exp <- "model{
    lambda <- 1/theta
    theta ~ dnorm(0,0.001)T(,minD)                      # standard deviation

    # Likelihood function
    for (i in 1:N){
        Y[i] ~ dexp(lambda)
    }
}"

# initial guess
inits <- function() {
  list(theta = runif(1, 0, 10))
}

# parameters
params <- c("theta")

# input data
Surv_fit <- jags(data = model.data,
                 inits = inits,
                 parameters.to.save  = params,
                 model.file = textConnection(Exp),
                 n.chains = 3,
                 n.iter = 5000,
                 n.thin = 1,
                 n.burnin = 1000)


S <- ggs(as.mcmc(Surv_fit)) %>%                     # format data for easy plotting
  filter(.,Parameter == "theta")


###  Plot via base R

hist(S$value)

### Plot chains via ggplot
ggplot(S,aes(x=value)) +
         geom_histogram(fill = "#f4901e",bins= 50,aes(y=..count../sum(..count..))) +
         # theme_xkcd() +
         xlab(expression(theta)) +
         ylab("Frequency")




