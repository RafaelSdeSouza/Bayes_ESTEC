#
# (c) 2017, Rafael S. de Souza 
# 
# Customized MCMC sampler
# Case: linear regression

# Generate data
set.seed(1056)                      # set seed to replicate example
nobs = 100                          # number of obs in model
x1 <- rnorm(nobs,0,5)               # random uniform variable
alpha = 1.5                         # intercept
beta = 4                            # angular coefficient
xb <- alpha + beta*x1              # linear predictor, xb
sd <- 1                             # Standard deviation
y <- rnorm(nobs, xb, sd = sd)         # create y as  random normal variate

# Fit with R function lm
summary(mod <- lm(y ~ x1))          # model of the synthetic data.

# Plot Output
ypred <- predict(mod,type="response")                # prediction from the model
plot(x1,y,pch=19,col="red")                          # plot scatter
lines(x1,ypred,col='grey40',lwd=2)                   # plot regression line
segments(x1,fitted(mod),x1,y,lwd=1,col="gray70")     # add the residuals




# MCMC solution

# Likelihood function
likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  xb <- a + b*x1
  LL = dnorm(y, mean = xb, sd = sd)
  return(sum(log(LL)))
}

# Fit via Maximum Likelihood

library(stats4)
LL <- function(a,b,sd){
  -likelihood(c(a,b,sd))
}


# Fit
fit <- mle(LL, start = list( a = 2, b = 2, sd = 1), method = "L-BFGS-B", lower = c(-Inf, -Inf, 0))
summary(fit)


# Prior distribution
prior <- function(param){
  a  = param[1]
  b =  param[2]
  sd = param[3]
  aprior = dnorm(a, mean = 0,sd = 10, log = T)
  bprior = dnorm(b, mean = 0,sd = 10, log = T)
  sdprior = dunif(sd, min = 0.001, max = 30, log = T)
  return(aprior + bprior + sdprior)
}


posterior <- function(param){
return(likelihood(param) + prior( param))
}

proposalfunction <- function(param,alpha=1){
  return(c(rnorm(2,mean = param[1:2],sd = 0.25),max(1e-3,runif(1,param[3]-alpha,param[3]+alpha))))
}


startvalue <- c(1,2,1)

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    
    # acceptance probability
    p_accept = min(1,probab)
    accept = rbinom(1,1,p_accept)
    
    
    if (accept == 1){
      chain[i + 1,] = proposal
    }else{
      chain[i + 1,] = chain[i,]
    }
  }
  return(chain)
}
N_it <- 5e4
burnIn <- 1e4
chain = run_metropolis_MCMC(startvalue, N_it)

chain2 <- data.frame(chain[,1],chain[,2],chain[,3])
chain2$burn <- as.factor(c(rep("burnin",burnIn + 1),rep("chain",N_it-burnIn)))
colnames(chain2) <- c("x","y","z","iteration")

plot(chain2[,1],chain2[,2],type ="l")


  d0 <- data.frame(a=1.5,b=4)
  ggplot(chain2,aes(x=x,y=y,colour=iteration)) +
  geom_path() + scale_color_tableau()  + 
  xlab(expression(alpha)) + ylab(expression(beta)) +
 geom_point(data=d0,mapping=aes(x=a,y=b),colour="red") +
  theme_xkcd() +
  theme(panel.background = element_rect(color = "black", fill = "gray85") ) +
  coord_cartesian(xlim=c(1,2.25),ylim=c(1.5,4.5))




par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],nclass=30,  main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = alpha , col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = beta, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main ="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = sd, col="red" )  
plot(chain[-(1:burnIn),1], type = "l", xlab = "True value = red line" , main = "Chain values of a", )
abline(h = alpha, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab = "True value = red line" , main = "Chain values of b", )
abline(h = beta, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab = "True value = red line" , main = "Chain values of sd", )
abline(h = sd, col="red" )








# Fit using MCMCpack
library(MCMCpack)
posteriors <- MCMCregress(y ~ x1, thin=1, seed=1056, burnin=1000,
                          mcmc=10000, verbose=1)
# Output
summary(posteriors)
plot(posteriors)


