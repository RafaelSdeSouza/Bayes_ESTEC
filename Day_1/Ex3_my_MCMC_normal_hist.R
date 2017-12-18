# ESTEC - 18-20 December 2017
# by Rafael S. de Souza
# 
# (c) 2017, Rafael S. de Souza 
# 
# Customized MCMC sampler
# Case: Normal distribution

library(fitdistrplus)
library(ggplot2)
library(ggthemes)
library(stats4)
require("plot3D")
# Generate data
set.seed(1056)                          # set seed to replicate example
nobs = 200                              # number of obs in model
y <- rnorm(nobs,5,1.5)                  # random uniform variable

gd <- data.frame(x=y)

hist(y)


ggplot(data = gd,aes(x = y)) +
  geom_histogram(binwidth=1,fill="#FF5733",color="black") + 
  scale_fill_tableau()  + 
  #theme_xkcd() +
  geom_rug(sides = "b", aes(y = 0), position = "jitter", colour = "cyan2")+
  theme(panel.background = element_rect(color = "black", fill = "gray85") )




## Fit a normal distribution to the 50 random data set
fitdist(y,"norm") 

# Likelihood function
likelihood <- function(param){
  "Normal likelihood"
  
  mu = param[1]
  sd = param[2]
  LL = dnorm(y, mean = mu, sd = sd)
  return(sum(log(LL)))
}

# 3D likelihood
grid.lines = 26                # predict values on regular xy grid
x1.pred <- seq(1,10, length.out = grid.lines)
x2.pred <- seq(0.5, 10, length.out = grid.lines)
x1x2 <- expand.grid( x1 = x1.pred, x2 = x2.pred)

y.pred <- matrix(apply(x1x2,1,likelihood), 
                 nrow = grid.lines, ncol = grid.lines)

# scatter plot with regression plane
contour3D(z = y.pred,                   
          cex = 0.5, cex.lab=1.5, colvar = z, 
          theta = 130, phi = 25, ticktype = "detailed",
          col="red2",bty = "b2",t="l",
          xlab="x1",
          ylab="x2",
          zlab="y")

slopevalues <- function(y){return(likelihood(c(y,1.5)))}
plot(seq(0,10,length.out = 100),sapply(seq(0,10,length.out = 100),slopevalues),type="l")
abline(v = 5)


x1 <- seq(-1,1,by=0.1)
plot(x1,dnorm(x1))


# Fit via ML
LL <- function(mu,sd){
  -likelihood(c(mu,sd))
}

fit <- mle(LL, start = list( mu = 3, sd = 1), method = "L-BFGS-B", lower = c(-Inf, 0))
summary(fit)

# Prior distribution
prior <- function(param){
  mu  = param[1]
  sd = param[2]
  muprior = dnorm(mu, mean = 0,sd = 10, log = T)
  sdprior = dunif(sd, min = 0.001, max = 10, log = T)
  return(muprior +  sdprior)
}

# Posterior
p_target <- function(param){
  return(likelihood(param) + prior(param))
}


# MCMC starts now

proposalfunction <- function(param,alpha=0.75){
  return(c(rnorm(1,mean = param[1],sd = 0.25),max(1e-3,runif(1,param[2] - alpha,param[2]+alpha))))
}

startvalue <- c(1,2)

# Define Bernoulli sampler as special case of binomial distribution
rbern <- function(p){
  rbinom(1,1,p)
  }

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations + 1,2))
  chain[1,] = startvalue
  
  for (i in 1:iterations){
    s_prop = proposalfunction(chain[i,])
    
    #Hasting ratio: exponential term accounts for loglikelihood and priors 
    H = exp(p_target(s_prop) - p_target(chain[i,]))
    
    # acceptance probability
    p_trans = min(1,H)
    Jump = rbern(p_trans)
    
    if (Jump == 1){
      chain[i + 1,] = s_prop
    }else{
      chain[i + 1,] = chain[i,]
    }
  }
  return(chain)
}
N_it <- 1000
burnIn <- 200
chain = run_metropolis_MCMC(startvalue, N_it)

chain2 <- data.frame(chain[,1],chain[,2])
chain2$burn <- as.factor(c(rep("burnin",burnIn + 1),rep("chain",N_it-burnIn)))
colnames(chain2) <- c("x","y","iteration")

# Plot
plot(chain2[,1],chain2[,2],type ="l")


ggplot(chain2[1:1000,],aes(x = x,fill = iteration)) +
  #  geom_point(size=0.2) +
  geom_bar(width=0.025) + scale_fill_tableau()  + 
  #  geom_density_2d(size=0.1) +
  xlab(expression(y)) + ylab("Counts") +
  # theme_xkcd() +
  theme(panel.background = element_rect(color = "black", fill = "gray85") )


d0 <- data.frame(x=5,y=1.5)
 ggplot(chain2[1:1000,],aes(x=x,y=y,colour=iteration)) +
 geom_path() + scale_color_tableau()  + 
  geom_point(data=d0,mapping=aes(x=x,y=y),colour="red") +
  #  geom_density_2d(size=0.1) +
  xlab(expression(mu)) + ylab(expression(sigma)) +
  # theme_xkcd() +
  theme(panel.background = element_rect(color = "black", fill = "gray85") ) 

