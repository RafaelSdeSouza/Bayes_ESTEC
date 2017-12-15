# period, luminosity, and color in early-type contact binaries

# my_MCMC

PLC <- read.csv("https://raw.githubusercontent.com/astrobayes/BMAD/master/data/Section_10p3/PLC.csv", header = T)

#PLC <- filter(PLC, type == "genuinely-contact ")


nobs = nrow(PLC)                        # number of data points
x1 <- PLC$logP                          # log period
x2 <- PLC$V_I                           # V-I color
y <- PLC$M_V                            # V magnitude


# Likelihood function
likelihood <- function(param){
  b0 = param[1]
  b1 = param[2]
  b2 = param[3]
  sd = param[4]
  xb <- b0 + b1*x1 + b2*x2
  LL = dnorm(y, mean = xb, sd = sd)
  return(sum(log(LL)))
}

# Fit via Maximum Likelihood

library(stats4)
LL <- function(b0,b1,b2,sd){
  -likelihood(c(b0,b1,b2,sd))
}

# Fit
fit <- mle(LL, start = list( b0 = 0, b1 = 0, b2 = 0, sd = 0.5), method = "L-BFGS-B", lower = c(-Inf, -Inf,-Inf, 0))
summary(fit)

# MCMC solution


# Prior distribution
prior <- function(param){
  b0  = param[1]
  b1 =  param[2]
  b2 =  param[3]
  sd = param[4]
  b0prior = dnorm(b0, mean = 0,sd = 2, log = T)
  b1prior = dnorm(b1, mean = 0,sd = 2, log = T)
  b2prior = dnorm(b2, mean = 0,sd = 2, log = T)
  sdprior = dunif(sd, min = 0.001, max = 10, log = T)
  return(b0prior + b1prior + b2prior + sdprior)
}

# Posterior
p_target <- function(param){
  return(likelihood(param) + prior(param))
}


proposalfunction <- function(param,alpha=1){
  return(c(rnorm(3,mean = param[1:3],sd = 0.25),max(1e-3,runif(1,param[3]-alpha,param[3]+alpha))))
}

startvalue <- c(1,-3,5,1)

# Define Bernoulli sampler as special case of binomial distribution
rbern <- function(p){
  rbinom(1,1,p)
}


run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations + 1,4))
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
N_it <- 30000
burnIn <- 5000
chain = run_metropolis_MCMC(startvalue, N_it)


chain2 <- data.frame(chain[,1],chain[,2],chain[,3],chain[,4])
chain2$burn <- as.factor(c(rep("burnin",burnIn + 1),rep("chain",N_it-burnIn)))
colnames(chain2) <- c("b0","b1","b2","sd","iteration")


ggplot(chain2,aes(x=b1,y=b2,colour=iteration)) +
  geom_path() + scale_color_tableau()  + 
  xlab(expression(b[1])) + ylab(expression(b[2])) +
  theme_xkcd() +
  theme(panel.background = element_rect(color = "black", fill = "gray85") ) 



