# From: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication
#
# Example of Bayesian normal linear regression in R using JAGS
# synthetic data
#
# 1 response (y) and 1 explanatory variable (x1)

require(R2jags)
source("../auxiliar_functions/jagsresults.R")
require(ggplot2)

set.seed(1056)                 # set seed to replicate example
nobs = 100                     # number of obs in model 
x1 <- rnorm(nobs,3,1)          # random normal variable
#x1 <- runif(nobs,2,5)         # random uniform variable

xb <- 2 + 3*x1                 # linear predictor
y <- rnorm(nobs, xb, sd=2)     # create y as adjusted random normal variate

plot(x1,y)

# Prepare data for prediction 
M=500
xx = seq(from =  min(x1), 
         to =  max(x1), 
         length.out = M)

jags_data <- list(Y = y,
                  X  = x1,
                  K  = 2,
                  N  = nobs,
                  M = M,
                  xx= xx)
    

# JAGS model
NORM <-" model{
    # Diffuse normal priors for predictors
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }
    
   # Uniform prior for standard deviation
    sigma ~ dunif(0, 100)       # standard deviation
    tau <- pow(sigma, -2)       # precision
   
   
    # Likelihood function 
    for (i in 1:N){
       Y[i]~dnorm(mu[i],tau)
       mu[i]  <- eta[i]
       eta[i] <- beta[1]+beta[2]*X[i]
    }

   # Prediction for new data
   for (j in 1:M){
   etax[j]<-beta[1]+beta[2]*xx[j]
   mux[j]  <- etax[j]
   Yx[j]~dnorm(mux[j],tau)
    }
}"

# set initial values
inits <- function () {
  list(
    beta = rnorm(2, 0, 0.01))
}

# define parameters
params <- c("beta", "sigma","Yx","mux")
#params <- c("beta", "sigma")


jagsfit <- jags(
           data       = jags_data,
           inits      = inits,
           parameters = params,
           model      = textConnection(NORM),
           n.chains   = 3,
           n.iter     = 5000,
           n.thin     = 1,
           n.burnin   = 2500)

# Plot
yx <- jagsresults(x=jagsfit, params=c('mux'))

normdata <- data.frame(x1,y)
gdata <- data.frame(x =xx, mean = yx[,"mean"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])


ggplot(normdata,aes(x=x1,y=y))+ geom_point(colour="#de2d26",size=1,alpha=0.35)+
  geom_point(size=1.5,colour="red3")+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("orange3"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.35, fill = c("orange"),show.legend=FALSE) +
  geom_line(data=gdata,aes(x=xx,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw()
 
