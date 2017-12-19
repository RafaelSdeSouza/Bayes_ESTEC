# From: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication

# Code 10.1 - Normal linear model in R using JAGS for accessing the relationship between
#             central black hole mass and bulge velocity dispersion

require(R2jags)
# Data
path_to_data = "data/M_sigma.csv"

# Read data
MS < -read.csv(path_to_data,header = T)

# Identify variables
N <- nrow(MS) # number of data points
obsx <- MS$obsx # log black hole mass
errx <- MS$errx # error on log black hole mass
obsy <- MS$obsy # log velocity dispersion
erry <- MS$erry # error on log velocity dispersion

# Prepare data for prediction 
M=500
xx = seq(from =  min(obsx),
         to =  max(obsx),
         length.out = M)



# Prepare data to JAGS
MS_data <- list(
  obsx = obsx,
  obsy = obsy,
  errx = errx,
  erry = erry,
  N = N,
  M = M,
  xx  = xx
)



# Fit
NORM_errors <-"model{
# Diffuse normal priors for predictors
alpha ~ dnorm(0,1e-3)
beta ~ dnorm(0,1e-3)

# Gamma prior for scatter
tau ~ dgamma(1e-3,1e-3) # precision
epsilon <- 1/sqrt(tau) # intrinsic scatter

# Diffuse normal priors for true x
for (i in 1:N){ x[i] ~ dnorm(0,1e-3) }

for (i in 1:N){
obsx[i] ~ dnorm(x[i], pow(errx[i], -2))
obsy[i] ~ dnorm(y[i], pow(erry[i], -2)) # likelihood function
y[i] ~ dnorm(mu[i], tau)
mu[i] <- alpha + beta*x[i] # linear predictor
}

# Prediction for new data
    for (j in 1:M){
    etax[j]<-alpha+beta*xx[j]
    mux[j]  <- etax[j]
    Yx[j]~dnorm(mux[j],tau)
    }



}"

# Define initial values
inits <- function () {
  list(alpha = runif(1,0,10),
       beta = runif(1,0,10))
}

# Identify parameters
params0 <- c("alpha","beta", "epsilon","mux")

# Fit
NORM_fit <- jags(data = MS_data,
                 inits = inits,
                 parameters.to.save  = params0,
                 model.file  = textConnection(NORM_errors),
                 n.chains = 3,
                 n.iter = 50000,
                 n.thin = 10,
                 n.burnin = 30000
)
# Output
#print(NORM_fit,justify = "left", digits=2)


# Plot
yx <- jagsresults(x = NORM_fit, params=c('mux'))

normdata <- data.frame(obsx,obsy)
gdata <- data.frame(x =xx, mean = yx[,"mean"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])


ggplot(normdata,aes(x=obsx,y=obsy))+ geom_point(colour="#de2d26",size=1,alpha=0.35)+
  geom_point(size=1.5,colour="red3") +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("orange3"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.35, fill = c("orange"),show.legend=FALSE) +
  geom_line(data=gdata,aes(x=xx,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  theme_bw()


