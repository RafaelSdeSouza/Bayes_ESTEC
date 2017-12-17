require(rjags)
require(R2jags)
require(latex2exp)
# Generate data
set.seed(1056)                      # set seed to replicate example
nobs = 100                          # number of obs in model
lambda <- runif(nobs,1e-1,2)               # random uniform variable

h <- 4.135e-15#eV.s
pi <- 3.141592654
c <- 2.997e10 * 1e4 #cm/s
k <- 8.6173303e-5#eV/K
Temp <-  5777#K 

sd <- 0.25   

mu <- (8*pi*h*c/lambda^5)*1/(exp(h*c/(k*Temp*lambda)) - 1)
muL <- log(mu) - 0.5*log(1 + (sd/mu)^2)

                       # Standard deviation
y <- rlnorm(nobs, muL, sdlog = sd)        # create y as  random normal variate, not quite right as I explain later
    
plot(lambda,y)


# MCMC solution via JAGS

model.data <- list(Y = y,               # Response variable
                   X = lambda,          # Predictors
                   N = nobs,             # Sample size
                   h =  4.135e-15, #erg.s
                   pi = 3.141592654,
                   c = 2.997e10 * 1e4,#cm/s
                   k = 8.6173303e-5
)

# Model set up
NORM <- "model{

# Uniform prior for standard deviation
tau ~ dgamma(0.01,0.01)

# Prior for temperature

Temp ~ dgamma(0.01,0.01)

# Likelihood function
for (i in 1:N){
Y[i] ~ dlnorm(muL[i],pow(tau,-2))
muL[i] <- log(eta[i]) - 0.5*log(1+(pow(tau/eta[i],2)))  
eta[i] <- (8*pi*h*c/pow(X[i],5))*1/(exp(h*c/(k*Temp*X[i])) - 1)
}
}"

# Initial values
inits <- function () {list(
  Temp = runif(1,1e2,1e4))
}

# Parameters to be displayed
params <- c("Temp", "tau","eta")

    
 
# MCMC
normfit <- jags(data = model.data,
                inits = inits,
                parameters = params,
                model = textConnection(NORM),
                n.chains = 3,
                n.iter = 10000,
                n.thin = 10,
                n.burnin = 5000)

print(normfit, intervals = c(0.025, 0.975), digits = 2)

as.mcmc(normfit)

jagsresults(x=normfit , params=c('Temp'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

y <- jagsresults(x=normfit , params=c('eta'),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))
xx <- lambda
gdata <- data.frame(x =xx, mean = y[,"mean"],lwr1=y[,"25%"],lwr2=y[,"2.5%"],lwr3=y[,"0.5%"],upr1=y[,"75%"],
                    upr2=y[,"97.5%"],upr3=y[,"99.5%"])

gobs <- data.frame(lambda,mu)



ggplot(gobs,aes(x=lambda,y=mu))+
  #  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),alpha=0.7,fill=c("gray70"),show.legend=FALSE)+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), fill = "#5DADE2",show.legend=FALSE) +
  #  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),alpha=0.4,fill=c("gray30"),show.legend=FALSE) +
  geom_point(color="#F1948A",size=2)+
  geom_line(data=gdata,aes(x=xx,y=mean),linetype="dashed",size=0.5,show.legend=FALSE) +
  theme_xkcd() +
  xlab(TeX('$\\mu m$')) +
  ylab(expression(B[lambda])) +
  theme(panel.background = element_rect(color = "black", fill = "gray85") ) 
 
  

