# Confidence vs Bayesian intervals 
# Truncated exponential 
require(ggplot2)
require(xkcd)
require(R2jags)
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

nobs <- 1e3
covs <- data.frame(id = 1:nobs)
y <- rexp(nobs) + 10 # Maybe there is a better way of doing this

# Check 
hist(y)



# Fit a survival model 
SurvObj <-  Surv(time=y)
model.par <- survreg(SurvObj~1, dist="exponential")
mu = exp(model.par$coeff)
CI <- exp(confint(model.par))
print(CI)

# Frequentist data
time <- seq(min(y),max(y),by=1)
S.t = exp(-(time)/mu)
survs.data <- data.frame(time=time, 
                         survs=S.t)





# Construct data dictionary
model.data <- list(Y = y,                               # Response variable
                   N = nobs,                           # Sample size
                   new.t = time,
                   new.N = length(time),
                   minD = min(y)
)

Exp <- "model{
lambda <- 1/mu
mu ~ dunif(1e-3, minD)                                        # standard deviation
# Likelihood function

for (i in 1:N){
Y[i] ~ dexp(lambda)
}

# Prediction
for (i in 1:new.N ){
S.t[i] <- exp(-(new.t[i])/mu)
}

}"

inits <- function() {
  list(mu = runif(1, 0, 10))
}
params <- c("mu","S.t")

fit <- jags(data = model.data,
                inits = inits,
                parameters.to.save  = params,
                model.file = textConnection(Exp),
                n.chains = 3,
                n.iter = 5000,
                n.thin = 1,
                n.burnin = 1000)


Sfit <- jagsresults(x = fit , params=c("S.t"),probs=c(0.005,0.025, 0.25, 0.5, 0.75, 0.975,0.995))

gdata <- data.frame(xx = time, mean = Sfit[,"mean"],lwr1=Sfit[,"25%"],lwr2=Sfit[,"2.5%"],lwr3=Sfit[,"0.5%"],upr1=Sfit[,"75%"],
                     upr2=Sfit[,"97.5%"],upr3=Sfit[,"99.5%"])



ggplot(gdata ,aes(x=xx)) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr3, ymax=upr3,y= NULL),alpha=0.8,fill=c("#f0f0f0"),show.legend=FALSE)+
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL),alpha=0.7,  fill = c("#bdbdbd"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL),alpha=0.6,fill=c("#636363"),show.legend=FALSE)  +
  #
theme_xkcd() +
  geom_line(survs.data ,mapping=aes(x=time,y=S.t))

hist(as.matrix(as.mcmc(fit)[,"mu"]))



