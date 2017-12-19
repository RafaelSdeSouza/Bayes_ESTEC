# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2016, Cambridge Univ. Press
#
# Example of Bayesian negative binomial regression in R using JAGS
# 
# real data from Harris et al.,2011
# complete catalog available at http://www.physics.mcmaster.ca/~harris/GCS_table.txt
#
# 1 response (y): nunmber of globular clusters
# 1 explanatory variable (x1): visual magnitude

require(R2jags)
require(ggplot2)
require(scales)
source("../../auxiliar_functions/jagsresults.R")

path_to_data = "../data/GCs.csv"

# Read data
NGC_dat = read.csv(file=path_to_data,header = T,dec=".")

# Look at the data
head(GC_dat)

# Prepare data for prediction 
M=500
xx = seq(from = 1.05 * min(NGC_dat$MV_T), 
         to = 0.875 * max(NGC_dat$MV_T), 
         length.out = M)


# prepare data for JAGS input
jags_data <- list(
  X = NGC_dat$MV_T,
  Y = NGC_dat$N_GC,
  N = nrow(NGC_dat),
  xx = xx,
  M = M, 
  K = 2
)

# JAGS model
model.NB <- "model{

# Diffuse normal priors betas
for (i in 1:K) { beta[i] ~ dnorm(0, 1e-5)}

# Prior for size
theta ~ dgamma(1e-3,1e-3)

# Likelihood
for (i in 1:N){
eta[i]<-beta[1]+beta[2]*X[i]
mu[i] <- exp(eta[i])
p[i]<-theta/(theta+mu[i])
Y[i]~dnegbin(p[i],theta)
}

# Prediction for new data
for (j in 1:M){
etax[j]<-beta[1]+beta[2]*xx[j]
mux[j] <- exp(etax[j])
px[j]<-theta/(theta+mux[j])
Yx[j]~dnegbin(px[j],theta)
}
}"

# set initial values
inits  <- function() {
  list(
    beta  = rnorm(2, 0, 0.1),
    theta = runif(1,0.1,5) 
  )  
}

# set parameters
params <- c("beta","theta","Yx","mux")


# fit
NB <- jags(data       = jags_data ,
           inits      = inits,
           parameters = params,
           model      = textConnection(model.NB),
           n.thin     = 1,
           n.chains   = 3,
           n.burnin   = 3500,
           n.iter     = 7000)

# numerical results
jagsresults(NB, params=c("beta","theta"))

# plot
yx <- jagsresults(NB, params=c('Yx'))
gdata <- data.frame(x=xx, mean = yx[,"50%"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

ggplot(NGC_dat,aes(x=MV_T,y=N_GC))+
  geom_ribbon(data=gdata,aes(x=x,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.45, fill=c("gray90"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=x,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.35, fill = c("gray45"),show.legend=FALSE) +
  geom_point(colour="red",size=1,alpha=0.85)+
  geom_line(data=gdata,aes(x=x,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),
                     labels=c("0",expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
  scale_shape_manual(name="",values=c(19,2,8,10))+
  scale_x_reverse()+
  theme_bw() +
  ylab(expression(N[GC]))+
  xlab(expression(M[V]))
