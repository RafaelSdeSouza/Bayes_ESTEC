# From: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication
#
# Code 10.3 - Gaussian linear mixed model (in R using JAGS) for modeling the relationship
#             between type Ia supernovae host galaxy mass and Hubble residuals
# 
# real data from Wolf et al., 2016 - http://adsabs.harvard.edu/abs/2016ApJ...821..115W
# 1 response (y): Hubble residuals 
# 1 explanatory variable (x1): host galaxy mass


require(R2jags)
require(ggplot2)
require(mcmcplots)
source("../auxiliar_functions/jagsresults.R")

library(R2jags)

# Data
path_to_data = "https://raw.githubusercontent.com/astrobayes/BMAD/master/data/Section_10p2/HR.csv"

dat <- read.csv(path_to_data, header = T)

# Prepare data to JAGS
nobs = nrow(dat)
obsx1 <- dat$LogMass
errx1 <- dat$e_LogMass
obsy <- dat$HR
erry <- dat$e_HR
type <- as.numeric(dat$Type) # convert class to numeric flag 1 or 2

jags_data <- list(
  obsx1 = obsx1,
  obsy = obsy,
  errx1 = errx1,
  erry = erry,
  K = 2,
  N = nobs,
  type = type)

# Fit
NORM_errors <-" model{
    tau0 ~ dunif(1e-1,5)
    mu0 ~ dnorm(0,1)

    # Diffuse normal priors for predictors
    for (j in 1:2){
        for (i in 1:K) {
            beta[i,j] ~  dnorm(mu0, tau0)
        }
    }

    # Gamma prior for standard deviation
    tau ~ dgamma(1e-3, 1e-3) # precision
    sigma <- 1 / sqrt(tau) # standard deviation

    # Diffuse normal priors for true x
    for (i in 1:N){
        x1[i] ~ dnorm(0,1e-3)
    }

    # Likelihood function
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tau)
        obsx1[i]  ~ dnorm(x1[i],pow(errx1[i],-2))
        mu[i] <- beta[1,type[i]] + beta[2,type[i]] * x1[i]
    }
}"

inits <- function () {
    list(beta = matrix(rnorm(4, 0, 0.01),ncol = 2))
}

params0 <- c("beta", "sigma")                        # define parameters

# Run MCMC
NORM <- jags(
             data = jags_data,
             inits = inits,
             parameters = params0,
             model = textConnection(NORM_errors),
             n.chains = 3,
             n.iter = 100000,
             n.thin = 1,
             n.burnin = 60000)

# Output
print(NORM,justify = "left", digits=3)

# Diagnostics 
# Try to increase n.iter and n.burnin is to improve mixing

# plot chains
traplot(NORM,c("beta", "sigma"))

# plot posteriors
denplot(NORM,c("beta", "sigma"))

# plot all coefficients
caterplot(NORM,c("beta", "sigma"))

# Look at output
jagsresults(NORM, params=c("beta", "sigma"),signif=2)


# Plot fitted values
yx <- jagsresults(NORM, params=c('Yx'))
gdata <- data.frame(x =xx, mean = yx[,"50%"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])
nmod <- data.frame(obsx1,obsy)

# plot  
ggplot(nmod,aes(x=obsx1,y=obsy))+
  geom_errorbar(aes(ymin=obsy-erry,ymax=obsy+erry),alpha=0.5,width=0.025,colour="gray70")+
  geom_errorbarh(aes(xmin=obsx1-errx1,
                     xmax=obsx1+errx1),alpha=0.5,height=0.025,colour="gray70")+
  geom_point(size=2,color="blue")+
  geom_ribbon(data=gdata,aes(x=xx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.15, fill=c("orange")) +
  geom_ribbon(data=gdata,aes(x=xx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("orange3")) +
  geom_line(data=gdata,aes(x=xx,y=mean),colour="red",linetype="dashed",size=1.5)+
  theme_bw()+
  ylab(expression(mu[SN]-mu[z]))+
  xlab(expression(log~M/M['\u0298']))+
  coord_cartesian(xlim=c(8,13),ylim=c(-1,1))


