# From: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication

# Code 4.3 - Normal linear model in Python using STAN
# 1 response (y) and 1 explanatory variable (x1)

import numpy as np
import statsmodels.api as sm
import pyjags
from scipy.stats import uniform

# Data
np.random.seed(1056)                 # set seed to replicate example
nobs= 5000                           # number of obs in model 
x1 = uniform.rvs(size=nobs)          # random uniform variable

x1.transpose()                   # create response matrix
X = sm.add_constant(x1)          # add intercept
beta = [2.0, 3.0]                # create vector of parameters

xb = np.dot(X, beta)                                  # linear predictor, xb
y = np.random.normal(loc=xb, scale=1.0, size=nobs)    # create y as adjusted
                                                      # random normal variate 

# Fit
toy_data = {}                  # build data dictionary
toy_data['N'] = nobs        # sample size
toy_data['X'] = x1             # explanatory variable
toy_data['Y'] = y              # response variable
toy_data['M'] = 5000
toy_data['K'] = 2
toy_data['xx'] = np.arange(min(x1), max(x1),((max(x1)-min(x1))/5000))


# JAGS code
NORM =""" model{
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
}"""


model = pyjags.Model(NORM, data=toy_data, chains=3)
samples = model.sample(5000, vars=['beta', 'sigma', 'Yx', 'mux'])


def summary(samples, varname, p=95):
    values = samples[varname]
    ci = np.percentile(values, [100-p, p])
    print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
      varname, np.mean(values), p, *ci))

for varname in ['beta', 'sigma', 'Yx', 'mux']:
    summary(samples, varname)
