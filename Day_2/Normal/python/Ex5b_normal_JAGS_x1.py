# Developed by Emille Ishida in Dec 2017
# ESTEC Bayesian course - Netherlands
#
# Adapted from: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# Normal linear model in Python using JAGS
# 1 response (y) and 1 explanatory variable (x1)

import numpy as np
import pyjags
from scipy.stats import uniform, norm
import pylab as plt

# Data
np.random.seed(1056)                                 # set seed to replicate example
nobs= 500                                            # number of obs in model 
x1 = uniform.rvs(loc=2, scale=3, size=nobs)          # random uniform variable

xb = 2 + 3 * x1                                       # linear predictor, xb
y = np.random.normal(loc=xb, scale=2.0, size=nobs)    # create y as adjusted
                                                      # random normal variate 

# plot data
plt.figure(figsize=(12,8))
plt.scatter(x1, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


# Fit
toy_data = {}                  # build data dictionary
toy_data['N'] = nobs           # sample size
toy_data['X'] = x1             # explanatory variable
toy_data['Y'] = y              # response variable
toy_data['M'] = 5000           # number of data for prediction

toy_data['xx'] = np.arange(min(x1), max(x1),((max(x1)-min(x1))/5000))      # grid for prediction


# JAGS code
NORM =""" model{
    # Diffuse normal priors for predictors
    beta0 ~ dnorm(0, 0.0001)
    beta1 ~ dnorm(0, 0.0001)    

   # Uniform prior for standard deviation
    sigma ~ dunif(0, 100)       # standard deviation
    tau <- pow(sigma, -2)       # precision
   
   
    # Likelihood function 
    for (i in 1:N){
       Y[i] ~ dnorm(mu[i], tau)
       mu[i]  <- eta[i]
       eta[i] <- beta0 + beta1 * X[i]
    }

   # Prediction for new data
   for (j in 1:M){
   etax[j] <- beta0 + beta1 * xx[j]
   mux[j]  <- etax[j]
   Yx[j] ~ dnorm(mux[j],tau)
    }
}"""


model = pyjags.Model(NORM, data=toy_data, chains=3)
samples = model.sample(5000, vars=['beta0', 'beta1', 'sigma', 'Yx', 'mux'])


def summary(samples, varname, p=95):
    values = samples[varname]
    ci = np.percentile(values, [100-p, p])
    print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
      varname, np.mean(values), p, *ci))

for varname in ['beta0', 'beta1', 'sigma']:
    summary(samples, varname)


# get Gaussian fit
beta0_mean, beta0_std = norm.fit(samples['beta0'][0][:,0])
beta1_mean, beta1_std = norm.fit(samples['beta1'][0][:,0])

beta0_axis = np.arange(min(samples['beta0'][0][:,0]), max(samples['beta0'][0][:,0]), 0.01)
beta1_axis = np.arange(min(samples['beta1'][0][:,0]), max(samples['beta1'][0][:,0]), 0.01)

# plot posteriors
plt.figure(figsize=(8,8))
plt.subplot(2,2,1)
plt.hist(samples['beta0'][0][:,0], normed=True)
plt.plot(beta0_axis, norm.pdf(beta0_axis, loc=beta0_mean, scale=beta0_std), color='red', lw=2)
plt.xlabel('beta0')

plt.subplot(2,2,3)
plt.scatter(samples['beta0'][0][:,0], samples['beta1'][0][:,0])
plt.xlabel('beta0')
plt.ylabel('beta1')

plt.subplot(2,2,4)
plt.hist(samples['beta1'][0][:,0], normed=True)
plt.plot(beta1_axis, norm.pdf(beta1_axis, loc=beta1_mean, scale=beta1_std), color='red', lw=2)
plt.xlabel('beta1')

plt.tight_layout()
plt.show()


