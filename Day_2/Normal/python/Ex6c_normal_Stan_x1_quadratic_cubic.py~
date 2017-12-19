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
x1 = uniform.rvs(loc=-10, scale=20, size=nobs)          # random uniform variable

xb = -3 + 6 * x1 - 1.75 * pow(x1, 2) + 0.25 * pow(x1, 3)                                       # linear predictor, xb
y = np.random.normal(loc=xb, scale=20.0, size=nobs)    # create y as adjusted
                                                      # random normal variate 

# plot data
plt.figure(figsize=(12,8))
plt.scatter(x1, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


# Fit
toy_data = {}                                                 # build data dictionary
toy_data['N'] = nobs                                          # sample size
toy_data['X'] = sm.add_constant(np.transpose(x1))             # explanatory variable
toy_data['Y'] = y                                             # response variable
toy_data['M'] = 500                                           # number of data for prediction
toy_data['K'] = 4                                             # number of parameters

toy_data['xx'] = np.arange(min(x1), max(x1),((max(x1)-min(x1))/500))      # grid for prediction


# JAGS code
NORM =""" model{
    # Diffuse normal priors for predictors
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001) }

    # Uniform prior for standard deviation
    tau <- pow(sigma, -2)       # precision
    sigma ~ dunif(0, 100)       # standard deviation

    # Likelihood function 
    for (i in 1:N){
        Y[i] ~ dnorm(mu[i],tau)
        mu[i]  <- eta[i]
        eta[i] <- beta[1] + beta[2] * X[i,2] + beta[3] * X[i,2]^2 + beta[4] * X[i,2]^3
    }

    # Prediction for new data
    for (j in 1:M){
        etax[j] <- beta[1] + beta[2] * xx[j] + beta[3] * xx[j]^2 + beta[4] * xx[j]^3
        mux[j]  <- etax[j]
        Yx[j] ~ dnorm(mux[j],tau)
    }
}"""

model = pyjags.Model(NORM, data=toy_data, chains=3)
samples = model.sample(5000, vars=['beta', 'sigma', 'Yx', 'mux'])


def summary(samples, varname, p=95):
    if varname == 'beta':
        for k in range(4):
            values = samples[varname][k]
            ci = np.percentile(values, [100-p, p])
            print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
                   varname[k], np.mean(values), p, *ci))

    else:
        values = samples[varname]
        ci = np.percentile(values, [100-p, p])
        print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
               varname, np.mean(values), p, *ci))

for varname in ['beta', 'sigma']:
    summary(samples, varname)


# get Gaussian fit
mean = []
std = []
beta_axis = []
for j in range(4):
    meanx, stdx = norm.fit(samples['beta'][j][:,0])
    mean.append(meanx)
    std.append(stdx)
    beta_axis.append(np.arange(min(samples['beta'][j][:,0]), max(samples['beta'][j][:,0]), 0.001))


# plot posteriors
plt.figure(figsize=(8,8))
plt.subplot(4,4,1)
plt.hist(samples['beta'][0][:,0], normed=True)
plt.plot(beta_axis[0], norm.pdf(beta_axis[0], loc=mean[0], scale=std[0]), color='red', lw=2)
plt.xlabel('beta[0]')

plt.subplot(4,4,5)
plt.scatter(samples['beta'][0][:,0], samples['beta'][1][:,0])
plt.xlabel('beta[0]')
plt.ylabel('beta[1]')

plt.subplot(4,4,6)
plt.hist(samples['beta'][1][:,0], normed=True)
plt.plot(beta_axis[1], norm.pdf(beta_axis[1], loc=mean[1], scale=std[1]), color='red', lw=2)
plt.xlabel('beta[1]')

plt.subplot(4,4,9)
plt.scatter(samples['beta'][0][:,0], samples['beta'][2][:,0])
plt.xlabel('beta[0]')
plt.ylabel('beta[2]')

plt.subplot(4,4,10)
plt.scatter(samples['beta'][1][:,0], samples['beta'][2][:,0])
plt.xlabel('beta[1]')
plt.ylabel('beta[2]')

plt.subplot(4,4,11)
plt.hist(samples['beta'][2][:,0], normed=True)
plt.plot(beta_axis[2], norm.pdf(beta_axis[2], loc=mean[2], scale=std[2]), color='red', lw=2)
plt.xlabel('beta[2]')

plt.subplot(4,4,13)
plt.scatter(samples['beta'][0][:,0], samples['beta'][3][:,0])
plt.xlabel('beta[0]')
plt.ylabel('beta[3]')

plt.subplot(4,4,14)
plt.scatter(samples['beta'][1][:,0], samples['beta'][3][:,0])
plt.xlabel('beta[1]')
plt.ylabel('beta[3]')

plt.subplot(4,4,15)
plt.scatter(samples['beta'][2][:,0], samples['beta'][3][:,0])
plt.xlabel('beta[2]')
plt.ylabel('beta[3]')

plt.subplot(4,4,16)
plt.hist(samples['beta'][3][:,0], normed=True)
plt.plot(beta_axis[3], norm.pdf(beta_axis[3], loc=mean[3], scale=std[3]), color='red', lw=2)
plt.xlabel('beta[3]')

plt.tight_layout()
plt.show()


