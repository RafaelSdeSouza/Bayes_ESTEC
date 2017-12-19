# Developed by Emille Ishida on Dec 2017
# ESTEC Bayesian course
#
# Adapted from: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication
#
# Normal linear model, in Python using Stan, for assessing the relationship
#           between central black hole mass and bulge velocity dispersion
#
# Statistical Model: Gaussian regression considering errors in variables
#                    in Python using Stan
#
# Astronomy case: Relation between mass of galaxy central supermassive black hole
#                 and its stelar bulge velocity dispersion
#                 taken from Harris, Poole and Harris, 2013, MNRAS, 438 (3), p.2117-2130
#
# 1 response (obsy - mass) and 1 explanatory variable (obsx - velocity dispersion)
#
# Data from: http://www.physics.mcmaster.ca/~harris/GCS_table.txt

import numpy as np
import pandas as pd
import pyjags
from scipy.stats import norm
import pylab as plt

path_to_data = '../../data/M_sigma.csv'

# read data
data_frame = dict(pd.read_csv(path_to_data))

# prepare data for Stan
data = {}
data['obsx'] = np.array(data_frame['obsx'])
data['errx'] = np.array(data_frame['errx'])
data['obsy'] = np.array(data_frame['obsy'])
data['erry'] = np.array(data_frame['erry'])
data['N'] = len(data['obsx'])
data['M'] = 500
data['xx'] = np.arange(min(data['obsx']), max(data['obsx']), (max(data['obsx']) - min(data['obsx']))/500)

# JAGS Gaussian model with errors
jags_code = """model{
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
}"""

# Run mcmc
model = pyjags.Model(jags_code, data=data, chains=3)
samples = model.sample(5000, vars=['alpha', 'beta', 'epsilon', 'mux'])
9

def summary(samples, varname, p=95):
    values = samples[varname]
    ci = np.percentile(values, [100-p, p])
    print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
               varname, np.mean(values), p, *ci))

for varname in ['alpha', 'beta', 'epsilon']:
    summary(samples, varname)



# get Gaussian fit
alpha_mean, alpha_std = norm.fit(samples['alpha'][0][:,0])
alpha_mean, alpha_std = norm.fit(samples['alpha'][0][:,0])

alpha_axis = np.arange(min(samples['alpha'][0][:,0]), max(samples['alpha'][0][:,0]), 0.01)
alpha_axis = np.arange(min(samples['alpha'][0][:,0]), max(samples['alpha'][0][:,0]), 0.01)


beta_mean, beta_std = norm.fit(samples['beta'][0][:,0])
beta_mean, beta_std = norm.fit(samples['beta'][0][:,0])

beta_axis = np.arange(min(samples['beta'][0][:,0]), max(samples['beta'][0][:,0]), 0.01)
beta_axis = np.arange(min(samples['beta'][0][:,0]), max(samples['beta'][0][:,0]), 0.01)

# plot posteriors
plt.figure(figsize=(8,8))
plt.subplot(2,2,1)
plt.hist(samples['alpha'][0][:,0], normed=True)
plt.plot(alpha_axis, norm.pdf(alpha_axis, loc=alpha_mean, scale=alpha_std), color='red', lw=2)
plt.xlabel('alpha')

plt.subplot(2,2,3)
plt.scatter(samples['alpha'][0][:,0], samples['beta'][0][:,0])
plt.xlabel('alpha')
plt.ylabel('beta')

plt.subplot(2,2,4)
plt.hist(samples['beta'][0][:,0], normed=True)
plt.plot(beta_axis, norm.pdf(beta_axis, loc=beta_mean, scale=beta_std), color='red', lw=2)
plt.xlabel('beta')

plt.tight_layout()
plt.show()


