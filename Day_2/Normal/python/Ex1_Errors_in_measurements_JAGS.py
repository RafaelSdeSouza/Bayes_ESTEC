# Developed by Emille Ishida on Dec 2017
# ESTEC Bayesian course
#
# Adapted from: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication

# Code 4.11 - Normal linear model in Python using JAGS and
#             including errors in variables

# 1 response (y) and 1 explanatory variable (x1)


import numpy as np
import statsmodels.api as sm
import pyjags

from scipy.stats import norm

############### Data
np.random.seed(1056)                      # set seed to replicate example
nobs = 750                               # number of obs in model 
sdobsx = 1.25
truex =  norm.rvs(0,2.5, size=nobs)     # normal variable
errx = norm.rvs(0, sdobsx, size=nobs)   # errors
obsx = truex + errx                     # observed

beta0 = -4
beta1 = 7           
sdy = 1.25
sdobsy = 2.5

erry = norm.rvs(0, sdobsy, size=nobs)
truey = norm.rvs(beta0 + beta1*truex, sdy, size=nobs)
obsy = truey + erry

# Fit
toy_data = {}                                # build data dictionary
toy_data['N'] = nobs                         # sample size
toy_data['obsx'] = obsx                      # explanatory variable       
toy_data['errx'] = errx                      # uncertainty in explanatory variable
toy_data['obsy'] = obsy                      # response variable
toy_data['erry'] = erry                      # uncertainty in response variable
toy_data['K'] = 2                            # number of parameters

# JAGS code
NORM_err = """ model{
# Diffuse normal priors for predictors
for (i in 1:K) { beta[i] ~ dnorm(0, 1e-3) }

# Uniform prior for standard deviation
tauy <- pow(sigma, -2)                               # precision
sigma ~ dunif(0, 100)                                # diffuse prior for standard deviation

# Diffuse normal priors for true x
for (i in 1:N){
x[i] ~ dnorm(0,1e-3)
}

# Likelihood
for (i in 1:N){
obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
y[i] ~ dnorm(mu[i],tauy)
obsx[i] ~ dnorm(x[i],pow(errx[i],-2))
mu[i] <- beta[1]+beta[2]*x[i]
}
}"""

# Run mcmc
model = pyjags.Model(NORM_err, data=toy_data, chains=3)
samples = model.sample(5000, vars=['beta', 'sigma'])


def summary(samples, varname, p=95):
    if varname == 'beta':
        for k in range(toy_data['K']):
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


