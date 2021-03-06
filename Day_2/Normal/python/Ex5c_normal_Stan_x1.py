# Developed by Emille Ishida on Dec 2017
# ESTEC Bayesian course, Netherlands
#
# Adapted from: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication

# Code 4.3 - Normal linear model in Python using STAN
# 1 response (y) and 1 explanatory variable (x1)

import numpy as np
import pystan
from scipy.stats import uniform
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
toy_data['nobs'] = nobs        # sample size
toy_data['x'] = x1             # explanatory variable
toy_data['y'] = y              # response variable
toy_data['xx'] = np.arange(min(x1), max(x1), (max(x1) - min(x1))/500)    # exploratory values for prediction
toy_data['M'] = 500                                                      # number of points for prediction

# STAN code
stan_code = """
data {
    int<lower=0> nobs;           
    int<lower=0> M;                      
    vector[nobs] x;                       
    vector[nobs] y;                       
    vector[M] xx;
}
parameters {
    real beta0;
    real beta1;                                                
    real<lower=0> sigma;               
}
model {
    vector[nobs] mu;

    mu = beta0 + beta1 * x;

    y ~ normal(mu, sigma);             # Likelihood function
}
generated quantities{
    vector[M] ypred;

    for (i in 1:M){
        ypred[i] = beta0 + beta1 * xx[i];
    }
}
"""

fit = pystan.stan(model_code=stan_code, data=toy_data, iter=5000, chains=3, verbose=False, n_jobs=3)

# Output
nlines = 9                    # number of lines in screen output

output = str(fit).split('\n')
for item in output[:nlines]:
    print(item)   


# Plot posteriors
fit.plot(['beta0', 'beta1', 'sigma'])
plt.tight_layout()
plt.show()

# plot prediction
plt.figure(figsize=(12,8))
plt.plot(toy_data['xx'], 2 + 3 * toy_data['xx'], color='black', lw=2)
plt.scatter(x1, y, color='blue')
plt.xlabel('x', fontsize=18)
plt.ylabel('y', fontsize=18)
plt.tight_layout()
plt.show()

