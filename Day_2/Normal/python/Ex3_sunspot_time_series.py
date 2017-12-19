# Developed by Emille Ishida on Dec 2017
# ESTEC Bayesian course
#
# Adapted from: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication
#
# Normal model in Python using Stan, 
#            for assessing the evolution of the number of sunspots 
#            through the years.
#
# Statistical Model: Time series in Python using Stan
#
# Astronomy case: Evolution of the number of sunspots with time
#                 from Lunn et al., 2012, 
#                 The BUGS Book: a Practical Introduction to
#                                Bayesian Analysis, CRC
#
# 1 response (Y - number of sunspots) 
# 1 explanatory variable (x - year)
#
# Data from: http://www.sidc.be/silso/DATA/EISN/EISN_current.csv

import numpy as np
import pystan 
import pylab as plt
import pandas as pd

# Data
path_to_data = "../../data/sunspot.csv"

# read data
data_frame = dict(pd.read_csv(path_to_data))

# prepare data for Stan
data = {}
data['Y'] = [round(item) for item in data_frame['nspots']]
data['N'] = len(data['Y'])
data['K'] = 2

# Fit
stan_code="""
data{
    int<lower=0> N;                # number of data points
    int<lower=0> K;                # number of coefficients
    real Y[N];                     # nuber of sunspots
}
parameters{
    vector[K] phi;                    # linear predictor coefficients
    real<lower=0> tau;                # noise parameter
}
model{
    real mu[N];
    
    mu[1] = Y[1];                     # set initial value

    # priors and likelihood
    tau ~ gamma(0.001, 0.001);
    for (i in 1:K) phi[i] ~ normal(0, 100);

    for (t in 2:N) mu[t] = phi[1] + phi[2] * Y[t - 1];
    Y ~ normal(mu, tau);
}
generated quantities{
    vector[N] new_mu;
    vector[N] ypred;

    new_mu[1] = Y[1];
    for (i in 2:N) new_mu[i] = phi[1] + phi[2] * Y[i - 1];
    for (j in 1:N) ypred[j] = normal_rng(new_mu[j], tau);
}
"""

# Run mcmc
fit = pystan.stan(model_code=stan_code, data=data, iter=7500, chains=3,
                  warmup=5000, thin=1, n_jobs=3)

# Output
nlines = 8                                   # number of lines in screen output

output = str(fit).split('\n')
for item in output[:nlines]:
    print(item) 


plt.figure(figsize=(15, 8))
plt.scatter(range(1700, 2016), data['Y'])
plt.plot(np.arange(1700, 2016), fit.extract(['ypred'])['ypred'][-1], color='black')
plt.xlabel('year', fontsize=18)
plt.ylabel('sunspots', fontsize=18)
plt.tight_layout()
plt.show()
