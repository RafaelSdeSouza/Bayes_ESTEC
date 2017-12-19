# Developed by Emille Ishida on Dec 2017
# ESTEC Bayesian course
#
# Adapted from Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# you are kindly asked to include the complete citation if you used this 
# material in a publication
#
# Code 10.13 Bernoulli model in Python using Stan, for assessing the 
#            relationship between bulge size and the fraction of red spirals.
#
# Statistical Model: Bernoulli model in Python using Stan
#
# Astronomy case: Relation between fraction of red spiral galaxies
#                 and  galaxy bulge size
#                 taken from Masters et al., 2010, MNRAS,  405 (2), 783-799
#
# 1 response (Y - galaxy type: red/blue spirals) 
# 1 explanatory variable (x - bulge size)
#
# Data from: http://data.galaxyzoo.org/data/redspirals/BlueSpiralsA2.txt
#            http://data.galaxyzoo.org/data/redspirals/RedSpiralsA1.txt

import numpy as np
import pandas as pd
import pystan 
import statsmodels.api as sm

# Data
path_to_data = '../data/Red_spirals.csv'

# read data
data_frame = dict(pd.read_csv(path_to_data))
x = np.array(data_frame['fracdeV'])

# for prediction
xx = np.arange(min(x), max(x), (max(x) - min(x))/500)

# prepare data for Stan
data = {}
data['X'] = sm.add_constant((x.transpose()))
data['Y'] = np.array(data_frame['type'])
data['nobs'] = data['X'].shape[0]
data['K'] = data['X'].shape[1]
data['M'] = 500
data['XX'] = sm.add_constant((xx.transpose()))



# Fit
stan_code="""
data{
    int<lower=0> nobs;                # number of data points
    int<lower=0> K;                   # number of coefficients
    int<lower=0> M;                   # number of points for prediction
    matrix[nobs, K] X;                # bulge size
    int Y[nobs];                      # galaxy type: 1 - red, 0 - blue
    matrix[M, K] XX;                     # exploratory variable for plotting
}
parameters{
    vector[K] beta;                   # linear predictor coefficients
}
model{
    # priors and likelihood
    for (i in 1:K) beta[i] ~ normal(0, 100);

    Y ~ bernoulli_logit(X * beta);
}
generated quantities{
    vector[M] ypred;

    for (j in 1:M){
        ypred[j] = bernoulli_logit_rng(XX[j] * beta);
    }
}
"""

# Run mcmc
fit = pystan.stan(model_code=stan_code, data=data, iter=6000, chains=3,
                  warmup=3000, thin=1, n_jobs=3)

# Output
nlines = 7                                   # number of lines in screen output

output = str(fit).split('\n')
for item in output[:nlines]:
    print(item)   



