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
import pyjags
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

# Fit
AR1_NORM = """model{
# Priors
sd2 ~ dgamma(1e-3,1e-3)
tau <- 1/sd2
sd <- sqrt(sd2)

for(i in 1:2){
    phi[i] ~ dnorm(0,1e-2)
}

mu[1] <- Y[1]

# Likelihood function
for (t in 2:N) {
    Y[t] ~ dnorm(mu[t],tau)
    mu[t] <- phi[1] + phi[2] * Y[t-1]
}

# Prediction
for (t in 1:N){
    Yx[t]~dnorm(mu[t],tau)
}
}"""

# Run mcmc
model = pyjags.Model(AR1_NORM, data=data, chains=3)
samples = model.sample(5000, vars=['sd', 'phi', 'Yx'])


def summary(samples, varname, p=95):
    if varname == 'phi':
        for k in range(2):
            values = samples[varname][k]
            ci = np.percentile(values, [100-p, p])
            print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
                   varname[k], np.mean(values), p, *ci))

    else:
        values = samples[varname]
        ci = np.percentile(values, [100-p, p])
        print('{:<6} mean = {:>5.1f}, {}% credible interval [{:>4.1f} {:>4.1f}]'.format(
               varname, np.mean(values), p, *ci))

for varname in ['sd', 'phi']:
    summary(samples, varname)

    print(item) 


plt.figure(figsize=(15, 8))
plt.scatter(range(1700, 2016), data['Y'])
plt.plot(np.arange(1700, 2016), samples['Yx'][:,0][:,0], color='black')
plt.xlabel('year', fontsize=18)
plt.ylabel('sunspots', fontsize=18)
plt.tight_layout()
plt.show()
