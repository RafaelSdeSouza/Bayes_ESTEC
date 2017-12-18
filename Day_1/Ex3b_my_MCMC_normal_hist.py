"""
Created by Emille Ishida in 12 November 2017. 

Example of simple MCMC algorithm.

Goal is to determine the mean of the original distribution.

Here we are assuming that:
	- we know that the correct model is Gaussian
 	- the uncertainty is known and constant

The first part of this script was taken from 
http://twiecki.github.io/blog/2015/11/10/mcmc-sampling/
"""



import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import norm, bernoulli, uniform, multivariate_normal

# fix random seed
np.random.seed(123)

# define fiducial parameters
true_mu = 2.0
std = 0.25
n = 50

# generate observations = true + gaussian noise
data = norm.rvs(loc=true_mu, scale=std, size=n)

# plot data
plt.figure()
plt.title('Data')
plt.hist(data)
plt.show()

# set starting point
mu_current = 10.0

# proposal size of step for random walk
proposal_width = 1.0

# number of samples to accept
n_samples = 500

# determine prior parameters
# here we consider  a gaussian prior
mu_prior_mean = 5 
mu_prior_std = 1.5

# store chain
chain = []

while len(chain) < n_samples:

    print str(len(chain)) + '  samples'

    # propose a new point
    mu_proposal = norm(mu_current, proposal_width).rvs()

    # calculate likelihood in current and proposed point
    likelihood_current = norm(mu_current, std).pdf(data).prod()
    likelihood_proposal = norm(mu_proposal, std).pdf(data).prod()

    # Compute prior probability of current and proposed mu        
    prior_current = norm(mu_prior_mean, mu_prior_std).pdf(mu_current)
    prior_proposal = norm(mu_prior_mean, mu_prior_std).pdf(mu_proposal)

    # Nominator of Bayes formula
    p_current = likelihood_current * prior_current
    p_proposal = likelihood_proposal * prior_proposal

    # acceptance probability
    p_accept = min([1, p_proposal / p_current])

    accept = bernoulli.rvs(p_accept)

    if accept:
        # store new point
        chain.append(mu_proposal)

        # Update position
        mu_current = mu_proposal

    else:
        chain.append(mu_proposal)


plt.figure()
plt.title('Chain convergence')
plt.plot(range(n_samples), chain)
plt.xlabel('sample')
plt.ylabel('mu')
plt.show()

# from the above plot you can see that the first 200 iterations where very
# far from the convergence point, so we can consider them burn-in

burn_in = 200

plt.figure()
plt.title('Chain after burn-in/warm-up')
plt.plot(range(n_samples - burn_in), chain[burn_in:])
plt.plot([0, n_samples - burn_in], [true_mu, true_mu], color='black', ls='--', lw=2.0)
plt.xlabel('sample')
plt.ylabel('mu')
plt.show()

# from this we can estimate
mean_posterior, std_posterior = norm.fit(chain[burn_in:])


x_axis = np.arange(min(chain[burn_in:]), max(chain[burn_in:]), 0.001)
y_axis = norm.pdf(x_axis, loc=mean_posterior, scale=std_posterior)

plt.figure()
plt.hist(chain[burn_in:], normed=True, alpha=0.5)
plt.plot(x_axis, y_axis, color='red', lw=2.0, label='posterior')
plt.legend()
plt.xlabel('mu')
plt.ylabel('PDF')
plt.show()


