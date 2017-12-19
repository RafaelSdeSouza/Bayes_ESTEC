# Developed by Emille Ishida in Dec 2017
# ESTEC Bayesian course - Netherlands
#
# Adapted from: Bayesian Models for Astrophysical Data, Cambridge Univ. Press
# (c) 2017,  Joseph M. Hilbe, Rafael S. de Souza and Emille E. O. Ishida 
# 
# Normal linear model in Python using Stan
# 1 response (y) and 1 explanatory variable (x1)

import numpy as np
import pystan
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
toy_data['xx'] = np.arange(min(x1), max(x1),((max(x1)-min(x1))/500))      # grid for prediction
toy_data['M'] = toy_data['xx'].shape[0]                                   # number of data for prediction
toy_data['K'] = 4                                             # number of parameters


# Stan code
stan_code =""" 
data {
    int<lower=0> N;           
    int<lower=0> M;     
    int<lower=0> K;                 
    matrix[N, 2] X;                       
    vector[N] Y;                       
    vector[M] xx;
}
parameters {
    vector[K] beta;                                              
    real<lower=0> sigma;               
}
model {
    vector[N] mu;
    vector[N] eta;

    for (i in 1:N){
        eta[i] = beta[1] + beta[2] * X[i,2] + beta[3] * X[i,2]^2 + beta[4] * X[i,2] ^3;
        mu[i] = eta[i];
    }

    Y ~ normal(mu, sigma);             # Likelihood function
}
generated quantities{
    vector[M] ypred;

    for (i in 1:M){
        ypred[i] = beta[1] + beta[2] * xx[i] + beta[3] * xx[i]^2 + beta[4] * xx[i]^3;
    }
}"""

fit = pystan.stan(model_code=stan_code, data=toy_data, iter=5000, chains=3, verbose=False, n_jobs=3)

# Output
nlines = 9                    # number of lines in screen output

output = str(fit).split('\n')
for item in output[:nlines]:
    print(item)   



# Plot posteriors
fit.plot(['beta', 'sigma'])
plt.tight_layout()
plt.show()


# get Gaussian fit
mean = []
std = []
beta_axis = []
for j in range(4):
    meanx, stdx = norm.fit(fit.extract(['beta'])['beta'][:,j])
    mean.append(meanx)
    std.append(stdx)
    beta_axis.append(np.arange(min(fit.extract(['beta'])['beta'][:,j]), max(fit.extract(['beta'])['beta'][:,j]), 0.001))


ypred = norm.fit(fit.extract(['ypred'])['ypred'][5000:])


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


# plot chains
plt.figure(figsize=(8, 12))

for k in range(4):
    plt.subplot(4,2, (2*k+1))
    plt.plot(beta_axis[k], norm.pdf(beta_axis[k], loc=mean[k], scale=std[k]), color='blue', lw=2)
    plt.xlabel('beta[' + str(k) + ']')

    plt.subplot(4,2, 2*k + 2)
    plt.plot(range(2500), fit.extract(['beta'])['beta'][:,k][5000:])
    plt.ylabel('beta[' + str(k) + ']')
    plt.xlabel('samples')

plt.tight_layout()
plt.show()

# plot prediction
plt.figure()
plt.scatter(x1, y, color='blue')
plt.plot(toy_data['xx'], fit.extract(['ypred'])['ypred'][-1], ls='--', color='black', lw=2)
plt.xlabel('x')
plt.ylabel('y')
plt.show()
