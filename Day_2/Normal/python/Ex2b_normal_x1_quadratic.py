# Developed by Emille E. O. Ishida on 18 Dec 2017
# ESTEC Bayesian course, Netherlands
#
# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2017, Cambridge Univ. Press
#
# Example of frequentist linear regression in Python
# synthetic data
# 1 response (y) and 1 explanatory variable (x)

import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pylab as plt

from scipy.stats import uniform, norm

# set random seed
np.random.seed(1056) 

nobs = 50                                             # number of observations
x = uniform.rvs(loc=0, scale=5, size=nobs)            # random uniform variable
mu = 1 + 5 * x + 0.75 * pow(x, 2)                     # linear predictor
y = norm.rvs(loc=mu, scale=1)                         # create y as adjusted random variate

# plot data
plt.figure()
plt.scatter(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# create data dictionary
my_data = {}
my_data['x'] = x
my_data['xx'] = pow(x, 2)
my_data['y'] = y

# fit a linear model using ordinary least squares
results = smf.glm(formula='y ~ 1 + x + xx', data=my_data, family=sm.families.Gaussian()).fit()

# output results
print(str(results.summary()))

# given fitted model, get predicted values
yfitted = results.predict()

# plot data + results
plt.figure()
plt.scatter(x, y, color='blue')
plt.plot(x, yfitted, color='red')
for i in range(len(x)):
    plt.plot([x[i], x[i]], [yfitted[i], y[i]], color='black')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
