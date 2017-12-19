# Developed by Emille E. O. Ishida on 18 Dec 2017
# ESTEC Bayesian course, Netherlands
#
# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2017, Cambridge Univ. Press
#
# Example of frequentist linear regression in Python
# synthetic data
# 1 response (y) and 1 explanatory variable (x1, x2)

import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.stats import uniform, norm

# set random seed
np.random.seed(1056) 

nobs = 150                                             # number of observations
x1 = uniform.rvs(loc=0, scale=2, size=nobs)            # random uniform variables
x2 = uniform.rvs(loc=0, scale=2, size=nobs)

xb = 2 + 3 * x1 + 1.5 * x2                            # linear predictor
y = norm.rvs(loc=xb, scale=1)                         # create y as adjusted random variate

# plot data
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1, x2, y)
ax.set_xlabel('x1', fontsize=18)
ax.set_ylabel('x2', fontsize=18)
ax.set_zlabel('y', fontsize=18)
plt.show()

# create data dictionary
my_data = {}
my_data['x1'] = x1
my_data['x2'] = x2
my_data['y'] = y

# fit a linear model using ordinary least squares
results = smf.glm(formula='y ~ 1 + x1 + x2', data=my_data, family=sm.families.Gaussian()).fit()

# output results
print(str(results.summary()))

# given fitted model, get predicted values
yfitted = results.predict()

# plot data + results
x1_axis = np.arange(0, 2, 0.01)
x2_axis = np.arange(0, 2, 0.01)
XX, YY = np.meshgrid(x1_axis, x2_axis)
Z = results.params[0] + results.params[1] * x1_axis + results.params[2] * x2_axis

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(XX, YY, Z, color='orange', alpha=0.5)
ax.scatter(x1, x2, y, color='blue')
ax.set_xlabel('x1', fontsize=18)
ax.set_ylabel('x2', fontsize=18)
ax.set_zlabel('y', fontsize=18)
plt.show()

