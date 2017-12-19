# Developed by Emille E. O. Ishida on 18 Dec 2017
# ESTEC Bayesian course, Netherlands
#
# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2017, Cambridge Univ. Press
#
# Example of frequentist linear regression in Python
# synthetic data
# 1 response (y) and 2 explanatory variable (x1, x2)

import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

from scipy.stats import uniform, norm

# set random seed
np.random.seed(1056) 

nobs = 150                                             # number of observations
x1 = uniform.rvs(loc=0, scale=2, size=nobs)            # random uniform variables
x2 = uniform.rvs(loc=0, scale=2, size=nobs)

xb = 2 + 3 * x1 - 4.5 * pow(x1, 2) + 1.5 * pow(x2, 2)    # linear predictor
y = norm.rvs(loc=xb, scale=0.5)                          # create y as adjusted random variate


# create data dictionary
my_data = {}
my_data['x1'] = x1
my_data['x12'] = pow(x1, 2)
my_data['x22'] = pow(x2, 2)
my_data['y'] = y

# fit a linear model using ordinary least squares
results = smf.ols(formula='y ~ 1 + x1 + x12 + x22', data=my_data).fit()

# output results
print(str(results.summary()))

# given fitted model, get predicted values
yfitted = results.predict()

# plot data + results
x1_axis = np.arange(0, 2, 0.01)
x2_axis = np.arange(0, 2, 0.01)
XX, YY = np.meshgrid(x1_axis, x2_axis)

exog = pd.core.frame.DataFrame({'x1':XX.ravel(), 'x12':pow(XX.ravel(),2), 'x22':pow(YY.ravel(),2)})
out = results.predict(exog)

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(XX, YY, out.reshape(XX.shape), color='orange', alpha=0.5)
ax.scatter(x1, x2, y, color='blue')
ax.set_xlabel('x1', fontsize=18)
ax.set_ylabel('x2', fontsize=18)
ax.set_zlabel('y', fontsize=18)
plt.show()

