"""
Created by Emille E. O. Ishida on 5 Dec 2017.
For the 11th TRR33 Winter School
10 - 16 December 2107, Passo del Tonale - Italy
"""


from cosmoabc.priors import flat_prior
from cosmoabc.ABC_sampler import ABC
from cosmoabc.ABC_functions import read_input
from cosmoabc.plots import plot_2p
import numpy as np

from my_functions import my_distance, my_simulation, my_prior

#user input file
filename = 'user.input'

#read  user input
Parameters = read_input(filename)

# update dictionary of user input parameters 
# with customized functions
Parameters['distance_func'] = my_distance
Parameters['simulation_func'] = my_simulation
Parameters['prior']['mean']['func'] = my_prior


# in case you want to generate a pseudo-observed data set
Parameters['dataset1'] = my_simulation(Parameters['simulation_input'])

#calculate distance between 2 catalogues
dtemp = my_distance(Parameters['dataset1'], Parameters)

#determine dimension of distance output
Parameters['dist_dim'] = len(dtemp)

#initiate ABC sampler
sampler_ABC = ABC(params=Parameters)

#build first particle system
sys1 = sampler_ABC.BuildFirstPSystem()

#update particle system until convergence
sampler_ABC.fullABC()

#plot results
plot_2p(sampler_ABC.T, 'results.pdf' , Parameters)
