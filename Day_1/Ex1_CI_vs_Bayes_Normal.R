# ESTEC - 18-20 December 2017
# by Rafael S. de Souza
# documentation by Emille E. O. Ishida

# Bayesian vs Frequentist
# taken from the Python implementation of 
# http://jakevdp.github.io/blog/2014/06/12/frequentism-and-bayesianism-3-confidence-credibility/

require(mdatools)
source("/Users/Rafael/Documents/GitHub/Bayes_ESTEC/auxiliar_functions/jagsresults.R")


# Example 1: check the outcomes from estimating the 
#            mean of 5 Gaussian distributed numbers given 
#            that we have 10^6 samples 


N = 5                                  # number of random numbers to be averaged
Nsamp = 1e5                          # number of sample for each random variable
sigma_x = 2                            # variance of all numbers


####  Frequentist approach
x <- matrix(rnorm(Nsamp*5, 0, 2), ncol = 5)          # generate Nsamp sets of 5 numbers each
mu_samp = rowMeans(x)                                # calculate the mean of each set (row)


# Formula
sig_samp = sigma_x/sqrt(N)                           # Standard error of the mean - given population value
                                                     # this expression is supposed to hold for large enough numbers

cat(sig_samp, "approximates", sd(mu_samp))           # compare expected  and calculated standard error of the mean


# Example on how to generate confidence intervals

# generate the data
true_B = 100                                        # mean                      
sigma_x = 10                                        # standard deviation

set.seed(1)
D = rnorm(3,true_B, sigma_x)                        # generate 3 data points
print(D)

# Create a function to calculate confidence intervals
freq_CI_mu <- function(D, sigma, frac=0.95){
  "Compute the confidence interval of the mean"
  
  Nsigma = sqrt(2)*erfinv(frac)
  mu = mean(D)
  sigma_mu = sigma/sqrt(length(D))
  
  return(c(mu - Nsigma * sigma_mu, mu + Nsigma * sigma_mu))
}
  
freq_CI_mu(D,sigma_x)                               # confidence intervals of 3 points

####   Bayesian approach

bayes_CR_mu <- function(D, sigma, frac=0.95){
    "Bayesian credbility interval"
  
    Nsigma = sqrt(2)*erfinv(frac)
    mu = mean(D)
    sigma_mu = sigma/sqrt(length(D))
    
    return(c(mu - Nsigma * sigma_mu, mu + Nsigma * sigma_mu))
}


bayes_CR_mu(D,sigma_x)                              # credible interval of 3 points
