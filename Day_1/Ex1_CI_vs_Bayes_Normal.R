# ESTEC - 18-20 December 2017
# by Rafael S. de Souza

# Bayesian vs Frequentist
# taken from the Python implementation of 
# http://jakevdp.github.io/blog/2014/06/12/frequentism-and-bayesianism-3-confidence-credibility/

require(mdatools)

# Example 1: check the outcomes from estimating the 
#            mean of 5 Gaussian distributed numbers given 
#            that we have 10^6 samples 

N = 5
Nsamp = 1000
sigma_x = 2

N = 5                                  # number of random numbers to be averaged
Nsamp = 1e6                            # number of sample for each random variable
sigma_x = 2                            # variance of all numbers

x <- matrix(rnorm(Nsamp*5, 0, 2), ncol = 5)          # generate 5 x Nsamp
mu_samp = rowMeans(x)                                # calculate the mean of each set


# Formula
sig_samp = sigma_x/sqrt(N)                           # Standard error of the mean

cat(sig_samp, "approximates", sd(mu_samp))           # 

true_B = 100
sigma_x = 10

set.seed(1)
D = rnorm(3,true_B, sigma_x)
print(D)

freq_CI_mu <- function(D,sigma,frac=0.95){
  
  Nsigma = sqrt(2)*erfinv(frac)
  mu = mean(D)
  sigma_mu = sigma/sqrt(length(D))
  return(c(mu - Nsigma * sigma_mu, mu + Nsigma * sigma_mu))
}
  
freq_CI_mu(D,sigma_x)

bayes_CR_mu <- function(D, sigma, frac=0.95){
    Nsigma = sqrt(2)*erfinv(frac)
    mu = mean(x)
    sigma_mu = sigma/sqrt(length(D))
    return(c(mu - Nsigma * sigma_mu, mu + Nsigma * sigma_mu))
}

type1 <- rnorm(50, mean = 0, sd = sqrt(2500))
type2 <- rnorm(50, mean = 0, sd = sqrt(1700))
var.test(type1, type2, alternative = "two.sided")


bayes_CR_mu <- function(D, sigma, frac=0.95){
Nsigma = sqrt(2)*erfinv(frac)
mu = mean(D)
sigma_mu = sigma/sqrt(length(D))
return(c(mu - Nsigma * sigma_mu, mu + Nsigma * sigma_mu))
}


bayes_CR_mu(D,sigma_x)
