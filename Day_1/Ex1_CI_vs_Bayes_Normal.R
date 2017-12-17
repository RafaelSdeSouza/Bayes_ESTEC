# Bayesian vs Frequentist
require(mdatools)

N = 5
Nsamp = 1000
sigma_x = 2


x <- matrix(rnorm(Nsamp*5, 0, 2), ncol = 5)
mu_samp = rowMeans(x)


# Formula
sig_samp = sigma_x/sqrt(N)

cat(sig_samp, "approximates", sd(mu_samp))



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
mu = mean(D)
sigma_mu = sigma/sqrt(length(D))
return(c(mu - Nsigma * sigma_mu, mu + Nsigma * sigma_mu))
}


bayes_CR_mu(D,sigma_x)



