#
# Bayesian tutorial by Rafael S. de Souza - UNC & COIN
#
# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2016, Cambridge Univ. Press
#
# Example of frequentist linear regression in R
# synthetic data
# 1 response (y) and 1 explanatory variable (x1)


set.seed(1056)                    # set seed to replicate example
nobs = 50                         # number of observations 
x <- runif(nobs,-2,5)                 # random uniform variable
mu <- -1 + 6*x  - 3*x^2                  # linear predictor
y <- rnorm(nobs, mu, sd = 1.5)     # create y as adjusted random normal variate 

d <- data.frame(x,y)
plot(x,y)

# Fit model
fit <- lm(y ~ x + I(x,2) , data=d)                 # Normal Fit of the synthetic data. 
summary(fit) 

d$predicted <- predict(fit)

ggplot(d, aes(x = x, y = y)) +
  geom_segment(aes(xend = x, yend = predicted)) +
  geom_point(color="purple",size=2) +
  geom_line(aes(x=x,y=predicted))+
  theme_xkcd() +
  theme(panel.background = element_rect(color = "black", fill = "gray85") )


