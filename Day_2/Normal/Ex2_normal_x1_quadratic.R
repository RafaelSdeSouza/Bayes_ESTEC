# ADA8 – Astronomical Data Analysis Summer School
# Bayesian tutorial by Rafael S. de Souza - ELTE, Hungary & COIN
#
# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2016, Cambridge Univ. Press
#
# Example of frequentist linear regression in R
# synthetic data
# 1 response (y) and 1 explanatory variable (x1) with 1 quadratic term


set.seed(1056)                    # set seed to replicate example
nobs = 50                         # number of obs in model 
x <- runif(nobs,0,5)             # random uniform variable
mu <- 1 + 5 * x - 0.75 * x ^ 2  # linear predictor, xb
y <- rnorm(nobs, mu, sd=1)      # create y as adjusted random normal variate 


d <- data.frame(x,y)
plot(x,y)



fit <- lm(y ~ x1+I(x^2))          # Normal Fit 
summary(fit) 

d$predicted <- predict(fit)

ggplot(d, aes(x = x, y = y)) +
  geom_segment(aes(xend = x, yend = predicted)) +
  geom_point(color="purple",size=2) +
  geom_line(aes(x=x,y=predicted))+
  theme_xkcd() +
  theme(panel.background = element_rect(color = "black", fill = "gray85") )


