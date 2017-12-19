# Bayesian tutorial by Rafael S. de Souza - UNC, USA & COIN
#
# Partial example from Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2017, Cambridge Univ. Press
#
# Example of frequentist Poisson regression in R
# synthetic data
# 1 response (y) and 1 explanatory variables (x1) with one quadratic term

require(ggplot2)

set.seed(2016)                                        # set seed to replicate example
nobs <- 75                                          # number of observations
x1 <- runif(nobs,-3,6)                                     # random uniform variable
xb <- 0.25 + 0.75*x1 - 0.075*x1^2                                        # linear predictor
py <- rpois(nobs, exp(xb))                            # create py as adjusted random Poisson variate
poismod <- data.frame(py,x1)


# Fit a Poisson model to the same data
fit_pois <- glm(py ~ x1+I(x1^2),data=poismod,family = 'poisson')
pred_pois <- predict(fit_pois, type = 'response')
data_pois <- data.frame(x1,py,pred_pois)


summary(fit_pois) 


ggplot(data_pois, aes(x = x1, y = py)) +
  geom_segment(aes(xend = x1, yend = pred_pois)) +
  geom_point(color="purple",size=2) +
  geom_line(aes(x=x1,y=pred_pois))+
  theme(panel.background = element_rect(color = "black", fill = "gray85"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white")) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
