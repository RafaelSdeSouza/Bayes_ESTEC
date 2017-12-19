# ESTEC - 18-20 December 2017
# (c) 2017, Rafael S. de Souza 
#
# synthetic data
# 1 response (y) and 1 explanatory variable (x1) with 1 quadratic term

require(ggplot2)

set.seed(1056)                    # set seed to replicate example
nobs = 50                         # number of obs in model 
x <- runif(nobs,-0.25,1) 
phi <- 0.05
r <-1/phi
beta1 <- 1
beta2 <- 0.25
beta3 <- -2.25
xb <- beta1 + beta2* x - beta3 * x^2  
exb <- exp(xb)
y <- rgamma(nobs, shape = r, rate= r/exb)
d <- data.frame(x,y)
plot(x,y)



fit <- glm(y ~ x+I(x^2),family = Gamma(link = "log"))          # Normal Fit 
summary(fit) 

d$predicted <- predict(fit,type="response")

ggplot(d, aes(x = x, y = y)) +
  geom_segment(aes(xend = x, yend = predicted)) +
  geom_point(color="purple",size=2) +
  geom_line(aes(x=x,y=predicted))+
  theme(panel.background = element_rect(color = "black", fill = "gray85") )
