# JAGS LOGIT  Model for Red Spirals
############### Required packages
require(R2jags)
require(jagstools)
library(extrafont)
require(ggplot2)
require(ggthemes)
############### Simulation 



path_to_data = '../data/Red_spirals.csv'

# Read data
Red <- read.csv(path_to_data,header=T)
# Prepare data to JAGS
N <- nrow(Red)
x <- Red$fracdeV
y <- Red$type

gr <- data.frame(frac=x,type=factor(y))
ggplot(data=gr,aes(x=type,y=frac)) +
  geom_boxplot() +
  xlab("red galaxy")



# Prepare data for prediction 
M=500
xx = seq(from =  min(x), 
         to =  max(x), 
         length.out = M)


############### Construct data dictionary
X <- model.matrix(~ x)
XX <- model.matrix(~ xx)

K <- ncol(X)
logit_data <- list(Y  = y, # Response variable
                   X  = X,           # Predictors
                   N  = N,        # Sample size 
                   XX = XX, 
                   K = K,
                   M = M
)
############### JAGS code
LOGIT<-"model{

# Diffuse normal priors for predictors
for(i in 1:K){
beta[i]  ~ dnorm(0, 1e-4) 
}
# Likelihood function 

for (i in 1:N){  
Y[i] ~ dbern(p[i])
logit(p[i]) <- eta[i]
eta[i]      <- inprod(beta[], X[i,])

}
# Prediction for new data
for (j in 1:M){
etax[j]<-inprod(beta[], XX[j,])
logit(px[j]) <- etax[j]
Yx[j]~dbern(px[j])
}

}
"
#A function to generate initial values for mcmc
inits  <- function () {
  list(
    beta  = rnorm(ncol(X), 0, 0.1)  )  }

params <- c("beta","px")

# Run mcmc
LOGIT_fit<- jags(data     = logit_data,
               inits      = inits,
               parameters = params,
               model      = textConnection(LOGIT),
               n.thin     = 1,
               n.chains   = 3,
               n.burnin   = 3000,
               n.iter     = 6000)

############### Output on screen
print(LOGIT_fit,intervals=c(0.025, 0.975),justify = "left", digits=2)



# Plot
yx <- jagsresults(x=LOGIT_fit, params=c('px'))
gdata <- data.frame(x =xx, mean = yx[,"mean"],lwr1=yx[,"25%"],lwr2=yx[,"2.5%"],upr1=yx[,"75%"],upr2=yx[,"97.5%"])
# Bin data for visualization
binx<-0.05
t.breaks <-cut(x, seq(0,0.5, by=binx))
means <-tapply(y, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(y, t.breaks, semean)

gbin<-data.frame(x=seq(binx,max(x), by=binx),y=means)




logitmod <-data.frame(y, x) 


ggplot(logitmod,aes(x=x,y=y))+ 
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr1, ymax=upr1,y=NULL), alpha=0.95, fill=c("#969696"),show.legend=FALSE) +
  geom_ribbon(data=gdata,aes(x=xx,ymin=lwr2, ymax=upr2,y=NULL), alpha=0.65, fill = c("#d9d9d9"),show.legend=FALSE) +
  geom_line(data=gdata,aes(x=xx,y=mean),colour="gray25",linetype="dashed",size=1,show.legend=FALSE)+
  geom_point(aes(x=x,y=y),size=2.5,data=gbin,colour="#de2d26")+
  geom_errorbar(data=gbin,aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.85,
                colour="#de2d26",width=0.005)+
  theme_bw() +
#  ylab(expression(N[red]/N[Spirals]))+
  ylab(expression(p[red]))+
  xlab(expression(FracDev[r]))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 25,family="serif"))+
  coord_cartesian(ylim=c(0,0.4),xlim=c(0,0.5))

quartz.save(type = 'pdf', file = 'LOGIT_red_spirals.pdf',width = 10, height = 9)



