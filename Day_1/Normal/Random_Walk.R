# A simple  Random Walk experiment


#Plot the position at N-th  step for nsim simulations  

plot_RW<-function(nsim,N){
for(i in 1:nsim){
plot(cumsum(runif(20,-1,1)),xlim=c(1,16),ylim=c(-6,6),type="l",col="red3",
     ylab="position",xlab="step")
par(new=TRUE)
}
}

# Plot the positon for each of 20 steps for 50 simulations
plot_RW(75,20)




# Plot a histogram of the final position
plot_hist<-function(nsim,N){
par(new=F)
pos <- replicate(nsim,sum(runif(N,-1,1)) ) 
hist(pos,main=paste("N=",nsim,sep=""),col = "green3")
}



par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(3, 3, 2, 2), oma = c(1, 1, 1, 1))
plot_hist(5,20)
plot_hist(10,20)
plot_hist(50,20)
plot_hist(100,20)
plot_hist(500,20)
plot_hist(2500,20)

