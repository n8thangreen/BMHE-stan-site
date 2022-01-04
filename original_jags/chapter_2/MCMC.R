## PROGRAMS THE GIBBS SAMPLING - pag 64-67

# Define the number of observations & simulations
n <- 30

# Load the data
#y <- rnorm(n,4,10)	# This would generate random values
y <- c(1.2697,7.7637,2.2532,3.4557,4.1776,6.4320,-3.6623,7.7567,5.9032,7.2671,
	-2.3447,8.0160,3.5013,2.8495,0.6467,3.2371,5.8573,-3.3749,4.1507,4.3092,
	11.7327,2.6174,9.4942,-2.7639,-1.5859,3.6986,2.4544,-0.3294,0.2329,5.2846)

# Defines the hyper-parameters to build the full conditionals
ybar <- mean(y)
mu_0 <- 0
sigma2_0 <- 10000
alpha_0 <- 0.01
beta_0 <- 0.01

# Initialises the parameters
mu <- tau <- numeric()
sigma2 <- 1/tau

mu[1] <- rnorm(1,0,3)
tau[1] <- runif(1,0,3)
sigma2[1] <- 1/tau[1]

# Gibbs sampling (samples from the full conditionals)
nsim <- 1000
for (i in 2:nsim) {
	sigma_n <- sqrt(1/(1/sigma2_0 + n/sigma2[i-1]))
	mu_n <- (mu_0/sigma2_0 + n*ybar/sigma2[i-1])*sigma_n^2
	mu[i] <- rnorm(1,mu_n,sigma_n)

	alpha_n <- alpha_0+n/2
	beta_n <- beta_0 + sum((y-mu[i])^2)/2
	tau[i] <- rgamma(1,alpha_n,beta_n)
	sigma2[i] <- 1/tau[i]
}

## Creates a bivariate contour, using the package mvtnorm
require(mvtnorm) 
theta <- c(mean(mu),mean(sqrt(sigma2)))
sigma <- c(var(mu),var(sqrt(sigma2)))
rho <- cor(mu,sqrt(sigma2))
ins <- c(-10,10)
x1 <- seq(theta[1]-abs(ins[1]),theta[1]+abs(ins[1]),length.out=1000)
x2 <- seq(theta[2]-abs(ins[2]),theta[2]+abs(ins[2]),length.out=1000)
all <- expand.grid(x1, x2)
Sigma<-matrix(c(sigma[1], sigma[1]*sigma[2]*rho, sigma[1]*sigma[2]*rho, sigma[2]), nrow=2, ncol=2)
f.x <- dmvnorm(all, mean = theta, sigma = Sigma)
f.x2 <- matrix(f.x, nrow=length(x1), ncol=length(x2))

##Plots of the results at different simulation lengths 
lw <- c(1,1.6)

# First 10 iterations
plot(mu,sqrt(sigma2),col="white",pch=20,cex=.3,xlim=range(mu),ylim=range(sqrt(sigma2)),xlab=expression(mu),ylab=expression(sigma),
	main="After 10 iterations")
for (i in 1:(9)) {
	points(c(mu[i],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i]),t="l",col="grey60")
	points(c(mu[i+1],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i+1]),t="l",col="grey60")
}
text(jitter(mu[1:10]),jitter(sqrt(sigma2[1:10])),1:10,col="grey20")
contour(x = x1, y = x2, z = f.x2, nlevels=5,add=T,col="black",drawlabels=FALSE,lwd=lw[1])


# First 30 iterations
plot(mu,sqrt(sigma2),col="white",pch=20,cex=.3,xlim=range(mu),ylim=range(sqrt(sigma2)),xlab=expression(mu),ylab=expression(sigma),
	main="After 30 iterations")
for (i in 1:29) {
	points(c(mu[i],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i]),t="l",col="grey60")
	points(c(mu[i+1],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i+1]),t="l",col="grey60")
}
text(jitter(mu[1:30]),jitter(sqrt(sigma2[1:30])),1:30,col="grey20")
contour(x = x1, y = x2, z = f.x2, nlevels=5,add=T,col="black",drawlabels=FALSE,lwd=lw[1])


# First 100 iterations
plot(mu,sqrt(sigma2),col="white",pch=20,cex=.3,xlim=range(mu),ylim=range(sqrt(sigma2)),xlab=expression(mu),ylab=expression(sigma),
	main="After 100 iterations")
for (i in 1:99) {
	points(c(mu[i],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i]),t="l",col="grey60")
	points(c(mu[i+1],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i+1]),t="l",col="grey60")
}
text(jitter(mu[1:100]),jitter(sqrt(sigma2[1:100])),1:100,col="grey20")
contour(x = x1, y = x2, z = f.x2, nlevels=5,add=T,col="black",drawlabels=FALSE,lwd=lw[1])


# All 1000 iterations
plot(mu,sqrt(sigma2),col="grey60",pch=20,cex=.3,xlim=range(mu),ylim=range(sqrt(sigma2)),xlab=expression(mu),ylab=expression(sigma),
	main="After 1000 iterations")
for (i in 1:(nsim-1)) {
	points(c(mu[i],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i]),t="l",col="grey60")
	points(c(mu[i+1],mu[i+1]),c(sqrt(sigma2)[i],sqrt(sigma2)[i+1]),t="l",col="grey60")
}
text(jitter(mu),jitter(sqrt(sigma2)),1:nsim,col="grey28")
contour(x = x1, y = x2, z = f.x2, nlevels=5,add=T,col="black",drawlabels=FALSE,lwd=lw[2])

