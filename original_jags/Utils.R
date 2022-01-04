### Set of utility functions to perform Bayesian estimations and health economic analysis

############################################################
## Produce a traceplot to assess convergence of a MCMC run
## node is a string with the name of the node to be plotted
## model is the name of the object containing the MCMC simulations
## Copyright Gianluca Baio 2012
mytraceplot <- function(node,model=m,title="",lab=""){
	xlab <- "Iteration"
	cmd <- ifelse(class(model)=="rjags",mdl <- model$BUGSoutput,mdl <- model)
	col <- colors()
	if (mdl$n.chains==1) {
		plot(mdl$sims.array[,1,node],t="l",col=col[which(col=="blue")],xlab=xlab,
			ylab=lab,main=title)
	}
	if (mdl$n.chains==2) {
		cols <- c("blue","red")
		plot(mdl$sims.array[,1,node],t="l",col=col[which(col==cols[1])],xlab=xlab,
			ylab=lab,main=title,ylim=range(mdl$sims.array[,1:2,node]))
		points(mdl$sims.array[,2,node],t="l",col=col[which(col==cols[2])])
	}
	if (mdl$n.chains>2) {
		cols <- c("blue","red","green","magenta","orange","brown","azure")
		plot(mdl$sims.array[,1,node],t="l",col=col[which(col==cols[1])],xlab=xlab,
			ylab=lab,main=title,
			ylim=range(mdl$sims.array[,1:mdl$n.chains,node]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,i,node],t="l",col=col[which(col==cols[i])])
		}
	}
}



############################################################
## Compute the value of parameters (a,b) for a Beta distribution to have mean and sd (m,s)
## Copyright Gianluca Baio 2012
betaPar <- function(m,s){
	a <- m*( (m*(1-m)/s^2) -1 )
	b <- (1-m)*( (m*(1-m)/s^2) -1 )
	list(a=a,b=b)
}



############################################################
## Compute the value of parameters (mulog,sigmalog) for a logNormal distribution to have mean and sd (m,s)
## Copyright Gianluca Baio 2012
lognPar <- function(m,s) {
	s2 <- s^2
	mulog <- log(m) - .5*log(1+s2/m^2)
	s2log <- log(1+(s2/m^2))
	sigmalog <- sqrt(s2log)
	list(mulog=mulog,sigmalog=sigmalog)
}



############################################################
## Compute the parameters of a Beta distribution, given a prior guess for:
##  mode = the mode of the distribution
##  upp  = an upper bound value for the distribution
##  prob = the estimated probability that (theta <= upp)
## Based on "Bayesian ideas and data analysis", page 100. 
## Optimisation method to identify the values of a,b that give required conditions on the Beta distribution
## Copyright Gianluca Baio 2012
betaPar2 <- function(mode,upp,prob){
N <- 10000
b <- 1:N
a <- (1+mode*(b-2))/(1-mode)
sim <- qbeta(prob,a,b)
m <- ifelse(prob>=.5,max(which(sim>=upp)),min(which(sim>=upp)))
M <- ifelse(prob>=.5,min(which(sim<=upp)),max(which(sim<=upp)))

b <- min(m,M)+(b/N)
a <- (1+mode*(b-2))/(1-mode)
sim <- qbeta(prob,a,b)
m <- ifelse(prob>=.5,max(which(sim>=upp)),min(which(sim>=upp)))
M <- ifelse(prob>=.5,min(which(sim<=upp)),max(which(sim<=upp)))
a <- ifelse(m==M,a[m],mean(a[m],a[M]))
b <- ifelse(m==M,b[m],mean(b[m],b[M]))

step <- 0.001
theta <- seq(0,1,step)
density <- dbeta(theta,a,b)

norm.dens <- density/sum(density)
cdf <- cumsum(norm.dens)
M <- min(which(cdf>=.5))
m <- max(which(cdf<=.5))

theta.mode <- theta[which(density==max(density))]
theta.mean <- a/(a+b)
theta.median <- mean(theta[m],theta[M])
theta.sd <- sqrt((a*b)/(((a+b)^2)*(a+b+1)))
beta.params <- c(a,b,theta.mode,theta.mean,theta.median,theta.sd)
res1 <- beta.params[1]
res2 <- beta.params[2]
theta.mode <- beta.params[3]
theta.mean <- beta.params[4]
theta.median <- beta.params[5]
theta.sd <- beta.params[6]
list(
res1=res1,res2=res2,theta.mode=theta.mode,theta.mean=theta.mean,theta.median=theta.median,theta.sd=theta.sd)
}



############################################################
## Posterior summary
## Creates summary statistics a' la WinBUGS --- currently works only for vectors
## Copyright Gianluca Baio 2012
stats <- function(x){
c(mean(x),sd(x),quantile(x,.025),median(x),quantile(x,.975))
}


