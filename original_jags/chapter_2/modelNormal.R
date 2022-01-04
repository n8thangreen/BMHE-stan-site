## ANALYSIS OF THE NORMAL MODEL - pages 69-73

# define the working directory assuming - the user is already in ~/WebMaterial/ch2
working.dir <- paste0(getwd(), "/")

# Load data from .dta format
library(foreign)
data <- read.dta("phbirths.dta")
attach(data)

# Arrange the data and define relevant variables
y <- grams
N <- length(y)
X <- gestate
k <- 5000

## Run JAGS
library(R2jags)

###########################################
## 1. Uncentred model

dataJags <- list("N","y","X","k") 				# define the data
filein <- paste(working.dir,"modelNormal.txt",sep="")		# define the model file
params <- c("alpha","beta","sigma")				# define the parameters to be monitored
inits <- function(){						# create a function to randomly simulate 
	list(alpha=rnorm(1),beta=rnorm(1),lsigma=rnorm(1))	#   initial values for the required nodes
}

n.iter <- 10000							# define the number of iterations
n.burnin <- 9500						# define the burn-in
n.thin <- floor((n.iter-n.burnin)/500)				# compute the thinning to have 1000 iterations monitored

# Launch JAGS and print the summary table
m.1 <- jags(dataJags, inits, params, model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(m.1,digits=3,intervals=c(0.025, 0.975))


#############################################################
# 2. Uncentred model with thinning (run for longer chains)

n.iter <- 50000							# increase number of iterations to improve convergence
n.burnin <- 9500						# same burn-in
n.thin <- floor((n.iter-n.burnin)/500)				# compute thinning, this time > 1

# Launch JAGS and print the summary table
m.2 <- jags(dataJags, inits, params, model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(m.2,digits=3,intervals=c(0.025, 0.975))


################################################
# 3. Centred model

X <- gestate-mean(gestate)					# create centred covariate to improve convergence
dataJags <- list("N","y","X","k") 				# need to re-define data list

n.iter <- 10000							# doesn't need so many iterations, now
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/500)				# thinning computed and set to 1 again

# Launch JAGS and print the summary table
m.3 <- jags(dataJags, inits, params, model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(m.3,digits=3,intervals=c(0.025, 0.975))


## Make plots
source("../Utils.R")						# load utility functions

# 1. Uncentred model
# Traceplot
par(mfrow=c(2,1))
mytraceplot("alpha",title="Uncentred model",lab=expression(alpha),model=m.1)
mytraceplot("beta",lab=expression(beta),model=m.1)

# Makes ACF plot
alpha.u <- c(m.1$BUGSoutput$sims.array[,1,1],m.1$BUGSoutput$sims.array[,2,1])
beta.u <- c(m.1$BUGSoutput$sims.array[,1,2],m.1$BUGSoutput$sims.array[,2,2])
txt <- substitute("Autocorrelation function for "*alpha* " - Uncentred model")
acf(alpha.u,main=txt,ylab="")

# Computes correlation for this model
cor(alpha.u,beta.u)
txt2 <- substitute("Scatterplot for "*alpha* " and "*beta* " - Uncentred model")
plot(m.1$BUGSoutput$sims.list$alpha,m.1$BUGSoutput$sims.list$beta,xlab=expression(alpha),
	ylab=expression(beta),main=txt2)


# 2. Uncentred model with thinning
# Traceplot
par(mfrow=c(2,1))
mytraceplot("alpha",title="Uncentred model with thinning",lab=expression(alpha),model=m.2)
mytraceplot("beta",lab=expression(beta),model=m.2)

# Makes ACF plot
alpha.l <- c(m.2$BUGSoutput$sims.array[,1,1],m.2$BUGSoutput$sims.array[,2,1])
beta.l <- c(m.2$BUGSoutput$sims.array[,1,2],m.2$BUGSoutput$sims.array[,2,2])
txt <- substitute("Autocorrelation function for "*alpha* " - Uncentred model with thinning")
acf(alpha.l,main=txt,ylab="")

# Computes correlation for this model
cor(alpha.l,beta.l)
txt4 <- substitute("Scatterplot for "*alpha* " and "*beta* " - Uncentred model with thinning")
plot(m.2$BUGSoutput$sims.list$alpha,m.2$BUGSoutput$sims.list$beta,xlab=expression(alpha),	
	ylab=expression(beta),main=txt4)


# 3. Centred model 
par(mfrow=c(2,1))
mytraceplot("alpha",model=m.3,title="Centred model",lab=expression(alpha))
mytraceplot("beta",model=m.3,lab=expression(beta))

# Makes ACF plot
alpha.c <- c(m.3$BUGSoutput$sims.array[,1,1],m.3$BUGSoutput$sims.array[,2,1])
beta.c <- c(m.3$BUGSoutput$sims.array[,1,2],m.3$BUGSoutput$sims.array[,2,2])
txt <- substitute("Autocorrelation function for "*alpha* " - Centred model")
acf(alpha.c,main=txt,ylab="")

# Plots Scatterplot for correlation
cor(alpha.c,beta.c)
txt3 <- substitute("Scatterplot for "*alpha* " and "*beta* " - Centred model")
plot(m.3$BUGSoutput$sims.list$alpha,m.3$BUGSoutput$sims.list$beta,xlab=expression(alpha),
	ylab=expression(beta),main=txt3)



