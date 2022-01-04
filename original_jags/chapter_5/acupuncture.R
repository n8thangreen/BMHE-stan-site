## Source: Wonderling et al (2004): http://www.bmj.com/highwire/filestream/347722/field_highwire_article_pdf/0.pdf

working.dir <- paste(getwd(),"/",sep="")
setwd(working.dir)

# Load and reformat the data on effects & costs for the two interventions
data <- read.csv("dataRCTacupuncture.csv",sep=";")
ctrl <- data.frame(pts=data[1:119,1],e=data[1:119,2],c=data[1:119,3])
trts <- data.frame(pts=data[1:136,4],e=data[1:136,5],c=data[1:136,6])
## e = QALYs 
## c = total cost


########################################################################################
## Model for cost-effectiveness 
# 1. logit/log normal model
# Transforms QALY values using logit transformation
e1 <- log(ctrl$e/(1-ctrl$e))
e2 <- log(trts$e/(1-trts$e))

# Transforms cost values using log transformation 
add <- 1		# Need to add a small value to avoid problems with 0s
c1 <- log(ctrl$c+add)
c2 <- log(trts$c+add)

# Sample sizes
n <- c(length(c1),length(c2))

# Run JAGS
library(R2jags)
dataJags <- list("n","c1","c2","e1","e2") 
filein <- "actptRCT.txt"
params <- c("mu.c","mu.e","beta","sigma.e","sigma.c")
inits <- function(){
	list(
		mu.c=rnorm(2,0,1),mu.e=rnorm(2,0,1),beta=runif(2,-.5,.5),
		logsigma.c=runif(2,0,1),logsigma.e=runif(2,0,1)
	)
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter-n.burnin)/500)
acupt <- jags(dataJags, inits, params, model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(acupt,digits=3,intervals=c(0.025, 0.975))
attach.bugs(acupt$BUGSoutput)

# Run the health economic analysis
# Cost-effectiveness parameters on the normal scale
# Numerical estimation of the mean for effectiveness on the natural scale
MC.sims <- 10000
logit.scale <- e.sim <- array(NA,c(MC.sims,n.sims,2))
for (t in 1:2) {
	for (i in 1:n.sims) {
		logit.scale[,i,t] <- rnorm(MC.sims,mu.e[i,t],sigma.e[i,t])
	# Samples of simulations from the naturale scale of the effectiveness measure
		e.sim[,i,t] <- exp(logit.scale[,i,t])/(1+exp(logit.scale[,i,t]))
	}
}

m.c <- m.e <- matrix(NA,acupt$BUGSoutput$n.sims,2)
for (t in 1:2) {
	m.c[,t] <- exp(mu.c[,t]+.5*(sigma.c[,t])^2)	# inverse lognormal transformation
	m.e[,t] <- apply(e.sim[,,t],2,mean)	 	# inverse logit transformation
}	
rm(e.sim,logit.scale)					# remove large objects from R workspace

ints <- c("Treatment as usual","Acupuncture")	# defines labels
e <- cbind(m.e[,1],m.e[,2])
c <- cbind(m.c[,1],m.c[,2])

library(BCEA)
m.nn <- bcea(e,c,ref=2,interventions=ints)

source("../Utils.R")
tab <- rbind(stats(m.c[,1])[c(-4)],stats(m.c[,2])[c(-4)],stats(m.e[,1])[c(-4)],stats(m.e[,2])[c(-4)])
rownames(tab) <- c("mean cost (ctrl)","mean cost (trts)","mean QALYs (ctrl)","mean QALYs (trts)")
colnames(tab) <- c("mean","sd","2.5%","97.5%")

# Saves the simulations in the object m1
m1 <- list(e=e,c=c)

########################################################################################
# 2. logit/normal-gamma model
# Transforms QALY values using logit transformation
e1 <- log(ctrl$e/(1-ctrl$e))
e2 <- log(trts$e/(1-trts$e))

# Transforms cost values using log transformation 
add <- 1		# Need to add a small value to avoid problems with 0s
c1 <- ctrl$c+add
c2 <- trts$c+add

# Sample sizes
n <- c(length(c1),length(c2))

# Extremes of mean costs
low <- 0
upp <- 2000

# Run JAGS
library(R2jags)
dataJags <- list("n","c1","c2","e1","e2","low","upp") 
filein <- "actptRCT_gamma.txt"
params <- c("mu.c","mu.e","beta","sigma.e")
inits <- function(){
	list(
		mu.c=runif(2,0,30),mu.e=rnorm(2,0,1),logsigma.e=runif(2,0,1),beta=runif(2,-.5,.5),eta=runif(2,0,1)
	)
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter-n.burnin)/500)
acupt_gam <- jags(dataJags, inits, params, model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(acupt_gam,digits=3,intervals=c(0.025, 0.975))
attach.bugs(acupt_gam$BUGSoutput)

# Run the health economic analysis
# Numerical estimation of the mean for effectiveness on the natural scale
MC.sims <- 10000
logit.scale <- e.sim <- array(NA,c(MC.sims,n.sims,2))
for (t in 1:2) {
	for (i in 1:n.sims) {
		logit.scale[,i,t] <- rnorm(MC.sims,mu.e[i,t],sigma.e[i,t])
	# Samples of simulations from the naturale scale of the effectiveness measure
		e.sim[,i,t] <- exp(logit.scale[,i,t])/(1+exp(logit.scale[,i,t]))
	}
}

# Cost-effectiveness parameters on the normal scale
ints <- c("Treatment as usual","Acupuncture")	# defines labels
m.c <- m.e <- matrix(NA,acupt_gam$BUGSoutput$n.sims,2)
for (t in 1:2) {
	m.c[,t] <- mu.c[,t]				# costs are already on the normal scale
	m.e[,t] <- apply(e.sim[,,t],2,mean)	 	# inverse logit transformation
}	
rm(e.sim,logit.scale)					# remove large objects from R workspace

tab <- rbind(stats(m.c[,1])[c(-4)],stats(m.c[,2])[c(-4)],stats(m.e[,1])[c(-4)],stats(m.e[,2])[c(-4)])
rownames(tab) <- c("mean cost (ctrl)","mean cost (trts)","mean QALYs (ctrl)","mean QALYs (trts)")
colnames(tab) <- c("mean","sd","2.5%","97.5%")

e <- cbind(m.e[,1],m.e[,2])
c <- cbind(m.c[,1],m.c[,2])

library(BCEA)
m.ng <- bcea(e,c,ref=2,interventions=ints)

# Saves the simulations in the object m2
m2 <- list(e=e,c=c)


########################################################################################
# 3. logit/normal-lognormal model
# Transforms QALY values using logit transformation
e1 <- log(ctrl$e/(1-ctrl$e))
e2 <- log(trts$e/(1-trts$e))

# Transforms cost values using log transformation 
add <- 1		# Need to add a small value to avoid problems with 0s
c1 <- ctrl$c+add
c2 <- trts$c+add

# Sample sizes
n <- c(length(c1),length(c2))

# Extremes of mean costs
low <- 0
upp <- 10000

# Run JAGS
library(R2jags)
dataJags <- list("n","c1","c2","e1","e2","low","upp") 
filein <- "actptRCT_logN.txt"
params <- c("mu.c","mu.e","beta","sigma.e")
inits <- function(){
	list(
		mu.c=runif(2,low,upp),mu.e=rnorm(2,0,1),beta=runif(2,-.5,.5),
		logsigma.c=runif(2,0,1),logsigma.e=runif(2,0,1)
	)
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter-n.burnin)/500)
acupt_ln <- jags(dataJags, inits, params, model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(acupt_ln,digits=3,intervals=c(0.025, 0.975))
attach.bugs(acupt_ln$BUGSoutput)

# Run the health economic analysis
# Numerical estimation of the mean for effectiveness on the natural scale
MC.sims <- 10000
logit.scale <- e.sim <- array(NA,c(MC.sims,n.sims,2))
for (t in 1:2) {
	for (i in 1:n.sims) {
		logit.scale[,i,t] <- rnorm(MC.sims,mu.e[i,t],sigma.e[i,t])
	# Samples of simulations from the naturale scale of the effectiveness measure
		e.sim[,i,t] <- exp(logit.scale[,i,t])/(1+exp(logit.scale[,i,t]))
	}
}

# Cost-effectiveness parameters on the normal scale
ints <- c("Treatment as usual","Acupuncture")	# defines labels
m.c <- m.e <- matrix(NA,acupt_ln$BUGSoutput$n.sims,2)
for (t in 1:2) {
	m.c[,t] <- mu.c[,t]				# costs are already on the normal scale
	m.e[,t] <- apply(e.sim[,,t],2,mean)	 	# inverse logit transformation
}	
rm(e.sim,logit.scale)					# remove large objects from R workspace

tab <- rbind(stats(m.c[,1])[c(-4)],stats(m.c[,2])[c(-4)],stats(m.e[,1])[c(-4)],stats(m.e[,2])[c(-4)])
rownames(tab) <- c("mean cost (ctrl)","mean cost (trts)","mean QALYs (ctrl)","mean QALYs (trts)")
colnames(tab) <- c("mean","sd","2.5%","97.5%")

e <- cbind(m.e[,1],m.e[,2])
c <- cbind(m.c[,1],m.c[,2])

library(BCEA)
m.nln <- bcea(e,c,ref=2,interventions=ints)

# Saves the simulations in the object m3
m3 <- list(e=e,c=c)


# Saves all JAGS objects
save(acupt,acupt_gam,acupt_ln,file="AcupunctureSimulations.Rdata")


