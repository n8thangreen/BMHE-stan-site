#source: Cooper et al (2004)

working.dir <- paste(getwd(),"/",sep="")
setwd(working.dir)
source("../Utils.R")

# Defines the data
# Number of interventions (t=0: control; t=1: prophylactic use of Neuramidase Inhibitors (NI) 
T <- 2 					

# Evidence synthesis on effectiveness of NIs prophylaxis vs placebo
r0 <- r1 <- n0 <- n1 <- numeric()	# defines observed cases & sample sizes
r0 <- c(34,40,9,19,6,34)
r1 <- c(11,7,3,3,3,4)
n0 <- c(554,423,144,268,251,462)
n1 <- c(553,414,144,268,252,493)
S <- length(r0)				# number of relevant studies

# Evidence synthesis on incidence of influenza in healthy adults (under t=0)
x <- m <- numeric()			# defines observed values for baseline risk
x <- c(0,6,5,6,25,18,14,3,27)
m <- c(23,241,159,137,519,298,137,24,132)
H <- length(x)

# Data on costs
unit.cost.drug <- 2.4			# unit (daily) cost of NI
length.treat <- 6*7			# 6 weeks course of treatment
c.gp <- 19				# cost of GP visit to administer prophylactic NI
vat <- 1.175				# VAT
c.ni <- unit.cost.drug*length.treat*vat 

# Informative prior on cost of influenza 
mu.inf <- 16.78				# mean cost of influenza episode
sigma.inf <- 2.34			# sd cost of influenza episode
tau.inf <- 1/sigma.inf^2		# precision cost of influenza episode

# Informative prior on length of influenza episodes
m.l <- 8.2					# original value in the paper: 8.2 ##2.092
s.l <- sqrt(2)					# original value in the paper: sqrt(2) ##sqrt(.172)
mu.l <- lognPar(m.l,s.l)$mulog			# mean time to recovery (log scale)
sigma.l <- lognPar(m.l,s.l)$sigmalog		# sd time to recovery (log scale)
tau.l <- 1/sigma.l^2				# precision time to recovery (log scale)

# Parameters of unstructured effects
mean.alpha <- 0
sd.alpha <- sqrt(10)
prec.alpha <- 1/sd.alpha^2
mean.mu.delta <- 0
sd.mu.delta <- sqrt(10)
prec.mu.delta <- 1/sd.mu.delta^2
mean.mu.gamma <- 0
sd.mu.gamma <- 1000
prec.mu.gamma <- 1/sd.mu.gamma^2


# Prepares to launch JAGS
library(R2jags)
dataJags <- list("S","H","r0","r1","n0","n1","x","m","mu.inf","tau.inf","mu.l","tau.l","mean.alpha","prec.alpha",
	"mean.mu.delta","prec.mu.delta","mean.mu.gamma","prec.mu.gamma") 
filein <- "EvSynth.txt"
params <- c("p1","p2","rho","l","c.inf","alpha","delta","gamma")
inits <- function(){
	list(alpha=rnorm(S,0,1),delta=rnorm(S,0,1),mu.delta=rnorm(1),sigma.delta=runif(1),
		gamma=rnorm(H,0,1),mu.gamma=rnorm(1),sigma.gamma=runif(1),c.inf=rnorm(1),l=runif(1))
}

n.iter <- 20000
n.burnin <- 10000
n.thin <- floor((n.iter-n.burnin)/500)
es <- jags(data=dataJags,inits=inits,parameters.to.save=params,model.file=filein,
	n.chains=2, n.iter, n.burnin, n.thin, 
	DIC=TRUE, working.directory=working.dir, progress.bar="text")
print(es,digits=3,intervals=c(0.025, 0.975))
attach.bugs(es$BUGSoutput)

# Runs economic analysis 
# cost of treatment
c <- e <- matrix(NA,n.sims,T)
c[,1] <- (1-p1)*(c.gp) + p1*(c.gp+c.inf)
c[,2] <- (1-p2)*(c.gp+c.ni) + p2*(c.gp+c.ni+c.inf)
e[,1] <- -l*p1
e[,2] <- -l*p2

library(BCEA)
treats <- c("status quo","prophylaxis with NIs")
m <- bcea(e,c,ref=2,treats,Kmax=10000)

