## HEALTH ECONOMIC EVALUATION --- EXAMPLE THROUGHOUT CHAPTER 3 AND CHAPTER 4

## Loosely based on Fox-Rushby & Cairns (2005). Economic Evaluation. Open University Press
## Parameters of the model
# pi = probability of side effects (specific to each treatment)
# rho = reduction in probability of side effects for t=2
# phi = probability of ambulatory visit, given side effects (constant across treatments)
# c.amb = cost of ambulatory visit after side effects
# c.hosp = cost of hospitalisation after side effects
# c.drug = cost of each treatment (known constant)

# define the working directory
working.dir <- paste0(getwd(),"/")
source("../Utils.R")

## Generates data for the number of side effects under t=0
N.studies <- 5
sample.size <- 32
prop.pi <- 0.24
n <- rpois(N.studies,sample.size)
se <- rbinom(N.studies,n,prop.pi)

## Generates the data for the number of patients with side effects requiring ambulatory care
prop.gamma <- 0.619
amb <- rbinom(N.studies,se,prop.gamma)

## Computes the hyperparameters for the informative prior distributions 
m.pi <- 0.5 
s.pi <- sqrt(0.125)
a.pi <- betaPar(m.pi,s.pi)$a
b.pi <- betaPar(m.pi,s.pi)$b

m.gamma <- 0.5
s.gamma <- sqrt(0.125)
a.gamma <- betaPar(m.gamma,s.gamma)$a
b.gamma <- betaPar(m.gamma,s.gamma)$b

m.rho <- 0.8
s.rho <- 0.2
tau.rho <- 1/s.rho^2

mu.amb <- 120
sd.amb <- 20
m.amb <- lognPar(mu.amb,sd.amb)$mulog
s.amb <- lognPar(mu.amb,sd.amb)$sigmalog
tau.amb <- 1/s.amb^2

mu.hosp <- 5500
sd.hosp <- 980
m.hosp <- lognPar(mu.hosp,sd.hosp)$mulog
s.hosp <- lognPar(mu.hosp,sd.hosp)$sigmalog
tau.hosp <- 1/s.hosp^2

c.drug <- c(110,520)

# Number of patients in the population
N <- 1000


library(R2jags)

dataJags <- list("se","amb","n","N.studies","a.pi","b.pi","a.gamma",
                 "b.gamma","m.amb","tau.amb","m.hosp",
                 "tau.hosp","m.rho","tau.rho","N") 
filein <- "model.txt"
params <- c("pi","gamma","c.amb","c.hosp","rho","SE","A","H")			# parameters list
inits <- function() {								# randomly initialise relevant variables
  SE=rbinom(2,N,.2)
  list(pi=c(runif(1),NA),
       gamma=runif(1),
       c.amb=rlnorm(1),
       c.hosp=rlnorm(1),
       rho=runif(1),
       SE=SE,
       A=rbinom(2,SE,.2))
}

n.iter <- 20000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/250)

chemo <-
  jags(dataJags,
       inits,
       params,
       model.file=filein,
       n.chains=2,
       n.iter,
       n.burnin,
       n.thin,
       DIC=TRUE,
       working.directory=working.dir,
       progress.bar="text")

## NB: If initial values for SE are not compatible with the rest of the model,
##     JAGS stops with an error.
##     The solution is to re-launch the whole model until it does work. 

print(chemo,digits=3,intervals=c(0.025, 0.975))
R2WinBUGS::attach.bugs(chemo$BUGSoutput)


#############################################################
## Cost-effectiveness analysis
# Defines the variables of cost and effectiveness
e <- c <- matrix(NA,chemo$BUGSoutput$n.sims,2)
e <- N - SE
for (t in 1:2) {
  c[,t] <- c.drug[t]*(N-SE[,t]) + (c.amb+c.drug[t])*A[,t] + (c.hosp+c.drug[t])*H[,t]
}

library(BCEA)

treats <- c("Old Chemotherapy","New Chemotherapy")
m <- bcea(e=e, c=c, ref=2, interventions=treats, Kmax=50000)
summary(m)


########################################################
## Compute EVPPI
# Analysis of Expected Value of Partial Information using
# Strong & Oakley's or Sadatsafavi et al.'s methods

inp <- CreateInputs(chemo)	# Create the inputs (parameters & matrix of MCMC simulations for all of them)
x.so <- evppi("rho", inp$mat, m, n.blocks=50)	# Use S&O method

# new BCEA version
# x.so <- evppi(he = m, param_idx = "rho", n.blocks = 50)	# Use S&O method


plot(x.so)
x.sal <- evppi("rho",inp$mat,m,n.seps=2)		# Use Sadatsafavi et al method
plot(x.sal)

# Compare the two methods
plot(x.so$k, x.so$evi,
     type="l",lwd=2,
     xlab="Willingness to pay", ylab="",
     main="Expected value of partial information")
points(x.so$k,x.so$evppi, type="l", col="blue")
points(x.sal$k,x.sal$evppi, type="l", col="red")
legend("topleft",
       c("EVPI","EVPPI for rho (S&O)","EVPPI for rho (S et al)"),
       cex = 0.7, bty = "n", lty = 1,
       col  =c("black","blue","red"))


# Analysis of Expected Value of Partial Information using 2stage approach
# Main parameters: 
# rho = reduction in side effects 
# gamma = chance of ambulatory care
rho.sim <- rho
gamma.sim <- gamma

detach.bugs()
inits <- NULL

###########################################################
# EVPPI for rho
# Replicates the C/E analysis for a fixed value of rho at each iteration

Ustar.phi <- matrix(NA,
                    nrow = chemo$BUGSoutput$n.sims,
                    ncol = length(m$k))
Rhat <- matrix(NA,
               nrow = chemo$BUGSoutput$n.sims,
               ncol = length(chemo$BUGSoutput$summary[,"Rhat"]) - 1)

for (i in 1:chemo$BUGSoutput$n.sims) {
  rho <- rho.sim[i]	
  dataJags <- list("se","amb","n","N.studies",
                   "a.pi","b.pi","a.gamma","b.gamma","m.amb",
                   "tau.amb","m.hosp","tau.hosp","N","rho") 
  filein <- "modelEVPPI_rho.txt"
  params <- c("pi","gamma","c.amb","c.hosp","SE","A","H")
  n.iter <- 10000
  n.burnin <- 9750
  n.thin <- floor((n.iter-n.burnin)/250)
  chemo.evppi <-
    jags(dataJags,
         inits,
         params,
         model.file=filein,
         n.chains=2,
         n.iter,
         n.burnin,
         n.thin,
         DIC=TRUE,
         working.directory=working.dir,
         progress.bar="text")
  
  Rhat[i,] <- chemo.evppi$BUGSoutput$summary[,"Rhat"]
  attach.bugs(chemo.evppi$BUGSoutput)
  
  # Performs the economic analysis given the simulations for all the
  # other parameters, conditionally on rho
  e.temp <- c.temp <- matrix(NA, nrow = chemo$BUGSoutput$n.sims, ncol = 2)
  e.temp <- N - SE
  for (t in 1:2) {
    c.temp[,t] <- c.drug[Rt]*(N-SE[,t]) + (c.amb+c.drug[t])*A[,t] + (c.hosp+c.drug[t])*H[,t]
  }
  
  ##TODO: what is Ktable?
  m.evppi <- bcea(e=e.temp,c=c.temp,ref=2,interventions=treats,Kmax=50000,Ktable=25000)
  
  # Computes the maximum expected utility for that iteration for each value of k
  Ustar.phi[i,] <- apply(m.evppi$Ustar,2,mean)
  rm(m.evppi)
  detach.bugs()
  print(i)
}

# Computes the average value of the maximum expected utility for each value of k
Vstar.phi <- apply(Ustar.phi,2,mean)
Umax <- apply(apply(m$U,c(2,3),mean),1,max)

# Computes the EVPPI for each value of k
EVPPI_rho <- Vstar.phi - Umax

# Plot the results
plot(x = m$k, y = m$evi,
     type = "l", xlab = "Willingness to pay", ylab = "") 
points(m$k,EVPPI_rho,type="l",lty=2)
text <- c("EVPI",expression(paste("EVPPI for ",rho, sep="")))
legend("topleft",text,lty=1:2,cex=0.9,box.lty=0)

# Checks for convergence in all (inner) simulations
par(mfrow=c(4,3))
for (i in 1:12) {
  plot(Rhat[,i],type="l",main=rownames(chemo.evppi$BUGSoutput$summary)[i],ylim=c(0.99,1.15))
  abline(h=1.1)
}

##########################################################
# EVPI for gamma
# Replicates the C/E analysis for a fixed value of rho at each iteration

Ustar.phi <- matrix(NA,
                    nrow = chemo$BUGSoutput$n.sims,
                    ncol = length(m$k))
Rhat <- matrix(NA,
               nrow = chemo$BUGSoutput$n.sims,
               ncol = length(chemo$BUGSoutput$summary[,"Rhat"]) - 1)
for (i in 1:chemo$BUGSoutput$n.sims) {
  gamma <- gamma.sim[i]
  
  dataJags <- list("se","amb","n","N.studies","a.pi","b.pi","m.rho","tau.rho","m.amb",
                   "tau.amb","m.hosp","tau.hosp","N","gamma") 
  filein <- "modelEVPPI_gamma.txt"
  params <- c("pi","rho","c.amb","c.hosp","SE","A","H")
  n.iter <- 10000
  n.burnin <- 9750
  n.thin <- floor((n.iter-n.burnin)/250)
  chemo.evppi <- jags(dataJags,
                      inits,
                      params,
                      model.file=filein,
                      n.chains=2,
                      n.iter,
                      n.burnin,
                      n.thin,
                      DIC=TRUE,
                      working.directory=working.dir,
                      progress.bar="text")
  Rhat[i,] <- chemo.evppi$BUGSoutput$summary[,"Rhat"]
  attach.bugs(chemo.evppi$BUGSoutput)
  
  # Performs the economic analysis given the simulations for all the
  # other parameters, conditionally on rho
  e.temp <- c.temp <- matrix(NA,chemo$BUGSoutput$n.sims,2)
  e.temp <- N - SE
  for (t in 1:2) {
    c.temp[,t] <- c.drug[t]*(N-SE[,t]) + (c.amb+c.drug[t])*A[,t] + (c.hosp+c.drug[t])*H[,t]
  }
  m.evppi <- bcea(e=e.temp,c=c.temp,ref=2,interventions=treats,Kmax=50000,Ktable=25000)
  
  # Computes the maximum expected utility for that iteration for each value of k
  Ustar.phi[i,] <- apply(m.evppi$Ustar,2,mean)
  rm(m.evppi); detach.bugs()
  print(i)
}

# Computes the average value of the maximum expected utility for each value of k
Vstar.phi <- apply(Ustar.phi,2,mean)
Umax <- apply(apply(m$U,c(2,3),mean),1,max)

# Computes the EVPPI for each value of k
EVPPI_gamma <- Vstar.phi - Umax

# Plot the results
plot(x = m$k, y = m$evi,type="l",xlab="Willingness to pay",ylab="") 
points(m$k,EVPPI_gamma,type="l",lty=2)
text <- c("EVPI",expression(paste("EVPPI for ",gamma,sep="")))
legend("topleft",text,lty=1:2,cex=0.9,box.lty=0)


# Checks for convergence in all (inner) simulations
par(mfrow=c(4,3))
for (i in 1:12) {
  plot(Rhat[,i],type="l",main=rownames(chemo.evppi$BUGSoutput$summary)[i],ylim=c(0.99,1.15))
  abline(h=1.1)
}



## Model averaging
# M1 = complete model
# M2 = rho =: 0.3
# M3 = rho =: 1

# analysis for M2
detach.bugs()
rho <- 1

dataJags <- list("se","n","amb","N.studies","a.pi","b.pi","a.gamma",
                 "b.gamma","m.amb","tau.amb","m.hosp","tau.hosp","N","rho") 
filein <- "modelEVPPI.txt"
params <- c("pi","gamma","c.amb","c.hosp","rho","SE","A","H")
inits <- function(){
  SE=rbinom(2,N,0.2)
  list(pi=c(runif(1),NA),
       gamma=runif(1),
       c.amb=rlnorm(1),
       c.hosp=rlnorm(1),
       SE=SE,
       A=rbinom(2,SE,0.2))
}

n.iter <- 20000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/250)
chemo2 <- jags(dataJags,
               inits,
               params,
               model.file=filein,
               n.chains=2,
               n.iter,
               n.burnin,
               n.thin,
               DIC=TRUE,
               working.directory=working.dir,
               progress.bar="text")
print(chemo2,digits=3,intervals=c(0.025, 0.975))
attach.bugs(chemo2$BUGSoutput,overwrite=TRUE)

## Cost-effectiveness analysis
# Defines the variables of cost and effectiveness
e <- c <- matrix(NA,
                 nrow = chemo$BUGSoutput$n.sims,
                 ncol = 2)
e <- N - SE
for (t in 1:2) {
  c[,t] <- c.drug[t]*(N-SE[,t]) + (c.amb+c.drug[t])*A[,t] + (c.hosp+c.drug[t])*H[,t]
}

m2 <- bcea(e=e,c=c,ref=2,interventions=treats,Kmax=50000)
fileout <- "Model_rho_fixed.ps"
postscript(fileout,paper="special",width=6,height=6,horizontal=FALSE)
plot(m2)
dev.off()

## Compares the models based on the DIC
d <- w <- numeric()
d[1] <- chemo$BUGSoutput$DIC
d[2] <- chemo2$BUGSoutput$DIC
dmin <- min(d)

# Computes the model weights
w <- exp(-0.5*(d-dmin))/sum(exp(-0.5*(d-dmin)))

# Computes the model average health economic analysis
SE <- w[1]*chemo$BUGSoutput$sims.list$SE + w[2]*chemo2$BUGSoutput$sims.list$SE
c.hosp <- w[1]*chemo$BUGSoutput$sims.list$c.hosp + w[2]*chemo2$BUGSoutput$sims.list$c.hosp
c.amb <- w[1]*chemo$BUGSoutput$sims.list$c.amb + w[2]*chemo2$BUGSoutput$sims.list$c.amb 
A <- w[1]*chemo$BUGSoutput$sims.list$A + w[2]*chemo2$BUGSoutput$sims.list$A
H <- w[1]*chemo$BUGSoutput$sims.list$H + w[2]*chemo2$BUGSoutput$sims.list$H

e <- c <- matrix(NA,chemo$BUGSoutput$n.sims,2)
e <- N - SE
for (t in 1:2) {
  c[,t] <- c.drug[t]*(N-SE[,t]) + (c.amb+c.drug[t])*A[,t] + (c.hosp+c.drug[t])*H[,t]
}
m.avg <- bcea(e=e,c=c,ref=2,interventions=treats,Kmax=50000)
plot(m.avg)

par(mfrow=c(2,2))
contour(m)
contour(m2)
contour(m.avg)

