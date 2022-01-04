## ANALYSIS OF THE NORMAL MODEL - pag 129-141

# define the working directory assuming - the user is already in ~/WebMaterial/ch2
working.dir <- paste0(getwd(), "/")

# Load data from STATA .dta format
data <- foreign::read.dta("phbirths.dta")
attach(data)

# Arrange the data and define relevant variables
y <- grams
N <- length(y)
X <- gestate
k <- 5000

library(R2jags)

## 4. Uncentred model with blocking (for chapter 4)
## NB Continues the analysis from chapter 2 (check modelNormal.R in the directory /ch2)
X <- gestate
m0 <- c(0,0)
prec <- 1/100000*matrix(c(1,0,0,1),nrow=2,ncol=2)
dataJags <- list("m0","prec","N","y","X","k") 
filein <- "modelNormalBlocking.txt"
params <- c("alpha","beta","sigma")
inits <- function(){
	list(lsigma=rnorm(1),coef=rnorm(2,0,1))
}

n.iter <- 10000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/500)
m.4 <- jags(dataJags,
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
print(m.4,digits=3,intervals=c(0.025, 0.975))


## 5. Predictive distribution (for chapter 4)
# 5a. Uncentered version with thinning
X.star <- 28
dataJags2 <- list("N","y","X","k","X.star") 
filein <- paste0(working.dir,"modelNormal2.txt")
params <- c("alpha","beta","sigma","y.star")
inits <- function(){
	list(alpha=rnorm(1),
	     beta=rnorm(1),
	     lsigma=rnorm(1),
	     y.star=runif(1))
}

n.iter <- 50000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/500)
m.5 <- jags(dataJags2,
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
print(m.5,digits=3,intervals=c(0.025, 0.975))

# 5b. Centred model with missing data
y <- c(y,NA)
N <- length(y)
X <- c(X,28)
# still call centred variable X so can still use same model file
#Xbar <- mean(X)
#X <- X - Xbar
dataJags3 <- list("N","y","X","k") 
filein <- "modelNormal.txt"
params3 <- c("alpha","beta","sigma","y[1116]")
inits <- function(){
list(alpha=rnorm(1),
     beta=rnorm(1),
     lsigma=rnorm(1),
     y=c(rep(NA,(N-1)),rnorm(1)))
}
n.iter <- 50000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/500)

m.6 <- jags(dataJags3, 
            inits, 
            params3, 
            model.file=filein,
            n.chains=2,
            n.iter=n.iter,
            n.burnin=n.burnin,
            n.thin=n.thin,
            DIC=TRUE)
print(m.6,digits=3,intervals=c(0.025, 0.975))
