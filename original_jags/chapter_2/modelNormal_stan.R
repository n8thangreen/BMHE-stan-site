## ANALYSIS OF THE NORMAL MODEL - pages 69-73

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

###########################################
## 1. Uncentred model

stan_dat <- list("N","y","X","k") 			         	# define the data
filein <- paste0(working.dir,"modelNormal.stan")	# define the model file
params <- c("alpha","beta","sigma")		         		# define the parameters to be monitored
inits <- function(){						# create a function to randomly simulate 
	list(alpha = rnorm(1),
	     beta = rnorm(1),
	     lsigma = rnorm(1))	      # initial values for the required nodes
}

n.iter <- 10000							# define the number of iterations
n.burnin <- 9500						# define the burn-in
n.thin <- floor((n.iter - n.burnin)/500)	# compute the thinning to have 1000 iterations monitored

# Launch Stan and print the summary table
m.1 <-
  rstan::stan(data = stan_dat,
              init = inits,
              pars = params,
              file = filein,
              chains = 2,
              iter = n.iter,
              warmup = n.burnin,
              thin = n.thin)

print(m.1, digits=3, intervals=c(0.025, 0.975))


#############################################################
# 2. Uncentred model with thinning (run for longer chains)

n.iter <- 50000							# increase number of iterations to improve convergence
n.burnin <- 9500						# same burn-in
n.thin <- floor((n.iter - n.burnin)/500)				# compute thinning, this time > 1

m.2 <- 
  rstan::stan(data = stan_dat,
              init = inits,
              pars = params,
              file = filein,
              chains = 2,
              iter = n.iter,
              warmup = n.burnin,
              thin = n.thin)

print(m.2, digits=3, intervals=c(0.025, 0.975))


################################################
# 3. Centred model

X <- gestate - mean(gestate)					# create centred covariate to improve convergence
stan_dat <- list("N","y","X","k") 		# need to re-define data list

n.iter <- 10000							            # doesn't need so many iterations, now
n.burnin <- 9500
n.thin <- floor((n.iter - n.burnin)/500)	# thinning computed and set to 1 again

m.3 <- 
  rstan::stan(data = stan_dat,
              init = inits,
              pars = params,
              file = filein,
              chains = 2,
              iter = n.iter,
              warmup = n.burnin,
              thin = n.thin)

print(m.3,digits=3, intervals=c(0.025, 0.975))


#########
# plots #
#########

### 1. Uncentred model
rstan::traceplot(m.1, pars = c("alpha", "beta"), nrow = 2)

# ACF plot
rstan::stan_ac(m.1, pars = c("alpha", "beta"))

# Computes correlation for this model
cor(alpha.u,beta.u)
txt2 <- substitute("Scatterplot for "*alpha* " and "*beta* " - Uncentred model")
rstan::stan_scat(m.1, pars = c("alpha", "beta"))

### 2. Uncentred model with thinning
# Traceplot
rstan::traceplot(m.2, pars = c("alpha", "beta"), nrow = 2)

# ACF plot
rstan::stan_ac(m.2, pars = c("alpha", "beta"))

# Computes correlation for this model
cor(alpha.l,beta.l)
txt4 <- substitute("Scatterplot for "*alpha* " and "*beta* " - Uncentred model with thinning")
rstan::stan_scat(m.2, pars = c("alpha", "beta"))

### 3. Centred model 
rstan::traceplot(m.3, pars = c("alpha", "beta"), nrow = 2)

# ACF plot
rstan::stan_ac(m.3, pars = c("alpha", "beta"))

# Plots Scatterplot for correlation
cor(alpha.c,beta.c)
txt3 <- substitute("Scatterplot for "*alpha* " and "*beta* " - Centred model")
rstan::stan_scat(m.3, pars = c("alpha", "beta"))



