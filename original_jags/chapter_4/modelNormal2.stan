//

data {
}
  
model {
	for (i in 1:N) {
	    y[i] ~ dnorm(mu[i],tau)
	    mu[i] = alpha + beta*X[i]
	}
	alpha ~ dnorm(0,0.00001)
	beta ~ dnorm(0,0.00001)

	lsigma ~ dunif(-k,k)
	sigma = exp(lsigma)
	tau = pow(sigma,-2)

	y.star ~ dnorm(mu.star,tau)
	mu.star = alpha + beta*X.star
}

