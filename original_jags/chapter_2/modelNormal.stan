// normal model

data {
  int<lower=1> N;   // sample size
  real<lower=0> k;  // prior uniform bounds
  real y[N];
  real X[N];
}

parameters {
  real beta;
  real alpha;
  real lsigma;
}

transformed parameters {
  real mu[N];
	real<lower=0> sigma = exp(lsigma);
	
	for (i in 1:N) {
	  mu[i] = alpha + beta*X[i];
	}
}

model {
  // likelihood
	for (i in 1:N) {
	  target += normal_lpdf(y[i] | mu[i], sigma);
	}
	
	// priors
  //	alpha ~ normal(0, 0.00001);
  //	beta ~ normal(0, 0.00001);
	alpha ~ uniform(-k, k);
	beta ~ uniform(-k, k);
	lsigma ~ uniform(-k, k);
}
