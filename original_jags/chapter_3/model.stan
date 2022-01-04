// Health economic evaluation --- example for Chapter 3
// Loosely based on Fox-Rushby & Cairns (2005)_ Economic Evaluation_ Open University Press

data {
int N;
int n;
int se;
int amb;
int N_studies;
real m_rho;
real tau_rho;
real a_gamma;
real b_gamma;
real m_amb;
real tau_amb;
real m_hosp;
real tau_hosp;

}

parameters {
  
}

transformed parameters {
  
}

model {
  // likelihood
	for (s in 1:N_studies) {
		se[s] ~ binomial(pi[1],n[s]);
		amb[s] ~ binomial(gamma,se[s]);
	}
	pi[1] ~ beta(a_pi,b_pi);
	pi[2] = pi[1]*rho;
	
	// priors
	rho ~ dnorm(m_rho,tau_rho);
	gamma ~ dbeta(a_gamma,b_gamma);
	c_amb ~ dlnorm(m_amb,tau_amb);
	c_hosp ~ dlnorm(m_hosp,tau_hosp);
}

generated quantities { 
  // Predictive distributions on the clinical outcomes
	for (t in 1:2) {
		SE[t] ~ binomial(pi[t],N)
		A[t] ~ binomial(gamma,SE[t])
		H[t] = SE[t] - A[t]
	}
}
