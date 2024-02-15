data {
  int<lower=0> N;                         // Number of peptides
  int<lower=0> M;                         // Number of methods
  array[N, M] real<lower=0, upper = 1> y; // Observed errors
  vector<lower=0, upper = 1>[M] tau;      // Estimated method contribution
}

parameters {
  real<lower=0> lambda; // Dataset specific precision prior
  vector<lower=0>[N] t_p; // Precision prior
  vector<lower=0>[N] t; // Precision prior
  vector<lower=0, upper = 1>[N] mu; // Mean prior
}

model {
  // Priors
  lambda ~ exponential(.001);
  t_p ~ exponential(lambda);
  t ~ exponential(t_p);
  // mu ~ uniform(0, 1)
  for (i in 1:N) {
    // Transformation of y
    vector[M] log_y = log(rep_vector(1, M) -  to_vector(y[i]));
    // Jacobian correction
    target += -log_y;
    // Likelihood
    log_y ~ normal(log(mu[i]), inv(t[i] * tau));
  }
}
