data {
  int<lower=1> N;
  vector[N] theta_hat;
  vector<lower=0>[N] V;
  real prior_mean;
  real<lower=0> prior_sd;
  real<lower=0> prior_tau_sd;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[N] z;
}
transformed parameters {
  vector[N] theta;
  theta = mu + tau * z;
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  z ~ normal(0, 1);
  for (i in 1:N)
    theta_hat[i] ~ normal(theta[i], 1 / sqrt(V[i]));}
