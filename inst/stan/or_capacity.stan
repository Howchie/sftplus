data {
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;
  matrix[J, K] B;

  int<lower=0> d_A[I, J];
  int<lower=0> d_B[I, J];
  int<lower=0> d_AB[I, J];
  matrix<lower=0>[I, J] exposure_A;
  matrix<lower=0>[I, J] exposure_B;
  matrix<lower=0>[I, J] exposure_AB;

  real log_hazard_mean;
  real<lower=0> log_hazard_sd;
  real delta_mean;
  real<lower=0> delta_sd;
  real<lower=0> speed_sd;
  real<lower=0> asymmetry_sd;
  real<lower=0> capacity_sd;
  real<lower=0> smooth_A_sd;
  real<lower=0> smooth_B_sd;
  real<lower=0> smooth_delta_sd;
}

parameters {
  vector[K] beta_A;
  vector[K] beta_B;
  vector[K] beta_delta;

  vector[I] z_speed;
  vector[I] z_asymmetry;
  vector[I] z_capacity;

  real<lower=0> sigma_speed;
  real<lower=0> sigma_asymmetry;
  real<lower=0> sigma_capacity;
  real<lower=0> sigma_smooth_A;
  real<lower=0> sigma_smooth_B;
  real<lower=0> sigma_smooth_delta;
}

transformed parameters {
  vector[I] speed = speed_sd * sigma_speed * z_speed;
  vector[I] asymmetry = asymmetry_sd * sigma_asymmetry * z_asymmetry;
  vector[I] capacity_shift = capacity_sd * sigma_capacity * z_capacity;

  vector[J] population_A = B * beta_A;
  vector[J] population_B = B * beta_B;
  vector[J] population_delta = B * beta_delta;

  matrix[I, J] eta_A;
  matrix[I, J] eta_B;
  matrix[I, J] eta_AB;

  for (i in 1:I) {
    for (j in 1:J) {
      eta_A[i, j] = population_A[j] + speed[i] + 0.5 * asymmetry[i];
      eta_B[i, j] = population_B[j] + speed[i] - 0.5 * asymmetry[i];
      eta_AB[i, j] = log_sum_exp(eta_A[i, j], eta_B[i, j]) +
        population_delta[j] + capacity_shift[i];
    }
  }
}

model {
  // Non-centred subject effects.  The half-normal scales are regularised by
  // the data-provided seconds-calibrated standard deviations.
  z_speed ~ std_normal();
  z_asymmetry ~ std_normal();
  z_capacity ~ std_normal();
  sigma_speed ~ normal(0, 1);
  sigma_asymmetry ~ normal(0, 1);
  sigma_capacity ~ normal(0, 1);

  // The first two basis coefficients anchor the broad log-hazard level and
  // put the population capacity process at the hazard-level UCIP centre.
  if (K == 1) {
    beta_A[1] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_B[1] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_delta[1] ~ normal(delta_mean, delta_sd);
  } else {
    beta_A[1:2] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_B[1:2] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_delta[1:2] ~ normal(delta_mean, delta_sd);
  }
  if (K > 2) {
    for (k in 3:K) {
      beta_A[k] - 2 * beta_A[k - 1] + beta_A[k - 2] ~ normal(0, sigma_smooth_A);
      beta_B[k] - 2 * beta_B[k - 1] + beta_B[k - 2] ~ normal(0, sigma_smooth_B);
      beta_delta[k] - 2 * beta_delta[k - 1] + beta_delta[k - 2] ~ normal(0, sigma_smooth_delta);
    }
  }
  sigma_smooth_A ~ normal(0, smooth_A_sd);
  sigma_smooth_B ~ normal(0, smooth_B_sd);
  sigma_smooth_delta ~ normal(0, smooth_delta_sd);

  // Aggregated piecewise-exponential likelihood.  Incorrect finite RTs are
  // already included in exposure and contribute no event count.
  for (i in 1:I) {
    for (j in 1:J) {
      if (exposure_A[i, j] > 0)
        d_A[i, j] ~ poisson_log(log(exposure_A[i, j]) + eta_A[i, j]);
      if (exposure_B[i, j] > 0)
        d_B[i, j] ~ poisson_log(log(exposure_B[i, j]) + eta_B[i, j]);
      if (exposure_AB[i, j] > 0)
        d_AB[i, j] ~ poisson_log(log(exposure_AB[i, j]) + eta_AB[i, j]);
    }
  }
}

generated quantities {
  // Pointwise log-likelihood over every subject-by-bin cell of each series, in
  // the order AB, A, B. Cells with no exposure contribute an exact zero and are
  // dropped before WAIC/LOO in R, so they affect neither lppd nor p_waic.
  vector[3 * I * J] log_lik;
  {
    int k = 0;
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_AB[i, j] > 0 ?
        poisson_log_lpmf(d_AB[i, j] | log(exposure_AB[i, j]) + eta_AB[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_A[i, j] > 0 ?
        poisson_log_lpmf(d_A[i, j] | log(exposure_A[i, j]) + eta_A[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_B[i, j] > 0 ?
        poisson_log_lpmf(d_B[i, j] | log(exposure_B[i, j]) + eta_B[i, j]) : 0;
    }
  }
}
