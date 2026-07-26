data {
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;
  matrix[J, K] B;

  int<lower=0> d_A_L[I, J];
  int<lower=0> d_A_H[I, J];
  int<lower=0> d_B_L[I, J];
  int<lower=0> d_B_H[I, J];
  int<lower=0> d_AB_LL[I, J];
  int<lower=0> d_AB_LH[I, J];
  int<lower=0> d_AB_HL[I, J];
  int<lower=0> d_AB_HH[I, J];

  matrix<lower=0>[I, J] exposure_A_L;
  matrix<lower=0>[I, J] exposure_A_H;
  matrix<lower=0>[I, J] exposure_B_L;
  matrix<lower=0>[I, J] exposure_B_H;
  matrix<lower=0>[I, J] exposure_AB_LL;
  matrix<lower=0>[I, J] exposure_AB_LH;
  matrix<lower=0>[I, J] exposure_AB_HL;
  matrix<lower=0>[I, J] exposure_AB_HH;

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
  real<lower=0> gamma_sd;
  real<lower=0> smooth_gamma_sd;
}

parameters {
  vector[K] beta_A_L;
  vector[K] beta_A_H;
  vector[K] beta_B_L;
  vector[K] beta_B_H;
  vector[K] beta_delta;
  matrix[K, 4] beta_gamma_raw;

  vector[I] z_speed;
  vector[I] z_asymmetry;
  vector[I] z_capacity;

  real<lower=0> sigma_speed;
  real<lower=0> sigma_asymmetry;
  real<lower=0> sigma_capacity;
  real<lower=0> sigma_smooth_A;
  real<lower=0> sigma_smooth_B;
  real<lower=0> sigma_smooth_delta;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_smooth_gamma;
}

transformed parameters {
  vector[I] speed = speed_sd * sigma_speed * z_speed;
  vector[I] asymmetry = asymmetry_sd * sigma_asymmetry * z_asymmetry;
  vector[I] capacity_shift = capacity_sd * sigma_capacity * z_capacity;

  vector[J] population_A_L = B * beta_A_L;
  vector[J] population_A_H = B * beta_A_H;
  vector[J] population_B_L = B * beta_B_L;
  vector[J] population_B_H = B * beta_B_H;
  vector[J] population_delta = B * beta_delta;
  matrix[K, 4] beta_gamma;
  matrix[J, 4] population_gamma;

  matrix[I, J] eta_A_L;
  matrix[I, J] eta_A_H;
  matrix[I, J] eta_B_L;
  matrix[I, J] eta_B_H;
  matrix[I, J] eta_AB_LL;
  matrix[I, J] eta_AB_LH;
  matrix[I, J] eta_AB_HL;
  matrix[I, J] eta_AB_HH;

  for (k in 1:K) {
    beta_gamma[k, 1] = beta_gamma_raw[k, 1] -
      (beta_gamma_raw[k, 1] + beta_gamma_raw[k, 2] +
       beta_gamma_raw[k, 3] + beta_gamma_raw[k, 4]) / 4;
    beta_gamma[k, 2] = beta_gamma_raw[k, 2] -
      (beta_gamma_raw[k, 1] + beta_gamma_raw[k, 2] +
       beta_gamma_raw[k, 3] + beta_gamma_raw[k, 4]) / 4;
    beta_gamma[k, 3] = beta_gamma_raw[k, 3] -
      (beta_gamma_raw[k, 1] + beta_gamma_raw[k, 2] +
       beta_gamma_raw[k, 3] + beta_gamma_raw[k, 4]) / 4;
    beta_gamma[k, 4] = beta_gamma_raw[k, 4] -
      (beta_gamma_raw[k, 1] + beta_gamma_raw[k, 2] +
       beta_gamma_raw[k, 3] + beta_gamma_raw[k, 4]) / 4;
  }
  for (j in 1:J) {
    for (c in 1:4) {
      population_gamma[j, c] = 0;
      for (k in 1:K)
        population_gamma[j, c] += B[j, k] * beta_gamma[k, c];
    }
  }

  for (i in 1:I) {
    for (j in 1:J) {
      eta_A_L[i, j] = population_A_L[j] + speed[i] + 0.5 * asymmetry[i];
      eta_A_H[i, j] = population_A_H[j] + speed[i] + 0.5 * asymmetry[i];
      eta_B_L[i, j] = population_B_L[j] + speed[i] - 0.5 * asymmetry[i];
      eta_B_H[i, j] = population_B_H[j] + speed[i] - 0.5 * asymmetry[i];

      eta_AB_LL[i, j] = log_sum_exp(eta_A_L[i, j], eta_B_L[i, j]) +
        population_delta[j] + capacity_shift[i] + population_gamma[j, 1];
      eta_AB_LH[i, j] = log_sum_exp(eta_A_L[i, j], eta_B_H[i, j]) +
        population_delta[j] + capacity_shift[i] + population_gamma[j, 2];
      eta_AB_HL[i, j] = log_sum_exp(eta_A_H[i, j], eta_B_L[i, j]) +
        population_delta[j] + capacity_shift[i] + population_gamma[j, 3];
      eta_AB_HH[i, j] = log_sum_exp(eta_A_H[i, j], eta_B_H[i, j]) +
        population_delta[j] + capacity_shift[i] + population_gamma[j, 4];
    }
  }
}

model {
  // Non-centred subject effects use seconds-calibrated log-hazard scales.
  z_speed ~ std_normal();
  z_asymmetry ~ std_normal();
  z_capacity ~ std_normal();
  sigma_speed ~ normal(0, 1);
  sigma_asymmetry ~ normal(0, 1);
  sigma_capacity ~ normal(0, 1);

  if (K == 1) {
    beta_A_L[1] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_A_H[1] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_B_L[1] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_B_H[1] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_delta[1] ~ normal(delta_mean, delta_sd);
  } else {
    beta_A_L[1:2] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_A_H[1:2] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_B_L[1:2] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_B_H[1:2] ~ normal(log_hazard_mean, log_hazard_sd);
    beta_delta[1:2] ~ normal(delta_mean, delta_sd);
  }
  if (K > 2) {
    for (k in 3:K) {
      beta_A_L[k] - 2 * beta_A_L[k - 1] + beta_A_L[k - 2] ~ normal(0, sigma_smooth_A);
      beta_A_H[k] - 2 * beta_A_H[k - 1] + beta_A_H[k - 2] ~ normal(0, sigma_smooth_A);
      beta_B_L[k] - 2 * beta_B_L[k - 1] + beta_B_L[k - 2] ~ normal(0, sigma_smooth_B);
      beta_B_H[k] - 2 * beta_B_H[k - 1] + beta_B_H[k - 2] ~ normal(0, sigma_smooth_B);
      beta_delta[k] - 2 * beta_delta[k - 1] + beta_delta[k - 2] ~ normal(0, sigma_smooth_delta);
      for (c in 1:4)
        beta_gamma_raw[k, c] - 2 * beta_gamma_raw[k - 1, c] +
          beta_gamma_raw[k - 2, c] ~ normal(0, sigma_smooth_gamma);
    }
  }
  sigma_smooth_A ~ normal(0, smooth_A_sd);
  sigma_smooth_B ~ normal(0, smooth_B_sd);
  sigma_smooth_delta ~ normal(0, smooth_delta_sd);
  sigma_gamma ~ normal(0, gamma_sd);
  sigma_smooth_gamma ~ normal(0, smooth_gamma_sd);
  for (k in 1:K) for (c in 1:4)
    beta_gamma_raw[k, c] ~ normal(0, sigma_gamma);

  for (i in 1:I) {
    for (j in 1:J) {
      if (exposure_A_L[i, j] > 0) d_A_L[i, j] ~ poisson_log(log(exposure_A_L[i, j]) + eta_A_L[i, j]);
      if (exposure_A_H[i, j] > 0) d_A_H[i, j] ~ poisson_log(log(exposure_A_H[i, j]) + eta_A_H[i, j]);
      if (exposure_B_L[i, j] > 0) d_B_L[i, j] ~ poisson_log(log(exposure_B_L[i, j]) + eta_B_L[i, j]);
      if (exposure_B_H[i, j] > 0) d_B_H[i, j] ~ poisson_log(log(exposure_B_H[i, j]) + eta_B_H[i, j]);
      if (exposure_AB_LL[i, j] > 0) d_AB_LL[i, j] ~ poisson_log(log(exposure_AB_LL[i, j]) + eta_AB_LL[i, j]);
      if (exposure_AB_LH[i, j] > 0) d_AB_LH[i, j] ~ poisson_log(log(exposure_AB_LH[i, j]) + eta_AB_LH[i, j]);
      if (exposure_AB_HL[i, j] > 0) d_AB_HL[i, j] ~ poisson_log(log(exposure_AB_HL[i, j]) + eta_AB_HL[i, j]);
      if (exposure_AB_HH[i, j] > 0) d_AB_HH[i, j] ~ poisson_log(log(exposure_AB_HH[i, j]) + eta_AB_HH[i, j]);
    }
  }
}

generated quantities {
  // Pointwise log-likelihood over every subject-by-bin cell of the eight
  // salience series, in the order A_L, A_H, B_L, B_H, AB_LL, AB_LH, AB_HL,
  // AB_HH. Cells with no exposure contribute an exact zero and are dropped
  // before WAIC/LOO in R.
  vector[8 * I * J] log_lik;
  {
    int k = 0;
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_A_L[i, j] > 0 ?
        poisson_log_lpmf(d_A_L[i, j] | log(exposure_A_L[i, j]) + eta_A_L[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_A_H[i, j] > 0 ?
        poisson_log_lpmf(d_A_H[i, j] | log(exposure_A_H[i, j]) + eta_A_H[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_B_L[i, j] > 0 ?
        poisson_log_lpmf(d_B_L[i, j] | log(exposure_B_L[i, j]) + eta_B_L[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_B_H[i, j] > 0 ?
        poisson_log_lpmf(d_B_H[i, j] | log(exposure_B_H[i, j]) + eta_B_H[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_AB_LL[i, j] > 0 ?
        poisson_log_lpmf(d_AB_LL[i, j] | log(exposure_AB_LL[i, j]) + eta_AB_LL[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_AB_LH[i, j] > 0 ?
        poisson_log_lpmf(d_AB_LH[i, j] | log(exposure_AB_LH[i, j]) + eta_AB_LH[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_AB_HL[i, j] > 0 ?
        poisson_log_lpmf(d_AB_HL[i, j] | log(exposure_AB_HL[i, j]) + eta_AB_HL[i, j]) : 0;
    }
    for (i in 1:I) for (j in 1:J) {
      k += 1;
      log_lik[k] = exposure_AB_HH[i, j] > 0 ?
        poisson_log_lpmf(d_AB_HH[i, j] | log(exposure_AB_HH[i, j]) + eta_AB_HH[i, j]) : 0;
    }
  }
}
