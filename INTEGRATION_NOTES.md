# sft_plus integration notes

The base repository is the upstream `jhoupt/sft` source. Our enhancements are integrated into the package's canonical implementations in `R/sft_core.R` and `R/simulation.R`. The upstream files that would otherwise redefine those public functions have been removed or reduced to their non-overlapping functionality.

Integrated local functionality:

- Rcpp-backed NAH, NAK, and UCIP evaluators with an R fallback;
- robust OR, AND, STST, ID, and Altieri capacity functions, including bounds;
- restored UCIP, SIC, stochastic-dominance, MIC, and resilience test outputs;
- assessment helpers, CDF/SIC/assessment data builders, and smoothing helpers;
- LBA and OU simulation with R and optional Rcpp backends, logical gates, and
  capacity/SIC/CDF/assessment summaries.

Bayesian routes:

- `sictestBayes(method = "DP")` implements the discretized Dirichlet-process
  Monte Carlo route in `R/sic.R`, with explicit Monte Carlo controls. Each of the
  four cells carries its own Dirichlet over the shared pooled-RT bins and the SIC
  is formed as the contrast of their CDFs before classification. `tolSIC` is
  applied to every sign test, so the all-negative and all-positive classes remain
  reachable for noisy posterior draws;
- `siDominance(method = "dp")` provides directional posterior probabilities
  and posterior/prior event-probability ratios for each dominance relation. Rows
  are labelled in *survivor* terms: `S.hh > S.hl` means the HH survivor dominates,
  which is the `less` alternative of `ks.test` on the CDFs;
- `capacityGroup.bayes()` is a hierarchical Bayesian companion to the UCIP/Cz z
  test (the Bayesian analogue of `capacityGroup()`; the single-subject kernel is
  `ucip.bayes()`). It extracts the score numerator and information for each participant
  and partially pools the resulting score effects with a conjugate
  Normal/inverse-Gamma hierarchy. It does not pool raw RTs.

  The pooled effect defaults to `score_method = "score"`: theta = U / V, the Peto
  one-step estimate of the log hazard ratio for the UCIP contrast, with precision
  V. theta is stable in location as trials accumulate while V grows linearly with
  them, so a fixed prior width is a genuine effect-size prior. `"standardized"`
  reports the identical fit on the reference-Cz scale phi = theta * sqrt(V_ref);
  `"multiplicative"` is the separate eta = log(A/B) estimand. A population
  point-null Bayes factor (mu = 0, tau free) is reported via a Rao-Blackwellised
  Savage-Dickey ratio; it is not a test of "every participant is UCIP".

Later Bayesian additions:

- `resilience.test()` was fixed to use the published estimator. Resilience
  (Little, Eidels, Fific & Wang, 2015, Eq. 11) has the *same functional form* as
  the OR capacity coefficient, R(t) = H_AB(t) / (H_AY(t) + H_XB(t)), differing
  only in that the denominator conditions carry distractors; Houpt & Little
  (2016) accordingly build its test on the Houpt-Townsend (2012) weighted-logrank
  score statistic. The old implementation instead read the Nelson-Aalen contrast
  off at the terminal time, which is degenerate: for complete data the cumulative
  hazard at a sample's own maximum is exactly the harmonic number sum(1/k)
  whatever the distribution, so all three series coincide and the statistic is
  the same constant for every dataset. Measured on three datasets (true UCIP-OR,
  strongly super-capacity, strongly limited) the old z was -2.5220 in all three
  and rejected the true null at p = 0.012; the score statistic gives -0.03, 5.73,
  and -6.11.
- `resilience.bayes()` / `resilienceGroup.bayes()` are consequently `ucip.bayes()`
  / `capacityGroup.bayes()` applied to the resilience conditions: same score
  machinery, same three `score_method` scales (`"score"` default), same priors.
  There is no evaluation horizon and no `at` / `at_quantile` argument.
  `theta_hat * sqrt(V)` reproduces `resilience.test()` exactly.
- The hierarchical group models accept several conditions. `Condition = NULL` now
  fits every condition present, giving one mu_c per condition over a shared tau
  plus all pairwise `condition_contrasts`. Subjects are independent given
  (mu, tau), so the contrasts are between-condition population comparisons, not
  within-subject differences.
- `prior_sensitivity()` refits over a grid of alternative-prior widths and reports
  both null Bayes factors, making the Jeffreys-Lindley dependence explicit.
- `spike_slab()` gives each participant an exact-zero spike or the population
  slab, returning per-participant inclusion probabilities and Bayes factors. The
  inclusion rate w is estimated, so participants shrink toward one another.
- Race model inequality (`R/rmi.R`). Miller's (1982) bound is a different
  hypothesis from the UCIP capacity coefficient: C_OR(t) > 1 rejects one specific
  model (independent, unlimited-capacity parallel), whereas the race bound holds
  for *every* race model regardless of dependence or capacity, so super capacity
  does not imply a Miller violation and a limited-capacity race violates UCIP
  while satisfying Miller. Nothing in the file reuses the UCIP score.
  - `rmi.test()` reproduces the classical Ulrich-Miller-Schroeter procedure on
    the millisecond scale, for replications and for the figure.
  - `rmiGroup.bayes()` tests the same inequality on the *probability* scale:
    delta_ik = F_A(t_ik) + F_B(t_ik) - F_AB(t_ik), where t_ik is the u_k-th
    quantile of participant i's own pooled correct single-target RTs. Evaluating
    the inequality within a participant makes overall speed cancel exactly, and
    the estimand is invariant to any strictly increasing reparameterisation of
    time applied to all three of a participant's conditions, which is what
    licenses pooling on the shared u index. The classical procedure instead
    averages millisecond quantiles across participants, where a slow
    participant's 40 ms violation outweighs a fast participant's equally severe
    15 ms one and real violations can cancel across differing base times.
  - The hierarchy is the existing normal-normal engine with the percentile level
    playing the role the Condition factor plays elsewhere: mu_k per level over a
    shared tau. The global summary is min_k mu_k, whose P(< 0) is the
    multiplicity-honest counterpart of the classical family of per-percentile
    t-tests. It concerns the tested grid only.
  - Percentile levels where the bound is already vacuous are not fitted: for the
    Miller bound F_A + F_B reaches 1 near the median of the single-target
    distributions, so a default grid covers the fast half of the distribution,
    which is where violations occur.
  - `errors = "defective"` counts errors in the denominator (each cell CDF is
    P(correct, RT <= t), asymptote = accuracy) instead of discarding them.
    Prefer it when accuracy differs across cells, since discarding renormalises
    away the extra correct responses of a condition that is both faster and more
    accurate. The classical time-scale route is correct-only by construction.
  - `build_rmi_bound_df()` / `build_rmi_cdf_df()` / `build_rmi_violation_df()`
    return tidy data frames in the same spirit as `build_cdf_df()`; ggplot2 stays
    in Suggests. The bound builder is the AB-against-(A + B) figure. Its
    quantiles are averaged across participants ("Vincentized"), which *is*
    base-time dependent -- it is the picture, and the probability-scale violation
    curve is the inference.
- `sft_waic()` returns WAIC (and PSIS-LOO when the optional `loo` package is
  installed) from the pointwise Poisson log-likelihood emitted by the Stan
  models' generated quantities blocks. The criteria are cell-wise: predictive
  accuracy for a new time bin of an observed participant, not a new participant.

# Model 2: A hierarchical semiparametric hazard model

The 2016 Bayesian SFT paper (sftbayes/Van Zandt et al (2016).pdf, corresponding code in sftbayes/parametricSIC) established the basic route: approximate each condition’s hazard using piecewise-constant rates, smooth adjacent rates with a stochastic process, and transform posterior hazard draws into posterior capacity and SIC curves. Their model used an AR(1) process over log hazards and hierarchical priors over observers and conditions.

The version I would implement now differs in two ways:

1. use a common time grid and spline or random-walk smoothing rather than subject-condition-specific AR processes;
2. parameterise the dual-target hazard around the UCIP prediction, so the hierarchy directly pools **capacity deviations**, not merely generic hazard parameters.

The second change is the important one.

## 2.1 Piecewise-exponential likelihood

Choose common boundaries

[
0=s_0<s_1<\cdots<s_J,
]

and assume that the hazard for subject (i), condition (c), and interval (j) is constant:

[
h_{ic}(t)=h_{icj},
\qquad
s_{j-1}<t\le s_j.
]

For each subject-condition-bin, calculate:

[
d_{icj}
=======

\text{number of correct responses terminating in bin }j,
]

and

[
E_{icj}
=======

\text{total time at risk contributed in bin }j.
]

For trial (k) with observed time (t_{ick}), its exposure in bin (j) is

[
e_{ickj}
========

\max\left[
0,;
\min(t_{ick},s_j)-s_{j-1}
\right].
]

Then

[
E_{icj}=\sum_k e_{ickj}.
]

The likelihood is

[
d_{icj}\sim
\operatorname{Poisson}(E_{icj}h_{icj}).
]

Equivalently,

[
d_{icj}\sim
\operatorname{PoissonLog}
\left(
\log E_{icj}+\eta_{icj}
\right),
\qquad
\eta_{icj}=\log h_{icj}.
]

This is the aggregated form of the piecewise-exponential survival likelihood. Incorrect responses or omissions can contribute exposure without an event, mirroring the censoring treatment in the Nelson–Aalen implementation. The original Bayesian SFT model uses the same underlying piecewise-exponential construction.

## 2.2 Use one common grid

The 2016 model used quantiles separately within samples so each bin had similar numbers of observations. That worked computationally, but it means different observer-condition hazards are defined over different intervals, which is awkward for direct pointwise pooling. The paper also acknowledges that the quantile boundaries use the data twice and that the first and final bins can produce systematic distortion. 
For a group model, I would instead use:

* RT measured in seconds;
* one set of boundaries shared across all subjects and conditions;
* approximately 20 to 30 pooled quantile bins;
* smoothing evaluated using the actual bin midpoints and widths.

Pooled quantiles give narrow bins where the data are dense and broad bins in the tail. Because the boundaries are common, (h_{icj}) always refers to the same time interval.

The final model should include all observations, but capacity summaries should generally be restricted to a central region where all conditions are informed, perhaps the interval spanning the pooled 5th to 95th percentiles.

---

# 2.3 A generic hierarchy

The direct analogue of my earlier suggestion is

[
\eta_{icj}
==========

\alpha_c(t_j)
+
a_i
+
r_{ic}
+
u_i(t_j).
]

Here:

* (\alpha_c(t)) is the population log-hazard curve for condition (c);
* (a_i) is a subject-wide speed effect;
* (r_{ic}) is a subject-by-condition shift;
* (u_i(t)) is an optional subject-specific shape deviation shared across conditions.

Represent the smooth functions with a basis matrix (B):

[
\alpha_c(t_j)=B_j^\mathsf T\boldsymbol\beta_c,
]

[
u_i(t_j)=\widetilde B_j^\mathsf T\mathbf b_i.
]

The hierarchy becomes

[
a_i\sim N(0,\sigma_a^2),
]

[
\mathbf r_i\sim
N_C(\mathbf0,\Sigma_r),
]

[
b_{ik}\sim N(0,\sigma_{b,k}^2).
]

This model works for OR, AND, and STST because it estimates arbitrary condition distributions, after which the appropriate SFT functional is calculated.

The catch is that it pools generic condition hazards. The prior does not have any particular centre on UCIP, limited capacity, or super capacity. An independent prior on three hazards can induce a surprisingly non-neutral prior over their nonlinear capacity contrast.

For OR capacity, there is a cleaner parameterisation.

---

# 2.4 Recommended OR model: centre the hierarchy on UCIP

Under OR UCIP,

[
h_{i,AB}(t)=h_{i,A}(t)+h_{i,B}(t).
]

This follows because

[
S_{AB}(t)=S_A(t)S_B(t),
]

so

[
H_{AB}(t)=H_A(t)+H_B(t),
]

and therefore their derivatives satisfy the hazard sum identity.

Define the dual-target hazard as

[
h_{i,AB}(t)
===========

\left[
h_{i,A}(t)+h_{i,B}(t)
\right]
\exp{\delta_i(t)}.
]

Equivalently,

[
\log h_{i,AB}(t)
================

\operatorname{logsumexp}
\left[
\log h_{i,A}(t),
\log h_{i,B}(t)
\right]
+
\delta_i(t).
]

This is not a restrictive processing model. For any three positive hazard functions, there exists

[
\delta_i(t)
===========

## \log h_{i,AB}(t)

\log{h_{i,A}(t)+h_{i,B}(t)}.
]

It is simply a reparameterisation in which:

[
\delta_i(t)=0
]

corresponds to hazard-level UCIP.

A parsimonious first hierarchy would be:

[
\log h_{i,A}(t_j)
=================

\alpha_A(t_j)+a_i+\frac{q_i}{2},
]

[
\log h_{i,B}(t_j)
=================

\alpha_B(t_j)+a_i-\frac{q_i}{2},
]

[
\delta_i(t_j)
=============

\delta_0(t_j)+d_i.
]

The random effects mean:

[
a_i\sim N(0,\sigma_a^2)
]

captures overall subject speed,

[
q_i\sim N(0,\sigma_q^2)
]

captures individual channel asymmetry, and

[
d_i\sim N(0,\sigma_d^2)
]

captures each subject’s tendency toward limited or super capacity across the curve.

The population functions are:

[
\alpha_A(t_j)=B_j^\mathsf T\boldsymbol\beta_A,
]

[
\alpha_B(t_j)=B_j^\mathsf T\boldsymbol\beta_B,
]

[
\delta_0(t_j)=B_j^\mathsf T\boldsymbol\beta_\delta.
]

Thus,

[
\log h_{i,AB}(t_j)
==================

\operatorname{logsumexp}
\left[
B_j^\mathsf T\boldsymbol\beta_A+a_i+q_i/2,;
B_j^\mathsf T\boldsymbol\beta_B+a_i-q_i/2
\right]
+
B_j^\mathsf T\boldsymbol\beta_\delta+d_i.
]

This model gives you:

* flexible population hazard shapes for both channels;
* a flexible population capacity trajectory;
* partial pooling of individual speed;
* partial pooling of channel asymmetry;
* partial pooling of individual capacity;
* a prior naturally centred on UCIP through (\boldsymbol\beta_\delta=\mathbf0).

That is much closer to the scientific estimand than independently pooling three hazard curves.

## Important interpretation

Although (\delta_i(t)) is a hazard-level capacity deviation,

[
C_i^{OR}(t)
\neq \exp{\delta_i(t)}
]

in general, because the capacity coefficient is based on integrated hazards:

[
H_{ic}(t_j)
===========

\sum_{\ell\le j}h_{ic\ell}\Delta_\ell.
]

You still derive:

[
C_i^{OR}(t_j)
=============

\frac{H_{i,AB}(t_j)}
{H_{i,A}(t_j)+H_{i,B}(t_j)},
]

and preferably the difference form

[
D_i^{OR}(t_j)
=============

## H_{i,AB}(t_j)

## H_{i,A}(t_j)

H_{i,B}(t_j).
]

The difference form is numerically much better behaved near the left tail.

---

# 2.5 Smoothness priors

A Bayesian P-spline representation would use second-difference penalties:

[
\beta_{c,k}
-----------

2\beta_{c,k-1}
+
\beta_{c,k-2}
\sim
N(0,\sigma_{\beta,c}^2),
]

and similarly

[
\beta_{\delta,k}
----------------

2\beta_{\delta,k-1}
+
\beta_{\delta,k-2}
\sim
N(0,\sigma_{\delta}^2).
]

Then:

[
\sigma_{\beta,c}\sim
\operatorname{HalfNormal}(0,s_\beta),
]

[
\sigma_\delta\sim
\operatorname{HalfNormal}(0,s_\delta).
]

The original paper used AR(1) smoothing and later noted that its autoregressive parameters tended toward the boundary, suggesting that the constant mean structure was inadequate. It explicitly proposed finer-grid Gaussian-process smoothing as a future direction. A low-rank spline gives most of that benefit without the computational weight of a full GP.

I would begin with roughly 8 to 12 spline coefficients over 20 to 30 exposure bins.

---

# 2.6 Stan skeleton for the OR-centred model

This is the central structure, omitting some generated quantities and prior calibration details:

```stan
data {
  int<lower=1> I;
  int<lower=1> J;
  int<lower=3> K;

  matrix[J, K] B;

  array[I, J] int<lower=0> d_A;
  array[I, J] int<lower=0> d_B;
  array[I, J] int<lower=0> d_AB;

  matrix<lower=0>[I, J] exposure_A;
  matrix<lower=0>[I, J] exposure_B;
  matrix<lower=0>[I, J] exposure_AB;
}

parameters {
  vector[K] beta_A;
  vector[K] beta_B;
  vector[K] beta_delta;

  real<lower=0> sigma_smooth_A;
  real<lower=0> sigma_smooth_B;
  real<lower=0> sigma_smooth_delta;

  vector[I] z_speed;
  vector[I] z_asymmetry;
  vector[I] z_capacity;

  real<lower=0> sigma_speed;
  real<lower=0> sigma_asymmetry;
  real<lower=0> sigma_capacity;
}

transformed parameters {
  vector[I] speed = sigma_speed * z_speed;
  vector[I] asymmetry = sigma_asymmetry * z_asymmetry;
  vector[I] capacity_shift = sigma_capacity * z_capacity;

  matrix[I, J] eta_A;
  matrix[I, J] eta_B;
  matrix[I, J] eta_AB;

  vector[J] population_A = B * beta_A;
  vector[J] population_B = B * beta_B;
  vector[J] population_delta = B * beta_delta;

  for (i in 1:I) {
    for (j in 1:J) {
      eta_A[i, j] =
        population_A[j]
        + speed[i]
        + 0.5 * asymmetry[i];

      eta_B[i, j] =
        population_B[j]
        + speed[i]
        - 0.5 * asymmetry[i];

      eta_AB[i, j] =
        log_sum_exp(eta_A[i, j], eta_B[i, j])
        + population_delta[j]
        + capacity_shift[i];
    }
  }
}

model {
  // Non-centred subject effects
  z_speed ~ std_normal();
  z_asymmetry ~ std_normal();
  z_capacity ~ std_normal();

  sigma_speed ~ normal(0, 0.5);
  sigma_asymmetry ~ normal(0, 0.5);
  sigma_capacity ~ normal(0, 0.5);

  // Anchor the first two spline coefficients
  beta_A[1:2] ~ normal(0, 2);
  beta_B[1:2] ~ normal(0, 2);

  // UCIP-centred capacity prior
  beta_delta[1:2] ~ normal(0, 0.5);

  sigma_smooth_A ~ normal(0, 0.5);
  sigma_smooth_B ~ normal(0, 0.5);
  sigma_smooth_delta ~ normal(0, 0.25);

  // Second-order difference penalties
  for (k in 3:K) {
    beta_A[k] - 2 * beta_A[k - 1] + beta_A[k - 2]
      ~ normal(0, sigma_smooth_A);

    beta_B[k] - 2 * beta_B[k - 1] + beta_B[k - 2]
      ~ normal(0, sigma_smooth_B);

    beta_delta[k] - 2 * beta_delta[k - 1] + beta_delta[k - 2]
      ~ normal(0, sigma_smooth_delta);
  }

  // Piecewise-exponential Poisson likelihood
  for (i in 1:I) {
    for (j in 1:J) {
      if (exposure_A[i, j] > 0) {
        d_A[i, j] ~ poisson_log(
          log(exposure_A[i, j]) + eta_A[i, j]
        );
      }

      if (exposure_B[i, j] > 0) {
        d_B[i, j] ~ poisson_log(
          log(exposure_B[i, j]) + eta_B[i, j]
        );
      }

      if (exposure_AB[i, j] > 0) {
        d_AB[i, j] ~ poisson_log(
          log(exposure_AB[i, j]) + eta_AB[i, j]
        );
      }
    }
  }
}
```

The numerical prior scales above are placeholders, not sacred tablets. Converting RTs to seconds makes log-hazard priors vastly easier to reason about.

---

# 2.7 Posterior capacity calculations

For every posterior draw (m), subject (i), and bin (j):

[
h^{(m)}_{icj}
=============

\exp{\eta^{(m)}_{icj}}.
]

With bin width

[
\Delta_j=s_j-s_{j-1},
]

calculate:

[
H^{(m)}_{ic}(t_j)
=================

\sum_{\ell=1}^j
h^{(m)}*{ic\ell}\Delta*\ell.
]

Then:

[
D_i^{OR,(m)}(t_j)
=================

## H^{(m)}_{i,AB}(t_j)

## H^{(m)}_{i,A}(t_j)

H^{(m)}_{i,B}(t_j),
]

[
C_i^{OR,(m)}(t_j)
=================

\frac{
H^{(m)}*{i,AB}(t_j)
}{
H^{(m)}*{i,A}(t_j)+H^{(m)}_{i,B}(t_j)
}.
]

The primary group-level curve should be calculated draw by draw as:

[
\overline D^{(m)}(t)
====================

\frac1I\sum_iD_i^{(m)}(t),
]

or

[
\overline{\log C}^{(m)}(t)
==========================

\frac1I\sum_i\log C_i^{(m)}(t).
]

Do **not** average the hazards first and then calculate capacity. Because capacity is nonlinear,

[
C\left(E[h_i]\right)
\neq
E[C(h_i)].
]

The former describes an artificial population-mixture distribution. The latter describes average within-person capacity.

I would return three levels:

1. population mean capacity trajectory;
2. shrunk subject-specific trajectories;
3. posterior predictive trajectory for a new subject.

For example:

[
P{\overline D(t)>0\mid y},
]

[
P{D_i(t)>0\mid y},
]

and

[
P{D_{\mathrm{new}}(t)>0\mid y}.
]

For a scalar summary, use something explicitly defined from the posterior curve rather than calling it Bayesian (C_z), perhaps:

[
A_i
===

\frac{
\int_{\mathcal T}w(t)D_i(t),dt
}{
\int_{\mathcal T}w(t),dt
},
]

or

[
L_i
===

\frac{
\int_{\mathcal T}w(t)\log C_i(t),dt
}{
\int_{\mathcal T}w(t),dt
}.
]

The weighting function (w(t)) should be fixed in advance, with zero weight in the unstable extreme tails.

---

# 2.8 AND and STST

For a general implementation supporting all stopping rules, fit three arbitrary ordinary-hazard curves:

[
\log h_{ic}(t_j)
================

B_j^\mathsf T\boldsymbol\beta_c
+
a_i+r_{ic},
]

then derive:

[
S_{ic}(t)=\exp{-H_{ic}(t)},
]

[
F_{ic}(t)=1-S_{ic}(t),
]

[
K_{ic}(t)=\log F_{ic}(t).
]

For AND capacity:

[
C_i^{AND}(t)
============

\frac{
K_{i,A}(t)+K_{i,B}(t)
}{
K_{i,AB}(t)
},
]

with signed difference

[
D_i^{AND}(t)
============

## K_{i,AB}(t)

## K_{i,A}(t)

K_{i,B}(t).
]

This sign convention gives:

[
D_i^{AND}(t)>0
]

for super capacity, matching the difference convention used in the current SFT code.

The OR-centred model is cleaner because ordinary hazards add under OR UCIP. There is no equally convenient forward-survival parameterisation for AND because its defining identity is in the CDF and reverse cumulative hazard. I would therefore implement:

* a generic three-hazard model for OR, AND, and STST;
* an optional UCIP-centred OR parameterisation for greater efficiency and interpretability.



