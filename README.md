# sftplus

Extended Systems Factorial Technology analysis in R.

`sftplus` preserves the public analysis and testing interfaces of Joseph Houpt's
`sft` package and adds the CapacityAND implementations: robust hazard/UCIP
estimators, capacity bounds and Altieri capacity, restored UCIP/SIC/MIC tests,
Bayesian distribution-free SIC model probabilities, and LBA/OU simulation
utilities with optional Rcpp backends.

The Bayesian SIC implementation is Monte Carlo integration over a discretized
Dirichlet-process-style prior. Its posterior model probabilities are not the
same quantity as a K-S p-value; priors, binning, Monte Carlo error, and the
model-class tolerances should be reported with results.

The UCIP Bayesian companion uses the score-test information directly. For each
participant it extracts the UCIP numerator and variance, models
`theta_hat = U / V` with its known score precision, and partially pools the
participant effects through a Normal/inverse-Gamma hierarchy. The
single-participant kernel is `ucip.bayes()` (the Bayesian analogue of
`ucip.test()`); the hierarchical, multi-subject estimator is
`capacityGroup.bayes()` (the Bayesian analogue of `capacityGroup()`). Both take
the canonical `inData` trial data frame (a pre-split RT list or nested list of
participant RT lists is still accepted); raw response-time distributions are
never pooled. The result includes
posterior HDIs, directional probabilities, prior/posterior predictive score
draws, and a basic between-chain R-hat diagnostic.

The OR-centred semiparametric hierarchical SFT model is available as
`semiparametricSFT.bayes()`. It accepts
the canonical `inData` trial data frame, uses a common piecewise-exponential grid, and
returns draw-by-draw difference and ratio capacity curves. The default
channel-presence mode pools positive salience into `AB`, `A`, and `B`; pass
`salience_split = TRUE` (or an explicit mapping) to retain `A_L`, `A_H`,
`B_L`, `B_H`, and `AB_LL`, `AB_LH`, `AB_HL`, `AB_HH`. In split mode capacity is
matched within salience and averaged over cells, and the fit additionally
identifies the hierarchical survivor interaction contrast and mean interaction
contrast.

`sic.bayes()` and `mic.bayes()` extract those interaction contrasts. Both are
derived draw-by-draw from the same posterior survivor curves at three
hierarchically distinct levels — the partially pooled subject, the
typical-subject population parameter, and the posterior-predictive new subject —
with sign tests against zero at every level and subject-versus-population
contrasts. This separates the population estimand from a mere average of
per-subject fits, improving on the original per-subject Van Zandt approach. The
MIC, `mean(LL) - mean(LH) - mean(HL) + mean(HH)`, is exactly the signed area
under the posterior SIC. `rstan` remains optional; `semiparametricSFT.bayes(...,
sample = FALSE)` performs validation, exposure construction, and prior-predictive
preparation without fitting.

The list-based UCIP, capacity, resilience, and assessment functions also
accept the canonical row-wise trial data frame. They use `sft_data_to_rt()`
internally to extract `AB`, `A`, and `B` (or STST context/target) cells; for
multiple subjects, `capacityGroup.bayes()` performs the corresponding hierarchical
conversion automatically.

All row-wise data-frame entry points use the same canonical names and aliases:
`RT` may be supplied as `rt`, `Subject` as `subjects`, and `Condition` as
`LogicalRule`. Canonical names take precedence if both a canonical column and
one of its aliases are present.

The LBA simulator also supports an optional shared-capacity fluctuation on
double-target (`AB`) trials. Set `kappa` and `tau` in `p_vec` (or pass them
directly to `simulate_sft()`):

```r
p_vec <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
           ve = .8, sv_c = .1, sv_e = .1, kappa = 1, tau = .25)
sim <- simulate_sft("lba", n = 1000, p_vec = p_vec,
                    design = c("AB", "AN", "NB", "NN"),
                    logical_rules = "AND")
```

For each AB trial, one `Z ~ N(0, 1)` is drawn and the two target drifts use
`V_i = (kappa + tau * Z) * v_i + epsilon_i`. Thus `kappa` controls the mean
target multiplier and `tau` controls shared trial-to-trial co-fluctuation;
the nontargets and non-AB cells are unchanged. The existing additive
`capacity_target` modifier can be used at the same time. Defaults are
`kappa = 1` and `tau = 0` (`rho` is accepted as an alias for `tau`).

`simulate_sft()` defaults to the parallel architecture. Serial simulations can
be requested with `architecture = "serial"`. Select
`serial_mode = "self-terminating"` (the default for serial calls) or
`serial_mode = "exhaustive"`, and choose `serial_order = "A"`, `"B"`, or
`"random"`. Random order checks channel A first on a proportion `p_A_first`
(default `0.5`) of trials. Self-terminating OR trials stop at the first target
or after both absences; explicitly requested self-terminating AND trials stop
at the first absence or after both targets. Serial exhaustive trials complete
both channels in sequence. Serial self-termination is currently available for
OR and AND (the default serial self-terminating logical rule is OR); use
exhaustive processing for ID data. Since serial channels do not run
concurrently, LBA capacity/shared-capacity modifiers (including `kappa` and
the `tau`/`rho` fluctuation parameter) and OU cross-channel interaction
parameters are set to zero in serial runs.
