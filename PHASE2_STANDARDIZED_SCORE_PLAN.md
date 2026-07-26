> **Status: superseded (2026-07-26).** Phase 2 shipped `"standardized"` as the
> *default* and dropped `"score"`. That was wrong in one specific way: a fixed
> `prior_sd` on `phi` implies a prior on `theta` of width
> `prior_sd / sqrt(V_ref)`, which tightens at exactly the rate the likelihood
> sharpens and therefore is never outweighed by data. On an eight-subject set
> with every `Cz` near 8.5 it drove the population mean to 0.39 and `BF10` to
> 1.36. The plan below verified invariance to *weight rescaling* but never to
> *sample size*, which is the property that actually matters.
>
> Resolution: `"score"` (theta, the Peto one-step log hazard ratio) is restored
> as the default effect scale; `"standardized"` is retained and its default
> priors are the theta priors pushed through the `sqrt(V_ref)` map, so the two
> are the same model reported in different units. See `.sft_default_priors()`.

# Phase 2 — `score_method = "standardized"` (reference-information Cz scale)

Follow-up to Phase 1 (Bayes factors + `"capacity"`→`"multiplicative"` rename). Phase 1
added `bayes_factor` outputs and matched each null test to the effect scale, but the raw
`score` effect `theta_i = U_i / V_i` is on an information-referenced scale whose alternative
prior is nearly uninterpretable: with a `Normal(0, prior_sd = 1)` prior, a typical
participant (`V_i ≈ 35`, `sqrt(V_i) ≈ 5.9`) has a prior SD of ~5.9 on the expected `Cz`
scale — a `±12` 95% range. That diffuseness inflates the point-null Bayes factor
(Jeffreys–Lindley) and leaves the ROPE undefinable, which is why Phase 1 gives the raw
`score` scale a point null only.

Phase 2 adds a third `score_method`, `"standardized"`, that rescales the score effect onto
a **reference-information** unit so a `prior_sd = 1` alternative means "one `Cz`-unit shift
for a participant with typical information." This makes both the interval null and the
alternative prior coherent, without changing the frequentist `Cz` statistic at all.

## The transformation

Let `V_ref = median(V_i)` over the participants entering the fit. Define

```
phi_hat_i        = theta_hat_i * sqrt(V_ref)          # = (U_i / V_i) * sqrt(V_ref)
precision(phi_i) = V_i / V_ref
```

Properties (all verified algebraically in the Phase-1 review):

- **`Cz` unchanged:** `Cz_i = phi_hat_i * sqrt(V_i / V_ref) = theta_hat_i * sqrt(V_i)`.
- **Point null invariant:** `mu_phi = 0` iff `mu_theta = 0`, so the UCIP point-null test is
  unaffected in *location*; only the alternative-prior *width* changes.
- **Score-weight invariant:** multiplying all score weights by `c` sends
  `theta_hat → theta_hat / c`, `V_ref → c^2 V_ref`, so `phi_hat` is unchanged. This removes
  the arbitrary-weight-scaling problem the raw `theta` scale has.
- **Interval null now meaningful:** one `phi`-unit is a `Cz` shift of 1 for a typical-
  information participant, so a `rope` on `phi` is interpretable and the interval Bayes
  factor applies.

`phi` is still an information-referenced *statistical* scale (not a physical RT scale), so it
is not a scientific ROPE the way the MIC's proportional region is — but it is a stable,
weight-invariant scale for alternative priors and Bayes factors. Document this caveat.

## Scope decision (ask the user before coding)

**Keep `"standardized"` opt-in; do not change the `"score"` default.** Changing the default
would silently move every existing `capacityGroup.bayes`/`ucip.bayes` score-effect posterior,
ROPE meaning, and Bayes factor — bad for CRAN stability and reproducibility. The three
choices become `c("score", "multiplicative", "standardized")` with `"score"` still first
(default). Revisit making it the default only if the user explicitly wants it.

`"standardized"` applies **only to the UCIP score model**. `"multiplicative"` (`eta`) and the
MIC `rho` are already dimensionless and already coherent — leave them alone.

## Implementation steps

### 1. `.sft_score_method` (`R/ucip.R`)
Add `standardized = "standardized"` (plus aliases: `std`, `referenceinformation`,
`phi`). Update the error string to list all three canonical methods.

### 2. `.sft_ucip_score_data` (`R/ucip_bayes_engine.R`)
This is where the scaling lives. Currently:

```r
theta_hat  <- U / V           # per subject
precision  <- V
effect_name<- "theta"         # for "score"
```

For `method == "standardized"`:
- Compute `V_ref <- stats::median(V)` **after** degenerate-subject dropping (so `V_ref`
  matches the participants actually fitted — important for the single-subject and
  post-drop hierarchy cases).
- `effect_hat <- theta_hat * sqrt(V_ref)`
- `precision  <- V / V_ref`
- `effect_name <- "phi"`
- Bootstrap path: the within-subject bootstrap already resamples `theta_hat = U/V`; scale
  each replicate’s implied effect by `sqrt(V_ref)`, or equivalently multiply the bootstrap
  precision by `V_ref` (since `var(phi_hat) = V_ref * var(theta_hat)` ⇒
  `precision_phi = precision_theta / V_ref`). Confirm which of `.sft_ucip_boot_effect` /
  `.sft_ucip_boot_precision` needs the factor and apply it once, consistently.
- Return `V_ref` in the result list and expose it in the `score` data frame (e.g. a
  `V_ref` column and a `phi_hat` column alongside the existing `theta_hat`, `Cz`, `se`).

Watch the single-subject route: `ucip.bayes` pulls `scored$effect_hat[[1L]]` /
`scored$precision[[1L]]`, so `V_ref = median(V)` of one subject is just that subject's `V`,
giving `phi_hat = Cz` and `precision = 1`. That is correct and desirable (the single-subject
standardized effect *is* `Cz`), but note it in the code and in a test.

### 3. Precision-name / Stan data (`R/ucip_bayes_engine.R`, `.sft_run_stan_hierarchy`)
`precision_name` is currently `if (effect_name == "eta") "P" else "V"`. Decide the Stan data
name for `phi` (e.g. `"W"`), and mirror it in `precision_name <- if (score_method == ...)`
in `ucip_bayes.R`. The Stan model code is generated from `effect_name`/`precision_name`, so
`"phi"` + a distinct precision name flows through `.sft_stan_code` and the model cache key
without further change — just verify the cache key (`paste0(method, "_", prior_tau_sd, "_",
effect_name)`) stays unique.

### 4. Orchestrators (`R/ucip_bayes.R`)
- `ucip.bayes`: `precision_name` branch, `interpretation` string for `phi`, and the BF
  gating. For `"standardized"`, **both** tests apply:
  `interval = score_method %in% c("multiplicative", "standardized")`, `point = TRUE`.
- `capacityGroup.bayes`: same gating change; nothing else structural since it reuses
  `.sft_normal_hierarchy_fit` and the Phase-1 assemblers.
- Add a `phi` interpretation branch to the `model$interpretation` blocks:
  "phi is the reference-information standardized score effect (Cz for a participant with
  median information)".

### 5. Bayes factors
No new helper needed — Phase 1's `.sft_analytic_bayes_factor` / `.sft_hierarchy_bayes_factor`
already take `point`/`interval` flags. Just flip `interval = TRUE` for `"standardized"` and
let `rope` supply δ (now on the `phi` scale).

### 6. Docs (`man/sft_plus-extensions.Rd`)
- `\usage`: `score_method = c("score", "multiplicative", "standardized")` in both signatures.
- Extend the `score_method` `\details` paragraph: define `phi`, `V_ref = median(V_i)`, the
  `Cz`-invariance and weight-invariance, and the coherent-prior/interval-null benefit; state
  the caveat that it is still an information-referenced (not physical-RT) scale.
- Note that `"standardized"` reports both `point_null` and `interval_null`.

### 7. Tests (`tests/testthat/test-extensions.R`)
- **`Cz` invariance:** `score_method = "standardized"` yields the same `score$Cz` as
  `"score"` (and the frequentist `ucip.test`), while `effect_hat`/posterior differ.
- **`V_ref` and definition:** `phi_hat == theta_hat * sqrt(median(V))` and
  `precision == V / median(V)` for a multi-subject fit.
- **Single subject:** `phi_hat == Cz` and `precision == 1` (median of one `V`).
- **Weight invariance:** scaling the score weights (or an RT unit change that scales `U`,`V`)
  leaves `phi_hat` unchanged.
- **Bayes factors present:** `bayes_factor` has both `point_null` and `interval_null` when a
  `rope` is supplied; interval δ equals `rope` on the `phi` scale.
- Confirm the point-null `BF01` for `"standardized"` differs from `"score"` at the *same*
  `prior_sd` (demonstrating the tighter, coherent alternative), and is closer to the
  Jeffreys-scale expectation.

### 8. Regression / verification
```
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-extensions.R")'
```
Then a smoke run of `ucip.bayes` and `capacityGroup.bayes` with `score_method =
"standardized"`, `rope = <small>`, confirming `Cz` matches the `"score"` fit and both BFs
appear. Update `INTEGRATION_NOTES.md` if it enumerates score methods.

## Out of scope (possible Phase 3)
- Spike-and-slab / per-subject model comparison for the individual UCIP null
  (`theta_i = 0` for each `i`), and the "everyone is UCIP" null (`mu = 0` **and** `tau = 0`).
  Phase 1's population point null answers only `mu = 0` with `tau` free.
- A physical-RT ROPE for the score model (would require choosing an RT reference like the
  MIC's grand-mean normalization; the MIC already has this, the score model does not).
