# Split monolithic R files into modules by top-level symbol, with strict
# line-coverage checking so no code can be dropped or duplicated.
# Usage: Rscript split_modules.R [--write]
WRITE <- "--write" %in% commandArgs(TRUE)

# ---- boundary detection -------------------------------------------------
# A top-level definition starts at column 0 with `symbol <- ` (all defs in
# these files use `<-`). Embedded Stan/C++ strings use `{`/`=`, never col-0 `<-`.
find_defs <- function(lines) {
  rx <- "^(`[^`]+`|[.A-Za-z][.A-Za-z0-9_]*)[[:space:]]*<-"
  rx_full <- "^(`[^`]+`|[.A-Za-z][.A-Za-z0-9_]*)[[:space:]]*<-.*$"
  idx <- grep(rx, lines, perl = TRUE)
  sym <- sub(rx_full, "\\1", lines[idx], perl = TRUE)
  sym <- gsub("`", "", sym)
  data.frame(line = idx, sym = sym, stringsAsFactors = FALSE)
}

# Extend each def's start upward over its contiguous preceding comment/blank
# block, so a function's leading comments travel with it.
comment_start <- function(lines, start, floor) {
  i <- start - 1L
  while (i >= floor && (grepl("^[[:space:]]*#", lines[i]) || !nzchar(trimws(lines[i])))) {
    i <- i - 1L
  }
  i + 1L
}

partition <- function(path) {
  lines <- readLines(path, warn = FALSE)
  defs <- find_defs(lines)
  n <- nrow(defs)
  # block start for each def = top of its attached comment block, but never
  # before the previous def's own line.
  cstart <- integer(n)
  for (k in seq_len(n)) {
    floor <- if (k == 1L) 1L else defs$line[k - 1L] + 1L
    cstart[k] <- comment_start(lines, defs$line[k], floor)
  }
  # block k spans [cstart[k], cstart[k+1]-1]; last to EOF.
  bstart <- cstart
  bend <- c(cstart[-1L] - 1L, length(lines))
  header <- if (cstart[1L] > 1L) seq_len(cstart[1L] - 1L) else integer(0)
  list(lines = lines, sym = defs$sym, bstart = bstart, bend = bend,
       header = header, total = length(lines))
}

get_block <- function(p, k) p$lines[p$bstart[k]:p$bend[k]]

# ---- mappings -----------------------------------------------------------
core_map <- list(
  sft_utils.R = c("%||%", ".sft_htest", ".sft_hdi", ".sft_bind_rows", ".to_correct_indicator"),
  data_prep.R = c(".sft_zero_curve", ".sft_curve", ".sft_finite_times", ".sft_correct",
                  ".sft_clean", ".sft_cr_list", ".sft_data_numeric", ".sft_data_correct",
                  "sft_data_to_rt", ".sft_as_rt_cr", ".sft_ecdf", "row_to_pvec", ".safe_ecdf"),
  hazard.R = c(".sft_hazard_env", ".ensure_hazard_rcpp", "estimateNAH", "estimateNAK",
               "estimateUCIPor", "estimateUCIPand"),
  capacity.R = c(".sft_capacity_bounds_or", ".sft_capacity_bounds_and", "capacity.or",
                 "capacity.and", "capacity.stst", "capacity.id", "capacity.altieri",
                 "capacityGroup"),
  ucip.R = c(".sft_ucip_components", ".sft_score_method", ".sft_ucip_score_components",
             ".sft_ucip_score", "ucip.test", "ucip.id.test"),
  ucip_bayes_engine.R = c(".sft_subject_ucip_input", ".sft_hierarchical_cz_chain",
                          ".sft_bayes_method", ".sft_validate_bayes_args", ".sft_ucip_score_data",
                          ".sft_bayes_seed", ".sft_restore_bayes_seed", ".sft_parameter_summary",
                          ".sft_split_chains", ".sft_fallback_ess", ".sft_manual_diagnostics",
                          ".sft_mcmc_diagnostics", ".sft_prior_predictive",
                          ".sft_posterior_predictive", ".sft_shrinkage_draws",
                          ".sft_shrinkage_summary", ".sft_stan_model_cache", ".sft_stan_code",
                          ".sft_stan_sampler_diagnostics", ".sft_run_stan_hierarchy"),
  ucip_bayes.R = c("ucip.bayes", "capacityGroup.bayes"),
  sic.R = c("sic.test", "mic.test", ".sft_bayes_pair", "siDominance", ".sft_dp_checkmods",
            "sicDPtest", "sictestBayes", "sic", "sicGroup", "sicGroupBF"),
  plot_builders.R = c("build_ct_df", "smooth_one_cdf", "build_sic_df", "build_at_df",
                      "build_cdf_df", "get_time_range", "plot_altieri"),
  assessment_ta.R = c(".at_clean", ".at_count", ".at_defective", ".at_log_ratio", ".at_pack",
                      "assessment_ta_or", "assessment_ta_and", ".at_weighted_survivor_integral",
                      "assessment_donkin_discrimination"),
  # appended to existing files
  `+resilience.R` = c("resilience.test", "resilience"),
  `+simulation.R` = c(".normalize_logical_rules", ".default_design", ".char_to_level",
                      ".expand_design", ".annotate_trials", ".split_by_rule",
                      "obj_max_vs_min", "find_equal_vc_yes")
)

bayes_map <- list(
  semiparametric_prep.R = c(".sft_bayes_series", ".sft_bayes_salience_series",
                            ".sft_bayes_salience_cells", ".sft_bayes_numeric", ".sft_bayes_correct",
                            ".sft_bayes_salience_labels", ".sft_bayes_resolve_salience",
                            ".sft_bayes_validate_input", ".sft_bayes_grid", ".sft_bayes_basis",
                            ".sft_bayes_units", ".sft_bayes_priors", ".sft_bayes_prepared",
                            ".sft_bayes_stan_path", ".sft_bayes_stan_cache", ".sft_bayes_stan_model",
                            ".sft_capacity_bayes_seed", ".sft_bayes_validate_sampling"),
  semiparametric_posterior.R = c(".sft_bayes_dim_array", ".sft_bayes_transform",
                                 ".sft_bayes_cumulative_series", ".sft_bayes_transform_salience",
                                 ".sft_bayes_summary_matrix", ".sft_bayes_split_summaries",
                                 ".sft_bayes_summaries", ".sft_bayes_survivor_dxJ",
                                 ".sft_bayes_population_sic", ".sft_bayes_integrate_sic",
                                 ".sft_bayes_scalar_sign", ".sft_bayes_mic_result",
                                 ".sft_bayes_sic_result", ".sft_bayes_effects_from_posterior",
                                 ".sft_bayes_effects_from_salience_posterior",
                                 ".sft_bayes_prior_draws", ".sft_bayes_prior_predictive",
                                 ".sft_bayes_prior_draws_salience",
                                 ".sft_bayes_predictive_from_series",
                                 ".sft_bayes_prior_predictive_salience", ".sft_bayes_diagnostics",
                                 ".sft_bayes_posterior_predictive", ".sft_bayes_empty_transformed",
                                 ".sft_bayes_empty_curve_data"),
  semiparametric_sft.R = c("semiparametricSFT.bayes", "print.sft_bayes", "mic.bayes", "sic.bayes")
)

headers <- c(
  sft_utils.R = "# Small cross-cutting helpers shared across the sftplus internals.",
  data_prep.R = "# Canonical SFT data adapters and response-time curve/ECDF utilities.",
  hazard.R = "# Nelson-Aalen hazard/cumulative-hazard estimators and the Rcpp kernel.",
  capacity.R = "# Workload capacity coefficients (OR/AND/STST/ID/Altieri) and the group orchestrator.",
  ucip.R = "# Houpt-Townsend UCIP frequentist tests and the score/information machinery.",
  ucip_bayes_engine.R = "# Hierarchical Bayesian UCIP estimation engine: samplers, diagnostics,\n# predictive checks, and shrinkage helpers. No public entry points.",
  ucip_bayes.R = "# Public Bayesian UCIP orchestrators: ucip.bayes() (single subject) and\n# capacityGroup.bayes() (hierarchical group). They score kernels and wrap the\n# estimation engine in ucip_bayes_engine.R.",
  sic.R = "# Survivor and mean interaction contrasts: frequentist and Bayes-factor SIC/MIC.",
  plot_builders.R = "# Tidy data-frame builders and plotting helpers for curve visualisation.",
  assessment_ta.R = "# Timed-assessment (trial-abundance) and Donkin discrimination extensions.",
  semiparametric_prep.R = "# Semiparametric hierarchical SFT model: input validation, salience resolution,\n# pooled RT grid, spline basis, priors, and Stan-data preparation.",
  semiparametric_posterior.R = "# Semiparametric hierarchical SFT model: posterior transforms, curve summaries,\n# and SIC/MIC/effects/predictive derivations.",
  semiparametric_sft.R = "# Semiparametric hierarchical SFT model: public API (semiparametricSFT.bayes(),\n# print method, and the sic.bayes()/mic.bayes() extractors)."
)

process <- function(path, map, dir = "R") {
  p <- partition(path)
  found <- p$sym
  wanted <- unlist(map, use.names = FALSE)
  cat("\n==== ", path, " ====\n", sep = "")
  cat("defs found:", length(found), " mapped:", length(wanted), "\n")
  miss <- setdiff(found, wanted)
  extra <- setdiff(wanted, found)
  if (length(miss)) cat("!! FOUND-BUT-UNMAPPED:", paste(miss, collapse = ", "), "\n")
  if (length(extra)) cat("!! MAPPED-BUT-NOT-FOUND:", paste(extra, collapse = ", "), "\n")
  if (any(duplicated(found))) cat("!! DUPLICATE SYMBOLS:", paste(found[duplicated(found)], collapse=", "), "\n")
  # coverage: every line from bstart[1] to EOF assigned exactly once; header separate
  covered <- integer(0)
  outputs <- list()
  for (fname in names(map)) {
    syms <- map[[fname]]
    ks <- match(syms, found)
    if (anyNA(ks)) next
    ks <- ks[order(ks)]  # preserve original order
    blk <- unlist(lapply(ks, function(k) c(get_block(p, k), "")), use.names = FALSE)
    for (k in ks) covered <- c(covered, p$bstart[k]:p$bend[k])
    outputs[[fname]] <- blk
  }
  body_lines <- if (length(p$bstart)) p$bstart[1]:p$total else integer(0)
  uncov <- setdiff(body_lines, covered)
  dup <- covered[duplicated(covered)]
  cat("body lines:", length(body_lines), " covered:", length(unique(covered)),
      " uncovered:", length(uncov), " double-covered:", length(dup), "\n")
  if (length(uncov)) cat("!! UNCOVERED LINES:", paste(head(uncov,20), collapse=","), "\n")
  if (length(dup)) cat("!! DOUBLE-COVERED:", paste(head(dup,20), collapse=","), "\n")
  ok <- !length(uncov) && !length(dup) && !length(miss) && !length(extra) && !anyNA(match(found, found))
  if (WRITE && ok) {
    for (fname in names(map)) {
      blk <- outputs[[fname]]
      if (is.null(blk)) next
      append <- startsWith(fname, "+")
      real <- sub("^\\+", "", fname)
      target <- file.path(dir, real)
      if (append) {
        con <- file(target, "a")
        writeLines(c("", headersFor(real), blk), con); close(con)
        cat("appended ->", target, "\n")
      } else {
        writeLines(c(headers[[fname]], "", blk), target)
        cat("wrote    ->", target, "\n")
      }
    }
  }
  ok
}
headersFor <- function(real) {
  # section banner when appending to an existing file
  paste0("# ---- moved from the former monolith ----")
}

ok1 <- process("R/sft_core.R", core_map)
ok2 <- process("R/sft_bayes.R", bayes_map)
cat("\nALL OK:", ok1 && ok2, " (WRITE=", WRITE, ")\n", sep = "")
