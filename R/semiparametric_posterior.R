# Semiparametric hierarchical SFT model: posterior transforms, curve summaries,
# and SIC/MIC/effects/predictive derivations.


# Lower-triangular indicator M[k, j] = 1 iff k <= j, cached per grid size. Right-
# multiplying by it turns a row of per-bin quantities into its running total.
.sft_bayes_cumsum_env <- new.env(parent = emptyenv())

.sft_bayes_cumsum_matrix <- function(J) {
  key <- as.character(J)
  if (!exists(key, envir = .sft_bayes_cumsum_env, inherits = FALSE)) {
    assign(key, upper.tri(matrix(0, J, J), diag = TRUE) * 1,
           envir = .sft_bayes_cumsum_env)
  }
  get(key, envir = .sft_bayes_cumsum_env, inherits = FALSE)
}


# Cumulative hazard from log-hazards on a common bin grid:
#   H[..., j] = sum_{k <= j} exp(eta[..., k]) * widths[k]
# eta may be draws-by-bin or draws-by-subject-by-bin; the bin index is always
# last. Flattening the leading margins and using one matrix product replaces the
# per-draw (and per-subject) R loops this file used to run, which dominated the
# cost of every posterior transform and predictive check.
.sft_bayes_cumhaz <- function(eta, widths) {
  d <- dim(eta)
  if (is.null(d)) d <- c(1L, length(eta))
  J <- d[[length(d)]]
  if (length(widths) != J) {
    stop("widths must have one entry per bin.", call. = FALSE)
  }
  x <- exp(as.numeric(eta))
  dim(x) <- c(length(x) / J, J)
  out <- (x * rep(widths, each = nrow(x))) %*% .sft_bayes_cumsum_matrix(J)
  dim(out) <- d
  dimnames(out) <- dimnames(eta)
  out
}


# Mean over one margin of an array, dropping it. Equivalent to
#   apply(x, setdiff(seq_along(dim(x)), margin), mean, na.rm = na.rm)
# but vectorised. apply() issues one scalar mean() call per retained cell, which
# on a draws-by-subject-by-bin-by-cell array is hundreds of thousands of calls
# and dominated the cost of every salience transform and predictive check.
.sft_bayes_mean_over <- function(x, margin, na.rm = TRUE) {
  d <- dim(x)
  if (is.null(d)) return(mean(x, na.rm = na.rm))
  keep <- setdiff(seq_along(d), margin)
  y <- if (identical(as.integer(margin), length(d))) x else aperm(x, c(keep, margin))
  flat <- matrix(as.numeric(y), ncol = prod(d[margin]))
  if (isTRUE(na.rm)) {
    ok <- !is.na(flat)
    cnt <- rowSums(ok)
    flat[!ok] <- 0
    out <- rowSums(flat) / cnt
    # apply(mean, na.rm = TRUE) yields NaN for an all-missing cell; match it.
    out[cnt == 0L] <- NaN
  } else {
    out <- rowMeans(flat)
  }
  if (length(keep) > 1L) dim(out) <- d[keep]
  out
}


.sft_bayes_dim_array <- function(x, draws, I, J, name) {
  d <- dim(x)
  if (is.null(d) || length(d) != 3L || !all(d == c(draws, I, J))) {
    stop("Unexpected posterior dimension for ", name, ".", call. = FALSE)
  }
  x
}


.sft_bayes_transform <- function(log_hazard, widths, times, subjects,
                                 new_subject = NULL) {
  required <- c("A", "B", "AB")
  if (!all(required %in% names(log_hazard))) stop("log_hazard needs A, B, and AB arrays.",
                                                     call. = FALSE)
  draws <- dim(log_hazard$A)[1L]; I <- dim(log_hazard$A)[2L]; J <- dim(log_hazard$A)[3L]
  H <- lapply(log_hazard, .sft_bayes_cumhaz, widths = widths)
  D <- H$AB - H$A - H$B
  denom <- H$A + H$B
  C <- H$AB / denom
  C[!is.finite(C)] <- NA_real_
  dimnames(D) <- list(draw = seq_len(draws), Subject = subjects, bin = seq_along(times))
  dimnames(C) <- dimnames(D)
  for (nm in names(H)) dimnames(H[[nm]]) <- dimnames(D)
  out <- list(H_A = H$A, H_B = H$B, H_AB = H$AB, D = D, C = C,
              subject = list(H_A = H$A, H_B = H$B, H_AB = H$AB, D = D, C = C),
              times = times, widths = widths, subjects = subjects,
              # identically D; kept as an alias rather than recomputed.
              identity_error = D)
  if (!is.null(new_subject)) {
    new_H <- lapply(new_subject, .sft_bayes_cumhaz, widths = widths)
    new_D <- new_H$AB - new_H$A - new_H$B
    new_C <- new_H$AB / (new_H$A + new_H$B); new_C[!is.finite(new_C)] <- NA_real_
    out$new_subject <- list(H_A = new_H$A, H_B = new_H$B, H_AB = new_H$AB,
                            D = new_D, C = new_C)
  }
  out
}


.sft_bayes_cumulative_series <- function(log_hazard, widths, subjects, times) {
  draws <- dim(log_hazard[[1L]])[1L]
  I <- dim(log_hazard[[1L]])[2L]
  J <- dim(log_hazard[[1L]])[3L]
  H <- lapply(log_hazard, function(x) {
    if (is.null(dim(x)) || !all(dim(x) == c(draws, I, J)))
      stop("All salience posterior arrays must have draw by subject by bin dimensions.", call. = FALSE)
    .sft_bayes_cumhaz(x, widths)
  })
  dn <- list(draw = seq_len(draws), Subject = subjects, bin = seq_along(times))
  for (nm in names(H)) dimnames(H[[nm]]) <- dn
  H
}


.sft_bayes_transform_salience <- function(log_hazard, widths, times, subjects,
                                          new_subject = NULL) {
  required <- .sft_bayes_salience_series
  if (!all(required %in% names(log_hazard))) {
    stop("Salience-split posterior draws must contain A_L, A_H, B_L, B_H, and AB_LL/AB_LH/AB_HL/AB_HH.",
         call. = FALSE)
  }
  H <- .sft_bayes_cumulative_series(log_hazard, widths, subjects, times)
  draws <- dim(H[[1L]])[1L]; I <- dim(H[[1L]])[2L]; J <- dim(H[[1L]])[3L]
  cells <- .sft_bayes_salience_cells
  Dcells <- array(NA_real_, dim = c(draws, I, J, 4L),
                  dimnames = list(draw = seq_len(draws), Subject = subjects,
                                  bin = seq_along(times), cell = cells))
  Ccells <- Dcells
  for (cc in seq_along(cells)) {
    cell <- cells[[cc]]
    u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
    ha <- H[[paste0("A_", u)]]; hb <- H[[paste0("B_", v)]]
    hab <- H[[paste0("AB_", cell)]]
    Dcells[, , , cc] <- hab - ha - hb
    denom <- ha + hb
    z <- hab / denom; z[!is.finite(z)] <- NA_real_
    Ccells[, , , cc] <- z
  }
  Davg <- .sft_bayes_mean_over(Dcells, 4L)
  logCavg <- .sft_bayes_mean_over(log(Ccells), 4L)
  Cavg <- exp(logCavg)
  sic <- exp(-H$AB_LL) - exp(-H$AB_LH) - exp(-H$AB_HL) + exp(-H$AB_HH)
  dn <- list(draw = seq_len(draws), Subject = subjects, bin = seq_along(times))
  dimnames(Davg) <- dn; dimnames(Cavg) <- dn; dimnames(logCavg) <- dn
  dimnames(sic) <- dn
  cap <- list(D = Dcells, C = Ccells, D_average = Davg,
              logC_average = logCavg, C_average = Cavg,
              cells = cells, weights = rep(1 / 4, 4L))
  out <- list(H = H, S = lapply(H, function(x) exp(-x)),
              H_A_L = H$A_L, H_A_H = H$A_H, H_B_L = H$B_L, H_B_H = H$B_H,
              H_AB_LL = H$AB_LL, H_AB_LH = H$AB_LH,
              H_AB_HL = H$AB_HL, H_AB_HH = H$AB_HH,
              D_cells = Dcells, C_cells = Ccells, D = Davg, C = Cavg,
              logC_average = logCavg, SIC = sic, capacity = cap,
              subject = list(D = Davg, C = Cavg, logC = logCavg,
                             D_cells = Dcells, C_cells = Ccells, SIC = sic),
              times = times, widths = widths, subjects = subjects,
              identity_error = Dcells)
  if (!is.null(new_subject)) {
    new_H <- lapply(new_subject, function(x) {
      if (is.null(dim(x)) || length(dim(x)) != 2L || dim(x)[1L] != draws ||
          dim(x)[2L] != J) stop("Unexpected new-subject salience posterior dimensions.", call. = FALSE)
      .sft_bayes_cumhaz(x, widths)
    })
    new_D <- array(NA_real_, dim = c(draws, J, 4L),
                   dimnames = list(draw = seq_len(draws), bin = seq_along(times), cell = cells))
    new_C <- new_D
    for (cc in seq_along(cells)) {
      cell <- cells[[cc]]; u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
      ha <- new_H[[paste0("A_", u)]]; hb <- new_H[[paste0("B_", v)]]
      hab <- new_H[[paste0("AB_", cell)]]
      new_D[, , cc] <- hab - ha - hb
      z <- hab / (ha + hb); z[!is.finite(z)] <- NA_real_; new_C[, , cc] <- z
    }
    new_Davg <- .sft_bayes_mean_over(new_D, 3L)
    new_logCavg <- .sft_bayes_mean_over(log(new_C), 3L)
    new_Cavg <- exp(new_logCavg)
    new_sic <- exp(-new_H$AB_LL) - exp(-new_H$AB_LH) -
      exp(-new_H$AB_HL) + exp(-new_H$AB_HH)
    out$new_subject <- list(H = new_H, D = new_Davg, C = new_Cavg,
                            logC = new_logCavg, D_cells = new_D, C_cells = new_C,
                            SIC = new_sic)
  }
  out
}


.sft_bayes_summary_matrix <- function(x, times, level, measure, subject = NA_character_,
                                      report, hdi = 0.94, rope = 0, series_label = "OR") {
  # x is draws x time for one curve.
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L, ncol = length(times))
  if (ncol(x) != length(times)) x <- matrix(as.numeric(x), nrow = dim(x)[1L], ncol = length(times))
  ans <- lapply(seq_along(times), function(j) {
    v <- x[, j]
    interval <- .sft_hdi(v, hdi)
    threshold <- if (measure %in% c("difference", "D", "SIC")) 0 else 1
    data.frame(Level = level, Subject = subject, Series = series_label, Measure = measure,
               Time = times[[j]], Ct = mean(v, na.rm = TRUE),
               Mean = mean(v, na.rm = TRUE), Lower = interval[[1L]],
               Upper = interval[[2L]],
               Prob_super = mean(v > threshold, na.rm = TRUE),
               Prob_limited = mean(v < threshold, na.rm = TRUE),
               Prob_rope = if (rope > 0) mean(abs(v - threshold) <= rope, na.rm = TRUE) else NA_real_,
               Report = report[[j]], stringsAsFactors = FALSE)
  })
  do.call(rbind, ans)
}


.sft_bayes_split_summaries <- function(transformed, grid, subjects, hdi = 0.94,
                                       rope = 0) {
  times <- grid$bins$upper; report <- grid$bins$report
  parts <- list()
  add_matrix <- function(x, level, measure, subject = NA_character_, series = "OR-average") {
    parts[[length(parts) + 1L]] <<- .sft_bayes_summary_matrix(
      x, times, level, measure, subject, report, hdi, rope, series)
  }
  for (measure in c("D", "C")) {
    for (i in seq_along(subjects))
      add_matrix(transformed$subject[[measure]][, i, , drop = FALSE][, 1L, ],
                 "subject", measure, subjects[[i]])
    add_matrix(transformed$population[[measure]], "population", measure)
    if (!is.null(transformed$new_subject)) add_matrix(transformed$new_subject[[measure]],
                                                     "new_subject", measure = measure,
                                                     series = "OR-average")
  }
  for (cc in seq_along(.sft_bayes_salience_cells)) {
    cell <- .sft_bayes_salience_cells[[cc]]
    for (measure in c("D", "C")) {
      for (i in seq_along(subjects))
        add_matrix(transformed$subject[[paste0(measure, "_cells")]][, i, , cc, drop = FALSE][, 1L, , 1L],
                   "subject", measure, subjects[[i]], cell)
      add_matrix(transformed$population[[paste0(measure, "_cells")]][, , cc],
                 "population", measure, series = cell)
      if (!is.null(transformed$new_subject))
        add_matrix(transformed$new_subject[[paste0(measure, "_cells")]][, , cc],
                   "new_subject", measure = measure, series = cell)
    }
  }
  for (i in seq_along(subjects))
    add_matrix(transformed$subject$SIC[, i, , drop = FALSE][, 1L, ],
               "subject", "SIC", subjects[[i]], "SIC")
  add_matrix(transformed$population$SIC, "population", "SIC", series = "SIC")
  if (!is.null(transformed$new_subject)) add_matrix(transformed$new_subject$SIC,
                                                   "new_subject", measure = "SIC",
                                                   series = "SIC")
  tidy <- do.call(rbind, parts); rownames(tidy) <- NULL
  list(tidy = tidy,
       subject = tidy[tidy$Level == "subject", , drop = FALSE],
       population = tidy[tidy$Level == "population", , drop = FALSE],
       new_subject = tidy[tidy$Level == "new_subject", , drop = FALSE])
}


.sft_bayes_summaries <- function(transformed, grid, subjects, hdi = 0.94, rope = 0) {
  times <- grid$bins$upper; report <- grid$bins$report
  parts <- list()
  for (measure in c("D", "C")) {
    for (i in seq_along(subjects)) {
      parts[[length(parts) + 1L]] <- .sft_bayes_summary_matrix(
        transformed$subject[[measure]][, i, , drop = FALSE][, 1L, ],
        times, "subject", measure, subjects[[i]], report, hdi, rope)
    }
    population <- .sft_bayes_mean_over(transformed$subject[[measure]], 2L,
                                       na.rm = FALSE)
    parts[[length(parts) + 1L]] <- .sft_bayes_summary_matrix(
      population, times, "population", measure, NA_character_, report, hdi, rope)
    if (!is.null(transformed$new_subject)) {
      parts[[length(parts) + 1L]] <- .sft_bayes_summary_matrix(
        transformed$new_subject[[measure]], times, "new_subject", measure,
        NA_character_, report, hdi, rope)
    }
  }
  tidy <- do.call(rbind, parts)
  rownames(tidy) <- NULL
  list(
    tidy = tidy,
    subject = tidy[tidy$Level == "subject", , drop = FALSE],
    population = tidy[tidy$Level == "population", , drop = FALSE],
    new_subject = tidy[tidy$Level == "new_subject", , drop = FALSE]
  )
}


# ---------------------------------------------------------------------------
# Hierarchical SIC and MIC
#
# The salience-split fit already produces a draw-by-draw survivor interaction
# contrast SIC(t) = S_AB_LL - S_AB_LH - S_AB_HL + S_AB_HH per subject.  The
# functions below add (1) a genuinely hierarchical group-level SIC/MIC that
# uses the model's population hyper-parameters rather than averaging the fitted
# subjects, and (2) the mean interaction contrast MIC = integral of SIC(t) at
# every level, with posterior sign tests against zero.
#
# Three group-level estimands are kept distinct, which the original Van Zandt
# fit-each-subject-then-average approach cannot separate:
#   * population            - the typical-subject parameter curve, obtained by
#                             evaluating the population log-hazards with the
#                             subject random effects at zero.  This is the
#                             hierarchical group estimand.
#   * population_finite_mean - the draw-by-draw mean of the (shrunken) observed
#                             subjects; conditions on this sample of subjects.
#   * new_subject           - the posterior-predictive curve for an unobserved
#                             participant; carries the between-subject variance.
# ---------------------------------------------------------------------------

# Cumulative hazard -> survivor for one draws-by-bin log-hazard matrix.
.sft_bayes_survivor_dxJ <- function(eta, widths) {
  eta <- matrix(as.numeric(eta), nrow = dim(eta)[1L])
  exp(-.sft_bayes_cumhaz(eta, widths))
}


# Typical-subject (random effects at zero) population SIC(t), draws by bin.
# Returns NULL when the population hyper-parameters were not extracted, e.g. a
# posterior_draws fit that supplied only the subject log-hazard arrays.
.sft_bayes_population_sic <- function(raw, widths) {
  need <- c("population_A_L", "population_A_H", "population_B_L", "population_B_H",
            "population_delta", "population_gamma")
  if (is.null(raw) || !all(need %in% names(raw))) return(NULL)
  paL <- raw$population_A_L; paH <- raw$population_A_H
  pbL <- raw$population_B_L; pbH <- raw$population_B_H
  pd <- raw$population_delta; g <- raw$population_gamma
  if (length(dim(g)) != 3L || dim(g)[3L] != 4L) return(NULL)
  lse <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
  # gamma columns follow .sft_bayes_salience_cells = c("LL","LH","HL","HH"),
  # matching the eta_AB_* assignments in the Stan model.
  eta_LL <- lse(paL, pbL) + pd + g[, , 1L]
  eta_LH <- lse(paL, pbH) + pd + g[, , 2L]
  eta_HL <- lse(paH, pbL) + pd + g[, , 3L]
  eta_HH <- lse(paH, pbH) + pd + g[, , 4L]
  .sft_bayes_survivor_dxJ(eta_LL, widths) - .sft_bayes_survivor_dxJ(eta_LH, widths) -
    .sft_bayes_survivor_dxJ(eta_HL, widths) + .sft_bayes_survivor_dxJ(eta_HH, widths)
}


# Reduce SIC(t) to the scalar MIC, draw by draw.  Because integration is linear,
# sum_j SIC(t_j) * width_j is identically the contrast of the four posterior cell
# means, mean(LL) - mean(LH) - mean(HL) + mean(HH): each cell mean is the area
# under that cell's survivor on the same grid, so no separate mean computation is
# needed and the MIC is provably the area under the SIC curve we report.
#
# Two approximations are worth stating plainly, because the identity above is
# exact only in the continuum:
#   * The pooled grid stops at the largest observed RT, so each cell mean is a
#     *truncated* mean and the MIC is the truncated-mean contrast.  Mass beyond
#     the last bin is ignored.
#   * SIC(t_j) is evaluated at bin upper edges, so this is a right-endpoint
#     Riemann sum; for a decreasing survivor it understates each cell's area.
# Both biases act in the same direction on all four cells and so largely cancel
# in the contrast, but they do not cancel exactly.  With a heavy right tail or a
# coarse grid, prefer more bins over reading the absolute MIC too literally.
#
# Accepts a draws-by-bin curve (returns a vector) or a draws-by-subject-by-bin
# array (returns a draws-by-subject matrix).
.sft_bayes_integrate_sic <- function(sic, widths) {
  d <- dim(sic)
  if (length(d) == 2L) return(as.numeric(sic %*% widths))
  draws <- d[1L]; I <- d[2L]
  out <- matrix(NA_real_, nrow = draws, ncol = I)
  for (i in seq_len(I)) out[, i] <- as.numeric(matrix(sic[, i, ], nrow = draws) %*% widths)
  out
}


# One-row posterior sign summary for a scalar contrast (MIC or MIC contrast).
.sft_bayes_scalar_sign <- function(v, level, subject, measure, hdi, rope) {
  v <- v[is.finite(v)]
  ci <- .sft_hdi(v, hdi)
  data.frame(Level = level, Subject = subject, Measure = measure,
             Mean = mean(v), Lower = ci[[1L]], Upper = ci[[2L]],
             Prob_positive = mean(v > 0), Prob_negative = mean(v < 0),
             Prob_rope = if (rope > 0) mean(abs(v) <= rope) else NA_real_,
             stringsAsFactors = FALSE)
}


# Assemble the hierarchical MIC posterior (salience-split fits only).
.sft_bayes_mic_result <- function(transformed, raw, widths, subjects, hdi, rope) {
  if (is.null(transformed$SIC)) return(NULL)
  draws <- dim(transformed$SIC)[1L]
  subject_mic <- .sft_bayes_integrate_sic(transformed$SIC, widths)  # draws x I
  colnames(subject_mic) <- subjects
  finite_mean_mic <- rowMeans(subject_mic)
  pop_sic <- .sft_bayes_population_sic(raw, widths)
  population_mic <- if (!is.null(pop_sic)) .sft_bayes_integrate_sic(pop_sic, widths) else NULL
  new_mic <- if (!is.null(transformed$new_subject) && !is.null(transformed$new_subject$SIC))
    .sft_bayes_integrate_sic(transformed$new_subject$SIC, widths) else NULL
  # Headline group draws: the typical-subject parameter when available, else the
  # finite-sample mean.  Contrasts pair each subject with the group within-draw.
  group_mic <- if (!is.null(population_mic)) population_mic else finite_mean_mic
  rows <- list(.sft_bayes_scalar_sign(group_mic, "population", NA_character_, "MIC", hdi, rope),
               .sft_bayes_scalar_sign(finite_mean_mic, "population_finite_mean",
                                      NA_character_, "MIC", hdi, rope))
  if (!is.null(new_mic))
    rows[[length(rows) + 1L]] <- .sft_bayes_scalar_sign(new_mic, "new_subject",
                                                        NA_character_, "MIC", hdi, rope)
  for (i in seq_along(subjects))
    rows[[length(rows) + 1L]] <- .sft_bayes_scalar_sign(subject_mic[, i], "subject",
                                                        subjects[[i]], "MIC", hdi, rope)
  summary <- do.call(rbind, rows); rownames(summary) <- NULL
  contrasts <- do.call(rbind, lapply(seq_along(subjects), function(i)
    .sft_bayes_scalar_sign(subject_mic[, i] - group_mic, "subject_vs_population",
                           subjects[[i]], "MIC_contrast", hdi, rope)))
  rownames(contrasts) <- NULL
  list(summary = summary,
       contrasts = contrasts,
       draws = list(population = population_mic,
                    population_finite_mean = finite_mean_mic,
                    new_subject = new_mic, subject = subject_mic),
       group_estimand = if (!is.null(population_mic))
         "population parameter (typical subject, random effects held at zero)" else
         "finite-sample mean of the observed subjects",
       definition = paste("MIC = integral of posterior SIC(t) over the pooled grid",
                          "= mean(LL) - mean(LH) - mean(HL) + mean(HH)"),
       interpretation = paste("Prob_positive = P(MIC > 0): over-additive",
                             "(parallel-OR / coactive); Prob_negative = P(MIC < 0):",
                             "under-additive (parallel-AND); mass near 0: additive (serial)."))
}


# Assemble the hierarchical SIC(t) curve summaries (salience-split fits only).
.sft_bayes_sic_result <- function(transformed, raw, grid, subjects, hdi, rope) {
  if (is.null(transformed$SIC)) return(NULL)
  draws <- dim(transformed$SIC)[1L]
  times <- grid$bins$upper; report <- grid$bins$report
  parts <- list()
  add_curve <- function(x, level, subject = NA_character_) {
    parts[[length(parts) + 1L]] <<- .sft_bayes_summary_matrix(
      x, times, level, "SIC", subject, report, hdi, rope, "SIC")
  }
  pop_sic <- .sft_bayes_population_sic(raw, grid$bins$width)
  if (!is.null(pop_sic)) add_curve(pop_sic, "population")
  add_curve(transformed$population$SIC, "population_finite_mean")
  if (!is.null(transformed$new_subject) && !is.null(transformed$new_subject$SIC))
    add_curve(transformed$new_subject$SIC, "new_subject")
  for (i in seq_along(subjects))
    add_curve(matrix(transformed$subject$SIC[, i, ], nrow = draws), "subject", subjects[[i]])
  tidy <- do.call(rbind, parts); rownames(tidy) <- NULL
  list(summary = tidy,
       subject = tidy[tidy$Level == "subject", , drop = FALSE],
       population = tidy[tidy$Level %in% c("population", "population_finite_mean"), , drop = FALSE],
       new_subject = tidy[tidy$Level == "new_subject", , drop = FALSE],
       group_estimand = if (!is.null(pop_sic))
         "population parameter (typical subject)" else
         "finite-sample mean of the observed subjects",
       note = paste("Hierarchical SIC: subject curves are partially pooled toward the",
                    "population; the population level is a model parameter (not an average",
                    "of point estimates), and new_subject is the posterior-predictive SIC",
                    "for an unobserved participant.  Prob_super = P(SIC>0), Prob_limited = P(SIC<0)."))
}


.sft_bayes_effects_from_posterior <- function(raw, priors, I, J, subjects, widths,
                                              seed, new_subject = TRUE) {
  needed <- c("eta_A", "eta_B", "eta_AB")
  if (!all(needed %in% names(raw))) stop("The Stan fit did not return all log-hazard arrays.",
                                          call. = FALSE)
  draws <- dim(raw$eta_A)[1L]
  log_hazard <- list(A = .sft_bayes_dim_array(raw$eta_A, draws, I, J, "eta_A"),
                     B = .sft_bayes_dim_array(raw$eta_B, draws, I, J, "eta_B"),
                     AB = .sft_bayes_dim_array(raw$eta_AB, draws, I, J, "eta_AB"))
  new_eta <- NULL
  if (isTRUE(new_subject)) {
    if (!all(c("population_A", "population_B", "population_delta",
               "sigma_speed", "sigma_asymmetry", "sigma_capacity") %in% names(raw))) {
      warning("The fit lacks population effect draws; new-subject prediction was omitted.", call. = FALSE)
    } else {
      set.seed(seed)
      # Draw a fresh subject with the same scaling as the Stan model, where each
      # effect is <fixed prior sd> * sigma_* * z with z ~ N(0, 1). Omitting the
      # fixed prior sd would over-disperse the predicted subject.
      sigma_speed <- as.numeric(raw$sigma_speed); sigma_asymmetry <- as.numeric(raw$sigma_asymmetry)
      sigma_capacity <- as.numeric(raw$sigma_capacity)
      speed <- stats::rnorm(draws, 0, priors$speed_sd * sigma_speed)
      asym <- stats::rnorm(draws, 0, priors$asymmetry_sd * sigma_asymmetry)
      cap <- stats::rnorm(draws, 0, priors$capacity_sd * sigma_capacity)
      pa <- raw$population_A; pb <- raw$population_B; pd <- raw$population_delta
      new_a <- pa + speed + asym / 2
      new_b <- pb + speed - asym / 2
      # AB centres on the UCIP log-sum-exp plus the population capacity deviation
      # (pd) and the subject capacity shift, mirroring eta_AB in the Stan model.
      new_eta <- list(A = new_a, B = new_b,
                      AB = pmax(new_a, new_b) + log1p(exp(-abs(new_a - new_b))) + pd + cap)
    }
  }
  list(log_hazard = log_hazard, new_log_hazard = new_eta)
}


.sft_bayes_effects_from_salience_posterior <- function(raw, priors, I, J,
                                                       subjects, seed,
                                                       new_subject = TRUE) {
  eta_names <- paste0("eta_", .sft_bayes_salience_series)
  if (!all(eta_names %in% names(raw))) {
    stop("The salience-split fit did not return all eight log-hazard arrays.", call. = FALSE)
  }
  first_name <- eta_names[[1L]]
  first_eta <- raw[[first_name]]
  draws <- dim(first_eta)[1L]
  log_hazard <- setNames(lapply(seq_along(.sft_bayes_salience_series), function(k) {
    nm <- eta_names[[k]]
    .sft_bayes_dim_array(raw[[nm]], draws, I, J, nm)
  }), .sft_bayes_salience_series)
  new_eta <- NULL
  if (isTRUE(new_subject)) {
    needed <- c("population_A_L", "population_A_H", "population_B_L", "population_B_H",
                "population_delta", "sigma_speed", "sigma_asymmetry", "sigma_capacity")
    if (!all(needed %in% names(raw))) {
      warning("The fit lacks population effect draws; new-subject prediction was omitted.",
              call. = FALSE)
    } else {
      set.seed(seed)
      sigma_speed <- as.numeric(raw$sigma_speed)
      sigma_asymmetry <- as.numeric(raw$sigma_asymmetry)
      sigma_capacity <- as.numeric(raw$sigma_capacity)
      speed <- stats::rnorm(draws, 0, priors$speed_sd * sigma_speed)
      asym <- stats::rnorm(draws, 0, priors$asymmetry_sd * sigma_asymmetry)
      cap <- stats::rnorm(draws, 0, priors$capacity_sd * sigma_capacity)
      new_eta <- list(
        A_L = raw$population_A_L + speed + asym / 2,
        A_H = raw$population_A_H + speed + asym / 2,
        B_L = raw$population_B_L + speed - asym / 2,
        B_H = raw$population_B_H + speed - asym / 2
      )
      gamma <- matrix(0, nrow = draws, ncol = J)
      if ("population_gamma" %in% names(raw)) {
        if (length(dim(raw$population_gamma)) != 3L ||
            !all(dim(raw$population_gamma) == c(draws, J, 4L))) {
          stop("Unexpected population_gamma dimensions in the Stan fit.", call. = FALSE)
        }
      }
      for (cell in .sft_bayes_salience_cells) {
        if ("population_gamma" %in% names(raw)) {
          cc <- match(cell, .sft_bayes_salience_cells)
          gamma <- matrix(raw$population_gamma[, , cc, drop = TRUE], nrow = draws, ncol = J)
        } else gamma[,] <- 0
        u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
        a <- new_eta[[paste0("A_", u)]]
        b <- new_eta[[paste0("B_", v)]]
        lse <- pmax(a, b) + log1p(exp(-abs(a - b)))
        new_eta[[paste0("AB_", cell)]] <- lse + raw$population_delta + cap + gamma
      }
    }
  }
  list(log_hazard = log_hazard, new_log_hazard = new_eta)
}


.sft_bayes_prior_draws <- function(prepared, n_draws, seed) {
  n_draws <- max(1L, as.integer(n_draws)); J <- prepared$n_bins; K <- ncol(prepared$basis)
  B <- prepared$basis; p <- prepared$priors
  set.seed(seed)
  pa <- pb <- pd <- matrix(NA_real_, nrow = n_draws, ncol = J)
  for (m in seq_len(n_draws)) {
    beta_a <- stats::rnorm(K, p$log_hazard_mean, p$log_hazard_sd)
    beta_b <- stats::rnorm(K, p$log_hazard_mean, p$log_hazard_sd)
    beta_d <- stats::rnorm(K, p$delta_mean, p$delta_sd)
    if (K >= 3L) {
      # Second-difference random walk: beta[k] - 2 beta[k-1] + beta[k-2] ~ N(0, sd).
      # Draw it as a cumulative recursion so the prior predictive matches the
      # smoothing prior used by the Stan model.
      for (k in 3:K) {
        beta_a[k] <- 2 * beta_a[k - 1L] - beta_a[k - 2L] + stats::rnorm(1L, 0, p$smooth_A_sd)
        beta_b[k] <- 2 * beta_b[k - 1L] - beta_b[k - 2L] + stats::rnorm(1L, 0, p$smooth_B_sd)
        beta_d[k] <- 2 * beta_d[k - 1L] - beta_d[k - 2L] + stats::rnorm(1L, 0, p$smooth_delta_sd)
      }
    }
    pa[m, ] <- as.numeric(B %*% beta_a)
    pb[m, ] <- as.numeric(B %*% beta_b)
    pd[m, ] <- as.numeric(B %*% beta_d)
  }
  speed <- stats::rnorm(n_draws, 0, p$speed_sd)
  asym <- stats::rnorm(n_draws, 0, p$asymmetry_sd)
  cap <- stats::rnorm(n_draws, 0, p$capacity_sd)
  a <- pa + speed + asym / 2; b <- pb + speed - asym / 2
  ab <- pmax(a, b) + log1p(exp(-abs(a - b))) + pd + cap
  list(A = a, B = b, AB = ab)
}


.sft_bayes_prior_predictive <- function(prepared, n_draws, seed) {
  raw <- .sft_bayes_prior_draws(prepared, n_draws, seed)
  widths <- prepared$grid$bins$width; times <- prepared$grid$bins$upper
  H <- lapply(raw, .sft_bayes_cumhaz, widths = widths)
  D <- H$AB - H$A - H$B; C <- H$AB / (H$A + H$B)
  C[!is.finite(C)] <- NA_real_
  counts <- lapply(seq_along(.sft_bayes_series), function(s) {
    E <- matrix(prepared$exposure[, s, ], nrow = prepared$n_subjects,
                ncol = prepared$n_bins)
    Hname <- .sft_bayes_series[[s]]
    rates <- raw[[Hname]]
    total_mean <- vapply(seq_len(nrow(rates)), function(m) {
      sum(E * matrix(exp(rates[m, ]) * prepared$grid$bins$width,
                     nrow = nrow(E), ncol = ncol(E), byrow = TRUE))
    }, numeric(1))
    stats::rpois(nrow(rates), lambda = pmax(total_mean, .Machine$double.eps))
  })
  names(counts) <- .sft_bayes_series
  list(log_hazard = raw, H_A = H$A, H_B = H$B, H_AB = H$AB, D = D, C = C,
       times = times, counts = counts)
}


.sft_bayes_prior_draws_salience <- function(prepared, n_draws, seed) {
  n_draws <- max(1L, as.integer(n_draws)); J <- prepared$n_bins
  I <- prepared$n_subjects; K <- ncol(prepared$basis); B <- prepared$basis; p <- prepared$priors
  set.seed(seed)
  out <- lapply(.sft_bayes_salience_series,
                function(nm) array(0, dim = c(n_draws, I, J)))
  names(out) <- .sft_bayes_salience_series
  for (m in seq_len(n_draws)) {
    draw_curve <- function(mean, sd, smooth_sd) {
      beta <- stats::rnorm(K, mean, sd)
      if (K >= 3L) for (k in 3:K)
        beta[k] <- 2 * beta[k - 1L] - beta[k - 2L] + stats::rnorm(1L, 0, smooth_sd)
      as.numeric(B %*% beta)
    }
    pa_l <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_A_sd)
    pa_h <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_A_sd)
    pb_l <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_B_sd)
    pb_h <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_B_sd)
    pd <- draw_curve(p$delta_mean, p$delta_sd, p$smooth_delta_sd)
    gamma_beta <- matrix(stats::rnorm(K * 4L, 0, p$gamma_sd), nrow = K, ncol = 4L)
    if (K >= 3L) for (k in 3:K) for (cc in 1:4)
      gamma_beta[k, cc] <- 2 * gamma_beta[k - 1L, cc] - gamma_beta[k - 2L, cc] +
        stats::rnorm(1L, 0, p$smooth_gamma_sd)
    gamma_beta <- gamma_beta - rowMeans(gamma_beta)
    gamma <- B %*% gamma_beta
    speed <- stats::rnorm(I, 0, p$speed_sd)
    asym <- stats::rnorm(I, 0, p$asymmetry_sd)
    cap <- stats::rnorm(I, 0, p$capacity_sd)
    for (i in seq_len(I)) {
      a_l <- pa_l + speed[[i]] + asym[[i]] / 2
      a_h <- pa_h + speed[[i]] + asym[[i]] / 2
      b_l <- pb_l + speed[[i]] - asym[[i]] / 2
      b_h <- pb_h + speed[[i]] - asym[[i]] / 2
      out$A_L[m, i, ] <- a_l; out$A_H[m, i, ] <- a_h
      out$B_L[m, i, ] <- b_l; out$B_H[m, i, ] <- b_h
      out$AB_LL[m, i, ] <- pmax(a_l, b_l) + log1p(exp(-abs(a_l - b_l))) + pd + cap[[i]] + gamma[, 1L]
      out$AB_LH[m, i, ] <- pmax(a_l, b_h) + log1p(exp(-abs(a_l - b_h))) + pd + cap[[i]] + gamma[, 2L]
      out$AB_HL[m, i, ] <- pmax(a_h, b_l) + log1p(exp(-abs(a_h - b_l))) + pd + cap[[i]] + gamma[, 3L]
      out$AB_HH[m, i, ] <- pmax(a_h, b_h) + log1p(exp(-abs(a_h - b_h))) + pd + cap[[i]] + gamma[, 4L]
    }
  }
  out
}


.sft_bayes_predictive_from_series <- function(prepared, log_hazard, seed,
                                              n_draws = NULL) {
  set.seed(seed)
  all_draws <- dim(log_hazard[[1L]])[1L]
  if (!is.null(n_draws) && length(n_draws) == 1L && is.finite(n_draws) && n_draws >= 1L) {
    n_draws <- min(all_draws, as.integer(n_draws))
    keep <- if (n_draws < all_draws) sample(seq_len(all_draws), n_draws) else seq_len(all_draws)
  } else keep <- seq_len(all_draws)
  predicted <- lapply(prepared$series, function(nm) {
    eta <- log_hazard[[nm]][keep, , , drop = FALSE]
    E <- prepared$exposure[, nm, , drop = FALSE]
    E <- matrix(E, nrow = prepared$n_subjects, ncol = prepared$n_bins)
    out <- array(0L, dim = c(length(keep), prepared$n_subjects, prepared$n_bins))
    for (m in seq_along(keep)) {
      rate <- exp(eta[m, , ])
      rate <- matrix(rate, nrow = prepared$n_subjects, ncol = prepared$n_bins)
      out[m, , ] <- matrix(stats::rpois(length(E), as.numeric(E * rate)),
                            nrow = prepared$n_subjects, ncol = prepared$n_bins)
    }
    out
  })
  names(predicted) <- prepared$series
  summary <- do.call(rbind, lapply(prepared$series, function(nm) {
    vals <- apply(predicted[[nm]], 1L, sum)
    interval <- .sft_hdi(vals)
    s <- match(nm, prepared$series)
    data.frame(Series = nm, Observed = sum(prepared$events[, s, ]),
               Mean = mean(vals), Lower = interval[[1L]], Upper = interval[[2L]],
               stringsAsFactors = FALSE)
  }))
  list(draws = predicted, summary = summary)
}


.sft_bayes_prior_predictive_salience <- function(prepared, n_draws, seed) {
  raw <- .sft_bayes_prior_draws_salience(prepared, n_draws, seed)
  H <- lapply(raw, .sft_bayes_cumhaz, widths = prepared$grid$bins$width)
  cells <- .sft_bayes_salience_cells
  Dcells <- array(NA_real_, dim = c(dim(H[[1L]]), 4L))
  Ccells <- Dcells
  for (cc in seq_along(cells)) {
    cell <- cells[[cc]]; u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
    Dcells[, , , cc] <- H[[paste0("AB_", cell)]] - H[[paste0("A_", u)]] - H[[paste0("B_", v)]]
    z <- H[[paste0("AB_", cell)]] /
      (H[[paste0("A_", u)]] + H[[paste0("B_", v)]])
    z[!is.finite(z)] <- NA_real_; Ccells[, , , cc] <- z
  }
  sic <- exp(-H$AB_LL) - exp(-H$AB_LH) - exp(-H$AB_HL) + exp(-H$AB_HH)
  list(log_hazard = raw, H = H, D_cells = Dcells, C_cells = Ccells,
       D = .sft_bayes_mean_over(Dcells, 4L),
       C = exp(.sft_bayes_mean_over(log(Ccells), 4L)),
       SIC = sic, times = prepared$grid$bins$upper,
       counts = .sft_bayes_predictive_from_series(prepared, raw, seed + 7L)$summary)
}


.sft_bayes_diagnostics <- function(fit) {
  sm <- rstan::summary(fit)$summary
  pars <- rownames(sm)
  convergence <- data.frame(parameter = pars, mean = sm[, "mean"], sd = sm[, "sd"],
                            rhat = sm[, "Rhat"], n_eff = sm[, "n_eff"], stringsAsFactors = FALSE)
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth <- max(vapply(sampler, function(x) max(x[, "treedepth__"]), numeric(1)))
  energy <- unlist(lapply(sampler, function(x) x[, "energy__"]))
  ebfmi <- if (length(energy) > 1L) mean(diff(energy)^2) / stats::var(energy) else NA_real_
  list(convergence = convergence, rhat = convergence$rhat, n_eff = convergence$n_eff,
       divergent = divergent, divergences = divergent, max_treedepth = treedepth,
       ebfmi = ebfmi,
       ok = !any(convergence$rhat > 1.01, na.rm = TRUE) && divergent == 0L)
}


.sft_bayes_posterior_predictive <- function(prepared, transformed, seed, n_draws = NULL) {
  set.seed(seed)
  all_draws <- dim(transformed$H_A)[1L]
  if (!is.null(n_draws) && length(n_draws) == 1L && is.finite(n_draws) && n_draws >= 1L) {
    n_draws <- min(all_draws, as.integer(n_draws))
    keep <- if (n_draws < all_draws) sample(seq_len(all_draws), n_draws) else seq_len(all_draws)
  } else {
    keep <- seq_len(all_draws)
  }
  draws <- length(keep)
  predicted <- lapply(seq_along(.sft_bayes_series), function(s) {
    # The likelihood is written in terms of hazards.  Recover interval rates
    # from cumulative hazards before simulating replicated event counts.
    hname <- .sft_bayes_series[[s]]
    Hs <- transformed[[paste0("H_", hname)]][keep, , , drop = FALSE]
    rates <- Hs
    rates[, , 1L] <- Hs[, , 1L] / prepared$grid$bins$width[[1L]]
    J <- dim(Hs)[3L]
    if (J > 1L) for (j in 2:J) {
      rates[, , j] <- (Hs[, , j] - Hs[, , j - 1L]) / prepared$grid$bins$width[[j]]
    }
    E <- matrix(prepared$exposure[, match(hname, .sft_bayes_series), ],
                nrow = prepared$n_subjects, ncol = J)
    out <- array(0L, dim = c(draws, nrow(E), J))
    for (m in seq_len(draws)) out[m, , ] <- matrix(
      stats::rpois(length(E), as.numeric(E * matrix(rates[m, , ],
                                                    nrow = nrow(E), ncol = J))),
      nrow = nrow(E), ncol = J)
    out
  })
  names(predicted) <- .sft_bayes_series
  summary <- do.call(rbind, lapply(seq_along(predicted), function(s) {
    vals <- apply(predicted[[s]], 1L, sum)
    interval <- .sft_hdi(vals)
    data.frame(Series = .sft_bayes_series[[s]], Observed = sum(prepared$events[, s, ]),
               Mean = mean(vals), Lower = interval[[1L]], Upper = interval[[2L]],
               stringsAsFactors = FALSE)
  }))
  list(draws = predicted, summary = summary)
}


.sft_bayes_empty_transformed <- function() {
  list(H_A = NULL, H_B = NULL, H_AB = NULL, D = NULL, C = NULL,
       subject = NULL, population = NULL, new_subject = NULL)
}


.sft_bayes_empty_curve_data <- function() {
  data.frame(Level = character(), Subject = character(), Series = character(),
             Measure = character(), Time = numeric(), Ct = numeric(),
             Mean = numeric(), Lower = numeric(), Upper = numeric(),
             Prob_super = numeric(), Prob_limited = numeric(), Prob_rope = numeric(),
             Report = logical(), stringsAsFactors = FALSE)
}

