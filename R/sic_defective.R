# Unconditional ("defective") interaction contrasts for correct and incorrect
# responses, following Gondan and Legner's generalisation of Townsend and
# Nozawa (1995) to tasks with non-negligible error rates.
#
# The whole extension turns on one estimator. The classical SIC route conditions
# on accuracy: errors are discarded and each cell's ECDF is normalised by its
# own count of correct trials. The unconditional route keeps every trial and
# normalises by the cell total, giving the two subdistributions
#
#   F+(t) = P(RT <= t and correct),   F-(t) = P(RT <= t and error),
#
# neither of which reaches 1; F+(Inf) is the cell accuracy. Their sum is the
# accuracy-agnostic distribution F = F+ + F-.
#
# Sign convention. Gondan and Legner state their results for
# D2F = F_AB + F_ab - F_Ab - F_aB with uppercase denoting high salience. This
# package's SIC is the survivor contrast, SIC = -D2F, so every sign in their
# tables flips here. A race model has D2F+ <= 0, hence SIC+ >= 0, recovering the
# familiar positive parallel-OR SIC; AND parallel-exhaustive has D2F+ >= 0,
# hence SIC+ <= 0, recovering the familiar negative parallel-AND SIC.


# Normalise the `errors` split. "discard" is the conditional, correct-only
# analysis the package has always done; "defective" is the unconditional route.
.sft_errors_method <- function(method) {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(method[[1L]])))
  switch(key,
         discard = "discard", conditional = "discard", correctonly = "discard",
         correct = "discard", drop = "discard", exclude = "discard",
         defective = "defective", unconditional = "defective",
         subdistribution = "defective", subdistributions = "defective",
         sub = "defective",
         stop("errors must be one of 'discard' or 'defective'.", call. = FALSE))
}


# Defective ECDF. Trials that did not produce the response type of interest are
# pushed to +Inf so they stay in the denominator but never contribute at any
# finite t -- the trick used in Gondan and Legner's own supplement. Non-finite
# RTs are treated as omissions and handled the same way, which is also what the
# go/no-go case requires.
.sft_subecdf <- function(rt, keep) {
  rt <- as.numeric(rt)
  rt[!is.finite(rt)] <- Inf
  rt[!keep] <- Inf
  if (!length(rt)) return(.sft_zero_curve(0))
  if (!any(is.finite(rt))) return(.sft_zero_curve(0))
  stats::ecdf(rt)
}


# Normalise the four SIC cells to a common (rt, correct, n) representation.
# `CR` may be NULL, a list of four correctness vectors in HH, HL, LH, LL order,
# or a list named with those cells. Under "discard" the result holds only the
# finite correct RTs, so downstream code sees exactly what it saw before.
.sft_sic_cells <- function(HH, HL, LH, LL, CR = NULL, errors = "discard") {
  rts <- lapply(list(HH = HH, HL = HL, LH = LH, LL = LL), as.numeric)
  if (is.null(CR)) {
    if (errors == "defective") {
      stop("errors = 'defective' needs correctness indicators; supply CR as a ",
           "list of four vectors in HH, HL, LH, LL order.", call. = FALSE)
    }
    cr <- lapply(rts, function(x) rep(TRUE, length(x)))
  } else {
    if (!is.list(CR) || length(CR) != 4L) {
      stop("CR must be a list of four correctness vectors in HH, HL, LH, LL order.",
           call. = FALSE)
    }
    if (!is.null(names(CR)) && all(c("HH", "HL", "LH", "LL") %in% names(CR))) {
      CR <- CR[c("HH", "HL", "LH", "LL")]
    }
    lens <- vapply(CR, length, integer(1))
    if (!identical(unname(lens), unname(vapply(rts, length, integer(1))))) {
      stop("Each CR element must have the same length as its response-time cell.",
           call. = FALSE)
    }
    cr <- Map(function(x, k) .to_correct_indicator(k, length(x)), rts, CR)
  }
  names(cr) <- names(rts)
  cells <- Map(function(rt, keep) {
    finite <- is.finite(rt)
    if (errors == "discard") {
      idx <- finite & keep
      list(rt = rt[idx], correct = rep(TRUE, sum(idx)), n = sum(idx))
    } else {
      rt[!finite] <- Inf
      list(rt = rt, correct = keep, n = length(rt))
    }
  }, rts, cr)
  names(cells) <- names(rts)
  if (any(vapply(cells, function(z) z$n, integer(1)) == 0L)) {
    stop("All four SIC cells must contain at least one usable trial.", call. = FALSE)
  }
  cells
}


# Per-cell subdistribution functions, plus the pooled finite support.
.sft_sic_subcdfs <- function(cells) {
  list(correct = lapply(cells, function(z) .sft_subecdf(z$rt, z$correct)),
       error = lapply(cells, function(z) .sft_subecdf(z$rt, !z$correct)),
       all = lapply(cells, function(z) .sft_subecdf(z$rt, rep(TRUE, z$n))))
}


.sft_sic_times <- function(cells) {
  sort(unique(unlist(lapply(cells, function(z) z$rt[is.finite(z$rt)]),
                     use.names = FALSE)))
}


# The interaction contrast in the package's survivor convention,
# SIC = F_LH + F_HL - F_HH - F_LL = -D2F.
.sft_sic_contrast <- function(fns, times) {
  fns$LH(times) + fns$HL(times) - fns$HH(times) - fns$LL(times)
}


# Assemble the three unconditional contrasts, the per-cell accuracies, and the
# large-t asymptotes. SIC = SIC.correct + SIC.error holds exactly at every t.
.sft_sic_defective <- function(cells, interpolate = FALSE) {
  times <- .sft_sic_times(cells)
  if (!length(times)) stop("No finite response times in the SIC cells.", call. = FALSE)
  cdfs <- .sft_sic_subcdfs(cells)
  values <- list(correct = .sft_sic_contrast(cdfs$correct, times),
                 error = .sft_sic_contrast(cdfs$error, times),
                 all = .sft_sic_contrast(cdfs$all, times))
  as_fun <- function(y) {
    if (interpolate) .sft_curve(times, y) else stats::stepfun(times, c(0, y))
  }
  accuracy <- vapply(cells, function(z) sum(z$correct & is.finite(z$rt)) / z$n, numeric(1))
  omissions <- vapply(cells, function(z) sum(!is.finite(z$rt)) / z$n, numeric(1))
  last <- length(times)
  list(times = times,
       SIC = as_fun(values$all), SIC.correct = as_fun(values$correct),
       SIC.error = as_fun(values$error),
       values = values, cdf = cdfs, accuracy = accuracy, omissions = omissions,
       # Asymptotes, evaluated at the largest observed finite time. In the
       # survivor convention these are -D2P+ and -D2P-; with no omissions they
       # sum to zero because the agnostic contrast must vanish as t grows.
       asymptote = c(correct = values$correct[[last]], error = values$error[[last]],
                     all = values$all[[last]]),
       n = vapply(cells, function(z) z$n, integer(1)))
}


# ---------------------------------------------------------------------------
# Stochastic dominance for the extended assumption set (Gondan & Legner 8A-8D)
# ---------------------------------------------------------------------------

# The four selective-influence cell pairs, expressed once so the survivor,
# correct-subdistribution, and error-subdistribution blocks stay aligned.
.sft_sic_dominance_pairs <- list(
  c("hh", "hl"), c("hh", "lh"), c("hl", "ll"), c("lh", "ll")
)


# Two-sample KS on defective samples. The shared atom at +Inf contributes
# nothing to the supremum, so the statistic is the sup over finite t of the
# subdistribution difference; the asymptotic p-value is computed from the
# continuous-case Kolmogorov law and is therefore conservative here.
.sft_sub_ks <- function(x, xkeep, y, ykeep, alternative) {
  ax <- as.numeric(x); ax[!is.finite(ax)] <- Inf; ax[!xkeep] <- Inf
  ay <- as.numeric(y); ay[!is.finite(ay)] <- Inf; ay[!ykeep] <- Inf
  suppressWarnings(stats::ks.test(ax, ay, alternative = alternative, exact = FALSE))
}


# Dirichlet dominance for a pair of subdistributions. Mass that is not a finite
# response of the requested type is collected in one extra category, so the two
# cells are compared on the same probability scale rather than after
# renormalisation.
.sft_bayes_pair_sub <- function(x, xkeep, y, ykeep, nx, ny, nbin = 20L,
                                nsamp = 10000L, tol = 0, seed = NULL) {
  x <- as.numeric(x)[xkeep & is.finite(x)]
  y <- as.numeric(y)[ykeep & is.finite(y)]
  pooled <- c(x, y)
  if (!length(pooled)) {
    return(list(BF_greater = 1, BF_less = 1, posterior_greater = .5,
                posterior_less = .5, prior_greater = .5, prior_less = .5,
                nbin = 0L, nsamp = nsamp))
  }
  nbin <- max(1L, min(as.integer(nbin), length(unique(pooled))))
  breaks <- unique(stats::quantile(pooled, probs = seq(0, 1, length.out = nbin + 1L),
                                   names = FALSE, type = 8))
  if (length(breaks) < 2L) breaks <- c(breaks[[1L]], breaks[[1L]] + 1)
  breaks[length(breaks)] <- max(breaks[length(breaks)], max(pooled) + .Machine$double.eps)
  nbin <- length(breaks) - 1L
  bin_counts <- function(v) {
    if (!length(v)) return(rep(0, nbin))
    as.numeric(hist(v, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts)
  }
  # The trailing category is "not a finite response of this type": errors,
  # omissions, or the other response alternative.
  cx <- c(bin_counts(x), nx - length(x))
  cy <- c(bin_counts(y), ny - length(y))
  alpha <- rep(1 / (nbin + 1L), nbin + 1L)
  old_seed <- .sft_stash_seed(seed)
  on.exit(.sft_unstash_seed(old_seed), add = TRUE)
  draw <- function(a) .sft_rdirichlet(nsamp, a)
  # Dominance is judged on the first nbin categories only, i.e. on the
  # subdistribution over finite response times.
  dominates <- function(a, b) {
    d <- (a - b)[, seq_len(nbin), drop = FALSE]
    cum <- if (nbin == 1L) d else t(apply(d, 1L, cumsum))
    rowSums(cum >= -tol) == nbin
  }
  prior_x <- draw(alpha); prior_y <- draw(alpha)
  post_x <- draw(cx + alpha); post_y <- draw(cy + alpha)
  prior_greater <- mean(dominates(prior_x, prior_y))
  post_greater <- mean(dominates(post_x, post_y))
  prior_less <- mean(dominates(prior_y, prior_x))
  post_less <- mean(dominates(post_y, post_x))
  eps <- .Machine$double.eps
  list(BF_greater = post_greater / max(prior_greater, eps),
       BF_less = post_less / max(prior_less, eps),
       posterior_greater = post_greater, posterior_less = post_less,
       prior_greater = prior_greater, prior_less = prior_less,
       nbin = nbin, nsamp = nsamp)
}


# The full 8A/8B/8D dominance battery on defective cells. Each block reports
# both directions of each pair, matching the layout of the correct-only route.
.sft_dominance_defective <- function(cells, method = "ks", nbin = 20L,
                                     nsamp = 10000L, seed = NULL) {
  # Every comparison is run on distribution functions, but 8B is *labelled* in
  # survivor terms to match the correct-only route. Because S = 1 - F, the
  # survivor label reads in the opposite direction to the F-test that produces
  # it: predicting F_first >= F_second is predicting S_first <= S_second.
  blocks <- list(
    list(assumption = "8B (survivor ordering)", label = "S",
         keep = function(z) rep(TRUE, z$n), flip_label = TRUE, predict_ge = TRUE),
    list(assumption = "8A (correct subdistribution)", label = "F+",
         keep = function(z) z$correct, flip_label = FALSE, predict_ge = TRUE),
    list(assumption = "8D (error subdistribution)", label = "F-",
         keep = function(z) !z$correct, flip_label = FALSE, predict_ge = FALSE)
  )
  cell_of <- c(hh = "HH", hl = "HL", lh = "LH", ll = "LL")
  rows <- list(); k <- 0L
  for (b in blocks) for (p in .sft_sic_dominance_pairs) {
    a <- cells[[cell_of[[p[[1L]]]]]]; d <- cells[[cell_of[[p[[2L]]]]]]
    ka <- b$keep(a); kd <- b$keep(d)
    # `alt_pred` is the ks.test alternative that matches the assumption: "greater"
    # means the first cell's distribution function lies above the second's.
    alt_pred <- if (b$predict_ge) "greater" else "less"
    alt_rev <- if (b$predict_ge) "less" else "greater"
    sym <- function(alt) {
      ge <- identical(alt, "greater")
      if (b$flip_label) ge <- !ge
      if (ge) " >= " else " <= "
    }
    lab <- function(alt) paste0(b$label, ".", p[[1L]], sym(alt), b$label, ".", p[[2L]])
    k <- k + 1L
    if (method == "ks") {
      pred <- .sft_sub_ks(a$rt, ka, d$rt, kd, alt_pred)
      rev <- .sft_sub_ks(a$rt, ka, d$rt, kd, alt_rev)
      rows[[k]] <- data.frame(
        Assumption = b$assumption, Test = c(lab(alt_pred), lab(alt_rev)),
        statistic = c(unname(pred$statistic), unname(rev$statistic)),
        p.value = c(pred$p.value, rev$p.value),
        Predicted = c(TRUE, FALSE), stringsAsFactors = FALSE)
    } else {
      bf <- .sft_bayes_pair_sub(a$rt, ka, d$rt, kd, a$n, d$n, nbin = nbin,
                                nsamp = nsamp,
                                seed = if (is.null(seed)) NULL else seed + k)
      # BF_greater is the evidence that the first cell's distribution function
      # dominates the second's, i.e. the "greater" alternative above.
      take <- function(alt) if (identical(alt, "greater")) "greater" else "less"
      rows[[k]] <- data.frame(
        Assumption = b$assumption, Test = c(lab(alt_pred), lab(alt_rev)),
        BF = c(bf[[paste0("BF_", take(alt_pred))]], bf[[paste0("BF_", take(alt_rev))]]),
        posterior = c(bf[[paste0("posterior_", take(alt_pred))]],
                      bf[[paste0("posterior_", take(alt_rev))]]),
        prior = c(bf[[paste0("prior_", take(alt_pred))]],
                  bf[[paste0("prior_", take(alt_rev))]]),
        p.value = c(bf[[paste0("posterior_", take(alt_pred))]],
                    bf[[paste0("posterior_", take(alt_rev))]]),
        Predicted = c(TRUE, FALSE), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  if (method != "ks") {
    attr(out, "inferential_type") <-
      "Bayesian posterior probability and posterior/prior event-probability ratio"
  }
  out
}


# ---------------------------------------------------------------------------
# Shared RNG and Dirichlet helpers
# ---------------------------------------------------------------------------

.sft_stash_seed <- function(seed) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  if (!is.null(seed)) set.seed(seed)
  old
}


.sft_unstash_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else assign(".Random.seed", old, envir = .GlobalEnv)
}


.sft_rdirichlet <- function(n, a) {
  g <- matrix(stats::rgamma(n * length(a), shape = rep(a, each = n)), nrow = n)
  g / rowSums(g)
}


# ---------------------------------------------------------------------------
# Dirichlet posterior over the three unconditional contrasts
# ---------------------------------------------------------------------------

# Sign class of one contrast curve. This is deliberately *not* the six-way SFT
# architecture classification: with errors retained, the architecture mapping
# needs the task type (detection / OR / AND) and admits substantial mimicry, so
# the defective route reports the sign of each curve and stops there.
.sft_sic_sign_class <- function(x, tol) {
  pos <- any(x > tol); neg <- any(x < -tol)
  if (!pos && !neg) return("Zero")
  if (!pos) return("Negative")
  if (!neg) return("Positive")
  "Crosses"
}


.sft_sic_sign_levels <- c("Zero", "Negative", "Positive", "Crosses")


# Dirichlet posterior for the defective SIC. Each cell carries one Dirichlet
# over 2 * (nbin + 1) categories: nbin correct bins, one "correct but never
# observed" catch-all, then the same for errors. Splitting the concentration
# evenly between the two halves makes the implied prior on cell accuracy
# uniform, which a single flat Dirichlet over all categories would not do.
.sft_sic_dp_defective <- function(cells, nbin = NULL, nsamp = 10000L,
                                  maxn = 500000L, tolSIC = 1e-2, ci = .95,
                                  seed = NULL) {
  pooled <- unlist(lapply(cells, function(z) z$rt[is.finite(z$rt)]), use.names = FALSE)
  if (!length(pooled)) stop("No finite response times in the SIC cells.", call. = FALSE)
  if (is.null(nbin)) nbin <- floor(min(vapply(cells, function(z) z$n, integer(1))) / 2)
  nbin <- max(1L, min(as.integer(nbin), length(unique(pooled))))
  bins <- unique(stats::quantile(pooled, probs = seq(0, 1, length.out = nbin + 1L),
                                 names = FALSE, type = 8))
  if (length(bins) < 2L) bins <- c(bins[[1L]], bins[[1L]] + 1)
  bins[length(bins)] <- max(bins[length(bins)], max(pooled) + .Machine$double.eps)
  nbin <- length(bins) - 1L
  K <- nbin + 1L
  cell_counts <- function(z) {
    ok <- is.finite(z$rt)
    cnt <- function(v) if (!length(v)) rep(0, nbin) else
      as.numeric(hist(v, breaks = bins, plot = FALSE, include.lowest = TRUE)$counts)
    correct <- cnt(z$rt[ok & z$correct])
    error <- cnt(z$rt[ok & !z$correct])
    c(correct, sum(z$correct) - sum(correct), error, sum(!z$correct) - sum(error))
  }
  counts <- vapply(cells, cell_counts, numeric(2L * K))
  # Half the unit concentration to each response type, spread over its K
  # categories: the marginal prior on accuracy is then Beta(1, 1).
  alpha <- rep(1 / K, 2L * K)
  old_seed <- .sft_stash_seed(seed)
  on.exit(.sft_unstash_seed(old_seed), add = TRUE)
  # `a` has one column per cell in HH, HL, LH, LL order. The SIC contrast
  # weights those columns c(-1, +1, +1, -1).
  contrasts_from <- function(n, a) {
    parts <- lapply(seq_len(4L), function(cc) {
      p <- .sft_rdirichlet(n, a[, cc])
      cum <- function(idx) {
        block <- p[, idx, drop = FALSE]
        if (ncol(block) == 1L) block else t(apply(block, 1L, cumsum))
      }
      list(correct = cum(seq_len(nbin)), error = cum(K + seq_len(nbin)))
    })
    w <- c(-1, 1, 1, -1)
    add <- function(what) Reduce(`+`, Map(function(pp, ww) ww * pp[[what]], parts, w))
    correct <- add("correct"); error <- add("error")
    list(correct = correct, error = error, all = correct + error)
  }
  classify <- function(n, a) {
    cur <- contrasts_from(n, a)
    lapply(cur, function(m) {
      cls <- apply(m, 1L, .sft_sic_sign_class, tol = tolSIC)
      table(factor(cls, levels = .sft_sic_sign_levels))
    })
  }
  zero <- stats::setNames(rep(0, length(.sft_sic_sign_levels)), .sft_sic_sign_levels)
  prior <- list(correct = zero, error = zero, all = zero)
  post <- prior; N <- 0L
  prior_alpha <- matrix(alpha, nrow = 2L * K, ncol = 4L)
  post_alpha <- counts + alpha
  repeat {
    pr <- classify(nsamp, prior_alpha)
    po <- classify(nsamp, post_alpha)
    prior <- Map(function(a, b) a + as.numeric(b), prior, pr)
    post <- Map(function(a, b) a + as.numeric(b), post, po)
    N <- N + nsamp
    if (N >= maxn) break
    # Stop once the most probable class of every curve is pinned down, matching
    # the adaptive rule sicDPtest() uses for the correct-only route.
    widest <- max(vapply(post, function(p) {
      hi <- which.max(p)
      diff(stats::qbeta(c((1 - ci) / 2, 1 - (1 - ci) / 2), p[[hi]] + 1,
                        N - p[[hi]] + 1))
    }, numeric(1)))
    if (widest < 0.01) break
  }
  prior <- lapply(prior, function(p) stats::setNames(p, .sft_sic_sign_levels))
  post <- lapply(post, function(p) stats::setNames(p, .sft_sic_sign_levels))
  bf <- Map(function(p, q) (p + 1) / (q + 1), post, prior)
  list(BF = bf, BFp = lapply(bf, function(b) b / sum(b)),
       posterior = lapply(post, function(p) p / N), N = N,
       counts = list(prior = prior, posterior = post),
       ci = ci, nbin = nbin, nsamp = nsamp, seed = seed,
       tolerances = c(SIC = tolSIC))
}


# ---------------------------------------------------------------------------
# Angelov-Ekstrom permutation test for the interaction contrast
# ---------------------------------------------------------------------------

# Test statistics. Each is computed on the pooled sample of the two mixtures,
# with `alternative = "greater"` meaning "F_X lies above F_Y". Sample sizes are
# fixed under the permutation scheme, so the unnormalised forms below are
# adequate; they follow Gondan and Legner's supplement.
.sft_ae_stat_cvm <- function(x, y, ...) {
  t12 <- c(x, y)
  n1 <- length(x); n2 <- length(y)
  d <- stats::ecdf(x)(t12) - stats::ecdf(y)(t12)
  sqrt((n1 * n2) / (n1 + n2)) * (1 / (n1 + n2)) * sum(d[d > 0])
}


.sft_ae_stat_ad <- function(x, y, gm = 1, ...) {
  t12 <- c(x, y)
  n1 <- length(x); n2 <- length(y)
  H <- stats::ecdf(t12)(t12)
  w <- (H * (1 - H))^(-1 / gm)
  w[!is.finite(w)] <- 0
  if (sum(w) > 0) w <- w / sum(w)
  d <- w * (stats::ecdf(x)(t12) - stats::ecdf(y)(t12))
  sqrt(n1 * n2 / (n1 + n2)) * sum(d[d > 0])
}


.sft_ae_stat_ks <- function(x, y, ...) {
  t12 <- c(x, y)
  n1 <- length(x); n2 <- length(y)
  ord <- ifelse(order(t12) <= n1, 1 / n1, -1 / n2)
  i_inf <- sum(is.finite(t12))
  if (!i_inf) return(0)
  max(cumsum(ord[seq_len(i_inf)]))
}


# Aly et al.'s slope statistic: the largest increase of psi(s) - psi(t) over
# s <= t, used for the monotonicity predictions on the subdensity contrasts.
# cummax gives this in linear time rather than forming the full outer matrix.
.sft_ae_stat_aly <- function(x, y, ...) {
  t <- sort(c(x, y))
  t <- t[is.finite(t)]
  if (!length(t)) return(0)
  psi <- stats::ecdf(x)(t) - stats::ecdf(y)(t)
  max(cummax(psi) - psi)
}


# What each statistic's two orientations actually detect. The dominance
# statistics compare the *level* of the two mixtures, so they resolve the sign
# of the interaction contrast. Aly's slope statistic compares its *increments*,
# so it resolves monotonicity instead -- that is the form the subdensity
# predictions take (e.g. D2f- >= 0 for all t). Slot "a" is the orientation the
# `less` statistic detects and slot "b" the one `greater` detects; for Aly the
# algebra pairs an upward level excursion with a downward slope excursion,
# because psi is formed from the opposite mixture difference in each case.
.sft_ae_semantics <- function(stat) {
  if (identical(stat, "Aly")) {
    return(list(p = c(a = "p.SICdec", b = "p.SICinc"),
                a = "SIC decreasing", b = "SIC increasing",
                mixed = "SIC non-monotone",
                alternative = "the unconditional interaction contrast is not monotone"))
  }
  list(p = c(a = "p.SICpos", b = "p.SICneg"),
       a = "SIC positive", b = "SIC negative", mixed = "SIC crosses",
       alternative = "the unconditional interaction contrast is not flat")
}


.sft_ae_stat <- function(stat) {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(stat[[1L]])))
  switch(key,
         cvm = .sft_ae_stat_cvm, cramervonmises = .sft_ae_stat_cvm,
         ks = .sft_ae_stat_ks, kolmogorovsmirnov = .sft_ae_stat_ks,
         ad = .sft_ae_stat_ad, andersondarling = .sft_ae_stat_ad,
         aly = .sft_ae_stat_aly, slope = .sft_ae_stat_aly,
         stop("stat must be one of 'CvM', 'KS', 'AD', or 'Aly'.", call. = FALSE))
}


# Reduce one participant's four cells to the two mixtures. X pools the two
# cells with matching salience (AB and ab, i.e. HH and LL), Y pools the mixed
# cells (Ab and aB, i.e. HL and LH), so that
#   F_X - F_Y = (F_HH + F_LL - F_HL - F_LH) / 2 = D2F / 2 = -SIC / 2.
# The half-and-half mixture assumes equal trial counts per cell; `balance`
# subsamples to the participant's smallest cell when they differ.
.sft_ae_balance <- function(cells, balance = TRUE) {
  ns <- vapply(cells, function(z) z$n, integer(1))
  if (!balance || length(unique(ns)) == 1L) return(lapply(ns, seq_len))
  m <- min(ns)
  lapply(cells, function(z) sort(sample.int(z$n, m)))
}


.sft_ae_mixture <- function(cells, idx, contrast) {
  pick <- switch(contrast,
                 correct = function(z) z$correct,
                 error = function(z) !z$correct,
                 all = function(z) rep(TRUE, z$n))
  vals <- Map(function(z, i) {
    rt <- z$rt[i]; keep <- pick(z)[i]
    rt[!is.finite(rt)] <- Inf
    rt[!keep] <- Inf
    rt
  }, cells, idx)
  list(X = c(vals$HH, vals$LL), Y = c(vals$HL, vals$LH))
}


#' Permutation test of the unconditional interaction contrast
#'
#' Applies Angelov and Ekstrom's (2023) repeated-measures test for stochastic
#' dominance to the interaction contrast of the unconditional response-time
#' subdistributions, following Gondan and Legner. Each participant contributes
#' two mixtures, \eqn{F_X = (F_{HH} + F_{LL})/2} and
#' \eqn{F_Y = (F_{HL} + F_{LH})/2}, whose difference is half the interaction
#' contrast. The null distribution is obtained by randomly swapping the two
#' mixture labels within each participant.
#'
#' @details
#' The test is run separately for the correct subdistribution, the error
#' subdistribution, and the accuracy-agnostic distribution. Responses of the
#' other type, and omissions, are set to \code{+Inf} so that they stay in the
#' denominator without contributing at any finite time.
#'
#' The Cramer-von Mises, Kolmogorov-Smirnov and Anderson-Darling statistics
#' compare the level of the two mixtures and so resolve the *sign* of the
#' interaction contrast, reporting \code{p.SICpos} and \code{p.SICneg} and a
#' decision among \code{"n.s."}, \code{"SIC positive"}, \code{"SIC negative"}
#' and \code{"SIC crosses"}. Aly's slope statistic compares their increments
#' instead and so resolves *monotonicity*, which is the form the subdensity
#' predictions take (for instance a non-negative error subdensity contrast at
#' every t). It reports \code{p.SICinc} and \code{p.SICdec} and decides among
#' \code{"n.s."}, \code{"SIC increasing"}, \code{"SIC decreasing"} and
#' \code{"SIC non-monotone"}.
#'
#' The four-way decision follows Angelov and Ekstrom's published rule, which is
#' deliberately asymmetric: \code{"SIC positive"} requires
#' \code{p.SICpos <= alpha} together with \code{p.SICneg > alpha}, whereas
#' \code{"SIC negative"} requires \code{p.SICneg <= alpha} together with
#' \code{p.SICpos > alpha.star}. A contrast with strong one-sided evidence can
#' therefore still be labelled \code{"SIC crosses"} when the opposing p-value
#' falls between \code{alpha} and \code{alpha.star}. Both one-sided p-values
#' are returned so a different rule can be applied.
#'
#' Signs are reported in this package's survivor convention, \code{SIC = -D2F};
#' \code{"SIC positive"} therefore corresponds to Gondan and Legner's
#' "IC negative".
#'
#' This is a group-level test: with a single participant the permutation
#' distribution has only two points and the p-values are uninformative.
#'
#' @param inData Canonical SFT trial data with \code{Subject}, \code{RT},
#'   \code{Correct}, \code{Channel1}, and \code{Channel2} columns. A
#'   \code{Condition} column is used to subset when \code{Condition} is given.
#' @param contrast Which subdistributions to test: any of \code{"correct"},
#'   \code{"error"}, and \code{"agnostic"}.
#' @param stat Test statistic: \code{"CvM"} (default), \code{"KS"},
#'   \code{"AD"}, or \code{"Aly"}.
#' @param nperm Number of permutations.
#' @param alpha Significance level.
#' @param alpha.star Tuning parameter separating crossing from dominating
#'   distribution functions.
#' @param gm Weight exponent for the Anderson-Darling statistic.
#' @param balance Subsample each participant's cells to a common size before
#'   forming the mixtures. Leave \code{TRUE} unless the design is already
#'   balanced.
#' @param min.trials Minimum trials per cell for a participant to be included.
#' @param Condition Optional condition to subset before testing.
#' @param seed Optional random seed.
#' @return A list with one element per requested contrast, each an \code{htest}
#'   carrying both one-sided p-values and the decision, plus a \code{summary}
#'   data frame and the participants used.
#' @references
#' Angelov, A. G., & Ekstrom, M. (2023). Tests of stochastic dominance with
#' repeated measurements data. \emph{AStA Advances in Statistical Analysis},
#' 107, 443-467.
#' @seealso [sic()], [sicGroup()], [siDominance()]
#' @export
sicPermTest <- function(inData, contrast = c("correct", "error", "agnostic"),
                        stat = c("CvM", "KS", "AD", "Aly"), nperm = 1001L,
                        alpha = .05, alpha.star = .8, gm = 1, balance = TRUE,
                        min.trials = 10L, Condition = NULL, seed = NULL) {
  contrast <- match.arg(contrast, several.ok = TRUE)
  stat <- match.arg(stat)
  stat_fn <- .sft_ae_stat(stat)
  sem <- .sft_ae_semantics(stat)
  nperm <- max(1L, as.integer(nperm))
  subjects <- .sft_sic_subject_cells(inData, errors = "defective",
                                     Condition = Condition,
                                     min.trials = min.trials)
  if (length(subjects) < 2L) {
    stop("sicPermTest() needs at least two participants with usable cells.",
         call. = FALSE)
  }
  # Pooling participants across conditions would test a contrast no architecture
  # predicts, so make the caller choose rather than doing it silently.
  conds <- unique(vapply(subjects, function(s) s$condition, character(1)))
  if (length(conds) > 1L) {
    stop("inData spans several conditions (", paste(conds, collapse = ", "),
         "); pass Condition to select one.", call. = FALSE)
  }
  old_seed <- .sft_stash_seed(seed)
  on.exit(.sft_unstash_seed(old_seed), add = TRUE)
  key <- c(correct = "correct", error = "error", agnostic = "all")
  # Balance once, so all three contrasts are tested on the same trials and
  # SIC = SIC.correct + SIC.error still holds for what the test sees.
  idx <- lapply(subjects, function(s) .sft_ae_balance(s$cells, balance))
  out <- list()
  for (cn in contrast) {
    mixes <- Map(function(s, i) .sft_ae_mixture(s$cells, i, key[[cn]]), subjects, idx)
    observed <- .sft_ae_pair(mixes, stat_fn, swap = rep(FALSE, length(mixes)), gm = gm)
    null <- replicate(nperm, {
      swap <- sample(c(TRUE, FALSE), length(mixes), replace = TRUE)
      .sft_ae_pair(mixes, stat_fn, swap = swap, gm = gm)
    })
    # For the dominance statistics `greater` detects F_X above F_Y (D2F > 0,
    # SIC < 0) and `less` the reverse; for Aly they detect the corresponding
    # slope excursions. See .sft_ae_semantics().
    p_b <- (1 + sum(null["greater", ] >= observed[["greater"]])) / (1 + nperm)
    p_a <- (1 + sum(null["less", ] >= observed[["less"]])) / (1 + nperm)
    # Angelov and Ekstrom's published decision rule, kept verbatim including its
    # asymmetry: the `a` verdict needs only alpha on the opposing side, the `b`
    # verdict needs alpha.star.
    decision <- if (p_a > alpha && p_b > alpha) "n.s." else
      if (p_a <= alpha && p_b > alpha) sem$a else
        if (p_a > alpha.star) sem$b else sem$mixed
    st <- c(observed[["greater"]], observed[["less"]])
    names(st) <- c("T(greater)", "T(less)")
    # The headline p.value combines the two one-sided tests with a Bonferroni
    # correction, so it remains a valid p-value for the exchangeability null.
    # The four-way `decision`, not this number, is what the theory is read off.
    h <- .sft_htest(st, min(1, 2 * min(p_a, p_b)), sem$alternative,
                    paste0("Angelov-Ekstrom permutation test of the ", cn,
                           " interaction contrast (", stat, ")"),
                    paste0(length(subjects), " participants"))
    h[[sem$p[["a"]]]] <- p_a; h[[sem$p[["b"]]]] <- p_b
    h$decision <- decision
    h$alpha <- alpha; h$alpha.star <- alpha.star; h$nperm <- nperm
    out[[cn]] <- h
  }
  out$summary <- data.frame(
    Contrast = contrast,
    a = vapply(contrast, function(cn) out[[cn]][[sem$p[["a"]]]], numeric(1)),
    b = vapply(contrast, function(cn) out[[cn]][[sem$p[["b"]]]], numeric(1)),
    Decision = vapply(contrast, function(cn) out[[cn]]$decision, character(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  names(out$summary)[2:3] <- unname(sem$p[c("a", "b")])
  out$subjects <- vapply(subjects, function(s) s$subject, character(1))
  out$stat <- stat
  out$nperm <- nperm
  out
}


# One permutation replicate: `swap` flips the X/Y labels for whole participants,
# which is exactly the exchangeability the repeated-measures test assumes.
.sft_ae_pair <- function(mixes, stat_fn, swap, gm) {
  xs <- unlist(Map(function(m, s) if (s) m$Y else m$X, mixes, swap), use.names = FALSE)
  ys <- unlist(Map(function(m, s) if (s) m$X else m$Y, mixes, swap), use.names = FALSE)
  c(greater = stat_fn(xs, ys, gm = gm), less = stat_fn(ys, xs, gm = gm))
}


# ---------------------------------------------------------------------------
# Data-frame plumbing shared by sicGroup() and sicPermTest()
# ---------------------------------------------------------------------------

# Split canonical trial data into per-subject (and optionally per-condition)
# SIC cells. Under "discard" only correct trials are kept, reproducing the
# historical sicGroup() behaviour; under "defective" every trial is retained
# with its correctness flag.
.sft_sic_subject_cells <- function(inData, errors = "discard", Condition = NULL,
                                   min.trials = 10L) {
  inData <- .sft_normalize_columns(inData)
  req <- c("Subject", "RT", "Correct", "Channel1", "Channel2")
  if (!all(req %in% names(inData))) {
    stop("inData is missing required SIC columns: ",
         paste(setdiff(req, names(inData)), collapse = ", "), call. = FALSE)
  }
  if (!"Condition" %in% names(inData)) inData$Condition <- "all"
  if (!is.null(Condition)) {
    inData <- inData[inData$Condition %in% Condition, , drop = FALSE]
  }
  correct <- as.logical(inData$Correct)
  correct[is.na(correct)] <- FALSE
  inData$.correct <- correct
  out <- list(); k <- 0L
  for (cond in unique(inData$Condition)) for (subj in unique(inData$Subject)) {
    d <- inData[inData$Condition == cond & inData$Subject == subj, , drop = FALSE]
    if (!nrow(d)) next
    grab <- function(c1, c2) d[d$Channel1 == c1 & d$Channel2 == c2, , drop = FALSE]
    parts <- list(HH = grab(2, 2), HL = grab(2, 1), LH = grab(1, 2), LL = grab(1, 1))
    usable <- vapply(parts, function(p) {
      if (errors == "discard") sum(p$.correct & is.finite(p$RT)) else nrow(p)
    }, numeric(1))
    if (min(usable) <= min.trials) next
    cells <- .sft_sic_cells(parts$HH$RT, parts$HL$RT, parts$LH$RT, parts$LL$RT,
                            CR = list(HH = parts$HH$.correct, HL = parts$HL$.correct,
                                      LH = parts$LH$.correct, LL = parts$LL$.correct),
                            errors = errors)
    k <- k + 1L
    out[[k]] <- list(subject = as.character(subj), condition = as.character(cond),
                     cells = cells)
  }
  out
}


# ---------------------------------------------------------------------------
# Plot builder
# ---------------------------------------------------------------------------

#' Tidy data frame of unconditional SIC curves
#'
#' Long-format companion to [build_sic_df()] for the three interaction
#' contrasts produced by `errors = "defective"`: correct, error, and
#' accuracy-agnostic.
#'
#' @param times Times at which to evaluate the curves.
#' @param series_specs List of series. Each element needs a `label` and an `fn`
#'   list with `correct` and `error` sublists holding the four cell
#'   subdistribution functions named `HH`, `HL`, `LH`, and `LL`.
#' @param trim_cdf_tails Drop times where every cell distribution is in its
#'   extreme tail.
#' @param cdf_tail_cut Tail cut used when trimming.
#' @param smooth_cdf Smoothing applied to each cell distribution before
#'   contrasting: `"none"`, `"mono"`, or `"spline"`.
#' @param smooth_spar Smoothing parameter for `smooth_cdf = "spline"`.
#' @return A data frame with `Time`, `SIC`, `Contrast`, and `Series` columns.
#' @seealso [build_sic_df()], [sicGroup()]
#' @export
build_sic_defective_df <- function(times, series_specs, trim_cdf_tails = TRUE,
                                   cdf_tail_cut = 5e-4,
                                   smooth_cdf = c("none", "mono", "spline"),
                                   smooth_spar = .65) {
  smooth_cdf <- match.arg(smooth_cdf)
  cells <- c("HH", "HL", "LH", "LL")
  .sft_bind_rows(lapply(series_specs, function(spec) {
    ev <- function(what) {
      y <- lapply(cells, function(nm) spec$fn[[what]][[nm]](times))
      names(y) <- cells
      if (smooth_cdf == "none") return(y)
      y <- lapply(y, function(v) smooth_one_cdf(times, v, smooth_cdf, smooth_spar))
      names(y) <- cells
      y
    }
    fc <- ev("correct"); fe <- ev("error")
    fa <- Map(`+`, fc, fe)
    sic <- function(y) y$LH + y$HL - y$HH - y$LL
    out <- .sft_bind_rows(list(
      data.frame(Time = times, SIC = sic(fc), Contrast = "Correct",
                 Series = spec$label, stringsAsFactors = FALSE),
      data.frame(Time = times, SIC = sic(fe), Contrast = "Error",
                 Series = spec$label, stringsAsFactors = FALSE),
      data.frame(Time = times, SIC = sic(fa), Contrast = "Agnostic",
                 Series = spec$label, stringsAsFactors = FALSE)))
    if (isTRUE(trim_cdf_tails)) {
      cut <- if (is.finite(cdf_tail_cut) && cdf_tail_cut > 0 && cdf_tail_cut < .5)
        cdf_tail_cut else 5e-4
      # Trim on the accuracy-agnostic distributions; the subdistributions never
      # reach 1, so the usual upper cut cannot be applied to them directly.
      lo <- Reduce(`&`, lapply(fa, function(v) v < cut))
      top <- max(vapply(fa, function(v) max(v), numeric(1)))
      hi <- Reduce(`&`, lapply(fa, function(v) v > top - cut))
      drop <- lo | hi
      out <- out[!rep(drop, times = 3L), , drop = FALSE]
    }
    out
  }))
}
