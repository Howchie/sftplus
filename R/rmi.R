# Race model inequality (RMI) scoring: the Miller upper bound and the Grice
# lower bound, evaluated on a per-participant percentile grid.
#
# Miller's (1982) inequality states that for a redundant-target condition AB and
# single-target conditions A and B,
#
#   F_AB(t) <= F_A(t) + F_B(t)   for all t,
#
# under any separate-activation (race) account.  Two assumptions are worth
# separating, because they are routinely confused (Miller, 2016).  Stochastic
# independence of the two channels' finishing times is *not* required: the
# channels may be correlated in an arbitrary way, which is exactly why only a
# bound, and not a point prediction, is available.  Context independence -- that
# each channel's finishing-time distribution is the same on redundant trials as
# on single-target trials -- *is* required, and Miller (2016) argues it is not a
# separate side condition but an inherent part of Raab's (1962) statistical
# facilitation account.  A model in which the channels run at different speeds
# when the other signal is present is therefore not a race model that happens to
# violate an auxiliary assumption; it is a species of coactivation, since both
# signals influence the response on a single trial.  This is what licenses the
# usual reading of a violation: RMI violations justify concluding coactivation.
# It also means the design must minimise the opportunity for context dependence
# -- mixed rather than blocked presentation of the stimulus conditions.
#
# The Grice bound F_AB(t) >= max(F_A(t), F_B(t)) is the corresponding lower
# limit.  Both are the two-channel OR case of the bounds in
# R/estimate.bounds.R, but this file estimates them from its own cumulative
# frequency polygons (see .sft_rmi_cdf below) rather than calling
# estimate.bounds(), so that the bound and the redundant-target quantile it is
# compared against come from one estimator with one tie rule.
#
# This is a different hypothesis from the UCIP capacity coefficient.  C_OR(t) > 1
# is violation of one *specific* model -- independent, unlimited-capacity
# parallel processing -- whereas Miller's bound holds for every context-
# independent race model regardless of stochastic dependence.  Super capacity
# therefore does not imply a Miller violation, and a race whose channels are
# slowed by load violates UCIP while satisfying Miller -- though note that under
# Miller (2016) such slowing is itself context dependence, so the two hypotheses
# partition the model space differently rather than nesting.  Nothing in this
# file reuses the UCIP score.
#
# The estimand is on the probability scale.  For participant i and percentile
# level u_k,
#
#   t_ik     = the u_k-th quantile of participant i's own pooled correct
#              single-target RTs
#   Miller:  delta_ik = F_A(t_ik) + F_B(t_ik) - F_AB(t_ik)
#   Grice:   delta_ik = F_AB(t_ik) - max(F_A(t_ik), F_B(t_ik))
#
# so that a *negative* delta is a violation under either bound and the shared
# sign / ROPE machinery reads identically for both.  Two properties make this
# free of participant base time.  First, the inequality is evaluated within a
# participant, so overall speed cancels exactly and no millisecond quantity
# enters the group statistic.  Second, the inequality is invariant to any
# strictly increasing reparameterisation of time applied to all three conditions
# of a participant, which is what licenses evaluating each participant on its own
# percentile grid and then pooling on the shared u index.  Miller's original
# procedure instead averaged millisecond quantiles across participants, so a slow
# participant's 40 ms gap outweighed a fast participant's 15 ms gap even when
# both violated the bound equally.


.sft_rmi_bound <- function(bound) {
  key <- tolower(gsub("[^a-z]", "", as.character(bound[[1L]])))
  switch(key,
         miller = "miller", upper = "miller", or = "miller", rmi = "miller",
         grice = "grice", lower = "grice",
         both = "both",
         stop("bound must be one of 'miller', 'grice', or 'both'."))
}


.sft_rmi_probs <- function(probs) {
  probs <- sort(unique(as.numeric(probs)))
  if (!length(probs) || any(!is.finite(probs)) || any(probs <= 0) ||
      any(probs >= 1)) {
    stop("probs must be finite percentile levels strictly inside (0, 1).",
         call. = FALSE)
  }
  if (length(probs) < 3L) {
    stop("probs must contain at least three percentile levels; the shape of the ",
         "violation curve and the min-over-grid global summary are not ",
         "meaningful below three.", call. = FALSE)
  }
  probs
}


# Label a percentile level for use as a population-mean name.
.sft_rmi_prob_label <- function(p) {
  paste0("p", format(round(p, 4), trim = TRUE, scientific = FALSE))
}


# ---- Per-block quantile estimation -----------------------------------------
#
# Miller (1982) estimated the quantiles separately within each participant x
# block and averaged over blocks.  That matters when overall speed drifts across
# the session: pooling a fast late block with a slow early block spreads the
# pooled distribution by the between-block shift, which flattens the CDFs of all
# three cells and moves the bound.  Averaging within-block quantiles instead
# estimates the block-average distribution, with the drift absorbed block by
# block.  It is the same argument that makes the estimand within-participant,
# applied one level down.
#
# Every scorer below therefore takes a *list of blocks* for one participant --
# each block a list(RT, CR) holding the AB / A / B cells -- and reduces over
# blocks by averaging.  With `by_block = FALSE` a participant has exactly one
# block covering the whole session, the reduction is the identity, and the
# default path is numerically unchanged.
#
# Because blocks are disjoint sets of trials they are independent, so the
# sampling covariance of the block-averaged delta is the sum of the per-block
# covariances divided by the squared number of blocks.  `min_n` applies per
# (participant, block, cell): a block short of it is dropped, and a participant
# with no surviving block is dropped as before.


.sft_rmi_block_levels <- function(x) {
  if (is.factor(x)) return(as.character(levels(droplevels(x))))
  as.character(sort(unique(x)))
}


# Condition-split input as .sft_condition_input() returns it, with an added
# `blocks` element: one list of per-block RT/CR cell lists per participant.
.sft_rmi_input <- function(inData, CR = NULL, Condition = NULL, Subject = NULL,
                           by_block = FALSE) {
  grouped <- .sft_condition_input(inData, CR, "OR", Condition, Subject)
  grouped$by_block <- isTRUE(by_block)
  if (!isTRUE(by_block)) return(grouped)
  if (!is.data.frame(inData)) {
    stop("by_block = TRUE needs a trial-level data frame; the pre-split RT/CR ",
         "list input carries no block index.", call. = FALSE)
  }
  data <- .sft_normalize_columns(inData)
  if (!"Block" %in% names(data)) {
    stop("by_block = TRUE needs a Block column in inData.", call. = FALSE)
  }
  if (anyNA(data$Block)) {
    stop("Block must not contain missing values.", call. = FALSE)
  }
  block <- as.character(data$Block)
  block_levels <- .sft_rmi_block_levels(data$Block)
  has_condition <- "Condition" %in% names(data)
  condition <- if (has_condition) as.character(data$Condition) else
    rep("default", nrow(data))
  subject <- as.character(data$Subject)

  conditions <- lapply(grouped$levels, function(level) {
    input <- grouped$conditions[[level]]
    blocks <- stats::setNames(vector("list", length(input$subject)),
                              input$subject)
    labels <- blocks
    for (b in block_levels) {
      rows <- block == b & condition == level
      if (!any(rows)) next
      present <- intersect(input$subject, unique(subject[rows]))
      if (!length(present)) next
      # A block with no usable channel cell for the selected participants is
      # simply not a block of this analysis.
      converted <- tryCatch(
        sft_data_to_rt(data[rows, , drop = FALSE],
                       Condition = if (has_condition) level else NULL,
                       Subject = present, stopping.rule = "OR",
                       by_subject = "always"),
        error = function(e) NULL)
      if (is.null(converted)) next
      for (id in present) {
        blocks[[id]] <- c(blocks[[id]], list(list(RT = converted$RT[[id]],
                                                  CR = converted$CR[[id]])))
        labels[[id]] <- c(labels[[id]], b)
      }
    }
    input$blocks <- unname(blocks)
    input$block_names <- unname(labels)
    input
  })
  names(conditions) <- grouped$levels
  grouped$conditions <- conditions
  grouped$block_levels <- block_levels
  grouped
}


# One participant's blocks.  Inputs that never went through .sft_rmi_input()
# (the list path, and direct internal callers) are treated as one block.
.sft_rmi_subject_blocks <- function(input, i) {
  if (!is.null(input$blocks)) return(input$blocks[[i]])
  list(list(RT = input$RT[[i]], CR = input$CR[[i]]))
}


# Average a per-block field across blocks, keeping the percentile index.  Levels
# no block could estimate stay NA rather than becoming NaN.
.sft_rmi_block_mean <- function(per, field, k, na.rm = FALSE) {
  m <- matrix(unlist(lapply(per, `[[`, field)), nrow = k)
  out <- rowMeans(m, na.rm = na.rm)
  out[!is.finite(out)] <- NA_real_
  out
}


# Response times per cell, sorted.
#
# errors = "discard" keeps correct trials only.  This matches the treatment in
# estimate.bounds() (R/estimate.bounds.R:45) and in Miller's own analyses:
# incorrect and omitted trials are dropped rather than censored.  The
# semiparametric hazard model in this package treats errors as censored exposure
# instead, but the race inequality is a statement about the observed correct-RT
# distributions, so discarding keeps the estimand the one Miller defined.
#
# errors = "defective" instead keeps every trial and records an error as Inf, so
# it is counted in the denominator but never satisfies RT <= t.  Each cell CDF is
# then the defective distribution F(t) = P(correct and RT <= t), whose asymptote
# is the accuracy rather than 1.  This is the right treatment when accuracy
# differs across cells, because under "discard" a redundant-target condition that
# is both faster and more accurate has its extra correct responses renormalised
# away.  Representing an error as Inf rather than as a separate denominator means
# every downstream estimator -- the ECDF, its covariance, and the bootstrap --
# needs no special case.
.sft_rmi_cells <- function(rt, cr = NULL, errors = c("discard", "defective")) {
  errors <- match.arg(errors)
  required <- c("AB", "A", "B")
  missing <- setdiff(required, names(rt))
  if (length(missing)) {
    stop("Race-model scoring needs the AB, A, and B cells; missing: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  out <- lapply(required, function(cell) {
    cleaned <- .sft_clean(rt[[cell]], if (is.null(cr)) NULL else cr[[cell]])
    if (errors == "discard") return(sort(cleaned$RT[cleaned$CR]))
    sort(ifelse(cleaned$CR, cleaned$RT, Inf))
  })
  names(out) <- required
  out
}


# Cells for each of one participant's blocks, keeping the blocks that carry
# enough correct trials in every cell and admit an estimable percentile grid.
# min_n counts observed responses: under errors = "defective" the Inf
# placeholders inflate the cell size without adding information about shape.
.sft_rmi_usable_blocks <- function(blocks, probs, min_n, errors = "discard",
                                   grid = TRUE) {
  cells <- lapply(blocks, function(b) {
    .sft_rmi_cells(b$RT, b$CR, errors = errors)
  })
  observed <- lapply(cells, function(cl) {
    vapply(cl, function(v) sum(is.finite(v)), integer(1))
  })
  ok <- vapply(seq_along(cells), function(j) {
    all(observed[[j]] >= min_n) &&
      (!grid || !anyNA(.sft_rmi_grid(cells[[j]], probs)))
  }, logical(1))
  # The reported cell counts are the trials that entered the estimate, so they
  # cover the surviving blocks; with none surviving they describe what was
  # available, which is what the drop message is about.
  n <- Reduce(`+`, lapply(cells[if (any(ok)) ok else TRUE], lengths))
  # Quantiles are estimated within block, so the estimable percentile range is
  # governed by the smallest cell of a surviving block, not by the pooled count.
  n_min <- if (any(ok)) min(unlist(observed[ok])) else NA_integer_
  # The probability-scale grid is read off the *pooled* single-target RTs, so
  # its estimable range is governed by n_A + n_B rather than by the worst cell.
  n_pool <- if (any(ok)) {
    min(vapply(observed[ok], function(v) v[["A"]] + v[["B"]], integer(1)))
  } else NA_integer_
  list(cells = cells[ok], n = n, n_blocks = sum(ok), n_min = n_min,
       n_pool = n_pool)
}


.sft_rmi_ecdf_at <- function(x, times) {
  if (!length(x)) return(rep(NA_real_, length(times)))
  colMeans(outer(x, times, "<="))
}


# ---- The Ulrich-Miller-Schroeter cumulative frequency polygon ---------------
#
# By default, everything on the *time* scale -- the reported cell quantiles, the
# bound curve, and the bound quantiles compared against them -- is estimated
# with the algorithm of Ulrich, Miller and Schroeter (2007), so that rmi.test()
# reproduces the published procedure exactly. The optional kernel path applies
# one Gaussian kernel-CDF estimator symmetrically to all of those quantities.
# (The probability-scale estimand used by rmiGroup.bayes() is a different
# quantity and uses the plain empirical CDF above; see the header.)
#
# Their Equation 2 is a cumulative frequency polygon rather than a step
# function: the ith of n ordered RTs is placed at (i - .5)/n and adjacent points
# are joined by a straight line.  UMS are explicit that the step function is not
# an acceptable substitute for percentile estimation.  `qtype` generalises the
# plotting-position rule to the others offered by stats::quantile(); qtype = 5
# is Equation 2 itself.
#
# Their Appendix A gives the correction for ties: RT is recorded to the nearest millisecond, so ties "could occur in
# any real data set".  A value repeated n_i times is a single vertical jump of
# the step function, from s_{i-1}/n to s_i/n, and the polygon must pass through
# that jump's *midpoint*, (s_{i-1} + s_i)/(2n) -- not its top, and not a flat
# plateau spanning the tied probabilities.  The midpoint is the arithmetic mean
# of the plotting positions of the tied observations, so averaging within a tied
# group implements Appendix A for any rule and collapses to the untied case when
# n_i = 1.  Both of the estimators previously used here got this wrong in
# opposite directions -- stats::quantile() produces the plateau and the
# quantile-rule CDF in estimate.bounds() took the top of the jump -- which
# biased the redundant-target quantile and the bound quantile against each
# other by up to a millisecond on the worked example of UMS Table 1.
.sft_rmi_positions <- function(x, qtype = 5) {
  x <- sort(as.numeric(x[is.finite(x)]))
  n <- length(x)
  if (!n) return(NULL)
  i <- seq_len(n)
  p <- switch(as.character(qtype),
              "1" = i / n,
              "2" = i / n,
              "3" = (i - 0.5) / n,
              "4" = i / n,
              "5" = (i - 0.5) / n,
              "6" = i / (n + 1),
              "7" = if (n > 1L) (i - 1) / (n - 1) else 1,
              "8" = (i - 1 / 3) / (n + 1 / 3),
              "9" = (i - 3 / 8) / (n + 0.25),
              stop("qtype must be an integer between 1 and 9.", call. = FALSE))
  ux <- unique(x)
  # x is sorted, so tied observations are contiguous and the group means are a
  # difference of cumulative sums.
  counts <- tabulate(match(x, ux), length(ux))
  upper <- cumsum(counts)
  lower <- c(0L, upper[-length(upper)])
  cp <- c(0, cumsum(p))
  list(x = ux, p = (cp[upper + 1L] - cp[lower + 1L]) / counts, n = n)
}


# The polygon itself.  Below the smallest observation UMS set G to zero, so the
# polygon carries a discontinuity there; that is what makes probabilities below
# the first plotting position inestimable rather than silently extrapolated.
.sft_rmi_cdf <- function(x, qtype = 5) {
  pos <- .sft_rmi_positions(x, qtype)
  if (is.null(pos)) return(function(t) rep(NA_real_, length(t)))
  function(t) {
    out <- numeric(length(t))
    out[t >= pos$x[[length(pos$x)]]] <- 1
    mid <- t >= pos$x[[1L]] & t < pos$x[[length(pos$x)]]
    if (any(mid)) {
      out[mid] <- if (length(pos$x) == 1L) pos$p else
        stats::approx(pos$x, pos$p, xout = t[mid])$y
    }
    out
  }
}


# Percentiles by inverting the polygon.  UMS note that there is no satisfactory
# estimate of a percentile below p = 1/(2n) or above p = (2n - 1)/(2n), since
# either would require extrapolating past the observed data; for qtype = 5 those
# limits are precisely the first and last plotting positions, so returning NA
# outside the polygon's range enforces the rule exactly rather than clamping to
# the sample minimum or maximum the way stats::quantile() does.
.sft_rmi_quantile <- function(x, probs, qtype = 5) {
  pos <- .sft_rmi_positions(x, qtype)
  if (is.null(pos) || length(pos$x) < 2L) return(rep(NA_real_, length(probs)))
  stats::approx(pos$p, pos$x, xout = probs, ties = "ordered")$y
}


# Invert an exact kernel CDF (or an algebraic bound made from kernel CDFs).
.sft_rmi_kernel_invert <- function(cdf, probs, lower, upper) {
  vapply(probs, function(p) {
    at_ends <- c(cdf(lower), cdf(upper))
    if (any(!is.finite(at_ends)) || at_ends[1L] > p || at_ends[2L] < p) {
      return(NA_real_)
    }
    stats::uniroot(function(t) cdf(t) - p, lower = lower, upper = upper,
                   tol = .Machine$double.eps^0.5)$root
  }, numeric(1))
}


# Symmetric in the statistical-procedure sense: every observed condition and
# every channel entering the bound uses the same Gaussian kernel-CDF estimator
# and bandwidth-selection rule.
.sft_rmi_kernel_estimates <- function(cells, probs, bound,
                                      kernel_bw = "nrd0") {
  cdfs <- lapply(cells[c("AB", "A", "B")], .make_kernel_cdf, bw = kernel_bw)
  lower <- min(vapply(cdfs, function(f) exp(log(attr(f, "range")[1L]) - 8 * attr(f, "bw")), numeric(1)))
  upper <- max(vapply(cdfs, function(f) exp(log(attr(f, "range")[2L]) + 8 * attr(f, "bw")), numeric(1)))
  bound_cdf <- if (bound == "miller") {
    function(t) cdfs$A(t) + cdfs$B(t)
  } else {
    function(t) pmax(cdfs$A(t), cdfs$B(t))
  }
  list(q_AB = .sft_rmi_kernel_invert(cdfs$AB, probs, lower, upper),
       q_A = .sft_rmi_kernel_invert(cdfs$A, probs, lower, upper),
       q_B = .sft_rmi_kernel_invert(cdfs$B, probs, lower, upper),
       q_bound = .sft_rmi_kernel_invert(bound_cdf, probs, lower, upper))
}


# The bound curve B(t) = G_A(t) + G_B(t) (Miller) or max(G_A, G_B) (Grice).
# For polygons, B is piecewise linear between their knots, so the knot union is
# exact and grid-free; points just below the cell minima preserve their jumps
# from zero. For kernel CDFs this helper returns a dense display grid. The RMI
# quantiles themselves use direct root finding in .sft_rmi_kernel_estimates().
.sft_rmi_bound_curve <- function(cells, bound, qtype = 5,
                                 cdf_method = c("polygon", "kernel"),
                                 kernel_bw = "nrd0") {
  cdf_method <- match.arg(cdf_method)
  a <- cells$A[is.finite(cells$A)]
  b <- cells$B[is.finite(cells$B)]
  if (length(a) < 2L || length(b) < 2L) return(NULL)
  if (cdf_method == "kernel") {
    cdf_a <- .make_kernel_cdf(a, bw = kernel_bw)
    cdf_b <- .make_kernel_cdf(b, bw = kernel_bw)
    lower <- min(exp(log(min(a)) - 8 * attr(cdf_a, "bw")),
                 exp(log(min(b)) - 8 * attr(cdf_b, "bw")))
    upper <- max(exp(log(max(a)) + 8 * attr(cdf_a, "bw")),
                 exp(log(max(b)) + 8 * attr(cdf_b, "bw")))
    grid <- seq(lower, upper, length.out = 4096L)
    g_a <- cdf_a(grid)
    g_b <- cdf_b(grid)
  } else {
    knots <- sort(unique(c(a, b)))
    eps <- max(diff(range(knots)) * 1e-9, 1e-9)
    grid <- sort(unique(c(knots, min(a) - eps, min(b) - eps)))
    g_a <- .sft_rmi_cdf(a, qtype)(grid)
    g_b <- .sft_rmi_cdf(b, qtype)(grid)
  }
  list(grid = grid, value = if (bound == "miller") g_a + g_b else pmax(g_a, g_b))
}


# Percentile levels at which at least one cell of at least one participant falls
# outside the UMS estimable range.  Reported once by the caller rather than per
# participant, which would be unreadable.
.sft_rmi_inestimable <- function(probs, n_min) {
  n_min <- n_min[is.finite(n_min) & n_min > 0]
  if (!length(n_min)) return(numeric(0))
  n <- min(n_min)
  probs[probs < 1 / (2 * n) | probs > (2 * n - 1) / (2 * n)]
}


.sft_rmi_warn_inestimable <- function(probs, n_min) {
  bad <- .sft_rmi_inestimable(probs, n_min)
  if (!length(bad)) return(invisible(NULL))
  warning("Percentile level(s) ", paste(format(bad), collapse = ", "),
          " lie outside the estimable range 1/(2n) to (2n-1)/(2n) for at ",
          "least one participant and cell (smallest cell n = ",
          min(n_min[is.finite(n_min) & n_min > 0]), "); Ulrich, Miller & ",
          "Schroeter (2007) note there is no satisfactory estimate there. ",
          "They are reported as NA.", call. = FALSE)
  invisible(NULL)
}


# The evaluation grid: quantiles of the participant's own pooled correct
# single-target RTs.  Pooling A and B (rather than using AB) keeps the grid
# independent of the condition whose CDF is being tested, and matches the
# distributions Miller inverted.  type = 8 follows the convention used elsewhere
# in the package (R/semiparametric_prep.R:263).  Under errors = "defective" the
# Inf placeholders are dropped here: the grid is a set of evaluation times, so it
# is read off the observed responses whichever denominator the CDFs use.
.sft_rmi_grid <- function(cells, probs) {
  pooled <- c(cells$A, cells$B)
  pooled <- pooled[is.finite(pooled)]
  if (length(pooled) < 2L) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(pooled, probs = probs, names = FALSE, type = 8))
}


# Bound values and the violation contrast at supplied times.  Reproduces the
# two-channel OR case of estimate.bounds(): upper = F_A + F_B (Miller), lower =
# max(F_A, F_B) (Grice).
.sft_rmi_delta_at <- function(cells, times, bound) {
  f_ab <- .sft_rmi_ecdf_at(cells$AB, times)
  f_a <- .sft_rmi_ecdf_at(cells$A, times)
  f_b <- .sft_rmi_ecdf_at(cells$B, times)
  bound_value <- if (bound == "miller") f_a + f_b else pmax(f_a, f_b)
  delta <- if (bound == "miller") bound_value - f_ab else f_ab - bound_value
  list(delta = delta, bound = bound_value, F_AB = f_ab, F_A = f_a, F_B = f_b)
}


# Analytic sampling covariance of the delta vector.
#
# Each cell CDF is an empirical proportion, so for one condition with n trials
# Cov(Fhat(s), Fhat(t)) = F(min(s,t)) (1 - F(max(s,t))) / n, and monotonicity of
# F lets that be written with pmin/pmax on the fitted values directly.  The three
# conditions come from disjoint trials and are independent, so their covariances
# add.
#
# For the Grice bound max(F_A, F_B) is not differentiable where the two channel
# CDFs cross; away from a crossing the covariance below is the delta-method
# covariance with the binding channel selected pointwise, and cross terms between
# grid points served by different channels vanish by independence.  Near a
# crossing it understates uncertainty, which is why the caller warns and
# recommends var_method = "bootstrap" there.
.sft_rmi_ecdf_cov <- function(x, times) {
  n <- length(x)
  k <- length(times)
  if (n < 2L) return(matrix(NA_real_, k, k))
  f <- .sft_rmi_ecdf_at(x, times)
  outer(f, f, pmin) * (1 - outer(f, f, pmax)) / n
}


.sft_rmi_analytic_cov <- function(cells, times, bound) {
  cov_ab <- .sft_rmi_ecdf_cov(cells$AB, times)
  cov_a <- .sft_rmi_ecdf_cov(cells$A, times)
  cov_b <- .sft_rmi_ecdf_cov(cells$B, times)
  if (bound == "miller") return(cov_ab + cov_a + cov_b)
  f_a <- .sft_rmi_ecdf_at(cells$A, times)
  f_b <- .sft_rmi_ecdf_at(cells$B, times)
  take_a <- f_a >= f_b
  k <- length(times)
  channel <- matrix(0, k, k)
  if (any(take_a)) channel[take_a, take_a] <- cov_a[take_a, take_a]
  if (any(!take_a)) channel[!take_a, !take_a] <- cov_b[!take_a, !take_a]
  cov_ab + channel
}


# Bootstrap sampling covariance.  Each condition is resampled with replacement
# and the whole estimator is recomputed, including the percentile grid: the grid
# is itself a statistic of the participant's data, so holding it fixed would
# understate uncertainty.  Because the estimand is indexed by the percentile
# level rather than by time, replicates remain comparable across a moving grid.
# `blocks` is the participant's list of per-block cell lists; resampling within
# block and re-averaging keeps the resampling scheme matched to the estimator.
.sft_rmi_boot_cov <- function(blocks, probs, bound, n_boot) {
  k <- length(probs)
  reps <- matrix(NA_real_, nrow = n_boot, ncol = k)
  for (b in seq_len(n_boot)) {
    per <- lapply(blocks, function(cells) {
      resampled <- lapply(cells, function(v) {
        if (length(v) < 2L) return(v)
        sort(v[sample.int(length(v), replace = TRUE)])
      })
      times <- .sft_rmi_grid(resampled, probs)
      if (anyNA(times)) return(NULL)
      list(delta = .sft_rmi_delta_at(resampled, times, bound)$delta)
    })
    if (any(vapply(per, is.null, logical(1)))) next
    reps[b, ] <- .sft_rmi_block_mean(per, "delta", k)
  }
  complete <- stats::complete.cases(reps)
  if (sum(complete) < max(10L, k + 1L)) return(matrix(NA_real_, k, k))
  stats::var(reps[complete, , drop = FALSE])
}


# One participant's delta vector, its covariance, and the grid points that carry
# information.
#
# Grid points where the Miller bound has already reached 1 are dropped: the
# inequality is vacuous there (F_AB <= 1 always) and the sampling variance
# degenerates.  Points with zero estimated variance are dropped for the same
# reason.  Dropped points are recorded per participant rather than removed from
# the shared percentile index, so participants with different informative ranges
# still contribute to the common population curve.
#
# With several blocks the reported delta, times, and cell CDFs are averages over
# the blocks that survive `min_n`, and the covariance is the covariance of that
# average.  A grid point is kept only where the average is estimable in every
# surviving block.
.sft_rmi_subject_score <- function(blocks, probs, bound,
                                   var_method = "analytic", n_boot = 2000L,
                                   min_n = 20L, errors = "discard") {
  k <- length(probs)
  empty <- list(delta = rep(NA_real_, k), cov = matrix(NA_real_, k, k),
                keep = rep(FALSE, k), times = rep(NA_real_, k),
                F_AB = rep(NA_real_, k), F_A = rep(NA_real_, k),
                F_B = rep(NA_real_, k), bound_value = rep(NA_real_, k),
                n = c(AB = 0L, A = 0L, B = 0L), n_blocks = 0L,
                n_pool = NA_integer_, ok = FALSE, near_crossing = FALSE)
  usable <- .sft_rmi_usable_blocks(blocks, probs, min_n, errors = errors)
  empty$n <- usable$n
  if (!usable$n_blocks) return(empty)

  cells <- usable$cells
  n_blocks <- usable$n_blocks
  per <- lapply(cells, function(cl) {
    times <- .sft_rmi_grid(cl, probs)
    c(.sft_rmi_delta_at(cl, times, bound), list(times = times))
  })
  delta <- .sft_rmi_block_mean(per, "delta", k)
  bound_value <- .sft_rmi_block_mean(per, "bound", k)
  covariance <- if (var_method == "bootstrap") {
    .sft_rmi_boot_cov(cells, probs, bound, n_boot)
  } else {
    # Blocks are disjoint sets of trials, so Var(mean of B block deltas) is the
    # sum of the block covariances over B^2.
    Reduce(`+`, lapply(seq_len(n_blocks), function(j) {
      .sft_rmi_analytic_cov(cells[[j]], per[[j]]$times, bound)
    })) / n_blocks^2
  }
  variance <- if (all(is.na(covariance))) rep(NA_real_, k) else diag(covariance)
  keep <- is.finite(delta) & is.finite(variance) & variance > 0
  if (bound == "miller") keep <- keep & bound_value < 1
  # max(F_A, F_B) is non-smooth where the channel CDFs cross; flag it so the
  # caller can recommend the bootstrap covariance.
  near_crossing <- bound == "grice" &&
    any(vapply(seq_len(n_blocks), function(j) {
      n_j <- lengths(cells[[j]])
      any(keep & abs(per[[j]]$F_A - per[[j]]$F_B) <
            1 / max(n_j[["A"]], n_j[["B"]]))
    }, logical(1)))
  list(delta = delta, cov = covariance, keep = keep,
       times = .sft_rmi_block_mean(per, "times", k),
       F_AB = .sft_rmi_block_mean(per, "F_AB", k),
       F_A = .sft_rmi_block_mean(per, "F_A", k),
       F_B = .sft_rmi_block_mean(per, "F_B", k),
       bound_value = bound_value, n = usable$n, n_blocks = n_blocks,
       n_pool = usable$n_pool, ok = sum(keep) >= 2L,
       near_crossing = near_crossing)
}


# Score every participant for one bound and assemble the vector-valued analogue
# of the scorer contract used by the scalar group models: instead of one
# effect_hat and one precision per participant, Y carries one delta vector per
# participant and S the matching sampling covariance.
.sft_rmi_score_data <- function(input, probs, bound, var_method = "analytic",
                                n_boot = 2000L, min_n = 20L,
                                errors = "discard") {
  var_method <- .sft_var_method(var_method)
  subject <- input$subject
  scores <- lapply(seq_along(subject), function(i) {
    .sft_rmi_subject_score(.sft_rmi_subject_blocks(input, i), probs, bound,
                           var_method = var_method, n_boot = n_boot,
                           min_n = min_n, errors = errors)
  })
  ok <- vapply(scores, `[[`, logical(1), "ok")
  if (!any(ok)) {
    stop("No participant has an estimable ", bound, "-bound violation curve: ",
         "each participant needs at least ", min_n, " correct trials in each of ",
         "the AB, A, and B cells and at least two informative percentile ",
         "levels.", call. = FALSE)
  }
  if (!all(ok)) {
    warning("Dropping ", sum(!ok), " participant(s) without an estimable ",
            bound, "-bound violation curve: ",
            paste(subject[!ok], collapse = ", "), ".", call. = FALSE)
    scores <- scores[ok]
    subject <- subject[ok]
  }
  .sft_rmi_warn_inestimable(probs,
                            vapply(scores, function(s) as.numeric(s$n_pool),
                                   numeric(1)))
  if (var_method == "analytic" && any(vapply(scores, `[[`, logical(1),
                                            "near_crossing"))) {
    warning("The channel CDFs cross within the tested range, where max(F_A, F_B) ",
            "is not differentiable and the analytic Grice covariance ",
            "understates uncertainty. Consider var_method = \"bootstrap\".",
            call. = FALSE)
  }

  y <- do.call(rbind, lapply(scores, `[[`, "delta"))
  keep <- do.call(rbind, lapply(scores, `[[`, "keep"))
  covariances <- lapply(scores, `[[`, "cov")
  score_df <- do.call(rbind, lapply(seq_along(scores), function(i) {
    s <- scores[[i]]
    data.frame(subject = subject[[i]], prob = probs, time = s$times,
               F_A = s$F_A, F_B = s$F_B, F_AB = s$F_AB, bound = s$bound_value,
               delta = s$delta,
               se = sqrt(ifelse(rep(all(is.na(s$cov)), length(probs)),
                                NA_real_, diag(s$cov))),
               used = s$keep, n_blocks = s$n_blocks, stringsAsFactors = FALSE)
  }))
  rownames(score_df) <- NULL
  cell_n <- do.call(rbind, lapply(scores, `[[`, "n"))

  list(subject = subject, Y = y, S = covariances, keep = keep, probs = probs,
       bound = bound, effect_name = "delta", var_method = var_method,
       score = score_df, n = cell_n,
       n_blocks = vapply(scores, `[[`, integer(1), "n_blocks"))
}


# Miller-comparable secondary estimand on the time scale: the signed distance
# between the redundant-target quantile and the bound quantile at percentile p,
# reported raw and divided by the participant's grand-mean correct RT.  The
# relative version follows the same reasoning as the relative MIC in
# R/mic_bayes.R (rho = MIC / M): dividing by a participant scale makes the
# quantity dimensionless so participant base time cannot drive the group mean.
# Negative values are violations, matching the probability-scale orientation.
.sft_rmi_invert <- function(times, values, probs) {
  keep <- is.finite(times) & is.finite(values)
  times <- times[keep]
  values <- values[keep]
  if (length(times) < 2L) return(rep(NA_real_, length(probs)))
  ord <- order(times)
  times <- times[ord]
  values <- cummax(values[ord])
  usable <- !duplicated(values)
  if (sum(usable) < 2L) return(rep(NA_real_, length(probs)))
  out <- stats::approx(values[usable], times[usable], xout = probs,
                       ties = "ordered")$y
  out[probs > max(values) | probs < min(values)] <- NA_real_
  out
}


#
# Under cdf_method = "polygon", `qtype` selects the plotting-position rule used
# for both the observed quantiles and the bound; see .sft_rmi_positions(). The
# default 5 is Equation 2 of Ulrich, Miller & Schroeter (2007), with ties
# handled by their Appendix A. Under cdf_method = "kernel", qtype is irrelevant:
# A, B, AB, and the bound all use exact Gaussian kernel CDFs.
#
# With several blocks each quantile is estimated within block and averaged over
# blocks, as Miller (1982) did, and the difference is then taken between the
# averaged quantiles.  A percentile the bound never reaches in a given block
# contributes nothing there, so q_bound averages over the blocks that reach it.
.sft_rmi_time_subject <- function(blocks, probs, bound, min_n = 20L, qtype = 5,
                                  cdf_method = c("polygon", "kernel"),
                                  kernel_bw = "nrd0") {
  cdf_method <- match.arg(cdf_method)
  k <- length(probs)
  empty <- list(diff = rep(NA_real_, k), relative = rep(NA_real_, k),
                q_AB = rep(NA_real_, k), q_A = rep(NA_real_, k),
                q_B = rep(NA_real_, k), q_bound = rep(NA_real_, k),
                mean_rt = NA_real_, n = c(AB = 0L, A = 0L, B = 0L),
                n_blocks = 0L, n_min = NA_integer_, ok = FALSE)
  usable <- .sft_rmi_usable_blocks(blocks, probs, min_n, grid = FALSE)
  empty$n <- usable$n
  if (!usable$n_blocks) return(empty)

  per <- lapply(usable$cells, function(cells) {
    estimates <- if (cdf_method == "kernel") {
      .sft_rmi_kernel_estimates(cells, probs, bound, kernel_bw = kernel_bw)
    } else {
      quantile_of <- function(x) .sft_rmi_quantile(x, probs, qtype = qtype)
      curve <- .sft_rmi_bound_curve(cells, bound, qtype = qtype)
      list(q_bound = if (is.null(curve)) rep(NA_real_, k) else
             .sft_rmi_invert(curve$grid, curve$value, probs),
           q_AB = quantile_of(cells$AB), q_A = quantile_of(cells$A),
           q_B = quantile_of(cells$B))
    }
    c(estimates, list(rt = c(cells$AB, cells$A, cells$B)))
  })
  q_bound <- .sft_rmi_block_mean(per, "q_bound", k, na.rm = TRUE)
  q_ab <- .sft_rmi_block_mean(per, "q_AB", k, na.rm = TRUE)
  mean_rt <- mean(unlist(lapply(per, `[[`, "rt")))
  difference <- if (bound == "miller") q_ab - q_bound else q_bound - q_ab
  ok <- is.finite(mean_rt) && mean_rt > 0 && sum(is.finite(difference)) >= 2L
  list(diff = difference,
       relative = if (isTRUE(ok)) difference / mean_rt else rep(NA_real_, k),
       q_AB = q_ab, q_A = .sft_rmi_block_mean(per, "q_A", k, na.rm = TRUE),
       q_B = .sft_rmi_block_mean(per, "q_B", k, na.rm = TRUE),
       q_bound = q_bound, mean_rt = mean_rt, n = usable$n,
       n_blocks = usable$n_blocks, n_min = usable$n_min, ok = ok)
}


.sft_rmi_time_table <- function(input, probs, bound, min_n = 20L, qtype = 5,
                                cdf_method = c("polygon", "kernel"),
                                kernel_bw = "nrd0") {
  cdf_method <- match.arg(cdf_method)
  parts <- lapply(seq_along(input$subject), function(i) {
    .sft_rmi_time_subject(.sft_rmi_subject_blocks(input, i), probs, bound,
                          min_n = min_n, qtype = qtype,
                          cdf_method = cdf_method, kernel_bw = kernel_bw)
  })
  ok <- vapply(parts, `[[`, logical(1), "ok")
  if (!any(ok)) return(NULL)
  # Raising min_n to the 20 that Ulrich et al. and Kiesel et al. recommend
  # excludes participants that a smaller threshold would have admitted, so say
  # which ones rather than shrinking the table quietly.
  if (!all(ok)) {
    warning("Dropping ", sum(!ok), " participant(s) without an estimable ",
            bound, "-bound quantile table: ",
            paste(input$subject[!ok], collapse = ", "),
            ". Each participant needs at least min_n = ", min_n,
            " correct trials in each of the AB, A, and B cells.",
            call. = FALSE)
  }
  subject <- input$subject[ok]
  parts <- parts[ok]
  .sft_rmi_warn_inestimable(probs,
                            vapply(parts, function(p) as.numeric(p$n_min),
                                   numeric(1)))
  long <- do.call(rbind, lapply(seq_along(parts), function(i) {
    p <- parts[[i]]
    data.frame(subject = subject[[i]], prob = probs, q_A = p$q_A, q_B = p$q_B,
               q_AB = p$q_AB, q_bound = p$q_bound, difference = p$diff,
               relative = p$relative, mean_RT = p$mean_rt,
               n_AB = p$n[["AB"]], n_A = p$n[["A"]], n_B = p$n[["B"]],
               n_blocks = p$n_blocks, stringsAsFactors = FALSE)
  }))
  rownames(long) <- NULL
  summarise <- function(values) {
    do.call(rbind, lapply(seq_along(probs), function(k) {
      x <- values[, k]
      x <- x[is.finite(x)]
      n <- length(x)
      se <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
      data.frame(prob = probs[[k]], n = n, mean = if (n) mean(x) else NA_real_,
                 se = se,
                 lower = if (n > 1L) mean(x) - stats::qt(.975, n - 1L) * se else NA_real_,
                 upper = if (n > 1L) mean(x) + stats::qt(.975, n - 1L) * se else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }
  diff_matrix <- do.call(rbind, lapply(parts, `[[`, "diff"))
  rel_matrix <- do.call(rbind, lapply(parts, `[[`, "relative"))
  # Quantile averaging across participants at each percentile level -- Miller's
  # (1982) "Vincentized" group CDF, which is what his Fig. 1 plots. It is a
  # display convention, not the estimand: the group mean of a millisecond
  # quantile is base-time dependent, which is exactly why the tested quantity is
  # the within-participant probability-scale delta.
  #
  vincentized <- function(column) {
    vapply(seq_along(probs), function(k) {
      x <- long[[column]][long$prob == probs[[k]]]
      x <- x[is.finite(x)]
      if (length(x)) mean(x) else NA_real_
    }, numeric(1))
  }
  list(subject = long, difference = summarise(diff_matrix),
       relative = summarise(rel_matrix),
       quantiles = data.frame(prob = probs, q_A = vincentized("q_A"),
                              q_B = vincentized("q_B"),
                              q_AB = vincentized("q_AB"),
                              q_bound = vincentized("q_bound"),
                              stringsAsFactors = FALSE),
       note = paste("Secondary time-scale summary for comparability with",
                    "Miller (1982). Negative values indicate violation.",
                    "Percentiles the bound never reaches are reported as NA."))
}


#' Race model inequality tests (Miller 1982; Ulrich, Miller & Schroeter 2007).
#'
#' \code{rmi.test} reproduces the classical procedure: per-participant quantiles
#' of the bound distribution and of the redundant-target distribution are
#' compared at each percentile level with a one-sided t-test across
#' participants.  It is provided so that a replication can run the classical and
#' the hierarchical Bayesian analyses (\code{\link{rmiGroup.bayes}}) on identical
#' data through one input path.
#'
#' @param inData canonical trial data frame with \code{Subject}, \code{RT},
#'   \code{Correct}, \code{Channel1}, \code{Channel2} and optionally
#'   \code{Condition}, or a pre-split RT/CR list.
#' @param bound \code{"miller"} (upper bound, the race model inequality),
#'   \code{"grice"} (lower bound), or \code{"both"}.
#' @param probs percentile levels to test.  The default is the 10\%-25\% range
#'   recommended by Kiesel, Miller and Ulrich (2007): testing the whole
#'   distribution accumulates Type I error to roughly 13\%, percentiles above
#'   the median are vacuous for the Miller bound, and 10\%-25\% is where
#'   violations are in practice found, so restricting the range costs little
#'   power.  Widen it deliberately, not by habit.
#' @param p.adjust.method multiplicity adjustment passed to
#'   \code{\link[stats]{p.adjust}}; the default \code{"holm"} is applied because
#'   the classical procedure applies none.  Kiesel et al. note that
#'   Bonferroni-type corrections are conservative here, because the tests at
#'   adjacent percentiles correlate .77 to .95, and prefer range restriction as
#'   the primary control; \code{p.adjust.method = "none"} recovers their
#'   procedure once \code{probs} is restricted.
#' @param Condition,Subject optional selectors.
#' @param min_n minimum correct trials required per cell; with
#'   \code{by_block = TRUE} it applies per participant, block, and cell.  The
#'   default of 10 is the sample size Ulrich et al. and Kiesel et al. recommend
#'   if, as in Miller (1982), there are at least 2 blocks per subject;
#'   below it the systematic bias in percentile estimation is large enough to
#'   produce spurious violations, and Kiesel et al. advise treating rejections
#'   obtained with fewer than 20 trials in any cell with suspicion.
#' @param by_block estimate the quantiles separately within each participant
#'   \emph{and} block, using the \code{Block} column of \code{inData}, and
#'   average over blocks, as Miller (1982) did.  \code{FALSE} (the default)
#'   pools the whole session and estimates one set of quantiles per participant
#'   (and condition).
#' @param qtype plotting-position rule for the cumulative frequency polygon
#'   from which both the reported quantiles and the bound are estimated.  The
#'   default 5 places the ith of n ordered RTs at \eqn{(i - .5)/n}, which is
#'   Equation 2 of Ulrich, Miller and Schroeter (2007) and the convention
#'   Miller (1982) used; ties follow their Appendix A whatever the rule.
#' @param cdf_method \code{"polygon"} (default) for the Ulrich et al.
#'   cumulative-frequency polygon, or \code{"kernel"} for a fully symmetric
#'   Gaussian kernel-CDF analysis. In the kernel analysis the A, B, and AB CDFs,
#'   their quantiles, and the bound all use the same estimator; \code{qtype} is
#'   then ignored.
#' @param kernel_bw bandwidth selector or one positive bandwidth for
#'   \code{cdf_method = "kernel"}. The default \code{"nrd0"} uses R's standard
#'   rule independently in each cell. Other selectors accepted by
#'   \code{\link[stats]{density}} may also be supplied.
#' @return A list with one element per requested bound, each containing the
#'   per-percentile \code{test} table, the per-participant \code{subject} table,
#'   the quantile-averaged \code{quantiles} table used by
#'   \code{\link{build_rmi_bound_df}}, and the relative (base-time free)
#'   summary.
#' @seealso \code{\link{rmiGroup.bayes}} for the hierarchical Bayesian analysis
#'   on the probability scale, and \code{\link{build_rmi_bound_df}} for the
#'   figure that accompanies this test.
#' @references
#' Miller, J. (1982). Divided attention: Evidence for coactivation with
#' redundant signals. \emph{Cognitive Psychology}, 14, 247-279.
#'
#' Ulrich, R., Miller, J., & Schroeter, H. (2007). Testing the race model
#' inequality: An algorithm and computer programs. \emph{Behavior Research
#' Methods}, 39, 291-302.
#'
#' Kiesel, A., Miller, J., & Ulrich, R. (2007). Systematic biases and Type I
#' error accumulation in tests of the race model inequality. \emph{Behavior
#' Research Methods}, 39, 539-551.
#' @export
rmi.test <- function(inData, bound = c("miller", "grice", "both"),
                     probs = seq(.10, .25, by = .05),
                     p.adjust.method = "holm",
                     Condition = NULL, Subject = NULL, min_n = 10L,
                     by_block = FALSE, qtype = 5,
                     cdf_method = c("polygon", "kernel"),
                     kernel_bw = "nrd0") {
  bound <- .sft_rmi_bound(if (missing(bound)) "miller" else bound)
  cdf_method <- match.arg(cdf_method)
  probs <- .sft_rmi_probs(probs)
  grouped <- .sft_rmi_input(inData, NULL, Condition, Subject, by_block)
  bounds <- if (bound == "both") c("miller", "grice") else bound

  one_bound <- function(which_bound) {
    per_condition <- lapply(grouped$levels, function(level) {
      table <- .sft_rmi_time_table(grouped$conditions[[level]], probs,
                                   which_bound, min_n = min_n, qtype = qtype,
                                   cdf_method = cdf_method,
                                   kernel_bw = kernel_bw)
      if (is.null(table)) return(NULL)
      values <- table$subject
      test <- do.call(rbind, lapply(probs, function(p) {
        x <- values$difference[values$prob == p]
        x <- x[is.finite(x)]
        n <- length(x)
        if (n < 2L || stats::sd(x) == 0) {
          return(data.frame(prob = p, n = n, mean_difference = if (n) mean(x) else NA_real_,
                            t = NA_real_, df = NA_real_, p.value = NA_real_,
                            stringsAsFactors = FALSE))
        }
        # One-sided: the alternative of interest is a violation, i.e. the
        # redundant-target quantile faster than the bound quantile.
        tt <- stats::t.test(x, alternative = "less")
        data.frame(prob = p, n = n, mean_difference = mean(x),
                   t = unname(tt$statistic), df = unname(tt$parameter),
                   p.value = tt$p.value, stringsAsFactors = FALSE)
      }))
      test$p.adjusted <- stats::p.adjust(test$p.value, method = p.adjust.method)
      test$violation <- !is.na(test$p.adjusted) & test$p.adjusted < .05
      rownames(test) <- NULL
      list(test = test, subject = table$subject,
           relative = table$relative, difference = table$difference,
           quantiles = table$quantiles, bound = which_bound,
           condition = level, by_block = isTRUE(grouped$by_block),
           cdf_method = cdf_method, kernel_bw = kernel_bw,
           method = paste(if (cdf_method == "kernel")
                            "Log-kernel CDF test of the"
                          else "Ulrich-Miller-Schroeter test of the",
                          which_bound, "bound"),
           p.adjust.method = p.adjust.method)
    })
    names(per_condition) <- grouped$levels
    per_condition <- per_condition[!vapply(per_condition, is.null, logical(1))]
    if (!length(per_condition)) {
      stop("No condition yielded an estimable ", which_bound, "-bound test.",
           call. = FALSE)
    }
    if (!grouped$multi_condition) per_condition[[1L]] else per_condition
  }

  out <- lapply(bounds, one_bound)
  names(out) <- bounds
  if (length(out) == 1L) out[[1L]] else out
}


# ---- Hierarchical Bayesian race-model inequality ---------------------------
#
# The scorer above is vector valued: one delta per percentile level per
# participant.  The hierarchy that consumes it treats the (participant,
# percentile) cell as the unit and gives each percentile level its own
# population mean over a shared between-participant scale,
#
#   delta_hat[i,k] ~ Normal(delta[i,k], 1 / sqrt(R[i,k]))
#   delta[i,k]     ~ Normal(mu[k], tau)
#   mu[k]          ~ Normal(prior_mean, prior_sd)
#
# which is the same normal-normal engine the other group models use, with the
# percentile level playing the role the Condition factor plays there.  Two
# consequences are worth stating plainly.  First, tau is a single
# between-participant scale assumed common across the grid; participants are
# taken to be independent given (mu, tau), so the within-participant covariance
# across percentile levels carried in S is used for the marginal variances only
# and does not induce dependence in the fit.  Second, a participant whose
# informative range covers only part of the grid contributes to the levels it
# covers and is silently absent elsewhere, which is why the levels are indexed by
# percentile rather than by time.
#
# Orientation is fixed by .sft_rmi_delta_at(): mu[k] < 0 is a violation of the
# requested bound at level k, for both the Miller and the Grice bound.
#
# The global summary is min_k mu[k].  P(min_k mu[k] < 0) is the posterior
# probability that the population violation curve dips below zero *somewhere on
# the tested grid*, which is the multiplicity-honest analogue of the classical
# procedure's family of per-percentile t-tests.  It is a statement about the
# tested levels only: a violation confined to percentiles outside `probs` is not
# in its scope.


# Flatten the vector-valued score into the scalar (effect, precision) contract
# used by .sft_normal_group_result(), one row per kept (participant, percentile)
# cell, with the group index running over (condition, percentile) combinations.
.sft_rmi_stack <- function(grouped, probs, bound, var_method = "analytic",
                           n_boot = 2000L, min_n = 20L, errors = "discard") {
  parts <- lapply(grouped$levels, function(level) {
    .sft_rmi_score_data(grouped$conditions[[level]], probs, bound,
                        var_method = var_method, n_boot = n_boot,
                        min_n = min_n, errors = errors)
  })
  names(parts) <- grouped$levels
  multi <- grouped$multi_condition
  labels <- .sft_rmi_prob_label(probs)

  cells <- do.call(rbind, lapply(seq_along(parts), function(ci) {
    p <- parts[[ci]]
    idx <- which(p$keep, arr.ind = TRUE)
    if (!nrow(idx)) return(NULL)
    variance <- vapply(seq_len(nrow(idx)), function(r) {
      p$S[[idx[r, 1L]]][idx[r, 2L], idx[r, 2L]]
    }, numeric(1))
    data.frame(condition = grouped$levels[[ci]],
               subject = p$subject[idx[, 1L]],
               prob = probs[idx[, 2L]],
               level = if (multi) {
                 paste0(grouped$levels[[ci]], "@", labels[idx[, 2L]])
               } else labels[idx[, 2L]],
               effect_hat = p$Y[idx], variance = variance,
               stringsAsFactors = FALSE)
  }))
  if (is.null(cells) || !nrow(cells)) {
    stop("No estimable (participant, percentile) cell remains.", call. = FALSE)
  }
  cells <- cells[is.finite(cells$effect_hat) & is.finite(cells$variance) &
                   cells$variance > 0, , drop = FALSE]
  ord <- order(match(cells$condition, grouped$levels), cells$prob)
  cells <- cells[ord, , drop = FALSE]
  rownames(cells) <- NULL

  levels <- unique(cells$level)
  group <- match(cells$level, levels)
  # A percentile level informed by a single participant cannot separate mu[k]
  # from that participant's delta; drop it rather than let the shared tau
  # silently carry it.
  counts <- tabulate(group, length(levels))
  if (any(counts < 2L)) {
    drop <- levels[counts < 2L]
    warning("Dropping percentile level(s) informed by fewer than two ",
            "participants: ", paste(drop, collapse = ", "), ".", call. = FALSE)
    keep <- !(cells$level %in% drop)
    cells <- cells[keep, , drop = FALSE]
    levels <- unique(cells$level)
    group <- match(cells$level, levels)
  }
  if (length(levels) < 2L) {
    stop("Fewer than two percentile levels have at least two participants; ",
         "the race-model hierarchy is not identified.", call. = FALSE)
  }

  map <- cells[!duplicated(cells$level), c("condition", "prob", "level"),
               drop = FALSE]
  map$column <- seq_len(nrow(map))
  rownames(map) <- NULL
  score <- do.call(rbind, lapply(seq_along(parts), function(ci) {
    df <- parts[[ci]]$score
    if (multi) cbind(Condition = grouped$levels[[ci]], df,
                     stringsAsFactors = FALSE) else df
  }))
  rownames(score) <- NULL

  list(effect_hat = cells$effect_hat, precision = 1 / cells$variance,
       subject = paste0(cells$subject, "@", cells$level),
       raw_subject = cells$subject, group = group, levels = levels,
       multi_condition = TRUE, effect_name = "delta", score = score,
       cells = cells, map = map, per_condition = parts,
       conditions = grouped$levels, bound = bound)
}


# Per-percentile population summary and the min-over-grid global statement.
.sft_rmi_population <- function(mu_matrix, map, hdi, rope) {
  rows <- lapply(seq_len(nrow(map)), function(k) {
    d <- mu_matrix[, map$column[[k]]]
    interval <- .sft_hdi(d, hdi)
    data.frame(condition = map$condition[[k]], prob = map$prob[[k]],
               mean = mean(d), median = stats::median(d),
               lower = interval[["lower"]], upper = interval[["upper"]],
               p_violation = mean(d < 0),
               p_rope = if (is.null(rope)) NA_real_ else mean(abs(d) <= rope),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


.sft_rmi_global <- function(mu_matrix, map, hdi) {
  conditions <- unique(map$condition)
  rows <- lapply(conditions, function(cond) {
    cols <- map$column[map$condition == cond]
    probs <- map$prob[map$condition == cond]
    m <- if (length(cols) == 1L) mu_matrix[, cols] else
      apply(mu_matrix[, cols, drop = FALSE], 1L, min)
    worst <- if (length(cols) == 1L) rep(1L, length(m)) else
      max.col(-mu_matrix[, cols, drop = FALSE], ties.method = "first")
    interval <- .sft_hdi(m, hdi)
    data.frame(condition = cond, mean_worst = mean(m),
               lower = interval[["lower"]], upper = interval[["upper"]],
               p_any_violation = mean(m < 0),
               p_bound_holds = mean(m >= 0),
               median_prob_of_worst = stats::median(probs[worst]),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# Between-condition contrasts evaluated separately at each percentile level.
# The generic pairwise-mu contrasts computed by .sft_normal_group_result() are
# not usable here: its levels index (condition, percentile) jointly, so most of
# its pairs compare different percentiles of different conditions.
.sft_rmi_contrasts <- function(mu_matrix, map, hdi, rope, prior_sd) {
  conditions <- unique(map$condition)
  if (length(conditions) < 2L) return(NULL)
  pairs <- utils::combn(conditions, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    shared <- intersect(map$prob[map$condition == pair[[1L]]],
                        map$prob[map$condition == pair[[2L]]])
    if (!length(shared)) return(NULL)
    do.call(rbind, lapply(shared, function(p) {
      a <- mu_matrix[, map$column[map$condition == pair[[1L]] & map$prob == p]]
      b <- mu_matrix[, map$column[map$condition == pair[[2L]] & map$prob == p]]
      d <- a - b
      interval <- .sft_hdi(d, hdi)
      bf <- if (is.null(rope)) NULL else
        .sft_interval_null_bf(d, rope, 0, sqrt(2) * prior_sd)
      data.frame(contrast = paste(pair[[1L]], "-", pair[[2L]]), prob = p,
                 mean = mean(d), lower = interval[["lower"]],
                 upper = interval[["upper"]],
                 posterior_positive = mean(d > 0),
                 posterior_negative = mean(d < 0),
                 BF01_interval = if (is.null(bf)) NA_real_ else bf$BF01,
                 BF10_interval = if (is.null(bf)) NA_real_ else bf$BF10,
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) return(NULL)
  rownames(out) <- NULL
  out
}


#' Hierarchical Bayesian race model inequality test.
#'
#' Fits the race-model inequality on the probability scale, where it is free of
#' participant base time, and reports a population violation curve indexed by
#' percentile level together with a single min-over-grid global summary.
#'
#' @details
#' For participant \eqn{i} and percentile level \eqn{u_k}, the evaluation time
#' \eqn{t_{ik}} is the \eqn{u_k}-th quantile of that participant's own pooled
#' correct single-target RTs, and the estimand is
#' \deqn{\delta_{ik} = F_A(t_{ik}) + F_B(t_{ik}) - F_{AB}(t_{ik})}
#' for the Miller bound, or \eqn{F_{AB}(t_{ik}) - \max(F_A(t_{ik}),
#' F_B(t_{ik}))} for the Grice bound.  Negative values are violations under
#' either bound.  Because the inequality is evaluated within a participant,
#' overall speed cancels exactly and no millisecond quantity enters the group
#' statistic; the estimand is also invariant to any strictly increasing
#' reparameterisation of time applied to all three of a participant's
#' conditions, which is what licenses pooling participants on the shared
#' percentile index.  Averaging millisecond quantiles across participants
#' instead -- Miller's original procedure, and the arrangement
#' \code{\link{rmi.test}} reproduces -- lets a slow participant's violation
#' outweigh a fast participant's equally severe one, and can cancel real
#' violations across participants with different base times.
#'
#' The hierarchy gives each percentile level its own population mean
#' \eqn{\mu_k} over a between-participant scale \eqn{\tau} shared across the
#' grid.  Participants are independent given \eqn{(\mu, \tau)}, so the
#' within-participant covariance across levels enters through the marginal
#' variances only.  The global summary is \eqn{\min_k \mu_k}, whose posterior
#' probability of being negative is the multiplicity-honest counterpart of the
#' classical family of per-percentile t-tests; it concerns the tested grid only.
#'
#' @param inData canonical trial data frame with \code{Subject}, \code{RT},
#'   \code{Correct}, \code{Channel1}, \code{Channel2} and optionally
#'   \code{Condition}, or a pre-split RT/CR list.
#' @param CR optional correct-response list when \code{inData} is a list.
#' @param bound \code{"miller"} (upper bound, the race model inequality) or
#'   \code{"grice"} (lower bound).
#' @param probs percentile levels to evaluate.  The default is the 10\%-25\%
#'   range recommended by Kiesel, Miller and Ulrich (2007); see
#'   \code{\link{rmi.test}}.  The multiplicity argument for restricting it is
#'   weaker here, since the global summary is a single min-over-grid statement
#'   rather than a family of tests, but the bias argument is not: percentile
#'   estimates outside that range are the ones most affected by the systematic
#'   bias they document, and the levels above the median are vacuous for the
#'   Miller bound in any case.
#' @param ndraws,burnin,thin,chains sampler settings.
#' @param seed optional RNG seed.
#' @param hdi mass of the reported highest-density intervals.
#' @param prior_mean,prior_sd Normal prior on each population mean \eqn{\mu_k}.
#'   \code{prior_sd} defaults to 0.25: \eqn{\delta} is a difference of
#'   probabilities, so the default puts most prior mass inside +/- 0.5, which is
#'   already a very large violation.
#' @param prior_shape,prior_rate inverse-Gamma prior on \eqn{\tau^2} for
#'   \code{method = "InvGamma"}.
#' @param prior_tau_sd half-Normal prior scale on \eqn{\tau} for the Stan
#'   methods.
#' @param rope optional half-width of a region of practical equivalence on the
#'   probability scale.
#' @param method \code{"InvGamma"} (conjugate Gibbs), or \code{"HalfNormal"} /
#'   \code{"Centered"} (Stan NUTS).
#' @param var_method \code{"analytic"} for the closed-form multinomial
#'   covariance of the cell CDFs, or \code{"bootstrap"} to resample within
#'   participant, which also propagates uncertainty in the percentile grid.
#' @param n_boot bootstrap replicates when \code{var_method = "bootstrap"}.
#' @param adapt_delta,max_treedepth,stan_control Stan sampler controls.
#' @param Condition,Subject optional selectors.
#' @param min_n minimum correct trials required per cell; with
#'   \code{by_block = TRUE} it applies per participant, block, and cell.  The
#'   default of 10 is half the sample size Ulrich et al. and Kiesel et al. recommend because 
#'   there are usually 2+ blocks per subject;
#'   see \code{\link{rmi.test}}.
#' @param by_block evaluate the inequality separately within each participant
#'   \emph{and} block, using the \code{Block} column of \code{inData}, and
#'   average the per-block \eqn{\delta} curves, as Miller (1982) estimated his
#'   quantiles.  Blocks are disjoint sets of trials, so the sampling covariance
#'   of the average is the sum of the per-block covariances over the squared
#'   number of blocks.  \code{FALSE} (the default) pools the whole session.
#' @param qtype plotting-position rule for the cumulative frequency polygon
#'   behind the reported time-scale quantiles; see \code{\link{rmi.test}}.  It
#'   affects the accompanying quantile tables only, never the
#'   probability-scale estimand, which uses the plain empirical CDF.
#' @param cdf_method,kernel_bw estimator and bandwidth for the accompanying
#'   time-scale quantile tables; see \code{\link{rmi.test}}. These do not alter
#'   the Bayesian probability-scale estimand.
#' @param errors how incorrect and omitted trials enter the cell CDFs.
#'   \code{"discard"} keeps correct trials only, as Miller did and as
#'   \code{\link{estimate.bounds}} does; \code{"defective"} counts an error in
#'   the denominator without ever satisfying RT <= t, giving each cell the
#'   defective CDF whose asymptote is the accuracy rather than 1.  Prefer
#'   \code{"defective"} when accuracy differs across cells.
#' @return An \code{sft_bayes} list.  Beyond the standard hierarchical fields it
#'   carries \code{percentile}, the population violation curve;
#'   \code{global}, the min-over-grid summary; \code{percentile_contrasts},
#'   between-condition differences at each level; and \code{quantiles}, the
#'   time-scale tables consumed by \code{\link{build_rmi_bound_df}}.
#' @seealso \code{\link{rmi.test}}, \code{\link{build_rmi_bound_df}},
#'   \code{\link{build_rmi_violation_df}}.
#' @references
#' Miller, J. (1982). Divided attention: Evidence for coactivation with
#' redundant signals. \emph{Cognitive Psychology}, 14, 247-279.
#'
#' Ulrich, R., Miller, J., & Schroeter, H. (2007). Testing the race model
#' inequality: An algorithm and computer programs. \emph{Behavior Research
#' Methods}, 39, 291-302.
#'
#' Kiesel, A., Miller, J., & Ulrich, R. (2007). Systematic biases and Type I
#' error accumulation in tests of the race model inequality. \emph{Behavior
#' Research Methods}, 39, 539-551.
#'
#' Miller, J. (2016). Statistical facilitation and the redundant signals
#' effect: What are race and coactivation models? \emph{Attention, Perception,
#' & Psychophysics}, 78, 516-519.
#' @export
rmiGroup.bayes <- function(inData, CR = NULL, bound = c("miller", "grice"),
                           probs = seq(.10, .25, by = .05),
                           ndraws = 10000L, prior_shape = NULL,
                           prior_rate = NULL, seed = NULL, hdi = .94,
                           prior_mean = 0, prior_sd = NULL, burnin = 1000L,
                           thin = 1L, chains = 4L, rope = NULL,
                           method = c("InvGamma", "HalfNormal", "Centered"),
                           prior_tau_sd = NULL, adapt_delta = .95,
                           max_treedepth = 12L, stan_control = list(),
                           Condition = NULL, Subject = NULL,
                           var_method = c("analytic", "bootstrap"),
                           n_boot = 2000L, min_n = 10L, by_block = FALSE,
                           qtype = 5, errors = c("discard", "defective"),
                           cdf_method = c("polygon", "kernel"),
                           kernel_bw = "nrd0") {
  bound <- .sft_rmi_bound(if (missing(bound)) "miller" else bound)
  errors <- match.arg(errors)
  cdf_method <- match.arg(cdf_method)
  if (bound == "both") {
    stop("bound must be a single bound; fit 'miller' and 'grice' separately.",
         call. = FALSE)
  }
  probs <- .sft_rmi_probs(probs)
  var_method <- .sft_var_method(var_method)
  method <- .sft_bayes_method(method)
  if (!is.list(stan_control)) stop("stan_control must be a list.")
  .sft_validate_bayes_args(ndraws, burnin, thin, chains, prior_shape,
                           prior_rate, prior_mean, prior_sd, prior_tau_sd,
                           hdi, rope, adapt_delta, max_treedepth)
  .sft_validate_boot(var_method, n_boot)

  grouped <- .sft_rmi_input(inData, CR, Condition, Subject, by_block)
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  scored <- .sft_rmi_stack(grouped, probs, bound, var_method = var_method,
                           n_boot = n_boot, min_n = min_n, errors = errors)
  if (length(unique(scored$raw_subject)) < 2L) {
    stop("The race-model hierarchy needs at least two participants; with one ",
         "participant tau would measure spread across percentile levels ",
         "rather than across participants.", call. = FALSE)
  }

  # delta is a difference of probabilities, so the scale is known a priori and
  # the priors do not key off a data-dependent reference information the way the
  # UCIP score priors do.
  priors <- list(
    prior_sd = if (is.null(prior_sd)) 0.25 else prior_sd,
    prior_tau_sd = if (is.null(prior_tau_sd)) 0.25 else prior_tau_sd,
    prior_shape = if (is.null(prior_shape)) 2 else prior_shape,
    prior_rate = if (is.null(prior_rate)) 0.25^2 else prior_rate)

  quantiles <- lapply(grouped$levels, function(level) {
    .sft_rmi_time_table(grouped$conditions[[level]], probs, bound,
                        min_n = min_n, qtype = qtype,
                        cdf_method = cdf_method, kernel_bw = kernel_bw)
  })
  names(quantiles) <- grouped$levels
  quantiles <- quantiles[!vapply(quantiles, is.null, logical(1))]

  out <- .sft_normal_group_result(
    scored = scored, priors = priors, prior_mean = prior_mean, method = method,
    ndraws = ndraws, burnin = burnin, thin = thin, chains = chains, seed = seed,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    stan_control = stan_control, hdi = hdi, rope = rope,
    label = paste("race-model", bound, "bound"),
    alternative = paste("the population", bound,
                        "bound is violated at some tested percentile"),
    interpretation = paste("mu[k] is the population", bound,
                           "bound margin at percentile level k, on the",
                           "probability scale; mu[k] < 0 is a violation"),
    precision_name = "R",
    extra = list(bound = bound, probs = probs, var_method = var_method,
                 quantiles = quantiles, cells = scored$cells,
                 conditions = scored$conditions, min_n = min_n,
                 by_block = isTRUE(grouped$by_block), errors = errors,
                 cdf_method = cdf_method, kernel_bw = kernel_bw))

  mu_matrix <- out$mu_matrix
  map <- scored$map
  out$percentile <- .sft_rmi_population(mu_matrix, map, hdi, rope)
  out$global <- .sft_rmi_global(mu_matrix, map, hdi)
  out$percentile_contrasts <- .sft_rmi_contrasts(mu_matrix, map, hdi, rope,
                                                 priors$prior_sd)
  # Replace the generic pairwise-mu contrasts: its levels index (condition,
  # percentile) jointly, so most of its pairs are not interpretable.
  out$condition_contrasts <- out$percentile_contrasts
  out$conditions <- scored$conditions
  out$statistic <- stats::setNames(out$global$p_any_violation,
                                   paste0("P(violation | ",
                                          out$global$condition, ")"))
  out$map <- map
  out$note <- paste(
    "Percentile levels where the bound is already vacuous are not fitted.",
    "For the Miller bound F_A + F_B reaches 1 near the median of the",
    "single-target distributions, above which F_AB <= 1 satisfies the",
    "inequality for free; the fit therefore covers the fast half of the",
    "distribution, which is where race-model violations occur.")
  out
}


# ---- Tidy builders for race-model figures ----------------------------------
#
# The package keeps ggplot2 in Suggests and never depends on it in R/, so these
# return tidy data frames in the same spirit as build_cdf_df() and friends
# rather than plot objects.


# Walk either a whole rmi.test()/rmiGroup.bayes() result or one of its nested
# per-bound / per-condition elements and collect every quantile table found,
# labelled by bound and condition.
.sft_rmi_collect <- function(x, bound = NULL, condition = NULL) {
  found <- list()
  visit <- function(node, bound_label, condition_label) {
    if (!is.list(node)) return(invisible(NULL))
    bound_label <- if (!is.null(node$bound) && is.character(node$bound))
      node$bound[[1L]] else bound_label
    if (!is.null(node$condition) && is.character(node$condition)) {
      condition_label <- node$condition[[1L]]
    }
    q <- node$quantiles
    if (is.data.frame(q) && all(c("prob", "q_AB", "q_bound") %in% names(q))) {
      found[[length(found) + 1L]] <<- data.frame(
        bound = bound_label %||% NA_character_,
        condition = condition_label %||% "default", q,
        stringsAsFactors = FALSE)
      return(invisible(NULL))
    }
    if (is.list(q) && !is.data.frame(q)) {
      for (nm in names(q)) visit(q[[nm]], bound_label, nm)
      return(invisible(NULL))
    }
    for (nm in names(node)) {
      if (nm %in% c("draws", "mu_matrix", "score", "cells", "subject",
                    "stan_fit", "shrinkage_draws")) next
      visit(node[[nm]], bound_label, condition_label)
    }
    invisible(NULL)
  }
  visit(x, NULL, NULL)
  if (!length(found)) {
    stop("No quantile table found; supply the result of rmi.test() or ",
         "rmiGroup.bayes().", call. = FALSE)
  }
  out <- .sft_bind_rows(found)
  out <- out[!duplicated(out[c("bound", "condition", "prob")]), , drop = FALSE]
  if (!is.null(bound)) out <- out[out$bound %in% bound, , drop = FALSE]
  if (!is.null(condition)) {
    out <- out[out$condition %in% condition, , drop = FALSE]
  }
  if (!nrow(out)) stop("No quantile table matched the requested bound / ",
                       "condition.", call. = FALSE)
  rownames(out) <- NULL
  out
}


# Keep an explicitly requested percentile window when a builder receives an
# already-computed result.  Raw trial data are passed through to rmi.test(),
# which performs the full percentile validation; this helper only needs to
# select levels that are present in a fitted result.
.sft_rmi_filter_probs <- function(x, probs) {
  if (is.null(probs)) return(x)
  requested <- sort(unique(as.numeric(probs)))
  if (!length(requested) || any(!is.finite(requested)) ||
      any(requested <= 0) || any(requested >= 1)) {
    stop("probs must be finite percentile levels strictly inside (0, 1).",
         call. = FALSE)
  }
  if (!"prob" %in% names(x)) return(x)
  available <- unique(x$prob[is.finite(x$prob)])
  matched <- vapply(requested, function(p) {
    any(abs(available - p) <= 1e-8)
  }, logical(1))
  if (any(!matched)) {
    stop("Requested percentile level(s) are not present in the supplied ",
         "result: ", paste(format(requested[!matched]), collapse = ", "),
         ". Pass raw trial data to recompute them.", call. = FALSE)
  }
  keep <- vapply(x$prob, function(p) {
    any(abs(p - requested) <= 1e-8)
  }, logical(1))
  x[keep, , drop = FALSE]
}


#' Tidy data frames for race model inequality figures.
#'
#' \code{build_rmi_bound_df} produces the classical race-model figure: the
#' redundant-target cumulative distribution plotted against the bound it must
#' not cross, so that a violation is visible as the AB curve lying to the left
#' of the \eqn{F_A + F_B} bound.  It is the visual complement to the
#' per-percentile t-tests of \code{\link{rmi.test}}, which report where that gap
#' is reliable but not how large it is.  \code{build_rmi_cdf_df} gives the three
#' condition CDFs behind it (Miller, 1982, Fig. 1), and
#' \code{build_rmi_violation_df} the violation curve with its interval.
#'
#' @details
#' All three consume the result of \code{\link{rmi.test}} or
#' \code{\link{rmiGroup.bayes}}, or any single per-bound / per-condition element
#' of one, and return a long data frame ready for \code{ggplot2}.
#'
#' The quantiles in \code{build_rmi_cdf_df} and \code{build_rmi_bound_df} are
#' averaged across participants at each percentile level -- Miller's
#' "Vincentized" group curves.  That average is on the millisecond scale and so
#' is base-time dependent; it is a display convention, and the quantity actually
#' tested by \code{rmiGroup.bayes} is the within-participant probability-scale
#' margin returned by \code{build_rmi_violation_df}.  Read the first two as the
#' picture and the third as the inference.
#'
#' \code{build_rmi_violation_df} reports the probability-scale population curve
#' with its highest-density interval and \code{p_violation} when given a
#' Bayesian fit, and the time-scale mean difference with its t-interval and
#' adjusted p-value when given \code{rmi.test} output.  Negative values are
#' violations on both scales.
#'
#' @param x result of \code{\link{rmi.test}} or \code{\link{rmiGroup.bayes}}.
#' @param bound optional bound to keep when the result holds both.
#' @param condition optional condition(s) to keep.
#' @param probs optional percentile levels. For raw trial data these are passed
#'   to \code{rmi.test}; for an existing result they select levels already
#'   present in that result. If omitted, the \code{rmi.test} default is used
#'   for raw data and all levels are retained from an existing result.
#' @param scale for \code{build_rmi_violation_df}, \code{"auto"} (probability
#'   scale if the fit is Bayesian, otherwise time), \code{"probability"}, or
#'   \code{"time"}.
#' @param by_block passed to \code{\link{rmi.test}} when \code{x} is a raw
#'   trial dataset: estimate the quantiles per participant and block rather than
#'   per participant.
#' @param qtype plotting-position rule passed to \code{\link{rmi.test}} when
#'   \code{x} is a raw trial dataset.
#' @param cdf_method,kernel_bw CDF estimator and kernel bandwidth passed to
#'   \code{\link{rmi.test}} when \code{x} is a raw trial dataset.
#' @return A long \code{data.frame}.
#' @references
#' Miller, J. (1982). Divided attention: Evidence for coactivation with
#' redundant signals. \emph{Cognitive Psychology}, 14, 247-279.
#' @name build_rmi_df
#' @export
build_rmi_cdf_df <- function(x, bound = NULL, condition = NULL,
                             by_block = FALSE, qtype = 5, probs = NULL,
                             cdf_method = c("polygon", "kernel"),
                             kernel_bw = "nrd0") {
  cdf_method <- match.arg(cdf_method)
  if (is.data.frame(x) && !all(c("prob", "q_AB", "q_bound") %in% names(x))) {
    args <- list(inData = x,
                 bound = if (is.null(bound)) "miller" else bound,
                 by_block = by_block, qtype = qtype,
                 cdf_method = cdf_method, kernel_bw = kernel_bw)
    if (!is.null(probs)) args$probs <- probs
    x <- do.call(rmi.test, args)
  }
  q <- .sft_rmi_collect(x, bound, condition)
  q <- .sft_rmi_filter_probs(q, probs)
  cells <- c(A = "q_A", B = "q_B", AB = "q_AB")
  out <- .sft_bind_rows(lapply(names(cells), function(cell) {
    data.frame(bound = q$bound, condition = q$condition, prob = q$prob,
               Condition = cell, RT = q[[cells[[cell]]]],
               stringsAsFactors = FALSE)
  }))
  out$Condition <- factor(out$Condition, levels = names(cells))
  out <- out[is.finite(out$RT), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' @rdname build_rmi_df
#' @export
build_rmi_bound_df <- function(x, bound = NULL, condition = NULL,
                               by_block = FALSE, qtype = 5, probs = NULL,
                               cdf_method = c("polygon", "kernel"),
                               kernel_bw = "nrd0") {
  cdf_method <- match.arg(cdf_method)
  if (is.data.frame(x) && !all(c("prob", "q_AB", "q_bound") %in% names(x))) {
    args <- list(inData = x,
                 bound = if (is.null(bound)) "miller" else bound,
                 by_block = by_block, qtype = qtype,
                 cdf_method = cdf_method, kernel_bw = kernel_bw)
    if (!is.null(probs)) args$probs <- probs
    x <- do.call(rmi.test, args)
  }
  q <- .sft_rmi_collect(x, bound, condition)
  q <- .sft_rmi_filter_probs(q, probs)
  bound_label <- ifelse(q$bound %in% "grice", "max(A, B) bound", "A + B bound")
  out <- rbind(
    data.frame(bound = q$bound, condition = q$condition, prob = q$prob,
               curve = "AB", RT = q$q_AB, stringsAsFactors = FALSE),
    data.frame(bound = q$bound, condition = q$condition, prob = q$prob,
               curve = bound_label, RT = q$q_bound, stringsAsFactors = FALSE))
  out <- out[is.finite(out$RT), , drop = FALSE]
  out <- out[order(out$bound, out$condition, out$curve, out$prob), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' @rdname build_rmi_df
#' @export
build_rmi_violation_df <- function(x, bound = NULL, condition = NULL,
                                   scale = c("auto", "probability", "time"),
                                   by_block = FALSE, qtype = 5, probs = NULL,
                                   cdf_method = c("polygon", "kernel"),
                                   kernel_bw = "nrd0") {
  scale <- match.arg(scale)
  cdf_method <- match.arg(cdf_method)
  if (is.data.frame(x) && !all(c("prob", "mean", "lower", "upper") %in% names(x))) {
    args <- list(inData = x,
                 bound = if (is.null(bound)) "miller" else bound,
                 by_block = by_block, qtype = qtype,
                 cdf_method = cdf_method, kernel_bw = kernel_bw)
    if (!is.null(probs)) args$probs <- probs
    x <- do.call(rmi.test, args)
  }
  bayesian <- is.data.frame(x$percentile)
  if (scale == "probability" && !bayesian) {
    stop("The probability-scale violation curve is produced by ",
         "rmiGroup.bayes(); rmi.test() reports the time scale.", call. = FALSE)
  }
  if (bayesian && scale != "time") {
    out <- x$percentile
    out <- .sft_rmi_filter_probs(out, probs)
    if (!is.null(condition)) out <- out[out$condition %in% condition, , drop = FALSE]
    out$scale <- "probability"
    out$bound <- x$bound %||% NA_character_
    rownames(out) <- NULL
    return(out)
  }
  collect <- function(node, bound_label, condition_label, acc) {
    if (!is.list(node)) return(acc)
    if (!is.null(node$bound) && is.character(node$bound)) bound_label <- node$bound[[1L]]
    if (!is.null(node$condition) && is.character(node$condition)) {
      condition_label <- node$condition[[1L]]
    }
    if (is.data.frame(node$test) && is.data.frame(node$difference)) {
      merged <- merge(node$difference[c("prob", "n", "mean", "lower", "upper")],
                      node$test[c("prob", "t", "df", "p.value", "p.adjusted")],
                      by = "prob", all = TRUE)
      acc[[length(acc) + 1L]] <- data.frame(
        bound = bound_label %||% NA_character_,
        condition = condition_label %||% "default", merged,
        scale = "time", stringsAsFactors = FALSE)
      return(acc)
    }
    for (nm in names(node)) {
      if (nm %in% c("subject", "quantiles", "relative", "difference", "test")) next
      acc <- collect(node[[nm]], bound_label, nm, acc)
    }
    acc
  }
  parts <- collect(x, NULL, NULL, list())
  if (!length(parts)) {
    stop("No per-percentile test table found; supply rmi.test() output.",
         call. = FALSE)
  }
  out <- .sft_bind_rows(parts)
  out <- .sft_rmi_filter_probs(out, probs)
  if (!is.null(bound)) out <- out[out$bound %in% bound, , drop = FALSE]
  if (!is.null(condition)) out <- out[out$condition %in% condition, , drop = FALSE]
  rownames(out) <- NULL
  out
}
