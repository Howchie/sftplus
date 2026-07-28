# Survivor and mean interaction contrasts: frequentist and Bayes-factor SIC/MIC.


sic.test <- function(HH, HL, LH, LL, method = "ks") {
  method <- match.arg(method, "ks")
  xs <- list(HH = HH, HL = HL, LH = LH, LL = LL)
  if (any(vapply(xs, function(x) !any(is.finite(x)), logical(1)))) {
    stop("All four SIC cells must contain at least one finite RT.")
  }
  times <- sort(unique(unlist(lapply(xs, function(x) x[is.finite(x)]), use.names = FALSE)))
  sic_values <- .sft_ecdf(LH)(times) + .sft_ecdf(HL)(times) -
    .sft_ecdf(HH)(times) - .sft_ecdf(LL)(times)
  n_eff <- 1 / sum(1 / vapply(xs, function(x) sum(is.finite(x)), numeric(1)))
  dplus <- max(0, sic_values)
  dminus <- abs(min(0, sic_values))
  pplus <- exp(-2 * n_eff * dplus^2)
  pminus <- exp(-2 * n_eff * dminus^2)
  dname <- paste("HH:", deparse(substitute(HH)), "HL:", deparse(substitute(HL)),
                 "LH:", deparse(substitute(LH)), "LL:", deparse(substitute(LL)))
  plus <- .sft_htest(setNames(dplus, "D+"), pplus, "the SIC is above 0 at some time",
                     "Houpt-Townsend KS-SIC test", dname)
  minus <- .sft_htest(setNames(dminus, "D-"), pminus, "the SIC is below 0 at some time",
                      "Houpt-Townsend KS-SIC test", dname)
  list(positive = plus, negative = minus)
}


mic.test <- function(HH, HL, LH, LL, method = c("art", "anova")) {
  method <- match.arg(method)
  xs <- lapply(list(HH, HL, LH, LL), function(x) as.numeric(x[is.finite(x)]))
  if (any(!lengths(xs))) stop("All four SIC cells must contain finite RTs.")
  statistic <- (mean(xs[[4]]) - mean(xs[[3]])) - (mean(xs[[2]]) - mean(xs[[1]]))
  names(statistic) <- "MIC"
  allrt <- unlist(xs, use.names = FALSE)
  h1 <- c(rep(1, length(xs[[1]]) + length(xs[[2]])), rep(0, length(xs[[3]]) + length(xs[[4]])))
  h2 <- c(rep(1, length(xs[[1]])), rep(0, length(xs[[2]])), rep(1, length(xs[[3]])), rep(0, length(xs[[4]])))
  if (method == "art") {
    mA0 <- mean(allrt[h1 == 0]); mA1 <- mean(allrt[h1 == 1])
    mB0 <- mean(allrt[h2 == 0]); mB1 <- mean(allrt[h2 == 1])
    adjusted <- round(allrt - (1 - h1) * mA0 - h1 * mA1 - (1 - h2) * mB0 - h2 * mB1, 15)
    fit <- stats::anova(stats::lm(rank(adjusted) ~ h1 * h2))
    p <- fit[["Pr(>F)"]][3]
    method_name <- "Adjusted Rank Transform test of the MIC"
  } else {
    old <- options(contrasts = c("contr.helmert", "contr.poly")); on.exit(options(old), add = TRUE)
    fit <- stats::anova(stats::lm(allrt ~ h1 * h2))
    p <- fit[["Pr(>F)"]][3]
    method_name <- "ANOVA test of the MIC"
  }
  .sft_htest(statistic, p, "the MIC is not zero", method_name,
             paste("HH:", deparse(substitute(HH)), "HL:", deparse(substitute(HL)),
                   "LH:", deparse(substitute(LH)), "LL:", deparse(substitute(LL))))
}


.sft_bayes_pair <- function(x, y, nbin = 20L, nsamp = 10000L, tol = 0, seed = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("Both distributions must contain finite RTs.")
  pooled <- c(x, y)
  nbin <- max(1L, min(as.integer(nbin), length(unique(pooled))))
  probs <- seq(0, 1, length.out = nbin + 1L)
  breaks <- unique(stats::quantile(pooled, probs = probs, names = FALSE, type = 8))
  if (length(breaks) < 2L) breaks <- c(breaks[[1L]], breaks[[1L]] + 1)
  breaks[length(breaks)] <- max(breaks[length(breaks)], max(pooled) + .Machine$double.eps)
  nbin <- length(breaks) - 1L
  alpha <- rep(1 / nbin, nbin)
  counts <- c(hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts,
              hist(y, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts)
  cx <- counts[seq_len(nbin)]; cy <- counts[nbin + seq_len(nbin)]
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  rdir <- function(n, a) {
    g <- matrix(stats::rgamma(n * length(a), shape = rep(a, each = n)), nrow = n)
    g / rowSums(g)
  }
  prior_x <- rdir(nsamp, alpha); prior_y <- rdir(nsamp, alpha)
  post_x <- rdir(nsamp, cx + alpha); post_y <- rdir(nsamp, cy + alpha)
  # dominates(a, b) is TRUE when the CDF of a is everywhere at or above that of
  # b, i.e. F_a >= F_b, which is S_a <= S_b. The reported relations are stated in
  # *survivor* terms, so "S_x > S_y" is dominates(y, x), not dominates(x, y).
  dominates <- function(a, b) rowSums(t(apply(a - b, 1, cumsum)) >= -tol) == ncol(a)
  prior_greater <- mean(dominates(prior_y, prior_x))
  post_greater <- mean(dominates(post_y, post_x))
  prior_less <- mean(dominates(prior_x, prior_y))
  post_less <- mean(dominates(post_x, post_y))
  list(BF_greater = post_greater / max(prior_greater, .Machine$double.eps),
       BF_less = post_less / max(prior_less, .Machine$double.eps),
       posterior_greater = post_greater, posterior_less = post_less,
       prior_greater = prior_greater, prior_less = prior_less,
       nbin = nbin, nsamp = nsamp)
}


siDominance <- function(HH, HL, LH, LL, method = c("ks", "dp"),
                        nbin = 20L, nsamp = 10000L, seed = NULL,
                        CR = NULL, errors = c("discard", "defective")) {
  method <- match.arg(method)
  errors <- .sft_errors_method(errors)
  cells <- .sft_sic_cells(HH, HL, LH, LL, CR = CR, errors = errors)
  if (errors == "defective") {
    # The unconditional route tests the extended assumption set: survivor
    # ordering (8B) as before, plus ordering of the correct subdistributions
    # (8A) and the reverse ordering of the error subdistributions (8D).
    return(.sft_dominance_defective(cells, method = method, nbin = nbin,
                                    nsamp = nsamp, seed = seed))
  }
  HH <- cells$HH$rt; HL <- cells$HL$rt; LH <- cells$LH$rt; LL <- cells$LL$rt
  pairs <- list(
    list(greater = "S.hh > S.hl", less = "S.hh < S.hl", x = HH, y = HL),
    list(greater = "S.hh > S.lh", less = "S.hh < S.lh", x = HH, y = LH),
    list(greater = "S.hl > S.ll", less = "S.hl < S.ll", x = HL, y = LL),
    list(greater = "S.lh > S.ll", less = "S.lh < S.ll", x = LH, y = LL)
  )
  if (method == "ks") {
    rows <- lapply(pairs, function(p) {
      # ks.test(alternative = "greater") tests whether the *CDF* of x lies above
      # that of y, which is S_x < S_y. The rows below are labelled in survivor
      # terms, so the alternatives are swapped relative to the label.
      greater <- suppressWarnings(stats::ks.test(p$x, p$y, alternative = "less", exact = FALSE))
      less <- suppressWarnings(stats::ks.test(p$x, p$y, alternative = "greater", exact = FALSE))
      list(test = c(p$greater, p$less), statistic = c(unname(greater$statistic), unname(less$statistic)),
           p = c(greater$p.value, less$p.value))
    })
    out <- data.frame(Test = unlist(lapply(rows, `[[`, "test")),
                      statistic = unlist(lapply(rows, `[[`, "statistic")),
                      p.value = unlist(lapply(rows, `[[`, "p")),
                      row.names = NULL)
    return(out)
  }
  rows <- lapply(seq_along(pairs), function(i) {
    p <- pairs[[i]]
    bf <- .sft_bayes_pair(p$x, p$y, nbin = nbin, nsamp = nsamp,
                          seed = if (is.null(seed)) NULL else seed + i)
    data.frame(Test = c(p$greater, p$less), BF = c(bf$BF_greater, bf$BF_less),
               posterior = c(bf$posterior_greater, bf$posterior_less),
               prior = c(bf$prior_greater, bf$prior_less),
               p.value = c(bf$posterior_greater, bf$posterior_less),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "inferential_type") <- "Bayesian posterior probability and posterior/prior event-probability ratio"
  out
}


# Classify one SIC curve into the six SFT shape classes. tolSIC is applied
# consistently to every sign test: an all-negative (N) or all-positive (P) curve
# means no excursion beyond tolSIC on the *other* side, not literal sign purity.
# Requiring strict purity makes N and P unreachable for any noisy curve -- a
# single bin a hair over zero disqualifies it -- which silently drains their
# posterior mass into the mixed classes.
.sft_dp_checkmods <- function(x, dx, tolSIC = 5e-2, tolMIC = 1e-2) {
  ans <- c(Z = FALSE, N = FALSE, P = FALSE, nP = FALSE, Np = FALSE, np = FALSE)
  pos <- x > tolSIC
  neg <- x < -tolSIC
  if (!any(pos) && !any(neg)) { ans["Z"] <- TRUE; return(ans) }
  if (!any(pos)) { ans["N"] <- TRUE; return(ans) }
  if (!any(neg)) { ans["P"] <- TRUE; return(ans) }
  # Both signs present: the only admissible mixed shapes cross once, negative
  # first (serial-AND and coactive).
  signed <- sign(x[pos | neg])
  if (signed[[1L]] > 0 || sum(diff(signed) != 0) > 1L) return(ans)
  neg_area <- -sum(x[neg] * dx[neg])
  pos_area <- sum(x[pos] * dx[pos])
  if (abs(neg_area - pos_area) < tolMIC) ans["np"] <- TRUE
  else if (neg_area < pos_area) ans["nP"] <- TRUE
  else ans["Np"] <- TRUE
  ans
}


sicDPtest <- function(dat, nbin = NULL, nsamp = 10000L, maxn = 500000L,
                      tolSIC = NULL, tolMIC = NULL, ci = .95, seed = NULL,
                      CR = NULL, errors = c("discard", "defective")) {
  if (!is.list(dat) || length(dat) != 4L) stop("dat must be a list in HH, HL, LH, LL order.")
  errors <- .sft_errors_method(errors)
  if (errors == "defective") {
    # Unconditional route: one Dirichlet per cell over correct and error bins,
    # yielding sign posteriors for each of the three contrasts. The default
    # tolerance is tighter than the correct-only default because the error and
    # agnostic contrasts are bounded by accuracy differences and so are an
    # order of magnitude smaller than a conditional SIC.
    cells <- .sft_sic_cells(dat[[1L]], dat[[2L]], dat[[3L]], dat[[4L]],
                            CR = CR, errors = errors)
    return(.sft_sic_dp_defective(cells, nbin = nbin, nsamp = nsamp, maxn = maxn,
                                 tolSIC = tolSIC %||% 1e-2, ci = ci, seed = seed))
  }
  tolSIC <- tolSIC %||% 5e-2
  if (!is.null(CR)) {
    dat <- lapply(.sft_sic_cells(dat[[1L]], dat[[2L]], dat[[3L]], dat[[4L]],
                                 CR = CR, errors = errors), function(z) z$rt)
  }
  dat <- lapply(dat, function(x) as.numeric(x[is.finite(x)]))
  if (any(!lengths(dat))) stop("All SIC cells must contain finite RTs.")
  pooled <- unlist(dat, use.names = FALSE)
  if (is.null(nbin)) nbin <- floor(min(lengths(dat)) / 2)
  nbin <- max(1L, min(as.integer(nbin), length(unique(pooled))))
  q <- stats::quantile(pooled, probs = seq(0, 1, length.out = nbin + 1L), names = FALSE, type = 8)
  bins <- unique(q)
  if (length(bins) < 2L) bins <- c(bins[[1L]], bins[[1L]] + 1)
  bins[length(bins)] <- max(bins[length(bins)], max(pooled) + .Machine$double.eps)
  nbin <- length(bins) - 1L
  dx <- diff(bins)
  counts <- matrix(vapply(dat, function(x)
    as.numeric(hist(x, breaks = bins, plot = FALSE, include.lowest = TRUE)$counts),
    numeric(nbin)), nrow = nbin, ncol = 4L,
    dimnames = list(NULL, c("HH", "HL", "LH", "LL")))
  alpha <- rep(1 / nbin, nbin)
  if (is.null(tolMIC)) tolMIC <- max(mean(pooled) / 700, 1e-8)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  rdir <- function(n, a) {
    g <- matrix(stats::rgamma(n * length(a), shape = rep(a, each = n)), nrow = n)
    g / rowSums(g)
  }
  # Each of the four cells carries its own Dirichlet over the shared bins, and a
  # draw of the SIC is the contrast of their CDFs,
  #   SIC = F_LH + F_HL - F_HH - F_LL,
  # the same definition sic.test() and sic() use. `dat` is in HH, HL, LH, LL
  # order, so the columns of `a` are weighted c(-1, +1, +1, -1).
  sic_draws <- function(n, a) {
    a <- matrix(a, nrow = nbin, ncol = 4L)
    cdf <- lapply(seq_len(4L), function(cc) {
      p <- rdir(n, a[, cc])
      if (nbin == 1L) matrix(p, nrow = n, ncol = 1L) else t(apply(p, 1L, cumsum))
    })
    cdf[[3L]] + cdf[[2L]] - cdf[[1L]] - cdf[[4L]]
  }
  classify <- function(n, a) {
    curves <- sic_draws(n, a)
    apply(curves, 1L, .sft_dp_checkmods, dx = dx, tolSIC = tolSIC, tolMIC = tolMIC)
  }
  prior <- matrix(0, nrow = 6, ncol = 1, dimnames = list(c("Z", "N", "P", "nP", "Np", "np"), "Prior"))
  post <- prior; N <- 0L
  repeat {
    prior <- prior + rowSums(classify(nsamp, matrix(alpha, nbin, 4L)))
    post <- post + rowSums(classify(nsamp, counts + alpha))
    N <- N + nsamp
    bf <- (post + 1) / (prior + 1)
    if (N >= maxn) break
    # The highest-BF class is the only one whose interval width is required for
    # the stopping rule, matching the original implementation's fast default.
    p_hi <- stats::qbeta(c((1 - ci) / 2, 1 - (1 - ci) / 2), post[which.max(bf)] + 1,
                         N - post[which.max(bf)] + 1)
    if (diff(p_hi) < 0.01) break
  }
  bf <- drop((post + 1) / (prior + 1)); names(bf) <- rownames(prior)
  bfp <- bf / sum(bf)
  list(BF = bf, BFp = bfp, N = rbind(Prior = drop(prior), Posterior = drop(post)),
       ci = ci, nbin = nbin, nsamp = nsamp, seed = seed,
       tolerances = c(SIC = tolSIC, MIC = tolMIC))
}


sictestBayes <- function(HH, HL, LH, LL, method = c("DP", "IG"), model = NULL,
                         nbin = NULL, nsamp = 10000L, maxn = 500000L,
                         seed = NULL, CR = NULL,
                         errors = c("discard", "defective"), ...) {
  method <- match.arg(method)
  errors <- .sft_errors_method(errors)
  if (method == "IG") {
    stop("The legacy IG/Stan SIC route is not bundled in sft_plus; use method='DP' or supply a separate model implementation.")
  }
  result <- sicDPtest(list(HH, HL, LH, LL), nbin = nbin, nsamp = nsamp,
                      maxn = maxn, seed = seed, CR = CR, errors = errors, ...)
  dname <- paste("HH:", deparse(substitute(HH)), "HL:", deparse(substitute(HL)),
                 "LH:", deparse(substitute(LH)), "LL:", deparse(substitute(LL)))
  if (errors == "defective") {
    # No architecture classification here: the sign of each contrast is
    # reported for the correct, error, and accuracy-agnostic curves separately.
    return(list(statistic = result$BF, BFp = result$BFp,
                posterior = result$posterior,
                method = "Nonparametric Bayesian unconditional SIC sign test",
                data.name = dname, errors = errors, details = result))
  }
  bf <- result$BF
  names(bf)[names(bf) == "Z"] <- "Zero"
  names(bf)[names(bf) == "N"] <- "Negative"
  names(bf)[names(bf) == "P"] <- "Positive"
  names(bf)[names(bf) == "np"] <- "NegPos.MIC0"
  names(bf)[names(bf) == "nP"] <- "NegPos.MICpos"
  names(bf)[names(bf) == "Np"] <- "NegPos.MICneg"
  list(statistic = bf, BFp = result$BFp, method = "Nonparametric Bayesian SIC test",
       data.name = dname, errors = errors, details = result)
}


sic <- function(HH, HL, LH, LL, domtest = "ks", sictest = "ks",
                mictest = c("art", "anova"), interpolate = FALSE,
                nbin = 20L, nsamp = 10000L, maxn = 500000L, seed = NULL,
                CR = NULL, errors = c("discard", "defective"), ...) {
  errors <- .sft_errors_method(errors)
  cells <- .sft_sic_cells(HH, HL, LH, LL, CR = CR, errors = errors)
  d <- siDominance(HH, HL, LH, LL, method = domtest, nbin = nbin, nsamp = nsamp,
                   seed = seed, CR = CR, errors = errors)
  if (errors == "defective") return(.sft_sic_defective_fit(
    cells, d, sictest = sictest, mictest = mictest, interpolate = interpolate,
    nbin = nbin, nsamp = nsamp, maxn = maxn, seed = seed))
  HH <- cells$HH$rt; HL <- cells$HL$rt; LH <- cells$LH$rt; LL <- cells$LL$rt
  xs <- list(HH = HH, HL = HL, LH = LH, LL = LL)
  times <- sort(unique(unlist(xs, use.names = FALSE)))
  sic_values <- .sft_ecdf(LH)(times) + .sft_ecdf(HL)(times) - .sft_ecdf(HH)(times) - .sft_ecdf(LL)(times)
  sic_fun <- if (interpolate) .sft_curve(times, sic_values) else stats::stepfun(times, c(0, sic_values))
  st <- if (sictest == "bf") sictestBayes(HH, HL, LH, LL, nbin = nbin, nsamp = nsamp, maxn = maxn, seed = seed) else sic.test(HH, HL, LH, LL, method = sictest)
  mt <- mic.test(HH, HL, LH, LL, method = mictest)
  n_eff <- 1 / sum(1 / vapply(xs, function(x) sum(is.finite(x)), numeric(1)))
  dvals <- if (sictest == "bf") NULL else rbind(c(unname(st$positive$statistic), st$positive$p.value),
                                                  c(unname(st$negative$statistic), st$negative$p.value))
  list(SIC = sic_fun, MIC = unname(mt$statistic), MICtest = mt,
       SICtest = st, Dominance = d, Dvals = dvals, N = n_eff,
       HH = .sft_ecdf(HH), HL = .sft_ecdf(HL), LH = .sft_ecdf(LH), LL = .sft_ecdf(LL),
       errors = errors)
}


# The unconditional counterpart of sic(). Returns the three interaction
# contrasts, per-cell accuracy, and the large-t asymptotes. No architecture
# classification is attempted: the mapping from these signs to an architecture
# needs the task type (simple detection, OR, or AND) and admits real mimicry
# between serial and parallel variants once error rates are non-negligible.
.sft_sic_defective_fit <- function(cells, dominance, sictest, mictest,
                                   interpolate, nbin, nsamp, maxn, seed) {
  fit <- .sft_sic_defective(cells, interpolate = interpolate)
  # The mean interaction contrast is taken over every observed response, not
  # just the correct ones, so that it estimates the same D2E the unconditional
  # predictions refer to.
  rts <- lapply(cells, function(z) z$rt[is.finite(z$rt)])
  mt <- mic.test(rts$HH, rts$HL, rts$LH, rts$LL, method = mictest)
  st <- if (sictest == "bf") {
    .sft_sic_dp_defective(cells, nbin = nbin, nsamp = nsamp, maxn = maxn,
                          seed = seed)
  } else {
    # The Houpt-Townsend KS route is not valid for subdistributions: its null
    # is a proper-CDF Kolmogorov law, whereas these curves carry mass < 1.
    # Group-level inference is available from sicPermTest().
    NULL
  }
  list(SIC = fit$SIC, SIC.correct = fit$SIC.correct, SIC.error = fit$SIC.error,
       MIC = unname(mt$statistic), MICtest = mt, SICtest = st, Dominance = dominance,
       Dvals = NULL, N = 1 / sum(1 / fit$n), times = fit$times, values = fit$values,
       accuracy = fit$accuracy, omissions = fit$omissions, asymptote = fit$asymptote,
       HH = fit$cdf$all$HH, HL = fit$cdf$all$HL, LH = fit$cdf$all$LH, LL = fit$cdf$all$LL,
       correct.cdf = fit$cdf$correct, error.cdf = fit$cdf$error,
       errors = "defective")
}


sicGroup <- function(inData, sictest = c("ks", "bf"), mictest = c("art", "anova"),
                     domtest = c("ks", "dp"), alpha.sic = .05, plotSIC = TRUE,
                     errors = c("discard", "defective"), min.trials = 10L,
                     perm = TRUE, nperm = 1001L, stat = c("CvM", "KS", "AD", "Aly"),
                     times = NULL, n.times = 1000L, seed = NULL, ...) {
  inData <- .sft_normalize_columns(inData)
  sictest <- match.arg(sictest); domtest <- match.arg(domtest); mictest <- match.arg(mictest)
  stat <- match.arg(stat)
  errors <- .sft_errors_method(errors)
  req <- c("Subject", "Condition", "RT", "Correct", "Channel1", "Channel2")
  if (!all(req %in% names(inData))) stop("inData is missing required SIC columns.")
  times <- if (is.null(times)) .sft_time_grid(inData$RT, n.times) else
    sort(unique(as.numeric(times)[is.finite(times)]))
  if (!length(times)) stop("times contains no finite values.")
  subjects <- .sft_sic_subject_cells(inData, errors = errors,
                                     min.trials = min.trials)
  records <- list(); curves <- list(); fits <- list()
  correct_curves <- list(); error_curves <- list()
  for (n in seq_along(subjects)) {
    s <- subjects[[n]]
    cell <- lapply(s$cells, function(z) z$rt)
    cr <- lapply(s$cells, function(z) z$correct)
    fit <- sic(cell$HH, cell$HL, cell$LH, cell$LL, domtest = domtest,
               sictest = sictest, mictest = mictest, CR = cr, errors = errors,
               seed = seed, ...)
    fits[[n]] <- fit; curves[[n]] <- fit$SIC(times)
    if (errors == "defective") {
      correct_curves[[n]] <- fit$SIC.correct(times)
      error_curves[[n]] <- fit$SIC.error(times)
    }
    records[[n]] <- data.frame(Subject = s$subject, Condition = s$condition,
                               stringsAsFactors = FALSE)
  }
  overview <- if (!length(records)) data.frame() else do.call(rbind, records)
  if (nrow(overview)) {
    if (errors == "defective") {
      # No Model column: architecture classification is not attempted on the
      # unconditional contrasts. What is reported instead is the accuracy each
      # cell reached and where each contrast settles as t grows.
      overview$Accuracy.HH <- vapply(fits, function(f) unname(f$accuracy[["HH"]]), numeric(1))
      overview$Accuracy.LL <- vapply(fits, function(f) unname(f$accuracy[["LL"]]), numeric(1))
      overview$SIC.correct.asymptote <- vapply(fits, function(f) unname(f$asymptote[["correct"]]), numeric(1))
      overview$SIC.error.asymptote <- vapply(fits, function(f) unname(f$asymptote[["error"]]), numeric(1))
      overview$MIC <- ifelse(vapply(fits, function(f) f$MICtest$p.value < alpha.sic, logical(1)),
                             "Significant", "Nonsignificant")
    } else if (sictest == "bf") {
      model <- vapply(fits, function(f) names(f$SICtest$statistic)[which.max(f$SICtest$statistic)], character(1))
      overview$Model <- model
    } else {
      positive <- vapply(fits, function(f) f$SICtest$positive$p.value < alpha.sic, logical(1))
      negative <- vapply(fits, function(f) f$SICtest$negative$p.value < alpha.sic, logical(1))
      overview$Positive.SIC <- ifelse(positive, "Significant", "Nonsignificant")
      overview$Negative.SIC <- ifelse(negative, "Significant", "Nonsignificant")
      overview$MIC <- ifelse(vapply(fits, function(f) f$MICtest$p.value < alpha.sic, logical(1)), "Significant", "Nonsignificant")
    }
  }
  empty <- matrix(numeric(), 0, length(times))
  out <- list(overview = overview,
              sic.fn = if (length(curves)) do.call(rbind, curves) else empty,
              sic = fits, times = times, errors = errors)
  if (errors == "defective") {
    out$sic.correct.fn <- if (length(correct_curves)) do.call(rbind, correct_curves) else empty
    out$sic.error.fn <- if (length(error_curves)) do.call(rbind, error_curves) else empty
    # The interaction contrast itself is tested at the group level: the
    # permutation scheme swaps whole participants between the two mixtures, so
    # it needs more than one participant to say anything.
    if (isTRUE(perm)) {
      conds <- unique(overview$Condition)
      out$perm <- stats::setNames(lapply(conds, function(cond) {
        if (sum(overview$Condition == cond) < 2L) return(NULL)
        sicPermTest(inData, stat = stat, nperm = nperm, alpha = alpha.sic,
                    min.trials = min.trials, Condition = cond, seed = seed)
      }), conds)
    }
  }
  out
}


sicGroupBF <- function(inData, domtest = "ks", plotSIC = TRUE, ...) {
  sicGroup(inData, sictest = "bf", domtest = domtest, plotSIC = plotSIC, ...)
}
