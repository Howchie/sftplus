# Unconditional (defective) interaction contrasts: Gondan & Legner's extension
# of Townsend & Nozawa to tasks with non-negligible error rates.

# Poisson race channel (their Example 1): correct and error accumulators race
# and salience raises only the correct rate, so assumptions 8A-8D all hold.
sicdef_channel <- function(n, kappa, K = 5, L = 5, lam = 6) {
  uc <- stats::rgamma(n, shape = K, rate = kappa)
  ue <- stats::rgamma(n, shape = L, rate = lam)
  list(d = pmin(uc, ue), r = uc <= ue)
}

sicdef_combine <- function(A, B, kind) {
  switch(kind,
    race = {
      first <- A$d <= B$d
      list(rt = ifelse(first, A$d, B$d), cr = ifelse(first, A$r, B$r))
    },
    or_pst = {
      both <- A$r & B$r
      list(rt = ifelse(both, pmin(A$d, B$d),
                ifelse(A$r, A$d, ifelse(B$r, B$d, pmax(A$d, B$d)))),
           cr = A$r | B$r)
    },
    and_pex = list(rt = pmax(A$d, B$d), cr = A$r & B$r))
}

sicdef_cells <- function(kind, n, kHi = 9, kLo = 5) {
  mk <- function(k1, k2) sicdef_combine(sicdef_channel(n, k1), sicdef_channel(n, k2), kind)
  list(HH = mk(kHi, kHi), HL = mk(kHi, kLo), LH = mk(kLo, kHi), LL = mk(kLo, kLo))
}

sicdef_df <- function(kind, nsub, n = 250L, seed = 1L) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(nsub), function(s) {
    d <- sicdef_cells(kind, n, kHi = stats::runif(1, 8, 10), kLo = stats::runif(1, 4.5, 5.5))
    do.call(rbind, Map(function(nm, x) data.frame(
      Subject = paste0("P", s), Condition = "C1", RT = x$rt, Correct = x$cr,
      Channel1 = if (substr(nm, 1L, 1L) == "H") 2 else 1,
      Channel2 = if (substr(nm, 2L, 2L) == "H") 2 else 1,
      stringsAsFactors = FALSE), names(d), d))
  }))
}


test_that("defective subdistributions decompose the agnostic contrast exactly", {
  set.seed(4)
  d <- sicdef_cells("or_pst", 400L)
  fit <- sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt,
             CR = lapply(d, `[[`, "cr"), errors = "defective")
  expect_equal(fit$values$correct + fit$values$error, fit$values$all)
  # F+ asymptotes at the cell accuracy, not at 1.
  last <- max(fit$times)
  expect_equal(unname(vapply(fit$correct.cdf, function(f) f(last), numeric(1))),
               unname(fit$accuracy))
  expect_true(all(vapply(fit$correct.cdf, function(f) f(last), numeric(1)) < 1))
  # The agnostic contrast vanishes at large t, so the two asymptotes cancel.
  expect_equal(unname(fit$asymptote[["correct"]] + fit$asymptote[["error"]]), 0)
  expect_null(fit$SICtest)
  expect_false("Model" %in% names(fit))
})


test_that("non-finite response times count as omissions in the denominator", {
  rt <- list(HH = c(1:9, NA), HL = 1:10, LH = 1:10, LL = 1:10)
  cr <- lapply(rt, function(x) rep(TRUE, length(x)))
  fit <- sic(rt$HH, rt$HL, rt$LH, rt$LL, CR = cr, errors = "defective")
  expect_equal(unname(fit$omissions[["HH"]]), 0.1)
  expect_equal(unname(fit$accuracy[["HH"]]), 0.9)
  # 9 correct responses out of a cell of 10 means F+ tops out at 0.9.
  expect_equal(fit$correct.cdf$HH(100), 0.9)
  expect_equal(fit$correct.cdf$HL(100), 1)
})


test_that("the three contrasts carry the signs the theory predicts", {
  set.seed(21)
  expectations <- list(
    race = c(correct = 1, error = -1, all = 1),
    or_pst = c(correct = 1, error = -1, all = 1),
    and_pex = c(correct = -1, error = 1, all = -1))
  for (kind in names(expectations)) {
    d <- sicdef_cells(kind, 20000L)
    fit <- sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt,
               CR = lapply(d, `[[`, "cr"), errors = "defective")
    # Judge on the interior of the support; the ECDF edges are pure noise.
    keep <- fit$times > stats::quantile(fit$times, .02) &
      fit$times < stats::quantile(fit$times, .98)
    for (nm in names(expectations[[kind]])) {
      v <- fit$values[[nm]][keep] * expectations[[kind]][[nm]]
      expect_gt(max(v), 0.02)
      expect_gt(min(v), -0.01)   # no meaningful excursion the other way
    }
  }
})


test_that("the dominance battery covers 8A, 8B and 8D and reads in the right direction", {
  set.seed(5)
  d <- sicdef_cells("or_pst", 3000L)
  dom <- siDominance(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt,
                     CR = lapply(d, `[[`, "cr"), errors = "defective")
  expect_equal(nrow(dom), 24L)
  expect_true(all(c("Assumption", "Test", "Predicted") %in% names(dom)))
  expect_equal(sort(unique(dom$Assumption)),
               c("8A (correct subdistribution)", "8B (survivor ordering)",
                 "8D (error subdistribution)"))
  # Salience makes responses faster, so the survivor rows are labelled S.hh <= S.hl.
  expect_true(all(grepl("<=", dom$Test[dom$Predicted & dom$Assumption == "8B (survivor ordering)"])))
  expect_true(all(grepl(">=", dom$Test[dom$Predicted & dom$Assumption == "8A (correct subdistribution)"])))
  expect_true(all(grepl("<=", dom$Test[dom$Predicted & dom$Assumption == "8D (error subdistribution)"])))
  # A model that satisfies 8A-8D confirms every prediction and no reversal.
  expect_true(all(dom$p.value[dom$Predicted] < .05))
  expect_true(all(dom$p.value[!dom$Predicted] > .05))
})


test_that("the Dirichlet route reports contrast signs rather than architectures", {
  set.seed(6)
  d <- sicdef_cells("or_pst", 1500L)
  bf <- sictestBayes(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt,
                     CR = lapply(d, `[[`, "cr"), errors = "defective",
                     nbin = 12L, nsamp = 1500L, maxn = 4500L, seed = 12L)
  expect_equal(names(bf$statistic), c("correct", "error", "all"))
  expect_equal(names(bf$statistic$correct),
               c("Zero", "Negative", "Positive", "Crosses"))
  best <- vapply(bf$statistic, function(b) names(b)[which.max(b)], character(1))
  expect_equal(unname(best[c("correct", "error")]), c("Positive", "Negative"))
  # None of the six architecture labels appear.
  expect_false(any(c("NegPos.MIC0", "NegPos.MICpos", "NegPos.MICneg") %in%
                     unlist(lapply(bf$statistic, names))))
})


test_that("sicPermTest tests all three contrasts across participants", {
  df <- sicdef_df("or_pst", nsub = 10L, seed = 31L)
  out <- sicPermTest(df, nperm = 99L, seed = 2L)
  expect_equal(out$summary$Contrast, c("correct", "error", "agnostic"))
  expect_true(all(out$summary$p.SICpos > 0 & out$summary$p.SICpos <= 1))
  expect_true(all(out$summary$p.SICneg > 0 & out$summary$p.SICneg <= 1))
  expect_s3_class(out$correct, "htest")
  expect_length(out$subjects, 10L)
  expect_identical(out$stat, "CvM")
  # A parallel self-terminating OR model has a positive correct-response SIC.
  expect_lt(out$correct$p.SICpos, .05)
  expect_gt(out$correct$p.SICneg, .05)
  expect_identical(out$correct$decision, "SIC positive")
  ks <- sicPermTest(df, contrast = "correct", stat = "KS", nperm = 99L, seed = 2L)
  expect_lt(ks$correct$p.SICpos, ks$correct$p.SICneg)
  expect_error(sicPermTest(df[df$Subject == "P1", , drop = FALSE]),
               "at least two participants")
})


test_that("the Aly statistic tests monotonicity, not sign", {
  df <- sicdef_df("or_pst", nsub = 10L, seed = 31L)
  aly <- sicPermTest(df, stat = "Aly", nperm = 99L, seed = 2L)
  # Different question, so different vocabulary.
  expect_equal(names(aly$summary)[2:3], c("p.SICdec", "p.SICinc"))
  expect_null(aly$correct$p.SICpos)
  expect_true(all(aly$summary$Decision %in%
                    c("n.s.", "SIC increasing", "SIC decreasing", "SIC non-monotone")))
  # A parallel self-terminating OR model has a monotone error contrast: the
  # error subdensity contrast has one sign, so SIC.error only ever decreases.
  expect_lt(aly$error$p.SICdec, .05)
  expect_gt(aly$error$p.SICinc, .05)
  # The correct contrast rises to a peak and then decays, so it is not monotone.
  expect_lt(aly$correct$p.SICinc, .05)
})


test_that("sicGroup gains the unconditional curves and drops classification", {
  df <- sicdef_df("or_pst", nsub = 6L, seed = 41L)
  g <- sicGroup(df, errors = "defective", nperm = 49L, seed = 3L)
  expect_identical(g$errors, "defective")
  expect_equal(nrow(g$overview), 6L)
  expect_false("Model" %in% names(g$overview))
  expect_false("Positive.SIC" %in% names(g$overview))
  expect_true(all(c("Accuracy.HH", "SIC.correct.asymptote", "SIC.error.asymptote")
                  %in% names(g$overview)))
  expect_equal(dim(g$sic.correct.fn), dim(g$sic.fn))
  expect_equal(dim(g$sic.error.fn), dim(g$sic.fn))
  expect_equal(g$sic.correct.fn + g$sic.error.fn, g$sic.fn)
  expect_named(g$perm, "C1")
  expect_equal(g$perm$C1$summary$Contrast, c("correct", "error", "agnostic"))
  # perm = FALSE skips the group test.
  expect_null(sicGroup(df, errors = "defective", perm = FALSE)$perm)
})


test_that("build_sic_defective_df returns one row per time and contrast", {
  set.seed(8)
  d <- sicdef_cells("or_pst", 300L)
  fit <- sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt,
             CR = lapply(d, `[[`, "cr"), errors = "defective")
  tt <- seq(stats::quantile(fit$times, .1), stats::quantile(fit$times, .9), length.out = 40)
  out <- build_sic_defective_df(tt, list(list(label = "S1",
    fn = list(correct = fit$correct.cdf, error = fit$error.cdf))),
    trim_cdf_tails = FALSE)
  expect_equal(nrow(out), 120L)
  expect_equal(sort(unique(out$Contrast)), c("Agnostic", "Correct", "Error"))
  expect_equal(out$SIC[out$Contrast == "Correct"] + out$SIC[out$Contrast == "Error"],
               out$SIC[out$Contrast == "Agnostic"])
})


test_that("errors = 'discard' is unchanged and 'defective' needs correctness", {
  set.seed(9)
  T1h <- stats::rexp(150, .2); T1l <- stats::rexp(150, .1)
  T2h <- stats::rexp(150, .21); T2l <- stats::rexp(150, .11)
  HH <- pmin(T1h, T2h); HL <- pmin(T1h, T2l)
  LH <- pmin(T1l, T2h); LL <- pmin(T1l, T2l)
  fit <- sic(HH, HL, LH, LL)
  tt <- seq(0.5, 30, by = .5)
  reference <- stats::ecdf(LH)(tt) + stats::ecdf(HL)(tt) -
    stats::ecdf(HH)(tt) - stats::ecdf(LL)(tt)
  expect_equal(fit$SIC(tt), reference)
  expect_identical(fit$errors, "discard")
  expect_equal(nrow(fit$Dominance), 8L)
  expect_false("Assumption" %in% names(fit$Dominance))
  expect_s3_class(fit$SICtest$positive, "htest")
  expect_equal(dim(fit$Dvals), c(2L, 2L))
  # CR is optional under "discard" but mandatory under "defective".
  expect_error(sic(HH, HL, LH, LL, errors = "defective"), "correctness indicators")
  cr <- list(rep(c(TRUE, FALSE), c(140, 10)), rep(TRUE, 150),
             rep(TRUE, 150), rep(TRUE, 150))
  filtered <- sic(HH, HL, LH, LL, CR = cr)
  expect_equal(filtered$HH(tt), stats::ecdf(HH[cr[[1]]])(tt))
  expect_error(sic(HH, HL, LH, LL, CR = cr[1:3], errors = "defective"),
               "list of four")
})


test_that("the errors argument accepts documented aliases", {
  set.seed(10)
  d <- sicdef_cells("race", 200L)
  cr <- lapply(d, `[[`, "cr")
  a <- sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt, CR = cr, errors = "defective")
  b <- sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt, CR = cr, errors = "unconditional")
  expect_equal(a$values, b$values)
  expect_identical(sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt, errors = "conditional")$errors,
                   "discard")
  expect_error(sic(d$HH$rt, d$HL$rt, d$LH$rt, d$LL$rt, errors = "nonsense"),
               "'discard' or 'defective'")
})


test_that("sicPermTest refuses to pool conditions and reports a valid p-value", {
  a <- sicdef_df("or_pst", nsub = 5L, seed = 51L)
  b <- sicdef_df("and_pex", nsub = 5L, seed = 52L)
  b$Condition <- "C2"; b$Subject <- paste0(b$Subject, "b")
  both <- rbind(a, b)
  expect_error(sicPermTest(both, nperm = 19L), "several conditions")
  out <- sicPermTest(both, contrast = "correct", nperm = 99L,
                     Condition = "C2", seed = 4L)
  expect_length(out$subjects, 5L)
  # The headline p.value is the Bonferroni combination of the two one-sided
  # tests, so it never dips below either of them.
  expect_equal(out$correct$p.value,
               min(1, 2 * min(out$correct$p.SICpos, out$correct$p.SICneg)))
  # An AND parallel exhaustive model gives a negative correct-response SIC.
  expect_lt(out$correct$p.SICneg, out$correct$p.SICpos)
  # sicGroup keeps the conditions apart for us.
  g <- sicGroup(both, errors = "defective", nperm = 19L, seed = 5L)
  expect_named(g$perm, c("C1", "C2"))
})
