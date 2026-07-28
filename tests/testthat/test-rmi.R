rmi_data <- function(n_sub = 12, n = 150, speedup = 4, condition = NULL,
                     tag = "", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- do.call(rbind, lapply(seq_len(n_sub), function(s) {
    # A large, participant-specific base time: this is what the classical
    # millisecond-scale group statistic is vulnerable to.
    base <- stats::runif(1, 200, 900)
    ab <- base + stats::rexp(n, 1 / 150) / speedup
    a <- base + stats::rexp(n, 1 / 150)
    b <- base + stats::rexp(n, 1 / 150)
    data.frame(Subject = paste0("s", s, tag), RT = c(ab, a, b), Correct = TRUE,
               Channel1 = rep(c(2, 2, 0), each = n),
               Channel2 = rep(c(2, 0, 2), each = n),
               stringsAsFactors = FALSE)
  }))
  if (!is.null(condition)) out$Condition <- condition
  out
}


# A genuine race: AB is the minimum of two independent single-target draws, so
# the Miller bound holds by construction.
rmi_race_data <- function(n_sub = 12, n = 150, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  do.call(rbind, lapply(seq_len(n_sub), function(s) {
    base <- stats::runif(1, 200, 900)
    ab <- base + pmin(stats::rexp(n, 1 / 150), stats::rexp(n, 1 / 150))
    a <- base + stats::rexp(n, 1 / 150)
    b <- base + stats::rexp(n, 1 / 150)
    data.frame(Subject = paste0("s", s), RT = c(ab, a, b), Correct = TRUE,
               Channel1 = rep(c(2, 2, 0), each = n),
               Channel2 = rep(c(2, 0, 2), each = n),
               stringsAsFactors = FALSE)
  }))
}


test_that("the probability-scale estimand is free of participant base time", {
  set.seed(2)
  d <- rmi_data(n_sub = 6, n = 120)
  shifted <- d
  # Shift every one of a participant's three conditions by the same amount: the
  # inequality is unchanged, so the estimand must be too.
  offsets <- stats::setNames(stats::runif(6, -100, 400), unique(d$Subject))
  shifted$RT <- shifted$RT + offsets[shifted$Subject]

  grouped <- sftplus:::.sft_condition_input(d, NULL, "OR", NULL, NULL)
  grouped_shifted <- sftplus:::.sft_condition_input(shifted, NULL, "OR", NULL, NULL)
  probs <- seq(.1, .4, by = .1)
  a <- sftplus:::.sft_rmi_score_data(grouped$conditions[[1L]], probs, "miller")
  b <- sftplus:::.sft_rmi_score_data(grouped_shifted$conditions[[1L]], probs,
                                     "miller")
  expect_equal(a$Y, b$Y, tolerance = 1e-10)

  # A strictly increasing reparameterisation of time applied to all three
  # conditions leaves it invariant as well.
  warped <- d
  warped$RT <- warped$RT^1.3
  grouped_warped <- sftplus:::.sft_condition_input(warped, NULL, "OR", NULL, NULL)
  w <- sftplus:::.sft_rmi_score_data(grouped_warped$conditions[[1L]], probs,
                                     "miller")
  expect_equal(a$Y, w$Y, tolerance = 1e-10)
})


test_that("rmiGroup.bayes separates a coactive violation from a true race", {
  probs <- seq(.1, .4, by = .1)
  violating <- rmiGroup.bayes(rmi_data(speedup = 5, seed = 3), probs = probs,
                              ndraws = 3000, burnin = 500, chains = 2, seed = 7)
  racing <- rmiGroup.bayes(rmi_race_data(seed = 4), probs = probs,
                           ndraws = 3000, burnin = 500, chains = 2, seed = 7)

  expect_equal(violating$global$p_any_violation, 1, tolerance = .01)
  expect_true(all(violating$percentile$mean < 0))
  expect_lt(racing$global$p_any_violation, .8)
  expect_gt(min(racing$percentile$mean), -.05)

  # mu is on the probability scale, so a violation is bounded by 1.
  expect_true(all(abs(violating$percentile$mean) < 1))
  expect_lt(max(violating$diagnostics$rhat), 1.05)
})


test_that("the population curve is indexed by percentile, not by time", {
  fit <- rmiGroup.bayes(rmi_data(speedup = 5, seed = 3),
                        probs = seq(.1, .4, by = .1), ndraws = 1500,
                        burnin = 300, chains = 2, seed = 7)
  expect_equal(fit$percentile$prob, seq(.1, .4, by = .1))
  expect_true(all(grepl("^mu\\[p", fit$mu_names %||% fit$population_summary$parameter[1:4])))
  expect_true(is.character(fit$note) && nzchar(fit$note))
})


test_that("vacuous Miller levels are not fitted", {
  # Above roughly the median of the single-target distributions F_A + F_B >= 1
  # and the inequality is satisfied for free.
  fit <- rmiGroup.bayes(rmi_data(speedup = 5, seed = 3),
                        probs = seq(.1, .9, by = .1), ndraws = 1500,
                        burnin = 300, chains = 2, seed = 7)
  expect_lt(max(fit$percentile$prob), .9)
  expect_gt(nrow(fit$percentile), 2L)
})


test_that("the hierarchy refuses inputs it cannot identify", {
  one <- rmi_data(n_sub = 1, n = 120, seed = 5)
  expect_error(suppressWarnings(
    rmiGroup.bayes(one, probs = seq(.1, .4, by = .1), ndraws = 500,
                   burnin = 100, chains = 2)),
    "at least two participants|fewer than two")
  # A percentile level informed by one participant cannot separate mu[k] from
  # that participant's own delta, so it is dropped with a warning first.
  expect_warning(try(rmiGroup.bayes(one, probs = seq(.1, .4, by = .1),
                                    ndraws = 500, burnin = 100, chains = 2),
                     silent = TRUE),
                 "fewer than two participants")
  expect_error(sftplus:::.sft_rmi_probs(c(.2, .5)), "at least three")
  expect_error(rmiGroup.bayes(rmi_data(n_sub = 4, seed = 5), bound = "both"),
               "single bound")
})


test_that("multiple conditions get one curve each and per-percentile contrasts", {
  d <- rbind(rmi_data(n_sub = 8, speedup = 5, condition = "coactive", seed = 6),
             rmi_data(n_sub = 8, speedup = 1.05, condition = "flat",
                      tag = "b", seed = 7))
  fit <- rmiGroup.bayes(d, probs = seq(.1, .4, by = .1), ndraws = 2000,
                        burnin = 400, chains = 2, seed = 9)
  expect_setequal(fit$conditions, c("coactive", "flat"))
  expect_setequal(unique(fit$percentile$condition), c("coactive", "flat"))
  expect_equal(nrow(fit$global), 2L)
  expect_gt(fit$global$p_any_violation[fit$global$condition == "coactive"], .95)
  expect_lt(fit$global$p_any_violation[fit$global$condition == "flat"], .5)

  # Contrasts compare the same percentile across conditions, never different
  # percentiles of different conditions.
  expect_true(all(fit$percentile_contrasts$contrast == "coactive - flat"))
  expect_setequal(fit$percentile_contrasts$prob, seq(.1, .4, by = .1))
  expect_true(all(fit$percentile_contrasts$posterior_negative > .95))
})


test_that("the bootstrap covariance agrees with the analytic one", {
  d <- rmi_data(n_sub = 5, n = 200, speedup = 3, seed = 8)
  grouped <- sftplus:::.sft_condition_input(d, NULL, "OR", NULL, NULL)
  probs <- seq(.1, .4, by = .1)
  analytic <- sftplus:::.sft_rmi_score_data(grouped$conditions[[1L]], probs,
                                            "miller")
  set.seed(1)
  boot <- sftplus:::.sft_rmi_score_data(grouped$conditions[[1L]], probs,
                                        "miller", var_method = "bootstrap",
                                        n_boot = 400)
  expect_equal(analytic$Y, boot$Y, tolerance = 1e-10)
  se_a <- sqrt(vapply(analytic$S, function(s) mean(diag(s)), numeric(1)))
  se_b <- sqrt(vapply(boot$S, function(s) mean(diag(s)), numeric(1)))
  # The bootstrap also propagates uncertainty in the percentile grid, so it is
  # the larger of the two, but within a factor of about two.
  expect_true(all(se_b > se_a * .8))
  expect_true(all(se_b < se_a * 3))
})


test_that("build_rmi_bound_df gives the AB-against-bound figure", {
  fit <- rmi.test(rmi_data(speedup = 5, seed = 3), probs = seq(.1, .9, by = .1))
  d <- build_rmi_bound_df(fit)
  expect_setequal(unique(d$curve), c("AB", "A + B bound"))
  expect_true(all(c("bound", "condition", "prob", "curve", "RT") %in% names(d)))
  # A coactive AB lies to the left of the bound: at each percentile the AB
  # quantile is the faster one.
  wide <- stats::reshape(d[c("prob", "curve", "RT")], idvar = "prob",
                         timevar = "curve", direction = "wide")
  expect_true(all(wide$`RT.AB` < wide$`RT.A + B bound`))

  cdfs <- build_rmi_cdf_df(fit)
  expect_setequal(levels(cdfs$Condition), c("A", "B", "AB"))
  expect_equal(nrow(cdfs), 3L * length(unique(d$prob)))
})


test_that("RMI builders accept and retain a requested percentile window", {
  probs <- seq(.05, .95, by = .1)
  d <- rmi_data(speedup = 5, seed = 12)

  bound <- build_rmi_bound_df(d, probs = probs)
  expect_setequal(unique(bound$prob), probs)

  cdfs <- build_rmi_cdf_df(d, probs = probs)
  expect_setequal(unique(cdfs$prob), probs)

  violation <- build_rmi_violation_df(d, probs = probs)
  expect_setequal(unique(violation$prob), probs)

  fit <- rmi.test(d, probs = probs)
  expect_equal(build_rmi_bound_df(d, probs = probs),
               build_rmi_bound_df(fit), tolerance = 1e-10)

  d_block <- d
  n_subject <- length(unique(d_block$Subject))
  n_per_cell <- sum(d_block$Subject == d_block$Subject[[1L]]) / 3
  d_block$Block <- rep(rep(c(1, 2), each = n_per_cell / 2),
                       times = 3 * n_subject)
  fit_block <- rmi.test(d_block, probs = probs, by_block = TRUE)
  expect_equal(build_rmi_bound_df(d_block, probs = probs, by_block = TRUE),
               build_rmi_bound_df(fit_block), tolerance = 1e-10)

  selected <- build_rmi_bound_df(fit, probs = c(.15, .35, .55))
  expect_equal(sort(unique(selected$prob)), c(.15, .35, .55),
               tolerance = 1e-10)
  expect_error(build_rmi_bound_df(fit, probs = .12),
               "not present in the supplied result")
})


test_that("the builders handle both bounds and several conditions", {
  d <- rbind(rmi_data(n_sub = 6, speedup = 5, condition = "coactive", seed = 6),
             rmi_data(n_sub = 6, speedup = 1.05, condition = "flat",
                      tag = "b", seed = 7))
  fit <- rmi.test(d, bound = "both", probs = seq(.1, .9, by = .1))
  b <- build_rmi_bound_df(fit)
  expect_setequal(unique(b$bound), c("miller", "grice"))
  expect_setequal(unique(b$condition), c("coactive", "flat"))
  expect_setequal(unique(b$curve), c("AB", "A + B bound", "max(A, B) bound"))
  expect_equal(nrow(build_rmi_bound_df(fit, bound = "miller",
                                       condition = "flat")),
               2L * 9L)
  expect_error(build_rmi_bound_df(fit, bound = "nonesuch"), "matched")
})


test_that("build_rmi_violation_df picks the scale to match the fit", {
  d <- rmi_data(speedup = 5, seed = 3)
  freq <- build_rmi_violation_df(rmi.test(d, probs = seq(.1, .9, by = .1)))
  expect_true(all(freq$scale == "time"))
  expect_true(all(c("t", "p.value", "p.adjusted") %in% names(freq)))

  fit <- rmiGroup.bayes(d, probs = seq(.1, .4, by = .1), ndraws = 1500,
                        burnin = 300, chains = 2, seed = 7)
  bayes <- build_rmi_violation_df(fit)
  expect_true(all(bayes$scale == "probability"))
  expect_true(all(c("lower", "upper", "p_violation") %in% names(bayes)))
  expect_equal(bayes$mean, fit$percentile$mean)

  expect_error(build_rmi_violation_df(rmi.test(d), scale = "probability"),
               "rmiGroup.bayes")
})


test_that("the time-scale estimator reproduces Ulrich, Miller & Schroeter (2007) Table 2", {
  # The worked example of UMS Table 1, with the percentiles they publish in
  # Table 2. This is exact ground truth for the whole time-scale path: the
  # cumulative frequency polygon of their Equation 2, the tie correction of
  # their Appendix A (condition Cz has two tied pairs), the bound
  # B(t) = Gx(t) + Gy(t), and the inversion of both to percentiles.
  Cx <- c(244, 249, 257, 260, 264, 268, 271, 274, 277, 291)
  Cy <- c(245, 246, 248, 250, 251, 252, 253, 254, 255, 259, 263, 265, 279,
          282, 284, 319)
  Cz <- c(234, 238, 240, 240, 243, 243, 245, 251, 254, 256, 259, 270, 280)
  probs <- seq(.05, .95, by = .10)
  published_z <- c(234.6, 238.6, 240.4, 242.3, 244.1, 248.9, 253.9, 256.8,
                   265.1, 278.5)
  published_b <- c(244.0, 245.6, 247.3, 249.3, 250.9, 252.3, 253.6, 254.9,
                   257.8, 259.8)

  expect_equal(sftplus:::.sft_rmi_quantile(Cz, probs, qtype = 5), published_z,
               tolerance = 5e-4)
  curve <- sftplus:::.sft_rmi_bound_curve(list(A = sort(Cx), B = sort(Cy)),
                                          "miller", qtype = 5)
  expect_equal(sftplus:::.sft_rmi_invert(curve$grid, curve$value, probs),
               published_b, tolerance = 5e-4)

  # Whole-pipeline check through the per-participant time-scale scorer.
  scored <- sftplus:::.sft_rmi_time_subject(
    list(list(RT = list(AB = Cz, A = Cx, B = Cy), CR = NULL)), probs, "miller",
    min_n = 10L, qtype = 5)
  expect_equal(scored$q_AB, published_z, tolerance = 5e-4)
  expect_equal(scored$q_bound, published_b, tolerance = 5e-4)
  # Negative is a violation, and this example violates at the fast percentiles.
  # UMS print Table 2 to one decimal, so the difference of two printed values
  # carries up to 0.1 ms of rounding.
  expect_equal(scored$diff, scored$q_AB - scored$q_bound, tolerance = 1e-12)
  expect_lt(max(abs(scored$diff - (published_z - published_b))), .11)
  expect_true(all(scored$diff[probs <= .45] < 0))
})


test_that("ties follow UMS Appendix A rather than a plateau or the jump top", {
  # A value repeated n_i times must sit at the midpoint of its vertical jump.
  x <- c(244, 244, 249, 249, 257, 260, 264, 268, 271, 291)
  cdf <- sftplus:::.sft_rmi_cdf(x, qtype = 5)
  # 244 occupies steps 1 and 2, so the polygon passes through (0 + 2)/(2 * 10);
  # 249 occupies steps 3 and 4, giving (2 + 4)/(2 * 10).
  expect_equal(cdf(c(244, 249)), c(.10, .30), tolerance = 1e-12)
  # The old behaviour took the top of the jump, .15 and .35.
  expect_false(isTRUE(all.equal(cdf(244), .15)))
  # stats::quantile(type = 5) instead spreads a plateau across the tied value
  # and so disagrees with the polygon inverse wherever ties bind.
  probs <- c(.15, .25, .35)
  expect_false(isTRUE(all.equal(
    sftplus:::.sft_rmi_quantile(x, probs, qtype = 5),
    as.numeric(stats::quantile(x, probs, type = 5, names = FALSE)))))
  # Without ties the polygon inverse *is* stats::quantile(type = 5).
  y <- c(244, 249, 257, 260, 264, 268, 271, 274, 277, 291)
  expect_equal(sftplus:::.sft_rmi_quantile(y, probs, qtype = 5),
               as.numeric(stats::quantile(y, probs, type = 5, names = FALSE)),
               tolerance = 1e-12)
  # ... and the same holds for the other interpolating rules.
  expect_equal(sftplus:::.sft_rmi_quantile(y, probs, qtype = 8),
               as.numeric(stats::quantile(y, probs, type = 8, names = FALSE)),
               tolerance = 1e-12)

  # estimate.bounds(qtype = ) carries the same correction.
  positions <- sftplus:::.make_quantile_cdf(x, qtype = 5)
  expect_equal(positions(c(244, 249)), c(.10, .30), tolerance = 1e-12)
})


test_that("percentiles outside the UMS estimable range are NA, with a warning", {
  # UMS: no satisfactory estimate below p = 1/(2n) or above p = (2n - 1)/(2n).
  x <- seq(200, 400, length.out = 10)
  expect_true(is.na(sftplus:::.sft_rmi_quantile(x, .04, qtype = 5)))
  expect_true(is.na(sftplus:::.sft_rmi_quantile(x, .96, qtype = 5)))
  expect_false(is.na(sftplus:::.sft_rmi_quantile(x, .05, qtype = 5)))
  # stats::quantile clamps to the sample minimum instead of declining to guess.
  expect_equal(as.numeric(stats::quantile(x, .01, type = 5, names = FALSE)),
               min(x))

  d <- rmi_data(n_sub = 4, n = 12, seed = 31)
  expect_warning(rmi.test(d, probs = c(.02, .10, .20, .30), min_n = 10L),
                 "estimable range")
})


test_that("the defaults follow the Kiesel et al. (2007) recommendations", {
  expect_equal(formals(rmi.test)$min_n, 20L)
  expect_equal(formals(rmiGroup.bayes)$min_n, 20L)
  expect_equal(eval(formals(rmi.test)$probs), seq(.10, .25, by = .05))
  expect_equal(eval(formals(rmiGroup.bayes)$probs), seq(.10, .25, by = .05))
  # The restricted range is the primary Type I control, so it must not reach
  # past the median, where the Miller bound is vacuous.
  expect_true(all(eval(formals(rmi.test)$probs) <= .5))
  expect_equal(formals(rmi.test)$qtype, 5)
  expect_equal(formals(rmiGroup.bayes)$qtype, 5)

  # min_n is enforced: a participant with a cell below it is dropped.
  d <- rmi_data(n_sub = 4, n = 150, seed = 32)
  s1_a <- which(d$Subject == "s1" & d$Channel1 == 2 & d$Channel2 == 0)
  thin <- d[-s1_a[-seq_len(15L)], ]
  expect_warning(out <- rmi.test(thin, probs = seq(.10, .25, by = .05)),
                 "Dropping")
  expect_false("s1" %in% out$subject$subject)
  # ... and admitted once min_n is lowered to match.
  expect_true("s1" %in%
                rmi.test(thin, probs = seq(.10, .25, by = .05),
                         min_n = 15L)$subject$subject)
})


test_that("qtype = 5 reproduces Miller's (i - .5)/n plotting positions for both data and bounds", {
  d <- rmi_data(n_sub = 6, n = 100, seed = 11)
  probs <- seq(.1, .9, by = .1)
  miller <- rmi.test(d, probs = probs, qtype = 5)
  # Miller placed the ith of n ordered RTs at (i - .5)/n, which is what
  # quantile(type = 5) interpolates.
  one <- d[d$Subject == "s1" & d$Channel1 == 2 & d$Channel2 == 2, "RT"]
  expect_equal(miller$subject$q_AB[miller$subject$subject == "s1"],
               as.numeric(stats::quantile(sort(one), probs, names = FALSE,
                                          type = 5)),
               tolerance = 1e-10)
  # The choice moves the reported quantiles and the bound quantiles ...
  hazen <- rmi.test(d, probs = probs, qtype = 8)
  expect_false(isTRUE(all.equal(hazen$subject$q_AB, miller$subject$q_AB)))
  expect_false(isTRUE(all.equal(hazen$subject$q_bound, miller$subject$q_bound)))

  # estimate.bounds with qtype passes qtype down to single-channel CDFs
  x_a <- d[d$Subject == "s1" & d$Channel1 == 2 & d$Channel2 == 0, "RT"]
  x_b <- d[d$Subject == "s1" & d$Channel1 == 0 & d$Channel2 == 2, "RT"]
  b5 <- estimate.bounds(list(x_a, x_b), stopping.rule = "OR", qtype = 5)
  b_step <- estimate.bounds(list(x_a, x_b), stopping.rule = "OR")
  grid <- sort(unique(c(x_a, x_b)))
  expect_false(isTRUE(all.equal(b5$Upper.Bound(grid), b_step$Upper.Bound(grid))))
})


test_that("kernel CDF estimation is exact and symmetric across RMI cells", {
  x <- c(210, 240, 275, 330)
  h <- 18
  at <- c(200, 250, 300)
  fitted <- sftplus:::.make_kernel_cdf(x, bw = h)
  expected <- vapply(at, function(t) mean(stats::pnorm((log(t) - log(x)) / h)),
                     numeric(1))
  expect_equal(fitted(at), expected, tolerance = 1e-14)

  d <- rmi_data(n_sub = 5, n = 80, speedup = 3, seed = 41)
  probs <- seq(.10, .40, by = .10)
  k5 <- rmi.test(d, probs = probs, cdf_method = "kernel",
                 kernel_bw = h, qtype = 5)
  k8 <- rmi.test(d, probs = probs, cdf_method = "kernel",
                 kernel_bw = h, qtype = 8)
  # qtype has no role in the fully kernel-smoothed analysis.
  expect_equal(k5$subject[c("q_A", "q_B", "q_AB", "q_bound")],
               k8$subject[c("q_A", "q_B", "q_AB", "q_bound")],
               tolerance = 1e-12)

  one <- d[d$Subject == "s1", ]
  cells <- list(
    AB = one$RT[one$Channel1 == 2 & one$Channel2 == 2],
    A = one$RT[one$Channel1 == 2 & one$Channel2 == 0],
    B = one$RT[one$Channel1 == 0 & one$Channel2 == 2])
  manual <- sftplus:::.sft_rmi_kernel_estimates(cells, probs, "miller",
                                                kernel_bw = h)
  observed <- k5$subject[k5$subject$subject == "s1", ]
  expect_equal(observed$q_AB, manual$q_AB, tolerance = 1e-8)
  expect_equal(observed$q_bound, manual$q_bound, tolerance = 1e-8)

  polygon <- rmi.test(d, probs = probs, cdf_method = "polygon")
  expect_false(isTRUE(all.equal(k5$subject$q_AB, polygon$subject$q_AB)))
  expect_equal(build_rmi_bound_df(d, probs = probs, cdf_method = "kernel",
                                  kernel_bw = h),
               build_rmi_bound_df(k5), tolerance = 1e-10)
})


test_that("estimate.bounds supports the same Gaussian kernel CDF", {
  a <- c(220, 245, 270, 310, 355)
  b <- c(210, 250, 285, 325, 370)
  h <- 20
  fit <- estimate.bounds(list(a, b), stopping.rule = "OR",
                         cdf_method = "kernel", kernel_bw = h)
  at <- seq(225, 350, by = 25)
  expected_upper <- vapply(at, function(t) {
    mean(stats::pnorm((log(t) - log(a)) / h)) + mean(stats::pnorm((log(t) - log(b)) / h))
  }, numeric(1))
  expected_lower <- vapply(at, function(t) {
    max(mean(stats::pnorm((log(t) - log(a)) / h)),
        mean(stats::pnorm((log(t) - log(b)) / h)))
  }, numeric(1))
  expect_equal(fit$Upper.Bound(at), expected_upper, tolerance = 1e-4)
  expect_equal(fit$Lower.Bound(at), expected_lower, tolerance = 1e-4)

  # Existing defaults retain the step-ECDF and qtype-triggered polygon paths.
  expect_equal(formals(estimate.bounds)$cdf_method, NULL)
  expect_error(estimate.bounds(list(a, b), cdf_method = "kernel",
                               kernel_bw = 0), "positive finite bandwidth")
})


test_that("the Grice bound is oriented so that negative is still a violation", {
  # AB slower than the faster single target violates the lower bound.
  set.seed(13)
  n <- 150
  d <- do.call(rbind, lapply(1:8, function(s) {
    base <- stats::runif(1, 200, 900)
    data.frame(Subject = paste0("s", s),
               RT = c(base + stats::rexp(n, 1 / 150) * 2.5,
                      base + stats::rexp(n, 1 / 150),
                      base + stats::rexp(n, 1 / 150)),
               Correct = TRUE, Channel1 = rep(c(2, 2, 0), each = n),
               Channel2 = rep(c(2, 0, 2), each = n), stringsAsFactors = FALSE)
  }))
  # max(F_A, F_B) is not differentiable where the channel CDFs cross, and the
  # analytic covariance says so.
  expect_warning(
    fit <- rmiGroup.bayes(d, bound = "grice", probs = seq(.2, .8, by = .2),
                          ndraws = 2000, burnin = 400, chains = 2, seed = 7),
    "not differentiable")
  expect_equal(fit$bound, "grice")
  expect_true(all(fit$percentile$mean < 0))
  expect_gt(fit$global$p_any_violation, .95)
})


test_that("errors = 'defective' counts errors in the denominator", {
  # AB is faster *and* more accurate. Under "discard" the extra correct AB
  # responses are renormalised away, which flatters the redundant-target CDF;
  # under "defective" they are counted against the same denominator.
  set.seed(17)
  n <- 200
  d <- do.call(rbind, lapply(1:10, function(s) {
    base <- stats::runif(1, 200, 900)
    ab <- base + stats::rexp(n, 1 / 150) / 3
    a <- base + stats::rexp(n, 1 / 150)
    b <- base + stats::rexp(n, 1 / 150)
    data.frame(Subject = paste0("s", s), RT = c(ab, a, b),
               Correct = c(stats::runif(n) < .95, stats::runif(2 * n) < .75),
               Channel1 = rep(c(2, 2, 0), each = n),
               Channel2 = rep(c(2, 0, 2), each = n), stringsAsFactors = FALSE)
  }))
  grouped <- sftplus:::.sft_condition_input(d, NULL, "OR", NULL, NULL)
  probs <- seq(.1, .4, by = .1)

  cells_keep <- sftplus:::.sft_rmi_cells(grouped$conditions[[1L]]$RT[[1L]],
                                         grouped$conditions[[1L]]$CR[[1L]],
                                         errors = "defective")
  # Every trial is retained, errors as Inf.
  expect_equal(unname(lengths(cells_keep)), rep(n, 3L))
  expect_true(all(vapply(cells_keep, function(v) any(!is.finite(v)), logical(1))))
  expect_equal(sftplus:::.sft_rmi_ecdf_at(cells_keep$A, Inf), 1)

  discard <- sftplus:::.sft_rmi_score_data(grouped$conditions[[1L]], probs,
                                           "miller")
  defective <- sftplus:::.sft_rmi_score_data(grouped$conditions[[1L]], probs,
                                             "miller", errors = "defective")
  # The single-target CDFs shrink by the error rate, so the Miller bound is
  # lower and the estimated margin is smaller under "defective".
  expect_true(all(defective$Y < discard$Y))
  expect_true(all(is.finite(defective$Y)))

  # The evaluation grid is read off the observed responses either way.
  expect_equal(discard$score$time, defective$score$time, tolerance = 1e-10)

  fit <- rmiGroup.bayes(d, probs = probs, errors = "defective", ndraws = 1500,
                        burnin = 300, chains = 2, seed = 7)
  expect_equal(fit$errors, "defective")
  expect_true(all(abs(fit$percentile$mean) < 1))
  expect_lt(max(fit$diagnostics$rhat), 1.05)
})


# Blocked data with an optional across-block drift in overall speed. The drift
# is common to all three cells within a block, so it changes nothing about the
# inequality within a block, but it spreads the session-pooled distributions.
rmi_block_data <- function(n_sub = 8, n = 240, n_block = 4, speedup = 4,
                           drift = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  do.call(rbind, lapply(seq_len(n_sub), function(s) {
    base <- stats::runif(1, 300, 700)
    block <- rep(seq_len(n_block), each = n / n_block)
    shift <- base + (block - 1) * drift
    data.frame(Subject = paste0("s", s), Block = rep(block, 3),
               RT = c(shift + stats::rexp(n, 1 / 150) / speedup,
                      shift + stats::rexp(n, 1 / 150),
                      shift + stats::rexp(n, 1 / 150)),
               Correct = TRUE, Channel1 = rep(c(2, 2, 0), each = n),
               Channel2 = rep(c(2, 0, 2), each = n),
               stringsAsFactors = FALSE)
  }))
}


test_that("by_block averages within-block quantiles, one block reproduces the default", {
  probs <- seq(.1, .4, by = .1)
  one <- rmi_block_data(n_sub = 6, n = 120, n_block = 1, seed = 21)
  # A single block is the whole session, so the block-averaging reduction is
  # the identity and the default path must be reproduced exactly.
  expect_equal(rmi.test(one, probs = probs, by_block = TRUE)$subject,
               rmi.test(one, probs = probs)$subject, tolerance = 1e-12)

  d <- rmi_block_data(n_sub = 6, n = 240, n_block = 4, seed = 22)
  blocked <- rmi.test(d, probs = probs, by_block = TRUE)
  expect_true(all(blocked$subject$n_blocks == 4L))
  # Each reported quantile is the mean of the four within-block quantiles.
  by_hand <- vapply(split(d[d$Subject == "s1" & d$Channel1 == 2 &
                              d$Channel2 == 2, ], d$Block[d$Subject == "s1" &
                              d$Channel1 == 2 & d$Channel2 == 2],
                          drop = TRUE),
                    function(part) stats::quantile(sort(part$RT), .2,
                                                   names = FALSE, type = 5),
                    numeric(1))
  expect_equal(blocked$subject$q_AB[blocked$subject$subject == "s1" &
                                      blocked$subject$prob == .2],
               mean(by_hand), tolerance = 1e-10)
})


test_that("by_block recovers a violation that a session-wide drift hides", {
  probs <- seq(.1, .4, by = .1)
  # The same coactive violation in every block, plus a large practice/fatigue
  # drift between blocks. Pooling the session spreads all three cells by the
  # drift, which flattens the CDFs and moves the bound; estimating within block
  # is Miller's remedy.
  d <- rmi_block_data(n_sub = 8, drift = 120, speedup = 4, seed = 23)
  pooled <- rmiGroup.bayes(d, probs = probs, ndraws = 2000, burnin = 400,
                           chains = 2, seed = 7)
  blocked <- rmiGroup.bayes(d, probs = probs, ndraws = 2000, burnin = 400,
                            chains = 2, seed = 7, by_block = TRUE)
  expect_true(all(blocked$percentile$mean < pooled$percentile$mean))
  expect_gt(blocked$global$p_any_violation, .95)
  expect_lt(pooled$global$p_any_violation, blocked$global$p_any_violation)
  expect_true(blocked$by_block)
  expect_true(all(blocked$score$n_blocks == 4L))

  # Without a drift the two agree on the sign and are close on the size.
  flat <- rmi_block_data(n_sub = 8, drift = 0, speedup = 4, seed = 24)
  a <- rmiGroup.bayes(flat, probs = probs, ndraws = 2000, burnin = 400,
                      chains = 2, seed = 7)
  b <- rmiGroup.bayes(flat, probs = probs, ndraws = 2000, burnin = 400,
                      chains = 2, seed = 7, by_block = TRUE)
  expect_equal(a$percentile$mean, b$percentile$mean, tolerance = .05)
})


test_that("blocks are pooled independently and short blocks are dropped", {
  probs <- seq(.1, .4, by = .1)
  d <- rmi_block_data(n_sub = 6, n = 240, n_block = 4, seed = 25)
  input <- sftplus:::.sft_rmi_input(d, NULL, NULL, NULL, TRUE)
  blocks <- sftplus:::.sft_rmi_subject_blocks(input$conditions[[1L]], 1L)
  expect_length(blocks, 4L)
  score <- sftplus:::.sft_rmi_subject_score(blocks, probs, "miller")
  singles <- lapply(blocks, function(b)
    sftplus:::.sft_rmi_subject_score(list(b), probs, "miller"))
  # The delta is the mean of the per-block deltas and, blocks being disjoint
  # trials, its variance is their summed variance over B^2.
  expect_equal(score$delta,
               rowMeans(vapply(singles, `[[`, numeric(length(probs)), "delta")),
               tolerance = 1e-12)
  expect_equal(diag(score$cov),
               rowSums(vapply(singles, function(s) diag(s$cov),
                              numeric(length(probs)))) / 16,
               tolerance = 1e-12)

  # A participant missing from a block, and a block too thin in one cell, lose
  # those blocks only.
  thinned <- d[!(d$Subject == "s2" & d$Block == 3L), ]
  thinned <- thinned[!(thinned$Subject == "s3" & thinned$Block == 4L &
                         thinned$Channel2 == 0), ]
  out <- rmi.test(thinned, probs = probs, by_block = TRUE)
  n_blocks <- out$subject$n_blocks[!duplicated(out$subject$subject)]
  names(n_blocks) <- unique(out$subject$subject)
  expect_equal(unname(n_blocks[c("s2", "s3")]), c(3L, 3L))
  expect_true(all(n_blocks[c("s1", "s4", "s5", "s6")] == 4L))
})


test_that("by_block needs a Block column and trial-level input", {
  d <- rmi_block_data(n_sub = 4, n = 120, n_block = 2, seed = 26)
  expect_error(rmi.test(d[setdiff(names(d), "Block")], by_block = TRUE),
               "Block column")
  expect_error(rmi.test(list(AB = 1:30, A = 1:30, B = 1:30), by_block = TRUE),
               "trial-level data frame")
  missing_block <- d
  missing_block$Block[1L] <- NA
  expect_error(rmi.test(missing_block, by_block = TRUE), "must not contain")
  # The Block column is found under the usual aliases.
  alias <- d
  names(alias)[names(alias) == "Block"] <- "block"
  expect_equal(rmi.test(alias, probs = seq(.1, .4, by = .1),
                        by_block = TRUE)$test,
               rmi.test(d, probs = seq(.1, .4, by = .1), by_block = TRUE)$test)
})
