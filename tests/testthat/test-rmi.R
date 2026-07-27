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


test_that("qtype = 5 reproduces Miller's (i - .5)/n plotting positions", {
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
  # The choice moves the reported quantiles ...
  hazen <- rmi.test(d, probs = probs, qtype = 8)
  expect_false(isTRUE(all.equal(hazen$subject$q_AB, miller$subject$q_AB)))
  # ... but it is a display convention only. The probability-scale estimand is
  # scored on its own grid and does not take qtype at all.
  expect_false("qtype" %in% names(formals(sftplus:::.sft_rmi_score_data)))
  expect_false("qtype" %in% names(formals(sftplus:::.sft_rmi_grid)))
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
