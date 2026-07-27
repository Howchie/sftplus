# Resilience Bayes, multi-condition hierarchies, prior sensitivity,
# spike-and-slab, and the WAIC/LOO helpers.

ucip_or_subject <- function(n = 150, rate = 1.6) {
  a <- stats::rexp(n, rate); b <- stats::rexp(n, rate)
  list(pmin(a, b), stats::rexp(n, rate), stats::rexp(n, rate))
}

super_subject <- function(n = 150, rate = 1.6, ab_rate = 8) {
  list(stats::rexp(n, ab_rate), stats::rexp(n, rate), stats::rexp(n, rate))
}


test_that("the resilience test uses the score statistic, not a terminal-time contrast", {
  # Resilience shares the OR capacity functional form, so its test is the
  # Houpt-Townsend weighted-logrank score statistic (Houpt & Little, 2016).
  # A terminal-time Nelson-Aalen contrast cannot work: for complete data the
  # cumulative hazard at a sample's own maximum is exactly sum(1/k) whatever
  # the distribution, so all three series coincide there.
  set.seed(9)
  n <- 120
  rt <- list(stats::rlnorm(n, log(.33), .25), stats::rlnorm(n, log(.5), .25),
             stats::rlnorm(n, log(.5), .25))
  h <- lapply(1:3, function(i) estimateNAH(rt[[i]], rep(TRUE, n)))
  tmax <- max(unlist(rt))
  harmonic <- sum(1 / seq_len(n))
  expect_equal(h[[1]]$H(tmax), harmonic, tolerance = 1e-8)
  expect_equal(h[[2]]$H(tmax), harmonic, tolerance = 1e-8)

  # The score statistic separates the three regimes the terminal contrast cannot.
  set.seed(11)
  z <- vapply(list(super_subject(), ucip_or_subject(),
                   list(stats::rexp(150, 1.8), stats::rexp(150, 1.6),
                        stats::rexp(150, 1.6))),
              function(d) unname(resilience.test(d)$statistic), numeric(1))
  expect_gt(z[[1]], 4)
  expect_lt(abs(z[[2]]), 2)
  expect_lt(z[[3]], -3)
})


test_that("resilience.bayes recovers super-capacity and the UCIP-OR null", {
  set.seed(4)
  sup <- lapply(1:6, function(i) super_subject())
  names(sup) <- paste0("s", 1:6)
  g <- resilienceGroup.bayes(sup, ndraws = 3000, burnin = 400, chains = 2, seed = 2)
  expect_gt(g$population_summary$lower[1], 0)
  expect_gt(g$bayes_factor$point_null$BF10, 50)

  set.seed(6)
  null <- lapply(1:6, function(i) ucip_or_subject())
  names(null) <- paste0("s", 1:6)
  g0 <- resilienceGroup.bayes(null, ndraws = 3000, burnin = 400, chains = 2, seed = 2)
  expect_lt(g0$bayes_factor$point_null$BF10, 1)
  expect_lt(abs(g0$population_summary$mean[1]), 0.25)
})


test_that("resilience.bayes is the score model on the resilience conditions", {
  set.seed(12)
  rt <- super_subject()
  fit <- resilience.bayes(rt, ndraws = 2000, chains = 2, seed = 3)
  # theta_hat * sqrt(V) is exactly the resilience.test() z statistic.
  expect_equal(fit$score$theta_hat * sqrt(fit$score$V),
               unname(resilience.test(rt)$statistic), tolerance = 1e-8)
  expect_equal(fit$score$Cz, unname(resilience.test(rt)$statistic),
               tolerance = 1e-8)
  # No evaluation horizon survives in the interface.
  expect_null(fit$at)
  expect_error(resilience.bayes(rt, at = 0.4), "unused argument")
})


test_that("the condition factor reduces to the single-mean model", {
  # With one condition the generalised sampler must be the old one exactly.
  set.seed(3)
  hat <- stats::rnorm(9, .4, .2); prec <- stats::runif(9, 15, 45)
  set.seed(77)
  one <- sftplus:::.sft_hierarchical_normal_chain(hat, prec, 400L, 100L, 1L,
                                                  0, .5, 2, .25)
  set.seed(77)
  grp <- sftplus:::.sft_hierarchical_normal_chain(hat, prec, 400L, 100L, 1L,
                                                  0, .5, 2, .25,
                                                  group = rep(1L, 9))
  expect_equal(as.numeric(one$mu), as.numeric(grp$mu))
  expect_equal(one$tau, grp$tau)
  expect_equal(one$theta, grp$theta)
})


test_that("multi-condition fits recover a known condition difference", {
  set.seed(15)
  mk <- function(i, ab_rate) list(stats::rexp(120, ab_rate), stats::rexp(120, 1.6),
                                  stats::rexp(120, 1.6))
  easy <- lapply(1:6, function(i) mk(i, 5.2))
  hard <- lapply(1:6, function(i) mk(i, 3.4))
  hat <- c(vapply(easy, function(rt) {
    s <- sftplus:::.sft_ucip_score(rt); s$numer / s$variance }, numeric(1)),
    vapply(hard, function(rt) {
      s <- sftplus:::.sft_ucip_score(rt); s$numer / s$variance }, numeric(1)))
  prec <- c(vapply(easy, function(rt) sftplus:::.sft_ucip_score(rt)$variance, numeric(1)),
            vapply(hard, function(rt) sftplus:::.sft_ucip_score(rt)$variance, numeric(1)))
  fit <- sftplus:::.sft_normal_hierarchy_fit(
    hat, prec, paste0("s", seq_along(hat)), "theta", "InvGamma",
    3000L, 500L, 1L, 2L, 0, .5, 2, .25, .5, 4, .95, 12L, list(), .94, NULL,
    group = rep(1:2, each = 6), group_levels = c("easy", "hard"))
  expect_equal(fit$mu_names, c("mu[easy]", "mu[hard]"))
  expect_equal(ncol(fit$mu_matrix), 2L)
  # The easy condition is further from the UCIP benchmark, so its mean is larger.
  expect_gt(mean(fit$mu_matrix[, "easy"]), mean(fit$mu_matrix[, "hard"]))
  contrasts <- sftplus:::.sft_condition_contrasts(fit$mu_matrix,
                                                  c("easy", "hard"), .94, .05, .5)
  expect_equal(nrow(contrasts), 1L)
  expect_equal(contrasts$contrast, "easy - hard")
  expect_gt(contrasts$posterior_positive, .9)
})


test_that("micGroup.bayes fits several conditions and contrasts them", {
  set.seed(3)
  n <- 60
  mk <- function(i, cond, shift) data.frame(
    Subject = paste0("s", i), Condition = cond,
    RT = c(stats::rlnorm(n, log(.42 * shift), .2), stats::rlnorm(n, log(.55), .2),
           stats::rlnorm(n, log(.55), .2), stats::rlnorm(n, log(.75), .2)),
    Correct = 1, Channel1 = rep(c(2, 2, 1, 1), each = n),
    Channel2 = rep(c(2, 1, 2, 1), each = n))
  d <- do.call(rbind, c(lapply(1:6, mk, cond = "easy", shift = 1),
                        lapply(1:6, mk, cond = "hard", shift = 1.22)))
  m <- micGroup.bayes(d, ndraws = 2500, burnin = 400, chains = 2, seed = 9,
                      rope = .025)
  expect_identical(m$conditions, c("easy", "hard"))
  expect_equal(m$population_summary$parameter, c("mu[easy]", "mu[hard]", "tau"))
  expect_equal(nrow(m$condition_contrasts), 1L)
  expect_true("Condition" %in% names(m$score))
  expect_named(m$bayes_factor, c("easy", "hard"))
  expect_true(all(c("population_overadditive", "subject_overadditive") %in%
                    names(m$posterior_probability)))
  # Selecting one condition restores the single-mean output.
  one <- micGroup.bayes(d, Condition = "easy", ndraws = 1500, burnin = 300,
                        chains = 2, seed = 9, rope = .025)
  expect_equal(one$population_summary$parameter, c("mu", "tau"))
  expect_null(one$condition_contrasts)
})


test_that("prior_sensitivity traces the Jeffreys-Lindley dependence", {
  set.seed(31)
  rt <- lapply(1:7, function(i) list(stats::rexp(120, 4.2), stats::rexp(120, 1.6),
                                     stats::rexp(120, 1.6)))
  names(rt) <- paste0("s", 1:7)
  g <- capacityGroup.bayes(rt, ndraws = 2500, burnin = 400, chains = 2,
                           seed = 3, rope = .1)
  ps <- prior_sensitivity(g, ndraws = 1500, burnin = 300, chains = 2, seed = 11)
  expect_s3_class(ps, "sft_prior_sensitivity")
  expect_true(any(ps$is_fitted))
  expect_equal(attr(ps, "reference_prior_sd"), g$prior$mu$sd)
  # A more diffuse alternative must not increase evidence against the point null.
  expect_lt(ps$BF10_point[nrow(ps)], ps$BF10_point[which(ps$is_fitted)])
  expect_true(all(is.finite(ps$BF10_interval)))
  # Widening the prior also widens the posterior toward the data estimate.
  expect_gt(ps$mean[nrow(ps)], ps$mean[1L])
})


test_that("spike_slab separates genuinely non-null participants", {
  set.seed(21)
  rt <- c(lapply(1:5, function(i) super_subject()),
          lapply(1:5, function(i) ucip_or_subject()))
  names(rt) <- c(paste0("super", 1:5), paste0("ucip", 1:5))
  g <- capacityGroup.bayes(rt, ndraws = 3000, burnin = 400, chains = 2, seed = 5)
  ss <- spike_slab(g, ndraws = 3000, burnin = 600, chains = 2, seed = 7)
  expect_equal(nrow(ss$summary), 10L)
  sup <- ss$summary$inclusion_probability[1:5]
  nul <- ss$summary$inclusion_probability[6:10]
  expect_true(all(sup > .75))
  expect_true(all(nul < .75))
  expect_gt(min(sup), max(nul))
  expect_true(all(ss$summary$BF10_inclusion[1:5] > 1))
  expect_equal(ss$prior_inclusion_probability, .5)
  expect_true(all(is.finite(ss$summary$BF10_inclusion)))
  expect_true(all(ss$population$parameter == c("w", "mu", "tau")))
  expect_true(all(ss$diagnostics$rhat < 1.1))
  # The spike is exact: excluded draws contribute a hard zero, not a small value.
  expect_true(any(ss$draws$theta[ss$draws$z == 0] == 0))
  expect_true(all(ss$draws$theta[ss$draws$z == 0] == 0))
})


test_that("spike_slab needs at least two participants", {
  set.seed(2)
  one <- ucip.bayes(ucip_or_subject(60), ndraws = 400, chains = 2, seed = 1)
  expect_error(spike_slab(one), "at least two participants")
})


test_that("WAIC is computed correctly from a pointwise log-likelihood", {
  set.seed(8)
  ll <- matrix(stats::rnorm(200 * 12, -2, .3), nrow = 200, ncol = 12)
  w <- sftplus:::.sft_bayes_waic(ll)
  lppd <- sum(apply(ll, 2, function(x) log(mean(exp(x)))))
  p <- sum(apply(ll, 2, stats::var))
  expect_equal(w$lppd, lppd)
  expect_equal(w$p_waic, p)
  expect_equal(w$elpd_waic, lppd - p)
  expect_equal(w$waic, -2 * (lppd - p))
  expect_equal(w$n_observations, 12L)
  # Zero-exposure columns are exact zeros on every draw and must be dropped
  # rather than counted as observations.
  ll2 <- cbind(ll, matrix(0, nrow = 200, ncol = 3))
  w2 <- sftplus:::.sft_bayes_waic(ll2)
  expect_equal(w2$n_observations, 12L)
  expect_equal(w2$dropped_zero_exposure, 3L)
  expect_equal(w2$waic, w$waic)
})


test_that("sft_waic errors helpfully when no likelihood was stored", {
  fake <- structure(list(waic = NULL), class = c("sft_bayes", "list"))
  expect_error(sft_waic(fake), "no pointwise log-likelihood")
  expect_error(sft_waic(list()), "fitted sft_bayes object")
})
