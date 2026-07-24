test_that("enhanced hazard and capacity APIs retain test outputs", {
  set.seed(101)
  rt <- list(rexp(30, 2), rexp(30, 1.8), rexp(30, 1.6))
  cr <- lapply(rt, function(x) rep(TRUE, length(x)))
  nah <- estimateNAH(rt[[1]])
  expect_length(nah$H(c(.1, .5)), 2)
  out <- capacity.or(rt, cr)
  expect_true(all(c("Ct", "Ctest", "Ct_upper", "Ct_lower") %in% names(out)))
  expect_s3_class(out$Ctest, "htest")
})

test_that("row-wise SFT data convert to list input and are accepted by UCIP APIs", {
  d <- data.frame(
    Subject = rep(c("P1", "P2"), each = 9L), Condition = "OR",
    RT = rep(seq_len(9L), 2L),
    Correct = rep(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE), 2L),
    Channel1 = rep(c(1, 2, 0, 1, 2, 0, 1, 2, 0), 2L),
    Channel2 = rep(c(1, 0, 2, 1, 0, 2, 1, 0, 2), 2L),
    stringsAsFactors = FALSE
  )
  converted <- sft_data_to_rt(d)
  expect_identical(names(converted$RT), c("P1", "P2"))
  expect_identical(names(converted$RT[[1L]]), c("AB", "A", "B"))
  expect_equal(unname(vapply(converted$RT[[1L]], length, integer(1))), c(3L, 3L, 3L))
  expect_equal(unname(vapply(converted$CR[[2L]], length, integer(1))), c(3L, 3L, 3L))
  expect_error(ucip.test(d), "single Subject")

  from_list <- ucip.test(converted$RT[[1L]], converted$CR[[1L]])
  from_data <- ucip.test(d, Subject = "P1")
  expect_equal(from_data$statistic, from_list$statistic, tolerance = 1e-12)

  ans <- capacityGroup.bayes(d, ndraws = 120, burnin = 30, chains = 2, seed = 111)
  expect_equal(nrow(ans$subject_summary), 2L)
  expect_equal(ans$score$subject, c("P1", "P2"))
  expect_true(is.function(assessment(d, Subject = "P1", stopping.rule = "OR",
                                     correct = TRUE, fast = TRUE)))
})

test_that("frequentist SIC tests and Bayesian SIC route are available", {
  set.seed(102)
  x <- list(rexp(24, 2), rexp(24, 1.8), rexp(24, 1.7), rexp(24, 1.5))
  fit <- sic(x[[1]], x[[2]], x[[3]], x[[4]], nsamp = 100, maxn = 200, seed = 4)
  expect_true(is.function(fit$SIC))
  expect_s3_class(fit$MICtest, "htest")
  expect_equal(nrow(fit$Dominance), 8)
  bf <- sictestBayes(x[[1]], x[[2]], x[[3]], x[[4]], nsamp = 100, maxn = 200, seed = 4)
  expect_length(bf$statistic, 6)
  dp <- siDominance(x[[1]], x[[2]], x[[3]], x[[4]], method = "dp", nbin = 4, nsamp = 100, seed = 4)
  expect_equal(nrow(dp), 8)
})

test_that("hierarchical UCIP Bayesian companion uses score information", {
  set.seed(103)
  rt <- lapply(1:3, function(i) list(rexp(20, 2), rexp(20, 1.8), rexp(20, 1.6)))
  names(rt) <- paste0("P", 1:3)
  ans <- capacityGroup.bayes(rt, ndraws = 200, burnin = 50, chains = 2,
                         rope = .05, seed = 5)
  expect_true(all(c("U", "V", "theta_hat", "Cz") %in% names(ans$score)))
  expect_equal(ans$score$Cz[1], ucip.test(rt[[1]])$statistic[[1]])
  expect_equal(nrow(ans$subject_summary), 3)
  expect_true(all(c("population_super", "population_limited", "population_rope") %in%
                    names(ans$posterior_probability)))
  expect_true(all(is.finite(ans$diagnostics$rhat)))
  expect_true(all(c("ess_bulk", "ess_tail") %in% names(ans$diagnostics)))
  expect_equal(nrow(ans$diagnostics), 5)
  expect_true(all(c("mean_weight_on_individual", "lower", "upper") %in%
                    names(ans$shrinkage_summary)))
  expect_equal(nrow(ans$posterior_predictive$draws), 400)
})

test_that("log-capacity and score scales agree near the UCIP model", {
  # For OR, independent exponential A/B times imply an AB minimum with rate
  # lambda_A + lambda_B, so this is a direct UCIP data-generating case.
  set.seed(20260723)
  rt <- list(rexp(80, 3.5), rexp(80, 1.5), rexp(80, 2))
  cr <- lapply(rt, function(x) rep(TRUE, length(x)))
  score <- sftplus:::.sft_ucip_score(rt, cr, stopping.rule = "OR")
  input <- list(RT = list(rt), CR = list(cr), subject = "P1")
  score_data <- sftplus:::.sft_ucip_score_data(
    input, "OR", method = "score"
  )
  capacity_data <- sftplus:::.sft_ucip_score_data(
    input, "OR", method = "capacity"
  )

  # Both paths reuse the same score decomposition; only the observation
  # scale and precision differ.
  expect_equal(score_data$theta_hat, score$numer / score$variance)
  expect_equal(score_data$precision, score$variance)
  expect_equal(capacity_data$eta_hat, score$log_capacity)
  expect_equal(capacity_data$precision, score$log_capacity_precision)

  eta_z <- capacity_data$eta_hat * sqrt(capacity_data$precision)
  expect_equal(eta_z, unname(score$statistic), tolerance = 0.25)
})

test_that("single-subject UCIP Bayesian companion uses the analytic posterior", {
  set.seed(104)
  rt <- list(rexp(30, 2), rexp(30, 1.8), rexp(30, 1.6))
  ans <- capacityGroup.bayes(rt, ndraws = 400, chains = 2, seed = 7,
                         prior_mean = .2, prior_sd = .8, rope = .05)
  score <- ans$score
  post_var <- 1 / (score$V + 1 / .8^2)
  post_mean <- post_var * (score$V * score$theta_hat + .2 / .8^2)
  expect_identical(ans$method_code, "Analytic")
  expect_null(ans$population_summary)
  expect_false(any(c("mu", "tau") %in% names(ans$draws)))
  expect_equal(ans$summary$parameter, "theta[1]")
  expect_equal(ans$summary$mean, post_mean, tolerance = .08)
  expect_equal(nrow(ans$draws), 800)
  expect_equal(ans$score$Cz, ucip.test(rt)$statistic[[1]])
  expect_equal(ans$diagnostics$rhat, 1)
})

test_that("hierarchical method aliases retain the Gibbs implementation", {
  rt <- lapply(1:2, function(i) list(rexp(20, 2), rexp(20, 1.8), rexp(20, 1.6)))
  ans <- capacityGroup.bayes(rt, ndraws = 120, burnin = 30, chains = 2,
                         method = "gibbs", seed = 8,
                         tau2_prior_shape = 2, tau2_prior_rate = 1)
  expect_identical(ans$method_code, "InvGamma")
  expect_true(all(is.finite(ans$diagnostics$ess_bulk)))
})

test_that("Stan UCIP alternatives expose centered and non-centred models", {
  skip_if_not_installed("rstan")
  rt <- lapply(1:2, function(i) list(rexp(20, 2), rexp(20, 1.8), rexp(20, 1.6)))
  for (method in c("HalfNormal", "Centered")) {
    ans <- capacityGroup.bayes(rt, ndraws = 100, burnin = 50, chains = 2,
                           method = method, seed = 9)
    expect_identical(ans$method_code, method)
    expect_equal(nrow(ans$draws), 200)
    expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(ans$diagnostics)))
    expect_true(!is.null(ans$stan_fit))
  }
})

test_that("LBA simulation is packaged", {
  skip_if_not_installed("rtdists")
  p <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
         ve = .8, sv_c = .1, sv_e = .1)
  ans <- simulate_sft(model = "lba", n = 3, p_vec = p,
                      design = c("AB", "AN", "NB", "NN"), logical_rules = "OR")
  expect_true(all(c("data", "by_rule", "metrics") %in% names(ans)))
  expect_gt(nrow(ans$data), 0)
})

test_that("simulate_sft_group generates per-subject LBA data", {
  skip_if_not_installed("rtdists")
  base <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
            ve = .8, sv_c = .1, sv_e = .1)

  # Vector params recycled across subjects, auto-named SubjectXXX.
  set.seed(11)
  grp <- simulate_sft_group(model = "lba", n = 5, params = base,
                            nSubjects = 3, design = c("AB", "AN", "NB", "NN"),
                            logical_rules = "OR")
  expect_s3_class(grp, "data.frame")
  expect_true(all(c("RT", "Channel1", "Channel2", "Correct", "Subject", "Condition")
                  %in% names(grp)))
  expect_identical(levels(grp$Subject), c("Subject001", "Subject002", "Subject003"))

  # Single-subject call collapses to one labelled subject.
  set.seed(11)
  one <- simulate_sft_group(model = "lba", n = 5, params = base,
                            design = c("AB", "AN", "NB", "NN"),
                            logical_rules = "OR")
  expect_identical(levels(one$Subject), "Subject001")

  # Per-subject parameter matrix: rows drive distinct subjects, names honoured.
  pmat <- rbind(s_fast = base, s_slow = base)
  pmat["s_slow", "vc_yes"] <- 1.0
  set.seed(12)
  grp2 <- simulate_sft_group(model = "lba", n = 20, params = pmat,
                             design = c("AB", "AN", "NB", "NN"),
                             logical_rules = "OR")
  expect_identical(levels(grp2$Subject), c("s_fast", "s_slow"))
  expect_equal(length(unique(grp2$Subject)), 2L)

  # Directly consumable by the group-level API.
  expect_silent(capacityGroup(grp2, OR = TRUE, plotCt = FALSE))
})

test_that("simulate_sft_group validates its parameter inputs", {
  base <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
            ve = .8, sv_c = .1, sv_e = .1)
  expect_error(simulate_sft_group("lba", n = 1, params = unname(base), nSubjects = 2),
               "named numeric vector")
  expect_error(simulate_sft_group("lba", n = 1, params = matrix(1:4, 2)),
               "named columns")
  expect_error(simulate_sft_group("lba", n = 1, params = rbind(base, base), nSubjects = 3),
               "does not match")
})

test_that("LBA shared capacity draws correlated AB target drifts", {
  skip_if_not_installed("rtdists")
  pars <- matrix(c(
    1, 1.5, .2, 1.4, .05,
    1, 1.5, .2, .7, .05,
    1, 1.5, .2, 1.2, .06,
    1, 1.5, .2, .7, .05
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("A", "n_A", "B", "n_B"), c("A", "b", "t0", "v", "sv")))
  shared <- list(kappa = 1.2, tau = .5, active = TRUE)
  set.seed(105)
  drifts <- sftplus:::.lba_shared_target_drifts(20000, pars, shared)
  expect_equal(mean(drifts[, "A"]), 1.2 * 1.4, tolerance = .03)
  expect_equal(mean(drifts[, "B"]), 1.2 * 1.2, tolerance = .03)
  expect_equal(stats::cov(drifts)["A", "B"], .5^2 * 1.4 * 1.2,
               tolerance = .03)
  expect_gt(stats::cor(drifts[, "A"], drifts[, "B"]), .9)
})

test_that("shared LBA capacity is AB-only and defaults preserve the ordinary path", {
  skip_if_not_installed("rtdists")
  p <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
         ve = .8, sv_c = .1, sv_e = .1)
  set.seed(106)
  ordinary <- simulate_sft("lba", n = 10, p_vec = p,
                           design = c("AB", "AN", "NB", "NN"), logical_rules = "OR")
  set.seed(106)
  explicit_default <- simulate_sft("lba", n = 10, p_vec = p,
                                   design = c("AB", "AN", "NB", "NN"),
                                   logical_rules = "OR", kappa = 1, tau = 0)
  expect_identical(ordinary$data, explicit_default$data)

  p2 <- c(p, kappa = 1.3, tau = .4)
  set.seed(107)
  ans <- simulate_sft("lba", n = 20, p_vec = p2,
                      design = c("AB", "AN", "NB", "NN"), logical_rules = "OR")
  expect_equal(as.integer(table(ans$data$S)), c(20L, 20L, 20L, 20L))
  expect_true(all(is.finite(ans$data$RT)))
  expect_error(simulate_sft("lba", n = 2, p_vec = p, kappa = 0,
                            design = "AB", logical_rules = "OR"),
               "greater than zero")
  expect_error(simulate_sft("lba", n = 2, p_vec = p, tau = -.1,
                            design = "AB", logical_rules = "OR"),
               "non-negative")
})

# Build a salience 2x2 (Channel1/Channel2 in {2 = High, 1 = Low}) data frame with
# a multiplicative per-subject speed factor so the relative MIC must divide it out.
make_mic_data <- function(subjects = paste0("P", 1:4), base = c(.4, .6, .5, .8),
                          n = 40L, seed = 202) {
  set.seed(seed)
  cells <- list(HH = c(2, 2), HL = c(2, 1), LH = c(1, 2), LL = c(1, 1))
  mult <- c(HH = 1.05, HL = 1.00, LH = 1.00, LL = 0.95)
  do.call(rbind, Map(function(subj, b) {
    do.call(rbind, lapply(names(cells), function(nm) {
      ch <- cells[[nm]]
      data.frame(Subject = subj, Condition = "A",
                 RT = abs(stats::rnorm(n, b * mult[[nm]], b * .15)),
                 Correct = TRUE, Channel1 = ch[1], Channel2 = ch[2],
                 stringsAsFactors = FALSE)
    }))
  }, subjects, base))
}

test_that("micGroup.bayes pools the dimensionless relative MIC hierarchically", {
  d <- make_mic_data()
  ans <- micGroup.bayes(d, ndraws = 400, burnin = 100, chains = 2, seed = 42,
                        rope = .025)
  expect_identical(ans$method_code, "InvGamma")
  expect_match(ans$estimand, "relative MIC")
  expect_equal(nrow(ans$subject_summary), 4L)
  expect_equal(ans$summary$parameter, c("mu", "tau", paste0("rho[P", 1:4, "]")))
  expect_true(all(c("population_overadditive", "population_underadditive",
                    "population_additive", "subject_overadditive") %in%
                    names(ans$posterior_probability)))
  # relative_MIC is exactly MIC / equal-weight grand-mean RT.
  expect_equal(ans$score$relative_MIC, ans$score$MIC / ans$score$mean_RT)
  # The MIC contrast matches the frequentist mic.test on the same cells.
  cells <- sftplus:::.sft_mic_cells_from_data(d, Subject = "P1")$cells[[1L]]
  mt <- mic.test(cells$HH, cells$HL, cells$LH, cells$LL)
  expect_equal(ans$score$MIC[ans$score$subject == "P1"], unname(mt$statistic))
})

test_that("relative MIC is invariant to multiplicative participant speed", {
  cells <- list(HH = rexp(60, 1 / .5), HL = rexp(60, 1 / .45),
                LH = rexp(60, 1 / .46), LL = rexp(60, 1 / .4))
  slow <- lapply(cells, function(x) 2 * x)
  s1 <- sftplus:::.sft_mic_relative_score(cells)
  s2 <- sftplus:::.sft_mic_relative_score(slow)
  expect_equal(s1$rho, s2$rho)
  expect_equal(2 * s1$mic, s2$mic)
})

test_that("micGroup.bayes single subject uses the analytic posterior", {
  d <- make_mic_data()
  ans <- micGroup.bayes(d, Subject = "P1", ndraws = 400, chains = 2, seed = 7,
                        rope = .025)
  expect_identical(ans$method_code, "Analytic")
  expect_null(ans$population_summary)
  expect_equal(ans$summary$parameter, "rho[P1]")
  expect_equal(nrow(ans$draws), 800L)
  expect_true(all(c("subject_overadditive", "subject_underadditive",
                    "subject_additive") %in% names(ans$posterior_probability)))
})

test_that("bootstrap var_method reports both variances and is close to analytic", {
  d <- make_mic_data()
  boot <- micGroup.bayes(d, ndraws = 300, burnin = 100, chains = 2, seed = 42,
                         var_method = "bootstrap", n_boot = 200)
  expect_identical(boot$var_method, "bootstrap")
  expect_true("se_analytic" %in% names(boot$score))
  # Bootstrap and analytic SEs agree to within a factor of two here.
  expect_true(all(boot$score$se / boot$score$se_analytic > .5 &
                    boot$score$se / boot$score$se_analytic < 2))

  rt <- list(rexp(60, 3.5), rexp(60, 1.5), rexp(60, 2))
  u <- ucip.bayes(rt, ndraws = 300, chains = 2, seed = 3,
                  var_method = "bootstrap", n_boot = 200)
  expect_identical(u$var_method, "bootstrap")
  expect_true("se_analytic" %in% names(u$score))
  expect_error(ucip.bayes(rt, var_method = "bootstrap", n_boot = 10),
               "n_boot")
})
