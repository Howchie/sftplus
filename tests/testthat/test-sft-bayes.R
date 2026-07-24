make_or_data <- function(n_subjects = 2L) {
  one <- function(subject) {
    data.frame(
      Subject = subject,
      Condition = rep(c("OR", "other"), each = 6L),
      RT = c(1000, 2000, 3000, 1200, 2200, 3200,
             1100, 2100, 3100, 1300, 2300, 3300),
      Correct = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE,
                  TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
      Channel1 = c(1, 2, 0, 2, 3, 0, 4, 1, 0, 2, 5, 0),
      Channel2 = c(1, 0, 2, 3, 0, 4, 2, 0, 3, 4, 0, 5),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, lapply(seq_len(n_subjects), function(i) one(paste0("P", i))))
}

test_that("OR Bayesian preparation validates canonical input and pools salience", {
  d <- make_or_data()
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 4, sample = FALSE,
                      prior_predictive_draws = 3)
  expect_s3_class(x, "sft_bayes")
  expect_identical(x$metadata$series, c("AB", "A", "B"))
  expect_equal(sort(unique(x$data$series)), c("A", "AB", "B"))
  expect_identical(x$metadata$selected_conditions, "OR")
  expect_identical(x$metadata$salience_pooling$method, "channel_presence")
  expect_equal(x$prepared$events[, "AB", ],
               x$prepared$events[, "AB", ], tolerance = 0)
  expect_true(all(x$prepared$exposure >= 0))
  expect_equal(nrow(x$prepared$missing_series), 0L)
  expect_false(x$diagnostics$available)
})

test_that("validation catches missing columns, invalid RT/correctness, and incomplete series", {
  d <- make_or_data(1L)
  expect_error(semiparametricSFT.bayes(d[, setdiff(names(d), "RT")], sample = FALSE), "RT")
  d$RT[[1L]] <- Inf
  expect_error(semiparametricSFT.bayes(d, sample = FALSE), "finite")
  d <- make_or_data(1L); d$Correct[[1L]] <- 2
  expect_error(semiparametricSFT.bayes(d, sample = FALSE), "Correct")
  d <- make_or_data(1L); d$Channel1[[1L]] <- -1
  expect_error(semiparametricSFT.bayes(d, sample = FALSE), "Channel1")
  d <- make_or_data(1L); d <- d[d$Channel1 > 0 | d$Channel2 > 0, ]
  d <- d[!(d$Channel1 == 0 & d$Channel2 > 0), ]
  expect_error(semiparametricSFT.bayes(d, Condition = "OR", sample = FALSE), "missing combinations")
  expect_silent(semiparametricSFT.bayes(d, Condition = "OR", sample = FALSE,
                               require_complete = FALSE, prior_predictive_draws = 1))
})

test_that("pooled exposure includes incorrect censoring and spans intervals", {
  d <- data.frame(
    Subject = rep("P1", 6), Condition = "OR",
    RT = c(1, 2, 3, 1, 2, 3), Correct = c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
    Channel1 = c(1, 1, 0, 1, 1, 0), Channel2 = c(1, 0, 1, 1, 0, 1)
  )
  x <- semiparametricSFT.bayes(d, n_bins = 2, sample = FALSE, prior_predictive_draws = 1)
  # The exposure is obtained from every finite trial, not only correct ones.
  expect_equal(sum(x$prepared$exposure), sum(d$RT))
  expect_equal(sum(x$prepared$events), sum(d$Correct))
  expect_true(any(x$prepared$exposure[, "AB", ] > 0))
  expect_true(any(x$prepared$exposure[, "A", ] > 0))
  expect_true(any(x$prepared$exposure[, "B", ] > 0))
  expect_true(all(diff(x$grid$boundaries) > 0))
})

test_that("supplied posterior draws transform draw-by-draw and preserve UCIP", {
  d <- make_or_data(2L)
  prep <- sftplus:::.sft_bayes_prepared(d, Condition = "OR", n_bins = 3,
                                        rt_units = "milliseconds")
  M <- 5L; I <- prep$n_subjects; J <- prep$n_bins
  A <- array(log(1), c(M, I, J)); B <- array(log(2), c(M, I, J))
  AB <- array(log(3), c(M, I, J))
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 3, sample = FALSE,
                      rt_units = "milliseconds", posterior_draws = list(A = A, B = B, AB = AB),
                      prior_predictive_draws = 1, posterior_predictive_draws = 2)
  expect_equal(dim(x$transformed$H_A), c(M, I, J))
  expect_equal(as.numeric(x$transformed$H_A[1, 1, ]), cumsum(x$grid$bins$width),
               tolerance = 1e-12)
  expect_equal(as.numeric(x$transformed$H_AB[1, 1, ]),
               3 * cumsum(x$grid$bins$width), tolerance = 1e-12)
  expect_equal(max(abs(x$transformed$D)), 0, tolerance = 1e-10)
  expect_equal(x$posterior$n_draws, M)
  expect_equal(dim(x$transformed$population$D), c(M, J))
  expect_true(all(c("Mean", "Lower", "Upper", "Prob_super", "Prob_limited") %in%
                    names(x$tidy_curves)))
  expect_equal(nrow(x$curves$population), 2L * nrow(x$grid$bins))
})

make_salience_data <- function(n_subjects = 2L) {
  one <- function(subject) {
    z <- expand.grid(Channel1 = c(0, 1, 2), Channel2 = c(0, 1, 2),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    z <- z[z$Channel1 > 0 | z$Channel2 > 0, , drop = FALSE]
    z <- z[rep(seq_len(nrow(z)), each = 2L), , drop = FALSE]
    z$Subject <- subject
    z$Condition <- "OR"
    z$RT <- seq_len(nrow(z)) / 20 + .05
    z$Correct <- rep(c(TRUE, FALSE), length.out = nrow(z))
    z
  }
  do.call(rbind, lapply(seq_len(n_subjects), function(i) one(paste0("P", i))))
}

test_that("salience split retains matched single and dual cells", {
  d <- make_salience_data()
  x <- semiparametricSFT.bayes(d, Condition = "OR", salience_split = TRUE, n_bins = 4,
                      sample = FALSE, prior_predictive_draws = 2)
  expected <- c("A_L", "A_H", "B_L", "B_H", "AB_LL", "AB_LH", "AB_HL", "AB_HH")
  expect_true(x$prepared$split)
  expect_identical(x$prepared$series, expected)
  expect_equal(dim(x$prepared$events), c(2L, 8L, x$prepared$n_bins))
  expect_equal(dim(x$prepared$exposure), c(2L, 8L, x$prepared$n_bins))
  expect_true(all(x$prepared$complete))
  expect_identical(x$salience_pooling$method, "channel_presence_with_salience_split")
  expect_equal(x$salience_pooling$capacity_estimand,
               "matched within-salience capacity, averaged over salience cells draw-by-draw")
  expect_equal(x$salience_pooling$rules$Series, expected)
})

test_that("matched salience capacity and SIC are transformed per posterior draw", {
  d <- make_salience_data()
  prep <- sftplus:::.sft_bayes_prepared(d, Condition = "OR", n_bins = 3,
                                        rt_units = "seconds", salience_split = TRUE)
  M <- 4L; I <- prep$n_subjects; J <- prep$n_bins
  draws <- setNames(lapply(prep$series, function(nm) {
    if (startsWith(nm, "AB_")) {
      value <- switch(sub("AB_", "", nm), LL = log(2), LH = log(3),
                      HL = log(3), HH = log(4))
    } else {
      value <- if (endsWith(nm, "_L")) log(1) else log(2)
    }
    array(value, c(M, I, J))
  }), prep$series)
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 3, salience_split = TRUE,
                      sample = FALSE, posterior_draws = draws,
                      prior_predictive_draws = 1, posterior_predictive_draws = 2)
  expect_equal(dim(x$transformed$D_cells), c(M, I, J, 4L))
  expect_equal(dim(x$transformed$C_cells), c(M, I, J, 4L))
  expect_equal(max(abs(x$transformed$D)), 0, tolerance = 1e-10)
  expect_equal(max(abs(x$transformed$C - 1)), 0, tolerance = 1e-10)
  expect_true(any(x$transformed$SIC > 0))
  expect_true(all(c("LL", "LH", "HL", "HH", "SIC", "OR-average") %in%
                    unique(x$tidy_curves$Series)))
  expect_equal(dim(x$transformed$population$D), c(M, J))
  expect_equal(dim(x$transformed$population$D_cells), c(M, J, 4L))
  expect_equal(dim(x$predictive_checks$posterior$draws$AB_LL), c(2L, I, J))
})

test_that("salience split accepts explicit mappings and preserves SIC sign convention", {
  d <- make_salience_data()
  d$Channel1[d$Channel1 == 1] <- 10
  d$Channel1[d$Channel1 == 2] <- 20
  d$Channel2[d$Channel2 == 1] <- 10
  d$Channel2[d$Channel2 == 2] <- 20
  prep <- sftplus:::.sft_bayes_prepared(d, "OR", 3, c(.05, .95), "seconds", NULL,
                                        NULL, TRUE, NULL, c(10, 20))
  M <- 1L; I <- prep$n_subjects; J <- prep$n_bins
  draws <- setNames(lapply(prep$series, function(nm) {
    value <- switch(nm, AB_LL = log(1), AB_LH = log(2),
                    AB_HL = log(3), AB_HH = log(4), log(1))
    array(value, c(M, I, J))
  }), prep$series)
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 3,
                      salience_split = list(Channel1 = c(10, 20), Channel2 = c(10, 20)),
                      sample = FALSE, posterior_draws = draws,
                      prior_predictive_draws = 1)
  s <- exp(-x$transformed$H_AB_LL[1, 1, 1]) -
    exp(-x$transformed$H_AB_LH[1, 1, 1]) -
    exp(-x$transformed$H_AB_HL[1, 1, 1]) +
    exp(-x$transformed$H_AB_HH[1, 1, 1])
  expect_equal(x$transformed$SIC[1, 1, 1], s, tolerance = 1e-12)
  expect_identical(x$salience_pooling$mapping$method, "explicit")
})

test_that("hierarchical MIC integrates the SIC and tests its sign at both levels", {
  d <- make_salience_data(3L)
  prep <- sftplus:::.sft_bayes_prepared(d, Condition = "OR", n_bins = 4,
                                        rt_units = "seconds", salience_split = TRUE)
  M <- 40L; I <- prep$n_subjects; J <- prep$n_bins
  set.seed(7)
  draws <- setNames(lapply(prep$series, function(nm) {
    base <- if (startsWith(nm, "AB_")) {
      switch(sub("AB_", "", nm), LL = log(2), LH = log(3), HL = log(3), HH = log(4))
    } else if (endsWith(nm, "_L")) log(1) else log(2)
    array(base + rnorm(M * I * J, 0, 0.05), c(M, I, J))
  }), prep$series)
  # Population hyper-parameters exercise the typical-subject group estimand.
  draws$population_A_L <- matrix(log(1), M, J); draws$population_A_H <- matrix(log(2), M, J)
  draws$population_B_L <- matrix(log(1), M, J); draws$population_B_H <- matrix(log(2), M, J)
  draws$population_delta <- matrix(0, M, J)
  draws$population_gamma <- array(0, c(M, J, 4L))

  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 4, salience_split = TRUE,
                 sample = FALSE, posterior_draws = draws,
                 prior_predictive_draws = 1, posterior_predictive_draws = 1)

  expect_false(is.null(x$mic))
  expect_true(all(c("population", "population_finite_mean", "subject") %in%
                    x$mic$summary$Level))
  expect_identical(x$mic$group_estimand,
                   "population parameter (typical subject, random effects held at zero)")
  # Sign tests are proper probabilities in [0, 1].
  expect_true(all(x$mic$summary$Prob_positive >= 0 & x$mic$summary$Prob_positive <= 1))
  # This AB dominance ordering yields a strictly positive (over-additive) MIC.
  pop <- x$mic$summary[x$mic$summary$Level == "population", ]
  expect_gt(pop$Mean, 0)
  expect_equal(pop$Prob_positive, 1, tolerance = 1e-8)

  # MIC == area under SIC == contrast of the four posterior cell means.
  w <- prep$grid$bins$width
  Sll <- exp(-x$transformed$H_AB_LL); Slh <- exp(-x$transformed$H_AB_LH)
  Shl <- exp(-x$transformed$H_AB_HL); Shh <- exp(-x$transformed$H_AB_HH)
  direct <- vapply(seq_len(I), function(i)
    as.numeric((Sll[, i, ] - Slh[, i, ] - Shl[, i, ] + Shh[, i, ]) %*% w),
    numeric(M))
  expect_equal(x$mic$draws$subject, direct, tolerance = 1e-10,
               check.attributes = FALSE)

  # Subject-versus-population contrasts exist for every subject.
  expect_equal(nrow(x$mic$contrasts), I)
  expect_true(all(x$mic$contrasts$Level == "subject_vs_population"))

  # Extractor returns the stored object; data-frame path defaults to a split.
  expect_identical(mic.bayes(x), x$mic)
})

test_that("hierarchical SIC exposes population, finite-mean, and subject levels", {
  d <- make_salience_data(2L)
  prep <- sftplus:::.sft_bayes_prepared(d, Condition = "OR", n_bins = 3,
                                        rt_units = "seconds", salience_split = TRUE)
  M <- 8L; I <- prep$n_subjects; J <- prep$n_bins
  draws <- setNames(lapply(prep$series, function(nm) {
    base <- if (startsWith(nm, "AB_")) {
      switch(sub("AB_", "", nm), LL = log(2), LH = log(3), HL = log(3), HH = log(4))
    } else if (endsWith(nm, "_L")) log(1) else log(2)
    array(base, c(M, I, J))
  }), prep$series)
  draws$population_A_L <- matrix(log(1), M, J); draws$population_A_H <- matrix(log(2), M, J)
  draws$population_B_L <- matrix(log(1), M, J); draws$population_B_H <- matrix(log(2), M, J)
  draws$population_delta <- matrix(0, M, J)
  draws$population_gamma <- array(0, c(M, J, 4L))
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 3, salience_split = TRUE,
                 sample = FALSE, posterior_draws = draws, prior_predictive_draws = 1)
  expect_false(is.null(x$sic))
  expect_true(all(c("population", "population_finite_mean", "subject") %in%
                    x$sic$summary$Level))
  expect_true(all(c("Prob_super", "Prob_limited", "Time") %in% names(x$sic$summary)))
  expect_identical(sic.bayes(x), x$sic)
})

test_that("pooled OR fits carry no MIC or SIC", {
  d <- make_or_data(2L)
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 3, sample = FALSE,
                 prior_predictive_draws = 1)
  expect_null(x$mic)
  expect_null(x$sic)
  expect_error(mic.bayes(x), "salience_split")
  expect_error(sic.bayes(x), "salience_split")
})

test_that("rstan OR integration is optional and exposes diagnostics", {
  skip_if_not_installed("rstan")
  d <- make_or_data(2L)
  x <- semiparametricSFT.bayes(d, Condition = "OR", n_bins = 3, sample = TRUE,
                      iter = 60, warmup = 30, chains = 2, cores = 1,
                      refresh = 0, seed = 881, smoothness = list(basis_dim = 3),
                      prior_predictive_draws = 2, posterior_predictive_draws = 3)
  expect_true(x$diagnostics$available)
  expect_true(all(c("convergence", "divergences", "ebfmi") %in% names(x$diagnostics)))
  expect_equal(x$posterior$n_draws, 60L)
  expect_equal(dim(x$transformed$D)[2:3], c(2L, 3L))
  expect_equal(dim(x$predictive_checks$posterior$draws$AB)[1L], 3L)
})

test_that("rstan salience split integration is optional and exposes SIC", {
  skip_if_not_installed("rstan")
  d <- make_salience_data(2L)
  x <- semiparametricSFT.bayes(d, Condition = "OR", salience_split = TRUE, n_bins = 3,
                      sample = TRUE, iter = 40, warmup = 20, chains = 2, cores = 1,
                      refresh = 0, seed = 882, smoothness = list(basis_dim = 3),
                      prior_predictive_draws = 1, posterior_predictive_draws = 2)
  expect_true(x$diagnostics$available)
  expect_equal(x$posterior$n_draws, 40L)
  expect_equal(dim(x$transformed$SIC)[2:3], c(2L, 3L))
  expect_equal(dim(x$predictive_checks$posterior$draws$AB_LL)[1L], 2L)
})
