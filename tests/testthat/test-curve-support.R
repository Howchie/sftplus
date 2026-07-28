test_that(".sft_curve returns NA outside the estimated support", {
  f <- sftplus:::.sft_curve(c(10, 20, 30), c(1, 2, 3))
  expect_equal(f(c(10, 20, 30)), c(1, 2, 3))
  expect_true(all(is.na(f(c(0, 5, 9.9, 30.1, 100)))))
})


test_that(".sft_curve keeps NA markers instead of interpolating across them", {
  # The callers write NA where an estimate is undefined.  Those points used to
  # be dropped and painted over by flat extrapolation.
  f <- sftplus:::.sft_curve(c(1, 2, 3, 4, 5), c(NA, NA, 1, 2, 3))
  expect_true(all(is.na(f(c(1, 1.5, 2, 2.5)))))
  expect_equal(f(c(3, 4, 5)), c(1, 2, 3))
})


test_that("capacity C(t) is not extrapolated below the participant's data", {
  set.seed(11)
  ab <- stats::rlnorm(120, log(300), .25)
  a <- stats::rlnorm(120, log(400), .25)
  b <- stats::rlnorm(120, log(400), .25)
  fit <- capacity.or(list(ab, a, b), list(rep(TRUE, 120), rep(TRUE, 120), rep(TRUE, 120)))
  early <- min(c(ab, a, b)) - c(1, 10, 100)
  expect_true(all(is.na(fit$Ct(early))))
  expect_true(all(is.na(fit$Ct(max(c(ab, a, b)) + 1))))
  mid <- stats::median(ab)
  expect_true(is.finite(fit$Ct(mid)))
})


test_that("capacityGroup trims each curve to its own support", {
  set.seed(12)
  mk <- function(subj, shift) {
    n <- 90
    data.frame(Subject = subj, Condition = "c1",
               RT = c(stats::rlnorm(n, log(300 * shift), .2),
                      stats::rlnorm(n, log(420 * shift), .2),
                      stats::rlnorm(n, log(420 * shift), .2)),
               Correct = TRUE,
               Channel1 = rep(c(1, 1, 0), each = n),
               Channel2 = rep(c(1, 0, 1), each = n),
               stringsAsFactors = FALSE)
  }
  # One fast and one slow participant, so the pooled time grid extends well
  # past either one's own data.
  dat <- rbind(mk("fast", 1), mk("slow", 2.5))
  cg <- capacityGroup(dat, plotCt = FALSE)
  expect_true("specs" %in% names(cg))

  fast <- cg$Ct.fn[cg$overview$Subject == "fast", ]
  slow <- cg$Ct.fn[cg$overview$Subject == "slow", ]
  fast_win <- range(cg$times[is.finite(fast)])
  slow_win <- range(cg$times[is.finite(slow)])
  expect_true(fast_win[2] < slow_win[2])
  expect_true(fast_win[1] < slow_win[1])
  # The grid now starts where the earliest curve does, and the slow
  # participant is not evaluated there.
  expect_equal(fast_win[1], cg$times[1])
  expect_true(is.na(slow[1]))
  expect_true(cg$times[1] > min(dat$RT))

  wide <- capacityGroup(dat, plotCt = FALSE, trim = NULL)
  expect_true(sum(is.finite(wide$Ct.fn[1, ])) >= sum(is.finite(cg$Ct.fn[1, ])))
  expect_error(capacityGroup(dat, plotCt = FALSE, trim = c(.5, 2)), "trim must be")
})


test_that("the sicGroup time grid survives second-scale response times", {
  set.seed(13)
  mk <- function(subj, scale) {
    n <- 40
    cells <- expand.grid(Channel1 = c(1, 2), Channel2 = c(1, 2))
    do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
      rate <- 1 / (cells$Channel1[i] + cells$Channel2[i])
      data.frame(Subject = subj, Condition = "c1",
                 RT = scale * (0.3 + stats::rexp(n, 1 / rate)),
                 Correct = TRUE, Channel1 = cells$Channel1[i],
                 Channel2 = cells$Channel2[i], stringsAsFactors = FALSE)
    }))
  }
  secs <- rbind(mk("s1", 1), mk("s2", 1))
  ms <- secs; ms$RT <- ms$RT * 1000

  g_secs <- sicGroup(secs, plotSIC = FALSE)
  g_ms <- sicGroup(ms, plotSIC = FALSE)
  # round() used to collapse second-scale times to a handful of knots.
  expect_gt(length(g_secs$times), 100L)
  expect_equal(length(g_secs$times), length(g_ms$times))
  expect_equal(g_secs$sic.fn, g_ms$sic.fn, tolerance = 1e-8)

  # An explicit grid is honoured, and the default thins to n.times.
  own <- sicGroup(secs, plotSIC = FALSE, times = seq(.3, 2, length.out = 25))
  expect_equal(length(own$times), 25L)
  expect_lte(length(sicGroup(secs, plotSIC = FALSE, n.times = 50L)$times), 50L)
})


test_that("probability-scale builders read every participant at their own position", {
  set.seed(14)
  specs <- lapply(c(1, 2.5), function(shift) {
    rt <- stats::rlnorm(200, log(300 * shift), .2)
    list(label = paste0("s", shift), Condition = "c1",
         fn = function(t) 1 + (t - stats::median(rt)) / 1000,
         rt = rt, cr = rep(TRUE, 200), benchmark = rt)
  })
  q <- build_ct_prob_df(specs, probs = seq(.1, .9, by = .1))
  expect_equal(nrow(q), 18L)
  expect_true(all(is.finite(q$Ct)))
  # The package indexes probability-scale curves by `prob`, as rmi.test() does.
  expect_true(all(c("prob", "Time", "Ct", "Series", "scale") %in% names(q)))
  expect_equal(unique(q$scale), "probability")

  # Both participants contribute at every prob, which is the point: on an
  # absolute time axis the fast one has no data where the slow one lives.
  s <- build_prob_summary_df(q, "Ct")
  expect_equal(s$n, rep(2L, 9L))
  expect_true(all(diff(s$Time) > 0))
  expect_true(all(s$lower <= s$mean & s$mean <= s$upper))
  # Column names match build_rmi_violation_df() so the plotting code is shared.
  expect_true(all(c("prob", "n", "mean", "lower", "upper", "scale") %in% names(s)))

  cdf <- build_cdf_prob_df(specs, probs = seq(.1, .9, by = .1))
  # Read at its own percentiles, each participant's CDF is the identity in prob.
  expect_equal(cdf$CDF, rep(seq(.1, .9, by = .1), 2), tolerance = .03)
  expect_equal(sort(unique(cdf$Series)), c("s1", "s2.5"))
  expect_error(build_ct_prob_df(specs, probs = c(0, .5)), "strictly inside")
  expect_equal(nrow(build_ct_prob_df(specs, probs = .5, condition = "nope")), 0L)
})


test_that("capacityGroup can plot the group overlay without the individual panels", {
  set.seed(15)
  mk <- function(subj) {
    n <- 90
    data.frame(Subject = subj, Condition = "c1",
               RT = c(stats::rlnorm(n, log(300), .2), stats::rlnorm(n, log(420), .2),
                      stats::rlnorm(n, log(420), .2)),
               Correct = TRUE, Channel1 = rep(c(1, 1, 0), each = n),
               Channel2 = rep(c(1, 0, 1), each = n), stringsAsFactors = FALSE)
  }
  dat <- rbind(mk("a"), mk("b"), mk("c"))
  count_pages <- function(mode) {
    f <- tempfile(fileext = ".pdf")
    grDevices::pdf(f, compress = FALSE, onefile = TRUE)
    capacityGroup(dat, plotCt = mode)
    grDevices::dev.off()
    length(grep("/Type /Page[^s]", readLines(f, warn = FALSE, skipNul = TRUE),
                useBytes = TRUE))
  }
  # Three participants: individual gives one panel each, group gives one.
  expect_gt(count_pages("individual"), count_pages("group"))
  expect_equal(count_pages("group"), 1L)
  expect_equal(count_pages("none"), 0L)
  expect_equal(count_pages("both"), count_pages("individual") + 1L)
  # The historical logical argument still works.
  expect_equal(count_pages(TRUE), count_pages("individual"))
  expect_equal(count_pages(FALSE), 0L)
})
