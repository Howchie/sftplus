test_that("row-wise aliases are normalized to canonical SFT columns", {
  d <- data.frame(
    subjects = rep(c("P1", "P2"), each = 9L),
    LogicalRule = "OR",
    rt = rep(seq_len(9L), 2L),
    Correct = TRUE,
    Channel1 = rep(c(1, 2, 0, 1, 2, 0, 1, 2, 0), 2L),
    Channel2 = rep(c(1, 0, 2, 1, 0, 2, 1, 0, 2), 2L),
    stringsAsFactors = FALSE
  )
  canonical <- d
  names(canonical)[match(c("subjects", "LogicalRule", "rt"), names(canonical))] <-
    c("Subject", "Condition", "RT")

  alias_names <- names(sftplus:::.sft_normalize_columns(d))
  expect_true(all(c("Subject", "Condition", "RT") %in% alias_names))

  alias_rt <- sft_data_to_rt(d)
  canonical_rt <- sft_data_to_rt(canonical)
  expect_equal(alias_rt$RT, canonical_rt$RT)
  expect_equal(alias_rt$CR, canonical_rt$CR)
  expect_equal(attr(alias_rt, "subjects"), c("P1", "P2"))
  expect_equal(attr(alias_rt, "Condition"), "OR")

  alias_capacity <- capacity.or(d[d$subjects == "P1", , drop = FALSE])
  canonical_capacity <- capacity.or(canonical[canonical$Subject == "P1", , drop = FALSE])
  expect_equal(alias_capacity$Ct(c(1, 2)), canonical_capacity$Ct(c(1, 2)))
})


test_that("aliases work through Bayesian and list-adapter entry points", {
  d <- data.frame(
    subjects = rep(c("P1", "P2"), each = 9L),
    LogicalRule = "OR",
    rt = rep(seq_len(9L), 2L),
    Correct = TRUE,
    Channel1 = rep(c(1, 2, 0, 1, 2, 0, 1, 2, 0), 2L),
    Channel2 = rep(c(1, 0, 2, 1, 0, 2, 1, 0, 2), 2L),
    stringsAsFactors = FALSE
  )

  prep <- semiparametricSFT.bayes(d, sample = FALSE, n_bins = 3L)
  expect_equal(prep$prepared$subjects, c("P1", "P2"))
  expect_equal(prep$prepared$selected, "OR")

  fit <- capacityGroup.bayes(d, ndraws = 100L, burnin = 20L,
                             chains = 2L, seed = 17L)
  expect_equal(fit$score$subject, c("P1", "P2"))
})
