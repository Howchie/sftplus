test_that("serial race combination implements self-termination and exhaustion", {
  rts <- matrix(c(
    1, 3, 2, 4, # both targets
    3, 1, 2, 4, # A absence, B target
    1, 3, 4, 2, # A target, B absence
    3, 1, 4, 2  # both absences
  ), nrow = 4, byrow = TRUE,
  dimnames = list(NULL, c("A", "n_A", "B", "n_B")))

  or_st <- sftplus:::.serial_combine_races(rts, "OR", "self-terminating", "A")
  expect_equal(or_st$RT, c(1, 3, 1, 3))
  and_st <- sftplus:::.serial_combine_races(rts, "AND", "self-terminating", "A")
  expect_equal(and_st$RT, c(3, 1, 3, 1))
  exhaustive <- sftplus:::.serial_combine_races(rts, "OR", "exhaustive", "A")
  expect_equal(exhaustive$RT, rep(3, 4))
})

test_that("serial simulation exposes order metadata and disables interactions", {
  skip_if_not_installed("rtdists")
  p <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
         ve = .8, sv_c = .1, sv_e = .1, capacity_target = 99,
         capacity_absence = 99, kappa = 4, tau = 2, rho = -9)
  set.seed(77)
  serial <- simulate_sft("lba", n = 40, p_vec = p, design = "AB",
                         logical_rules = "OR", architecture = "serial",
                         serial_mode = "exhaustive", serial_order = "B")
  expect_true(all(serial$data$Architecture == "serial"))
  expect_true(all(serial$data$SerialMode == "exhaustive"))
  expect_true(all(serial$data$SerialOrder == "B"))

  p0 <- p; p0[c("capacity_target", "capacity_absence")] <- 0
  p0["kappa"] <- 1; p0["tau"] <- 0; p0["rho"] <- 0
  set.seed(77)
  serial0 <- simulate_sft("lba", n = 40, p_vec = p0, design = "AB",
                          logical_rules = "OR", architecture = "serial",
                          serial_mode = "exhaustive", serial_order = "B")
  expect_equal(serial$data$RT, serial0$data$RT)

  set.seed(78)
  random <- simulate_sft("lba", n = 1000, p_vec = p, design = "AB",
                         logical_rules = "OR", architecture = "serial",
                         serial_mode = "exhaustive", serial_order = "random",
                         p_A_first = .8)
  expect_equal(mean(random$data$SerialOrder == "A"), .8, tolerance = .06)
})

test_that("parallel remains the default architecture", {
  skip_if_not_installed("rtdists")
  p <- c(A = 1, b = 1.5, t0 = .2, vc_yes = 1.5, vc_no = .5,
         ve = .8, sv_c = .1, sv_e = .1)
  set.seed(79)
  implicit <- simulate_sft("lba", n = 3, p_vec = p, design = "AB", logical_rules = "OR")
  set.seed(79)
  explicit <- simulate_sft("lba", n = 3, p_vec = p, design = "AB", logical_rules = "OR",
                           architecture = "parallel")
  expect_equal(implicit$data, explicit$data)
  expect_false(any(c("Architecture", "SerialMode", "SerialOrder", "p_A_first") %in%
                   names(implicit$data)))
})
