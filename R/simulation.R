# Simulation utilities adapted from the CapacityAND Sim_* modules.
# This is the package's canonical simulation implementation.

.param_get <- function(p_vec, name, default = NA_real_) {
  if (name %in% names(p_vec)) return(unname(as.numeric(p_vec[[name]])))
  default
}

correctfun <- function(d) {
  cond <- as.character(d$Condition); s <- as.character(d$S); r <- as.character(d$R)
  out <- logical(length(cond))
  i <- cond == "OR"; out[i] <- (s[i] == "NN" & r[i] == "no") | (s[i] != "NN" & r[i] == "yes")
  i <- cond == "AND"; out[i] <- (s[i] == "AB" & r[i] == "yes") | (s[i] != "AB" & r[i] == "no")
  i <- cond == "XOR"; out[i] <- (s[i] %in% c("AN", "NB") & r[i] == "yes") | (s[i] %in% c("AB", "NN") & r[i] == "no")
  i <- cond == "ID"; out[i] <- s[i] == r[i]
  out
}

apply_logic_gates <- function(LogicalRule, A_t, nA_t, B_t, nB_t) {
  if (!LogicalRule %in% c("OR", "AND", "XOR", "ID")) stop("Unknown LogicalRule.")
  A_yes <- A_t < nA_t; B_yes <- B_t < nB_t
  tieA <- is.finite(A_t) & A_t == nA_t; tieB <- is.finite(B_t) & B_t == nB_t
  if (any(tieA)) A_yes[tieA] <- stats::runif(sum(tieA)) < .5
  if (any(tieB)) B_yes[tieB] <- stats::runif(sum(tieB)) < .5
  tA <- pmin(A_t, nA_t); tB <- pmin(B_t, nB_t)
  if (LogicalRule == "OR") {
    R <- ifelse(A_yes | B_yes, "yes", "no")
    RT <- ifelse(R == "yes", pmin(ifelse(A_yes, A_t, Inf), ifelse(B_yes, B_t, Inf)), pmax(nA_t, nB_t))
  } else if (LogicalRule == "AND") {
    R <- ifelse(A_yes & B_yes, "yes", "no")
    RT <- ifelse(R == "yes", pmax(A_t, B_t), pmin(ifelse(!A_yes, nA_t, Inf), ifelse(!B_yes, nB_t, Inf)))
  } else if (LogicalRule == "XOR") {
    R <- ifelse(xor(A_yes, B_yes), "yes", "no"); RT <- pmax(tA, tB)
  } else {
    R <- ifelse(!A_yes & !B_yes, "NN", ifelse(A_yes & !B_yes, "AN", ifelse(!A_yes & B_yes, "NB", "AB")))
    RT <- pmax(tA, tB)
  }
  keep <- is.finite(RT)
  data.frame(R = R[keep], RT = RT[keep], A_yes = A_yes[keep], B_yes = B_yes[keep],
             Channel1_RT = tA[keep], Channel2_RT = tB[keep], stringsAsFactors = FALSE)
}

# Assemble sequential trials from independently generated channel finish
# times. Both channels are drawn with the existing RT-generating machinery;
# this function controls only order and logical stopping.
.serial_combine_races <- function(rts, LogicalRule,
                                  serial_mode = c("self-terminating", "exhaustive"),
                                  serial_order = c("random", "A", "B"),
                                  p_A_first = .5) {
  serial_mode <- match.arg(serial_mode); serial_order <- match.arg(serial_order)
  if (!is.matrix(rts) || !all(c("A", "n_A", "B", "n_B") %in% colnames(rts))) {
    stop("Serial simulation requires finish times for A, n_A, B, and n_B.", call. = FALSE)
  }
  if (!LogicalRule %in% c("OR", "AND", "XOR", "ID")) stop("Unknown LogicalRule.")
  if (!is.numeric(p_A_first) || length(p_A_first) != 1L || is.na(p_A_first) || !is.finite(p_A_first) ||
      p_A_first < 0 || p_A_first > 1) stop("p_A_first must be one finite value in [0, 1].", call. = FALSE)
  n <- nrow(rts)
  A_t <- as.numeric(rts[, "A"]); nA_t <- as.numeric(rts[, "n_A"])
  B_t <- as.numeric(rts[, "B"]); nB_t <- as.numeric(rts[, "n_B"])
  A_yes <- A_t < nA_t; B_yes <- B_t < nB_t
  tieA <- is.finite(A_t) & A_t == nA_t; tieB <- is.finite(B_t) & B_t == nB_t
  if (any(tieA)) A_yes[tieA] <- stats::runif(sum(tieA)) < .5
  if (any(tieB)) B_yes[tieB] <- stats::runif(sum(tieB)) < .5
  A_rt <- pmin(A_t, nA_t); B_rt <- pmin(B_t, nB_t)
  A_first <- switch(serial_order, A = rep(TRUE, n), B = rep(FALSE, n), random = stats::runif(n) < p_A_first)
  first_yes <- ifelse(A_first, A_yes, B_yes); second_yes <- ifelse(A_first, B_yes, A_yes)
  first_rt <- ifelse(A_first, A_rt, B_rt); second_rt <- ifelse(A_first, B_rt, A_rt)
  exhaustive_rt <- first_rt + second_rt
  if (serial_mode == "self-terminating" && LogicalRule == "ID") {
    stop("Serial self-terminating simulation currently supports OR and AND, not ID.", call. = FALSE)
  }
  if (LogicalRule == "OR") {
    R <- ifelse(first_yes | second_yes, "yes", "no")
    RT <- if (serial_mode == "self-terminating") ifelse(first_yes, first_rt, exhaustive_rt) else exhaustive_rt
  } else if (LogicalRule == "AND") {
    R <- ifelse(first_yes & second_yes, "yes", "no")
    RT <- if (serial_mode == "self-terminating") ifelse(!first_yes, first_rt, exhaustive_rt) else exhaustive_rt
  } else if (LogicalRule == "XOR") {
    R <- ifelse(xor(A_yes, B_yes), "yes", "no"); RT <- exhaustive_rt
  } else {
    R <- ifelse(!A_yes & !B_yes, "NN", ifelse(A_yes & !B_yes, "AN", ifelse(!A_yes & B_yes, "NB", "AB")))
    RT <- exhaustive_rt
  }
  keep <- is.finite(RT)
  data.frame(R = R[keep], RT = RT[keep], A_yes = A_yes[keep], B_yes = B_yes[keep],
             Channel1_RT = A_rt[keep], Channel2_RT = B_rt[keep],
             SerialOrder = ifelse(A_first, "A", "B")[keep], stringsAsFactors = FALSE)
}

LogicalRules_rfun_lba <- function(n, pars, LogicalRule, shared_capacity = NULL,
                                   kappa = NULL, tau = NULL, rho = NULL, return_finish = FALSE) {
  if (!requireNamespace("rtdists", quietly = TRUE)) stop("LBA simulation requires rtdists.")
  racers <- c("A", "n_A", "B", "n_B")
  if (!all(racers %in% rownames(pars))) stop("pars must have rows A, n_A, B, n_B.")
  req <- c("A", "b", "t0", "v", "sv")
  if (!all(req %in% colnames(pars))) stop("pars must have columns A, b, t0, v, sv.")
  if (is.null(tau) && !is.null(rho)) tau <- rho
  shared <- .lba_resolve_shared_capacity(shared_capacity, kappa, tau)
  rts <- .lba_draw_rts(n, pars, racers, shared)
  if (is.null(dim(rts))) rts <- matrix(rts, nrow = n, ncol = length(racers))
  colnames(rts) <- racers
  if (isTRUE(return_finish)) return(rts)
  apply_logic_gates(LogicalRule, rts[, "A"], rts[, "n_A"], rts[, "B"], rts[, "n_B"])
}

.lba_resolve_shared_capacity <- function(shared_capacity = NULL, kappa = NULL, tau = NULL) {
  # A lower-level caller must opt in explicitly because `pars` alone does not
  # identify whether it represents an AB, AN, NB, or NN design cell.  The
  # higher-level simulator supplies active = TRUE only for AB cells.
  if (is.null(shared_capacity)) {
    if (is.null(kappa) && is.null(tau)) {
      return(list(kappa = 1, tau = 0, active = FALSE))
    }
    shared_capacity <- list()
  } else if (!is.list(shared_capacity)) {
    if (is.null(names(shared_capacity))) {
      stop("shared_capacity must be a list or a named vector.", call. = FALSE)
    }
    shared_capacity <- as.list(shared_capacity)
  }
  if (!is.null(kappa)) shared_capacity$kappa <- kappa
  if (!is.null(tau)) shared_capacity$tau <- tau

  get_one <- function(name, default) {
    value <- shared_capacity[[name]]
    if (is.null(value)) default else value
  }
  kappa <- get_one("kappa", 1)
  tau <- get_one("tau", get_one("rho", 0))
  active <- get_one("active", TRUE)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("shared-capacity kappa must be one finite value greater than zero.", call. = FALSE)
  }
  if (length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("shared-capacity tau must be one finite non-negative value.", call. = FALSE)
  }
  if (length(active) != 1L || is.na(active)) {
    stop("shared-capacity active must be one non-missing logical value.", call. = FALSE)
  }
  kappa <- as.numeric(kappa); tau <- as.numeric(tau)
  list(kappa = kappa, tau = tau,
       active = isTRUE(active) && (kappa != 1 || tau != 0))
}

.lba_shared_target_drifts <- function(n, pars, shared_capacity) {
  # This is the requested bivariate-normal construction.  The returned rates
  # are the raw normal draws; the LBA sampler clips negative rates to zero in
  # the same way as rtdists' internal LBA generator.
  z <- stats::rnorm(n)
  multiplier <- shared_capacity$kappa + shared_capacity$tau * z
  cbind(
    A = multiplier * pars["A", "v"] + stats::rnorm(n, mean = 0, sd = pars["A", "sv"]),
    B = multiplier * pars["B", "v"] + stats::rnorm(n, mean = 0, sd = pars["B", "sv"])
  )
}

.lba_finish_times <- function(n, pars, racers, drifts) {
  # Equivalent to the st0 = 0 portion of rtdists::rlba_norm, retaining a
  # separate finish time for each accumulator after supplying correlated
  # target drift draws.
  A <- as.numeric(pars[racers, "A"])
  b <- as.numeric(pars[racers, "b"])
  t0 <- as.numeric(pars[racers, "t0"])
  if (any(!is.finite(A)) || any(!is.finite(b)) || any(!is.finite(t0)) || any(b < A)) {
    stop("LBA parameters must be finite and satisfy b >= A.", call. = FALSE)
  }
  starts <- matrix(stats::runif(n * length(racers)), nrow = n,
                   ncol = length(racers))
  starts <- sweep(starts, 2L, A, `*`)
  drifts <- pmax(0, drifts)
  ttf <- sweep(-starts, 2L, b, `+`)
  ttf <- ttf / drifts
  out <- sweep(ttf, 2L, t0, `+`)
  colnames(out) <- racers
  out
}

.lba_draw_rts <- function(n, pars, racers = c("A", "n_A", "B", "n_B"),
                          shared_capacity = list(kappa = 1, tau = 0, active = FALSE)) {
  if (isTRUE(shared_capacity$active)) {
    target <- .lba_shared_target_drifts(n, pars, shared_capacity)
    target_finish <- .lba_finish_times(n, pars, c("A", "B"), target)
    rts <- cbind(
      A = target_finish[, "A"],
      n_A = rtdists::rlba_norm(n = n, A = pars["n_A", "A"], b = pars["n_A", "b"],
                               t0 = pars["n_A", "t0"], mean_v = pars["n_A", "v"],
                               sd_v = pars["n_A", "sv"], st0 = 0, posdrift = TRUE)[, "rt"],
      B = target_finish[, "B"],
      n_B = rtdists::rlba_norm(n = n, A = pars["n_B", "A"], b = pars["n_B", "b"],
                               t0 = pars["n_B", "t0"], mean_v = pars["n_B", "v"],
                               sd_v = pars["n_B", "sv"], st0 = 0, posdrift = TRUE)[, "rt"]
    )
    return(rts[, racers, drop = FALSE])
  }
  vapply(racers, function(i) {
    p <- pars[i, c("A", "b", "t0", "v", "sv"), drop = TRUE]
    rtdists::rlba_norm(n = n, A = p[["A"]], b = p[["b"]], t0 = p[["t0"]],
                       mean_v = p[["v"]], sd_v = p[["sv"]], st0 = 0, posdrift = TRUE)[, "rt"]
  }, numeric(n))
}

.lba_init_state <- function(p_vec, kappa = NULL, tau = NULL) {
  vc <- .param_get(p_vec, "vc_yes")
  b <- .param_get(p_vec, "b")
  get <- function(...) {
    nms <- c(...); hit <- nms[nms %in% names(p_vec)]
    if (length(hit)) .param_get(p_vec, hit[[1]]) else vc
  }
  kappa_value <- if (is.null(kappa)) {
    hit <- c("capacity_kappa", "capacity_target_kappa", "shared_capacity_kappa",
             "shared_kappa", "kappa")
    .param_get(p_vec, hit[hit %in% names(p_vec)][1L], 1)
  } else kappa
  tau_value <- if (is.null(tau)) {
    # `rho` is accepted as a backwards-compatible alias for `tau`.
    hit <- c("capacity_tau", "capacity_target_tau", "shared_capacity_tau",
             "shared_tau", "tau", "capacity_rho", "shared_capacity_rho",
             "shared_rho", "rho")
    .param_get(p_vec, hit[hit %in% names(p_vec)][1L], 0)
  } else tau
  list(vc_yesA_H = get("vc_yesA_H", "vc_yesA", "vc_yesH"),
       vc_yesA_L = get("vc_yesA_L", "vc_yesL", "vc_yesA"),
       vc_yesB_H = get("vc_yesB_H", "vc_yesB", "vc_yesH"),
       vc_yesB_L = get("vc_yesB_L", "vc_yesL", "vc_yesB"),
       b_yes = .param_get(p_vec, "b_yes", b), b_no = .param_get(p_vec, "b_no", b),
       cap_target = .param_get(p_vec, "capacity_target", .param_get(p_vec, "capacity", 0)),
       cap_absence = .param_get(p_vec, "capacity_absence", 0),
       shared_capacity = .lba_resolve_shared_capacity(
         list(kappa = kappa_value, tau = tau_value, active = TRUE)))
}

.lba_cell_pars <- function(p_vec, meta_row, state) {
  pars <- matrix(NA_real_, 4, 5, dimnames = list(c("A", "n_A", "B", "n_B"), c("A", "b", "t0", "v", "sv")))
  l1 <- as.integer(meta_row$Channel1); l2 <- as.integer(meta_row$Channel2)
  pars[, "A"] <- .param_get(p_vec, "A", .param_get(p_vec, "A_A", 1)); pars[, "t0"] <- .param_get(p_vec, "t0", .param_get(p_vec, "A_t0", 0))
  pars[c("A", "B"), "b"] <- state$b_yes; pars[c("n_A", "n_B"), "b"] <- state$b_no
  pars["A", "v"] <- if (l1 == 2) state$vc_yesA_H else if (l1 == 1) state$vc_yesA_L else .param_get(p_vec, "ve", 0)
  pars["B", "v"] <- if (l2 == 2) state$vc_yesB_H else if (l2 == 1) state$vc_yesB_L else .param_get(p_vec, "ve", 0)
  pars["n_A", "v"] <- if (l1 > 0) .param_get(p_vec, "ve", 0) else .param_get(p_vec, "vc_no", 0)
  pars["n_B", "v"] <- if (l2 > 0) .param_get(p_vec, "ve", 0) else .param_get(p_vec, "vc_no", 0)
  if (l1 > 0 && l2 > 0) pars[c("A", "B"), "v"] <- pars[c("A", "B"), "v"] + state$cap_target
  if (l1 == 0 && l2 == 0) pars[c("n_A", "n_B"), "v"] <- pars[c("n_A", "n_B"), "v"] + state$cap_absence
  pars["A", "sv"] <- if (l1 > 0) .param_get(p_vec, "sv_c", .param_get(p_vec, "A_sv", .1)) else .param_get(p_vec, "sv_e", .1)
  pars["n_A", "sv"] <- if (l1 > 0) .param_get(p_vec, "sv_e", .1) else .param_get(p_vec, "sv_c", .1)
  pars["B", "sv"] <- if (l2 > 0) .param_get(p_vec, "sv_c", .1) else .param_get(p_vec, "sv_e", .1)
  pars["n_B", "sv"] <- if (l2 > 0) .param_get(p_vec, "sv_e", .1) else .param_get(p_vec, "sv_c", .1)
  pars
}

.lba_simulate_rule_cell <- function(n, pars, logical_rule, shared_capacity = NULL) {
  LogicalRules_rfun_lba(n, pars, logical_rule, shared_capacity = shared_capacity)
}

.lba_simulate_serial_cell <- function(n, pars, logical_rule, serial_mode, serial_order, p_A_first) {
  rts <- LogicalRules_rfun_lba(n, pars, logical_rule, return_finish = TRUE)
  .serial_combine_races(rts, logical_rule, serial_mode, serial_order, p_A_first)
}

.ou_init_state <- function(p_vec) {
  if ("uc_yes" %in% names(p_vec) && !"uc_H" %in% names(p_vec)) p_vec["uc_H"] <- .param_get(p_vec, "uc_yes")
  if (!"uc_L" %in% names(p_vec)) p_vec["uc_L"] <- .param_get(p_vec, "uc_H")
  if (!"sv_e" %in% names(p_vec) && any(grepl("sv$", names(p_vec)))) p_vec["sv_e"] <- .param_get(p_vec, "sv")
  if (!"sv_c" %in% names(p_vec) && any(grepl("sv$", names(p_vec)))) p_vec["sv_c"] <- .param_get(p_vec, "sv")
  list(p_vec = p_vec, cap_target = .param_get(p_vec, "capacity_target", 0), cap_absence = .param_get(p_vec, "capacity_absence", 0))
}

.ou_cell_pars <- function(p_vec, meta_row, state) {
  pv <- state$p_vec; pars <- matrix(NA_real_, 4, 5, dimnames = list(c("A", "n_A", "B", "n_B"), c("a", "u", "sv", "b", "t0")))
  pars[, "a"] <- .param_get(pv, "a", -1); pars[, "b"] <- .param_get(pv, "b", 1); pars[, "t0"] <- .param_get(pv, "t0", 0)
  set_channel <- function(ch, level) {
    if (level == 2) c(u = .param_get(pv, "uc_H"), sv = .param_get(pv, "sv_c"))
    else if (level == 1) c(u = .param_get(pv, "uc_L"), sv = .param_get(pv, "sv_c"))
    else c(u = .param_get(pv, "ue_yes"), sv = .param_get(pv, "sv_e"))
  }
  l1 <- as.integer(meta_row$Channel1); l2 <- as.integer(meta_row$Channel2)
  a <- set_channel("A", l1); b <- set_channel("B", l2)
  pars["A", c("u", "sv")] <- a; pars["B", c("u", "sv")] <- b
  pars["n_A", "u"] <- if (l1 > 0) .param_get(pv, "ue_no") else .param_get(pv, "uc_no")
  pars["n_B", "u"] <- if (l2 > 0) .param_get(pv, "ue_no") else .param_get(pv, "uc_no")
  pars["n_A", "sv"] <- if (l1 > 0) .param_get(pv, "sv_e") else .param_get(pv, "sv_c")
  pars["n_B", "sv"] <- if (l2 > 0) .param_get(pv, "sv_e") else .param_get(pv, "sv_c")
  if (l1 > 0 && l2 > 0) { pars[c("A", "B"), "u"] <- pars[c("A", "B"), "u"] + state$cap_target * pars[c("A", "B"), "u"] }
  if (l1 == 0 && l2 == 0) { pars[c("n_A", "n_B"), "u"] <- pars[c("n_A", "n_B"), "u"] + state$cap_absence * pars[c("n_A", "n_B"), "u"] }
  pars
}

LogicalRules_OU_rfun_R <- function(n, pars, LogicalRule, cross_talk = c(pos = 0, neg = 0), dt = .001, max_time = 3,
                                   bridge = TRUE, bridge_exp_cutoff = -20, return_finish = FALSE) {
  racers <- c("A", "n_A", "B", "n_B"); a <- pars[racers, "a"]; u <- pars[racers, "u"]; sv <- pars[racers, "sv"]; b <- pars[racers, "b"]; t0 <- pars[racers, "t0"]
  threshold <- b; threshold[c(2, 4)] <- -threshold[c(2, 4)]
  x <- matrix(0, n, 4); finish <- matrix(Inf, n, 4); done <- matrix(FALSE, n, 4); colnames(finish) <- racers
  A_res <- B_res <- rep(FALSE, n); active <- seq_len(n); steps <- ceiling(max_time / dt)
  for (step in seq_len(steps)) {
    if (!length(active)) break
    keep <- !(A_res[active] & B_res[active]); active <- active[keep]; if (!length(active)) break
    old <- x[active, , drop = FALSE]; z <- matrix(stats::rnorm(length(active) * 4), nrow = length(active), ncol = 4)
    drift <- old; drift[, 1] <- a[1] * old[, 1] + cross_talk[["pos"]] * old[, 3] + u[1]; drift[, 3] <- cross_talk[["pos"]] * old[, 1] + a[3] * old[, 3] + u[3]
    drift[, 2] <- a[2] * old[, 2] + cross_talk[["neg"]] * old[, 4] + u[2]; drift[, 4] <- cross_talk[["neg"]] * old[, 2] + a[4] * old[, 4] + u[4]
    now <- old + drift * dt + z * rep(sv, each = nrow(old)) * sqrt(dt); x[active, ] <- now
    for (j in 1:4) {
      pos <- which(!done[active, j]); if (!length(pos)) next
      neg <- j %in% c(2, 4); hit <- if (neg) now[pos, j] <= threshold[j] else now[pos, j] >= threshold[j]
      bridge_hit <- rep(FALSE, length(pos))
      if (bridge) {
        candidates <- which(!hit); if (length(candidates)) {
          d0 <- threshold[j] - old[pos[candidates], j]; d1 <- threshold[j] - now[pos[candidates], j]
          prob <- exp(pmax(bridge_exp_cutoff, -2 * d0 * d1 / max(sv[j]^2 * dt, 1e-12)))
          bridge_hit[candidates] <- stats::runif(length(candidates)) < prob; hit[candidates] <- bridge_hit[candidates]
        }
      }
      if (any(hit)) {
        hit_time <- rep(step * dt, length(pos)); ii <- which(hit & !bridge_hit)
        if (length(ii)) { frac <- (threshold[j] - old[pos[ii], j]) / (now[pos[ii], j] - old[pos[ii], j]); frac[!is.finite(frac)] <- 1; hit_time[ii] <- (step - 1) * dt + pmin(1, pmax(0, frac)) * dt }
        ii <- which(hit & bridge_hit); if (length(ii)) hit_time[ii] <- (step - 1) * dt + stats::runif(length(ii)) * dt
        rows <- active[pos[hit]]; finish[rows, j] <- hit_time[hit] + t0[j]; done[rows, j] <- TRUE
        if (j <= 2) A_res[rows] <- TRUE else B_res[rows] <- TRUE
      }
    }
  }
  if (isTRUE(return_finish)) return(finish)
  apply_logic_gates(LogicalRule, finish[, 1], finish[, 2], finish[, 3], finish[, 4])
}

.sft_ou_env <- new.env(parent = baseenv())
.ensure_ou_rcpp <- local({
  compiled <- FALSE
  function() {
    if (compiled && exists("LogicalRules_OU_finish_rcpp", envir = .sft_ou_env, mode = "function")) return(TRUE)
    if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
    ok <- tryCatch({
      Rcpp::sourceCpp(code = '
        #include <Rcpp.h>
        #include <algorithm>
        #include <cmath>
        using namespace Rcpp;
        // [[Rcpp::export]]
        NumericMatrix LogicalRules_OU_finish_rcpp(int n, NumericVector a, NumericVector u,
          NumericVector sv, NumericVector b, NumericVector t0, double cp, double cn,
          double dt, double max_time, bool bridge, double bridge_cut, bool adaptive,
          double adapt_factor, double adapt_eps) {
          NumericMatrix out(n, 4); std::fill(out.begin(), out.end(), R_PosInf);
          double dmin = std::max(dt / 32.0, 1e-8), dmax = std::max(dt, dt * adapt_factor);
          for (int tr = 0; tr < n; ++tr) {
            double x[4] = {0,0,0,0}; double time = 0; bool done[4] = {false,false,false,false}; bool ar = false, br = false;
            while (time < max_time && !(ar && br)) {
              double h = std::min(dt, max_time - time);
              if (adaptive) { double md = R_PosInf; for (int j=0;j<4;++j) if (!done[j]) md = std::min(md, std::fabs((j==1||j==3 ? -b[j] : b[j])-x[j])); double smax = std::max(sv[0], std::max(sv[1], std::max(sv[2], sv[3]))); h = std::min(dmax, std::max(dmin, dt * std::pow(std::max(md / std::max(adapt_eps * smax, 1e-8), 1.0), 2.0))); h = std::min(h, max_time-time); }
              double old[4] = {x[0],x[1],x[2],x[3]}; double sq = std::sqrt(h);
              x[0] += (a[0]*old[0] + cp*old[2] + u[0])*h + sv[0]*sq*R::rnorm(0,1);
              x[1] += (a[1]*old[1] + cn*old[3] + u[1])*h + sv[1]*sq*R::rnorm(0,1);
              x[2] += (cp*old[0] + a[2]*old[2] + u[2])*h + sv[2]*sq*R::rnorm(0,1);
              x[3] += (cn*old[1] + a[3]*old[3] + u[3])*h + sv[3]*sq*R::rnorm(0,1);
              for (int j=0;j<4;++j) if (!done[j]) { bool neg=(j==1||j==3); double th=neg?-b[j]:b[j]; bool str=neg?x[j]<=th:x[j]>=th; bool hit=str; bool bridge_hit=false; if (bridge && !str) { double d0=th-old[j], d1=th-x[j]; double p=std::exp(std::max(bridge_cut, -2*d0*d1/std::max(sv[j]*sv[j]*h,1e-12))); if (R::runif(0,1)<p) {hit=true; bridge_hit=true;} } if (hit) { double ht=time+h; if (str) {double f=(th-old[j])/(x[j]-old[j]); if (!R_finite(f)) f=1; ht=time+std::min(1.0,std::max(0.0,f))*h;} else if (bridge_hit) ht=time+R::runif(0,1)*h; out(tr,j)=ht+t0[j]; done[j]=true; if (j<2) ar=true; else br=true; } }
              time += h;
            }
          }
          return out;
        }', env = .sft_ou_env)
      TRUE
    }, error = function(e) {
      assign("compile_error", conditionMessage(e), envir = .sft_ou_env)
      FALSE
    })
    compiled <<- isTRUE(ok) && exists("LogicalRules_OU_finish_rcpp", envir = .sft_ou_env, mode = "function")
    compiled
  }
})

LogicalRules_OU_rfun_Rcpp <- function(n, pars, LogicalRule, cross_talk = c(pos = 0, neg = 0), dt = .001, max_time = 3,
                                      bridge = TRUE, bridge_exp_cutoff = -20, adaptive = FALSE,
                                      adaptive_mode = c("distance", "ptol"), adapt_factor = 32, adapt_eps = 1.5,
                                      p_tol = 1e-9, curv_tol = .1, return_finish = FALSE) {
  if (!.ensure_ou_rcpp()) stop("Rcpp backend not available; use backend='R'.")
  adaptive_mode <- match.arg(adaptive_mode); racers <- c("A", "n_A", "B", "n_B")
  finish <- .sft_ou_env$LogicalRules_OU_finish_rcpp(as.integer(n), as.numeric(pars[racers, "a"]), as.numeric(pars[racers, "u"]),
                                                    as.numeric(pars[racers, "sv"]), as.numeric(pars[racers, "b"]), as.numeric(pars[racers, "t0"]),
                                                    as.numeric(cross_talk[["pos"]]), as.numeric(cross_talk[["neg"]]), as.numeric(dt), as.numeric(max_time),
                                                    isTRUE(bridge), as.numeric(bridge_exp_cutoff), isTRUE(adaptive), as.numeric(adapt_factor), as.numeric(adapt_eps))
  colnames(finish) <- racers
  if (isTRUE(return_finish)) return(finish)
  apply_logic_gates(LogicalRule, finish[, 1], finish[, 2], finish[, 3], finish[, 4])
}

LogicalRules_OU_rfun <- function(n, pars, LogicalRule, cross_talk = c(pos = 0, neg = 0), dt = .001, max_time = 3,
                                 bridge = TRUE, bridge_exp_cutoff = -20, backend = c("auto", "Rcpp", "R"),
                                 adaptive = FALSE, adaptive_mode = c("distance", "ptol"), adapt_factor = 32,
                                 adapt_eps = 1.5, p_tol = 1e-9, curv_tol = .1, return_finish = FALSE) {
  backend <- match.arg(backend); adaptive_mode <- match.arg(adaptive_mode)
  use <- if (backend == "Rcpp") TRUE else if (backend == "R") FALSE else .ensure_ou_rcpp()
  if (use) LogicalRules_OU_rfun_Rcpp(n, pars, LogicalRule, cross_talk, dt, max_time, bridge, bridge_exp_cutoff,
                                      adaptive, adaptive_mode, adapt_factor, adapt_eps, p_tol, curv_tol,
                                      return_finish = return_finish)
  else { if (adaptive) warning("adaptive=TRUE is only available in the Rcpp backend; using fixed dt.");
    LogicalRules_OU_rfun_R(n, pars, LogicalRule, cross_talk, dt, max_time, bridge, bridge_exp_cutoff,
                           return_finish = return_finish) }
}

.ou_simulate_rule_cell <- function(n, pars, logical_rule, p_vec, dt = .001, max_time = 2,
                                   bridge = TRUE, bridge_exp_cutoff = -20, backend = "auto", adaptive = FALSE,
                                   adaptive_mode = "distance", adapt_factor = 32, adapt_eps = 1.5,
                                   p_tol = 1e-9, curv_tol = .1) {
  LogicalRules_OU_rfun(n, pars, logical_rule,
                       cross_talk = c(pos = .param_get(p_vec, "cross_talk_pos", 0), neg = .param_get(p_vec, "cross_talk_neg", 0)),
                       dt = dt, max_time = max_time, bridge = bridge, bridge_exp_cutoff = bridge_exp_cutoff,
                       backend = backend, adaptive = adaptive, adaptive_mode = adaptive_mode,
                       adapt_factor = adapt_factor, adapt_eps = adapt_eps, p_tol = p_tol, curv_tol = curv_tol)
}

.ou_simulate_serial_cell <- function(n, pars, logical_rule, serial_mode, serial_order, p_A_first,
                                     dt = .001, max_time = 2, bridge = TRUE, bridge_exp_cutoff = -20,
                                     backend = "auto", adaptive = FALSE, adaptive_mode = "distance",
                                     adapt_factor = 32, adapt_eps = 1.5, p_tol = 1e-9, curv_tol = .1) {
  finish <- LogicalRules_OU_rfun(n, pars, logical_rule, cross_talk = c(pos = 0, neg = 0),
                                  dt = dt, max_time = max_time, bridge = bridge,
                                  bridge_exp_cutoff = bridge_exp_cutoff, backend = backend,
                                  adaptive = adaptive, adaptive_mode = adaptive_mode,
                                  adapt_factor = adapt_factor, adapt_eps = adapt_eps,
                                  p_tol = p_tol, curv_tol = curv_tol, return_finish = TRUE)
  .serial_combine_races(finish, logical_rule, serial_mode, serial_order, p_A_first)
}

.compute_capacity <- function(data, logical_rules) {
  out <- list(); by <- split(data, data$Condition)
  if ("OR" %in% logical_rules && !is.null(by$OR)) {
    d <- by$OR; rt <- lapply(c("AB", "AN", "NB"), function(s) d$RT[d$S == s]); cr <- lapply(c("AB", "AN", "NB"), function(s) d$Correct[d$S == s])
    if (all(lengths(rt) > 0L)) out$OR_Ct <- capacity.or(rt, cr)
  }
  if ("AND" %in% logical_rules && !is.null(by$AND)) {
    d <- by$AND; rt <- lapply(c("AB", "AN", "NB"), function(s) d$RT[d$S == s]); cr <- lapply(c("AB", "AN", "NB"), function(s) d$Correct[d$S == s])
    if (all(lengths(rt) > 0L)) out$AND_Ct <- capacity.and(rt, cr)
    absence_rt <- list(d$RT[d$S == "NN"], rt[[2]], rt[[3]]); absence_cr <- list(d$Correct[d$S == "NN"], cr[[2]], cr[[3]])
    if (all(lengths(absence_rt) > 0L)) out$AND_Absence_Ct <- capacity.or(absence_rt, absence_cr)
    alt_rt <- c(rt, list(d$RT[d$S == "NN"])); alt_cr <- c(cr, list(d$Correct[d$S == "NN"]))
    if (all(lengths(alt_rt) > 0L)) out$Altieri_Ct <- capacity.altieri(alt_rt, alt_cr)
  }
  if ("ID" %in% logical_rules && !is.null(by$ID)) {
    d <- by$ID; id_rt <- list(d$RT[d$S == "AB"], d$RT[d$S == "NN"], d$RT[d$S == "AN"], d$RT[d$S == "NB"])
    if (all(lengths(id_rt) > 0L)) out$ID_Ct <- capacity.id(id_rt[[1]], id_rt[[2]], id_rt[3:4], d$Correct[d$S == "AB"], d$Correct[d$S == "NN"], list(d$Correct[d$S == "AN"], d$Correct[d$S == "NB"]))
  }
  out
}

.safe_sic <- function(d) {
  if (!is.data.frame(d) || !nrow(d)) return(list())
  cells <- list(d$RT[d$Channel1 == 2 & d$Channel2 == 2 & d$Correct == 1], d$RT[d$Channel1 == 2 & d$Channel2 == 1 & d$Correct == 1], d$RT[d$Channel1 == 1 & d$Channel2 == 2 & d$Correct == 1], d$RT[d$Channel1 == 1 & d$Channel2 == 1 & d$Correct == 1])
  if (any(!lengths(cells))) return(list())
  sic(cells[[1]], cells[[2]], cells[[3]], cells[[4]])
}

.compute_sics <- function(data, logical_rules) {
  by <- split(data, data$Condition); out <- list()
  for (r in logical_rules) if (!is.null(by[[r]])) { z <- .safe_sic(by[[r]]); if (length(z)) out[[paste0(r, "_SIC")]] <- z }
  out
}

.compute_cdfs <- function(data, logical_rules) {
  by <- split(data, data$Condition); out <- list(); cells <- c("AB", "AN", "NB", "NN")
  for (r in logical_rules) if (!is.null(by[[r]])) { d <- split(by[[r]], by[[r]]$S); out[[r]] <- lapply(intersect(cells, names(d)), function(s) .sft_ecdf(d[[s]]$RT[d[[s]]$Correct == 1])); names(out[[r]]) <- intersect(cells, names(d)) }
  out
}

.compute_at <- function(data, logical_rules) {
  by <- split(data, data$Condition); out <- list()
  for (r in intersect(c("OR", "AND"), logical_rules)) if (!is.null(by[[r]])) { d <- split(by[[r]], by[[r]]$S); if (all(c("AB", "AN", "NB") %in% names(d))) out[[paste0(r, "_At")]] <- if (r == "OR") assessment_ta_or(lapply(c("AB", "AN", "NB"), function(s) d[[s]]$RT), lapply(c("AB", "AN", "NB"), function(s) d[[s]]$Correct)) else assessment_ta_and(lapply(c("AB", "AN", "NB"), function(s) d[[s]]$RT), lapply(c("AB", "AN", "NB"), function(s) d[[s]]$Correct)) }
  out
}

.resolve_serial_options <- function(architecture, serial_mode, serial_order, p_A_first,
                                    serial_first = NULL, serial_stopping = NULL,
                                    serial_p_A_first = NULL, serial_type = NULL) {
  architecture <- match.arg(architecture, c("parallel", "serial"))
  if (!is.null(serial_stopping)) serial_mode <- serial_stopping
  if (!is.null(serial_type)) serial_mode <- serial_type
  serial_mode <- sub("_", "-", as.character(serial_mode)[1L])
  serial_mode <- match.arg(serial_mode, c("self-terminating", "exhaustive"))
  serial_order <- match.arg(serial_order, c("random", "A", "B", "fixed"))
  if (!is.null(serial_p_A_first)) p_A_first <- serial_p_A_first
  if (!is.null(serial_first)) {
    serial_first <- match.arg(as.character(serial_first)[1L], c("A", "B"))
    if (serial_order %in% c("random", "fixed")) serial_order <- serial_first
    else if (!identical(serial_order, serial_first)) stop("serial_first conflicts with serial_order.", call. = FALSE)
  }
  if (serial_order == "fixed") serial_order <- "A"
  if (!is.numeric(p_A_first) || length(p_A_first) != 1L || is.na(p_A_first) || !is.finite(p_A_first) || p_A_first < 0 || p_A_first > 1) {
    stop("p_A_first must be one finite value in [0, 1].", call. = FALSE)
  }
  list(architecture = architecture, serial_mode = serial_mode,
       serial_order = serial_order, p_A_first = as.numeric(p_A_first))
}

simulate_sft <- function(model = c("lba", "ou"), n, p_vec, design = NULL, salience = FALSE,
                          logical_rules = NULL, sics = FALSE, cdfs = FALSE, a_t = FALSE,
                          subject = 1, backend = c("auto", "Rcpp", "R"), dt = .001, max_time = 2,
                          bridge = TRUE, bridge_exp_cutoff = -20, adaptive = TRUE,
                          adaptive_mode = c("ptol", "distance"), adapt_factor = 32, adapt_eps = 1.5,
                          p_tol = 1e-8, curv_tol = .05, kappa = NULL, tau = NULL, rho = NULL,
                          architecture = c("parallel", "serial"),
                          serial_mode = c("self-terminating", "exhaustive"),
                          serial_order = c("random", "A", "B", "fixed"), p_A_first = .5,
                          serial_first = NULL, serial_stopping = NULL, serial_p_A_first = NULL,
                          serial_type = NULL) {
  model <- match.arg(model); backend <- match.arg(backend); adaptive_mode <- match.arg(adaptive_mode)
  if (length(n) != 1L || n < 1 || n != as.integer(n)) stop("n must be a positive integer.")
  serial <- .resolve_serial_options(architecture, serial_mode, serial_order, p_A_first,
                                     serial_first, serial_stopping, serial_p_A_first, serial_type)
  architecture <- serial$architecture; serial_mode <- serial$serial_mode
  serial_order <- serial$serial_order; p_A_first <- serial$p_A_first
  if (is.null(tau) && !is.null(rho)) tau <- rho
  if (is.null(logical_rules)) {
    # ID is not defined for serial self-termination. Keep the default focused
    # on the requested OR case; AND remains available when explicitly asked
    # for and uses the analogous first-absence rule.
    logical_rules <- if (architecture == "serial" && serial_mode == "self-terminating") c("OR") else c("OR", "AND", "ID")
  }
  logical_rules <- .normalize_logical_rules(logical_rules); design_df <- .expand_design(design, salience)
  if (architecture == "serial" && serial_mode == "self-terminating" && "ID" %in% logical_rules) {
    stop("Serial self-terminating simulation currently supports OR and AND, not ID.", call. = FALSE)
  }
  sim_p_vec <- p_vec; init_p_vec <- p_vec
  if (architecture == "serial" && model == "lba") {
    # Ignore shared-capacity values before validation as well as after it.
    if ("kappa" %in% names(init_p_vec)) init_p_vec[["kappa"]] <- 1
    for (nm in c("tau", "rho")) if (nm %in% names(init_p_vec)) init_p_vec[[nm]] <- 0
    state <- .lba_init_state(init_p_vec, kappa = 1, tau = 0)
  } else state <- if (model == "lba") .lba_init_state(p_vec, kappa = kappa, tau = tau) else .ou_init_state(p_vec)
  if (architecture == "serial") {
    state$cap_target <- 0; state$cap_absence <- 0
    if (model == "lba") state$shared_capacity <- .lba_resolve_shared_capacity(list(kappa = 1, tau = 0, active = FALSE))
    else {
      for (nm in c("cross_talk_pos", "cross_talk_neg")) if (nm %in% names(sim_p_vec)) sim_p_vec[[nm]] <- 0
      state$p_vec <- sim_p_vec
    }
  }
  chunks <- lapply(logical_rules, function(rule) lapply(seq_len(nrow(design_df)), function(i) {
    meta <- design_df[i, , drop = FALSE]
    pars <- if (model == "lba") .lba_cell_pars(p_vec, meta, state) else .ou_cell_pars(p_vec, meta, state)
    shared <- if (architecture == "parallel" && model == "lba" && meta$Channel1 > 0 && meta$Channel2 > 0) state$shared_capacity else NULL
    sim <- if (architecture == "serial") {
      if (model == "lba") .lba_simulate_serial_cell(n, pars, rule, serial_mode, serial_order, p_A_first)
      else .ou_simulate_serial_cell(n, pars, rule, serial_mode, serial_order, p_A_first, dt, max_time,
                                    bridge, bridge_exp_cutoff, backend, adaptive, adaptive_mode,
                                    adapt_factor, adapt_eps, p_tol, curv_tol)
    } else if (model == "lba") .lba_simulate_rule_cell(n, pars, rule, shared)
    else .ou_simulate_rule_cell(n, pars, rule, p_vec, dt, max_time, bridge, bridge_exp_cutoff,
                                backend, adaptive, adaptive_mode, adapt_factor, adapt_eps, p_tol, curv_tol)
    .annotate_trials(sim, rule, meta)
  }))
  data <- .sft_bind_rows(unlist(chunks, recursive = FALSE)); data$Correct <- as.numeric(correctfun(data)); data$Subject <- rep(subject, length.out = nrow(data)); data$Model <- model; data$Backend <- if (model == "ou") backend else "R"
  if (architecture == "serial") {
    data$Architecture <- rep(architecture, length.out = nrow(data))
    data$SerialMode <- rep(serial_mode, length.out = nrow(data))
    if (!"SerialOrder" %in% names(data)) data$SerialOrder <- character(nrow(data))
    data$p_A_first <- rep(p_A_first, length.out = nrow(data))
  }
  by_rule <- .split_by_rule(data, logical_rules); out <- list(data = data, by_rule = by_rule, metrics = .compute_capacity(data, logical_rules))
  if (sics) out$sics <- .compute_sics(data, logical_rules); if (cdfs) out$cdfs <- .compute_cdfs(data, logical_rules); if (a_t) out$a_t <- .compute_at(data, logical_rules)
  out
}

# ---- Factorial design expansion, logical-rule labelling, and LBA search helpers ----

.normalize_logical_rules <- function(logical_rules) {
  if (is.null(logical_rules)) logical_rules <- c("OR", "AND", "ID")
  logical_rules <- unique(as.character(logical_rules))
  bad <- setdiff(logical_rules, c("OR", "AND", "XOR", "ID"))
  if (length(bad)) stop("Unsupported logical rules: ", paste(bad, collapse = ", "))
  logical_rules
}


.default_design <- function(salience = FALSE) {
  if (salience) c("HH", "HL", "LH", "LL", "HN", "LN", "NH", "NL", "NN") else c("AB", "AN", "NB", "NN")
}


.char_to_level <- function(ch) {
  if (ch %in% c("A", "B", "H")) return(2L)
  if (ch == "L") return(1L)
  if (ch %in% c("N", "_")) return(0L)
  stop("Unknown design character '", ch, "'.")
}


.expand_design <- function(design = NULL, salience = FALSE) {
  if (is.null(design)) design <- .default_design(salience)
  .sft_bind_rows(lapply(as.character(design), function(cell) {
    if (nchar(cell) != 2L) stop("Design cells must contain two characters.")
    l1 <- .char_to_level(substr(cell, 1, 1)); l2 <- .char_to_level(substr(cell, 2, 2))
    s <- if (l1 > 0 && l2 > 0) "AB" else if (l1 > 0) "AN" else if (l2 > 0) "NB" else "NN"
    data.frame(cell_raw = cell, S = s, Channel1 = l1, Channel2 = l2,
               A_present = l1 > 0, B_present = l2 > 0, stringsAsFactors = FALSE)
  }))
}


.annotate_trials <- function(tmp, logical_rule, meta_row) {
  n <- nrow(tmp)
  tmp$Condition <- rep(logical_rule, n); tmp$S <- rep(meta_row$S, n)
  tmp$Channel1 <- rep(meta_row$Channel1, n); tmp$Channel2 <- rep(meta_row$Channel2, n)
  tmp$Channel1_acc <- as.numeric((tmp$Channel1 > 0) == tmp$A_yes)
  tmp$Channel2_acc <- as.numeric((tmp$Channel2 > 0) == tmp$B_yes)
  tmp
}


.split_by_rule <- function(data, logical_rules) {
  parts <- split(data, data$Condition)
  setNames(lapply(logical_rules, function(r) if (is.null(parts[[r]])) data[0, , drop = FALSE] else parts[[r]]), logical_rules)
}


obj_max_vs_min <- function(vy, vno, t, A, b, t0, sv) {
  if (!requireNamespace("rtdists", quietly = TRUE)) stop("This helper requires rtdists.")
  fy <- rtdists::plba_norm(t, A = A, b = b, t0 = t0, mean_v = vy, sd_v = sv, posdrift = TRUE)
  fn <- rtdists::plba_norm(t, A = A, b = b, t0 = t0, mean_v = vno, sd_v = sv, posdrift = TRUE)
  mean((fy^2 - (1 - (1 - fn)^2))^2)
}


find_equal_vc_yes <- function(vno, t, A, b, t0, sv, interval = c(1, 4.5)) {
  stats::optimize(obj_max_vs_min, interval = interval, vno = vno, t = t, A = A, b = b, t0 = t0, sv = sv)$minimum
}
