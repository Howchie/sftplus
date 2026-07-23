# Core capacity, UCIP, SIC, resilience, Bayesian, plotting, and assessment
# implementations for sft.plus.  Public functions in this file replace the
# corresponding upstream implementations; each public name has one canonical
# definition in the package.

.sft_zero_curve <- function(value = 0) {
  force(value)
  function(t) rep(value, length(t))
}

.sft_curve <- function(x, y, default = NA_real_) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y) & !is.na(y)
  if (!any(keep)) return(.sft_zero_curve(default))
  x <- x[keep]
  y <- y[keep]
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (length(x) == 1L) return(.sft_zero_curve(y[[1L]]))
  stats::approxfun(x, y, rule = 2, ties = "ordered")
}

.sft_finite_times <- function(RT) {
  x <- as.numeric(unlist(RT, recursive = TRUE, use.names = FALSE))
  x <- sort(unique(x[is.finite(x)]))
  if (!length(x)) stop("RT contains no finite response times.")
  x
}

.sft_correct <- function(CR, n) {
  if (is.null(CR) || length(CR) != n) return(rep(TRUE, n))
  out <- as.logical(CR)
  out[is.na(out)] <- FALSE
  out
}

.sft_clean <- function(RT, CR = NULL) {
  RT <- as.numeric(RT)
  CR <- .sft_correct(CR, length(RT))
  keep <- is.finite(RT)
  list(RT = RT[keep], CR = CR[keep])
}

.sft_cr_list <- function(RT, CR = NULL) {
  if (is.null(CR) || length(CR) != length(RT)) {
    return(lapply(RT, function(x) rep(TRUE, length(x))))
  }
  lapply(seq_along(RT), function(i) .sft_correct(CR[[i]], length(RT[[i]])))
}

.sft_data_numeric <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) stop(name, " must contain numeric values.", call. = FALSE)
  out
}

.sft_data_correct <- function(x) {
  if (is.logical(x)) {
    out <- x
  } else {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      key <- tolower(trimws(x))
      bad <- !is.na(key) & !key %in% c("0", "1", "true", "false")
      if (any(bad)) stop("Correct must contain only logical or 0/1 values.",
                         call. = FALSE)
      out <- key %in% c("1", "true")
      out[is.na(key)] <- NA
    } else {
      y <- suppressWarnings(as.numeric(x))
      bad <- is.na(y) & !is.na(x)
      if (any(bad) || any(!is.na(y) & !y %in% c(0, 1))) {
        stop("Correct must contain only logical or 0/1 values.", call. = FALSE)
      }
      out <- y == 1
    }
  }
  out[is.na(out)] <- FALSE
  out
}

#' Convert canonical row-wise SFT data to RT/CR list input.
#'
#' The list-based SFT estimators historically accept one response-time vector
#' per experimental cell.  This adapter extracts the usual two-channel cells
#' from a trial-level data frame.  For OR and AND it returns AB, A, and B,
#' using positive channel values as target-present.  For STST it returns
#' context and target cells using the package's signed-channel convention.
#' With multiple subjects, the default result is nested by subject so it can
#' be passed directly to `ucip.bayes.test()`; select one subject or use
#' `by_subject = "never"` for single-subject functions.
#'
#' @param data A data frame containing `Subject`, `RT`, `Correct`, `Channel1`,
#'   and `Channel2`; `Condition` is optional when only one condition is used.
#' @param Condition Optional single condition value to retain.  If omitted,
#'   the data must contain one condition (or no `Condition` column).
#' @param Subject Optional subject value or values to retain.
#' @param stopping.rule One of `"OR"`, `"AND"`, or `"STST"`.
#' @param by_subject `"auto"` returns a flat list for one subject and nested
#'   lists for multiple subjects; `"always"` or `"never"` forces either form.
#' @param include_nn Include the both-off `NN` cell for OR/AND conversion.
#' @return A list with `RT` and matching `CR` lists.  Metadata are stored as
#'   attributes, including selected subjects, condition, and dropped rows.
#' @export
sft_data_to_rt <- function(data, Condition = NULL, Subject = NULL,
                           stopping.rule = c("OR", "AND", "STST"),
                           by_subject = c("auto", "always", "never"),
                           include_nn = FALSE) {
  if (!is.data.frame(data)) stop("data must be a data.frame.", call. = FALSE)
  rule <- match.arg(stopping.rule)
  by_subject <- match.arg(by_subject)
  required <- c("Subject", "RT", "Correct", "Channel1", "Channel2")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("data is missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!nrow(data)) stop("data must contain at least one trial.", call. = FALSE)

  subject <- as.character(data$Subject)
  if (anyNA(subject) || any(!nzchar(subject))) {
    stop("Subject must not contain missing or empty values.", call. = FALSE)
  }
  if ("Condition" %in% names(data)) {
    condition <- as.character(data$Condition)
    if (anyNA(condition)) stop("Condition must not contain missing values.", call. = FALSE)
  } else {
    condition <- rep("default", nrow(data))
  }
  rt <- .sft_data_numeric(data$RT, "RT")
  correct <- .sft_data_correct(data$Correct)
  channel1 <- .sft_data_numeric(data$Channel1, "Channel1")
  channel2 <- .sft_data_numeric(data$Channel2, "Channel2")
  if (any(!is.finite(channel1)) || any(!is.finite(channel2))) {
    stop("Channel1 and Channel2 must contain finite numeric values.", call. = FALSE)
  }
  if (rule != "STST" && (any(channel1 < 0) || any(channel2 < 0))) {
    stop("OR/AND conversion requires non-negative Channel1 and Channel2 values.",
         call. = FALSE)
  }

  if (is.null(Condition)) {
    selected_condition <- unique(condition)
    if (length(selected_condition) != 1L) {
      stop("data contains multiple Condition values; supply Condition to select one.",
           call. = FALSE)
    }
  } else {
    selected_condition <- as.character(Condition)
    if (length(selected_condition) != 1L || is.na(selected_condition) ||
        !selected_condition %in% condition) {
      stop("Condition must select exactly one value present in data.", call. = FALSE)
    }
  }
  keep <- condition == selected_condition
  if (!any(keep)) stop("No rows remain after Condition selection.", call. = FALSE)

  if (is.null(Subject)) {
    selected_subject <- unique(subject[keep])
  } else {
    selected_subject <- as.character(Subject)
    if (!length(selected_subject) || anyNA(selected_subject) ||
        any(!selected_subject %in% subject[keep])) {
      stop("Subject must select values present in the selected condition.", call. = FALSE)
    }
    selected_subject <- unique(selected_subject)
  }
  if (by_subject == "never" && length(selected_subject) != 1L) {
    stop("A single Subject is required when by_subject = 'never'.", call. = FALSE)
  }
  if (by_subject == "auto") by_subject <- if (length(selected_subject) > 1L) "always" else "never"

  positive1 <- channel1 > 0
  positive2 <- channel2 > 0
  if (rule == "STST") {
    context <- positive1 + positive2 == 1 &
      (channel1 < 0 | channel2 < 0)
    target <- channel1 >= 0 & channel2 >= 0 & (channel1 != 0) + (channel2 != 0) == 1
    cell <- ifelse(context, "context", ifelse(target, "target", NA_character_))
    cells <- c("context", "target")
  } else {
    cell <- ifelse(positive1 & positive2, "AB",
                   ifelse(positive1, "A",
                          ifelse(positive2, "B", if (include_nn) "NN" else NA_character_)))
    cells <- c("AB", "A", "B")
    if (isTRUE(include_nn)) cells <- c(cells, "NN")
  }
  if (!any(!is.na(cell) & keep)) {
    stop("No usable channel cells remain after conversion.", call. = FALSE)
  }

  make_one <- function(id) {
    take <- keep & subject == id & !is.na(cell)
    rt_out <- setNames(lapply(cells, function(label) rt[take & cell == label]), cells)
    cr_out <- setNames(lapply(cells, function(label) correct[take & cell == label]), cells)
    list(RT = rt_out, CR = cr_out)
  }
  parts <- lapply(selected_subject, make_one)
  if (by_subject == "always") {
    rt_out <- lapply(parts, `[[`, "RT")
    cr_out <- lapply(parts, `[[`, "CR")
    names(rt_out) <- names(cr_out) <- selected_subject
  } else {
    rt_out <- parts[[1L]]$RT
    cr_out <- parts[[1L]]$CR
  }
  out <- list(RT = rt_out, CR = cr_out)
  attr(out, "subjects") <- selected_subject
  attr(out, "Condition") <- selected_condition
  attr(out, "stopping.rule") <- rule
  attr(out, "cells") <- cells
  attr(out, "dropped_rows") <- sum(keep & is.na(cell))
  class(out) <- c("sft_rt_data", "list")
  out
}

.sft_as_rt_cr <- function(RT, CR = NULL, stopping.rule = "OR",
                          Condition = NULL, Subject = NULL,
                          by_subject = "never", include_nn = FALSE) {
  if (!is.data.frame(RT)) return(list(RT = RT, CR = CR, converted = FALSE))
  if (!is.null(CR)) {
    stop("When RT is a data.frame, correctness must come from its Correct column; omit CR.",
         call. = FALSE)
  }
  converted <- sft_data_to_rt(RT, Condition = Condition, Subject = Subject,
                              stopping.rule = stopping.rule,
                              by_subject = by_subject, include_nn = include_nn)
  list(RT = converted$RT, CR = converted$CR, converted = TRUE)
}

.sft_ecdf <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(.sft_zero_curve(0))
  stats::ecdf(x)
}

# Helper_Functions.R compatibility helpers used by the simulation and plotting
# workflows. They remain internal because they are data-preparation utilities,
# not part of the original sft public API.
row_to_pvec <- function(df_row) {
  v <- as.numeric(df_row[1, , drop = TRUE])
  names(v) <- names(df_row)
  v
}

.safe_ecdf <- .sft_ecdf

.sft_hazard_env <- new.env(parent = baseenv())

.ensure_hazard_rcpp <- local({
  compiled <- FALSE
  function() {
    if (compiled && exists("nah_nak_eval_rcpp", envir = .sft_hazard_env,
                           mode = "function")) return(TRUE)
    if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
    ok <- tryCatch({
      Rcpp::sourceCpp(code = '
        #include <Rcpp.h>
        #include <vector>
        #include <algorithm>
        using namespace Rcpp;

        struct RTCR { double rt; bool cr; };

        static NumericMatrix eval_hazard(NumericVector rt, LogicalVector cr,
                                         NumericVector q, bool reverse) {
          NumericMatrix out(q.size(), 2);
          std::vector<RTCR> v;
          for (int i = 0; i < rt.size(); ++i) {
            if (!R_finite(rt[i])) continue;
            bool event = i < cr.size() && cr[i] != NA_LOGICAL && cr[i] == TRUE;
            v.push_back(RTCR{rt[i], event});
          }
          std::sort(v.begin(), v.end(), [](const RTCR& a, const RTCR& b) {
            return a.rt < b.rt;
          });
          if (v.empty()) return out;
          std::vector<double> event_rt, cum_h, cum_v;
          double h = 0.0, vv = 0.0;
          for (int i = 0; i < v.size(); ++i) {
            if (!v[i].cr) continue;
            double risk = reverse ? (i + 1.0) : (v.size() - i);
            double inc = 1.0 / risk;
            h += inc; vv += inc * inc;
            event_rt.push_back(v[i].rt);
            cum_h.push_back(h); cum_v.push_back(vv);
          }
          if (event_rt.empty()) return out;
          if (!reverse) {
            for (int j = 0; j < q.size(); ++j) {
              if (!R_finite(q[j])) continue;
              int k = std::upper_bound(event_rt.begin(), event_rt.end(), q[j]) - event_rt.begin() - 1;
              if (k >= 0) { out(j, 0) = cum_h[k]; out(j, 1) = cum_v[k]; }
            }
          } else {
            std::vector<double> tail_h(event_rt.size()), tail_v(event_rt.size());
            double th = 0.0, tv = 0.0;
            for (int i = event_rt.size() - 1; i >= 0; --i) {
              double inc_h = i == 0 ? cum_h[0] : cum_h[i] - cum_h[i - 1];
              double inc_v = i == 0 ? cum_v[0] : cum_v[i] - cum_v[i - 1];
              th += inc_h; tv += inc_v;
              tail_h[i] = -th; tail_v[i] = tv;
            }
            for (int j = 0; j < q.size(); ++j) {
              if (!R_finite(q[j])) continue;
              int k = std::lower_bound(event_rt.begin(), event_rt.end(), q[j]) - event_rt.begin();
              if (k < event_rt.size()) { out(j, 0) = tail_h[k]; out(j, 1) = tail_v[k]; }
            }
          }
          return out;
        }

        // [[Rcpp::export]]
        NumericMatrix nah_nak_eval_rcpp(NumericVector rt, LogicalVector cr,
                                        NumericVector q, bool reverse) {
          return eval_hazard(rt, cr, q, reverse);
        }

        // [[Rcpp::export]]
        NumericMatrix ucip_eval_rcpp(List rt_list, List cr_list,
                                     NumericVector q, bool reverse) {
          NumericMatrix out(q.size(), 2);
          for (int i = 0; i < rt_list.size(); ++i) {
            NumericVector rt = as<NumericVector>(rt_list[i]);
            LogicalVector cr = i < cr_list.size() && !Rf_isNull(cr_list[i])
              ? as<LogicalVector>(cr_list[i]) : LogicalVector(rt.size(), true);
            NumericMatrix one = eval_hazard(rt, cr, q, reverse);
            for (int j = 0; j < q.size(); ++j) {
              out(j, 0) += one(j, 0); out(j, 1) += one(j, 1);
            }
          }
          return out;
        }', env = .sft_hazard_env)
      TRUE
    }, error = function(e) FALSE)
    compiled <<- isTRUE(ok) &&
      exists("nah_nak_eval_rcpp", envir = .sft_hazard_env, mode = "function")
    compiled
  }
})

estimateNAH <- function(RT, CR = NULL) {
  x <- .sft_clean(RT, CR)
  if (!length(x$RT)) return(list(H = .sft_zero_curve(), Var = .sft_zero_curve()))
  if (.ensure_hazard_rcpp()) {
    rt <- x$RT; cr <- x$CR
    H <- function(t) .sft_hazard_env$nah_nak_eval_rcpp(rt, cr, as.numeric(t), FALSE)[, 1]
    V <- function(t) .sft_hazard_env$nah_nak_eval_rcpp(rt, cr, as.numeric(t), FALSE)[, 2]
    return(list(H = H, Var = V))
  }
  ord <- order(x$RT)
  rt <- x$RT[ord]; cr <- x$CR[ord]
  risk <- rev(seq_along(rt))
  hit <- which(cr)
  if (!length(hit)) return(list(H = .sft_zero_curve(), Var = .sft_zero_curve()))
  list(H = stats::stepfun(rt[hit], c(0, cumsum(1 / risk[hit]))),
       Var = stats::stepfun(rt[hit], c(0, cumsum(1 / risk[hit]^2))))
}

estimateNAK <- function(RT, CR = NULL) {
  x <- .sft_clean(RT, CR)
  if (!length(x$RT)) return(list(K = .sft_zero_curve(), Var = .sft_zero_curve()))
  if (.ensure_hazard_rcpp()) {
    rt <- x$RT; cr <- x$CR
    K <- function(t) .sft_hazard_env$nah_nak_eval_rcpp(rt, cr, as.numeric(t), TRUE)[, 1]
    V <- function(t) .sft_hazard_env$nah_nak_eval_rcpp(rt, cr, as.numeric(t), TRUE)[, 2]
    return(list(K = K, Var = V))
  }
  ord <- order(x$RT)
  rt <- x$RT[ord]; cr <- x$CR[ord]
  risk <- seq_along(rt)
  hit <- which(cr)
  if (!length(hit)) return(list(K = .sft_zero_curve(), Var = .sft_zero_curve()))
  list(K = stats::stepfun(rt[hit], c(rev(-cumsum(rev(1 / risk[hit]))), 0), right = TRUE),
       Var = stats::stepfun(rt[hit], c(rev(cumsum(rev(1 / risk[hit]^2))), 0), right = TRUE))
}

estimateUCIPor <- function(RT, CR = NULL, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (!is.list(RT) || !length(RT)) stop("RT must be a non-empty list.")
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  if (.ensure_hazard_rcpp()) {
    rt <- lapply(RT, as.numeric)
    cr <- lapply(seq_along(rt), function(i) .sft_correct(CR[[i]], length(rt[[i]])))
    H <- function(t) .sft_hazard_env$ucip_eval_rcpp(rt, cr, as.numeric(t), FALSE)[, 1]
    V <- function(t) .sft_hazard_env$ucip_eval_rcpp(rt, cr, as.numeric(t), FALSE)[, 2]
    return(list(H = H, Var = V))
  }
  hs <- lapply(seq_along(RT), function(i) estimateNAH(RT[[i]], CR[[i]]))
  h <- rowSums(vapply(hs, function(z) z$H(times), numeric(length(times))))
  v <- rowSums(vapply(hs, function(z) z$Var(times), numeric(length(times))))
  list(H = stats::stepfun(times, c(0, h)), Var = stats::stepfun(times, c(0, v)))
}

estimateUCIPand <- function(RT, CR = NULL, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "AND",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (!is.list(RT) || !length(RT)) stop("RT must be a non-empty list.")
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  if (.ensure_hazard_rcpp()) {
    rt <- lapply(RT, as.numeric)
    cr <- lapply(seq_along(rt), function(i) .sft_correct(CR[[i]], length(rt[[i]])))
    K <- function(t) .sft_hazard_env$ucip_eval_rcpp(rt, cr, as.numeric(t), TRUE)[, 1]
    V <- function(t) .sft_hazard_env$ucip_eval_rcpp(rt, cr, as.numeric(t), TRUE)[, 2]
    return(list(K = K, Var = V))
  }
  ks <- lapply(seq_along(RT), function(i) estimateNAK(RT[[i]], CR[[i]]))
  k <- rowSums(vapply(ks, function(z) z$K(times), numeric(length(times))))
  v <- rowSums(vapply(ks, function(z) z$Var(times), numeric(length(times))))
  list(K = stats::stepfun(times, c(k, 0), right = TRUE),
       Var = stats::stepfun(times, c(v, 0), right = TRUE))
}

.sft_capacity_bounds_or <- function(times, RT, CR, denom, ratio) {
  f1 <- .sft_ecdf(RT[[2]][.sft_correct(CR[[2]], length(RT[[2]]))])
  f2 <- .sft_ecdf(RT[[3]][.sft_correct(CR[[3]], length(RT[[3]]))])
  upper <- pmin(1, pmax(0, f1(times) + f2(times)))
  lower <- pmin(1, pmax(0, pmax(f1(times), f2(times))))
  hu <- -log1p(-upper); hl <- -log1p(-lower); dh <- denom(times)
  if (ratio) {
    cu <- hu / dh; cl <- hl / dh
  } else {
    cu <- hu - dh; cl <- hl - dh
  }
  cu[!is.finite(cu)] <- NA_real_; cl[!is.finite(cl)] <- NA_real_
  list(Ct_upper = .sft_curve(times, cu), Ct_lower = .sft_curve(times, cl),
       H_upper = .sft_curve(times, hu), H_lower = .sft_curve(times, hl),
       F_upper = .sft_curve(times, upper), F_lower = .sft_curve(times, lower))
}

.sft_capacity_bounds_and <- function(times, RT, CR, denom, ratio) {
  f1 <- .sft_ecdf(RT[[2]][.sft_correct(CR[[2]], length(RT[[2]]))])
  f2 <- .sft_ecdf(RT[[3]][.sft_correct(CR[[3]], length(RT[[3]]))])
  upper <- pmin(1, pmax(0, pmin(f1(times), f2(times))))
  lower <- pmin(1, pmax(0, f1(times) + f2(times) - 1))
  ku <- log(pmax(upper, .Machine$double.eps))
  kl <- log(pmax(lower, .Machine$double.eps))
  d <- denom(times)
  if (ratio) {
    a <- kl / d; b <- ku / d
    cu <- pmax(a, b); cl <- pmin(a, b)
  } else {
    cu <- d - kl; cl <- d - ku
  }
  cu[!is.finite(cu)] <- NA_real_; cl[!is.finite(cl)] <- NA_real_
  list(Ct_upper = .sft_curve(times, cu), Ct_lower = .sft_curve(times, cl),
       K_upper = .sft_curve(times, ku), K_lower = .sft_curve(times, kl),
       F_upper = .sft_curve(times, upper), F_lower = .sft_curve(times, lower))
}

capacity.or <- function(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (length(RT) < 3L) stop("capacity.or() needs redundant, A-alone, and B-alone RTs.")
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  numerator <- estimateNAH(RT[[1]], CR[[1]])
  denominator <- estimateUCIPor(RT[-1], CR[-1])
  ctest <- ucip.test(RT, CR, OR = TRUE)
  h <- numerator$H(times); d <- denominator$H(times)
  if (ratio) {
    ct <- h / d; ct[!is.finite(ct)] <- NA_real_
    out <- list(Ct = .sft_curve(times, ct), Ctest = ctest, times = times)
  } else {
    ct <- h - d; v <- numerator$Var(times) + denominator$Var(times)
    out <- list(Ct = .sft_curve(base::c(0, times), base::c(0, ct)), Var = .sft_curve(base::c(0, times), base::c(0, v)),
                Ctest = ctest, p.val = ctest$p.value, times = times)
  }
  b <- .sft_capacity_bounds_or(times, RT, CR, denominator$H, ratio)
  out[c("Ct_upper", "Ct_lower")] <- b[c("Ct_upper", "Ct_lower")]
  out
}

capacity.and <- function(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "AND",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (length(RT) < 3L) stop("capacity.and() needs redundant, A-alone, and B-alone RTs.")
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  denominator <- estimateNAK(RT[[1]], CR[[1]])
  numerator <- estimateUCIPand(RT[-1], CR[-1])
  ctest <- ucip.test(RT, CR, OR = FALSE)
  d <- denominator$K(times); k <- numerator$K(times)
  if (ratio) {
    ct <- k / d; ct[!is.finite(ct)] <- NA_real_
    out <- list(Ct = .sft_curve(times, ct), Ctest = ctest, times = times)
  } else {
    ct <- d - k; v <- denominator$Var(times) + numerator$Var(times)
    out <- list(Ct = .sft_curve(base::c(times, Inf), base::c(ct, 0)), Var = .sft_curve(base::c(times, Inf), base::c(v, 0)),
                Ctest = ctest, times = times)
  }
  b <- .sft_capacity_bounds_and(times, RT, CR, denominator$K, ratio)
  out[c("Ct_upper", "Ct_lower")] <- b[c("Ct_upper", "Ct_lower")]
  out
}

capacity.stst <- function(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "STST",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (length(RT) < 2L) stop("capacity.stst() needs context and target-alone RTs.")
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  denominator <- estimateNAK(RT[[1]], CR[[1]])
  numerator <- estimateNAK(RT[[2]], CR[[2]])
  ctest <- ucip.test(RT, CR, OR = FALSE)
  d <- denominator$K(times); k <- numerator$K(times)
  if (ratio) {
    ct <- k / d; ct[!is.finite(ct)] <- NA_real_
    list(Ct = .sft_curve(times, ct), Ctest = ctest, times = times)
  } else {
    ct <- d - k; v <- denominator$Var(times) + numerator$Var(times)
    list(Ct = .sft_curve(base::c(0, times), base::c(0, ct)), Var = .sft_curve(base::c(0, times), base::c(0, v)),
         Ctest = ctest, p.val = ctest$p.value, times = times)
  }
}

capacity.id <- function(dt.rt, nt.rt, st.rts, dt.cr = NULL, nt.cr = NULL,
                        st.crs = NULL, ratio = TRUE) {
  if (!is.list(st.rts) || !length(st.rts)) stop("st.rts must be a non-empty list.")
  RT <- c(list(dt.rt, nt.rt), st.rts)
  CR <- c(list(dt.cr, nt.cr), st.crs)
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  ks <- lapply(seq_along(RT), function(i) estimateNAK(RT[[i]], CR[[i]]))
  denom <- ks[[1]]$K(times) + ks[[2]]$K(times)
  numer <- rowSums(vapply(ks[-c(1, 2)], function(z) z$K(times), numeric(length(times))))
  ctest <- ucip.id.test(dt.rt, nt.rt, st.rts, CR[[1]], CR[[2]], CR[-c(1, 2)])
  if (ratio) {
    ct <- numer / denom; ct[!is.finite(ct)] <- NA_real_
    return(list(Ct = .sft_curve(times, ct), Ctest = ctest, times = times))
  }
  v <- ks[[1]]$Var(times) + ks[[2]]$Var(times) +
    rowSums(vapply(ks[-c(1, 2)], function(z) z$Var(times), numeric(length(times))))
  list(Ct = .sft_curve(c(times, Inf), c(denom - numer, 0)),
       Var = .sft_curve(c(times, Inf), c(v, 0)), Ctest = ctest, times = times)
}

capacity.altieri <- function(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject,
                             include_nn = TRUE)
  RT <- converted$RT; CR <- converted$CR
  if (length(RT) < 4L) stop("capacity.altieri() needs AB, A, B, and NN RTs.")
  CR <- .sft_cr_list(RT, CR)
  times <- .sft_finite_times(RT)
  num <- estimateNAH(RT[[1]], CR[[1]])
  d1 <- estimateUCIPor(RT[2:3], CR[2:3])
  d2 <- estimateNAH(RT[[4]], CR[[4]])
  h <- num$H(times)
  if (ratio) {
    a1 <- h / d1$H(times); a2 <- h / d2$H(times)
  } else {
    a1 <- h - d1$H(times); a2 <- h - d2$H(times)
  }
  a1[!is.finite(a1)] <- NA_real_; a2[!is.finite(a2)] <- NA_real_
  list(Ct_a1 = .sft_curve(times, a1), Ct_a2 = .sft_curve(times, a2), times = times)
}

.sft_htest <- function(statistic, p.value, alternative, method, data.name = NULL) {
  ans <- list(statistic = statistic, p.value = p.value,
              alternative = alternative, method = method)
  if (!is.null(data.name)) ans$data.name <- data.name
  class(ans) <- "htest"
  ans
}

.sft_ucip_components <- function(RT, CR, OR = NULL, stopping.rule = NULL) {
  if (!is.null(OR)) {
    rule <- if (isTRUE(OR)) "OR" else "AND"
  } else {
    rule <- match.arg(stopping.rule %||% "OR", c("OR", "AND", "STST"))
  }
  if (!is.list(RT) || length(RT) < 2L) stop("RT must be a list of at least two conditions.")
  CR <- .sft_cr_list(RT, CR)
  rt <- lapply(RT, as.numeric)
  all_rt <- unlist(rt, use.names = FALSE)
  keep <- is.finite(all_rt)
  if (!any(keep)) stop("RT contains no finite response times.")
  tvec <- sort(all_rt[keep])
  Y <- t(vapply(rt, function(x) vapply(tvec, function(t) sum(is.finite(x) & x >= t), numeric(1)),
               numeric(length(tvec))))
  G <- t(vapply(rt, function(x) vapply(tvec, function(t) sum(is.finite(x) & x <= t), numeric(1)),
               numeric(length(tvec))))
  index <- rep(seq_along(rt), lengths(rt))
  flat_rt <- unlist(rt, use.names = FALSE)
  flat_cr <- unlist(CR, use.names = FALSE)
  ord <- order(flat_rt)
  list(rule = rule, rt = rt, cr = CR, tvec = tvec, Y = Y, G = G,
       index = index[ord], cr_sorted = flat_cr[ord], rt_sorted = flat_rt[ord])
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.sft_score_method <- function(method) {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(method[[1L]])))
  switch(key,
         score = "score", theta = "score",
         capacity = "capacity", logcapacity = "capacity", eta = "capacity",
         stop("score_method must be one of 'score' or 'capacity'."))
}

.sft_ucip_score <- function(RT, CR = NULL, OR = NULL,
                            stopping.rule = c("OR", "AND", "STST")) {
  z <- .sft_ucip_components(RT, CR, OR, stopping.rule)
  ncond <- length(z$rt)
  if (z$rule == "OR") {
    risk <- z$Y
    weight <- risk[1, ] * colSums(risk[-1, , drop = FALSE]) / colSums(risk)
    signs <- c(1, rep(-1, ncond - 1L))
    alternative <- "response times are different than those predicted by the UCIP-OR model"
  } else {
    risk <- z$G
    weight <- risk[1, ] * colSums(risk[-1, , drop = FALSE]) / colSums(risk)
    signs <- c(-1, rep(1, ncond - 1L))
    alternative <- paste0("response times are different than those predicted by the UCIP-",
                          z$rule, " model")
  }
  weight[!is.finite(weight)] <- 0
  event <- z$cr_sorted & is.finite(z$rt_sorted)
  cond <- z$index

  capacity_numerator <- 0
  capacity_denominator <- 0
  numerator_variance <- 0
  denominator_variance <- 0

  for (i in seq_len(ncond)) {
    idx <- event & cond == i

    time_index <- findInterval(
      z$rt_sorted[idx], z$tvec, rightmost.closed = TRUE
    )

    denom_i <- risk[i, time_index]
    denom_i[denom_i <= 0] <- NA_real_

    terms <- weight[time_index] / denom_i
    terms <- terms[is.finite(terms)]

    component <- sum(terms)
    component_variance <- sum(terms^2)

    if (signs[i] > 0) {
      capacity_numerator <- capacity_numerator + component
      numerator_variance <- numerator_variance + component_variance
    } else {
      capacity_denominator <- capacity_denominator + component
      denominator_variance <- denominator_variance + component_variance
    }
  }

  numer <- capacity_numerator - capacity_denominator
  variance <- numerator_variance + denominator_variance
  denom <- sqrt(variance)

  statistic <- if (denom > 0) numer / denom else NA_real_
  names(statistic) <- "z"
  p <- if (is.finite(statistic)) 2 * min(stats::pnorm(statistic), 1 - stats::pnorm(statistic)) else NA_real_

  capacity_ok <- is.finite(capacity_numerator) &&
    is.finite(capacity_denominator) &&
    capacity_numerator > 0 && capacity_denominator > 0

  log_capacity <- if (capacity_ok) {
    log(capacity_numerator / capacity_denominator)
  } else {
    NA_real_
  }

  log_capacity_variance <- if (capacity_ok) {
    numerator_variance / capacity_numerator^2 +
      denominator_variance / capacity_denominator^2
  } else {
    NA_real_
  }

  log_capacity_precision <- if (is.finite(log_capacity_variance) &&
                                log_capacity_variance > 0) {
    1 / log_capacity_variance
  } else {
    NA_real_
  }

  list(numer = numer, variance = variance, denom = denom,
       statistic = statistic, p.value = p, alternative = alternative,
       rule = z$rule,
       capacity_numerator = capacity_numerator,
       capacity_denominator = capacity_denominator,
       numerator_variance = numerator_variance,
       denominator_variance = denominator_variance,
       capacity_ratio = if (capacity_ok) {
         capacity_numerator / capacity_denominator
       } else {
         NA_real_
       },
       log_capacity = log_capacity,
       log_capacity_variance = log_capacity_variance,
       log_capacity_precision = log_capacity_precision)
}

ucip.test <- function(RT, CR = NULL, OR = NULL, stopping.rule = c("OR", "AND", "STST"),
                      Condition = NULL, Subject = NULL) {
  rule <- if (!is.null(OR)) if (isTRUE(OR)) "OR" else "AND" else match.arg(stopping.rule)
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = rule,
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  score <- .sft_ucip_score(RT, CR, OR, stopping.rule)
  out <- .sft_htest(score$statistic, score$p.value, score$alternative,
                    "Houpt-Townsend UCIP test", deparse(substitute(RT)))
  # Keep the score components available to the hierarchical Bayesian companion
  # without changing the established htest statistic or p-value interface.
  out$numer <- score$numer
  out$variance <- score$variance
  out$denom <- score$denom
  out
}

ucip.id.test <- function(dt.rt, nt.rt, st.rts, dt.cr = NULL, nt.cr = NULL, st.crs = NULL) {
  if (!is.list(st.rts) || !length(st.rts)) stop("st.rts must be a non-empty list.")
  RT <- c(list(dt.rt, nt.rt), st.rts)
  CR <- .sft_cr_list(RT, c(list(dt.cr, nt.cr), st.crs))
  z <- .sft_ucip_components(RT, CR, OR = FALSE)
  risk <- z$G
  weight <- (risk[1, ] + risk[2, ]) * colSums(risk[-c(1, 2), , drop = FALSE]) / colSums(risk)
  event <- z$cr_sorted & is.finite(z$rt_sorted)
  cond <- z$index
  numer <- 0; denom_terms <- numeric(0)
  for (i in seq_len(nrow(risk))) {
    idx <- event & cond == i
    ti <- findInterval(z$rt_sorted[idx], z$tvec, rightmost.closed = TRUE)
    ri <- risk[i, ti]; ri[ri <= 0] <- NA_real_
    terms <- weight[ti] / ri
    terms <- terms[is.finite(terms)]
    numer <- numer + if (i <= 2L) -sum(terms) else sum(terms)
    denom_terms <- c(denom_terms, terms^2)
  }
  denom <- sqrt(sum(denom_terms))
  statistic <- if (denom > 0) numer / denom else NA_real_
  names(statistic) <- "z"
  p <- if (is.finite(statistic)) 2 * min(stats::pnorm(statistic), 1 - stats::pnorm(statistic)) else NA_real_
  .sft_htest(statistic, p,
             "response times are different than those predicted by the adjusted UCIP-AND model",
             "Houpt-Townsend UCIP test")
}

.sft_hdi <- function(x, mass = .94) {
  x <- sort(as.numeric(x[is.finite(x)]))
  if (!length(x)) return(c(lower = NA_real_, upper = NA_real_))
  n <- length(x)
  width <- max(0L, ceiling(mass * n) - 1L)
  if (width == 0L) return(c(lower = x[1L], upper = x[1L]))
  left <- seq_len(n - width)
  j <- left[which.min(x[left + width] - x[left])]
  c(lower = x[j], upper = x[j + width])
}

.sft_subject_ucip_input <- function(RT, CR = NULL) {
  if (!is.list(RT) || !length(RT)) stop("RT must be a non-empty list.")
  nested <- all(vapply(RT, is.list, logical(1)))
  if (!nested) {
    return(list(RT = list(RT), CR = list(CR), subject = "1"))
  }
  if (is.null(CR)) CR <- vector("list", length(RT))
  if (!is.list(CR) || length(CR) != length(RT) ||
      any(vapply(CR, function(x) !is.null(x) && !is.list(x), logical(1)))) {
    stop("For nested RT input, CR must be a matching list of participant CR lists.")
  }
  subject <- names(RT)
  if (is.null(subject) || any(!nzchar(subject))) subject <- as.character(seq_along(RT))
  list(RT = RT, CR = CR, subject = subject)
}

.sft_hierarchical_cz_chain <- function(theta_hat, precision, ndraws, burnin,
                                       thin, prior_mean, prior_sd,
                                       tau2_prior_shape, tau2_prior_rate) {
  n <- length(theta_hat)
  niter <- burnin + ndraws * thin
  mu_draws <- numeric(ndraws)
  tau2_draws <- numeric(ndraws)
  theta_draws <- matrix(NA_real_, nrow = ndraws, ncol = n)
  mu <- stats::weighted.mean(theta_hat, precision)
  theta <- theta_hat
  tau2 <- 1 / stats::rgamma(1L, shape = tau2_prior_shape, rate = tau2_prior_rate)
  save_i <- 0L

  for (iter in seq_len(niter)) {
    theta_var <- 1 / (precision + 1 / tau2)
    theta_mean <- theta_var * (precision * theta_hat + mu / tau2)
    theta <- stats::rnorm(n, theta_mean, sqrt(theta_var))

    mu_var <- 1 / (n / tau2 + 1 / prior_sd^2)
    mu_mean <- mu_var * (sum(theta) / tau2 + prior_mean / prior_sd^2)
    mu <- stats::rnorm(1L, mu_mean, sqrt(mu_var))

    # tau^2 has an inverse-Gamma full conditional. Sampling its reciprocal
    # as Gamma keeps the update stable and avoids a parameter at zero.
    tau2 <- 1 / stats::rgamma(
      1L, shape = tau2_prior_shape + n / 2,
      rate = tau2_prior_rate + sum((theta - mu)^2) / 2
    )

    if (iter > burnin && (iter - burnin) %% thin == 0L) {
      save_i <- save_i + 1L
      mu_draws[save_i] <- mu
      tau2_draws[save_i] <- tau2
      theta_draws[save_i, ] <- theta
    
      if (save_i == ndraws) break
    }
  }
  list(mu = mu_draws, tau = sqrt(tau2_draws), theta = theta_draws)
}

.sft_bayes_method <- function(method) {
  if (length(method) == 0L || is.null(method)) return("InvGamma")
  key <- tolower(gsub("[^A-Za-z]", "", as.character(method[[1L]])))
  switch(key,
         invgamma = "InvGamma", invgammagibbs = "InvGamma",
         ig = "InvGamma", gibbs = "InvGamma",
         halfnormal = "HalfNormal", hnorm = "HalfNormal",
         halfnormalstan = "HalfNormal", stanhalfnormal = "HalfNormal",
         centered = "Centered", centred = "Centered",
         centeredstan = "Centered", centredstan = "Centered",
         stancentered = "Centered", stancentred = "Centered",
         stop("method must be one of 'InvGamma', 'HalfNormal', or 'Centered'."))
}

.sft_validate_bayes_args <- function(ndraws, burnin, thin, chains,
                                     prior_shape, prior_rate, prior_mean,
                                     prior_sd, prior_tau_sd, hdi, rope,
                                     adapt_delta, max_treedepth) {
  if (ndraws < 100L || ndraws != as.integer(ndraws)) {
    stop("ndraws must be an integer >= 100.")
  }
  if (burnin < 0L || burnin != as.integer(burnin)) {
    stop("burnin must be a non-negative integer.")
  }
  if (thin < 1L || thin != as.integer(thin)) {
    stop("thin must be a positive integer.")
  }
  if (chains < 2L || chains != as.integer(chains)) {
    stop("chains must be an integer >= 2.")
  }
  if (!is.finite(prior_shape) || prior_shape <= 0 ||
      !is.finite(prior_rate) || prior_rate <= 0) {
    stop("prior_shape and prior_rate must be positive finite values.")
  }
  if (!is.finite(prior_mean) || !is.finite(prior_sd) || prior_sd <= 0) {
    stop("prior_mean must be finite and prior_sd must be positive.")
  }
  if (!is.finite(prior_tau_sd) || prior_tau_sd <= 0) {
    stop("prior_tau_sd must be positive and finite.")
  }
  if (!is.finite(hdi) || hdi <= 0 || hdi >= 1) stop("hdi must lie in (0, 1).")
  if (!is.null(rope) && (!is.finite(rope) || rope <= 0)) {
    stop("rope must be NULL or a positive finite practical-equivalence threshold.")
  }
  if (!is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
    stop("adapt_delta must lie in (0, 1).")
  }
  if (max_treedepth < 1L || max_treedepth != as.integer(max_treedepth)) {
    stop("max_treedepth must be a positive integer.")
  }
  invisible(TRUE)
}

.sft_ucip_score_data <- function(input, stopping.rule,
                                 method = c("score", "capacity")) {
  method <- .sft_score_method(method)
  scores <- lapply(seq_along(input$RT), function(i) {
    score <- .sft_ucip_score(input$RT[[i]], input$CR[[i]],
                              stopping.rule = stopping.rule)
    if (!is.finite(score$numer) || !is.finite(score$variance) || score$variance <= 0) {
      stop("Each participant must have a finite UCIP score with positive variance.")
    }
    if (method == "capacity" &&
        (!is.finite(score$log_capacity) ||
         !is.finite(score$log_capacity_precision) ||
         score$log_capacity_precision <= 0)) {
      stop("Each participant must have an estimable positive weighted capacity numerator and denominator.")
    }
    score
  })

  theta_hat <- vapply(scores, function(x) x$numer / x$variance, numeric(1))
  theta_precision <- vapply(scores, function(x) x$variance, numeric(1))
  eta_hat <- vapply(scores, function(x) x$log_capacity, numeric(1))
  eta_precision <- vapply(scores, function(x) x$log_capacity_precision, numeric(1))
  U <- vapply(scores, function(x) x$numer, numeric(1))
  V <- vapply(scores, function(x) x$variance, numeric(1))

  effect_hat <- if (method == "capacity") eta_hat else theta_hat
  precision <- if (method == "capacity") eta_precision else theta_precision
  effect_name <- if (method == "capacity") "eta" else "theta"

  score_df <- if (method == "score") {
    data.frame(
      subject = input$subject,
      U = U,
      V = V,
      theta_hat = theta_hat,
      se = 1 / sqrt(theta_precision),
      Cz = vapply(scores, function(x) unname(x$statistic), numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      subject = input$subject,
      U = U,
      V = V,
      Cz = vapply(scores, function(x) unname(x$statistic), numeric(1)),
      capacity_numerator = vapply(
        scores, function(x) x$capacity_numerator, numeric(1)
      ),
      capacity_denominator = vapply(
        scores, function(x) x$capacity_denominator, numeric(1)
      ),
      capacity_ratio_hat = exp(eta_hat),
      eta_hat = eta_hat,
      eta_se = 1 / sqrt(eta_precision),
      score_effect_hat = U / V,
      score_effect_se = 1 / sqrt(V),
      stringsAsFactors = FALSE
    )
  }

  list(scores = scores, theta_hat = theta_hat, eta_hat = eta_hat,
       precision = precision, theta_precision = theta_precision,
       eta_precision = eta_precision, effect_hat = effect_hat,
       effect_name = effect_name, method = method, score = score_df)
}

.sft_bayes_seed <- function(seed) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  if (!is.null(seed)) set.seed(seed)
  old_seed
}

.sft_restore_bayes_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
  invisible(NULL)
}

.sft_parameter_summary <- function(parameter_draws, parameter_names, hdi, rope) {
  summaries <- lapply(seq_len(ncol(parameter_draws)), function(j) {
    x <- parameter_draws[, j]
    interval <- .sft_hdi(x, hdi)
    rope_probability <- if (is.null(rope)) NA_real_ else mean(abs(x) <= rope)
    data.frame(parameter = parameter_names[j], mean = mean(x),
               lower = interval[["lower"]], upper = interval[["upper"]],
               posterior_positive = mean(x > 0),
               posterior_negative = mean(x < 0),
               posterior_rope = rope_probability,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, summaries)
}

.sft_split_chains <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n < 4L) return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  half <- floor(n / 2L)
  out <- matrix(NA_real_, nrow = half, ncol = ncol(x) * 2L)
  for (j in seq_len(ncol(x))) {
    out[, 2L * j - 1L] <- x[seq_len(half), j]
    out[, 2L * j] <- x[n - half + seq_len(half), j]
  }
  out
}

.sft_fallback_ess <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (n < 3L || stats::sd(x) == 0) return(as.numeric(n))
  lag_max <- min(n - 1L, 1000L)
  ac <- stats::acf(x, plot = FALSE, lag.max = lag_max, demean = TRUE)$acf[, 1L, 1L]
  rho <- ac[-1L]
  pairs <- floor(length(rho) / 2L)
  if (!pairs) return(as.numeric(n))
  gamma <- rho[2L * seq_len(pairs) - 1L] + rho[2L * seq_len(pairs)]
  keep <- which(gamma > 0)
  if (!length(keep)) return(as.numeric(n))
  first_nonpositive <- which(gamma <= 0)
  last <- if (length(first_nonpositive)) min(first_nonpositive) - 1L else length(gamma)
  if (last < 1L) return(as.numeric(n))
  tau_int <- max(1, 1 + 2 * sum(gamma[seq_len(last)]))
  min(as.numeric(n), as.numeric(n) / tau_int)
}

.sft_manual_diagnostics <- function(chain_array, parameter_names) {
  out <- lapply(seq_along(parameter_names), function(j) {
    chains <- chain_array[, , j, drop = TRUE]
    split <- .sft_split_chains(chains)
    if (!length(split)) {
      return(data.frame(parameter = parameter_names[j], rhat = NA_real_,
                        ess_bulk = NA_real_, ess_tail = NA_real_,
                        stringsAsFactors = FALSE))
    }
    within <- mean(apply(split, 2L, stats::var))
    between <- nrow(split) * stats::var(colMeans(split))
    rhat <- if (is.finite(within) && within > 0) {
      sqrt(((nrow(split) - 1) * within + between) /
             (nrow(split) * within))
    } else if (is.finite(between) && between == 0) 1 else NA_real_
    flat <- as.vector(split)
    bulk <- stats::qnorm((rank(flat, ties.method = "average") - .5) / length(flat))
    q <- stats::quantile(flat, c(.05, .95), names = FALSE, type = 8)
    tail <- as.numeric(flat <= q[1L] | flat >= q[2L])
    data.frame(parameter = parameter_names[j], rhat = rhat,
               ess_bulk = .sft_fallback_ess(bulk),
               ess_tail = .sft_fallback_ess(tail),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.sft_mcmc_diagnostics <- function(chain_array, parameter_names) {
  if (requireNamespace("posterior", quietly = TRUE)) {
    dimnames(chain_array) <- list(NULL, paste0("chain", seq_len(dim(chain_array)[2L])),
                                  parameter_names)
    ans <- tryCatch({
      draws <- posterior::as_draws_array(chain_array)
      stats <- posterior::summarise_draws(
        draws, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
      )
      data.frame(parameter = as.character(stats$variable),
                 rhat = as.numeric(stats[[2L]]),
                 ess_bulk = as.numeric(stats[[3L]]),
                 ess_tail = as.numeric(stats[[4L]]),
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(ans)) return(ans)
  }
  .sft_manual_diagnostics(chain_array, parameter_names)
}

.sft_prior_predictive <- function(effect_hat, precision, ndraws, method,
                                  prior_shape, prior_rate, prior_mean,
                                  prior_sd, prior_tau_sd,
                                  effect_name = "theta") {
  prior_mu <- stats::rnorm(ndraws, prior_mean, prior_sd)
  if (method == "InvGamma") {
    prior_tau2 <- 1 / stats::rgamma(ndraws, shape = prior_shape, rate = prior_rate)
    prior_tau <- sqrt(prior_tau2)
  } else {
    prior_tau <- abs(stats::rnorm(ndraws, 0, prior_tau_sd))
    prior_tau2 <- prior_tau^2
  }
  n <- length(effect_hat)
  prior_effect <- matrix(
    stats::rnorm(ndraws * n,
                 mean = rep(prior_mu, each = n),
                 sd = rep(prior_tau, each = n)),
    nrow = ndraws, ncol = n, byrow = TRUE
  )
  prior_score <- prior_effect + sweep(
    matrix(stats::rnorm(length(prior_effect)), nrow = ndraws),
    2L, 1 / sqrt(precision), "*"
  )
  out <- list(mu = prior_mu, tau = prior_tau, score = prior_score)
  out[[effect_name]] <- prior_effect
  out[c("mu", "tau", effect_name, "score")]
}

.sft_posterior_predictive <- function(effect_draws, effect_hat, precision) {
  posterior_score <- effect_draws + sweep(
    matrix(stats::rnorm(length(effect_draws)), nrow = nrow(effect_draws)),
    2L, 1 / sqrt(precision), "*"
  )
  observed_mean <- mean(effect_hat)
  observed_sd <- if (length(effect_hat) > 1L) stats::sd(effect_hat) else 0
  predictive_mean <- rowMeans(posterior_score)
  predictive_sd <- if (ncol(posterior_score) > 1L) {
    apply(posterior_score, 1L, stats::sd)
  } else {
    rep(0, nrow(posterior_score))
  }
  list(
    draws = posterior_score,
    observed = effect_hat,
    discrepancy = data.frame(
      measure = c("mean", "sd"),
      observed = c(observed_mean, observed_sd),
      predictive_mean = c(mean(predictive_mean), mean(predictive_sd)),
      posterior_tail_probability = c(
        mean(predictive_mean >= observed_mean),
        mean(predictive_sd >= observed_sd)
      ),
      stringsAsFactors = FALSE
    )
  )
}

.sft_shrinkage_draws <- function(tau_draws, precision) {
  sweep(
    matrix(tau_draws^2, nrow = length(tau_draws), ncol = length(precision)),
    2L, precision,
    function(tau2, information) information * tau2 / (1 + information * tau2)
  )
}

.sft_shrinkage_summary <- function(shrinkage_draws, subjects, hdi) {
  probs <- c((1 - hdi) / 2, 1 - (1 - hdi) / 2)
  data.frame(
    subject = subjects,
    mean_weight_on_individual = colMeans(shrinkage_draws),
    lower = apply(shrinkage_draws, 2L, stats::quantile, probs = probs[1L], names = FALSE),
    upper = apply(shrinkage_draws, 2L, stats::quantile, probs = probs[2L], names = FALSE),
    mean_weight_on_population = 1 - colMeans(shrinkage_draws),
    stringsAsFactors = FALSE
  )
}

.sft_stan_model_cache <- new.env(parent = emptyenv())

.sft_stan_code <- function(method, effect_name = "theta",
                           precision_name = "V") {
  effect_hat_name <- paste0(effect_name, "_hat")
  if (method == "HalfNormal") {
    return(paste0("data {
  int<lower=1> N;
  vector[N] ", effect_hat_name, ";
  vector<lower=0>[N] ", precision_name, ";
  real prior_mean;
  real<lower=0> prior_sd;
  real<lower=0> prior_tau_sd;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[N] z;
}
transformed parameters {
  vector[N] ", effect_name, ";
  ", effect_name, " = mu + tau * z;
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  z ~ normal(0, 1);
  for (i in 1:N)
    ", effect_hat_name, "[i] ~ normal(", effect_name, "[i], 1 / sqrt(", precision_name, "[i]));", "}"))
  }
  paste0("data {
  int<lower=1> N;
  vector[N] ", effect_hat_name, ";
  vector<lower=0>[N] ", precision_name, ";
  real prior_mean;
  real<lower=0> prior_sd;
  real<lower=0> prior_tau_sd;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[N] ", effect_name, ";
}
model {
  mu ~ normal(prior_mean, prior_sd);
  tau ~ normal(0, prior_tau_sd);
  ", effect_name, " ~ normal(mu, tau);
  for (i in 1:N)
    ", effect_hat_name, "[i] ~ normal(", effect_name, "[i], 1 / sqrt(", precision_name, "[i]));", "}")
}

.sft_stan_sampler_diagnostics <- function(fit, max_treedepth) {
  params <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                     error = function(e) NULL)
  if (is.null(params)) return(NULL)
  per_chain <- lapply(seq_along(params), function(i) {
    x <- params[[i]]
    divergent <- if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"] > 0) else NA_real_
    treedepth <- if ("treedepth__" %in% colnames(x)) sum(x[, "treedepth__"] >= max_treedepth) else NA_real_
    data.frame(chain = i, divergences = divergent,
               max_treedepth_hits = treedepth, stringsAsFactors = FALSE)
  })
  per_chain <- do.call(rbind, per_chain)
  list(per_chain = per_chain,
       divergences = sum(per_chain$divergences, na.rm = TRUE),
       max_treedepth_hits = sum(per_chain$max_treedepth_hits, na.rm = TRUE))
}

.sft_run_stan_hierarchy <- function(method, effect_hat, precision, ndraws,
                                    burnin, thin, chains, prior_mean,
                                    prior_sd, prior_tau_sd, seed,
                                    adapt_delta, max_treedepth, stan_control,
                                    effect_name = "theta") {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("method='", method, "' requires the optional 'rstan' package.")
  }
  precision_name <- if (effect_name == "eta") "P" else "V"
  code <- .sft_stan_code(method, effect_name, precision_name)
  key <- paste0(method, "_", as.character(prior_tau_sd), "_", effect_name)
  model <- if (exists(key, envir = .sft_stan_model_cache, inherits = FALSE)) {
    get(key, envir = .sft_stan_model_cache, inherits = FALSE)
  } else {
    fitted_model <- rstan::stan_model(model_code = code,
                                       model_name = paste0("ucip_bayes_", method),
                                       verbose = FALSE)
    assign(key, fitted_model, envir = .sft_stan_model_cache)
    fitted_model
  }
  control <- utils::modifyList(
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    stan_control
  )
  stan_seed <- if (is.null(seed)) -1L else seed
  stan_data <- list(N = length(effect_hat),
                    prior_mean = prior_mean, prior_sd = prior_sd,
                    prior_tau_sd = prior_tau_sd)
  stan_data[[paste0(effect_name, "_hat")]] <- effect_hat
  stan_data[[precision_name]] <- precision
  fit <- rstan::sampling(
    model,
    data = stan_data,
    chains = chains, iter = burnin + ndraws * thin, warmup = burnin,
    thin = thin, seed = stan_seed, refresh = 0, cores = 1L,
    control = control
  )
  raw <- rstan::extract(fit, pars = c("mu", "tau", effect_name),
                        permuted = FALSE, inc_warmup = FALSE)
  variables <- dimnames(raw)[[3L]]
  effect_idx <- grep(paste0("^", effect_name, "\\["), variables)
  if (length(effect_idx) != length(effect_hat)) {
    stop("Stan did not return the expected participant effects.")
  }
  effect_idx <- effect_idx[order(as.integer(gsub("[^0-9]", "", variables[effect_idx])))]
  chain_array <- array(NA_real_,
                       dim = c(dim(raw)[1L], dim(raw)[2L], 2L + length(effect_hat)))
  chain_array[, , 1L] <- raw[, , match("mu", variables)]
  chain_array[, , 2L] <- raw[, , match("tau", variables)]
  for (j in seq_along(effect_idx)) chain_array[, , 2L + j] <- raw[, , effect_idx[j]]
  dimnames(chain_array) <- list(NULL, NULL,
                                c("mu", "tau", paste0(effect_name, "[", seq_along(effect_hat), "]")))
  list(chain_array = chain_array, fit = fit,
       sampler_diagnostics = .sft_stan_sampler_diagnostics(fit, max_treedepth))
}

# Analytic normal-normal posterior for one participant. There is no
# population mean or between-participant scale in this model.
ucip.bayes.single_subject.test <- function(RT, CR = NULL,
                                           stopping.rule = c("OR", "AND", "STST"),
                                           ndraws = 10000L, seed = NULL,
                                           hdi = .94, prior_mean = 0, prior_sd = 1,
                                           chains = 4L, rope = NULL,
                                           Condition = NULL, Subject = NULL,
                                           score_method = c("score", "capacity")) {
  stopping.rule <- match.arg(stopping.rule)
  score_method <- .sft_score_method(score_method)
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = stopping.rule,
                             Condition = Condition, Subject = Subject,
                             by_subject = "never")
  RT <- converted$RT; CR <- converted$CR
  .sft_validate_bayes_args(ndraws, 0L, 1L, chains, 1, 1, prior_mean,
                           prior_sd, 1, hdi, rope, .95, 12L)
  input <- .sft_subject_ucip_input(RT, CR)
  if (length(input$RT) != 1L) {
    stop("ucip.bayes.single_subject.test() requires exactly one subject.")
  }
  scored <- .sft_ucip_score_data(input, stopping.rule, score_method)
  effect_hat <- scored$effect_hat[[1L]]
  precision <- scored$precision[[1L]]
  effect_name <- scored$effect_name
  observed_name <- paste0(effect_name, "_hat")
  precision_name <- if (score_method == "capacity") "P" else "V"
  post_var <- 1 / (precision + 1 / prior_sd^2)
  post_mean <- post_var * (precision * effect_hat + prior_mean / prior_sd^2)

  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  effect_chain <- matrix(stats::rnorm(ndraws * chains, mean = post_mean,
                                      sd = sqrt(post_var)),
                         nrow = ndraws, ncol = chains)
  effect_draws <- matrix(as.vector(effect_chain), ncol = 1L)
  subject_parameter <- paste0(effect_name, "[", input$subject[[1L]], "]")
  parameter_names <- subject_parameter
  summary <- .sft_parameter_summary(effect_draws, parameter_names, hdi, rope)
  rownames(summary) <- NULL
  subject_summary <- summary
  rownames(subject_summary) <- input$subject[[1L]]
  prior_effect <- stats::rnorm(ndraws, prior_mean, prior_sd)
  prior_score <- prior_effect + stats::rnorm(ndraws, 0, 1 / sqrt(precision))
  posterior_predictive <- .sft_posterior_predictive(effect_draws,
                                                     scored$effect_hat,
                                                     scored$precision)
  prior_predictive <- list(score = prior_score)
  prior_predictive[[effect_name]] <- prior_effect
  prior_predictive <- prior_predictive[c(effect_name, "score")]
  observation_weight <- precision / (precision + 1 / prior_sd^2)
  shrinkage_summary <- data.frame(
    subject = input$subject[[1L]],
    mean_weight_on_individual = observation_weight,
    lower = observation_weight,
    upper = observation_weight,
    mean_weight_on_population = NA_real_,
    stringsAsFactors = FALSE
  )
  draws <- data.frame(effect_draws, check.names = FALSE)
  names(draws) <- subject_parameter
  prior <- list(rope = rope)
  prior[[effect_name]] <- list(mean = prior_mean, sd = prior_sd)
  list(
    statistic = setNames(mean(effect_draws[, 1L]), subject_parameter),
    posterior_probability = list(
      subject_super = setNames(mean(effect_draws[, 1L] > 0), input$subject[[1L]]),
      subject_limited = setNames(mean(effect_draws[, 1L] < 0), input$subject[[1L]]),
      subject_rope = setNames(if (is.null(rope)) NA_real_ else
                                mean(abs(effect_draws[, 1L]) <= rope),
                              input$subject[[1L]])
    ),
    summary = summary, population_summary = NULL,
    subject_summary = subject_summary, score = scored$score, draws = draws,
    prior_predictive = prior_predictive,
    posterior_predictive = posterior_predictive,
    diagnostics = data.frame(parameter = subject_parameter, rhat = 1,
                              ess_bulk = ndraws * chains,
                              ess_tail = ndraws * chains,
                              stringsAsFactors = FALSE),
    shrinkage_summary = shrinkage_summary,
    shrinkage_draws = matrix(observation_weight, nrow = ndraws * chains, ncol = 1L),
    method = "Analytic single-subject normal-normal posterior",
    method_code = "Analytic",
    alternative = paste("the participant UCIP", stopping.rule, "effect is nonzero"),
    prior = prior,
    model = list(
      observation = sprintf("%s ~ Normal(mean = %s, sd = 1 / sqrt(%s))",
                            observed_name, effect_name, precision_name),
      subject = sprintf("%s ~ Normal(mean = prior_mean, sd = prior_sd)",
                        effect_name),
      interpretation = if (score_method == "capacity") {
        "exp(eta) is the participant weighted capacity ratio"
      } else {
        "theta is the unstandardised score-effect parameter"
      },
      score_method = score_method,
      hierarchy = FALSE
    ),
    score_method = score_method,
    stopping.rule = stopping.rule, hdi = hdi, ndraws = ndraws,
    burnin = 0L, thin = 1L, chains = chains, seed = seed
  )
}

# Hierarchical Bayesian companion to the UCIP/Cz score test. The default
# score_method = "score" preserves the original U/V score-effect model;
# score_method = "capacity" uses eta = log(A/B) with its delta-method
# precision. The InvGamma method is the exact conjugate Gibbs sampler.
# HalfNormal is a non-centred Stan model and Centered is its centred Stan
# counterpart; both use a half-normal prior on tau.
ucip.bayes.test <- function(RT, CR = NULL, stopping.rule = c("OR", "AND", "STST"),
                            ndraws = 10000L, prior_shape = 2, prior_rate = 1,
                            seed = NULL, hdi = .94, prior_mean = 0,
                            prior_sd = 1, burnin = 1000L, thin = 1L,
                            chains = 4L, rope = NULL,
                            method = c("InvGamma", "HalfNormal", "Centered"),
                            prior_tau_sd = 1, adapt_delta = .95,
                            max_treedepth = 12L, stan_control = list(),
                            tau2_prior_shape = NULL, tau2_prior_rate = NULL,
                            Condition = NULL, Subject = NULL,
                            score_method = c("score", "capacity")) {
  stopping.rule <- match.arg(stopping.rule)
  if (length(method) == 1L &&
      tolower(as.character(method)) %in% c("capacity", "log_capacity", "eta")) {
    score_method <- "capacity"
    method <- "InvGamma"
  }
  score_method <- .sft_score_method(score_method)
  method <- .sft_bayes_method(method)
  if (!is.list(stan_control)) stop("stan_control must be a list.")
  if (!is.null(tau2_prior_shape)) prior_shape <- tau2_prior_shape
  if (!is.null(tau2_prior_rate)) prior_rate <- tau2_prior_rate
  .sft_validate_bayes_args(ndraws, burnin, thin, chains, prior_shape,
                           prior_rate, prior_mean, prior_sd, prior_tau_sd,
                           hdi, rope, adapt_delta, max_treedepth)
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = stopping.rule,
                             Condition = Condition, Subject = Subject,
                             by_subject = "always")
  if (converted$converted) {
    RT <- converted$RT; CR <- converted$CR
  }
  input <- .sft_subject_ucip_input(RT, CR)

  # A one-subject input is deliberately routed to the analytic posterior.
  # The hierarchy cannot identify a population mean and between-subject scale
  # from one participant, even though priors would make the Gibbs model proper.
  if (length(input$RT) == 1L) {
    return(ucip.bayes.single_subject.test(
      RT = RT, CR = CR, stopping.rule = stopping.rule, ndraws = ndraws,
      seed = seed, hdi = hdi, prior_mean = prior_mean, prior_sd = prior_sd,
      chains = chains, rope = rope, score_method = score_method
    ))
  }

  scored <- .sft_ucip_score_data(input, stopping.rule, score_method)
  effect_hat <- scored$effect_hat
  precision <- scored$precision
  effect_name <- scored$effect_name
  observed_name <- paste0(effect_name, "_hat")
  precision_name <- if (score_method == "capacity") "P" else "V"
  old_seed <- .sft_bayes_seed(seed)
  on.exit(.sft_restore_bayes_seed(old_seed), add = TRUE)
  prior_predictive <- .sft_prior_predictive(
    effect_hat, precision, ndraws, method, prior_shape, prior_rate,
    prior_mean, prior_sd, prior_tau_sd, effect_name = effect_name
  )

  stan_fit <- NULL
  sampler_diagnostics <- NULL
  if (method == "InvGamma") {
    chain_draws <- lapply(seq_len(chains), function(chain) {
      .sft_hierarchical_cz_chain(effect_hat, precision, ndraws, burnin, thin,
                                 prior_mean, prior_sd, prior_shape, prior_rate)
    })
    chain_array <- array(NA_real_, dim = c(ndraws, chains, 2L + length(effect_hat)))
    for (j in seq_len(chains)) {
      chain_array[, j, 1L] <- chain_draws[[j]]$mu
      chain_array[, j, 2L] <- chain_draws[[j]]$tau
      chain_array[, j, seq_len(length(effect_hat)) + 2L] <- chain_draws[[j]]$theta
    }
  } else {
    stan <- .sft_run_stan_hierarchy(
      method, effect_hat, precision, ndraws, burnin, thin, chains,
      prior_mean, prior_sd, prior_tau_sd, seed, adapt_delta,
      max_treedepth, stan_control, effect_name = effect_name
    )
    chain_array <- stan$chain_array
    stan_fit <- stan$fit
    sampler_diagnostics <- stan$sampler_diagnostics
  }
  parameter_names <- c("mu", "tau", paste0(effect_name, "[", input$subject, "]"))
  dimnames(chain_array) <- list(NULL, paste0("chain", seq_len(chains)), parameter_names)
  mu_chain <- chain_array[, , 1L, drop = TRUE]
  tau_chain <- chain_array[, , 2L, drop = TRUE]
  effect_array <- chain_array[, , seq_len(length(effect_hat)) + 2L, drop = FALSE]
  mu_draws <- as.vector(mu_chain)
  tau_draws <- as.vector(tau_chain)
  effect_draws <- matrix(effect_array, nrow = ndraws * chains,
                         ncol = length(effect_hat))
  parameter_draws <- cbind(mu = mu_draws, tau = tau_draws, effect_draws)
  summary <- .sft_parameter_summary(parameter_draws, parameter_names, hdi, rope)
  population <- summary[summary$parameter %in% c("mu", "tau"), , drop = FALSE]
  subject_summary <- summary[grepl(paste0("^", effect_name, "\\["), summary$parameter), , drop = FALSE]
  rownames(population) <- NULL
  rownames(subject_summary) <- input$subject
  posterior_predictive <- .sft_posterior_predictive(effect_draws, effect_hat, precision)
  shrinkage_draws <- .sft_shrinkage_draws(tau_draws, precision)
  shrinkage_summary <- .sft_shrinkage_summary(shrinkage_draws, input$subject, hdi)
  draws <- data.frame(mu = mu_draws, tau = tau_draws, effect_draws,
                      check.names = FALSE)
  names(draws)[-(1:2)] <- paste0(effect_name, "[", input$subject, "]")
  diagnostics <- .sft_mcmc_diagnostics(chain_array, parameter_names)
  prior <- if (method == "InvGamma") {
    list(mu = list(mean = prior_mean, sd = prior_sd),
         tau2 = list(family = "inverse-Gamma", shape = prior_shape,
                     rate = prior_rate), rope = rope)
  } else {
    list(mu = list(mean = prior_mean, sd = prior_sd),
         tau = list(family = "half-Normal", sd = prior_tau_sd), rope = rope)
  }
  model <- if (method == "InvGamma") {
    list(observation = sprintf("%s[i] ~ Normal(mean = %s[i], sd = 1 / sqrt(%s[i]))",
                               observed_name, effect_name, precision_name),
         subject = sprintf("%s[i] ~ Normal(mean = mu, sd = tau)", effect_name),
         population = "mu ~ Normal(mean = prior_mean, sd = prior_sd)",
         interpretation = if (score_method == "capacity") {
           "exp(eta[i]) is the participant weighted capacity ratio"
         } else {
           "theta[i] is the participant unstandardised score-effect parameter"
         },
         score_method = score_method,
         parameterization = "centered", sampler = "conjugate Gibbs")
  } else {
    list(observation = sprintf("%s[i] ~ Normal(mean = %s[i], sd = 1 / sqrt(%s[i]))",
                               observed_name, effect_name, precision_name),
         subject = sprintf("%s[i] ~ Normal(mean = mu, sd = tau)", effect_name),
         population = "mu ~ Normal(mean = prior_mean, sd = prior_sd)",
         interpretation = if (score_method == "capacity") {
           "exp(eta[i]) is the participant weighted capacity ratio"
         } else {
           "theta[i] is the participant unstandardised score-effect parameter"
         },
         score_method = score_method,
         tau = "tau ~ HalfNormal(prior_tau_sd)",
         parameterization = if (method == "Centered") "centered" else "non-centred",
         sampler = "Stan NUTS")
  }
  list(
    statistic = setNames(mean(mu_draws), "mu"),
    posterior_probability = list(
      population_super = mean(mu_draws > 0),
      population_limited = mean(mu_draws < 0),
      population_rope = if (is.null(rope)) NA_real_ else mean(abs(mu_draws) <= rope),
      subject_super = setNames(colMeans(effect_draws > 0), input$subject),
      subject_limited = setNames(colMeans(effect_draws < 0), input$subject),
      subject_rope = setNames(if (is.null(rope)) rep(NA_real_, length(effect_hat)) else
                                colMeans(abs(effect_draws) <= rope), input$subject)
    ),
    summary = summary, population_summary = population,
    subject_summary = subject_summary, score = scored$score, draws = draws,
    prior_predictive = prior_predictive,
    posterior_predictive = posterior_predictive,
    diagnostics = diagnostics,
    sampler_diagnostics = sampler_diagnostics,
    shrinkage_summary = shrinkage_summary,
    shrinkage_draws = shrinkage_draws,
    stan_fit = stan_fit,
    method = if (method == "InvGamma") {
      "Hierarchical Bayesian UCIP score/information model"
    } else {
      paste(method, "Bayesian UCIP score/information model")
    },
    method_code = method,
    score_method = score_method,
    alternative = paste("the population UCIP", stopping.rule, "effect is nonzero"),
    prior = prior, model = model,
    stopping.rule = stopping.rule, hdi = hdi, ndraws = ndraws,
    burnin = burnin, thin = thin, chains = chains, seed = seed
  )
}

resilience.test <- function(RT, CR = NULL, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (length(RT) < 3L) stop("RT must contain AB, A, and B conditions.")
  CR <- .sft_cr_list(RT[1:3], CR)
  times <- .sft_finite_times(RT[1:3])
  h <- lapply(1:3, function(i) estimateNAH(RT[[i]], CR[[i]]))
  contrast <- h[[1]]$H(times) - h[[2]]$H(times) - h[[3]]$H(times)
  v <- h[[1]]$Var(times) + h[[2]]$Var(times) + h[[3]]$Var(times)
  z <- contrast[length(contrast)] / sqrt(v[length(v)])
  p <- if (is.finite(z)) 2 * min(stats::pnorm(z), 1 - stats::pnorm(z)) else NA_real_
  .sft_htest(setNames(z, "z"), p,
             "response times are different than those predicted by the UCIP-OR model",
             "Resilience test", deparse(substitute(RT)))
}

resilience <- function(RT, CR = NULL, ratio = TRUE, Condition = NULL, Subject = NULL) {
  converted <- .sft_as_rt_cr(RT, CR, stopping.rule = "OR",
                             Condition = Condition, Subject = Subject)
  RT <- converted$RT; CR <- converted$CR
  if (length(RT) < 3L) stop("RT must contain AB, A, and B conditions.")
  CR <- .sft_cr_list(RT[1:3], CR)
  times <- .sft_finite_times(RT[1:3])
  h <- lapply(1:3, function(i) estimateNAH(RT[[i]], CR[[i]]))
  if (ratio) {
    r <- h[[1]]$H(times) / (h[[2]]$H(times) + h[[3]]$H(times))
    r[!is.finite(r)] <- NA_real_
    list(Rt = .sft_curve(times, r), Rtest = resilience.test(RT, CR))
  } else {
    r <- h[[1]]$H(times) - h[[2]]$H(times) - h[[3]]$H(times)
    v <- h[[1]]$Var(times) + h[[2]]$Var(times) + h[[3]]$Var(times)
    list(Rt = .sft_curve(times, r), Var = .sft_curve(times, v),
         Rtest = resilience.test(RT, CR))
  }
}

capacityGroup <- function(inData, acc.cutoff = .9, ratio = TRUE, OR = NULL,
                          stopping.rule = c("OR", "AND", "STST"), plotCt = TRUE, ...) {
  required <- c("Subject", "Condition", "RT", "Correct")
  if (!is.data.frame(inData) || !all(required %in% names(inData))) {
    stop("inData must contain Subject, Condition, RT, and Correct columns.")
  }
  channels <- grep("^Channel", names(inData), value = TRUE)
  if (length(channels) < 2L) stop("Not enough channels for capacity analysis.")
  rule <- if (!is.null(OR)) if (isTRUE(OR)) "OR" else "AND" else match.arg(stopping.rule)
  subjects <- unique(inData$Subject)
  conditions <- unique(inData$Condition)
  rt <- as.numeric(inData$RT)
  finite_rt <- rt[is.finite(rt)]
  if (!length(finite_rt)) stop("inData$RT contains no finite response times.")
  times <- seq(stats::quantile(finite_rt, .001, names = FALSE),
               stats::quantile(finite_rt, .999, names = FALSE), length.out = 1000)
  fits <- list(); curves <- list(); rows <- list(); z_by_condition <- list(); n <- 0L

  for (cond in conditions) {
    z <- numeric()
    for (subj in subjects) {
      d <- inData[inData$Condition == cond & inData$Subject == subj, , drop = FALSE]
      if (!nrow(d)) next
      if (rule == "STST") {
        context <- rowSums(d[channels] > 0) == 1L & rowSums(d[channels] < 0) > 0L
        target <- rowSums(d[channels] >= 0) == length(channels) & rowSums(d[channels] != 0) == 1L
        rts <- list(d$RT[context], d$RT[target])
        cr <- list(d$Correct[context], d$Correct[target])
      } else {
        context <- rowSums(d[channels] > 0) == length(channels)
        singles <- lapply(seq_along(channels), function(j) {
          other <- setdiff(seq_along(channels), j)
          d[["RT"]][d[[channels[j]]] > 0 & rowSums(d[channels[other]] == 0) == length(other)]
        })
        single_cr <- lapply(seq_along(channels), function(j) {
          other <- setdiff(seq_along(channels), j)
          d[["Correct"]][d[[channels[j]]] > 0 & rowSums(d[channels[other]] == 0) == length(other)]
        })
        rts <- c(list(d$RT[context]), singles)
        cr <- c(list(d$Correct[context]), single_cr)
      }
      accuracy <- vapply(cr, function(x) {
        x <- as.logical(x)
        if (!length(x)) return(FALSE)
        mean(x, na.rm = TRUE) >= acc.cutoff && sum(x, na.rm = TRUE) >= 10L
      }, logical(1))
      fit <- if (all(accuracy)) {
        if (rule == "OR") capacity.or(rts, cr, ratio)
        else if (rule == "AND") capacity.and(rts, cr, ratio)
        else capacity.stst(rts, cr, ratio)
      } else NULL

      n <- n + 1L
      fits[[n]] <- fit
      curves[[n]] <- if (is.null(fit)) rep(NA_real_, length(times)) else fit$Ct(times)
      if (!is.null(fit) && !is.null(fit$Ctest)) z <- c(z, unname(fit$Ctest$statistic))
      label <- if (is.null(fit)) NA_character_ else if (!is.null(fit$Ctest) &&
          is.finite(fit$Ctest$p.value) && fit$Ctest$p.value < .05) {
        if (unname(fit$Ctest$statistic) < 0) "Limited" else "Super"
      } else "Nonsignificant"
      rows[[n]] <- data.frame(Subject = as.character(subj), Condition = as.character(cond),
                              Capacity = label, stringsAsFactors = FALSE)
      if (plotCt && !is.null(fit)) {
        plot(times, fit$Ct(times), type = "l", xlab = "Time", ylab = "C(t)",
             main = paste(cond, "\nParticipant", subj), ...)
        abline(if (ratio) 1 else 0, 0, lty = 2)
      }
    }
    z_by_condition[[as.character(cond)]] <- z
    if (length(z) > 2L) {
      p <- stats::t.test(z)$p.value
      group_label <- if (p < .05) if (mean(z) < 0) "Limited" else "Super" else "Nonsignificant"
    } else group_label <- "Unknown"
    n <- n + 1L
    fits[[n]] <- list()
    curves[[n]] <- rep(NA_real_, length(times))
    rows[[n]] <- data.frame(Subject = "Group", Condition = as.character(cond),
                            Capacity = group_label, stringsAsFactors = FALSE)
  }
  overview <- if (length(rows)) do.call(rbind, rows) else data.frame()
  curve_mat <- if (length(curves)) do.call(rbind, curves) else matrix(numeric(), 0, length(times))
  list(overview = overview, Ct.fn = curve_mat, capacity = fits,
       times = times, z = z_by_condition, rule = rule)
}

sic.test <- function(HH, HL, LH, LL, method = "ks") {
  method <- match.arg(method, "ks")
  xs <- list(HH = HH, HL = HL, LH = LH, LL = LL)
  if (any(vapply(xs, function(x) !any(is.finite(x)), logical(1)))) {
    stop("All four SIC cells must contain at least one finite RT.")
  }
  times <- sort(unique(unlist(lapply(xs, function(x) x[is.finite(x)]), use.names = FALSE)))
  sic_values <- .sft_ecdf(LH)(times) + .sft_ecdf(HL)(times) -
    .sft_ecdf(HH)(times) - .sft_ecdf(LL)(times)
  n_eff <- 1 / sum(1 / vapply(xs, function(x) sum(is.finite(x)), numeric(1)))
  dplus <- max(0, sic_values)
  dminus <- abs(min(0, sic_values))
  pplus <- exp(-2 * n_eff * dplus^2)
  pminus <- exp(-2 * n_eff * dminus^2)
  dname <- paste("HH:", deparse(substitute(HH)), "HL:", deparse(substitute(HL)),
                 "LH:", deparse(substitute(LH)), "LL:", deparse(substitute(LL)))
  plus <- .sft_htest(setNames(dplus, "D+"), pplus, "the SIC is above 0 at some time",
                     "Houpt-Townsend KS-SIC test", dname)
  minus <- .sft_htest(setNames(dminus, "D-"), pminus, "the SIC is below 0 at some time",
                      "Houpt-Townsend KS-SIC test", dname)
  list(positive = plus, negative = minus)
}

mic.test <- function(HH, HL, LH, LL, method = c("art", "anova")) {
  method <- match.arg(method)
  xs <- lapply(list(HH, HL, LH, LL), function(x) as.numeric(x[is.finite(x)]))
  if (any(!lengths(xs))) stop("All four SIC cells must contain finite RTs.")
  statistic <- (mean(xs[[4]]) - mean(xs[[3]])) - (mean(xs[[2]]) - mean(xs[[1]]))
  names(statistic) <- "MIC"
  allrt <- unlist(xs, use.names = FALSE)
  h1 <- c(rep(1, length(xs[[1]]) + length(xs[[2]])), rep(0, length(xs[[3]]) + length(xs[[4]])))
  h2 <- c(rep(1, length(xs[[1]])), rep(0, length(xs[[2]])), rep(1, length(xs[[3]])), rep(0, length(xs[[4]])))
  if (method == "art") {
    mA0 <- mean(allrt[h1 == 0]); mA1 <- mean(allrt[h1 == 1])
    mB0 <- mean(allrt[h2 == 0]); mB1 <- mean(allrt[h2 == 1])
    adjusted <- round(allrt - (1 - h1) * mA0 - h1 * mA1 - (1 - h2) * mB0 - h2 * mB1, 15)
    fit <- stats::anova(stats::lm(rank(adjusted) ~ h1 * h2))
    p <- fit[["Pr(>F)"]][3]
    method_name <- "Adjusted Rank Transform test of the MIC"
  } else {
    old <- options(contrasts = c("contr.helmert", "contr.poly")); on.exit(options(old), add = TRUE)
    fit <- stats::anova(stats::lm(allrt ~ h1 * h2))
    p <- fit[["Pr(>F)"]][3]
    method_name <- "ANOVA test of the MIC"
  }
  .sft_htest(statistic, p, "the MIC is not zero", method_name,
             paste("HH:", deparse(substitute(HH)), "HL:", deparse(substitute(HL)),
                   "LH:", deparse(substitute(LH)), "LL:", deparse(substitute(LL))))
}

.sft_bayes_pair <- function(x, y, nbin = 20L, nsamp = 10000L, tol = 0, seed = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("Both distributions must contain finite RTs.")
  pooled <- c(x, y)
  nbin <- max(1L, min(as.integer(nbin), length(unique(pooled))))
  probs <- seq(0, 1, length.out = nbin + 1L)
  breaks <- unique(stats::quantile(pooled, probs = probs, names = FALSE, type = 8))
  if (length(breaks) < 2L) breaks <- c(breaks[[1L]], breaks[[1L]] + 1)
  breaks[length(breaks)] <- max(breaks[length(breaks)], max(pooled) + .Machine$double.eps)
  nbin <- length(breaks) - 1L
  alpha <- rep(1 / nbin, nbin)
  counts <- c(hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts,
              hist(y, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts)
  cx <- counts[seq_len(nbin)]; cy <- counts[nbin + seq_len(nbin)]
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  rdir <- function(n, a) {
    g <- matrix(stats::rgamma(n * length(a), shape = rep(a, each = n)), nrow = n)
    g / rowSums(g)
  }
  prior_x <- rdir(nsamp, alpha); prior_y <- rdir(nsamp, alpha)
  post_x <- rdir(nsamp, cx + alpha); post_y <- rdir(nsamp, cy + alpha)
  dominates <- function(a, b) rowSums(t(apply(a - b, 1, cumsum)) >= -tol) == ncol(a)
  prior_greater <- mean(dominates(prior_x, prior_y))
  post_greater <- mean(dominates(post_x, post_y))
  prior_less <- mean(dominates(prior_y, prior_x))
  post_less <- mean(dominates(post_y, post_x))
  list(BF_greater = post_greater / max(prior_greater, .Machine$double.eps),
       BF_less = post_less / max(prior_less, .Machine$double.eps),
       posterior_greater = post_greater, posterior_less = post_less,
       prior_greater = prior_greater, prior_less = prior_less,
       nbin = nbin, nsamp = nsamp)
}

siDominance <- function(HH, HL, LH, LL, method = c("ks", "dp"),
                        nbin = 20L, nsamp = 10000L, seed = NULL) {
  method <- match.arg(method)
  pairs <- list(
    list(greater = "S.hh > S.hl", less = "S.hh < S.hl", x = HH, y = HL),
    list(greater = "S.hh > S.lh", less = "S.hh < S.lh", x = HH, y = LH),
    list(greater = "S.hl > S.ll", less = "S.hl < S.ll", x = HL, y = LL),
    list(greater = "S.lh > S.ll", less = "S.lh < S.ll", x = LH, y = LL)
  )
  if (method == "ks") {
    rows <- lapply(pairs, function(p) {
      greater <- suppressWarnings(stats::ks.test(p$x, p$y, alternative = "greater", exact = FALSE))
      less <- suppressWarnings(stats::ks.test(p$x, p$y, alternative = "less", exact = FALSE))
      list(test = c(p$greater, p$less), statistic = c(unname(greater$statistic), unname(less$statistic)),
           p = c(greater$p.value, less$p.value))
    })
    out <- data.frame(Test = unlist(lapply(rows, `[[`, "test")),
                      statistic = unlist(lapply(rows, `[[`, "statistic")),
                      p.value = unlist(lapply(rows, `[[`, "p")),
                      row.names = NULL)
    return(out)
  }
  rows <- lapply(seq_along(pairs), function(i) {
    p <- pairs[[i]]
    bf <- .sft_bayes_pair(p$x, p$y, nbin = nbin, nsamp = nsamp,
                          seed = if (is.null(seed)) NULL else seed + i)
    data.frame(Test = c(p$greater, p$less), BF = c(bf$BF_greater, bf$BF_less),
               posterior = c(bf$posterior_greater, bf$posterior_less),
               prior = c(bf$prior_greater, bf$prior_less),
               p.value = c(bf$posterior_greater, bf$posterior_less),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "inferential_type") <- "Bayesian posterior probability and posterior/prior event-probability ratio"
  out
}

.sft_dp_checkmods <- function(x, dx, tolSIC = 5e-2, tolMIC = 1e-2) {
  ans <- c(Z = FALSE, N = FALSE, P = FALSE, nP = FALSE, Np = FALSE, np = FALSE)
  nonzero <- abs(x) > tolSIC
  if (!any(nonzero)) { ans["Z"] <- TRUE; return(ans) }
  if (all(x <= 0) && any(x < 0)) { ans["N"] <- TRUE; return(ans) }
  if (all(x >= 0) && any(x > 0)) { ans["P"] <- TRUE; return(ans) }
  signed <- x[nonzero]
  if (signed[[1L]] > 0 || sum(diff(sign(signed)) != 0) > 1L) return(ans)
  neg_area <- -sum(x[x < -tolSIC] * dx[x < -tolSIC])
  pos_area <- sum(x[x > tolSIC] * dx[x > tolSIC])
  if (abs(neg_area - pos_area) < tolMIC) ans["np"] <- TRUE
  else if (neg_area < pos_area) ans["nP"] <- TRUE
  else ans["Np"] <- TRUE
  ans
}

sicDPtest <- function(dat, nbin = NULL, nsamp = 10000L, maxn = 500000L,
                      tolSIC = 5e-2, tolMIC = NULL, ci = .95, seed = NULL) {
  if (!is.list(dat) || length(dat) != 4L) stop("dat must be a list in HH, HL, LH, LL order.")
  dat <- lapply(dat, function(x) as.numeric(x[is.finite(x)]))
  if (any(!lengths(dat))) stop("All SIC cells must contain finite RTs.")
  pooled <- unlist(dat, use.names = FALSE)
  if (is.null(nbin)) nbin <- floor(min(lengths(dat)) / 2)
  nbin <- max(1L, min(as.integer(nbin), length(unique(pooled))))
  q <- stats::quantile(pooled, probs = seq(0, 1, length.out = nbin + 1L), names = FALSE, type = 8)
  bins <- unique(q)
  if (length(bins) < 2L) bins <- c(bins[[1L]], bins[[1L]] + 1)
  bins[length(bins)] <- max(bins[length(bins)], max(pooled) + .Machine$double.eps)
  nbin <- length(bins) - 1L
  dx <- diff(bins)
  counts <- sapply(dat, function(x) hist(x, breaks = bins, plot = FALSE, include.lowest = TRUE)$counts)
  alpha <- rep(1 / nbin, nbin)
  if (is.null(tolMIC)) tolMIC <- max(mean(pooled) / 700, 1e-8)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  rdir <- function(n, a) {
    g <- matrix(stats::rgamma(n * length(a), shape = rep(a, each = n)), nrow = n)
    g / rowSums(g)
  }
  classify <- function(n, a) {
    p <- rdir(n, a)
    cdfs <- t(apply(p, 1, cumsum))
    apply(cdfs, 1, .sft_dp_checkmods, dx = dx, tolSIC = tolSIC, tolMIC = tolMIC)
  }
  prior <- matrix(0, nrow = 6, ncol = 1, dimnames = list(c("Z", "N", "P", "nP", "Np", "np"), "Prior"))
  post <- prior; N <- 0L
  repeat {
    prior <- prior + rowSums(classify(nsamp, alpha))
    post <- post + rowSums(classify(nsamp, counts + alpha))
    N <- N + nsamp
    bf <- (post + 1) / (prior + 1)
    if (N >= maxn) break
    # The highest-BF class is the only one whose interval width is required for
    # the stopping rule, matching the original implementation's fast default.
    p_hi <- stats::qbeta(c((1 - ci) / 2, 1 - (1 - ci) / 2), post[which.max(bf)] + 1,
                         N - post[which.max(bf)] + 1)
    if (diff(p_hi) * N < 0.01 * N) break
  }
  bf <- drop((post + 1) / (prior + 1)); names(bf) <- rownames(prior)
  bfp <- bf / sum(bf)
  list(BF = bf, BFp = bfp, N = rbind(Prior = drop(prior), Posterior = drop(post)),
       ci = ci, nbin = nbin, nsamp = nsamp, seed = seed,
       tolerances = c(SIC = tolSIC, MIC = tolMIC))
}

sictestBayes <- function(HH, HL, LH, LL, method = c("DP", "IG"), model = NULL,
                         nbin = NULL, nsamp = 10000L, maxn = 500000L,
                         seed = NULL, ...) {
  method <- match.arg(method)
  if (method == "IG") {
    stop("The legacy IG/Stan SIC route is not bundled in sft_plus; use method='DP' or supply a separate model implementation.")
  }
  result <- sicDPtest(list(HH, HL, LH, LL), nbin = nbin, nsamp = nsamp,
                      maxn = maxn, seed = seed, ...)
  bf <- result$BF
  names(bf)[names(bf) == "Z"] <- "Zero"
  names(bf)[names(bf) == "N"] <- "Negative"
  names(bf)[names(bf) == "P"] <- "Positive"
  names(bf)[names(bf) == "np"] <- "NegPos.MIC0"
  names(bf)[names(bf) == "nP"] <- "NegPos.MICpos"
  names(bf)[names(bf) == "Np"] <- "NegPos.MICneg"
  list(statistic = bf, BFp = result$BFp, method = "Nonparametric Bayesian SIC test",
       data.name = paste("HH:", deparse(substitute(HH)), "HL:", deparse(substitute(HL)),
                         "LH:", deparse(substitute(LH)), "LL:", deparse(substitute(LL))),
       details = result)
}

sic <- function(HH, HL, LH, LL, domtest = "ks", sictest = "ks",
                mictest = c("art", "anova"), interpolate = FALSE,
                nbin = 20L, nsamp = 10000L, maxn = 500000L, seed = NULL, ...) {
  xs <- list(HH = HH, HL = HL, LH = LH, LL = LL)
  if (any(!vapply(xs, function(x) any(is.finite(x)), logical(1)))) stop("All SIC cells need finite RTs.")
  times <- sort(unique(unlist(lapply(xs, function(x) x[is.finite(x)]), use.names = FALSE)))
  sic_values <- .sft_ecdf(LH)(times) + .sft_ecdf(HL)(times) - .sft_ecdf(HH)(times) - .sft_ecdf(LL)(times)
  sic_fun <- if (interpolate) .sft_curve(times, sic_values) else stats::stepfun(times, c(0, sic_values))
  d <- siDominance(HH, HL, LH, LL, method = domtest, nbin = nbin, nsamp = nsamp, seed = seed)
  st <- if (sictest == "bf") sictestBayes(HH, HL, LH, LL, nbin = nbin, nsamp = nsamp, maxn = maxn, seed = seed) else sic.test(HH, HL, LH, LL, method = sictest)
  mt <- mic.test(HH, HL, LH, LL, method = mictest)
  n_eff <- 1 / sum(1 / vapply(xs, function(x) sum(is.finite(x)), numeric(1)))
  dvals <- if (sictest == "bf") NULL else rbind(c(unname(st$positive$statistic), st$positive$p.value),
                                                  c(unname(st$negative$statistic), st$negative$p.value))
  list(SIC = sic_fun, MIC = unname(mt$statistic), MICtest = mt,
       SICtest = st, Dominance = d, Dvals = dvals, N = n_eff,
       HH = .sft_ecdf(HH), HL = .sft_ecdf(HL), LH = .sft_ecdf(LH), LL = .sft_ecdf(LL))
}

sicGroup <- function(inData, sictest = c("ks", "bf"), mictest = c("art", "anova"),
                     domtest = c("ks", "dp"), alpha.sic = .05, plotSIC = TRUE,
                     ...) {
  sictest <- match.arg(sictest); domtest <- match.arg(domtest); mictest <- match.arg(mictest)
  req <- c("Subject", "Condition", "RT", "Correct", "Channel1", "Channel2")
  if (!all(req %in% names(inData))) stop("inData is missing required SIC columns.")
  times <- sort(unique(round(inData$RT[is.finite(inData$RT)])))
  records <- list(); curves <- list(); fits <- list(); n <- 0L
  for (cond in unique(inData$Condition)) for (subj in unique(inData$Subject)) {
    d <- inData[inData$Condition == cond & inData$Subject == subj & inData$Correct == 1, , drop = FALSE]
    cell <- list(d$RT[d$Channel1 == 2 & d$Channel2 == 2], d$RT[d$Channel1 == 2 & d$Channel2 == 1],
                  d$RT[d$Channel1 == 1 & d$Channel2 == 2], d$RT[d$Channel1 == 1 & d$Channel2 == 1])
    if (min(lengths(cell)) <= 10L) next
    n <- n + 1L; fit <- sic(cell[[1]], cell[[2]], cell[[3]], cell[[4]],
                             domtest = domtest, sictest = sictest, mictest = mictest, ...)
    fits[[n]] <- fit; curves[[n]] <- fit$SIC(times)
    records[[n]] <- data.frame(Subject = as.character(subj), Condition = as.character(cond),
                                stringsAsFactors = FALSE)
  }
  overview <- if (!length(records)) data.frame() else do.call(rbind, records)
  if (nrow(overview)) {
    if (sictest == "bf") {
      model <- vapply(fits, function(f) names(f$SICtest$statistic)[which.max(f$SICtest$statistic)], character(1))
      overview$Model <- model
    } else {
      positive <- vapply(fits, function(f) f$SICtest$positive$p.value < alpha.sic, logical(1))
      negative <- vapply(fits, function(f) f$SICtest$negative$p.value < alpha.sic, logical(1))
      overview$Positive.SIC <- ifelse(positive, "Significant", "Nonsignificant")
      overview$Negative.SIC <- ifelse(negative, "Significant", "Nonsignificant")
      overview$MIC <- ifelse(vapply(fits, function(f) f$MICtest$p.value < alpha.sic, logical(1)), "Significant", "Nonsignificant")
    }
  }
  list(overview = overview, sic.fn = if (length(curves)) do.call(rbind, curves) else matrix(numeric(), 0, length(times)),
       sic = fits, times = times)
}

sicGroupBF <- function(inData, domtest = "ks", plotSIC = TRUE, ...) {
  sicGroup(inData, sictest = "bf", domtest = domtest, plotSIC = plotSIC, ...)
}

.to_correct_indicator <- function(cr, n) {
  if (is.null(cr)) return(rep(TRUE, n))
  out <- rep(FALSE, n); m <- min(length(cr), n)
  if (m) {
    v <- as.logical(cr[seq_len(m)]); v[is.na(v)] <- FALSE; out[seq_len(m)] <- v
  }
  out
}

.sft_bind_rows <- function(parts) {
  parts <- parts[vapply(parts, is.data.frame, logical(1))]
  if (!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(cols, names(x)); x[missing] <- NA
    x[cols]
  })
  do.call(rbind, parts)
}

build_ct_df <- function(times, series_specs, probs = c(.05, .95)) {
  .sft_bind_rows(lapply(series_specs, function(spec) {
    values <- spec$fn(times)
    if (!is.null(spec$rt) && length(spec$rt)) {
      q <- stats::quantile(spec$rt, probs = probs, names = FALSE, na.rm = TRUE)
      values[times < q[1] | times > q[2]] <- NA_real_
    }
    data.frame(Time = times, Ct = values, Series = rep(spec$label, length(times)),
               stringsAsFactors = FALSE)
  }))
}

smooth_one_cdf <- function(t, y, smooth_cdf = "none", smooth_spar = 0.65) {
  smooth_cdf <- match.arg(smooth_cdf, c("none", "mono", "spline"))
  if (smooth_cdf == "none") return(y)
  if (!length(t) || !length(y)) return(y)
  ord <- order(t); to <- t[ord]; yo <- y[ord]
  keep <- is.finite(to) & is.finite(yo)
  out <- y; if (!any(keep)) return(out)
  tu <- unique(to[keep]); yu <- yo[match(tu, to)]
  yu <- pmin(1, pmax(0, cummax(yu)))
  if (length(tu) < 3L) ys <- yo else if (smooth_cdf == "mono") {
    ys <- stats::splinefun(tu, yu, method = "monoH.FC")(to)
  } else {
    spar <- max(.4, min(1, if (is.finite(smooth_spar)) smooth_spar else .65))
    fit <- tryCatch(stats::smooth.spline(to, yo, spar = spar), error = function(e) NULL)
    ys <- if (is.null(fit)) stats::splinefun(tu, yu, method = "monoH.FC")(to) else stats::predict(fit, x = to)$y
  }
  ys <- cummax(pmin(1, pmax(0, ys))); out[ord] <- ys; out
}

build_sic_df <- function(times, series_specs, trim_cdf_tails = TRUE,
                         cdf_tail_cut = 5e-4,
                         smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65) {
  smooth_cdf <- match.arg(smooth_cdf)
  .sft_bind_rows(lapply(series_specs, function(spec) {
    raw <- lapply(c("HH", "HL", "LH", "LL"), function(nm) spec$fn[[nm]](times))
    names(raw) <- c("HH", "HL", "LH", "LL")
    sm <- lapply(raw, function(y) smooth_one_cdf(times, y, smooth_cdf = smooth_cdf, smooth_spar = smooth_spar))
    names(sm) <- names(raw)
    sic_raw <- raw$LH + raw$HL - raw$HH - raw$LL
    sic_smooth <- sm$LH + sm$HL - sm$HH - sm$LL
    out <- data.frame(Time = times, SICt = if (smooth_cdf == "none") sic_raw else sic_smooth,
                      SICt_raw = sic_raw, SICt_smooth = sic_smooth,
                      HH_s = 1 - sm$HH, HL_s = 1 - sm$HL, LH_s = 1 - sm$LH, LL_s = 1 - sm$LL,
                      Series = rep(spec$label, length(times)), stringsAsFactors = FALSE)
    if (isTRUE(trim_cdf_tails)) {
      cut <- if (is.finite(cdf_tail_cut) && cdf_tail_cut > 0 && cdf_tail_cut < .5) cdf_tail_cut else 5e-4
      left <- sm$HH < cut & sm$HL < cut & sm$LH < cut & sm$LL < cut
      right <- sm$HH > 1 - cut & sm$HL > 1 - cut & sm$LH > 1 - cut & sm$LL > 1 - cut
      out <- out[!(left | right), , drop = FALSE]
    }
    out
  }))
}

build_at_df <- function(times, series_specs, probs = c(.05, .95)) {
  .sft_bind_rows(lapply(series_specs, function(spec) {
    a <- spec$fn
    values <- lapply(c("A_CF", "A_IF", "A_CS", "A_IS"), function(nm) a[[nm]](times))
    names(values) <- c("A_CF", "A_IF", "A_CS", "A_IS")
    if (!is.null(spec$rt) && length(spec$rt)) {
      q <- stats::quantile(spec$rt, probs = probs, names = FALSE, na.rm = TRUE)
      values <- lapply(values, function(x) { x[times < q[1] | times > q[2]] <- NA_real_; x })
    }
    data.frame(Time = times, as.data.frame(values), Series = rep(spec$label, length(times)),
               stringsAsFactors = FALSE)
  }))
}

build_cdf_df <- function(times, series_specs, probs = c(.05, .95),
                         smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65) {
  smooth_cdf <- match.arg(smooth_cdf)
  .sft_bind_rows(lapply(series_specs, function(spec) {
    if (is.null(spec$rt)) stop("Each series spec must include an `rt` vector.")
    keep <- .to_correct_indicator(spec$cr, length(spec$rt))
    rt <- as.numeric(spec$rt)[keep & is.finite(spec$rt)]
    raw <- .sft_ecdf(rt)(times)
    sm <- smooth_one_cdf(times, raw, smooth_cdf, smooth_spar)
    if (length(rt)) {
      q <- stats::quantile(rt, probs = probs, names = FALSE)
      raw[times < q[1] | times > q[2]] <- NA_real_
      sm[times < q[1] | times > q[2]] <- NA_real_
    } else raw[] <- sm[] <- NA_real_
    data.frame(Time = times, CDF = if (smooth_cdf == "none") raw else sm,
               CDF_raw = raw, CDF_smooth = sm, Series = rep(spec$label, length(times)),
               stringsAsFactors = FALSE)
  }))
}

get_time_range <- function(rt_list, probs = c(.05, .95), by = .001) {
  if (!is.list(rt_list)) rt_list <- list(rt_list)
  q <- lapply(rt_list, function(x) stats::quantile(x[is.finite(x)], probs = probs, names = FALSE, na.rm = TRUE))
  if (!length(q) || any(vapply(q, length, integer(1)) < 2L)) stop("No usable RT quantiles.")
  lo <- max(vapply(q, `[[`, numeric(1), 1)); hi <- min(vapply(q, `[[`, numeric(1), 2))
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) stop("No usable quantile window for the provided RTs.")
  seq(lo, hi, by = by)
}

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

plot_altieri <- function(n, p_vec, smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_altieri() requires ggplot2.")
  smooth_cdf <- match.arg(smooth_cdf)
  sims <- simulate_sft(model = "lba", n = n, p_vec = p_vec, design = c("AB", "AN", "NB", "NN"), logical_rules = "AND")
  d <- sims$by_rule$AND; cells <- c("AB", "AN", "NB", "NN")
  rts <- lapply(cells, function(s) d$RT[d$S == s & d$Correct == 1])
  times <- get_time_range(rts)
  ct <- build_ct_df(times, list(
    list(label = "C_AND(t)", fn = sims$metrics$AND_Ct$Ct, rt = d$RT),
    list(label = "C_Absence(t)", fn = sims$metrics$AND_Absence_Ct$Ct, rt = d$RT),
    list(label = "C_Altieri(t)", fn = sims$metrics$Altieri_Ct$Ct_a1, rt = d$RT),
    list(label = "C_Bound(t)", fn = sims$metrics$Altieri_Ct$Ct_a2, rt = d$RT)))
  cdf <- build_cdf_df(times, list(list(label = "AB", rt = d$RT[d$S == "AB"], cr = d$Correct[d$S == "AB"]),
                                  list(label = "NN", rt = d$RT[d$S == "NN"], cr = d$Correct[d$S == "NN"])),
                      smooth_cdf = smooth_cdf, smooth_spar = smooth_spar)
  p1 <- ggplot2::ggplot(ct, ggplot2::aes_string(x = "Time", y = "Ct", color = "Series", linetype = "Series")) + ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") + ggplot2::theme_classic()
  p2 <- ggplot2::ggplot(cdf, ggplot2::aes_string(x = "Time", y = "CDF", color = "Series")) + ggplot2::geom_line() + ggplot2::theme_classic()
  list(sims = sims, plot_df = ct, cdf_df = cdf, plot = p1, cdf_plot = p2, p_vec = p_vec)
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

.at_clean <- function(RT, CR = NULL) {
  RT <- as.numeric(RT); CR <- .sft_correct(CR, length(RT)); keep <- is.finite(RT)
  list(RT = RT[keep], CR = CR[keep])
}
.at_count <- function(x, t) if (length(x)) findInterval(t, sort(x), rightmost.closed = TRUE) else rep(0, length(t))
.at_defective <- function(RT, CR = NULL) {
  z <- .at_clean(RT, CR); n <- length(z$RT); if (!n) stop("No finite RTs.")
  c_rt <- sort(z$RT[z$CR]); i_rt <- sort(z$RT[!z$CR]); pc <- length(c_rt) / n; pi <- length(i_rt) / n
  cdf <- function(x) function(t) .at_count(x, t) / n
  list(n = n, pC = pc, pI = pi, CF = cdf(c_rt), IF = cdf(i_rt),
       CS = function(t) pmax(0, pmin(1, pc - cdf(c_rt)(t))),
       IS = function(t) pmax(0, pmin(1, pi - cdf(i_rt)(t))))
}
.at_log_ratio <- function(num, den) {
  out <- log(pmax(num, 0)) / log(pmax(den, 0)); out[!is.finite(out) | num <= 0 | den <= 0] <- NA_real_; out
}
.at_pack <- function(times, values) {
  values <- as.data.frame(values); funs <- lapply(values[setdiff(names(values), "time")], function(x) .sft_curve(times, x))
  c(funs, list(values = values, times = times))
}

assessment_ta_or <- function(RT, CR = NULL, times = NULL) {
  if (length(RT) < 3L) stop("RT must be list(AB, A, B).")
  if (is.null(CR) || length(CR) != length(RT)) CR <- lapply(RT, function(x) rep(TRUE, length(x)))
  if (is.null(times)) times <- .sft_finite_times(RT[1:3])
  ab <- .at_defective(RT[[1]], CR[[1]]); a <- .at_defective(RT[[2]], CR[[2]]); b <- .at_defective(RT[[3]], CR[[3]])
  cfA <- a$CF(times); csA <- a$CS(times); ifA <- a$IF(times); isA <- a$IS(times)
  cfB <- b$CF(times); csB <- b$CS(times); ifB <- b$IF(times); isB <- b$IS(times)
  pred_IF <- ifA * ifB
  pred_CF <- cfA * b$pI + cfB * a$pI + cfA * csB + cfB * csA + cfA * cfB
  pred_IS <- pmax(0, pmin(1, isA * b$pI + isB * a$pI - isA * isB))
  pred_CS <- csA * b$pI + csB * a$pI + csA * csB
  obs <- list(IF = ab$IF(times), CF = ab$CF(times), IS = ab$IS(times), CS = ab$CS(times))
  .at_pack(times, list(time = times, A_CF = .at_log_ratio(pred_CF, obs$CF), A_IF = .at_log_ratio(pred_IF, obs$IF),
                       A_CS = .at_log_ratio(pred_CS, obs$CS), A_IS = .at_log_ratio(pred_IS, obs$IS),
                       pred_CF = pred_CF, obs_CF = obs$CF, pred_IF = pred_IF, obs_IF = obs$IF,
                       pred_CS = pred_CS, obs_CS = obs$CS, pred_IS = pred_IS, obs_IS = obs$IS))
}

assessment_ta_and <- function(RT, CR = NULL, times = NULL) {
  if (length(RT) < 3L) stop("RT must be list(AB, A, B).")
  if (is.null(CR) || length(CR) != length(RT)) CR <- lapply(RT, function(x) rep(TRUE, length(x)))
  if (is.null(times)) times <- .sft_finite_times(RT[1:3])
  ab <- .at_defective(RT[[1]], CR[[1]]); a <- .at_defective(RT[[2]], CR[[2]]); b <- .at_defective(RT[[3]], CR[[3]])
  cfA <- a$CF(times); csA <- a$CS(times); ifA <- a$IF(times); isA <- a$IS(times)
  cfB <- b$CF(times); csB <- b$CS(times); ifB <- b$IF(times); isB <- b$IS(times)
  pred_IF <- ifA * b$pC + ifB * a$pC + ifA * isB + ifB * isA + ifA * ifB
  pred_CF <- cfA * cfB
  pred_IS <- isA * b$pC + isB * a$pC + isA * isB
  pred_CS <- pmax(0, pmin(1, csA * b$pC + csB * a$pC - csA * csB))
  obs <- list(IF = ab$IF(times), CF = ab$CF(times), IS = ab$IS(times), CS = ab$CS(times))
  .at_pack(times, list(time = times, A_CF = .at_log_ratio(pred_CF, obs$CF), A_IF = .at_log_ratio(pred_IF, obs$IF),
                       A_CS = .at_log_ratio(pred_CS, obs$CS), A_IS = .at_log_ratio(pred_IS, obs$IS),
                       pred_CF = pred_CF, obs_CF = obs$CF, pred_IF = pred_IF, obs_IF = obs$IF,
                       pred_CS = pred_CS, obs_CS = obs$CS, pred_IS = pred_IS, obs_IS = obs$IS))
}

.at_weighted_survivor_integral <- function(event_rt, other_rt, n_self) {
  event_rt <- sort(as.numeric(event_rt[is.finite(event_rt)])); other_rt <- sort(as.numeric(other_rt[is.finite(other_rt)]))
  if (!length(event_rt) || n_self <= 0 || !length(other_rt)) {
    zero <- function(t) rep(0, length(t)); return(list(fast = zero, slow = zero, total = 0))
  }
  weights <- (1 - .at_count(other_rt, event_rt) / length(other_rt)) / n_self
  cumulative <- cumsum(weights); total <- sum(weights)
  fast <- function(t) {
    idx <- .at_count(event_rt, t); out <- numeric(length(t)); ok <- idx > 0
    out[ok] <- cumulative[idx[ok]]; pmax(0, pmin(1, out))
  }
  list(fast = fast, slow = function(t) pmax(0, pmin(1, total - fast(t))), total = total)
}

assessment_donkin_discrimination <- function(RT, CR = NULL, times = NULL) {
  if (length(RT) < 3L) stop("RT must be list(AB, A, B).")
  if (is.null(CR) || length(CR) != length(RT)) CR <- lapply(RT, function(x) rep(TRUE, length(x)))
  if (is.null(times)) times <- .sft_finite_times(RT[1:3])
  ab <- .at_clean(RT[[1]], CR[[1]]); aa <- .at_clean(RT[[2]], CR[[2]]); bb <- .at_clean(RT[[3]], CR[[3]])
  obs <- .at_defective(ab$RT, ab$CR)
  ac <- .at_weighted_survivor_integral(aa$RT[aa$CR], bb$RT, length(aa$RT))
  bc <- .at_weighted_survivor_integral(bb$RT[bb$CR], aa$RT, length(bb$RT))
  ai <- .at_weighted_survivor_integral(aa$RT[!aa$CR], bb$RT, length(aa$RT))
  bi <- .at_weighted_survivor_integral(bb$RT[!bb$CR], aa$RT, length(bb$RT))
  pred <- list(CF = ac$fast(times) + bc$fast(times), CS = ac$slow(times) + bc$slow(times),
               IF = ai$fast(times) + bi$fast(times), IS = ai$slow(times) + bi$slow(times))
  ob <- list(CF = obs$CF(times), CS = obs$CS(times), IF = obs$IF(times), IS = obs$IS(times))
  .at_pack(times, list(time = times, A_IF = .at_log_ratio(pred$IF, ob$IF), A_CF = .at_log_ratio(pred$CF, ob$CF),
                       A_IS = .at_log_ratio(pred$IS, ob$IS), A_CS = .at_log_ratio(pred$CS, ob$CS),
                       pred_IF = pmax(0, pmin(1, pred$IF)), obs_IF = pmax(0, pmin(1, ob$IF)),
                       pred_CF = pmax(0, pmin(1, pred$CF)), obs_CF = pmax(0, pmin(1, ob$CF)),
                       pred_IS = pmax(0, pmin(1, pred$IS)), obs_IS = pmax(0, pmin(1, ob$IS)),
                       pred_CS = pmax(0, pmin(1, pred$CS)), obs_CS = pmax(0, pmin(1, ob$CS))))
}
