# Tidy data-frame builders and plotting helpers for curve visualisation.


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

