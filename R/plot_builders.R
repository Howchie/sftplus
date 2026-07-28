# Tidy data-frame builders and plotting helpers for curve visualisation.


build_ct_df <- function(times, series_specs, probs = c(.05, .95), qtype = 5) {
  .sft_bind_rows(lapply(series_specs, function(spec) {
    values <- spec$fn(times)
    if (!is.null(spec$rt) && length(spec$rt)) {
      q <- stats::quantile(spec$rt, probs = probs, names = FALSE, na.rm = TRUE, type = qtype)
      values[times < q[1] | times > q[2]] <- NA_real_
    }
    data.frame(Time = times, Ct = values, Series = rep(spec$label, length(times)),
               stringsAsFactors = FALSE)
  }))
}


.sft_prob_specs <- function(x) {
  # Accept either a list of series specs or a capacityGroup() result, matching
  # the dual dispatch the build_rmi_*_df() builders use.
  if (is.list(x) && !is.null(x$specs) && is.list(x$specs)) return(x$specs)
  if (!is.list(x)) stop("x must be a list of series specs or a capacityGroup() result.")
  x
}


.sft_spec_labels <- function(spec, n) {
  # Carry any Subject/Condition labels through to the long frame, so the
  # population summary can be grouped by them.
  extra <- intersect(c("Subject", "Condition"), names(spec))
  if (!length(extra)) return(NULL)
  as.data.frame(lapply(spec[extra], function(v) rep(as.character(v)[1L], n)),
                stringsAsFactors = FALSE)
}


.sft_spec_benchmark <- function(spec, benchmark) {
  from_pool <- identical(benchmark, "all") && !is.null(spec$benchmark)
  x <- as.numeric(if (from_pool) spec$benchmark else spec$rt)
  keep <- if (from_pool) rep(TRUE, length(x)) else .to_correct_indicator(spec$cr, length(x))
  x <- x[keep & is.finite(x)]
  if (length(x) < 2L) NULL else x
}


.sft_prob_frame <- function(spec, probs, value_name, values, times) {
  out <- data.frame(prob = probs, Time = times, values,
                    Series = rep(spec$label, length(probs)),
                    scale = "probability", stringsAsFactors = FALSE)
  names(out)[3L] <- value_name
  labs <- .sft_spec_labels(spec, length(probs))
  if (is.null(labs)) out else cbind(out, labs, row.names = NULL)
}


build_ct_prob_df <- function(x, probs = seq(.05, .95, by = .05), condition = NULL,
                             benchmark = c("series", "all"), qtype = 5) {
  probs <- .sft_probs(probs)
  benchmark <- match.arg(benchmark)
  series_specs <- .sft_prob_specs(x)
  out <- .sft_bind_rows(lapply(series_specs, function(spec) {
    rt <- .sft_spec_benchmark(spec, benchmark)
    if (is.null(rt)) {
      empty <- rep(NA_real_, length(probs))
      return(.sft_prob_frame(spec, probs, "Ct", empty, empty))
    }
    # Each participant is read at their own percentiles, so prob indexes a
    # common position in the distribution, not a common clock time.
    t_p <- .sft_rmi_quantile(rt, probs, qtype = qtype)
    .sft_prob_frame(spec, probs, "Ct", as.numeric(spec$fn(t_p)), as.numeric(t_p))
  }))
  .sft_filter_condition(out, condition)
}


build_cdf_prob_df <- function(x, probs = seq(.05, .95, by = .05), condition = NULL,
                              benchmark = c("all", "series"),
                              smooth_cdf = c("none", "mono", "spline"),
                              smooth_spar = .65, qtype = 5) {
  probs <- .sft_probs(probs)
  benchmark <- match.arg(benchmark)
  smooth_cdf <- match.arg(smooth_cdf)
  series_specs <- .sft_prob_specs(x)
  out <- .sft_bind_rows(lapply(series_specs, function(spec) {
    empty <- rep(NA_real_, length(probs))
    bench <- .sft_spec_benchmark(spec, benchmark)
    rt <- as.numeric(spec$rt)
    keep <- .to_correct_indicator(spec$cr, length(rt))
    rt <- rt[keep & is.finite(rt)]
    if (is.null(bench) || !length(rt)) {
      return(.sft_prob_frame(spec, probs, "CDF", empty, empty))
    }
    # The benchmark times come from the participant's pooled trials, so every
    # condition of that participant is read on one common time axis and the
    # ordering between conditions is preserved.
    t_p <- .sft_rmi_quantile(bench, probs, qtype = qtype)
    y <- .sft_ecdf(rt)(t_p)
    y <- smooth_one_cdf(t_p, y, smooth_cdf = smooth_cdf, smooth_spar = smooth_spar)
    .sft_prob_frame(spec, probs, "CDF", as.numeric(y), as.numeric(t_p))
  }))
  .sft_filter_condition(out, condition)
}


.sft_filter_condition <- function(df, condition) {
  if (is.null(condition) || !"Condition" %in% names(df)) return(df)
  out <- df[df$Condition %in% condition, , drop = FALSE]
  rownames(out) <- NULL
  out
}


build_prob_summary_df <- function(df, value = c("Ct", "CDF"),
                                  center = c("median", "mean"),
                                  ci_probs = c(.1, .9), by = "Condition") {
  center <- match.arg(center)
  value <- if (length(value) > 1L) intersect(value, names(df))[1L] else value
  if (!is.data.frame(df) || is.na(value) || !all(c("prob", value) %in% names(df))) {
    stop("df must be a build_ct_prob_df() or build_cdf_prob_df() result.",
         call. = FALSE)
  }
  ci_probs <- sort(as.numeric(ci_probs))
  keys <- c("prob", intersect(by, names(df)))
  parts <- split(seq_len(nrow(df)),
                 do.call(interaction, c(lapply(keys, function(k) df[[k]]),
                                        list(drop = TRUE, lex.order = TRUE))))
  # Column names follow build_rmi_violation_df(), so the same plotting code
  # serves both probability-scale summaries.
  out <- .sft_bind_rows(lapply(parts, function(idx) {
    v <- as.numeric(df[[value]][idx]); v <- v[is.finite(v)]
    tm <- as.numeric(df$Time[idx]); tm <- tm[is.finite(tm)]
    cbind(df[idx[1L], keys, drop = FALSE],
          data.frame(n = length(v),
                     Time = if (length(tm)) mean(tm) else NA_real_,
                     mean = if (!length(v)) NA_real_ else
                       if (center == "mean") mean(v) else stats::median(v),
                     lower = if (length(v)) stats::quantile(v, ci_probs[1], names = FALSE) else NA_real_,
                     upper = if (length(v)) stats::quantile(v, ci_probs[2], names = FALSE) else NA_real_,
                     center = center, scale = "probability",
                     stringsAsFactors = FALSE),
          row.names = NULL)
  }))
  out <- out[order(out$prob), , drop = FALSE]
  rownames(out) <- NULL
  out
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


build_at_df <- function(times, series_specs, probs = c(.05, .95), qtype = 5) {
  .sft_bind_rows(lapply(series_specs, function(spec) {
    a <- spec$fn
    values <- lapply(c("A_CF", "A_IF", "A_CS", "A_IS"), function(nm) a[[nm]](times))
    names(values) <- c("A_CF", "A_IF", "A_CS", "A_IS")
    if (!is.null(spec$rt) && length(spec$rt)) {
      q <- stats::quantile(spec$rt, probs = probs, names = FALSE, na.rm = TRUE, type = qtype)
      values <- lapply(values, function(x) { x[times < q[1] | times > q[2]] <- NA_real_; x })
    }
    data.frame(Time = times, as.data.frame(values), Series = rep(spec$label, length(times)),
               stringsAsFactors = FALSE)
  }))
}


build_cdf_df <- function(times, series_specs, probs = c(.05, .95),
                         smooth_cdf = c("none", "mono", "spline"), smooth_spar = .65,
                         qtype = 5) {
  smooth_cdf <- match.arg(smooth_cdf)
  .sft_bind_rows(lapply(series_specs, function(spec) {
    if (is.null(spec$rt)) stop("Each series spec must include an `rt` vector.")
    keep <- .to_correct_indicator(spec$cr, length(spec$rt))
    rt <- as.numeric(spec$rt)[keep & is.finite(spec$rt)]
    raw <- .sft_ecdf(rt)(times)
    sm <- smooth_one_cdf(times, raw, smooth_cdf, smooth_spar)
    if (length(rt)) {
      q <- stats::quantile(rt, probs = probs, names = FALSE, type = qtype)
      raw[times < q[1] | times > q[2]] <- NA_real_
      sm[times < q[1] | times > q[2]] <- NA_real_
    } else raw[] <- sm[] <- NA_real_
    data.frame(Time = times, CDF = if (smooth_cdf == "none") raw else sm,
               CDF_raw = raw, CDF_smooth = sm, Series = rep(spec$label, length(times)),
               stringsAsFactors = FALSE)
  }))
}


get_time_range <- function(rt_list, probs = c(.05, .95), by = .001, qtype = 5) {
  if (!is.list(rt_list)) rt_list <- list(rt_list)
  q <- lapply(rt_list, function(x) stats::quantile(x[is.finite(x)], probs = probs,
                                                   names = FALSE, na.rm = TRUE, type = qtype))
  if (!length(q) || any(vapply(q, length, integer(1)) < 2L)) stop("No usable RT quantiles.")
  lo <- max(vapply(q, `[[`, numeric(1), 1)); hi <- min(vapply(q, `[[`, numeric(1), 2))
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) stop("No usable quantile window for the provided RTs.")
  seq(lo, hi, by = by)
}

