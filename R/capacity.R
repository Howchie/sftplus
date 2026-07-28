# Workload capacity coefficients (OR/AND/STST/ID/Altieri) and the group orchestrator.


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


.sft_capacity_window <- function(rts, crs, trim, qtype = 5) {
  # The window over which every contributing estimate is supported.  C(t) is a
  # ratio of cumulative hazards, so near the edge of a cell's own correct RTs
  # one of the risk sets has seen almost no events and the ratio is dominated
  # by a single observation.  Intersecting the per-cell quantile windows keeps
  # only the times where all of them have something to say.
  if (is.null(trim)) return(c(-Inf, Inf))
  trim <- sort(as.numeric(trim))
  if (length(trim) != 2L || any(!is.finite(trim)) || trim[1] < 0 || trim[2] > 1) {
    stop("trim must be NULL or two probabilities in [0, 1].")
  }
  lo <- -Inf; hi <- Inf
  for (i in seq_along(rts)) {
    x <- as.numeric(rts[[i]])[.sft_correct(crs[[i]], length(rts[[i]]))]
    x <- x[is.finite(x)]
    if (!length(x)) return(c(NA_real_, NA_real_))
    q <- stats::quantile(x, probs = trim, names = FALSE, type = qtype)
    lo <- max(lo, q[1]); hi <- min(hi, q[2])
  }
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(c(NA_real_, NA_real_))
  c(lo, hi)
}


.sft_capacity_cells <- function(d, channels, rule) {
  if (rule == "STST") {
    context <- rowSums(d[channels] > 0) == 1L & rowSums(d[channels] < 0) > 0L
    target <- rowSums(d[channels] >= 0) == length(channels) & rowSums(d[channels] != 0) == 1L
    return(list(rt = list(d$RT[context], d$RT[target]),
                cr = list(d$Correct[context], d$Correct[target])))
  }
  context <- rowSums(d[channels] > 0) == length(channels)
  pick <- function(col, j) {
    other <- setdiff(seq_along(channels), j)
    d[[col]][d[[channels[j]]] > 0 & rowSums(d[channels[other]] == 0) == length(other)]
  }
  list(rt = c(list(d$RT[context]), lapply(seq_along(channels), function(j) pick("RT", j))),
       cr = c(list(d$Correct[context]), lapply(seq_along(channels), function(j) pick("Correct", j))))
}


.sft_capacity_plot_mode <- function(plotCt) {
  if (is.logical(plotCt)) return(if (isTRUE(plotCt)) "individual" else "none")
  match.arg(as.character(plotCt)[1L], c("individual", "group", "both", "none"))
}


.sft_capacity_group_plot <- function(times, mat, cond, ratio, col = NULL, lwd = NULL, ...) {
  # Every participant on one axis, with the pointwise median over those who
  # are defined at that time.  Each curve stops where its own support does, so
  # the median is taken over a shrinking set towards the edges; `n` is drawn
  # alongside so that thinning is visible rather than implied.
  if (is.null(mat) || !length(mat) || !any(is.finite(mat))) return(invisible(NULL))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1L)
  keep <- which(colSums(is.finite(mat)) > 0L)
  if (!length(keep)) return(invisible(NULL))
  tt <- times[keep]; mm <- mat[, keep, drop = FALSE]
  ylim <- range(mm[is.finite(mm)])
  graphics::plot(range(tt), ylim, type = "n", xlab = "Time", ylab = "C(t)",
                 main = paste(cond, "\nAll participants"), ...)
  faint <- if (is.null(col)) grDevices::rgb(0, 0, 0, .25) else col
  for (i in seq_len(nrow(mm))) {
    ok <- is.finite(mm[i, ])
    if (any(ok)) graphics::lines(tt[ok], mm[i, ok], col = faint)
  }
  med <- apply(mm, 2L, function(v) if (any(is.finite(v))) stats::median(v[is.finite(v)]) else NA_real_)
  ok <- is.finite(med)
  if (any(ok)) graphics::lines(tt[ok], med[ok], lwd = if (is.null(lwd)) 2.5 else lwd)
  graphics::abline(if (ratio) 1 else 0, 0, lty = 2)
  n_ok <- colSums(is.finite(mm))
  graphics::mtext(paste0("n = ", min(n_ok), " to ", max(n_ok)), side = 3, line = .2, cex = .8)
  invisible(NULL)
}


capacityGroup <- function(inData, acc.cutoff = .9, ratio = TRUE, OR = NULL,
                          stopping.rule = c("OR", "AND", "STST"),
                          plotCt = c("individual", "group", "both", "none"),
                          trim = c(.05, .95), qtype = 5, ...) {
  inData <- .sft_normalize_columns(inData)
  required <- c("Subject", "Condition", "RT", "Correct")
  if (!is.data.frame(inData) || !all(required %in% names(inData))) {
    stop("inData must contain Subject, Condition, RT, and Correct columns.")
  }
  channels <- grep("^Channel", names(inData), value = TRUE)
  if (length(channels) < 2L) stop("Not enough channels for capacity analysis.")
  rule <- if (!is.null(OR)) if (isTRUE(OR)) "OR" else "AND" else match.arg(stopping.rule)
  subjects <- unique(inData$Subject)
  conditions <- unique(inData$Condition)
  plotCt <- .sft_capacity_plot_mode(plotCt)
  if (!any(is.finite(as.numeric(inData$RT)))) {
    stop("inData$RT contains no finite response times.")
  }
  # The grid spans the union of the windows the curves are actually reported
  # on.  Taking the .001 and .999 quantiles of the pooled response times
  # instead let a handful of anticipations from any one participant stretch the
  # grid far below the range where any capacity coefficient is defined.
  windows <- list()
  for (cond in conditions) {
    for (subj in subjects) {
      d <- inData[inData$Condition == cond & inData$Subject == subj, , drop = FALSE]
      if (!nrow(d)) next
      cells <- .sft_capacity_cells(d, channels, rule)
      w <- .sft_capacity_window(cells$rt, cells$cr, trim, qtype)
      if (all(is.finite(w))) windows[[length(windows) + 1L]] <- w
    }
  }
  if (!length(windows)) {
    finite_rt <- as.numeric(inData$RT); finite_rt <- finite_rt[is.finite(finite_rt)]
    span <- stats::quantile(finite_rt, c(.001, .999), names = FALSE)
  } else {
    span <- range(unlist(windows))
  }
  times <- seq(span[1L], span[2L], length.out = 1000)
  fits <- list(); curves <- list(); rows <- list(); z_by_condition <- list()
  specs <- list(); n <- 0L

  for (cond in conditions) {
    z <- numeric()
    for (subj in subjects) {
      d <- inData[inData$Condition == cond & inData$Subject == subj, , drop = FALSE]
      if (!nrow(d)) next
      cells <- .sft_capacity_cells(d, channels, rule)
      rts <- cells$rt; cr <- cells$cr
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
      curve <- if (is.null(fit)) rep(NA_real_, length(times)) else fit$Ct(times)
      if (!is.null(fit)) {
        win <- .sft_capacity_window(rts, cr, trim, qtype)
        curve[is.na(win[1]) | times < win[1] | times > win[2]] <- NA_real_
      }
      curves[[n]] <- curve
      # A ready-made series spec, so the tidy builders can be driven straight
      # from the group fit without re-deriving each participant's cells.  `rt`
      # is the redundant-target cell, which is what C(t) is a statement about.
      if (!is.null(fit)) {
        specs[[length(specs) + 1L]] <-
          list(label = as.character(subj), Subject = as.character(subj),
               Condition = as.character(cond), fn = fit$Ct,
               rt = rts[[1L]], cr = cr[[1L]],
               benchmark = d$RT[.sft_correct(d$Correct, nrow(d))])
      }
      if (!is.null(fit) && !is.null(fit$Ctest)) z <- c(z, unname(fit$Ctest$statistic))
      label <- if (is.null(fit)) NA_character_ else if (!is.null(fit$Ctest) &&
          is.finite(fit$Ctest$p.value) && fit$Ctest$p.value < .05) {
        if (unname(fit$Ctest$statistic) < 0) "Limited" else "Super"
      } else "Nonsignificant"
      rows[[n]] <- data.frame(Subject = as.character(subj), Condition = as.character(cond),
                              Capacity = label, stringsAsFactors = FALSE)
      if (plotCt %in% c("individual", "both") && !is.null(fit) && any(is.finite(curve))) {
        keep <- is.finite(curve)
        plot(times[keep], curve[keep], type = "l", xlab = "Time", ylab = "C(t)",
             main = paste(cond, "\nParticipant", subj), ...)
        abline(if (ratio) 1 else 0, 0, lty = 2)
      }
    }
    if (plotCt %in% c("group", "both")) {
      cond_rows <- which(vapply(rows, function(r) identical(r$Condition, as.character(cond)) &&
                                  !identical(r$Subject, "Group"), logical(1)))
      mat <- do.call(rbind, curves[cond_rows])
      .sft_capacity_group_plot(times, mat, cond, ratio, ...)
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
       times = times, z = z_by_condition, rule = rule, specs = specs,
       trim = trim, qtype = qtype)
}
