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


capacityGroup <- function(inData, acc.cutoff = .9, ratio = TRUE, OR = NULL,
                          stopping.rule = c("OR", "AND", "STST"), plotCt = TRUE, ...) {
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
