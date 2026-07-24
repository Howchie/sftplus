# Houpt-Townsend UCIP frequentist tests and the score/information machinery.


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


.sft_score_method <- function(method) {
  key <- tolower(gsub("[^A-Za-z]", "", as.character(method[[1L]])))
  switch(key,
         score = "score", theta = "score",
         capacity = "capacity", logcapacity = "capacity", eta = "capacity",
         stop("score_method must be one of 'score' or 'capacity'."))
}


# Accumulate the four weighted score sums (capacity numerator/denominator and
# their martingale variances) that every UCIP score quantity is built from.
# The Rcpp kernel and the vectorised R fallback are numerically identical: both
# use the pooled event grid with duplicates and the rightmost.closed tie rule.
.sft_ucip_score_components <- function(rt, CR, signs_positive, reverse) {
  if (.ensure_hazard_rcpp() &&
      exists("ucip_score_rcpp", envir = .sft_hazard_env, mode = "function")) {
    cr <- lapply(seq_along(rt), function(i) .sft_correct(CR[[i]], length(rt[[i]])))
    acc <- .sft_hazard_env$ucip_score_rcpp(rt, cr, signs_positive, reverse)
    return(list(capacity_numerator = acc[["capacity_numerator"]],
                capacity_denominator = acc[["capacity_denominator"]],
                numerator_variance = acc[["numerator_variance"]],
                denominator_variance = acc[["denominator_variance"]]))
  }
  all_rt <- unlist(rt, use.names = FALSE)
  tvec <- sort(all_rt[is.finite(all_rt)])
  count_at_least <- function(x) {
    xs <- sort(x[is.finite(x)])
    length(xs) - findInterval(tvec, xs, left.open = TRUE)
  }
  count_at_most <- function(x) findInterval(tvec, sort(x[is.finite(x)]))
  risk <- if (reverse) {
    t(vapply(rt, count_at_most, numeric(length(tvec))))
  } else {
    t(vapply(rt, count_at_least, numeric(length(tvec))))
  }
  weight <- risk[1, ] * colSums(risk[-1, , drop = FALSE]) / colSums(risk)
  weight[!is.finite(weight)] <- 0
  index <- rep(seq_along(rt), lengths(rt))
  flat_cr <- unlist(CR, use.names = FALSE)
  ord <- order(all_rt)
  index <- index[ord]; cr_sorted <- flat_cr[ord]; rt_sorted <- all_rt[ord]
  event <- as.logical(cr_sorted) & is.finite(rt_sorted)
  event[is.na(event)] <- FALSE
  num <- 0; den <- 0; num_v <- 0; den_v <- 0
  for (i in seq_along(rt)) {
    idx <- event & index == i
    ti <- findInterval(rt_sorted[idx], tvec, rightmost.closed = TRUE)
    ri <- risk[i, ti]; ri[ri <= 0] <- NA_real_
    terms <- weight[ti] / ri
    terms <- terms[is.finite(terms)]
    if (signs_positive[i]) {
      num <- num + sum(terms); num_v <- num_v + sum(terms^2)
    } else {
      den <- den + sum(terms); den_v <- den_v + sum(terms^2)
    }
  }
  list(capacity_numerator = num, capacity_denominator = den,
       numerator_variance = num_v, denominator_variance = den_v)
}


.sft_ucip_score <- function(RT, CR = NULL, OR = NULL,
                            stopping.rule = c("OR", "AND", "STST")) {
  if (!is.null(OR)) {
    rule <- if (isTRUE(OR)) "OR" else "AND"
  } else {
    rule <- match.arg(stopping.rule %||% "OR", c("OR", "AND", "STST"))
  }
  if (!is.list(RT) || length(RT) < 2L) stop("RT must be a list of at least two conditions.")
  CR <- .sft_cr_list(RT, CR)
  rt <- lapply(RT, as.numeric)
  if (!any(is.finite(unlist(rt, use.names = FALSE)))) {
    stop("RT contains no finite response times.")
  }
  ncond <- length(rt)
  reverse <- rule != "OR"
  signs <- if (rule == "OR") c(1, rep(-1, ncond - 1L)) else c(-1, rep(1, ncond - 1L))
  alternative <- if (rule == "OR") {
    "response times are different than those predicted by the UCIP-OR model"
  } else {
    paste0("response times are different than those predicted by the UCIP-",
           rule, " model")
  }

  acc <- .sft_ucip_score_components(rt, CR, signs > 0, reverse)
  capacity_numerator <- acc$capacity_numerator
  capacity_denominator <- acc$capacity_denominator
  numerator_variance <- acc$numerator_variance
  denominator_variance <- acc$denominator_variance

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
       rule = rule,
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

