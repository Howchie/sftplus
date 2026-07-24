# Timed-assessment (trial-abundance) and Donkin discrimination extensions.


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

