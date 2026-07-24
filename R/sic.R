# Survivor and mean interaction contrasts: frequentist and Bayes-factor SIC/MIC.


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

