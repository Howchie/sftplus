# Semiparametric hierarchical SFT model: input validation, salience resolution,
# pooled RT grid, spline basis, priors, and Stan-data preparation.

# OR-centred hierarchical Bayesian capacity model
#
# The data preparation and posterior transformation functions in this file are
# deliberately generic in their series labels.  A later SIC implementation can
# replace c("AB", "A", "B") with c("HH", "HL", "LH", "LL") without changing
# the pooled grid, exposure likelihood, or curve-summary machinery.

.sft_bayes_series <- c("AB", "A", "B")

.sft_bayes_salience_series <- c(
  "A_L", "A_H", "B_L", "B_H",
  "AB_LL", "AB_LH", "AB_HL", "AB_HH"
)

.sft_bayes_salience_cells <- c("LL", "LH", "HL", "HH")


# Highest-density intervals use the shared empirical HDI helper .sft_hdi()
# defined in sft_core.R (default mass 0.94).

.sft_bayes_numeric <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  if (length(out) != length(x) || anyNA(out)) {
    stop(name, " must contain finite numeric values.", call. = FALSE)
  }
  out
}


.sft_bayes_correct <- function(x) {
  if (is.logical(x)) {
    if (anyNA(x)) stop("Correct must not contain missing values.", call. = FALSE)
    return(x)
  }
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    key <- tolower(trimws(x))
    if (any(!key %in% c("0", "1", "true", "false"))) {
      stop("Correct must contain only logical or 0/1 values.", call. = FALSE)
    }
    return(key %in% c("1", "true"))
  }
  y <- suppressWarnings(as.numeric(x))
  if (anyNA(y) || any(!y %in% c(0, 1))) {
    stop("Correct must contain only logical or 0/1 values.", call. = FALSE)
  }
  y == 1
}


.sft_bayes_salience_labels <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  key <- toupper(trimws(as.character(x)))
  key[key %in% c("LOW", "L")] <- "L"
  key[key %in% c("HIGH", "H")] <- "H"
  if (anyNA(key) || any(!key %in% c("L", "H"))) {
    stop(name, " must return only low/high labels (L/H).", call. = FALSE)
  }
  key
}


.sft_bayes_resolve_salience <- function(channel1, channel2, salience_split) {
  if (is.null(salience_split) || identical(salience_split, FALSE)) {
    return(list(split = FALSE, series = ifelse(channel1 > 0 & channel2 > 0, "AB",
                                                ifelse(channel1 > 0, "A",
                                                       ifelse(channel2 > 0, "B", NA_character_))),
                labels = NULL, mapping = NULL))
  }

  map_one <- function(x, mapper, channel_name) {
    out <- rep(NA_character_, length(x)); positive <- x > 0
    if (!any(positive)) return(out)
    if (is.function(mapper)) {
      out[positive] <- .sft_bayes_salience_labels(mapper(x[positive]), channel_name)
      return(out)
    }
    vals <- as.numeric(mapper)
    if (length(vals) != 2L || any(!is.finite(vals)) || vals[1L] == vals[2L]) {
      stop(channel_name, " salience mapping must contain two distinct finite values.",
           call. = FALSE)
    }
    vals <- sort(vals)
    out[x == vals[1L]] <- "L"
    out[x == vals[2L]] <- "H"
    if (any(positive & is.na(out))) {
      stop(channel_name, " contains positive salience values outside its low/high mapping.",
           call. = FALSE)
    }
    out
  }

  auto_map <- function(x, channel_name) {
    vals <- sort(unique(x[x > 0]))
    if (length(vals) != 2L) {
      stop("salience_split = TRUE/'auto' requires exactly two positive salience levels in ",
           channel_name, "; supply an explicit mapping or function.", call. = FALSE)
    }
    map_one(x, vals, channel_name)
  }

  if (isTRUE(salience_split) ||
      (is.character(salience_split) && length(salience_split) == 1L &&
       tolower(salience_split) == "auto")) {
    l1 <- auto_map(channel1, "Channel1")
    l2 <- auto_map(channel2, "Channel2")
    mapping <- list(Channel1 = sort(unique(channel1[channel1 > 0])),
                    Channel2 = sort(unique(channel2[channel2 > 0])), method = "auto")
  } else if (is.function(salience_split)) {
    l1 <- map_one(channel1, salience_split, "salience_split(Channel1)")
    l2 <- map_one(channel2, salience_split, "salience_split(Channel2)")
    mapping <- list(method = "function")
  } else if (is.list(salience_split)) {
    if (!all(c("Channel1", "Channel2") %in% names(salience_split))) {
      stop("A salience_split list must contain Channel1 and Channel2 mappings.", call. = FALSE)
    }
    l1 <- map_one(channel1, salience_split$Channel1, "Channel1")
    l2 <- map_one(channel2, salience_split$Channel2, "Channel2")
    mapping <- list(Channel1 = salience_split$Channel1, Channel2 = salience_split$Channel2,
                    method = "explicit")
  } else if (is.numeric(salience_split) && length(salience_split) == 2L) {
    l1 <- map_one(channel1, salience_split, "Channel1")
    l2 <- map_one(channel2, salience_split, "Channel2")
    mapping <- list(Channel1 = sort(salience_split), Channel2 = sort(salience_split),
                    method = "explicit")
  } else {
    stop("salience_split must be NULL, TRUE/'auto', a two-value numeric mapping, ",
         "a Channel1/Channel2 mapping list, or a function.", call. = FALSE)
  }

  if (any(channel1 > 0 & is.na(l1)) || any(channel2 > 0 & is.na(l2))) {
    stop("Every positive channel value must receive a low/high salience label.", call. = FALSE)
  }
  series <- ifelse(channel1 > 0 & channel2 > 0, paste0("AB_", l1, l2),
                   ifelse(channel1 > 0, paste0("A_", l1), paste0("B_", l2)))
  list(split = TRUE, series = series, labels = list(Channel1 = l1, Channel2 = l2),
       mapping = mapping)
}


.sft_bayes_validate_input <- function(data, Condition = NULL, cell_mapping = NULL,
                                      salience_split = NULL) {
  if (!is.data.frame(data)) stop("data must be a data.frame.", call. = FALSE)
  required <- c("Subject", "Condition", "RT", "Correct", "Channel1", "Channel2")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("data is missing canonical column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!nrow(data)) stop("data must contain at least one trial.", call. = FALSE)

  subject <- as.character(data$Subject)
  condition <- as.character(data$Condition)
  if (anyNA(subject) || any(!nzchar(subject))) {
    stop("Subject must not contain missing or empty values.", call. = FALSE)
  }
  if (anyNA(condition)) {
    stop("Condition must not contain missing values.", call. = FALSE)
  }
  rt <- .sft_bayes_numeric(data$RT, "RT")
  if (any(!is.finite(rt))) {
    stop("RT must contain only finite values.", call. = FALSE)
  }
  if (any(rt < 0)) stop("RT must be non-negative.", call. = FALSE)
  correct <- .sft_bayes_correct(data$Correct)
  channel1 <- .sft_bayes_numeric(data$Channel1, "Channel1")
  channel2 <- .sft_bayes_numeric(data$Channel2, "Channel2")
  if (any(!is.finite(channel1)) || any(!is.finite(channel2)) ||
      any(channel1 < 0) || any(channel2 < 0)) {
    stop("Channel1 and Channel2 must contain finite, non-negative salience values.",
         call. = FALSE)
  }

  if (!is.null(Condition)) {
    selected <- as.character(Condition)
    if (!length(selected) || anyNA(selected)) {
      stop("Condition selection must contain at least one non-missing value.",
           call. = FALSE)
    }
    keep <- condition %in% selected
    if (!any(keep)) {
      stop("None of the requested Condition values are present in data.", call. = FALSE)
    }
  } else {
    selected <- sort(unique(condition))
    keep <- rep(TRUE, nrow(data))
  }

  data <- data[keep, , drop = FALSE]
  data$.sft_row <- which(keep)
  data$.sft_subject <- subject[keep]
  data$.sft_condition <- condition[keep]
  data$.sft_rt_input <- rt[keep]
  data$.sft_correct <- correct[keep]
  data$.sft_channel1 <- channel1[keep]
  data$.sft_channel2 <- channel2[keep]
  default_series <- ifelse(data$.sft_channel1 > 0 & data$.sft_channel2 > 0, "AB",
                           ifelse(data$.sft_channel1 > 0, "A",
                                  ifelse(data$.sft_channel2 > 0, "B", NA_character_)))
  if (is.null(cell_mapping)) {
    salience <- .sft_bayes_resolve_salience(data$.sft_channel1, data$.sft_channel2,
                                             salience_split)
    data$.sft_series <- salience$series
  } else if (is.function(cell_mapping)) {
    data$.sft_series <- as.character(cell_mapping(data$.sft_channel1, data$.sft_channel2))
    if (length(data$.sft_series) != nrow(data)) {
      stop("cell_mapping must return one series label per trial.", call. = FALSE)
    }
  } else {
    stop("cell_mapping must be NULL or a function of Channel1 and Channel2.", call. = FALSE)
  }
  if (anyNA(data$.sft_series)) {
    stop("Rows with Channel1 == 0 and Channel2 == 0 have no OR capacity series.",
         call. = FALSE)
  }
  split_labels <- all(data$.sft_series %in% .sft_bayes_salience_series)
  pooled_labels <- all(data$.sft_series %in% .sft_bayes_series)
  if (!split_labels && !pooled_labels) {
    stop("cell_mapping must return either pooled labels AB/A/B or salience labels ",
         paste(.sft_bayes_salience_series, collapse = ", "), ".", call. = FALSE)
  }
  if (split_labels && !is.null(salience_split) && is.null(cell_mapping)) {
    # This branch is documentary: the resolver has already validated the
    # labels.  Keeping the flag here makes the prepared object explicit.
    NULL
  }
  if (split_labels && pooled_labels) {
    stop("cell_mapping cannot mix pooled and salience-split series labels.", call. = FALSE)
  }
  data <- data[order(data$.sft_subject, data$.sft_row), , drop = FALSE]
  rownames(data) <- NULL
  list(data = data, selected = selected, original_rows = sum(keep),
       subjects = unique(data$.sft_subject), split = split_labels,
       series = if (split_labels) .sft_bayes_salience_series else .sft_bayes_series,
       salience = if (is.null(cell_mapping)) salience else
         list(split = split_labels, labels = NULL, mapping = list(method = "cell_mapping")))
}


.sft_bayes_grid <- function(rt, n_bins, report_quantiles) {
  if (length(n_bins) != 1L || !is.finite(n_bins) || n_bins < 1 ||
      n_bins != as.integer(n_bins)) {
    stop("n_bins must be a positive integer.", call. = FALSE)
  }
  if (length(report_quantiles) != 2L || any(!is.finite(report_quantiles)) ||
      report_quantiles[1L] <= 0 || report_quantiles[2L] >= 1 ||
      report_quantiles[1L] >= report_quantiles[2L]) {
    stop("report_quantiles must be two ordered probabilities strictly inside (0, 1).",
         call. = FALSE)
  }
  rt <- as.numeric(rt)
  if (!length(rt) || any(!is.finite(rt)) || max(rt) <= 0) {
    stop("At least one strictly positive finite RT is required to construct a grid.",
         call. = FALSE)
  }

  probs <- seq(0, 1, length.out = as.integer(n_bins) + 1L)
  cuts <- as.numeric(stats::quantile(rt, probs = probs, names = FALSE, type = 8))
  cuts[1L] <- 0
  cuts[length(cuts)] <- max(rt)
  boundaries <- sort(unique(c(0, cuts, max(rt))))
  if (length(boundaries) < 2L) stop("The pooled RT grid has no positive-width interval.",
                                    call. = FALSE)
  lower <- boundaries[-length(boundaries)]
  upper <- boundaries[-1L]
  width <- upper - lower
  midpoint <- (lower + upper) / 2
  q <- as.numeric(stats::quantile(rt, probs = report_quantiles, names = FALSE, type = 8))
  report <- midpoint >= q[1L] & midpoint <= q[2L]
  if (!any(report)) {
    # A very small sample or a tied RT distribution can leave no midpoint in
    # the requested window.  Retain the closest available grid point so a
    # valid summary is still returned, while preserving the requested window
    # in metadata.
    distance <- ifelse(midpoint < q[1L], q[1L] - midpoint,
                       ifelse(midpoint > q[2L], midpoint - q[2L], 0))
    report[which.min(distance)] <- TRUE
  }
  grid <- data.frame(bin = seq_along(width), lower = lower, upper = upper,
                     midpoint = midpoint, width = width, report = report,
                     stringsAsFactors = FALSE)
  list(boundaries = boundaries, bins = grid, requested_bins = as.integer(n_bins),
       actual_bins = nrow(grid), report_quantiles = report_quantiles,
       report_range = q)
}


.sft_bayes_basis <- function(grid, basis_dim) {
  J <- nrow(grid$bins)
  if (length(basis_dim) != 1L || !is.finite(basis_dim) || basis_dim < 1 ||
      basis_dim != as.integer(basis_dim)) {
    stop("smoothness$basis_dim must be a positive integer.", call. = FALSE)
  }
  K <- min(as.integer(basis_dim), J)
  x <- if (J == 1L) 0.5 else (grid$bins$midpoint - min(grid$bins$midpoint)) /
    (max(grid$bins$midpoint) - min(grid$bins$midpoint))
  if (K == 1L) {
    B <- matrix(1, nrow = J, ncol = 1L)
  } else {
    degree <- min(3L, K - 1L)
    B <- splines::bs(x, df = K, degree = degree, intercept = TRUE)
    B <- as.matrix(B)
  }
  colnames(B) <- paste0("basis", seq_len(ncol(B)))
  rownames(B) <- paste0("bin", seq_len(nrow(B)))
  list(B = B, dimension = ncol(B), x = x)
}


.sft_bayes_units <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    if (length(x) > 1L) x <- x[[1L]] else stop("rt_units must be supplied.", call. = FALSE)
  }
  key <- tolower(as.character(x))
  if (key %in% c("s", "sec", "secs", "second", "seconds")) return("seconds")
  if (key %in% c("ms", "msec", "msecs", "millisecond", "milliseconds")) return("milliseconds")
  stop("rt_units must be 'seconds' or 'milliseconds'.", call. = FALSE)
}


.sft_bayes_priors <- function(priors = NULL, smoothness = NULL) {
  if (is.null(priors)) priors <- list()
  if (!is.list(priors)) stop("priors must be a list.", call. = FALSE)
  if (is.null(smoothness)) smoothness <- list()
  if (!is.list(smoothness)) stop("smoothness must be a list.", call. = FALSE)
  defaults <- list(
    log_hazard_mean = 0,       # log(1 per second): a neutral seconds-scale centre
    log_hazard_sd = 2,         # allows hazards from much slower to much faster than 1/s
    delta_mean = 0,             # exact hazard-level UCIP centre
    delta_sd = 0.5,              # mildly informative capacity-deviation prior
    speed_sd = 0.5,              # subject speed variation on the log-hazard scale
    asymmetry_sd = 0.5,          # A/B subject asymmetry on the log-hazard scale
    capacity_sd = 0.5,            # subject capacity shift on the log-hazard scale
    smooth_A_sd = 0.5,
    smooth_B_sd = 0.5,
    smooth_delta_sd = 0.25,
    gamma_sd = 0.15,             # shrink salience departures strongly toward zero
    smooth_gamma_sd = 0.15,
    basis_dim = 8L
  )
  # Friendly aliases keep the public controls readable without making the Stan
  # data structure part of the API.
  aliases <- list(
    baseline_log_hazard = "log_hazard_mean",
    hazard_log_sd = "log_hazard_sd",
    sigma_speed = "speed_sd",
    sigma_asymmetry = "asymmetry_sd",
    sigma_capacity = "capacity_sd",
    sigma_A = "smooth_A_sd",
    sigma_B = "smooth_B_sd",
    sigma_delta = "smooth_delta_sd",
    sigma_gamma = "gamma_sd",
    smooth_salience = "smooth_gamma_sd"
  )
  for (nm in names(aliases)) if (!is.null(priors[[nm]]) && is.null(priors[[aliases[[nm]]]])) {
    priors[[aliases[[nm]]]] <- priors[[nm]]
  }
  if (!is.null(priors$baseline_hazard) && is.null(priors$log_hazard_mean)) {
    if (length(priors$baseline_hazard) != 1L || !is.finite(priors$baseline_hazard) ||
        priors$baseline_hazard <= 0) stop("priors$baseline_hazard must be positive.", call. = FALSE)
    priors$log_hazard_mean <- log(priors$baseline_hazard)
  }
  if (!is.null(priors$subject_sd) && length(priors$subject_sd)) {
    ss <- priors$subject_sd
    if (is.null(names(ss))) names(ss) <- c("speed", "asymmetry", "capacity")[seq_along(ss)]
    if (!is.null(ss["speed"])) priors$speed_sd <- ss[["speed"]]
    if (!is.null(ss["asymmetry"])) priors$asymmetry_sd <- ss[["asymmetry"]]
    if (!is.null(ss["capacity"])) priors$capacity_sd <- ss[["capacity"]]
  }
  if (!is.null(smoothness$basis_dim)) priors$basis_dim <- smoothness$basis_dim
  if (!is.null(smoothness$A)) priors$smooth_A_sd <- smoothness$A
  if (!is.null(smoothness$B)) priors$smooth_B_sd <- smoothness$B
  if (!is.null(smoothness$delta)) priors$smooth_delta_sd <- smoothness$delta
  if (!is.null(smoothness$delta_population)) priors$smooth_delta_sd <- smoothness$delta_population
  if (!is.null(smoothness$gamma)) priors$smooth_gamma_sd <- smoothness$gamma
  if (!is.null(smoothness$salience)) priors$smooth_gamma_sd <- smoothness$salience
  out <- defaults
  for (nm in intersect(names(priors), names(out))) out[[nm]] <- priors[[nm]]
  positive <- c("log_hazard_sd", "delta_sd", "speed_sd", "asymmetry_sd", "capacity_sd",
                "smooth_A_sd", "smooth_B_sd", "smooth_delta_sd", "gamma_sd",
                "smooth_gamma_sd")
  for (nm in c("log_hazard_mean", "delta_mean", positive, "basis_dim")) {
    if (length(out[[nm]]) != 1L || !is.finite(out[[nm]])) {
      stop("Prior/smoothness value `", nm, "` must be one finite value.", call. = FALSE)
    }
  }
  if (any(vapply(out[positive], function(x) x <= 0, logical(1)))) {
    stop("Prior and smoothness standard deviations must be positive.", call. = FALSE)
  }
  if (out$basis_dim < 1 || out$basis_dim != as.integer(out$basis_dim)) {
    stop("smoothness$basis_dim must be a positive integer.", call. = FALSE)
  }
  out$basis_dim <- as.integer(out$basis_dim)
  out
}


.sft_bayes_prepared <- function(data, Condition = NULL, n_bins = 25L,
                                report_quantiles = c(0.05, 0.95),
                                rt_units = "seconds", priors = NULL,
                                smoothness = NULL, require_complete = TRUE,
                                cell_mapping = NULL, salience_split = NULL) {
  valid <- .sft_bayes_validate_input(data, Condition, cell_mapping, salience_split)
  d <- valid$data
  unit <- .sft_bayes_units(rt_units)
  scale <- if (unit == "milliseconds") 1 / 1000 else 1
  rt <- d$.sft_rt_input * scale
  grid <- .sft_bayes_grid(rt, n_bins, report_quantiles)
  prior <- .sft_bayes_priors(priors, smoothness)
  basis <- .sft_bayes_basis(grid, prior$basis_dim)
  subjects <- unique(d$.sft_subject)
  I <- length(subjects); J <- grid$actual_bins
  subject_index <- match(d$.sft_subject, subjects)
  series <- valid$series
  series_index <- match(d$.sft_series, series)
  events <- array(0L, dim = c(I, length(series), J),
                  dimnames = list(subjects, series, paste0("bin", seq_len(J))))
  exposure <- array(0, dim = c(I, length(series), J), dimnames = dimnames(events))
  for (k in seq_len(nrow(d))) {
    i <- subject_index[[k]]; s <- series_index[[k]]; t <- rt[[k]]
    e <- pmax(0, pmin(t, grid$bins$upper) - grid$bins$lower)
    exposure[i, s, ] <- exposure[i, s, ] + e
    # Boundaries define (lower, upper] intervals.  Subtracting one from
    # findInterval keeps an event exactly on an upper boundary in the interval
    # whose exposure reaches that boundary.
    b <- findInterval(t, grid$boundaries) - 1L
    b <- max(1L, min(J, b))
    if (isTRUE(d$.sft_correct[[k]])) events[i, s, b] <- events[i, s, b] + 1L
  }
  complete <- matrix(FALSE, nrow = I, ncol = length(series),
                      dimnames = list(subjects, series))
  for (i in seq_len(I)) for (s in seq_along(series)) {
    complete[i, s] <- any(d$.sft_subject == subjects[[i]] &
                            d$.sft_series == series[[s]])
  }
  missing <- which(!complete, arr.ind = TRUE)
  missing_table <- if (nrow(missing)) {
    data.frame(Subject = rownames(complete)[missing[, 1L]],
               Series = colnames(complete)[missing[, 2L]],
               stringsAsFactors = FALSE)
  } else data.frame(Subject = character(), Series = character(), stringsAsFactors = FALSE)
  if (isTRUE(require_complete) && nrow(missing_table)) {
    stop("Each subject must contain every requested series; missing combinations: ",
         paste(paste(missing_table$Subject, missing_table$Series, sep = "/"), collapse = ", "),
         call. = FALSE)
  }
  zero_event_exposure <- sum(events > 0 & exposure == 0)
  if (zero_event_exposure) warning(zero_event_exposure,
                                   " correct event(s) occur at zero exposure on the pooled grid.",
                                   call. = FALSE)
  d$rt_seconds <- rt
  d$Subject <- d$.sft_subject
  d$Condition <- d$.sft_condition
  d$Correct <- d$.sft_correct
  d$Channel1 <- d$.sft_channel1
  d$Channel2 <- d$.sft_channel2
  d$series <- d$.sft_series
  d$subject_index <- subject_index
  d$series_index <- series_index
  d$.sft_row <- NULL; d$.sft_subject <- NULL; d$.sft_condition <- NULL
  d$.sft_rt_input <- NULL; d$.sft_correct <- NULL; d$.sft_channel1 <- NULL
  d$.sft_channel2 <- NULL; d$.sft_series <- NULL
  pooling <- list(
    method = if (valid$split) "channel_presence_with_salience_split" else
      if (is.null(cell_mapping)) "channel_presence" else "custom_function",
    split = valid$split, series = series,
    salience_levels = if (valid$split) c("L", "H") else NULL,
    mapping = valid$salience$mapping, positive_salience = "> 0",
    rules = if (valid$split) data.frame(
      Channel1 = c("> 0, L", "> 0, H", "== 0", "== 0", "> 0, L", "> 0, L", "> 0, H", "> 0, H"),
      Channel2 = c("== 0", "== 0", "> 0, L", "> 0, H", "> 0, L", "> 0, H", "> 0, L", "> 0, H"),
      Series = .sft_bayes_salience_series, stringsAsFactors = FALSE) else data.frame(
        Channel1 = c("> 0", "> 0", "== 0"),
        Channel2 = c("> 0", "== 0", "> 0"),
        Series = .sft_bayes_series, stringsAsFactors = FALSE),
    zero_zero = if (is.null(cell_mapping)) "rejected: no target-present OR series" else
      "handled by the supplied cell_mapping function",
    future_cell_mapping = "supply cell_mapping = function(Channel1, Channel2) ...",
    capacity_estimand = if (valid$split)
      "matched within-salience capacity, averaged over salience cells draw-by-draw" else
      "capacity of the channel-presence mixture"
  )
  stan_data <- list(
    I = I, J = J, K = ncol(basis$B), B = basis$B,
    log_hazard_mean = prior$log_hazard_mean,
    log_hazard_sd = prior$log_hazard_sd,
    delta_mean = prior$delta_mean,
    delta_sd = prior$delta_sd,
    speed_sd = prior$speed_sd,
    asymmetry_sd = prior$asymmetry_sd,
    capacity_sd = prior$capacity_sd,
    smooth_A_sd = prior$smooth_A_sd,
    smooth_B_sd = prior$smooth_B_sd,
    smooth_delta_sd = prior$smooth_delta_sd
  )
  if (!valid$split) {
    stan_data$d_A <- matrix(unname(events[, match("A", series), ]), nrow = I, ncol = J)
    stan_data$d_B <- matrix(unname(events[, match("B", series), ]), nrow = I, ncol = J)
    stan_data$d_AB <- matrix(unname(events[, match("AB", series), ]), nrow = I, ncol = J)
    stan_data$exposure_A <- matrix(unname(exposure[, match("A", series), ]), nrow = I, ncol = J)
    stan_data$exposure_B <- matrix(unname(exposure[, match("B", series), ]), nrow = I, ncol = J)
    stan_data$exposure_AB <- matrix(unname(exposure[, match("AB", series), ]), nrow = I, ncol = J)
  } else {
    for (nm in series) {
      stan_data[[paste0("d_", nm)]] <- matrix(unname(events[, match(nm, series), ]),
                                               nrow = I, ncol = J)
      stan_data[[paste0("exposure_", nm)]] <- matrix(unname(exposure[, match(nm, series), ]),
                                                      nrow = I, ncol = J)
    }
    stan_data$gamma_sd <- prior$gamma_sd
    stan_data$smooth_gamma_sd <- prior$smooth_gamma_sd
  }
  list(
    trials = d, subjects = subjects, series = series,
    events = events, exposure = exposure, complete = complete,
    missing_series = missing_table, zero_event_exposure = zero_event_exposure,
    grid = grid, basis = basis$B, basis_x = basis$x, priors = prior,
    stan_data = stan_data, n_subjects = I, n_bins = J,
    rt_units = unit, rt_scale = scale, selected_conditions = valid$selected,
    pooling = pooling, require_complete = isTRUE(require_complete),
    source_rows = valid$original_rows, split = valid$split,
    salience_labels = valid$salience$labels,
    stan_file = if (valid$split) "or_capacity_salience.stan" else "or_capacity.stan"
  )
}


.sft_bayes_stan_path <- function(split = FALSE) {
  filename <- if (isTRUE(split)) "or_capacity_salience.stan" else "or_capacity.stan"
  paths <- c(system.file("stan", filename, package = "sftplus"),
             file.path("inst", "stan", filename),
             file.path(getwd(), "inst", "stan", filename))
  paths <- paths[nzchar(paths)]
  paths <- paths[file.exists(paths)]
  if (!length(paths)) stop("Cannot locate the packaged OR capacity Stan model.", call. = FALSE)
  normalizePath(paths[[1L]], mustWork = TRUE)
}


.sft_bayes_stan_cache <- new.env(parent = emptyenv())

.sft_bayes_stan_model <- function(split = FALSE) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("semiparametricSFT.bayes() fitting requires the optional rstan package.", call. = FALSE)
  }
  path <- .sft_bayes_stan_path(split)
  key <- normalizePath(path, mustWork = TRUE)
  if (!exists(key, envir = .sft_bayes_stan_cache, inherits = FALSE)) {
    # Keep compiled artifacts in the in-memory cache rather than writing a
    # platform-specific .rds beside the installed Stan source.
    assign(key, rstan::stan_model(file = key, auto_write = FALSE),
           envir = .sft_bayes_stan_cache)
  }
  get(key, envir = .sft_bayes_stan_cache, inherits = FALSE)
}


.sft_capacity_bayes_seed <- function(seed) {
  if (is.null(seed)) seed <- sum(utf8ToInt("capacity-bayes-or-v1"))
  if (length(seed) != 1L || !is.finite(seed) || seed < 1 || seed != as.integer(seed)) {
    stop("seed must be one positive integer.", call. = FALSE)
  }
  as.integer(seed)
}


.sft_bayes_validate_sampling <- function(iter, warmup, chains, cores, refresh, control) {
  scalar_int <- function(x, name, min = 1L) {
    if (length(x) != 1L || !is.finite(x) || x < min || x != as.integer(x))
      stop(name, " must be an integer >= ", min, ".", call. = FALSE)
    as.integer(x)
  }
  iter <- scalar_int(iter, "iter", 2L)
  if (is.null(warmup)) warmup <- floor(iter / 2L)
  warmup <- scalar_int(warmup, "warmup", 0L)
  if (warmup >= iter) stop("warmup must be less than iter.", call. = FALSE)
  chains <- scalar_int(chains, "chains", 1L)
  cores <- scalar_int(cores, "cores", 1L)
  refresh <- scalar_int(refresh, "refresh", 0L)
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list.", call. = FALSE)
  adapt_delta <- if (is.null(control$adapt_delta)) 0.95 else control$adapt_delta
  max_treedepth <- if (is.null(control$max_treedepth)) 12L else control$max_treedepth
  if (length(adapt_delta) != 1L || !is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1)
    stop("control$adapt_delta must lie in (0, 1).", call. = FALSE)
  max_treedepth <- scalar_int(max_treedepth, "control$max_treedepth", 1L)
  list(iter = iter, warmup = warmup, chains = chains, cores = cores,
       refresh = refresh, control = list(adapt_delta = adapt_delta,
                                         max_treedepth = max_treedepth))
}

