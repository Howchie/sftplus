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
    stop("capacity.bayes() fitting requires the optional rstan package.", call. = FALSE)
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

.sft_bayes_dim_array <- function(x, draws, I, J, name) {
  d <- dim(x)
  if (is.null(d) || length(d) != 3L || !all(d == c(draws, I, J))) {
    stop("Unexpected posterior dimension for ", name, ".", call. = FALSE)
  }
  x
}

.sft_bayes_transform <- function(log_hazard, widths, times, subjects,
                                 new_subject = NULL) {
  required <- c("A", "B", "AB")
  if (!all(required %in% names(log_hazard))) stop("log_hazard needs A, B, and AB arrays.",
                                                     call. = FALSE)
  draws <- dim(log_hazard$A)[1L]; I <- dim(log_hazard$A)[2L]; J <- dim(log_hazard$A)[3L]
  h <- lapply(log_hazard, exp)
  H <- lapply(h, function(x) {
    out <- x
    for (m in seq_len(draws)) for (i in seq_len(I)) out[m, i, ] <- cumsum(x[m, i, ] * widths)
    out
  })
  D <- H$AB - H$A - H$B
  denom <- H$A + H$B
  C <- H$AB / denom
  C[!is.finite(C)] <- NA_real_
  dimnames(D) <- list(draw = seq_len(draws), Subject = subjects, bin = seq_along(times))
  dimnames(C) <- dimnames(D)
  for (nm in names(H)) dimnames(H[[nm]]) <- dimnames(D)
  out <- list(H_A = H$A, H_B = H$B, H_AB = H$AB, D = D, C = C,
              subject = list(H_A = H$A, H_B = H$B, H_AB = H$AB, D = D, C = C),
              times = times, widths = widths, subjects = subjects,
              identity_error = H$AB - H$A - H$B)
  if (!is.null(new_subject)) {
    new_h <- lapply(new_subject, exp)
    new_H <- lapply(new_h, function(x) {
      out_x <- x
      for (m in seq_len(draws)) out_x[m, ] <- cumsum(x[m, ] * widths)
      out_x
    })
    new_D <- new_H$AB - new_H$A - new_H$B
    new_C <- new_H$AB / (new_H$A + new_H$B); new_C[!is.finite(new_C)] <- NA_real_
    out$new_subject <- list(H_A = new_H$A, H_B = new_H$B, H_AB = new_H$AB,
                            D = new_D, C = new_C)
  }
  out
}

.sft_bayes_cumulative_series <- function(log_hazard, widths, subjects, times) {
  draws <- dim(log_hazard[[1L]])[1L]
  I <- dim(log_hazard[[1L]])[2L]
  J <- dim(log_hazard[[1L]])[3L]
  H <- lapply(log_hazard, function(x) {
    if (is.null(dim(x)) || !all(dim(x) == c(draws, I, J)))
      stop("All salience posterior arrays must have draw by subject by bin dimensions.", call. = FALSE)
    ans <- exp(x)
    for (m in seq_len(draws)) for (i in seq_len(I))
      ans[m, i, ] <- cumsum(ans[m, i, ] * widths)
    ans
  })
  dn <- list(draw = seq_len(draws), Subject = subjects, bin = seq_along(times))
  for (nm in names(H)) dimnames(H[[nm]]) <- dn
  H
}

.sft_bayes_transform_salience <- function(log_hazard, widths, times, subjects,
                                          new_subject = NULL) {
  required <- .sft_bayes_salience_series
  if (!all(required %in% names(log_hazard))) {
    stop("Salience-split posterior draws must contain A_L, A_H, B_L, B_H, and AB_LL/AB_LH/AB_HL/AB_HH.",
         call. = FALSE)
  }
  H <- .sft_bayes_cumulative_series(log_hazard, widths, subjects, times)
  draws <- dim(H[[1L]])[1L]; I <- dim(H[[1L]])[2L]; J <- dim(H[[1L]])[3L]
  cells <- .sft_bayes_salience_cells
  Dcells <- array(NA_real_, dim = c(draws, I, J, 4L),
                  dimnames = list(draw = seq_len(draws), Subject = subjects,
                                  bin = seq_along(times), cell = cells))
  Ccells <- Dcells
  for (cc in seq_along(cells)) {
    cell <- cells[[cc]]
    u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
    ha <- H[[paste0("A_", u)]]; hb <- H[[paste0("B_", v)]]
    hab <- H[[paste0("AB_", cell)]]
    Dcells[, , , cc] <- hab - ha - hb
    denom <- ha + hb
    z <- hab / denom; z[!is.finite(z)] <- NA_real_
    Ccells[, , , cc] <- z
  }
  Davg <- apply(Dcells, c(1L, 2L, 3L), mean, na.rm = TRUE)
  logC <- log(Ccells)
  logCavg <- apply(logC, c(1L, 2L, 3L), mean, na.rm = TRUE)
  Cavg <- exp(logCavg)
  sic <- exp(-H$AB_LL) - exp(-H$AB_LH) - exp(-H$AB_HL) + exp(-H$AB_HH)
  dn <- list(draw = seq_len(draws), Subject = subjects, bin = seq_along(times))
  dimnames(Davg) <- dn; dimnames(Cavg) <- dn; dimnames(logCavg) <- dn
  dimnames(sic) <- dn
  cap <- list(D = Dcells, C = Ccells, D_average = Davg,
              logC_average = logCavg, C_average = Cavg,
              cells = cells, weights = rep(1 / 4, 4L))
  out <- list(H = H, S = lapply(H, function(x) exp(-x)),
              H_A_L = H$A_L, H_A_H = H$A_H, H_B_L = H$B_L, H_B_H = H$B_H,
              H_AB_LL = H$AB_LL, H_AB_LH = H$AB_LH,
              H_AB_HL = H$AB_HL, H_AB_HH = H$AB_HH,
              D_cells = Dcells, C_cells = Ccells, D = Davg, C = Cavg,
              logC_average = logCavg, SIC = sic, capacity = cap,
              subject = list(D = Davg, C = Cavg, logC = logCavg,
                             D_cells = Dcells, C_cells = Ccells, SIC = sic),
              times = times, widths = widths, subjects = subjects,
              identity_error = Dcells)
  if (!is.null(new_subject)) {
    new_H <- lapply(new_subject, function(x) {
      if (is.null(dim(x)) || length(dim(x)) != 2L || dim(x)[1L] != draws ||
          dim(x)[2L] != J) stop("Unexpected new-subject salience posterior dimensions.", call. = FALSE)
      ans <- exp(x)
      for (m in seq_len(draws)) ans[m, ] <- cumsum(ans[m, ] * widths)
      ans
    })
    new_D <- array(NA_real_, dim = c(draws, J, 4L),
                   dimnames = list(draw = seq_len(draws), bin = seq_along(times), cell = cells))
    new_C <- new_D
    for (cc in seq_along(cells)) {
      cell <- cells[[cc]]; u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
      ha <- new_H[[paste0("A_", u)]]; hb <- new_H[[paste0("B_", v)]]
      hab <- new_H[[paste0("AB_", cell)]]
      new_D[, , cc] <- hab - ha - hb
      z <- hab / (ha + hb); z[!is.finite(z)] <- NA_real_; new_C[, , cc] <- z
    }
    new_Davg <- apply(new_D, c(1L, 2L), mean, na.rm = TRUE)
    new_logCavg <- apply(log(new_C), c(1L, 2L), mean, na.rm = TRUE)
    new_Cavg <- exp(new_logCavg)
    new_sic <- exp(-new_H$AB_LL) - exp(-new_H$AB_LH) -
      exp(-new_H$AB_HL) + exp(-new_H$AB_HH)
    out$new_subject <- list(H = new_H, D = new_Davg, C = new_Cavg,
                            logC = new_logCavg, D_cells = new_D, C_cells = new_C,
                            SIC = new_sic)
  }
  out
}

.sft_bayes_summary_matrix <- function(x, times, level, measure, subject = NA_character_,
                                      report, hdi = 0.94, rope = 0, series_label = "OR") {
  # x is draws x time for one curve.
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L, ncol = length(times))
  if (ncol(x) != length(times)) x <- matrix(as.numeric(x), nrow = dim(x)[1L], ncol = length(times))
  ans <- lapply(seq_along(times), function(j) {
    v <- x[, j]
    interval <- .sft_hdi(v, hdi)
    threshold <- if (measure %in% c("difference", "D", "SIC")) 0 else 1
    data.frame(Level = level, Subject = subject, Series = series_label, Measure = measure,
               Time = times[[j]], Ct = mean(v, na.rm = TRUE),
               Mean = mean(v, na.rm = TRUE), Lower = interval[[1L]],
               Upper = interval[[2L]],
               Prob_super = mean(v > threshold, na.rm = TRUE),
               Prob_limited = mean(v < threshold, na.rm = TRUE),
               Prob_rope = if (rope > 0) mean(abs(v - threshold) <= rope, na.rm = TRUE) else NA_real_,
               Report = report[[j]], stringsAsFactors = FALSE)
  })
  do.call(rbind, ans)
}

.sft_bayes_split_summaries <- function(transformed, grid, subjects, hdi = 0.94,
                                       rope = 0) {
  times <- grid$bins$upper; report <- grid$bins$report
  parts <- list()
  add_matrix <- function(x, level, measure, subject = NA_character_, series = "OR-average") {
    parts[[length(parts) + 1L]] <<- .sft_bayes_summary_matrix(
      x, times, level, measure, subject, report, hdi, rope, series)
  }
  for (measure in c("D", "C")) {
    for (i in seq_along(subjects))
      add_matrix(transformed$subject[[measure]][, i, , drop = FALSE][, 1L, ],
                 "subject", measure, subjects[[i]])
    add_matrix(transformed$population[[measure]], "population", measure)
    if (!is.null(transformed$new_subject)) add_matrix(transformed$new_subject[[measure]],
                                                     "new_subject", measure = measure,
                                                     series = "OR-average")
  }
  for (cc in seq_along(.sft_bayes_salience_cells)) {
    cell <- .sft_bayes_salience_cells[[cc]]
    for (measure in c("D", "C")) {
      for (i in seq_along(subjects))
        add_matrix(transformed$subject[[paste0(measure, "_cells")]][, i, , cc, drop = FALSE][, 1L, , 1L],
                   "subject", measure, subjects[[i]], cell)
      add_matrix(transformed$population[[paste0(measure, "_cells")]][, , cc],
                 "population", measure, series = cell)
      if (!is.null(transformed$new_subject))
        add_matrix(transformed$new_subject[[paste0(measure, "_cells")]][, , cc],
                   "new_subject", measure = measure, series = cell)
    }
  }
  for (i in seq_along(subjects))
    add_matrix(transformed$subject$SIC[, i, , drop = FALSE][, 1L, ],
               "subject", "SIC", subjects[[i]], "SIC")
  add_matrix(transformed$population$SIC, "population", "SIC", series = "SIC")
  if (!is.null(transformed$new_subject)) add_matrix(transformed$new_subject$SIC,
                                                   "new_subject", measure = "SIC",
                                                   series = "SIC")
  tidy <- do.call(rbind, parts); rownames(tidy) <- NULL
  list(tidy = tidy,
       subject = tidy[tidy$Level == "subject", , drop = FALSE],
       population = tidy[tidy$Level == "population", , drop = FALSE],
       new_subject = tidy[tidy$Level == "new_subject", , drop = FALSE])
}

.sft_bayes_summaries <- function(transformed, grid, subjects, hdi = 0.94, rope = 0) {
  times <- grid$bins$upper; report <- grid$bins$report
  parts <- list()
  for (measure in c("D", "C")) {
    for (i in seq_along(subjects)) {
      parts[[length(parts) + 1L]] <- .sft_bayes_summary_matrix(
        transformed$subject[[measure]][, i, , drop = FALSE][, 1L, ],
        times, "subject", measure, subjects[[i]], report, hdi, rope)
    }
    population <- apply(transformed$subject[[measure]], c(1L, 3L), mean)
    parts[[length(parts) + 1L]] <- .sft_bayes_summary_matrix(
      population, times, "population", measure, NA_character_, report, hdi, rope)
    if (!is.null(transformed$new_subject)) {
      parts[[length(parts) + 1L]] <- .sft_bayes_summary_matrix(
        transformed$new_subject[[measure]], times, "new_subject", measure,
        NA_character_, report, hdi, rope)
    }
  }
  tidy <- do.call(rbind, parts)
  rownames(tidy) <- NULL
  list(
    tidy = tidy,
    subject = tidy[tidy$Level == "subject", , drop = FALSE],
    population = tidy[tidy$Level == "population", , drop = FALSE],
    new_subject = tidy[tidy$Level == "new_subject", , drop = FALSE]
  )
}

.sft_bayes_effects_from_posterior <- function(raw, priors, I, J, subjects, widths,
                                              seed, new_subject = TRUE) {
  needed <- c("eta_A", "eta_B", "eta_AB")
  if (!all(needed %in% names(raw))) stop("The Stan fit did not return all log-hazard arrays.",
                                          call. = FALSE)
  draws <- dim(raw$eta_A)[1L]
  log_hazard <- list(A = .sft_bayes_dim_array(raw$eta_A, draws, I, J, "eta_A"),
                     B = .sft_bayes_dim_array(raw$eta_B, draws, I, J, "eta_B"),
                     AB = .sft_bayes_dim_array(raw$eta_AB, draws, I, J, "eta_AB"))
  new_eta <- NULL
  if (isTRUE(new_subject)) {
    if (!all(c("population_A", "population_B", "population_delta",
               "sigma_speed", "sigma_asymmetry", "sigma_capacity") %in% names(raw))) {
      warning("The fit lacks population effect draws; new-subject prediction was omitted.", call. = FALSE)
    } else {
      set.seed(seed)
      # Draw a fresh subject with the same scaling as the Stan model, where each
      # effect is <fixed prior sd> * sigma_* * z with z ~ N(0, 1). Omitting the
      # fixed prior sd would over-disperse the predicted subject.
      sigma_speed <- as.numeric(raw$sigma_speed); sigma_asymmetry <- as.numeric(raw$sigma_asymmetry)
      sigma_capacity <- as.numeric(raw$sigma_capacity)
      speed <- stats::rnorm(draws, 0, priors$speed_sd * sigma_speed)
      asym <- stats::rnorm(draws, 0, priors$asymmetry_sd * sigma_asymmetry)
      cap <- stats::rnorm(draws, 0, priors$capacity_sd * sigma_capacity)
      pa <- raw$population_A; pb <- raw$population_B; pd <- raw$population_delta
      new_a <- pa + speed + asym / 2
      new_b <- pb + speed - asym / 2
      # AB centres on the UCIP log-sum-exp plus the population capacity deviation
      # (pd) and the subject capacity shift, mirroring eta_AB in the Stan model.
      new_eta <- list(A = new_a, B = new_b,
                      AB = pmax(new_a, new_b) + log1p(exp(-abs(new_a - new_b))) + pd + cap)
    }
  }
  list(log_hazard = log_hazard, new_log_hazard = new_eta)
}

.sft_bayes_effects_from_salience_posterior <- function(raw, priors, I, J,
                                                       subjects, seed,
                                                       new_subject = TRUE) {
  eta_names <- paste0("eta_", .sft_bayes_salience_series)
  if (!all(eta_names %in% names(raw))) {
    stop("The salience-split fit did not return all eight log-hazard arrays.", call. = FALSE)
  }
  first_name <- eta_names[[1L]]
  first_eta <- raw[[first_name]]
  draws <- dim(first_eta)[1L]
  log_hazard <- setNames(lapply(seq_along(.sft_bayes_salience_series), function(k) {
    nm <- eta_names[[k]]
    .sft_bayes_dim_array(raw[[nm]], draws, I, J, nm)
  }), .sft_bayes_salience_series)
  new_eta <- NULL
  if (isTRUE(new_subject)) {
    needed <- c("population_A_L", "population_A_H", "population_B_L", "population_B_H",
                "population_delta", "sigma_speed", "sigma_asymmetry", "sigma_capacity")
    if (!all(needed %in% names(raw))) {
      warning("The fit lacks population effect draws; new-subject prediction was omitted.",
              call. = FALSE)
    } else {
      set.seed(seed)
      sigma_speed <- as.numeric(raw$sigma_speed)
      sigma_asymmetry <- as.numeric(raw$sigma_asymmetry)
      sigma_capacity <- as.numeric(raw$sigma_capacity)
      speed <- stats::rnorm(draws, 0, priors$speed_sd * sigma_speed)
      asym <- stats::rnorm(draws, 0, priors$asymmetry_sd * sigma_asymmetry)
      cap <- stats::rnorm(draws, 0, priors$capacity_sd * sigma_capacity)
      new_eta <- list(
        A_L = raw$population_A_L + speed + asym / 2,
        A_H = raw$population_A_H + speed + asym / 2,
        B_L = raw$population_B_L + speed - asym / 2,
        B_H = raw$population_B_H + speed - asym / 2
      )
      gamma <- matrix(0, nrow = draws, ncol = J)
      if ("population_gamma" %in% names(raw)) {
        if (length(dim(raw$population_gamma)) != 3L ||
            !all(dim(raw$population_gamma) == c(draws, J, 4L))) {
          stop("Unexpected population_gamma dimensions in the Stan fit.", call. = FALSE)
        }
      }
      for (cell in .sft_bayes_salience_cells) {
        if ("population_gamma" %in% names(raw)) {
          cc <- match(cell, .sft_bayes_salience_cells)
          gamma <- matrix(raw$population_gamma[, , cc, drop = TRUE], nrow = draws, ncol = J)
        } else gamma[,] <- 0
        u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
        a <- new_eta[[paste0("A_", u)]]
        b <- new_eta[[paste0("B_", v)]]
        lse <- pmax(a, b) + log1p(exp(-abs(a - b)))
        new_eta[[paste0("AB_", cell)]] <- lse + raw$population_delta + cap + gamma
      }
    }
  }
  list(log_hazard = log_hazard, new_log_hazard = new_eta)
}

.sft_bayes_prior_draws <- function(prepared, n_draws, seed) {
  n_draws <- max(1L, as.integer(n_draws)); J <- prepared$n_bins; K <- ncol(prepared$basis)
  B <- prepared$basis; p <- prepared$priors
  set.seed(seed)
  pa <- pb <- pd <- matrix(NA_real_, nrow = n_draws, ncol = J)
  for (m in seq_len(n_draws)) {
    beta_a <- stats::rnorm(K, p$log_hazard_mean, p$log_hazard_sd)
    beta_b <- stats::rnorm(K, p$log_hazard_mean, p$log_hazard_sd)
    beta_d <- stats::rnorm(K, p$delta_mean, p$delta_sd)
    if (K >= 3L) {
      # Second-difference random walk: beta[k] - 2 beta[k-1] + beta[k-2] ~ N(0, sd).
      # Draw it as a cumulative recursion so the prior predictive matches the
      # smoothing prior used by the Stan model.
      for (k in 3:K) {
        beta_a[k] <- 2 * beta_a[k - 1L] - beta_a[k - 2L] + stats::rnorm(1L, 0, p$smooth_A_sd)
        beta_b[k] <- 2 * beta_b[k - 1L] - beta_b[k - 2L] + stats::rnorm(1L, 0, p$smooth_B_sd)
        beta_d[k] <- 2 * beta_d[k - 1L] - beta_d[k - 2L] + stats::rnorm(1L, 0, p$smooth_delta_sd)
      }
    }
    pa[m, ] <- as.numeric(B %*% beta_a)
    pb[m, ] <- as.numeric(B %*% beta_b)
    pd[m, ] <- as.numeric(B %*% beta_d)
  }
  speed <- stats::rnorm(n_draws, 0, p$speed_sd)
  asym <- stats::rnorm(n_draws, 0, p$asymmetry_sd)
  cap <- stats::rnorm(n_draws, 0, p$capacity_sd)
  a <- pa + speed + asym / 2; b <- pb + speed - asym / 2
  ab <- pmax(a, b) + log1p(exp(-abs(a - b))) + pd + cap
  list(A = a, B = b, AB = ab)
}

.sft_bayes_prior_predictive <- function(prepared, n_draws, seed) {
  raw <- .sft_bayes_prior_draws(prepared, n_draws, seed)
  widths <- prepared$grid$bins$width; times <- prepared$grid$bins$upper
  H <- lapply(raw, function(x) {
    ans <- x
    for (m in seq_len(nrow(x))) ans[m, ] <- cumsum(exp(x[m, ]) * widths)
    ans
  })
  D <- H$AB - H$A - H$B; C <- H$AB / (H$A + H$B)
  C[!is.finite(C)] <- NA_real_
  counts <- lapply(seq_along(.sft_bayes_series), function(s) {
    E <- matrix(prepared$exposure[, s, ], nrow = prepared$n_subjects,
                ncol = prepared$n_bins)
    Hname <- .sft_bayes_series[[s]]
    rates <- raw[[Hname]]
    total_mean <- vapply(seq_len(nrow(rates)), function(m) {
      sum(E * matrix(exp(rates[m, ]) * prepared$grid$bins$width,
                     nrow = nrow(E), ncol = ncol(E), byrow = TRUE))
    }, numeric(1))
    stats::rpois(nrow(rates), lambda = pmax(total_mean, .Machine$double.eps))
  })
  names(counts) <- .sft_bayes_series
  list(log_hazard = raw, H_A = H$A, H_B = H$B, H_AB = H$AB, D = D, C = C,
       times = times, counts = counts)
}

.sft_bayes_prior_draws_salience <- function(prepared, n_draws, seed) {
  n_draws <- max(1L, as.integer(n_draws)); J <- prepared$n_bins
  I <- prepared$n_subjects; K <- ncol(prepared$basis); B <- prepared$basis; p <- prepared$priors
  set.seed(seed)
  out <- lapply(.sft_bayes_salience_series,
                function(nm) array(0, dim = c(n_draws, I, J)))
  names(out) <- .sft_bayes_salience_series
  for (m in seq_len(n_draws)) {
    draw_curve <- function(mean, sd, smooth_sd) {
      beta <- stats::rnorm(K, mean, sd)
      if (K >= 3L) for (k in 3:K)
        beta[k] <- 2 * beta[k - 1L] - beta[k - 2L] + stats::rnorm(1L, 0, smooth_sd)
      as.numeric(B %*% beta)
    }
    pa_l <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_A_sd)
    pa_h <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_A_sd)
    pb_l <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_B_sd)
    pb_h <- draw_curve(p$log_hazard_mean, p$log_hazard_sd, p$smooth_B_sd)
    pd <- draw_curve(p$delta_mean, p$delta_sd, p$smooth_delta_sd)
    gamma_beta <- matrix(stats::rnorm(K * 4L, 0, p$gamma_sd), nrow = K, ncol = 4L)
    if (K >= 3L) for (k in 3:K) for (cc in 1:4)
      gamma_beta[k, cc] <- 2 * gamma_beta[k - 1L, cc] - gamma_beta[k - 2L, cc] +
        stats::rnorm(1L, 0, p$smooth_gamma_sd)
    gamma_beta <- gamma_beta - rowMeans(gamma_beta)
    gamma <- B %*% gamma_beta
    speed <- stats::rnorm(I, 0, p$speed_sd)
    asym <- stats::rnorm(I, 0, p$asymmetry_sd)
    cap <- stats::rnorm(I, 0, p$capacity_sd)
    for (i in seq_len(I)) {
      a_l <- pa_l + speed[[i]] + asym[[i]] / 2
      a_h <- pa_h + speed[[i]] + asym[[i]] / 2
      b_l <- pb_l + speed[[i]] - asym[[i]] / 2
      b_h <- pb_h + speed[[i]] - asym[[i]] / 2
      out$A_L[m, i, ] <- a_l; out$A_H[m, i, ] <- a_h
      out$B_L[m, i, ] <- b_l; out$B_H[m, i, ] <- b_h
      out$AB_LL[m, i, ] <- pmax(a_l, b_l) + log1p(exp(-abs(a_l - b_l))) + pd + cap[[i]] + gamma[, 1L]
      out$AB_LH[m, i, ] <- pmax(a_l, b_h) + log1p(exp(-abs(a_l - b_h))) + pd + cap[[i]] + gamma[, 2L]
      out$AB_HL[m, i, ] <- pmax(a_h, b_l) + log1p(exp(-abs(a_h - b_l))) + pd + cap[[i]] + gamma[, 3L]
      out$AB_HH[m, i, ] <- pmax(a_h, b_h) + log1p(exp(-abs(a_h - b_h))) + pd + cap[[i]] + gamma[, 4L]
    }
  }
  out
}

.sft_bayes_predictive_from_series <- function(prepared, log_hazard, seed,
                                              n_draws = NULL) {
  set.seed(seed)
  all_draws <- dim(log_hazard[[1L]])[1L]
  if (!is.null(n_draws) && length(n_draws) == 1L && is.finite(n_draws) && n_draws >= 1L) {
    n_draws <- min(all_draws, as.integer(n_draws))
    keep <- if (n_draws < all_draws) sample(seq_len(all_draws), n_draws) else seq_len(all_draws)
  } else keep <- seq_len(all_draws)
  predicted <- lapply(prepared$series, function(nm) {
    eta <- log_hazard[[nm]][keep, , , drop = FALSE]
    E <- prepared$exposure[, nm, , drop = FALSE]
    E <- matrix(E, nrow = prepared$n_subjects, ncol = prepared$n_bins)
    out <- array(0L, dim = c(length(keep), prepared$n_subjects, prepared$n_bins))
    for (m in seq_along(keep)) {
      rate <- exp(eta[m, , ])
      rate <- matrix(rate, nrow = prepared$n_subjects, ncol = prepared$n_bins)
      out[m, , ] <- matrix(stats::rpois(length(E), as.numeric(E * rate)),
                            nrow = prepared$n_subjects, ncol = prepared$n_bins)
    }
    out
  })
  names(predicted) <- prepared$series
  summary <- do.call(rbind, lapply(prepared$series, function(nm) {
    vals <- apply(predicted[[nm]], 1L, sum)
    interval <- .sft_hdi(vals)
    s <- match(nm, prepared$series)
    data.frame(Series = nm, Observed = sum(prepared$events[, s, ]),
               Mean = mean(vals), Lower = interval[[1L]], Upper = interval[[2L]],
               stringsAsFactors = FALSE)
  }))
  list(draws = predicted, summary = summary)
}

.sft_bayes_prior_predictive_salience <- function(prepared, n_draws, seed) {
  raw <- .sft_bayes_prior_draws_salience(prepared, n_draws, seed)
  H <- lapply(raw, function(x) {
    ans <- exp(x); for (m in seq_len(dim(ans)[1L])) for (i in seq_len(dim(ans)[2L]))
      ans[m, i, ] <- cumsum(ans[m, i, ] * prepared$grid$bins$width)
    ans
  })
  cells <- .sft_bayes_salience_cells
  Dcells <- array(NA_real_, dim = c(dim(H[[1L]]), 4L))
  Ccells <- Dcells
  for (cc in seq_along(cells)) {
    cell <- cells[[cc]]; u <- substr(cell, 1L, 1L); v <- substr(cell, 2L, 2L)
    Dcells[, , , cc] <- H[[paste0("AB_", cell)]] - H[[paste0("A_", u)]] - H[[paste0("B_", v)]]
    z <- H[[paste0("AB_", cell)]] /
      (H[[paste0("A_", u)]] + H[[paste0("B_", v)]])
    z[!is.finite(z)] <- NA_real_; Ccells[, , , cc] <- z
  }
  sic <- exp(-H$AB_LL) - exp(-H$AB_LH) - exp(-H$AB_HL) + exp(-H$AB_HH)
  list(log_hazard = raw, H = H, D_cells = Dcells, C_cells = Ccells,
       D = apply(Dcells, c(1L, 2L, 3L), mean, na.rm = TRUE),
       C = exp(apply(log(Ccells), c(1L, 2L, 3L), mean, na.rm = TRUE)),
       SIC = sic, times = prepared$grid$bins$upper,
       counts = .sft_bayes_predictive_from_series(prepared, raw, seed + 7L)$summary)
}

.sft_bayes_diagnostics <- function(fit) {
  sm <- rstan::summary(fit)$summary
  pars <- rownames(sm)
  convergence <- data.frame(parameter = pars, mean = sm[, "mean"], sd = sm[, "sd"],
                            rhat = sm[, "Rhat"], n_eff = sm[, "n_eff"], stringsAsFactors = FALSE)
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth <- max(vapply(sampler, function(x) max(x[, "treedepth__"]), numeric(1)))
  energy <- unlist(lapply(sampler, function(x) x[, "energy__"]))
  ebfmi <- if (length(energy) > 1L) mean(diff(energy)^2) / stats::var(energy) else NA_real_
  list(convergence = convergence, rhat = convergence$rhat, n_eff = convergence$n_eff,
       divergent = divergent, divergences = divergent, max_treedepth = treedepth,
       ebfmi = ebfmi,
       ok = !any(convergence$rhat > 1.01, na.rm = TRUE) && divergent == 0L)
}

.sft_bayes_posterior_predictive <- function(prepared, transformed, seed, n_draws = NULL) {
  set.seed(seed)
  all_draws <- dim(transformed$H_A)[1L]
  if (!is.null(n_draws) && length(n_draws) == 1L && is.finite(n_draws) && n_draws >= 1L) {
    n_draws <- min(all_draws, as.integer(n_draws))
    keep <- if (n_draws < all_draws) sample(seq_len(all_draws), n_draws) else seq_len(all_draws)
  } else {
    keep <- seq_len(all_draws)
  }
  draws <- length(keep)
  predicted <- lapply(seq_along(.sft_bayes_series), function(s) {
    # The likelihood is written in terms of hazards.  Recover interval rates
    # from cumulative hazards before simulating replicated event counts.
    hname <- .sft_bayes_series[[s]]
    Hs <- transformed[[paste0("H_", hname)]][keep, , , drop = FALSE]
    rates <- Hs
    rates[, , 1L] <- Hs[, , 1L] / prepared$grid$bins$width[[1L]]
    J <- dim(Hs)[3L]
    if (J > 1L) for (j in 2:J) {
      rates[, , j] <- (Hs[, , j] - Hs[, , j - 1L]) / prepared$grid$bins$width[[j]]
    }
    E <- matrix(prepared$exposure[, match(hname, .sft_bayes_series), ],
                nrow = prepared$n_subjects, ncol = J)
    out <- array(0L, dim = c(draws, nrow(E), J))
    for (m in seq_len(draws)) out[m, , ] <- matrix(
      stats::rpois(length(E), as.numeric(E * matrix(rates[m, , ],
                                                    nrow = nrow(E), ncol = J))),
      nrow = nrow(E), ncol = J)
    out
  })
  names(predicted) <- .sft_bayes_series
  summary <- do.call(rbind, lapply(seq_along(predicted), function(s) {
    vals <- apply(predicted[[s]], 1L, sum)
    interval <- .sft_hdi(vals)
    data.frame(Series = .sft_bayes_series[[s]], Observed = sum(prepared$events[, s, ]),
               Mean = mean(vals), Lower = interval[[1L]], Upper = interval[[2L]],
               stringsAsFactors = FALSE)
  }))
  list(draws = predicted, summary = summary)
}

.sft_bayes_empty_transformed <- function() {
  list(H_A = NULL, H_B = NULL, H_AB = NULL, D = NULL, C = NULL,
       subject = NULL, population = NULL, new_subject = NULL)
}

.sft_bayes_empty_curve_data <- function() {
  data.frame(Level = character(), Subject = character(), Series = character(),
             Measure = character(), Time = numeric(), Ct = numeric(),
             Mean = numeric(), Lower = numeric(), Upper = numeric(),
             Prob_super = numeric(), Prob_limited = numeric(), Prob_rope = numeric(),
             Report = logical(), stringsAsFactors = FALSE)
}

#' Fit an OR-centred hierarchical Bayesian capacity model.
#'
#' The input must contain Subject, Condition, RT, Correct, Channel1, and
#' Channel2.  By default positive salience values are pooled by channel
#' presence into AB, A, and B.  With `salience_split`, the model retains
#' A_L/A_H, B_L/B_H, and AB_LL/AB_LH/AB_HL/AB_HH, then hierarchically pools
#' matched capacity deviations across those cells.  Correct finite trials are
#' events; incorrect finite trials are censored exposure.  Fitting requires
#' the optional rstan package, but `sample = FALSE` still performs and returns
#' the complete preparation stage.
#'
#' @param data Canonical SFT data frame. `sftData` is accepted as an alias.
#' @param Condition Optional condition value or values to retain.
#' @param n_bins Requested number of pooled RT intervals.
#' @param rt_units Input RT units, either seconds or milliseconds.
#' @param report_quantiles Central pooled RT quantiles retained in tidy curves.
#' @param priors Named prior overrides; see the returned `prior_metadata`.
#' @param smoothness Named smoothness overrides, including `basis_dim`, `A`,
#'   `B`, and `delta`.
#' @param require_complete Require every subject to have AB, A, and B.
#' @param cell_mapping Optional function of Channel1 and Channel2 returning an
#'   OR series label per trial; it may return pooled labels or the salience-split
#'   labels documented below.
#' @param salience_split Retain low/high salience cells.  Use `TRUE` or
#'   `"auto"` when each channel has exactly two positive levels, a numeric
#'   vector of the two levels, a list with `Channel1` and `Channel2` mappings,
#'   or a function returning `L`/`H` labels.  `NULL` preserves pooled mode.
#' @param sample If TRUE, run the optional rstan sampler.
#' @param posterior_draws Optional precomputed log-hazard arrays for testing or
#'   downstream shared-engine use; bypasses rstan when supplied.
#' @param ... Compatibility aliases: `bins`, `nbin`, `grid_bins`, `rt.units`,
#'   `central_quantiles`, `prior`, `smooth`, `fit`, and `salience`.
#' @return An object of class `sft_bayes_capacity`.
#' @export
capacity.bayes <- function(data = NULL, Condition = NULL, n_bins = 25L,
                           rt_units = c("seconds", "milliseconds"),
                           report_quantiles = c(0.05, 0.95), priors = NULL,
                           smoothness = NULL, require_complete = TRUE,
                           cell_mapping = NULL, salience_split = NULL,
                           sample = TRUE, posterior_draws = NULL,
                           chains = 4L, iter = 2000L, warmup = NULL,
                           seed = NULL, cores = 1L, refresh = 0L,
                           control = list(adapt_delta = 0.95, max_treedepth = 12L),
                           hdi = 0.94, rope = 0, prior_predictive_draws = 200L,
                           posterior_predictive_draws = 200L, return_fit = TRUE,
                           sftData = NULL, ...) {
  dots <- list(...)
  aliases <- c("bins", "nbin", "grid_bins")
  for (nm in aliases) if (!is.null(dots[[nm]])) n_bins <- dots[[nm]]
  if (!is.null(dots$rt.units)) rt_units <- dots$rt.units
  if (!is.null(dots$central_quantiles)) report_quantiles <- dots$central_quantiles
  if (!is.null(dots$prior)) priors <- dots$prior
  if (!is.null(dots$smooth)) smoothness <- dots$smooth
  if (!is.null(dots$fit)) sample <- dots$fit
  if (!is.null(dots$salience)) salience_split <- dots$salience
  if (!is.null(dots$salience_mapping)) salience_split <- dots$salience_mapping
  unknown <- setdiff(names(dots), c(aliases, "rt.units", "central_quantiles", "prior", "smooth",
                                    "fit", "salience", "salience_mapping"))
  if (length(unknown)) stop("Unused argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  if (!is.null(sftData)) {
    if (!is.null(data)) stop("Supply either data or sftData, not both.", call. = FALSE)
    data <- sftData
  }
  if (is.null(data)) stop("data must be supplied.", call. = FALSE)
  rt_units <- .sft_bayes_units(rt_units)
  if (length(hdi) != 1L || !is.finite(hdi) || hdi <= 0 || hdi >= 1)
    stop("hdi must lie strictly between 0 and 1.", call. = FALSE)
  if (length(rope) != 1L || !is.finite(rope) || rope < 0)
    stop("rope must be a non-negative finite value.", call. = FALSE)
  seed <- .sft_capacity_bayes_seed(seed)
  prepared <- .sft_bayes_prepared(data, Condition, n_bins, report_quantiles,
                                  rt_units, priors, smoothness, require_complete,
                                  cell_mapping, salience_split)
  ppc_seed <- seed + 101L
  prior_ppc <- if (prepared$split)
    .sft_bayes_prior_predictive_salience(prepared, prior_predictive_draws, ppc_seed) else
    .sft_bayes_prior_predictive(prepared, prior_predictive_draws, ppc_seed)
  fit <- NULL; raw <- NULL; effects <- NULL
  transformed <- .sft_bayes_empty_transformed()
  diagnostics <- list(available = FALSE, reason = "No posterior fit requested.")
  posterior_ppc <- NULL
  if (!is.null(posterior_draws)) {
    if (!is.list(posterior_draws)) stop("posterior_draws must be a named list of log-hazard arrays.",
                                       call. = FALSE)
    if (prepared$split) {
      if (!all(prepared$series %in% names(posterior_draws)))
        stop("Salience-split posterior_draws must include all eight requested series.", call. = FALSE)
      raw <- setNames(lapply(prepared$series, function(nm) posterior_draws[[nm]]), prepared$series)
      names(raw) <- paste0("eta_", names(raw))
      effects <- .sft_bayes_effects_from_salience_posterior(raw, prepared$priors,
                                                             prepared$n_subjects, prepared$n_bins,
                                                             prepared$subjects, seed + 17L,
                                                             new_subject = FALSE)
      transformed <- .sft_bayes_transform_salience(effects$log_hazard,
                                                   prepared$grid$bins$width,
                                                   prepared$grid$bins$upper,
                                                   prepared$subjects)
      posterior_ppc <- .sft_bayes_predictive_from_series(prepared, effects$log_hazard,
                                                         seed + 29L, posterior_predictive_draws)
    } else {
      if (!all(c("A", "B", "AB") %in% names(posterior_draws)))
        stop("posterior_draws must be a list with A, B, and AB log-hazard arrays.", call. = FALSE)
      raw <- list(eta_A = posterior_draws$A, eta_B = posterior_draws$B,
                  eta_AB = posterior_draws$AB)
      if (isTRUE(all(c("population_A", "population_B", "population_delta") %in% names(posterior_draws)))) {
        raw$population_A <- posterior_draws$population_A
        raw$population_B <- posterior_draws$population_B
        raw$population_delta <- posterior_draws$population_delta
      }
      effects <- .sft_bayes_effects_from_posterior(raw, prepared$priors,
                                                     prepared$n_subjects, prepared$n_bins,
                                                     prepared$subjects, prepared$grid$bins$width,
                                                     seed + 17L, new_subject = FALSE)
      transformed <- .sft_bayes_transform(effects$log_hazard, prepared$grid$bins$width,
                                           prepared$grid$bins$upper, prepared$subjects)
      posterior_ppc <- .sft_bayes_posterior_predictive(prepared, transformed, seed + 29L,
                                                       posterior_predictive_draws)
    }
  } else if (isTRUE(sample)) {
    sampling <- .sft_bayes_validate_sampling(iter, warmup, chains, cores, refresh, control)
    model <- .sft_bayes_stan_model(prepared$split)
    fit <- rstan::sampling(model, data = prepared$stan_data, iter = sampling$iter,
                           warmup = sampling$warmup, chains = sampling$chains,
                           cores = sampling$cores, seed = seed, refresh = sampling$refresh,
                           control = sampling$control)
    pars <- if (prepared$split) c(
      "beta_A_L", "beta_A_H", "beta_B_L", "beta_B_H", "beta_delta", "beta_gamma_raw",
      "z_speed", "z_asymmetry", "z_capacity", "speed", "asymmetry", "capacity_shift",
      "eta_A_L", "eta_A_H", "eta_B_L", "eta_B_H", "eta_AB_LL", "eta_AB_LH",
      "eta_AB_HL", "eta_AB_HH", "population_A_L", "population_A_H",
      "population_B_L", "population_B_H", "population_delta", "population_gamma",
      "sigma_speed", "sigma_asymmetry", "sigma_capacity", "sigma_smooth_A",
      "sigma_smooth_B", "sigma_smooth_delta", "sigma_gamma", "sigma_smooth_gamma") else c(
        "beta_A", "beta_B", "beta_delta", "z_speed", "z_asymmetry",
        "z_capacity", "speed", "asymmetry", "capacity_shift", "eta_A",
        "eta_B", "eta_AB", "population_A", "population_B",
        "population_delta", "sigma_speed", "sigma_asymmetry", "sigma_capacity",
        "sigma_smooth_A", "sigma_smooth_B", "sigma_smooth_delta")
    raw <- rstan::extract(fit, pars = pars, permuted = TRUE, inc_warmup = FALSE)
    if (prepared$split) {
      effects <- .sft_bayes_effects_from_salience_posterior(raw, prepared$priors,
                                                             prepared$n_subjects, prepared$n_bins,
                                                             prepared$subjects, seed + 17L,
                                                             new_subject = TRUE)
      transformed <- .sft_bayes_transform_salience(effects$log_hazard,
                                                   prepared$grid$bins$width,
                                                   prepared$grid$bins$upper,
                                                   prepared$subjects, effects$new_log_hazard)
      posterior_ppc <- .sft_bayes_predictive_from_series(prepared, effects$log_hazard,
                                                         seed + 29L, posterior_predictive_draws)
    } else {
      effects <- .sft_bayes_effects_from_posterior(raw, prepared$priors,
                                                     prepared$n_subjects, prepared$n_bins,
                                                     prepared$subjects, prepared$grid$bins$width,
                                                     seed + 17L, new_subject = TRUE)
      transformed <- .sft_bayes_transform(effects$log_hazard, prepared$grid$bins$width,
                                           prepared$grid$bins$upper, prepared$subjects,
                                           effects$new_log_hazard)
      posterior_ppc <- .sft_bayes_posterior_predictive(prepared, transformed, seed + 29L,
                                                       posterior_predictive_draws)
    }
    diagnostics <- .sft_bayes_diagnostics(fit)
    diagnostics$available <- TRUE
  }
  if (!is.null(transformed$D)) {
    if (prepared$split) {
      transformed$population <- list(
        D = apply(transformed$D, c(1L, 3L), mean),
        C = exp(apply(log(transformed$C), c(1L, 3L), mean, na.rm = TRUE)),
        logC = apply(transformed$logC_average, c(1L, 3L), mean, na.rm = TRUE),
        D_cells = apply(transformed$D_cells, c(1L, 3L, 4L), mean),
        C_cells = apply(transformed$C_cells, c(1L, 3L, 4L), function(x)
          exp(mean(log(x), na.rm = TRUE))),
        SIC = apply(transformed$SIC, c(1L, 3L), mean)
      )
      summaries <- .sft_bayes_split_summaries(transformed, prepared$grid,
                                               prepared$subjects, hdi, rope)
    } else {
      transformed$population <- list(
        D = apply(transformed$D, c(1L, 3L), mean),
        C = apply(transformed$C, c(1L, 3L), mean),
        H_A = apply(transformed$H_A, c(1L, 3L), mean),
        H_B = apply(transformed$H_B, c(1L, 3L), mean),
        H_AB = apply(transformed$H_AB, c(1L, 3L), mean)
      )
      summaries <- .sft_bayes_summaries(transformed, prepared$grid, prepared$subjects, hdi, rope)
    }
  } else {
    empty <- .sft_bayes_empty_curve_data()
    summaries <- list(tidy = empty, subject = empty, population = empty,
                      new_subject = empty)
  }
  posterior <- if (is.null(raw)) NULL else list(
    draws = raw,
    log_hazard = if (!is.null(effects)) effects$log_hazard else NULL,
    n_draws = if (!is.null(transformed$D)) dim(transformed$D)[1L] else 0L
  )
  # The direct log-hazard arrays are more useful and exact than attempting to
  # reconstruct them from cumulative hazards.  Keep them in a compact alias.
  if (!is.null(effects)) posterior$log_hazard <- effects$log_hazard
  if (!is.null(effects)) {
    transformed$log_hazard <- effects$log_hazard
    transformed$hazard <- lapply(effects$log_hazard, exp)
  }
  tidy <- summaries$tidy
  out <- list(
    call = match.call(), method = "OR-centred hierarchical piecewise-exponential",
    stopping_rule = "OR", data = prepared$trials, prepared = prepared,
    grid = prepared$grid, posterior = posterior, fit = if (isTRUE(return_fit)) fit else NULL,
    transformed = transformed, curves = summaries, tidy_curves = tidy,
    curve_data = tidy, diagnostics = diagnostics,
    predictive_checks = list(prior = prior_ppc, posterior = posterior_ppc),
    prior_metadata = list(priors = prepared$priors,
                          units = "log hazards are calibrated in seconds",
                          smoothness = list(type = "second-difference random walk",
                                            basis_dimension = ncol(prepared$basis),
                                            salience_effect = if (prepared$split)
                                              "sum-to-zero gamma cells with strong shrinkage" else NULL),
                          central_interval = prepared$grid$report_range,
                          hdi = hdi, rope = rope),
    salience_pooling = prepared$pooling,
    metadata = list(series = prepared$series, subjects = prepared$subjects,
                    selected_conditions = prepared$selected_conditions,
                    salience_pooling = prepared$pooling,
                    incorrect_trials = "finite censored exposure with no event",
                    complete_subject_series = prepared$complete,
                    missing_subject_series = prepared$missing_series,
                    report_quantiles = report_quantiles,
                    rt_units = prepared$rt_units,
                    posterior_summary = "all capacity summaries are draw-by-draw",
                    sic = if (prepared$split)
                      "posterior SIC = S_AB_LL - S_AB_LH - S_AB_HL + S_AB_HH" else
                      "not identified without an explicit salience split"),
    stan = list(path = tryCatch(.sft_bayes_stan_path(prepared$split), error = function(e) NA_character_),
                data = prepared$stan_data),
    settings = list(n_bins = n_bins, iter = iter, warmup = warmup,
                    chains = chains, seed = seed, sample = isTRUE(sample),
                    return_fit = isTRUE(return_fit))
  )
  class(out) <- c("sft_bayes_capacity", "list")
  out
}

print.sft_bayes_capacity <- function(x, ...) {
  cat("OR-centred hierarchical Bayesian capacity model\n")
  cat("Subjects:", length(x$prepared$subjects), " Series:", paste(x$prepared$series, collapse = ", "),
      " Bins:", x$prepared$n_bins, "\n")
  cat("RT units:", x$prepared$rt_units, "  Conditions:",
      paste(x$prepared$selected_conditions, collapse = ", "), "\n")
  if (isTRUE(x$diagnostics$available)) {
    cat("Posterior draws:", x$posterior$n_draws, "  Divergences:", x$diagnostics$divergences, "\n")
  } else cat("No posterior fit stored; preparation/prior predictive results are available.\n")
  invisible(x)
}
