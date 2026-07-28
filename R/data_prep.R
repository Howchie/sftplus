# Canonical SFT data adapters and response-time curve/ECDF utilities.

# Core capacity, UCIP, SIC, resilience, Bayesian, plotting, and assessment
# implementations for sft.plus.  Public functions in this file replace the
# corresponding upstream implementations; each public name has one canonical
# definition in the package.

.sft_zero_curve <- function(value = 0) {
  force(value)
  function(t) rep(value, length(t))
}


.sft_curve <- function(x, y, default = NA_real_, rule = 1L) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  # Only non-finite times are dropped.  Non-finite y values are kept as NA
  # markers: the callers use them to record where an estimate is undefined
  # (a zero denominator, a log of zero), and interpolating across them would
  # invent a curve in exactly the region the caller flagged as unusable.
  ok <- is.finite(x)
  if (!any(ok)) return(.sft_zero_curve(default))
  x <- x[ok]
  y <- y[ok]
  y[!is.finite(y)] <- NA_real_
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (!any(!is.na(y))) return(.sft_zero_curve(default))
  if (length(x) == 1L) return(.sft_zero_curve(y[[1L]]))
  # rule = 1 returns NA outside the estimated support; na.rm = FALSE lets the
  # NA markers propagate into the intervals adjoining them.
  stats::approxfun(x, y, rule = rule, ties = "ordered", na.rm = FALSE)
}


.sft_probs <- function(probs) {
  probs <- sort(unique(as.numeric(probs)))
  if (!length(probs) || any(!is.finite(probs)) || any(probs <= 0) ||
      any(probs >= 1)) {
    stop("probs must be finite percentile levels strictly inside (0, 1).",
         call. = FALSE)
  }
  probs
}


.sft_time_grid <- function(rt, n.times = 1000L) {
  # A scale-free replacement for round(rt): keep the observed times when there
  # are few enough of them, otherwise thin to an evenly spaced grid spanning
  # the same range.  Rounding to the nearest integer collapsed second-scale
  # response times to a handful of distinct values.
  x <- as.numeric(rt)
  x <- sort(unique(x[is.finite(x)]))
  if (!length(x)) stop("RT contains no finite response times.")
  n.times <- as.integer(n.times)
  if (is.na(n.times) || n.times < 2L || length(x) <= n.times) return(x)
  seq(x[1L], x[length(x)], length.out = n.times)
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


# Canonical names used by the row-wise SFT APIs and the aliases accepted when
# data come from external workflows.  Keep this mapping in one place so every
# data-frame entry point (frequentist, Bayesian, and conversion helpers) has
# identical column-name behaviour.
.sft_column_aliases <- list(
  RT = c("rt"),
  Subject = c("subjects"),
  Condition = c("LogicalRule"),
  Block = c("block", "Blocks", "blocks", "Session", "session")
)


# Rename aliases to their canonical names when the canonical column is absent.
# Canonical names take precedence when both names are supplied; retaining the
# alias column in that case is intentional so this helper never discards user
# data.  Multiple aliases for the same missing canonical column are rejected
# rather than selected arbitrarily.
.sft_normalize_columns <- function(data, aliases = .sft_column_aliases) {
  if (!is.data.frame(data)) return(data)
  if (is.null(names(data)) || anyNA(names(data)) || any(!nzchar(names(data)))) {
    stop("data must have non-empty column names.", call. = FALSE)
  }
  for (canonical in names(aliases)) {
    if (canonical %in% names(data)) next
    hits <- intersect(aliases[[canonical]], names(data))
    if (length(hits) > 1L) {
      stop("data contains multiple aliases for the '", canonical,
           "' column: ", paste(hits, collapse = ", "), call. = FALSE)
    }
    if (length(hits) == 1L) {
      names(data)[match(hits[[1L]], names(data))] <- canonical
    }
  }
  data
}


#' Convert canonical row-wise SFT data to RT/CR list input.
#'
#' Converts a trial-level data frame to one response-time vector and one
#' correctness vector per experimental cell. For OR and AND stopping rules it
#' returns AB, A, and B, using positive channel values as target-present. For
#' STST it returns
#' context and target cells using the package's signed-channel convention.
#' With multiple subjects, the default result is nested by subject so it can
#' be passed to `capacityGroup.bayes()`; select one subject or use
#' `by_subject = "never"` for single-subject functions.
#'
#' Rows with both channels off are omitted from OR and AND output unless
#' `include_nn = TRUE`. Under the STST convention, context rows have one
#' positive and one negative channel; target rows have one nonzero,
#' non-negative channel.
#'
#' @param data A data frame containing `Subject`, `RT`, `Correct`, `Channel1`,
#'   and `Channel2`; `Condition` is optional when only one condition is used.
#'   The aliases `subjects`, `rt`, and `LogicalRule` are accepted for
#'   `Subject`, `RT`, and `Condition`, respectively.
#' @param Condition Optional single condition value to retain.  If omitted,
#'   the data must contain one condition (or no `Condition` column).
#' @param Subject Optional subject value or values to retain.
#' @param stopping.rule One of `"OR"`, `"AND"`, or `"STST"`.
#' @param by_subject `"auto"` returns a flat list for one subject and nested
#'   lists for multiple subjects; `"always"` or `"never"` forces either form.
#' @param include_nn Include the both-off `NN` cell for OR/AND conversion.
#' @return A list with `RT` and matching `CR` lists.  Metadata are stored as
#'   attributes, including selected subjects, condition, and dropped rows.
#' @examples
#' one <- sft_data_to_rt(dots, Condition = "OR", Subject = "S1")
#' names(one$RT)
#'
#' group <- sft_data_to_rt(dots, Condition = "OR")
#' length(group$RT)
#' @export
sft_data_to_rt <- function(data, Condition = NULL, Subject = NULL,
                           stopping.rule = c("OR", "AND", "STST"),
                           by_subject = c("auto", "always", "never"),
                           include_nn = FALSE) {
  data <- .sft_normalize_columns(data)
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
