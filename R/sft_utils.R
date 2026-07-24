# Small cross-cutting helpers shared across the sftplus internals.


.sft_htest <- function(statistic, p.value, alternative, method, data.name = NULL) {
  ans <- list(statistic = statistic, p.value = p.value,
              alternative = alternative, method = method)
  if (!is.null(data.name)) ans$data.name <- data.name
  class(ans) <- "htest"
  ans
}


`%||%` <- function(x, y) if (is.null(x)) y else x


.sft_hdi <- function(x, mass = .94) {
  x <- sort(as.numeric(x[is.finite(x)]))
  if (!length(x)) return(c(lower = NA_real_, upper = NA_real_))
  n <- length(x)
  width <- max(0L, ceiling(mass * n) - 1L)
  if (width == 0L) return(c(lower = x[1L], upper = x[1L]))
  left <- seq_len(n - width)
  j <- left[which.min(x[left + width] - x[left])]
  c(lower = x[j], upper = x[j + width])
}


.to_correct_indicator <- function(cr, n) {
  if (is.null(cr)) return(rep(TRUE, n))
  out <- rep(FALSE, n); m <- min(length(cr), n)
  if (m) {
    v <- as.logical(cr[seq_len(m)]); v[is.na(v)] <- FALSE; out[seq_len(m)] <- v
  }
  out
}


.sft_bind_rows <- function(parts) {
  parts <- parts[vapply(parts, is.data.frame, logical(1))]
  if (!length(parts)) return(data.frame())
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(cols, names(x)); x[missing] <- NA
    x[cols]
  })
  do.call(rbind, parts)
}

