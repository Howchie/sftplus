generateData <- function(nSubjects, nTrials, arch, stoprule, highlowratio=3) {
  nu.l <- 0.1
  nu.h <- nu.l * highlowratio
  alpha <- 30
  diff <- 1
  sop <- .5

  rData <- matrix(NA, 4 * nSubjects * nTrials, 6)
  rData[,1] <- rep(1:nSubjects, each=4*nTrials)
  rData[,2] <- 1
  rData[,3] <- rep(c(rep(2,2*nTrials), rep(1, 2*nTrials)), nSubjects)
  rData[,4] <- rep(c(rep(2,nTrials), rep(1, nTrials)), 2*nSubjects)

  allrt <- numeric()
  for( subj in 1:nSubjects) {
    if(arch=="coactive") {
      rt.hh <- rinvGauss(nTrials, nu=alpha/(nu.h+nu.h), lambda=.5*(alpha/diff)^2)
      rt.hl <- rinvGauss(nTrials, nu=alpha/(nu.h+nu.l), lambda=.5*(alpha/diff)^2)
      rt.lh <- rinvGauss(nTrials, nu=alpha/(nu.l+nu.h), lambda=.5*(alpha/diff)^2)
      rt.ll <- rinvGauss(nTrials, nu=alpha/(nu.l+nu.l), lambda=.5*(alpha/diff)^2)

    } else {
      # Channel1 High
      x1h1 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)
      x1h2 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)

      # Channel1 Low
      x1l1 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)
      x1l2 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)

      # Channel2 High
      x2h1 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)
      x2h2 <- rinvGauss(nTrials, nu=alpha/nu.h, lambda=(alpha/diff)^2)

      # Channel2 Low
      x2l1 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)
      x2l2 <- rinvGauss(nTrials, nu=alpha/nu.l, lambda=(alpha/diff)^2)

      if (arch == "parallel") {
        if (stoprule == "and") {
          rt.hh <- pmax(x1h1, x2h1)
          rt.hl <- pmax(x1h2, x2l1)
          rt.lh <- pmax(x1l1, x2h2)
          rt.ll <- pmax(x1l2, x2l2)
        }else if (stoprule == "or") {
          rt.hh <- pmin(x1h1, x2h1)
          rt.hl <- pmin(x1h2, x2l1)
          rt.lh <- pmin(x1l1, x2h2)
          rt.ll <- pmin(x1l2, x2l2)
        } else {
          cat("Unknown stopping rule!\n")
          return(NULL)
        }
      } else if (arch == "serial") {
        if (stoprule == "and") {
          rt.hh <- x1h1 + x2h1
          rt.hl <- x1h2 + x2l1
          rt.lh <- x1l1 + x2h2
          rt.ll <- x1l2 + x2l2
        }else if (stoprule == "or") {
          oneFirst <- runif(nTrials) < .5
          rt.hh <- oneFirst * x1h1 + (1-oneFirst)*x2h1
          rt.hl <- oneFirst * x1h2 + (1-oneFirst)*x2l1
          rt.lh <- oneFirst * x1l1 + (1-oneFirst)*x2h2
          rt.ll <- oneFirst * x1l2 + (1-oneFirst)*x2l2
        } else {
          cat("Unknown stopping rule!\n")
          return(NULL)
        }
      } else {
        cat("Unknown architecture!\n")
        return(NULL)
      }
    }
    allrt <- c(allrt, rt.hh, rt.hl, rt.lh, rt.ll)
  }
  rData[,5] <- allrt
  rData[,6] <- TRUE
  rData <- as.data.frame(rData)
  rData[,2] <- paste(arch,stoprule, sep="-")
  names(rData) <- c("Subject", "Condition", "Channel1", "Channel2", "RT", "Correct")
  rData$Subject <- as.factor(rData$Subject)
  rData$Condition <- as.factor(rData$Condition)
  return(rData)
}


# ---------------------------------------------------------------------------
# Multi-subject LBA/OU generation
#
# `simulate_sft()` (see simulation.R) is the canonical single-subject LBA/OU
# simulator. `simulate_sft_group()` wraps it to generate `nSubjects` data sets
# from per-subject process parameters, then row-binds them into one data frame
# that is directly compatible with the package's group-level functions
# (capacityGroup(), sicGroup(), assessmentGroup(), ...). Single-subject use of
# simulate_sft() is unchanged.
# ---------------------------------------------------------------------------

# Coerce the user-supplied `params` into a subjects x parameters matrix with
# named columns (the parameter names expected in a simulate_sft() `p_vec`).
#   * a matrix / data.frame  -> one row per subject (names taken from columns)
#   * a named numeric vector -> a single p_vec, recycled to `nSubjects` rows
.sft_subject_param_matrix <- function(params, nSubjects = NULL) {
  if (is.data.frame(params)) params <- as.matrix(params)
  if (is.matrix(params)) {
    if (is.null(colnames(params))) {
      stop("`params` matrix must have named columns (the simulation parameter names).")
    }
    if (!is.null(nSubjects) && as.integer(nSubjects) != nrow(params)) {
      stop("nSubjects (", nSubjects, ") does not match the number of rows in `params` (",
           nrow(params), ").")
    }
    storage.mode(params) <- "double"
    return(params)
  }
  if (!is.numeric(params) || is.null(names(params)) || any(!nzchar(names(params)))) {
    stop("`params` must be a named numeric vector or a matrix/data.frame with named columns.")
  }
  N <- if (is.null(nSubjects)) 1L else as.integer(nSubjects)
  if (length(N) != 1L || is.na(N) || N < 1L) stop("nSubjects must be a positive integer.")
  matrix(as.numeric(params), nrow = N, ncol = length(params), byrow = TRUE,
         dimnames = list(NULL, names(params)))
}

# Resolve subject identifiers. Explicit `subject_ids` win, then row names of the
# parameter matrix, then auto-generated zero-padded "SubjectXXX" labels (1..N).
.sft_subject_ids <- function(subject_ids, N, row_names = NULL) {
  if (!is.null(subject_ids)) {
    if (length(subject_ids) != N) {
      stop("length(subject_ids) must equal the number of subjects (", N, ").")
    }
    return(as.character(subject_ids))
  }
  if (!is.null(row_names) && !anyNA(row_names) && all(nzchar(row_names))) {
    return(as.character(row_names))
  }
  width <- max(3L, nchar(as.character(N)))
  paste0("Subject", formatC(seq_len(N), width = width, flag = "0"))
}

simulate_sft_group <- function(model = c("lba", "ou"), n, params,
                               nSubjects = NULL, subject_ids = NULL,
                               combine = TRUE, ...) {
  model <- match.arg(model)
  pmat <- .sft_subject_param_matrix(params, nSubjects)
  N <- nrow(pmat)
  ids <- .sft_subject_ids(subject_ids, N, rownames(pmat))
  if (anyDuplicated(ids)) stop("Subject identifiers must be unique.")

  results <- vector("list", N)
  for (i in seq_len(N)) {
    p_vec <- stats::setNames(as.numeric(pmat[i, ]), colnames(pmat))
    results[[i]] <- simulate_sft(model = model, n = n, p_vec = p_vec,
                                 subject = ids[[i]], ...)
  }
  names(results) <- ids
  if (!combine) return(results)

  data <- .sft_bind_rows(lapply(results, function(r) r$data))
  data$Subject <- factor(as.character(data$Subject), levels = ids)
  data$Condition <- as.factor(data$Condition)
  rownames(data) <- NULL
  data
}
