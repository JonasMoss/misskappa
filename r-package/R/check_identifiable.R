# Identifiability / estimability guard for the observed missing-data pattern.
#
# A saturated agreement coefficient -- Cronbach's alpha, or Conger / Cohen /
# Brennan-Prediger kappa -- is a functional of *every* pairwise second moment, so
# it is identified from the observed data only if every column pair is jointly
# observed by at least one subject: the co-observation graph must be COMPLETE, not
# merely connected. Connectivity suffices only for structured / exchangeable
# coefficients (a one-factor covariance, or Fleiss under exchangeable raters).
# Companion-paper notes keep the population-identification discussion separate
# from this package-level estimator guard.
#
# FIML does not fail on a violation: the saturated EM returns an
# initialisation-pinned value in a flat, rank-deficient direction -- more
# dangerous than the NA pairwise deletion returns -- and its pseudo-inverse
# silently absorbs that direction. So this check gates *before* the estimator
# runs and errors for every estimator. The counts interface (kappa_counts) is the
# exchangeable representation -- rater identity is discarded -- so the
# complete-graph condition is vacuous there and this guard is not applied to it.

# Logical subjects-by-columns observation matrix, or NULL when x is not a usable
# numeric matrix / data frame (let the downstream backend raise the type error).
.pattern_observed <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) return(NULL)
  xm <- suppressWarnings(as.matrix(x))
  if (!is.numeric(xm)) return(NULL)
  is.finite(xm)
}

# Format up to `cap` offending column tuples as "a-b, c-d" using column names.
.fmt_tuples <- function(idx, nm, cap = 6L) {
  if (is.null(dim(idx))) idx <- matrix(idx, nrow = 1L)
  labs <- apply(idx, 1L, function(row) paste(nm[row], collapse = "-"))
  if (length(labs) > cap) {
    paste0(paste(labs[seq_len(cap)], collapse = ", "),
           sprintf(", (and %d more)", length(labs) - cap))
  } else {
    paste(labs, collapse = ", ")
  }
}

# Gate a fit on the identifiability of the requested coefficient given `observed`
# (logical subjects-by-columns, from .pattern_observed). `require = "complete"`
# enforces complete arity-wise co-observation for the saturated coefficients;
# `"each_subject_2"` is the exchangeable relaxation (each subject has >= 2
# observed columns). For `arity = 2`, returns the pairwise co-observation count
# matrix invisibly. Errors with call. = FALSE so the message reads as a
# user-facing precondition, not an internal trace.
.check_pattern_identifiable <- function(observed,
                                        unit = c("item", "rater"),
                                        require = c("complete", "each_subject_2"),
                                        coefficient = "the coefficient",
                                        arity = 2L) {
  unit <- match.arg(unit)
  require <- match.arg(require)
  if (is.null(observed)) return(invisible(NULL))

  obs <- matrix(as.numeric(observed), nrow(observed), ncol(observed),
                dimnames = dimnames(observed))
  p <- ncol(obs)
  nm <- colnames(obs)
  if (is.null(nm)) nm <- as.character(seq_len(p))

  # A pattern with fewer than two columns carries no co-observation constraint;
  # defer to the estimator's own "needs at least two columns" validation rather
  # than pre-empting it with the arity-range error below.
  if (p < 2L) return(invisible(NULL))

  if (!is.numeric(arity) || length(arity) != 1L || !is.finite(arity) ||
      arity != round(arity) || arity < 2L || arity > p) {
    stop("'arity' must be an integer between 2 and the number of columns.",
         call. = FALSE)
  }
  arity <- as.integer(arity)

  # Every column observed at least once (subsumes, with a clearer message, the
  # case of a never-seen column whose every pair would otherwise be flagged).
  if (any(colSums(obs) == 0)) {
    stop(sprintf("every %s must be observed for at least one subject.", unit),
         call. = FALSE)
  }

  if (require == "each_subject_2") {
    if (any(rowSums(obs) < 2)) {
      stop(sprintf("every subject must have at least two observed %ss.", unit),
           call. = FALSE)
    }
    return(invisible(crossprod(obs)))
  }

  # Saturated pairwise functional: complete co-observation graph
  # (every n_jk >= 1).
  N <- crossprod(obs)                 # p-by-p co-observation counts n_jk
  if (arity > 2L) {
    tuples <- t(utils::combn(p, arity))
    counts <- apply(tuples, 1L, function(cols) {
      sum(rowSums(obs[, cols, drop = FALSE]) == arity)
    })
    zero <- which(counts == 0L)
    if (length(zero) > 0L) {
      stop(sprintf(paste0(
        "%s %s-tuple(s) %s never jointly observed, so %s is not identified ",
        "from this missing-data pattern. Drop an offending %s, or provide ",
        "data in which every %s %s-tuple is co-observed."),
        unit, arity, .fmt_tuples(tuples[zero, , drop = FALSE], nm),
        coefficient, unit, unit, arity), call. = FALSE)
    }
    thin <- which(counts == 1L)
    if (length(thin) > 0L) {
      warning(sprintf(paste0(
        "%s %s-tuple(s) %s co-observed by only one subject; the corresponding ",
        "observed-disagreement term is degenerate and the standard error ",
        "unreliable."),
        unit, arity, .fmt_tuples(tuples[thin, , drop = FALSE], nm)),
        call. = FALSE)
    }
    return(invisible(counts))
  }

  ut <- upper.tri(N)
  zero <- which(N == 0 & ut, arr.ind = TRUE)
  if (nrow(zero) > 0L) {
    stop(sprintf(paste0(
      "%s pair(s) %s never jointly observed, so %s is not identified from this ",
      "missing-data pattern (the co-observation graph is incomplete). Drop an ",
      "offending %s, or provide data in which every %s pair is co-observed."),
      unit, .fmt_tuples(zero, nm), coefficient, unit, unit), call. = FALSE)
  }
  thin <- which(N == 1 & ut, arr.ind = TRUE)
  if (nrow(thin) > 0L) {
    warning(sprintf(paste0(
      "%s pair(s) %s co-observed by only one subject; the corresponding pairwise ",
      "covariance is degenerate and the standard error unreliable."),
      unit, .fmt_tuples(thin, nm)), call. = FALSE)
  }
  invisible(N)
}
