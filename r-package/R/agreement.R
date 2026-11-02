# Attach exact complete-data delete-one evaluations without changing the
# estimator object when the requested backend has a different point
# functional. The latter guard matters for structured likelihood estimators.
.attach_delete_one <- function(fit, out, coefficient_names) {
  if (is.null(out) || length(out$estimate) != length(coefficient_names)) {
    return(fit)
  }
  names(out$estimate) <- coefficient_names
  colnames(out$leave_one_out) <- coefficient_names
  if (!is.null(out$chance) && length(out$chance) > 0L) {
    names(out$chance) <- coefficient_names
    colnames(out$chance_leave_one_out) <- coefficient_names
  }
  target <- unname(fit$estimates[coefficient_names])
  scale <- max(1, abs(target), abs(out$estimate), na.rm = TRUE)
  if (any(!is.finite(target)) || any(!is.finite(out$estimate)) ||
      max(abs(target - out$estimate)) > 1e-7 * scale) {
    return(fit)
  }
  fit$delete_one <- out
  fit
}

.attach_raw_kappa_delete_one <- function(fit, x, weight, values,
                                         em_options = list(),
                                         categories = NULL) {
  if (nrow(x) < 3L || anyNA(x) || any(!is.finite(x))) return(fit)
  flatten <- em_options[["flatten"]]
  if (identical(fit$method, "fiml") && !is.null(flatten) && flatten > 0) {
    return(fit)
  }
  out <- rcpp_kappa_delete_one(x, weight, values, categories)
  .attach_delete_one(fit, out, c("Conger", "Fleiss", "Brennan-Prediger"))
}

.attach_continuous_kappa_delete_one <- function(fit, x, weight) {
  if (nrow(x) < 3L || !all(is.finite(x))) return(fit)
  out <- rcpp_kappa_continuous_delete_one(x, weight)
  .attach_delete_one(fit, out, c("Conger", "Fleiss"))
}

.attach_counts_kappa_delete_one <- function(fit, x, weight, values,
                                            em_options = list()) {
  if (nrow(x) < 3L || anyNA(x) || any(!is.finite(x))) return(fit)
  flatten <- em_options[["flatten"]]
  if (identical(fit$method, "cat_fiml") && !is.null(flatten) && flatten > 0) {
    return(fit)
  }
  out <- rcpp_kappa_counts_delete_one(x, weight, values)
  .attach_delete_one(fit, out, c("Fleiss", "Brennan-Prediger"))
}

# Resolve the rating scale of a subjects-by-raters matrix of labels and code
# the ratings as integers 1..C. The scale is, in order of precedence, the
# declared `categories`, the shared levels of an all-factor data frame, or the
# sorted observed labels. Labels may be numeric, character, logical, or factor;
# numeric labels are matched exactly, so 0.5 and 1 are distinct categories.
# `values` scores the categories (default: numeric labels themselves, else
# 1..C). `declared` records whether the user fixed the scale, which is what
# Brennan-Prediger needs for the scored estimators.
.code_ratings <- function(x, categories = NULL, values = NULL) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  factor_levels <- NULL
  if (is.data.frame(x)) {
    is_fac <- vapply(x, is.factor, logical(1))
    if (length(is_fac) && all(is_fac)) {
      lv <- lapply(x, levels)
      if (all(vapply(lv, identical, logical(1), lv[[1L]]))) {
        factor_levels <- lv[[1L]]
      }
    }
    x <- as.matrix(x)
  }
  # Numeric and logical labels carry their own order; text labels are ordered
  # only through declared categories, values, or shared factor levels.
  ordered <- is.numeric(x) || is.logical(x) || !is.null(categories) ||
    !is.null(values) || !is.null(factor_levels)
  if (!is.numeric(x) && !is.character(x) && !is.logical(x)) {
    stop("'x' must hold numeric, character, logical, or factor ratings.",
         call. = FALSE)
  }
  if (is.numeric(x)) x[!is.finite(x)] <- NA
  observed <- unique(x[!is.na(x)])
  if (!length(observed)) stop("all ratings are missing.", call. = FALSE)

  declared <- !is.null(categories) || !is.null(values) ||
    !is.null(factor_levels)
  if (is.null(categories)) {
    categories <- if (!is.null(factor_levels)) factor_levels else sort(observed)
  } else {
    if (is.factor(categories)) categories <- as.character(categories)
    if (!is.atomic(categories) || length(categories) < 2L ||
        anyNA(categories) || anyDuplicated(categories)) {
      stop("'categories' must list at least two distinct, non-missing labels.",
           call. = FALSE)
    }
    extra <- observed[is.na(match(observed, categories))]
    if (length(extra)) {
      shown <- extra[seq_len(min(5L, length(extra)))]
      stop("Observed ratings absent from 'categories': ",
           paste(shown, collapse = ", "), if (length(extra) > 5L) ", ...",
           call. = FALSE)
    }
  }
  C <- length(categories)
  if (is.null(values)) {
    values <- if (is.numeric(categories)) as.numeric(categories) else seq_len(C)
  } else {
    if (!is.numeric(values) || any(!is.finite(values))) {
      stop("'values' must be a finite numeric vector.", call. = FALSE)
    }
    if (length(values) != C) {
      stop(sprintf(paste0(
        "Length of 'values' (%d) must equal the number of categories (%d). ",
        "Without 'categories' these are the distinct observed ratings; pass ",
        "'categories' to declare the full scale."), length(values), C),
        call. = FALSE)
    }
    values <- as.numeric(values)
  }
  codes <- matrix(match(x, categories), nrow(x), ncol(x),
                  dimnames = dimnames(x))
  storage.mode(codes) <- "integer"
  list(codes = codes, categories = categories, values = values,
       declared = declared, ordered = ordered)
}

# Linear and quadratic weights, and the scored estimators, depend on the order
# of the categories. Refuse to invent one (alphabetical) for text labels.
.require_ordered <- function(r, what) {
  if (!isTRUE(r$ordered)) {
    stop(what, " needs ordered categories, but the ratings are text with no ",
         "declared order. Pass 'categories' in scale order, or factor columns ",
         "with shared levels.", call. = FALSE)
  }
  invisible(r)
}

# Numeric scores for the scored (quadratic) estimators. Plain numeric input
# with no declared scale is used as is; otherwise the ratings are coded and
# replaced by their category scores, and the scale counts as declared.
.scored_ratings <- function(x, categories = NULL, values = NULL) {
  if (is.null(categories) && is.null(values)) {
    X <- if (is.data.frame(x) &&
             all(vapply(x, is.numeric, logical(1)))) as.matrix(x) else x
    if (is.matrix(X) && is.numeric(X)) {
      storage.mode(X) <- "double"
      return(list(X = X, values = NULL, declared = FALSE))
    }
  }
  r <- .code_ratings(x, categories, values)
  .require_ordered(r, "A scored (quadratic) estimator")
  X <- matrix(r$values[r$codes], nrow(r$codes), ncol(r$codes),
              dimnames = dimnames(r$codes))
  list(X = X, values = r$values, declared = r$declared,
       categories = r$categories)
}

# Categorical raw-rating fit for any C++ method ("ipw", "fiml", and the
# research methods "available", "gwet"). agreement() calls it with "ipw" or
# "fiml"; research_agreement_raw() exposes the rest.
.kappa_raw_fit <- function(x, method, weight = "identity", values = NULL,
                           em_options = list(), categories = NULL) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.", call. = FALSE)
  if (any(is.finite(x_mat) & x_mat != round(x_mat))) {
    stop("category codes must be whole numbers.", call. = FALSE)
  }
  storage.mode(x_mat) <- "integer"
  if (!is.null(categories)) categories <- as.integer(categories)

  out <- rcpp_kappa_raw(x = x_mat, method = method, weight_type = weight,
                        values = values, em_options = em_options,
                        categories = categories)
  fit <- .new_estimate(out, c("Conger", "Fleiss", "Brennan-Prediger"),
                       method = method, weight = weight)
  .attach_raw_kappa_delete_one(fit, x_mat, weight, values, em_options,
                               categories)
}

#' Weighted agreement coefficients with missing data
#'
#' @description
#' Estimates Conger, Fleiss, and Brennan-Prediger agreement coefficients (with
#' two raters, Conger is Cohen's kappa) from a subjects-by-raters matrix in
#' which some ratings may be missing. Choose the estimator with `estimator`:
#' \itemize{
#'   \item `"ipw"` (default, recommended): inverse-probability weighting.
#'     Consistent under MCAR, also with raters observed at different rates;
#'     works with any `weight`.
#'   \item `"cat_fiml"`: categorical full-information maximum likelihood.
#'     Consistent under MAR; feasible for few raters and categories.
#'   \item `"pairwise"`: moment estimator for quadratic weights, MCAR.
#'   \item `"nt_fiml"`: robust normal-theory FIML for quadratic weights,
#'     suited to MAR with roughly normal scores. Falls back to `"pairwise"`
#'     when the fit fails; see the section on the fallback.
#'   \item `"nt_fiml_cs"`: `"nt_fiml"` under the compound-symmetry covariance
#'     \eqn{\Sigma = \psi J + e I}; needs only a connected rater co-observation
#'     graph, so it handles planned-incomplete designs.
#' }
#' The last three treat ratings as numeric scores and use the quadratic
#' weight.
#'
#' @section Rating scale:
#' Ratings may be numbers, strings, logicals, or factors. The scale defaults
#' to the shared factor levels, or else the sorted observed ratings. Text
#' ratings have no order of their own, so linear and quadratic weights and the
#' scored estimators need `categories` (or factors with shared levels) for
#' them. Pass
#' `categories` to declare the full scale including categories nobody used
#' (`categories = 1:5` on a 1-5 scale where no one chose 5). The scale sets
#' Brennan-Prediger's chance baseline and the range that scales the linear and
#' quadratic weights. `"nt_fiml"` and `"nt_fiml_cs"` report Brennan-Prediger
#' only when `categories`, `values`, or factor levels declare the scale;
#' `"pairwise"` reports Conger and Fleiss only.
#'
#' @section Missing data:
#' Subjects with no ratings are dropped. Every rater must be observed, and
#' every rater pair must rate at least one subject together (every g-tuple for
#' `g > 2`), otherwise the coefficient is not identified and the call errors.
#' `"nt_fiml_cs"` needs only a connected co-observation graph.
#'
#' @section Fallback for `"nt_fiml"`:
#' The normal-theory likelihood breaks down when the fitted covariance of the
#' ratings is singular or nearly so, for example when two raters agree
#' perfectly on the subjects they share. Kappa itself is still well defined
#' there, but the maximum-likelihood fit and its sandwich standard error are
#' not. This is common in small samples with high agreement and in planned
#' designs where each rater pair shares few subjects. On complete data
#' `"nt_fiml"` equals `"pairwise"` and never fails.
#'
#' By default (`fallback = "pairwise"`), a failed or non-converged `"nt_fiml"`
#' fit is replaced by the `"pairwise"` estimate of the same coefficient, with
#' a warning and the reason recorded in `fit$fallback`. `"pairwise"` is the
#' moment estimator that `"nt_fiml"` reduces to on complete data, so the
#' combination is consistent under MCAR everywhere and keeps `"nt_fiml"`'s
#' advantage under MAR whenever the fit succeeds. In the simulation study of
#' Moss (in preparation) the combination had smaller bias at n = 40 than
#' either estimator alone, because `"nt_fiml"` fails selectively on
#' high-agreement samples. Use `fallback = "none"` to get the error instead.
#'
#' @section Generalizations:
#' `g > 2` uses g-wise kernels that compare `g` raters at once (with `"ipw"`
#' or `"cat_fiml"`). A subjects-by-raters-by-features array gives
#' vector-valued agreement; see the beyond-pairwise article.
#'
#' @param x A subjects-by-raters matrix or data frame of ratings, or a
#'   subjects-by-raters-by-features array of vector-valued ratings. `NA` marks
#'   missing entries.
#' @param estimator `"ipw"` (default), `"cat_fiml"`, `"pairwise"`,
#'   `"nt_fiml"`, or `"nt_fiml_cs"`.
#' @param weight `"nominal"` (default), `"linear"`, or `"quadratic"`. The
#'   scored estimators always use `"quadratic"`. `"hubert"` is an
#'   all-raters-equal kernel for `g > 2`. For vector-valued ratings it selects
#'   the component loss (Hamming, absolute, or squared).
#' @param g Number of raters the disagreement kernel compares at once. `2`
#'   (default) is ordinary kappa. Quadratic agreement does not depend on `g`.
#' @param values Optional numeric scores, one per category in scale order.
#'   Defaults to the categories themselves when numeric, else `1:C`.
#' @param categories Optional vector declaring the full rating scale in order,
#'   for example `1:5` or `c("low", "mid", "high")`.
#' @param em_options Named list of EM settings for the likelihood estimators:
#'   `tol` and `max_iter`, and for `"cat_fiml"` also `prune_tol`,
#'   `start_alpha`, `info_rcond`, and `flatten`. `info_rcond` is the relative
#'   eigenvalue cutoff of the information matrix used for the standard
#'   errors. `flatten > 0` adds a small Dirichlet pseudo-count that makes the
#'   fitted table unique when the likelihood is flat; the estimate barely
#'   moves but standard errors become conservative. Settings an estimator
#'   does not use are an error.
#' @param fallback For `estimator = "nt_fiml"` only: `"pairwise"` (the
#'   default) returns the `"pairwise"` estimate when the fit fails or does
#'   not converge; `"none"` raises the error. See the section on the fallback.
#' @param ... Arguments for vector-valued ratings: `W` (a feature metric for
#'   the squared loss), `ranges` (a features-by-2 matrix of score bounds, for
#'   Brennan-Prediger), `feature_weights`, and `loss = "rms"`. For `g > 2`:
#'   `max_chance_tuples` (default `5e6`). Anything else is an error.
#'
#' @return A `misskappa_estimate`: the estimates, their covariance matrix,
#'   the number of subjects `n`, the rating scale `categories`, and the
#'   per-subject influence rows `psi`, with `vcov = crossprod(psi) / n^2`.
#'   `"cat_fiml"` fits also carry two diagnostics per coefficient:
#'   `null_frac`, the share of the coefficient's gradient that lies in
#'   directions the likelihood does not identify, and `se_cutoff_ratio`, the
#'   standard error recomputed with `info_rcond` lowered 100-fold divided by
#'   the reported one (values well above 1 mean the standard error leans on
#'   weakly identified directions of a sparse table).
#'   Methods: `print`, `coef`, `vcov`, `confint`, `as.data.frame`, and
#'   `stats::influence`.
#'
#' @references
#' Moss, J. (in preparation). *Agreement coefficients with missing ratings.*
#'
#' Moss, J. (in preparation). *Quadratically weighted agreement coefficients.*
#'
#' @seealso [agreement_counts()] for counts-format data, [agreement_test()]
#'   to compare coefficients, [confint.misskappa_estimate()].
#'
#' @examples
#' # Five raters, three categories, some ratings missing.
#' agreement(dat.klein2018)
#'
#' # Ordered categories with linear weights, on a declared 1-4 scale.
#' agreement(dat.klein2018, weight = "linear", categories = 1:4)
#'
#' # Quadratic agreement on a 1-5 scale under MAR, normal theory.
#' agreement(dat.zapf2016, estimator = "nt_fiml", categories = 1:5)
#'
#' # A planned-incomplete design: some rater pairs never meet, so only the
#' # compound-symmetry estimator is identified.
#' agreement(dat.tenhove2025, estimator = "nt_fiml_cs")
#'
#' @export
agreement <- function(x,
                  estimator = c("ipw", "cat_fiml", "pairwise", "nt_fiml",
                                "nt_fiml_cs"),
                  weight = c("nominal", "linear", "quadratic", "hubert"),
                  g = 2L,
                  values = NULL,
                  categories = NULL,
                  em_options = list(),
                  fallback = NULL,
                  ...) {
  estimator <- match.arg(estimator)
  weight_supplied <- !missing(weight)
  g_supplied <- !missing(g)
  weight <- match.arg(weight)
  fallback <- .resolve_fallback(fallback, estimator)
  dots <- list(...)
  kept <- .drop_empty_subjects(x)
  x <- kept$x
  fit_one <- function(estimator, em_options) {
    .agreement_dispatch(x, estimator, weight, weight_supplied, g, g_supplied,
                        values, categories, em_options, dots, kept$ids)
  }
  if (fallback == "none") return(fit_one(estimator, em_options))
  .with_fallback(fit_one(estimator, em_options), fit_one(fallback, list()),
                 from = estimator, to = fallback)
}

# fallback = NULL means the estimator's default: "pairwise" for "nt_fiml",
# none otherwise. An explicit fallback for another estimator is an error.
.resolve_fallback <- function(fallback, estimator) {
  if (is.null(fallback)) return(if (estimator == "nt_fiml") "pairwise" else "none")
  if (!is.character(fallback) || length(fallback) != 1L ||
      !fallback %in% c("pairwise", "none")) {
    stop('\'fallback\' must be "pairwise" or "none".', call. = FALSE)
  }
  if (estimator != "nt_fiml") {
    stop('\'fallback\' applies only to estimator = "nt_fiml".', call. = FALSE)
  }
  fallback
}

# Evaluate `primary`; if its estimation fails (singular or nearly singular
# covariance, numerical failure, or no EM convergence), return `secondary`
# instead, recording why in fit$fallback and warning. Any other error (bad
# arguments, an unidentified design) is raised unchanged, as is the primary
# error when the fallback fails too. Both arguments are lazy.
.fallback_triggers <- "covariance of the ratings is singular|Numerical error|did not converge"

.with_fallback <- function(primary, secondary, from, to) {
  reason <- NULL
  fit <- tryCatch(
    withCallingHandlers(primary, warning = function(w) {
      if (grepl("did not converge", conditionMessage(w))) {
        reason <<- "its EM did not converge"
        invokeRestart("muffleWarning")
      }
    }),
    error = function(e) e)
  if (inherits(fit, "error")) {
    if (!grepl(.fallback_triggers, conditionMessage(fit))) stop(fit)
    reason <- conditionMessage(fit)
  } else if (is.null(reason)) {
    return(fit)
  }
  primary_error <- if (inherits(fit, "error")) fit else NULL
  alt <- tryCatch(secondary, error = function(e) e)
  if (inherits(alt, "error")) {
    if (!is.null(primary_error)) stop(conditionMessage(primary_error), call. = FALSE)
    stop(conditionMessage(alt), call. = FALSE)
  }
  alt$fallback <- list(from = from, reason = reason)
  warning(sprintf(paste0(
    'estimator = "%s" failed (%s); returning the "%s" estimate, which assumes ',
    'ratings are missing completely at random. Use fallback = "none" to get ',
    'the error instead.'), from, sub("[.]$", "", reason), to), call. = FALSE)
  alt
}

# The body of agreement() for one estimator, after argument matching.
.agreement_dispatch <- function(x, estimator, weight, weight_supplied, g,
                                g_supplied, values, categories, em_options,
                                dots, ids) {
  em_allowed <- switch(estimator,
    cat_fiml = c(.em_options_categorical, "flatten"),
    nt_fiml = , nt_fiml_cs = .em_options_normal,
    character(0))

  ## Vector-valued ratings: subjects-by-raters-by-features array.
  if (is.array(x) && length(dim(x)) == 3L) {
    if (g_supplied) {
      stop("'g' does not apply to vector-valued ratings.", call. = FALSE)
    }
    if (!is.null(values) || !is.null(categories)) {
      stop("'values' and 'categories' do not apply to vector-valued ratings; ",
           "use 'ranges' to declare the score range.", call. = FALSE)
    }
    em_options <- .check_em_options(em_options, em_allowed, estimator)
    fit <- .kappa_vector_dispatch(x, estimator, weight, weight_supplied,
                                  em_options, dots)
    return(.label_subjects(fit, ids))
  }

  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g < 2 ||
      g != round(g)) {
    stop("'g' must be an integer >= 2.", call. = FALSE)
  }
  g <- as.integer(g)
  if (weight == "hubert" && g == 2L) {
    stop('weight = "hubert" is only defined for g > 2; use weight = ',
         '"nominal" for g = 2.', call. = FALSE)
  }

  ## Quadratic weighting is a moment functional and therefore g-invariant: any
  ## g collapses to the same coefficient, so always take the closed form
  ## (the g = 2 path) and never the combinatorial g-wise enumeration.
  if (weight == "quadratic" || g == 2L) {
    .check_dots(dots, character(0), "scalar ratings with g = 2")
    em_options <- .check_em_options(em_options, em_allowed, estimator)
    ## Every public scalar g = 2 path surfaces a fixed-rater coefficient, so
    ## the saturated functional needs a complete rater co-observation graph.
    ## The compound-symmetry normal FIML only needs a connected one; the
    ## reduced information check in C++ remains the final guard.
    is_cs <- estimator == "nt_fiml_cs"
    .check_pattern_identifiable(
      .pattern_observed(x), unit = "rater",
      require = if (is_cs) "connected" else "complete",
      coefficient = if (is_cs) {
        "Conger's kappa under the compound-symmetry covariance model"
      } else {
        "Conger's kappa (and Brennan-Prediger)"
      },
      arity = 2L)
    fit <- .kappa_scalar_g2(x, estimator, weight, weight_supplied, values,
                            categories, em_options)
    return(.label_subjects(fit, ids))
  }

  ## g > 2, non-quadratic: the Frechet / Hubert g-wise family.
  .check_dots(dots, "max_chance_tuples", "g-wise kappa")
  if (estimator == "cat_fiml") em_allowed <- .em_options_categorical
  em_options <- .check_em_options(em_options, em_allowed, estimator)
  .check_pattern_identifiable(
    .pattern_observed(x), unit = "rater",
    coefficient = sprintf("%s-wise kappa", g),
    arity = g)
  fit <- .kappa_scalar_gwise(x, estimator, weight, g, values, categories,
                             em_options, dots)
  .label_subjects(fit, ids)
}

# Subjects with no observed rating carry no information about agreement, so
# every estimator drops them before fitting. Returns the reduced input and the
# identifiers of the kept subjects (row names, or row numbers of the input).
.drop_empty_subjects <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x) &&
      !(is.array(x) && length(dim(x)) == 3L)) {
    return(list(x = x, ids = NULL))
  }
  n <- if (is.data.frame(x)) nrow(x) else dim(x)[1L]
  ids <- if (is.data.frame(x)) {
    # Automatic data-frame row names count as none: use row numbers.
    if (.row_names_info(x) < 0L) NULL else row.names(x)
  } else {
    dimnames(x)[[1L]]
  }
  if (is.null(ids)) ids <- seq_len(n)
  ids <- as.character(ids)
  m <- if (is.data.frame(x)) as.matrix(x) else x
  observed <- !is.na(m)
  if (is.numeric(m)) observed <- observed & is.finite(m)
  keep <- if (length(dim(m)) == 3L) {
    apply(observed, 1L, any)
  } else {
    rowSums(observed) > 0L
  }
  if (!all(keep)) {
    x <- if (length(dim(x)) == 3L) x[keep, , , drop = FALSE] else
      x[keep, , drop = FALSE]
  }
  list(x = x, ids = ids[keep])
}

# Name the influence rows by subject so paired tests can align fits whose
# subject sets differ (for example when a subject is unrated on one attribute).
.label_subjects <- function(fit, ids) {
  if (!is.null(ids) && !is.null(fit$psi) && is.matrix(fit$psi) &&
      nrow(fit$psi) == length(ids)) {
    rownames(fit$psi) <- ids
  }
  fit
}

# Ordinary pairwise (g = 2) kappa, also reached by any quadratic request
# because quadratic kappa does not depend on g.
.kappa_scalar_g2 <- function(x, estimator, weight, weight_supplied, values,
                             categories, em_options) {
  if (estimator %in% c("pairwise", "nt_fiml", "nt_fiml_cs")) {
    if (weight_supplied && weight != "quadratic") {
      stop(sprintf('estimator = "%s" requires weight = "quadratic".', estimator),
           call. = FALSE)
    }
    sc <- .scored_ratings(x, categories, values)
    if (estimator == "pairwise") {
      fit <- .kappa_continuous_fit(sc$X, method = "available",
                                   weight = "quadratic")
    } else {
      em_options$covariance <- if (estimator == "nt_fiml_cs") {
        "variance_components"
      } else {
        "saturated"
      }
      fit <- .kappa_quadratic_fiml(sc$X, values = sc$values,
                                  em_options = em_options)
      fit$covariance <- NULL    # the estimator name carries the structure
      fit$mean <- NULL
      # Brennan-Prediger needs a declared scale (see ?agreement).
      if (!sc$declared) fit <- .drop_coefficients(fit, "Brennan-Prediger")
    }
    fit$method <- estimator
    fit$weight <- "quadratic"
    fit$categories <- sc$categories
    return(fit)
  }

  r <- .code_ratings(x, categories, values)
  if (weight != "nominal") {
    .require_ordered(r, sprintf('weight = "%s"', weight))
  }
  cpp_weight <- if (weight == "nominal") "identity" else weight
  fit <- .kappa_raw_fit(
    r$codes,
    method = if (estimator == "ipw") "ipw" else "fiml",
    weight = cpp_weight,
    values = r$values,
    em_options = em_options,
    categories = seq_along(r$categories)
  )
  fit$method <- estimator
  fit$weight <- weight
  fit$categories <- r$categories
  fit
}

# g > 2, non-quadratic. Maps the public weight to a g-wise kernel and the
# public estimator to the backend method, then delegates to .kappa_gwise().
.kappa_scalar_gwise <- function(x, estimator, weight, g, values, categories,
                                em_options, dots) {
  map <- switch(weight,
    nominal = list(distance = "nominal", family = "categorical"),
    hubert  = list(distance = "hubert",  family = "categorical"),
    linear  = list(distance = "absolute", family = "continuous"),
    stop(sprintf(
      'weight = "%s" has no g-wise (g > 2) kernel; use "nominal", "linear", "hubert", or "quadratic".',
      weight), call. = FALSE)
  )

  method <- switch(estimator,
    ipw      = "ipw",
    cat_fiml = if (map$family == "categorical") "fiml" else stop(
      'estimator = "cat_fiml" has no continuous g-wise variant; use ',
      'weight = "nominal"/"hubert", or estimator = "ipw".', call. = FALSE),
    stop(sprintf(paste0(
      'estimator = "%s" is a quadratic estimator; for g > 2 with weight = ',
      '"%s" use "ipw" (equal to the complete-data estimator when no ratings ',
      'are missing) or "cat_fiml".'), estimator, weight), call. = FALSE)
  )

  max_chance_tuples <- if (!is.null(dots$max_chance_tuples)) {
    dots$max_chance_tuples
  } else {
    5000000L
  }
  r <- .code_ratings(x, categories, values)
  if (map$family == "continuous") {
    .require_ordered(r, sprintf('weight = "%s"', weight))
  }
  xx <- if (map$family == "continuous") {
    matrix(r$values[r$codes], nrow(r$codes), ncol(r$codes))
  } else {
    r$codes
  }

  fit <- .kappa_gwise(xx, distance = map$distance, method = method, g = g,
                     em_options = em_options,
                     max_chance_tuples = max_chance_tuples)
  fit$method <- estimator
  fit$weight <- weight
  fit
}

# Vector-valued ratings (subjects-by-raters-by-features array).
.kappa_vector_dispatch <- function(x, estimator, weight, weight_supplied,
                                   em_options, dots) {
  W <- dots$W
  ranges <- dots$ranges
  # Full-matrix-W quadratic family: the squared vector loss (x - y)' W (x - y).
  # The normal-FIML covariance models ("nt_fiml" saturated, "nt_fiml_cs" block
  # compound symmetry) report Conger and Fleiss, plus Brennan-Prediger when
  # `ranges` declares the score range.
  if (estimator %in% c("nt_fiml", "nt_fiml_cs")) {
    .check_dots(dots, c("W", "ranges"),
                sprintf('vector-valued estimator = "%s"', estimator))
    if (weight_supplied && weight != "quadratic") {
      stop(sprintf('estimator = "%s" requires weight = "quadratic".', estimator),
           call. = FALSE)
    }
    fit <- .kappa_vector_quadratic(x, method = estimator, W = W,
                                  em_options = em_options, ranges = ranges)
    if (is.null(ranges)) fit <- .drop_coefficients(fit, "Brennan-Prediger")
    fit$method <- estimator
    fit$weight <- "quadratic"
    return(fit)
  }
  if (estimator == "pairwise" && !is.null(W)) {
    .check_dots(dots, "W",
                'vector-valued estimator = "pairwise" with a feature metric W')
    if (weight_supplied && weight != "quadratic") {
      stop('a full feature metric W requires weight = "quadratic".',
           call. = FALSE)
    }
    fit <- .kappa_vector_quadratic(x, method = "pairwise", W = W)
    fit$method <- estimator
    fit$weight <- "quadratic"
    return(fit)
  }
  if (!estimator %in% c("pairwise", "ipw")) {
    stop(sprintf(
      'vector-valued ratings support estimator = "ipw", "pairwise", "nt_fiml", or "nt_fiml_cs"; got "%s".',
      estimator), call. = FALSE)
  }
  # Component-separable (diagonal-weight) path.
  .check_dots(dots, c("feature_weights", "loss"),
              sprintf('component-separable vector estimator = "%s"', estimator))
  loss <- if (!is.null(dots$loss)) {
    dots$loss
  } else {
    switch(weight,
      nominal   = "hamming",
      linear    = "absolute",
      quadratic = "squared",
      stop(sprintf(
        'weight = "%s" has no vector component loss; use "nominal", "linear", or "quadratic", or pass loss = "rms".',
        weight), call. = FALSE)
    )
  }
  fit <- .kappa_vector(x, method = estimator, loss = loss,
                      feature_weights = dots$feature_weights)
  fit$method <- estimator
  fit$weight <- if (!is.null(dots$loss)) dots$loss else weight
  fit
}

#' Convert rater-identified ratings to counts format
#'
#' @description
#' Turns a subjects-by-raters matrix of ratings into the subjects-by-categories
#' counts used by [agreement_counts()]. Missing ratings are not counted, so a
#' subject's row total is the number of raters who rated it. This discards
#' rater identity: the raters become exchangeable.
#'
#' @param x A subjects-by-raters matrix or data frame of category labels
#'   (integer codes, characters, or factors). `NA` (and, for numeric
#'   input, any non-finite value) marks a missing rating.
#' @param categories Optional vector of all categories, in column order.
#'   Unused categories become all-zero columns. Defaults to the shared factor
#'   levels, or else the sorted observed labels. Pass it for text ratings
#'   that have an order, since the column order sets the linear and
#'   quadratic weights in [agreement_counts()].
#'
#' @return An integer subjects-by-categories matrix whose columns are named by
#'   `categories` and whose row names are carried over from `x`. Suitable as
#'   the `x` argument of [agreement_counts()].
#'
#'
#' @examples
#' counts <- ratings_to_counts(dat.gwet2014)
#' head(counts)
#' agreement_counts(counts)
#'
#' @export
ratings_to_counts <- function(x, categories = NULL) {
  if (is.data.frame(x)) {
    is_fac <- vapply(x, is.factor, logical(1))
    if (is.null(categories) && length(is_fac) && all(is_fac)) {
      lv <- lapply(x, levels)
      if (all(vapply(lv, identical, logical(1), lv[[1L]]))) categories <- lv[[1L]]
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("'x' must be a subjects-by-raters matrix or data frame.",
         call. = FALSE)
  }
  if (is.numeric(x)) x[!is.finite(x)] <- NA

  observed <- x[!is.na(x)]
  if (length(observed) == 0L) {
    stop("'x' has no observed (non-missing) ratings.", call. = FALSE)
  }

  if (is.null(categories)) {
    categories <- sort(unique(as.vector(observed)))
  } else {
    extra <- setdiff(unique(as.vector(observed)), categories)
    if (length(extra) > 0L) {
      shown <- extra[seq_len(min(5L, length(extra)))]
      stop("Observed labels absent from 'categories': ",
           paste(shown, collapse = ", "),
           if (length(extra) > 5L) ", ..." else "", call. = FALSE)
    }
  }
  C <- length(categories)

  idx <- matrix(match(x, categories), nrow = nrow(x))  # NA ratings stay NA
  counts <- matrix(0L, nrow = nrow(x), ncol = C,
                   dimnames = list(rownames(x), as.character(categories)))
  for (i in seq_len(nrow(x))) {
    row_idx <- idx[i, !is.na(idx[i, ])]
    if (length(row_idx) > 0L) counts[i, ] <- tabulate(row_idx, nbins = C)
  }
  counts
}

#' Agreement coefficients for counts-format data
#'
#' @description
#' Estimates Fleiss and Brennan-Prediger coefficients when each row is a
#' subject and each column a category, with `x[i, k]` the number of raters who
#' put subject `i` in category `k`. Rows may have different totals, which is
#' how missing ratings show up in this format. Raters are anonymous, so they
#' are treated as exchangeable and Conger is not reported. For rater-identified
#' data use [agreement()]; [ratings_to_counts()] converts it to this format.
#'
#' \itemize{
#'   \item `"pairwise"` (default, recommended): the moment estimator of
#'     Fleiss and Cuzick (1979), extended to linear and quadratic weights.
#'     It averages agreement over the rater pairs within each subject, the
#'     counts-format analogue of `"pairwise"` in [agreement()]. Consistent
#'     under MCAR, any `weight`.
#'   \item `"cat_fiml"`: full-information maximum likelihood over the
#'     category compositions, consistent under MAR. It is the
#'     `"cat_fiml"` of [agreement()] with the raters constrained to be
#'     exchangeable (the same likelihood), which is why it shares the name.
#'   \item `"nt_fiml_cs"`: normal-theory compound-symmetry FIML, quadratic
#'     weight only.
#' }
#' Rows with fewer than two ratings carry no agreement information; at least
#' one row needs two or more.
#'
#' @param x A subjects-by-categories matrix of non-negative counts.
#' @param estimator `"pairwise"` (default), `"cat_fiml"`, or
#'   `"nt_fiml_cs"`.
#' @param weight `"nominal"` (default), `"linear"`, or `"quadratic"`.
#'   `"nt_fiml_cs"` always uses `"quadratic"`.
#' @param values Optional numeric category scores for the linear and
#'   quadratic weights. Defaults to `1:C`.
#' @param r_total For `"cat_fiml"`: the number of raters each subject would
#'   have with no missing ratings. Defaults to the largest row total.
#' @param em_options Named list of EM settings for the likelihood estimators:
#'   `tol` and `max_iter`, and for `"cat_fiml"` also `prune_tol`,
#'   `start_alpha`, and `info_rcond`.
#'
#' @return A `misskappa_estimate` with `Fleiss` and `Brennan-Prediger`; see
#'   [agreement()].
#'
#' @references
#' Fleiss, J. L., & Cuzick, J. (1979). The reliability of dichotomous
#' judgments: Unequal numbers of judges per subject. *Applied Psychological
#' Measurement*, 3(4), 537-542.
#'
#' Moss, J. (in preparation). *Quadratically weighted agreement coefficients.*
#'
#' @examples
#' # Fleiss (1971): 30 subjects, 6 raters, 5 diagnostic categories.
#' agreement_counts(dat.fleiss1971)
#'
#' # Ordered categories with quadratic weights.
#' agreement_counts(dat.fleiss1971, weight = "quadratic")
#'
#' @export
agreement_counts <- function(x,
                         estimator = c("pairwise", "cat_fiml", "nt_fiml_cs"),
                         weight = c("nominal", "linear", "quadratic"),
                         values = NULL,
                         r_total = NULL,
                         em_options = list()) {
  estimator <- match.arg(estimator)
  weight_supplied <- !missing(weight)
  weight <- match.arg(weight)
  if (estimator == "nt_fiml_cs") {
    if (weight_supplied && weight != "quadratic") {
      stop('estimator = "nt_fiml_cs" requires weight = "quadratic".',
           call. = FALSE)
    }
    weight <- "quadratic"
  }
  em_options <- .check_em_options(em_options, switch(estimator,
    cat_fiml = .em_options_categorical,
    nt_fiml_cs = .em_options_normal,
    character(0)), estimator)
  cpp_weight <- if (weight == "nominal") "identity" else weight
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.", call. = FALSE)
  if (anyNA(x_mat) || any(x_mat < 0) || any(x_mat != round(x_mat))) {
    stop("'x' must hold non-negative whole-number counts.", call. = FALSE)
  }
  ids <- rownames(x_mat)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(x_mat)))
  # Subjects nobody rated carry no information; drop them as agreement() does.
  keep <- rowSums(x_mat) > 0
  if (!any(keep)) stop("'x' has no observed ratings.", call. = FALSE)
  x_mat <- x_mat[keep, , drop = FALSE]
  ids <- ids[keep]
  storage.mode(x_mat) <- "integer"

  if (estimator == "nt_fiml_cs") {
    # Direct-on-counts exchangeable CS quadratic kappa (no realizing table;
    # see kappa_counts_cs.R).
    return(.label_subjects(.kappa_counts_cs(x_mat, values, em_options), ids))
  }
  out <- if (estimator == "pairwise") {
    rcpp_kappa_counts(x = x_mat, weight_type = cpp_weight, values = values)
  } else {
    if (is.null(r_total)) r_total <- as.integer(max(rowSums(x_mat)))
    rcpp_kappa_fiml_counts(
      x = x_mat,
      weight_type = cpp_weight,
      values = values,
      r_total = as.integer(r_total),
      em_options = em_options
    )
  }
  fit <- .new_estimate(out, c("Fleiss", "Brennan-Prediger"),
                       method = estimator, weight = weight)
  fit <- .attach_counts_kappa_delete_one(fit, x_mat, cpp_weight, values,
                                        em_options)
  .label_subjects(fit, ids)
}

# Continuous-rating moment fit (plumbing). agreement(estimator = "pairwise")
# calls it with method "available" and the quadratic kernel;
# research_agreement_continuous() exposes the other methods and kernels.
.kappa_continuous_fit <- function(x_mat, method, weight) {
  out <- rcpp_kappa_continuous(x = x_mat, method = method, weight_type = weight)
  fit <- .new_estimate(out, c("Conger", "Fleiss"), method = method,
                       weight = weight)
  .attach_continuous_kappa_delete_one(fit, x_mat, weight)
}

#' Component-separable vector agreement coefficients
#'
#' @description
#' Internal pilot for pairwise vector-valued ratings with component-wise
#' missingness. Input is a rectangular subjects-by-raters-by-features array.
#' The estimator forms diagonal-weighted component-loss moments and reports
#' Conger and Fleiss coefficients.
#'
#' With IPW, each positive-weight rater-feature cell must be observed at
#' least once. A vector coefficient is undefined when no positive-weight rater
#' pairs are jointly observed within subjects.
#'
#' @param x Numeric array with dimensions subjects, raters, features.
#' @param method `"pairwise"` for observed component-pair moments or `"ipw"`
#'   for inverse-probability weighting by rater-feature and
#'   rater-pair-feature observation rates.
#' @param loss Component-separable loss: `"hamming"`, `"absolute"`,
#'   `"squared"`, or `"rms"`.
#' @param feature_weights Non-negative diagonal feature weights. Defaults to
#'   equal weights.
#'
#' @return A `misskappa_estimate` object with `Conger` and `Fleiss`
#'   coefficients and a 2x2 influence-function covariance matrix.
#'
#' @noRd
.kappa_vector <- function(x,
                         method = c("pairwise", "ipw"),
                         loss = c("hamming", "absolute", "squared", "rms"),
                         feature_weights = NULL) {
  method <- match.arg(method)
  loss <- match.arg(loss)
  dims <- dim(x)
  if (length(dims) != 3L) {
    stop("'x' must be a subjects-by-raters-by-features array.", call. = FALSE)
  }
  if (dims[1L] < 1L || dims[2L] < 2L || dims[3L] < 1L) {
    stop("'x' must have at least one subject, two raters, and one feature.",
         call. = FALSE)
  }
  if (!is.numeric(x)) stop("'x' must be numeric.", call. = FALSE)
  storage.mode(x) <- "double"

  n <- dims[1L]
  R <- dims[2L]
  p <- dims[3L]
  if (is.null(feature_weights)) {
    feature_weights <- rep(1, p)
  }
  if (!is.numeric(feature_weights) || length(feature_weights) != p ||
      any(!is.finite(feature_weights)) || any(feature_weights < 0) ||
      !any(feature_weights > 0)) {
    stop("'feature_weights' must be a non-negative numeric vector of length ",
         "dim(x)[3] with at least one positive entry.", call. = FALSE)
  }

  flat <- matrix(NA_real_, nrow = n, ncol = R * p)
  for (r in seq_len(R)) {
    for (l in seq_len(p)) {
      flat[, (r - 1L) * p + l] <- x[, r, l]
    }
  }

  out <- rcpp_kappa_vector(
    x = flat,
    features = as.integer(p),
    method = method,
    loss_type = loss,
    feature_weights = as.numeric(feature_weights)
  )

  .new_estimate(out, c("Conger", "Fleiss"),
                method = paste0("vector_", method), weight = loss,
                loss = loss, feature_weights = as.numeric(feature_weights))
}

#' G-wise agreement coefficients for rectangular ratings
#'
#' @description
#' Estimates g-wise Conger-type and Fleiss-type agreement coefficients using
#' multirater disagreement kernels (the Frechet / Hubert family). The
#' `"complete"` estimator requires every entry observed. `"ipw"` and `"fiml"`
#' admit missing entries (FIML is categorical only) when every requested rater
#' g-tuple is jointly observed by at least one subject.
#'
#' @param x A subjects-by-raters matrix or data frame.
#' @param distance Multirater disagreement kernel. `"nominal"` uses Frechet
#'   mode disagreement for categorical ratings. `"absolute"` uses median
#'   absolute deviation. `"quadratic"` uses mean squared deviation. `"hubert"`
#'   uses all-raters-equal disagreement.
#' @param method `"complete"` (complete-data), `"ipw"` (MCAR inverse-probability
#'   weighting), or `"fiml"` (categorical full-information ML under ignorable
#'   missingness, not defined for the continuous kernels).
#' @param g Arity of the multirater distance. Defaults to all raters
#'   (`ncol(x)`). Must be between 2 and `ncol(x)`.
#' @param em_options Named list of EM options used by `method = "fiml"`.
#' @param max_chance_tuples Maximum number of direct `n^g` item tuples to
#'   evaluate before stopping for continuous distances, or finite category
#'   tuples for categorical distances.
#'
#' @return A `misskappa_estimate` object with `Conger` and `Fleiss`
#'   coefficients and a 2x2 influence-function covariance matrix.
#'
#' @noRd
.kappa_gwise <- function(x,
                        distance = c("nominal", "absolute", "quadratic", "hubert"),
                        method = c("complete", "ipw", "fiml"),
                        g = NULL,
                        em_options = list(),
                        max_chance_tuples = 5000000L) {
  distance <- match.arg(distance)
  method <- match.arg(method)
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.", call. = FALSE)
  }
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat)) stop("'x' must be numeric.", call. = FALSE)
  if (is.null(g)) g <- ncol(x_mat)
  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) ||
      g < 2 || g > ncol(x_mat)) {
    stop("'g' must be an integer between 2 and ncol(x).", call. = FALSE)
  }
  if (!is.numeric(max_chance_tuples) || length(max_chance_tuples) != 1L ||
      !is.finite(max_chance_tuples) || max_chance_tuples < 1 ||
      max_chance_tuples > .Machine$integer.max) {
    stop("'max_chance_tuples' must be a positive integer.", call. = FALSE)
  }

  categorical <- distance %in% c("nominal", "hubert")
  if (method == "fiml" && !categorical) {
    stop('method = "fiml" is only available for categorical g-wise distances.',
         call. = FALSE)
  }
  if (method == "complete" && any(!is.finite(x_mat))) {
    stop("the complete-data g-wise estimator requires complete, finite ratings.",
         call. = FALSE)
  }

  if (categorical) {
    cats <- sort(unique(c(x_mat)))
    x_indexed <- matrix(match(x_mat, cats) - 1L, nrow = nrow(x_mat), ncol = ncol(x_mat))
    storage.mode(x_indexed) <- "integer"
    out <- rcpp_kappa_gwise_categorical(
      x = x_indexed,
      distance_type = distance,
      method = method,
      g = as.integer(g),
      max_chance_tuples = as.integer(max_chance_tuples),
      em_options = em_options
    )
  } else {
    storage.mode(x_mat) <- "double"
    out <- rcpp_kappa_gwise_continuous(
      x = x_mat,
      distance_type = distance,
      method = method,
      g = as.integer(g),
      max_chance_tuples = as.integer(max_chance_tuples)
    )
  }

  .new_estimate(out, c("Conger", "Fleiss"),
                method = paste0("gwise_", method), weight = distance,
                distance = distance, g = as.integer(g))
}
