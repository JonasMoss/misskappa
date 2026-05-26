#' Estimation and inference for the quadratically weighted Fleiss' kappa and
#'    Conger's (Cohen's) kappa
#'
#' Confidence intervals for quadratic agreement coefficients with optional
#'    bootstrapping and transforms. Based on the formulas of Moss and van Oest
#'    (wip) and Moss (wip) along with standard asymptotic theory
#'    (Magnus, Neudecker, 2019) and the missing data theory of
#'    van Praag et al. (1985).
#'
#' There are two kinds of functions. The functions ending in `aggr` should be
#'    applied to data on aggregated form, where each row contains the number
#'    of selected ratings for each category. The data set `dat.fleiss1971`
#'    provides an example. Missing data is not supported for the `aggr`
#'    functions. The other functions should be applied to data on long form,
#'    where  each row contains the ratings of every rater. The data sets
#'    `dat.zapf2016` and `dat.klein2018` are examples. Missing data, and
#'    continuous data, is supported. See the usage vignette for
#'    more information.
#'
#' For data on long form, the methods handle missing data using pairwise
#'    available information, i.e., the option `use = "pairwise.complete.obs"`
#'    in [stats::cov()] along with the asymptotic theory of
#'    van Praag et al. (1985). The bootstrap option
#'    uses the studentized bootstrap (Efron, B. 1987), which is second order
#'    correct. Both functions makes use of [`future.apply`] when bootstrapping.
#'
#' The `type` variables defaults to `adf`, asymptotically distribution-free,
#'    which is consistent when the fourth moment is finite. The `normal` option
#'     assumes normality, and is not consistent for models with excess
#'    kurtosis unequal to `0`. The `elliptical` option assumes an
#'    elliptical or pseudo-elliptical distribution of the data. The resulting
#'    confidence intervals are corrected variants of the normal theory
#'    intervals with a kurtosis correction (Yuan & Bentler 2002). The
#'    common kurtosis parameter is calculated using the unbiased sample
#'    kurtosis (Joanes, 1998).
#'
#' Conger's (1980) kappa is a multi-rater generalization of Cohen's kappa.
#'    All functions in this package work for multiple raters, so functions
#'    starting with `cohen` or `conger` are aliases. The quadratically
#'    weighted Cohen's kappa is also known as Lin's concordance coefficient.
#'
#' The only difference between Cohen's kappa and Fleiss' kappa lies on how they
#'    measure disagreement due to chance. Here Fleiss' marginalizes the rating
#'    distribution across raters, essentially assuming there is no difference in
#'    the rating distribution across raters, while Cohen's kappa does not.
#'    There is a large literature comparing Fleiss' kappa to Cohen's kappa, and
#'    there is no consensus on which to prefer.
#'
#' The aggregated functions takes an argument `values`, which specifies what
#'    numerical value to attach to each category. The default value for `values`
#'    is `1...C`, where `C` is the number of categories.
#'
#' The Brennan-Prediger coefficients take an argument `kind`. If equal to `1`,
#'    it returns the traditional Brennan-Prediger coefficient. If `kind`
#'    equals `2`, it returns the new Brennan-Prediger coefficient of Moss (wip).
#'
#' @references
#'
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. Journal of the
#' American Statistical Association, 82(397), 171-185.
#' https://doi.org/10.2307/2289144
#'
#' Van Praag, B. M. S., Dijkstra, T. K., & Van Velzen, J. (1985).
#' Least-squares theory based on general distributional assumptions with
#' an application to the incomplete observations problem.
#' Psychometrika, 50(1), 25-36. https://doi.org/10.1007/BF02294145
#'
#' Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample skewness
#' and kurtosis. Journal of the Royal Statistical Society: Series D
#' (The Statistician), 47(1), 183-189. https://doi.org/10.1111/1467-9884.00122
#'
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement with provision for
#' scaled disagreement or partial credit. Psychological Bulletin, 70(4),
#' 213-220. https://doi.org/10.1037/h0026256
#'
#' Fleiss, J. L. (1975). Measuring agreement between two judges on the presence
#' or absence of a trait. Biometrics, 31(3), 651-659.
#' https://www.ncbi.nlm.nih.gov/pubmed/1174623
#'
#' Conger, A. J. (1980). Integration and generalization of kappas for multiple
#' raters. Psychological Bulletin, 88(2), 322-328.
#' https://doi.org/10.1037/0033-2909.88.2.322
#'
#' Lin, L. I. (1989). A concordance correlation coefficient to evaluate
#' reproducibility. Biometrics, 45(1), 255-268.
#' https://www.ncbi.nlm.nih.gov/pubmed/2720055
#'
#' Moss, van Oest (work in progress). Inference for quadratically weighted
#' multi-rater kappas with missing raters.
#'
#' Moss (work in progress). On the Brennan-Prediger coefficients.
#'
#' Magnus, J. R., & Neudecker, H. (2019). Matrix Differential Calculus with
#' Applications in Statistics and Econometrics. John Wiley & Sons.
#' https://doi.org/10.1002/9781119541219
#'
#' @export
#' @param x Input data data can be converted to a matrix using `as.matrix`.
#' @param values to attach to each column on the Fleiss form data.
#'    Defaults to `1:C`, where `C` is the number of categories. Only used
#'    in `fleiss_aggr` and `bp_aggr`.
#' @param type Type of confidence interval. Either `adf`, `elliptical`, or
#'   `normal`. Ignored in `fleiss_aggrci`.
#' @param kind The kind of Brennan-Prediger coefficient used, `1` for the
#'   classical kind and `2` for the kind introduced in Moss (2023). Only
#'   relevant for `bp_aggr` and `bp`.
#' @param transform One of `"none"`, `"log"`, `"fisher"`, and `"arcsin`.
#'   Defaults to `"none"`.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of `"two.sided"` (default), `"greater"` or `"less"`.
#' @param conf_level Confidence level. Defaults to `0.95`.
#' @param bootstrap If `TRUE`, performs a studentized bootstrap with `n_reps`
#'   repetitions. Defaults to `FALSE`.
#' @param n_reps Number of bootstrap samples if `bootstrap = TRUE`. Ignored if
#'   `bootstrap = FALSE`. Defaults to `1000`.
#' @return A vector of class `quadagree` containing the confidence end points.
#'   The arguments of the function call are included as attributes.
#' @name quadagree
#' @examples
#' library("quadagree")
#' # Fleiss' kappa for data on long form
#' fleissci(dat.zapf2016)
#'
#' # Brennan-Prediger for data on aggregated form
#' bpci_aggr(dat.fleiss1971)
#'
#' # Conger's (Cohen's) kappa for data on long form with missing values
#' congerci(dat.klein2018)
fleissci <- function(x,
                     type = c("adf", "elliptical", "normal", "unbiased"),
                     transform = "none",
                     conf_level = 0.95,
                     alternative = c("two.sided", "greater", "less"),
                     bootstrap = FALSE,
                     n_reps = 1000) {
  call <- match.call()
  type <- match.arg(type)
  fun <- fleiss_fun
  calc <- fleiss_prepare(x, type)

  args <- c(
    calc = list(calc),
    utils::tail(sapply(names(formals()), str2lang), -1),
    call = quote(call),
    fun = fun
  )
  do.call(what = quadagree_internal, args = args)
}

#' @export
#' @rdname quadagree
congerci <- function(x,
                     type = c("adf", "elliptical", "normal", "unbiased"),
                     transform = "none",
                     conf_level = 0.95,
                     alternative = c("two.sided", "greater", "less"),
                     bootstrap = FALSE,
                     n_reps = 1000) {
  call <- match.call()
  type <- match.arg(type)

  fun <- conger_fun
  calc <- conger_prepare(x, type)

  args <- c(
    calc = list(calc),
    utils::tail(sapply(names(formals()), str2lang), -1),
    call = quote(call),
    fun = fun
  )
  do.call(what = quadagree_internal, args = args)
}

#' @export
#' @rdname quadagree
cohenci <- congerci

#' @export
#' @rdname quadagree
bpci <- function(x,
                 values = NULL,
                 kind = 1,
                 type = c("adf", "elliptical", "normal", "unbiased"),
                 transform = "none",
                 conf_level = 0.95,
                 alternative = c("two.sided", "greater", "less"),
                 bootstrap = FALSE,
                 n_reps = 1000) {
  stopifnot(kind == 1 || kind == 2)
  call <- match.call()
  type <- match.arg(type)
  calc <- bp_prepare(x, values, kind, type)
  fun <- bp_fun

  args <- c(
    calc = list(calc),
    utils::tail(sapply(names(formals()), str2lang), -3),
    call = quote(call),
    fun = fun
  )
  do.call(what = quadagree_internal, args = args)
}

#' @export
#' @rdname quadagree
fleissci_aggr <- function(x,
                          values = seq_len(ncol(x)),
                          transform = "none",
                          conf_level = 0.95,
                          alternative = c("two.sided", "greater", "less"),
                          bootstrap = FALSE,
                          n_reps = 1000) {
  stopifnot(ncol(x) == length(values))
  call <- match.call()
  calc <- fleiss_aggr_prepare(x, values)
  fun <- fleiss_aggr_fun

  args <- c(
    calc = list(calc),
    utils::tail(sapply(names(formals()), str2lang), -1),
    call = quote(call),
    fun = fun
  )

  do.call(what = quadagree_internal, args = args)
}

#' @export
#' @rdname quadagree
bpci_aggr <- function(x,
                      values = seq_len(ncol(x)),
                      kind = 1,
                      transform = "none",
                      conf_level = 0.95,
                      alternative = c("two.sided", "greater", "less"),
                      bootstrap = FALSE,
                      n_reps = 1000) {
  stopifnot(kind == 1 | kind == 2)
  stopifnot(ncol(x) == length(values))
  call <- match.call()
  calc <- bp_aggr_prepare(x, values, kind)
  fun <- bp_aggr_fun

  args <- c(
    calc = list(calc),
    utils::tail(sapply(names(formals()), str2lang), -1),
    call = quote(call),
    fun = fun
  )
  do.call(what = quadagree_internal, args = args)
}

#' @export
#' @rdname quadagree
fleiss <- function(x, variant = c("normal", "unbiased", "unbiased2")) {
  variant <- match.arg(variant)
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs")
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")

  if (variant == "normal") {
    sigma <- sigma * (n - 1) / n
  }

  corr <- if (variant == "unbiased2") {
    1 / (n * r^2) * sum(sigma) - 1 / (r * n) * tr(sigma)
  } else {
    0
  }
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(x, na.rm = TRUE)
  fleiss_pop(mu, sigma, corr)
}

#' @export
#' @rdname quadagree
bp <- function(x, values = stats::na.omit(unique(c(x))), kind = 1) {
  x <- as.matrix(x)
  c1 <- bp_get_c1(values, kind)
  bp_calc(x, c1)
}

#' @export
#' @rdname quadagree
conger <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(x, na.rm = TRUE)
  conger_pop(mu, sigma)
}

#' @export
#' @rdname quadagree
cohen <- function(x) conger(x)

#' @export
#' @rdname quadagree
fleiss_aggr <- \(x, values = seq_len(ncol(x))) {
  r <- sum(x[1, ])
  stopifnot(ncol(x) == length(values))

  y <- as.matrix(x)
  xtx <- tcrossprod(values^2, y)
  xt1 <- tcrossprod(values, y)

  extx <- mean(xtx)
  ext1 <- mean(xt1)
  ext2 <- mean(xt1^2)

  1 / (r - 1) * ((ext2 - ext1^2) / (extx - ext1^2 / r) - 1)
}

#' @export
#' @rdname quadagree
bp_aggr <- function(x, values = seq_len(ncol(x)), kind = 1) {
  stopifnot(kind == 1 | kind == 2)
  stopifnot(ncol(x) == length(values))
  calc <- bp_aggr_prepare(x, values, kind)
  bp_aggr_fun(calc)$est
}
