#' @description
#' Package for quadratically weighted agreement coefficients. See `[quadagree]`
#' for basic documentation and browse the vignettes for more details.
#' `browseVignettes(package = "quadagree")`
"_PACKAGE"
#' Calculate limits of a confidence interval.
#'
#' @param alternative Alternative choosen.
#' @param conf_level Confidence level.
#' @keywords internal
limits <- function(alternative, conf_level) {
  half <- (1 - conf_level) / 2
  if (alternative == "two.sided") {
    return(c(half, 1 - half))
  }
  if (alternative == "greater") {
    return(c(2 * half, 1))
  }
  if (alternative == "less") {
    return(c(0, conf_level))
  }
}

#' Half-vectorize matrix.
#'
#' @param x Matrix to vectorize.
#' @keywords internal
vech <- function(x) x[row(x) >= col(x)]

#' Get the indices of the diagonal for the variances.
#' @param r Number of raters.
#' @return A vector of `TRUE` where the variances are.
#' @keywords internal
get_diag_indices <- function(r) {
  indices <- rep(0, choose(r + 1, 2))
  indices[c(1, 1 + cumsum(r:2))] <- 1
  indices
}

#' Gamma matrix
#'
#' Calculate the gamma matrix from a matrix of observations.
#' @param x A numeric matrix of observations.
#' @param sigma Covariance matrix of the data.
#' @param type One of `adf`, `normal`, `elliptical` or `unbiased`.
#' @return The sample estimate of the gamma matrix.
#' @keywords internal
gamma_est <- function(x, sigma, type = "adf") {
  if (type == "adf") {
    gamma_est_adf(x, sigma)
  } else if (type == "unbiased") {
    gamma_est_unbiased(x, sigma)
  } else {
    gamma <- gamma_est_nt(sigma)
    gamma_est_nt(sigma) * kurtosis_correction(x, type = type)
  }
}

#' Asymptotically distribution free covariance matrix.
#' @param x Data.
#' @param sigma Covariance of the data.
#' @return Estimate of the ADF covariance matrix.
#' @keywords Internal.
gamma_est_adf <- function(x, sigma) {
  i_row <- \(n) unlist(lapply(seq_len(n), seq.int, n))
  i_col <- \(n) rep.int(seq_len(n), times = rev(seq_len(n)))
  rows <- i_row(ncol(x))
  cols <- i_col(ncol(x))
  y <- t(x) - colMeans(x, na.rm = TRUE)
  z <- y[cols, ] * y[rows, ]
  mat <- z - rowMeans(z, na.rm = TRUE)

  if (!anyNA(mat)) {
    div <- nrow(x)
  } else {
    nas <- is.na(mat)
    mat[nas] <- 0
    div <- tcrossprod(!nas)
  }
  tcrossprod(mat) / div
}

#' Normal theory gamma matrix
#'
#' Code obtained from `lavaan`:
#' https://github.com/yrosseel/lavaan/blob/6f047c800206d23f246d484b9522295257614222/R/lav_matrix.R
#'
#' Calculate the gamma matrix from a matrix of observations.
#' @param sigma Covariance matrix of the data.
#' @return Normal theory gamma matrix.
#' @keywords internal
gamma_est_nt <- function(sigma) {
  n <- ncol(sigma)

  lower <- lower_vec_indices(n)
  upper <- upper_vec_indices(n)

  y <- sigma %x% sigma
  out <- (y[lower, , drop = FALSE] + y[upper, , drop = FALSE]) / 2
  out[, lower, drop = FALSE] + out[, upper, drop = FALSE]
}

#' Unbiased asymptotic covariance matrix.
#'
#' @param x Data.
#' @param sigma Covariance matrix of the data.
#' @return Unbiased asymptotic covariance matrix.
#' @keywords internal
gamma_est_unbiased <- function(x, sigma) {
  gamma_adf <- gamma_est_adf(x, sigma)
  gamma_nt <- gamma_est_nt(sigma)
  gamma_rem <- tcrossprod(vech(sigma))
  n <- nrow(x)
  mult <- n / ((n - 2) * (n - 3))
  mult * ((n - 1) * gamma_adf - (gamma_nt - 2 / (n - 1) * gamma_rem))
}


#' Obtain indices of lower or upper triangular matrix in vec indices.
#'
#' Code obtained from `lavaan`:
#' https://github.com/yrosseel/lavaan/blob/6f047c800206d23f246d484b9522295257614222/R/lav_matrix.R
#'
#' @param n Dimension of square matrix.
#' @param diagonal If `TRUE`, includes the diagonal elements.
#' @returns Indices `x` so that `a[x] = c(a)[x]` returns the elements
#'    of the lower (upper) diagonal matrix in row-wise (column-wise)
#'    order.
#' @keywords internal
#'
lower_vec_indices <- function(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  if (diagonal) which(rows >= cols) else which(rows > cols)
}

upper_vec_indices <- function(n = 1L, diagonal = TRUE) {
  rows <- matrix(seq_len(n), n, n)
  cols <- matrix(seq_len(n), n, n, byrow = TRUE)
  tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
  if (diagonal) tmp[rows >= cols] else tmp[rows > cols]
}

#' Calculate unbiased sample kurtosis.
#' @param x Matrix of valus.
#' @return Unbiased sample kurtosis.
#' @keywords internal
kurtosis <- function(x) {
  n <- nrow(x)
  g2 <- \(x) mean((x - mean(x, na.rm = TRUE))^4) / stats::var(x, na.rm = TRUE)^2
  kurtosis <- \(x) (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2(x) + 6)
  mean(apply(x, 2, kurtosis), na.rm = TRUE) - 3
}

#' Calculate kurtosis correction
#' @param x Matrix of values
#' @param type The type of correction, either "normal" or "elliptical".
#' @keywords internal
kurtosis_correction <- function(x, type) {
  kurt <- if (type == "normal") 0 else kurtosis(x)
  1 + kurt / 3
}

#' Calculates the empirical capital pi matrix.
#'
#' @param x Vector of probabilities for being missing.
#' @return The capital pi matrix.
#' @keywords internal
pi_mat_empirical <- \(x) {
  r <- ncol(x)
  if (!anyNA(x)) {
    return(matrix(1, choose(r + 1, 2), choose(r + 1, 2)))
  }

  ind2 <- arrangements::combinations(r, 2, replace = TRUE)
  ind4 <- arrangements::combinations(seq_len(nrow(ind2)), 2, replace = TRUE)

  nisna <- !is.na(x)
  combs <- apply(ind2, 1, \(i) nisna[, i[1]] & nisna[, i[2]])
  p2_hats <- colMeans(combs)

  hats <- apply(ind4, 1, \(i) {
    mean(combs[, i[1]] & combs[, i[2]]) / (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r + 1, 2), choose(r + 1, 2))
  new_mat[ind4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat, uplo = "U"))
}

#' Get the required value of c1.
#' @param values Vector of values.
#' @param kind The kind of c1 requested.
#' @return The value of c1.
#' @keywords internal
bp_get_c1 <- \(values, kind) {
  if (kind == 1) {
    n_cat <- length(values)
    (2 * n_cat * sum(values^2) - 2 * sum(values)^2) / n_cat^2
  } else {
    0.5 * (max(values) - min(values))^2
  }
}

#' @export
print.quadagree <- function(x, digits = getOption("digits"), ...) {
  at <- \(y) attr(x, y)
  cat("Call: ", paste(deparse(at("call")),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  if (!is.null(x)) {
    cat(format(100 * at("conf_level")),
      "% confidence interval (n = ", at("n"), ").\n",
      sep = ""
    )
    print(x[1:2], digits = digits)
    cat("\n")
  }

  if (!is.null(at("estimate"))) {
    cat("Sample estimates.\n")
    print(
      c(
        kappa = at("estimate"),
        sd = at("sd")
      ),
      digits = digits
    )
  }
  invisible(x)
}

get_transformer <- function(transform) {
  transformers <- list(
    fisher = transformer_fisher,
    none = transformer_none,
    arcsin = transformer_arcsin,
    log = transformer_log
  )

  if (transform %in% names(transformers)) {
    transformers[[transform]]
  } else {
    stop(paste0("`transformer = ", transform, "` not supported."))
  }
}

transformer_fisher <- c(
  est = \(est) atanh(est),
  sd = \(est, sd) sd / (1 - est^2),
  inv = tanh
)

transformer_log <- c(
  est = \(est) log(1 - est),
  sd = \(est, sd) sd / abs(1 - est),
  inv = \(x) 1 - exp(x)
)

transformer_none <- c(
  est = \(est) est,
  sd = \(est, sd) sd,
  inv = \(x) x
)

transformer_arcsin <- c(
  est = asin,
  sd = \(est, sd) sd / sqrt(1 - est^2),
  inv = sin
)

tr <- \(x) {
  sum(x[1L + 0L:(dim(x)[1L] - 1L) * (dim(x)[1L] + 1L)])
}
avar <- \(x, sigma, mu, type, fleiss, pi) {
  p <- colSums(!is.na(x)) / nrow(x)
  gamma <- gamma_est(x, sigma, type)
  mat <- cov_mat_kappa(p, mu, sigma, gamma, x, pi)
  avar_(mat, mu, sigma, fleiss)
}

avar_ <- function(mat, mu, sigma, fleiss) {
  r <- nrow(sigma)
  yy <- sum(diag(sigma))
  xx <- sum(sigma) - yy
  zz <- mean(mu^2) - mean(mu)^2

  vec <- if (fleiss) {
    mult <- 1 / ((r - 1) * (r * zz + yy)^2)
    c(r * zz + yy, r * zz - xx, -r * (xx + yy))
  } else {
    mult <- 1 / (r^2 * zz + (r - 1) * yy)^2
    c(r^2 * zz + (r - 1) * yy, -(r - 1) * xx, -r^2 * xx)
  }

  c(t(vec) %*% mat %*% vec) * mult^2
}

cov_mat_kappa <- function(p, mu, sigma, gamma, x, pi) {
  r <- length(p)

  # The covariances involving s.
  gamma_pi <- gamma * pi
  d <- get_diag_indices(r)
  ones <- rep(1, length(d))
  ones_minus_d <- ones - d
  cov_ss_ss <- t(ones_minus_d) %*% (gamma_pi) %*% ones_minus_d
  cov_ss_tr <- t(ones_minus_d) %*% (gamma_pi) %*% d
  cov_tr_tr <- t(d) %*% gamma_pi %*% d

  # The variance of mean(mu^2) - mean(mu)^2
  p_mat <- diag(1 / p)
  mu_middle <- (sigma + (p_mat - diag(r)) * diag(sigma))
  mu_vec <- (diag(r) - matrix(1 / r, r, r)) %*% mu
  cov_mu_mu <- c(4 / r^2 * t(mu_vec) %*% mu_middle %*% mu_vec)
  matrix(c(
    4 * cov_ss_ss, 2 * cov_ss_tr, 0,
    2 * cov_ss_tr, cov_tr_tr, 0,
    0, 0, cov_mu_mu
  ), nrow = 3, byrow = FALSE)
}

avar_bp <- \(x, type, c1, pi) {
  p <- colSums(!is.na(x)) / nrow(x)
  mu <- colMeans(x, na.rm = TRUE)
  sigma <- stats::cov(x, use = "pairwise.complete.obs")
  gamma <- gamma_est(x, sigma, type)
  var <- cov_mat_bp(p, mu, sigma, gamma, x, pi)
  avar_bp_(var, mu, sigma, c1)
}

cov_mat_bp <- function(p, mu, sigma, gamma, x, pi) {
  r <- length(p)

  # The variance of mean(mu^2) - mean(mu)^2
  p_mat <- diag(1 / p)
  mu_middle <- (sigma + (p_mat - diag(r)) * diag(sigma))
  mu_vec <- (diag(r) - matrix(1, r, r) / r) %*% mu
  var_a <- c(4 * t(mu_vec) %*% mu_middle %*% mu_vec)

  # The covariances involving Sigma.
  gamma_pi <- gamma * pi
  d <- get_diag_indices(r)
  ones <- rep(1, length(d))
  d_minus_ones <- d * (1 + 1 / r) - 2 / r * ones
  var_b <- c(t(d_minus_ones) %*% gamma_pi %*% d_minus_ones)
  var_a + var_b
}

avar_bp_ <- function(vars, mu, sigma, c1) {
  r <- nrow(sigma)
  4 / (r - 1)^2 / c1^2 * vars
}
fleiss_prepare <- \(x, type) {
  x <- as.matrix(x)
  pi <- pi_mat_empirical(x)
  list(xx = as.matrix(x), type = type, n = nrow(x), pi = pi)
}

fleiss_fun <- \(calc) {
  n <- calc$n
  sigma <- stats::cov(calc$xx, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(calc$xx, na.rm = TRUE)
  est <- fleiss_pop(mu, sigma)
  var <- avar(calc$xx, sigma, mu, calc$type, TRUE, calc$pi)
  list(est = est, var = var)
}

conger_prepare <- \(x, type) {
  x <- as.matrix(x)
  pi <- pi_mat_empirical(x)
  list(xx = x, type = type, n = nrow(x), pi = pi)
}

conger_fun <- \(calc) {
  n <- calc$n
  sigma <- stats::cov(calc$xx, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(calc$xx, na.rm = TRUE)
  est <- conger_pop(mu, sigma)
  var <- avar(calc$xx, sigma, mu, calc$type, FALSE, calc$pi)
  list(est = est, var = var)
}

bp_prepare <- \(x, values, kind, type) {
  x <- as.matrix(x)
  pi <- pi_mat_empirical(x)
  if (is.null(values)) values <- stats::na.omit(unique(c(x)))
  c1 <- bp_get_c1(values, kind)
  list(xx = x, c1 = c1, type = type, n = nrow(x), pi = pi)
}

bp_fun <- \(calc) {
  est <- bp_calc(calc$xx, calc$c1)
  var <- avar_bp(calc$xx, calc$type, calc$c1, calc$pi)
  list(est = est, var = var)
}
#' @keywords internal
quadagree_internal <- function(calc,
                               transform,
                               conf_level,
                               alternative = c("two.sided", "greater", "less"),
                               bootstrap,
                               n_reps,
                               call,
                               fun,
                               ...) {
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)
  quants <- limits(alternative, conf_level)
  est_var <- fun(calc)
  est <- est_var$est
  sd <- sqrt(max(est_var$var, 0))

  out <- if (!bootstrap) {
    ci_asymptotic(est, sd, calc$n, transformer, quants)
  } else {
    ci_boot(calc, fun, transformer, quants, n_reps)
  }

  out[2] <- min(out[2], 1)
  names(out) <- quants
  structure(out,
    conf_level = conf_level,
    alternative = alternative,
    type = calc$type,
    n = calc$n,
    transform = transform,
    bootstrap = bootstrap,
    n_reps = n_reps,
    estimate = est,
    sd = sd,
    call = call,
    class = "quadagree"
  )
}

#' @keywords internal
ci_asymptotic <- function(est, sd, n, transformer, quants) {
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  multiplier <- stats::qt(quants, n - 1) / sqrt(n - 1)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
ci_boot <- function(calc, fun, transformer, quants, n_reps) {
  est_var <- fun(calc)
  est_t <- transformer$est(est_var$est)
  sd_t <- transformer$sd(est_var$est, sqrt(est_var$var))
  boots <- bootstrapper(calc, est_t, fun, transformer, n_reps)
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
bootstrapper <- function(calc, est_t, fun, transformer, n_reps) {
  calc_new <- calc
  trans_est <- transformer$est
  trans_sd <- transformer$sd
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample.int(calc$n, replace = TRUE)
      calc_new$xx <- calc$xx[indices_star, ]
      est_var <- fun(calc_new)
      est_star <- est_var$est
      sd_star <- sqrt(est_var$var)
      (trans_est(est_star) - est_t) / trans_sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}
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
#' Agreement study from Zapf et. al (2016)
#'
#' Agreement study (n = 200) from Zapf et al. (2016) in wide format. There are
#'   `50` items `4` judges, and ratings from `1` to `5`. It is
#'    the case that the same set four judges rated every item, hence this
#'    data is suitable for Cohen's kappa.
#'
#' @usage dat.zapf2016
#'
#' @format A `n` times `R` matrix. There are `n = 50` row corresponding to the
#'   different items. Each of the `R = 4` columns contains the ratings of
#'   the `j`th judge.
#'
#' @keywords datasets
#'
#' @references
#' Zapf, A., Castell, S., Morawietz, L. et al. Measuring inter-rater
#' reliability for nominal data <U+2013> which coefficients and confidence
#' intervals are appropriate?. BMC Med Res Methodol 16, 93 (2016).
#' https://doi.org/10.1186/s12874-016-0200-9
#'
#' @source
#' <https://biomedcentral.com/articles/10.1186/s12874-016-0200-9#Sec14>
"dat.zapf2016"
#' Agreement study from Klein (2018)
#'
#' Agreement study (n = 10) from Klein (2018). This data set contains
#'    missing ratings.
#'
#' @usage dat.zapf2016
#'
#' @format A `n` times `R` matrix. There are `n = 20` row corresponding to the
#'   different items. Each of the `R = 5` columns contains the ratings of
#'   the `j`th judge.
#'
#' @keywords datasets
#'
#' @source
#' Klein, D. (2018). Implementing a General Framework for Assessing Interrater
#' Agreement in Stata. The Stata Journal, 18(4), 871-901.
#' https://doi.org/10.1177/1536867X1801800408
"dat.klein2018"
#' Agreement study from Gwet (2014)
#'
#' Agreement study (n = 20) from Gwet (2014), p. 125. This data set contains
#'    missing ratings.
#'
#' @usage dat.zapf2016
#'
#' @format A `n` times `R` matrix. There are `n = 20` row corresponding to the
#'   different items. Each of the `R = 5` columns contains the ratings of
#'   the `j`th judge.
#'
#' @keywords datasets
#'
#' @source
#' Gwet, K. L. (2014). Handbook of Inter-Rater Reliability.
#' Advanced Analytics, LLC.
"dat.gwet2014"
#' Agreement study from Fleiss, J. L. (1971)
#'
#' Fleiss, J. L. (1971) is best known for introducing Fleiss' kappa. In Table 1
#'    he presents diagnosis data of psychiatric patients. The data is on Fleiss
#'    form, with `n = 30` patients diagnosed by `6` psychiatrists each. It is
#'    not the case that the same set of six psychiatrists diagnosed every
#'    patient.
#'
#' @usage dat.fleiss1971
#'
#' @format The tibble contains the five columns "depression",
#'     "personality disorder", "schizophrenia" "neurosis", and "other".
#'     The content of the `ij`th cell is the number of raters who diagnosed
#'     the `i`th patient as having mental illness `j`.
#'
#' @keywords datasets
#'
#' @references
#' Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
#' Psychological Bulletin, 76(5), 378<U+2013>382.
#' https://doi.org/10.1037/h0031619
#'
#' @source
#' The data was scraped by hand from Fleiss (1971).
"dat.fleiss1971"
bp_calc <- function(x, c1) {
  y <- as.matrix(x)
  n <- nrow(y)
  sigma <- stats::cov(y, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(y, na.rm = TRUE)
  bp_pop_c1(mu, sigma, c1)
}

bp_pop_c1 <- function(mu, sigma, c1) {
  r <- ncol(sigma)
  trace <- tr(sigma)
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r
  cov_sum <- sum(sigma) / r
  d <- 2 / (r - 1) * (trace + mean_diff - cov_sum)
  1 - d / c1
}

bp_pop <- function(mu, sigma, values, type = 1) {
  c1 <- bp_get_c1(values, type)
  bp_pop_c1(mu, sigma, c1)
}

conger_pop <- function(mu, sigma, corr = 0) {
  r <- ncol(sigma)
  trace <- tr(sigma)
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r^2 + corr
  top <- sum(sigma) - trace
  bottom <- (r - 1) * trace + mean_diff
  top / bottom
}

cohen_pop <- function(mu, sigma) conger_pop(mu, sigma)

fleiss_pop <- function(mu, sigma, corr = 0) {
  n <- length(mu)
  r <- ncol(sigma)
  trace <- tr(sigma)
  mean_diff <- (sum(mu^2) / n - (sum(mu) / n)^2 + corr) * r
  top <- sum(sigma) - trace - mean_diff
  bottom <- (r - 1) * (trace + mean_diff)
  top / bottom
}
bp_aggr_prepare <- \(x, values, kind) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  c1 <- bp_get_c1(values, kind)
  xtx <- c(tcrossprod(values^2, y))
  xt12 <- c(tcrossprod(values, y))^2
  list(
    xx = cbind(xtx, xt12),
    n = nrow(y),
    r = r,
    c1 = c1
  )
}

bp_aggr_fun <- \(calc) {
  means <- colMeans(calc$xx)
  theta <- tcrossprod(t(calc$xx) - means) / calc$n

  calcr_inv_1 <- 1 / (calc$r - 1)
  calcr_inv <- 1 / calc$r
  calcc1_inv <- 1 / calc$c1
  k <- theta[1, 1] - 2 * theta[1, 2] * calcr_inv + theta[2, 2] * calcr_inv^2
  disagreement <- 2 * calcr_inv_1 * (means[1] - calcr_inv * means[2])

  est <- unname(1 - disagreement * calcc1_inv)
  var <- calcc1_inv^2 * 4 * calcr_inv_1^2 * k
  list(est = est, var = max(var, 0))
}

fleiss_aggr_prepare <- \(x, values) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  xtx <- c(tcrossprod(values^2, y))
  xt1 <- c(tcrossprod(values, y))
  xt12 <- xt1^2

  list(
    xx = cbind(xt1, xt12, xtx),
    n = nrow(y),
    r = r
  )
}

fleiss_aggr_fun <- \(calc) {
  means <- colMeans(calc$xx)
  theta <- tcrossprod(t(calc$xx) - means) / calc$n
  k <- calc$r / (means[3] * calc$r - means[1]^2)
  km <- k * (means[2] - means[1]^2)
  calc_r_inv <- 1 / (calc$r - 1)

  grad <- k * calc_r_inv * c(2 * means[1] * (km / calc$r - 1), 1, -km)
  est <- calc_r_inv * (km - 1)
  var <- c(crossprod(grad, theta %*% grad))

  list(est = est, var = var)
}
