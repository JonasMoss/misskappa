# Estimators of coefficient alpha under item-level missingness. Each returns a
# uniform record: estimate, se, lwr, upr (95% Wald unless the estimator carries
# its own interval), plus flags. Failures degrade to NA rather than erroring,
# so a single bad cell never aborts the overnight run.

Z975 <- stats::qnorm(0.975)

.wald <- function(est, se) {
  if (is.finite(est) && is.finite(se)) c(est - Z975 * se, est + Z975 * se) else c(NA_real_, NA_real_)
}

# Errors -> NA (e.g. a singular information matrix when alpha is unidentified);
# warnings (e.g. "did not fully converge") are NOT fatal -- suppressing rather
# than catching them keeps valid-but-warned estimates in the hard cells.
.fit_continuous <- function(X, se_type) {
  fit <- tryCatch(suppressWarnings(misskappa::alpha_continuous(X, se_type = se_type)),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(estimate = NA_real_, se = NA_real_))
  est <- as.numeric(stats::coef(fit))
  se <- tryCatch(sqrt(as.numeric(stats::vcov(fit))), error = function(e) NA_real_)
  list(estimate = est, se = se)
}

.fit_cat_fiml <- function(X) {
  Xi <- matrix(as.integer(round(X)), nrow(X), ncol(X))
  fit <- tryCatch(suppressWarnings(misskappa::alpha(Xi, method = "fiml")),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(estimate = NA_real_, se = NA_real_))
  list(estimate = as.numeric(stats::coef(fit)),
       se = tryCatch(sqrt(as.numeric(stats::vcov(fit))), error = function(e) NA_real_))
}

# Case (row) bootstrap SE -- the oracle for "is the analytic SE right",
# independent of coverage. `point_fun` returns a scalar alpha given a matrix.
bootstrap_se <- function(X, point_fun, B = 200L) {
  n <- nrow(X)
  bs <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(point_fun(X[idx, , drop = FALSE]), error = function(e) NA_real_)
  }, numeric(1))
  stats::sd(bs[is.finite(bs)])
}

# Dispatch one estimator by name. Returns a one-row data.frame.
run_estimator <- function(name, X, do_boot = FALSE, boot_B = 200L) {
  t0 <- proc.time()[["elapsed"]]
  undefined <- NA; npd <- NA; se_boot <- NA_real_

  if (name == "pairwise") {
    r <- pairwise_alpha(X)
    est <- r$estimate; se <- r$se; undefined <- r$undefined; npd <- r$npd
    if (do_boot && !isTRUE(undefined)) {
      se_boot <- bootstrap_se(X, function(M) alpha_point_from_cov(
        stats::cov(M, use = "pairwise.complete.obs")), boot_B)
    }
    ci <- .wald(est, se)

  } else if (name == "fiml_normal") {
    r <- .fit_continuous(X, "normal"); est <- r$estimate; se <- r$se; ci <- .wald(est, se)

  } else if (name == "fiml_sandwich") {
    r <- .fit_continuous(X, "sandwich"); est <- r$estimate; se <- r$se; ci <- .wald(est, se)

  } else if (name == "cat_fiml") {
    r <- .fit_cat_fiml(X); est <- r$estimate; se <- r$se; ci <- .wald(est, se)

  } else if (name == "feldt") {
    r <- feldt_alpha(X); est <- r$estimate; se <- r$se; ci <- c(r$lwr, r$upr)

  } else {
    stop("unknown estimator: ", name)
  }

  data.frame(
    estimator = name, estimate = est, se = se, se_boot = se_boot,
    lwr = ci[1], upr = ci[2],
    undefined = undefined, npd = npd,
    time = proc.time()[["elapsed"]] - t0,
    stringsAsFactors = FALSE
  )
}
