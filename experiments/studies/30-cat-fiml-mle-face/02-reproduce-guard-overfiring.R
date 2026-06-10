library(misskappa)
run <- function(n, C, R, seed) {
  set.seed(seed)
  z <- sample(1:C, n, TRUE)
  rate <- function(z) ifelse(runif(n) < 0.7, z, sample(1:C, n, TRUE))
  x <- sapply(1:R, function(j) rate(z))
  # planned missingness: each subject sees a random pair of raters
  pairs <- t(combn(R, 2))
  for (i in 1:n) {
    p <- pairs[sample(nrow(pairs), 1), ]
    x[i, setdiff(1:R, p)] <- NA
  }
  res <- try(misskappa::kappa(x, estimator = "cat_fiml", weight = "nominal"),
             silent = TRUE)
  if (inherits(res, "try-error")) {
    sprintf("n=%d C=%d R=%d seed=%d: ERROR: %s", n, C, R, seed,
            trimws(attr(res, "condition")$message))
  } else {
    sprintf("n=%d C=%d R=%d seed=%d: conger=%.4f se=%.4f", n, C, R, seed,
            res$estimates[1], sqrt(res$vcov[1,1]))
  }
}
for (cfg in list(c(40,3,3), c(40,4,4), c(40,5,4), c(100,4,4), c(200,5,5)))
  for (s in 1:3)
    cat(run(cfg[1], cfg[2], cfg[3], s), "\n")
