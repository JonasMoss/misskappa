# The rating scale (declared categories, labels, factors), subject handling,
# and degenerate inputs.

test_that("declared categories keep unused categories in the scale", {
  set.seed(2)
  x <- matrix(sample(1:4, 120, TRUE, prob = c(.4, .3, .2, .1)), 40)
  obs <- agreement(x, estimator = "ipw")
  full <- agreement(x, estimator = "ipw", categories = 1:5)
  # Conger and Fleiss under nominal weights do not depend on unused
  # categories; Brennan-Prediger's uniform chance baseline does.
  expect_equal(coef(full)[c("Conger", "Fleiss")], coef(obs)[c("Conger", "Fleiss")])
  p_o <- 1 - (1 - coef(obs)[["Brennan-Prediger"]]) * (1 - 1 / 4)
  expect_equal(coef(full)[["Brennan-Prediger"]], (p_o - 1 / 5) / (1 - 1 / 5))
  expect_equal(full$categories, 1:5)

  # Linear weights are scaled by the range of the declared scale.
  lin4 <- agreement(x, estimator = "ipw", weight = "linear")
  lin5 <- agreement(x, estimator = "ipw", weight = "linear", categories = 1:5)
  expect_false(isTRUE(all.equal(coef(lin4), coef(lin5))))

  # cat_fiml accepts the declared scale too.
  fiml <- agreement(x, estimator = "cat_fiml", categories = 1:5)
  expect_equal(coef(fiml), coef(full), tolerance = 1e-6)

  expect_error(agreement(x, estimator = "ipw", categories = 1:3),
               "absent from 'categories'")
  expect_error(agreement(x, estimator = "ipw", categories = c(1, 1, 2)),
               "distinct")
  expect_error(agreement(x, estimator = "ipw", categories = 1:5, values = 1:4),
               "number of categories")
})

test_that("scored estimators report Brennan-Prediger on a declared scale", {
  set.seed(3)
  x <- matrix(sample(1:4, 160, TRUE), 40)
  undeclared <- agreement(x, estimator = "nt_fiml")
  expect_false("Brennan-Prediger" %in% names(coef(undeclared)))
  declared <- agreement(x, estimator = "nt_fiml", categories = 1:5)
  by_values <- agreement(x, estimator = "nt_fiml", categories = 1:5,
                         values = 1:5)
  expect_true("Brennan-Prediger" %in% names(coef(declared)))
  expect_equal(coef(declared), coef(by_values))
  # Declaring an unused fifth category widens the BP chance baseline.
  four <- agreement(x, estimator = "nt_fiml", values = 1:4)
  expect_gt(coef(declared)[["Brennan-Prediger"]],
            coef(four)[["Brennan-Prediger"]])
  expect_equal(coef(declared)[c("Conger", "Fleiss")],
               coef(four)[c("Conger", "Fleiss")])
})

test_that("labels can be strings, factors, or non-integer numbers", {
  set.seed(4)
  codes <- matrix(sample(1:3, 90, TRUE), 30)
  codes[sample(90, 10)] <- NA
  ref <- agreement(codes, estimator = "ipw", weight = "linear")

  chr <- matrix(c("lo", "mid", "hi")[codes], 30)
  fit_chr <- agreement(chr, estimator = "ipw", weight = "linear",
                       categories = c("lo", "mid", "hi"))
  expect_equal(coef(fit_chr), coef(ref))

  df <- as.data.frame(lapply(seq_len(3), function(j)
    factor(chr[, j], levels = c("lo", "mid", "hi"))))
  fit_fac <- agreement(df, estimator = "ipw", weight = "linear")
  expect_equal(coef(fit_fac), coef(ref))
  expect_equal(fit_fac$categories, c("lo", "mid", "hi"))

  # Non-integer numeric labels are distinct categories, never truncated.
  halves <- codes / 2
  fit_half <- agreement(halves, estimator = "ipw", weight = "linear")
  expect_equal(coef(fit_half), coef(ref))
  expect_equal(fit_half$categories, c(0.5, 1, 1.5))
  merged <- matrix(c(1.1, 1.2, 1.1, 1.2, 2, 2), 3)
  expect_equal(length(agreement(merged, estimator = "ipw")$categories), 3L)
})

test_that("subjects with no ratings are dropped by every estimator", {
  set.seed(5)
  X <- matrix(sample(1:5, 200, TRUE), 50)
  X[1, ] <- NA
  X[2:10, 1] <- NA
  ns <- vapply(c("ipw", "cat_fiml", "pairwise", "nt_fiml", "nt_fiml_cs"),
               function(e) suppressWarnings(agreement(X, estimator = e))$n,
               numeric(1))
  expect_true(all(ns == 49))
  fit <- agreement(X, estimator = "ipw")
  expect_equal(rownames(fit$psi), as.character(2:50))
  expect_equal(agreement_counts(ratings_to_counts(X))$n, 49L)

  # Paired tests line up subjects by id, so estimators with different
  # missing rows still compare.
  a <- agreement(X, estimator = "pairwise")
  b <- agreement(X[-(11:12), ], estimator = "pairwise")
  expect_s3_class(agreement_test(a, b, paired = TRUE), "htest")
})

test_that("paired covariance matches joint_vcov on aligned fits", {
  set.seed(6)
  X <- matrix(rnorm(300), 100)
  a <- agreement(X[, 1:2], estimator = "pairwise")
  b <- agreement(X[, 2:3], estimator = "pairwise")
  V <- research_joint_vcov(a = a, b = b)
  t1 <- agreement_test(a, b, paired = TRUE)
  d <- coef(a)[["Conger"]] - coef(b)[["Conger"]]
  v <- V["a.Conger", "a.Conger"] + V["b.Conger", "b.Conger"] -
    2 * V["a.Conger", "b.Conger"]
  expect_equal(unname(t1$statistic), d^2 / v, tolerance = 1e-10)
})

test_that("constant ratings give undefined coefficients with a warning", {
  x <- matrix(1, 6, 3)
  expect_warning(fit <- agreement(x, estimator = "ipw"), "undefined")
  expect_true(all(is.nan(coef(fit))))
  expect_true(all(is.na(vcov(fit))))
})

test_that("errors do not name internal functions", {
  err <- tryCatch(agreement(matrix(1:6, 3), estimator = "ipw", g = 0),
                  error = identity)
  expect_null(conditionCall(err))
  err <- tryCatch(agreement(matrix(1:6, 3), estimator = "ipw", values = 1:5),
                  error = identity)
  expect_null(conditionCall(err))
})

test_that("the exported surface is the four public functions", {
  expect_setequal(getNamespaceExports("misskappa"),
                  c("agreement", "agreement_counts", "agreement_test",
                    "ratings_to_counts"))
  ns <- ls(asNamespace("misskappa"), all.names = TRUE)
  research <- grep("^research_", ns, value = TRUE)
  expect_true(length(research) > 0L)
  expect_length(intersect(research, getNamespaceExports("misskappa")), 0L)
})

test_that("text ratings need a declared order for order-dependent weights", {
  set.seed(7)
  x <- matrix(sample(c("low", "mid", "high"), 60, TRUE), 20)
  expect_s3_class(agreement(x), "misskappa_estimate")          # nominal is fine
  expect_error(agreement(x, weight = "linear"), "no declared order")
  expect_error(agreement(x, estimator = "nt_fiml"), "no declared order")
  expect_error(agreement(x, weight = "linear", g = 3), "no declared order")
  fit <- agreement(x, weight = "linear", categories = c("low", "mid", "high"))
  expect_equal(fit$categories, c("low", "mid", "high"))

  # Factor columns with different level sets give no shared order either.
  df <- data.frame(a = factor(x[, 1], levels = c("low", "mid", "high")),
                   b = factor(x[, 2], levels = c("low", "mid", "high", "x")),
                   c = factor(x[, 3], levels = c("low", "mid", "high")))
  expect_error(agreement(df, weight = "linear"), "no declared order")
})

test_that("ratings_to_counts keeps shared factor levels as column order", {
  df <- data.frame(a = factor(c("low", "high", "mid"), levels = c("low", "mid", "high")),
                   b = factor(c("mid", "high", "mid"), levels = c("low", "mid", "high")))
  expect_equal(colnames(ratings_to_counts(df)), c("low", "mid", "high"))
})

test_that("subject ids come from explicit row names, and must be unique", {
  set.seed(8)
  d <- as.data.frame(matrix(sample(1:3, 18, TRUE), 6))
  row.names(d) <- c(10L, 20L, 30L, 40L, 50L, 60L)
  expect_equal(rownames(agreement(d)$psi), as.character(c(10, 20, 30, 40, 50, 60)))
  expect_equal(rownames(agreement(as.data.frame(unname(as.matrix(d))))$psi),
               as.character(1:6))

  x <- matrix(sample(1:3, 18, TRUE), 6)
  rownames(x) <- rep("same", 6)
  a <- agreement(x)
  expect_error(agreement_test(a, a, paired = TRUE), "duplicated row names")
})
