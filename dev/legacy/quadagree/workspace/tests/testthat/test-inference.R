x <- dat.zapf2016
n <- nrow(x)
sigma <- stats::cov(x) * (n - 1) / n
mu <- colMeans(x)

alleq <- \(x, y) isTRUE(all.equal(x, y))

fleisscis <- list(
  fleissci(x, type = "adf"),
  fleissci(x, type = "elliptical"),
  fleissci(x, type = "normal")
)

fleissci_fishers <- list(
  fleissci(x, type = "adf", transform = "fisher"),
  fleissci(x, type = "elliptical", transform = "fisher"),
  fleissci(x, type = "normal", transform = "fisher")
)

congercis <- list(
  congerci(x, type = "adf"),
  congerci(x, type = "elliptical"),
  congerci(x, type = "normal")
)

congerci_fishers <- list(
  congerci(x, type = "adf", transform = "fisher"),
  congerci(x, type = "elliptical", transform = "fisher"),
  congerci(x, type = "normal", transform = "fisher")
)

bpcis <- list(
  bpci(x, type = "adf"),
  bpci(x, type = "elliptical"),
  bpci(x, type = "normal")
)

bpci_fishers <- list(
  bpci(x, type = "adf", transform = "fisher"),
  bpci(x, type = "elliptical", transform = "fisher"),
  bpci(x, type = "normal", transform = "fisher")
)

y <- dat.fleiss1971

fleissci_aggrs <- list(
  fleissci_aggr(y, transform = "none"),
  fleissci_aggr(y, transform = "fisher"),
  fleissci_aggr(y, transform = "arcsin")
)

bpci_aggrs <- list(
  bpci_aggr(y, transform = "none"),
  bpci_aggr(y, transform = "fisher"),
  bpci_aggr(y, transform = "arcsin")
)


test_that("fleissci yield different results", {
  for (i in seq(length(fleisscis) - 1)) {
    expect_false(alleq(fleisscis[[i]], fleisscis[[i + 1]]))
  }
})

test_that("bpci yield different results", {
  for (i in seq(length(bpcis) - 1)) {
    expect_false(alleq(bpcis[[i]], bpcis[[i + 1]]))
  }
})

test_that("bpci_aggrs yield different results", {
  for (i in seq(length(bpci_aggrs) - 1)) {
    expect_false(alleq(bpci_aggrs[[i]], bpci_aggrs[[i + 1]]))
  }
})

test_that("fleissci_aggrs yield different results", {
  for (i in seq(length(fleissci_aggrs) - 1)) {
    expect_false(alleq(fleissci_aggrs[[i]], fleissci_aggrs[[i + 1]]))
  }
})

test_that("fleissci fisher transform does something", {
  for (i in seq(length(fleisscis) - 1)) {
    expect_false(alleq(fleisscis[[i]], fleissci_fishers[[i]]))
  }
})

test_that("congerci fisher transform does something", {
  for (i in seq(length(fleisscis) - 1)) {
    expect_false(alleq(congercis[[i]], congerci_fishers[[i]]))
  }
})

test_that("bpci fisher transform does something", {
  for (i in seq(length(fleisscis) - 1)) {
    expect_false(alleq(bpcis[[i]], bpci_fishers[[i]]))
  }
})

test_that("congerci yield different results", {
  for (i in seq(length(congercis) - 1)) {
    expect_false(alleq(congercis[[i]], congercis[[i + 1]]))
  }
})

test_that("fleissci and congerci yield different results", {
  for (i in seq(length(congercis) - 1)) {
    expect_false(alleq(congercis[[i]], fleisscis[[i]]))
  }
})

test_that("bpci and congerci yield different results", {
  for (i in seq(length(congercis) - 1)) {
    expect_false(alleq(congercis[[i]], bpcis[[i]]))
  }
})

test_that("bpci_aggr and fleissci_aggr yield different results", {
  for (i in seq(length(fleissci_aggrs) - 1)) {
    expect_false(alleq(fleissci_aggrs[[i]], bpci_aggrs[[i]]))
  }
})


test_that("fleissci bootstrap does something", {
  expect_false(alleq(
    fleisscis[[1]],
    fleissci(x, type = "adf", bootstrap = TRUE, n_reps = 10)
  ))
})

test_that("congerci bootstrap does something", {
  expect_false(alleq(
    congercis[[1]],
    congerci(x, type = "adf", bootstrap = TRUE, n_reps = 10)
  ))
})

test_that("bpci bootstrap does something", {
  expect_false(alleq(
    bpcis[[1]],
    bpci(x, type = "adf", bootstrap = TRUE, n_reps = 10)
  ))
})

test_that("fleissci_aggr bootstrap does something", {
  expect_false(alleq(
    fleissci_aggrs[[1]],
    fleissci_aggr(x, bootstrap = TRUE, n_reps = 10)
  ))
})

test_that("bpci_aggr bootstrap does something", {
  expect_false(alleq(
    bpci_aggrs[[1]],
    bpci_aggr(x, bootstrap = TRUE, n_reps = 10)
  ))
})

test_that("print is invisible", {
  capture.output(expect_invisible(print(congercis[[1]])))
})

test_that("bp returns error for data with too many missing values", {
  expect_error(bpci(dat.gwet2014))
})

test_that("bp does not return error for data with missing values", {
  expect_no_error(bpci(dat.klein2018))
})

test_that("fleiss returns error for data with too many missing values", {
  expect_error(fleissci(dat.gwet2014))
})


test_that("conger returns error for data with too many missing values", {
  expect_error(congerci(dat.gwet2014))
})
