#!/usr/bin/env Rscript
# R package test coverage via covr.
#
# Covers r-package/R/ (the R layer) and, because covr recompiles the package
# with --coverage, the vendored C++ in r-package/src/ as exercised *through*
# Rcpp -- a second view of the C++ from the R side, complementing the standalone
# doctest coverage from dev/cpp-coverage.sh.
#
# Run `just vendor` first (the `just r-cov` recipe does) so r-package/src is
# in sync with the canonical src/ + include/.
#
# Noble binary mirror so curl/httr (covr's deps) resolve without compiling.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))
options(HTTPUserAgent = sprintf(
  "R/%s R (%s)", getRversion(),
  paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
))

suppressMessages(library(covr))

# type = "tests" runs only the testthat suite (not the heavy quarto vignettes).
cov <- package_coverage("r-package", type = "tests", quiet = FALSE)

cat("\n===================== R coverage (testthat) =====================\n")
print(cov)
cat(sprintf("\nTotal line coverage: %.2f%%\n", percent_coverage(cov)))

zero <- covr::zero_coverage(cov)
if (nrow(zero) > 0) {
  cat("\nUncovered lines, counted per file:\n")
  print(sort(table(zero$filename), decreasing = TRUE))
}

# HTML report needs DT; emit it only when available.
if (requireNamespace("DT", quietly = TRUE)) {
  covr::report(cov, file = "build-cov/r-coverage.html", browse = FALSE)
  cat("\nHTML report: build-cov/r-coverage.html\n")
} else {
  cat("\n(install.packages(\"DT\") for an HTML report)\n")
}
