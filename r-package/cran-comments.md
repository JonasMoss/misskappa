## Resubmission

This is a resubmission. In this version I have:

* Replaced the relative links to the pkgdown articles in the vignette by
  absolute URLs, fixing the "(possibly) invalid file URIs" NOTE.

Since the first submission I also improved the numerics of the
likelihood estimators (a higher default EM iteration limit, clearer handling
of singular covariances, and a fallback for the normal-theory estimator).

## First submission

This is the first CRAN submission of misskappa. Versions up to 2.0.0 were
released on GitHub only.

## Test environments

* local: Ubuntu 26.04, R 4.6.1 (`R CMD check --as-cran`)
* GitHub Actions (`R CMD check --as-cran`): windows-latest (R release),
  macos-latest (R release), ubuntu-latest (R devel), ubuntu-latest
  (R oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

The check also prints two INFO lines:

* "specified C++17": the C++ code needs at least C++17.
* On Linux, the installed size (about 87 Mb, nearly all in `libs/`): the
  package compiles a C++ library using Eigen (via RcppEigen), and the
  shared library carries debugging symbols. It is small on Windows and macOS.

## Notes for the reviewers

* The references in the Description field have DOIs. Two methods papers on
  the missing-data estimators are in preparation; `inst/CITATION` lists them
  as working papers, and they will be updated in a later release.
* Examples, tests, and the vignette do not write files.
