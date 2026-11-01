## New submission

This is the first CRAN submission of misskappa. Versions up to 2.0.0 were
released on GitHub only.

## Test environments

* local: Ubuntu 26.04, R 4.6.1
* GitHub Actions: windows-latest (R release), macos-latest (R release),
  ubuntu-latest (R devel), ubuntu-latest (R oldrel-1), ubuntu-24.04
  (R release)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

On some platforms the check may also report the installed size (the shared
library is large because it is compiled C++ using Eigen via RcppEigen, with
debugging symbols). The size is dominated by the compiled code in `libs/`.

## Notes for the reviewers

* The two references in the Description field have DOIs. Two methods papers
  describing the missing-data estimators are in preparation and are listed
  in `inst/CITATION` as unpublished working papers; they will be updated in
  a later release.
* Examples, tests, and the vignette do not write files.
