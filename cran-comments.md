## Submission of OWEAdesign 0.6.0 (new package)

### Test environments

* Windows 11, R 4.3.2 (local): `R CMD check --as-cran`
* win-builder (R-devel and R-release) and macOS builder: see the submission
  comments for the results at the time of submission

### R CMD check results

0 errors | 0 warnings | 1 note

* "New submission": this is the package's first submission.

### Comments for the CRAN team

* The package implements the Optimal Weights Exchange Algorithm of Yang,
  Biedermann and Tang (2013) <doi:10.1080/01621459.2013.806268> for general
  linear, nonlinear and generalized linear models. It is unrelated to the
  existing CRAN package OWEA (crossover designs), apart from the shared
  algorithm name; that is why this package is called OWEAdesign.
* The compiled code is C++ via Rcpp / RcppArmadillo; no system requirements
  beyond a C++ compiler.
* The slower tests (a few Shiny server tests and the alternative solver's
  tests) are skipped on CRAN with `skip_on_cran()`; the remaining tests run
  in well under a minute.
