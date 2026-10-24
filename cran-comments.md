## Resubmission of OWEAdesign 0.6.0 (new package)

This is a resubmission. Compared with the version submitted on 2026-10-18 it
addresses the two points raised in the CRAN reply:

* Software names in the Description field are now single-quoted: 'C++',
  'RcppArmadillo' and 'Shiny'.
* README.md no longer contains the relative file URI `../deploy/deploy.R`,
  which pointed outside the package. The README's other relative links (to
  the `slides/` folder, which is not part of the package) were replaced by
  absolute URLs to the package's GitHub repository.

### Test environments

* Windows 11, R 4.3.2 (local): `R CMD check --as-cran`
* GitHub Actions (R-devel, R-release and R-oldrel on Ubuntu; R-release on
  macOS and Windows): https://github.com/minyangstat-eng/OWEAdesign/actions

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
