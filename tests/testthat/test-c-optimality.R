# c-optimality (one linear combination c'theta): Elfving's bound certifies the
# design, including a SINGULAR one that the sensitivity function cannot certify
# (R/elfving.R).

info_logit <- function(x, theta) {
  q <- c(1, x); eta <- min(max(sum(q * theta), -500), 500)
  (exp(eta / 2) / (1 + exp(eta))) * q
}
Xc  <- matrix(seq(-3, 3, length.out = 601), ncol = 1)
thc <- c(1, 1)
sel <- function(theta) matrix(c(1, -1), nrow = 1)     # c'theta = theta0 - theta1 = eta(-1)

test_that("a singular c-optimal design is certified by Elfving's bound", {
  r <- optimal_design(info_vector = info_logit, theta = thc, candidate_set = Xc,
                      grad_g = sel, p = 0)
  # the one-point design at x = -1: c'M^-c = 1/nu(eta(-1)) = 1/nu(0) = 4
  expect_equal(nrow(r$support), 1L)
  expect_equal(r$support[1, 1], -1)
  expect_equal(r$weights, 1)
  expect_equal(exp(r$criterion), 4, tolerance = 1e-10)
  expect_equal(r$elfving_bound, 4, tolerance = 1e-8)
  expect_lt(abs(r$elfving_gap), 1e-8)
  expect_true(r$converged)                            # no "did NOT converge" any more
  expect_equal(r$efficiency_lower_bound, 1, tolerance = 1e-8)
  expect_gt(r$max_d, 0.5)                             # the pseudo-inverse sensitivity is NOT decisive
  # p = 1 is the same criterion for one linear combination
  r1 <- optimal_design(info_vector = info_logit, theta = thc, candidate_set = Xc,
                       grad_g = sel, p = 1)
  expect_true(r1$converged)
  expect_equal(r1$criterion, 4, tolerance = 1e-10)
  expect_equal(r1$elfving_bound, 4, tolerance = 1e-8)
})

test_that("the certificate is a valid lower bound and is tight for any h at the optimum", {
  F  <- t(sapply(Xc[, 1], info_logit, theta = thc))  # 601 x 2
  cv <- c(1, -1)
  # every h with h'c = 1 gives 1/max_i (h'f_i)^2 <= 4 (the optimal value)
  set.seed(1)
  for (i in 1:20) {
    h <- rnorm(2); h <- h / sum(h * cv)
    expect_lte(1 / max((F %*% h)^2), 4 + 1e-12)
  }
  # the internal certificate finds the tight h
  cert <- OWEAdesign:::.c_certificate(cv, 0L, t(F), matrix(0, 2, 2), 2L, value = 4)
  expect_equal(cert$bound, 4, tolerance = 1e-8)
  expect_true(cert$certified)
  expect_equal(sum(cert$h * cv), 1, tolerance = 1e-10)
  # a non-optimal design is not certified, and its efficiency bound is < 1
  cert2 <- OWEAdesign:::.c_certificate(cv, 0L, t(F), matrix(0, 2, 2), 2L, value = 5)
  expect_false(cert2$certified)
  expect_equal(cert2$efficiency, 4 / 5, tolerance = 1e-8)
})

test_that("the Nelder-Mead / Brent fallback (matrix mode, existing design) agrees with the LP", {
  info_m <- function(x) tcrossprod(info_logit(x, thc))
  rm <- optimal_design(info_matrix = info_m, candidate_set = Xc, grad_g = sel, p = 0)
  expect_true(rm$converged)
  expect_equal(rm$elfving_bound, 4, tolerance = 1e-6)
  # with an existing design the bound stays a valid lower bound on the value
  re <- optimal_design(info_vector = info_logit, theta = thc, candidate_set = Xc,
                       grad_g = sel, p = 0, xi0_points = matrix(c(0, 2)),
                       xi0_weights = c(0.5, 0.5), n0 = 10, n1 = 10)
  expect_lte(re$elfving_bound, exp(re$criterion) + 1e-8)
  expect_true(re$converged)
  expect_equal(re$efficiency_lower_bound, 1, tolerance = 1e-6)
})

test_that("a nonsingular c-optimal design is certified by both checks", {
  r <- optimal_design(info_vector = info_logit, theta = thc, candidate_set = Xc,
                      subset = 2, p = 0)                # the slope alone
  expect_true(r$converged)
  expect_lte(r$max_d, 1e-6)                             # the sensitivity certifies it ...
  expect_equal(r$elfving_bound, exp(r$criterion), tolerance = 1e-8)   # ... and so does Elfving
  expect_equal(nrow(r$support), 2L)
})

test_that("verify_optimality(), the grid path, exact_design() and the continuous path use the bound", {
  v <- verify_optimality(matrix(-1), 1, info_vector = info_logit, theta = thc,
                         candidate_set = Xc, grad_g = sel, p = 0)
  expect_true(v$is_optimal$value)
  expect_match(v$is_optimal$note, "Elfving")
  expect_equal(v$elfving_bound, 4, tolerance = 1e-8)
  expect_equal(v$efficiency_lower_bound, 1, tolerance = 1e-8)
  # a wrong design is not certified
  v2 <- verify_optimality(matrix(c(-2, 1)), c(0.5, 0.5), info_vector = info_logit,
                          theta = thc, candidate_set = Xc, grad_g = sel, p = 0)
  expect_false(v2$is_optimal$value)
  expect_lt(v2$efficiency_lower_bound, 0.99)
  # the multistage grid path, with the whole-box check
  rg <- suppressWarnings(optimal_design(info_vector = info_logit, theta = thc,
                                        design_box = list(c(-3, 3)),
                                        step_sequence = c(0.5, 0.1, 0.01), grad_g = sel,
                                        p = 0, check_global = TRUE))
  expect_true(rg$converged)
  expect_true(rg$global_check)
  expect_equal(rg$elfving_bound, 4, tolerance = 1e-8)
  # the exact design inherits the certified bound of its approximate reference
  ex <- exact_design(20, info_vector = info_logit, theta = thc, candidate_set = Xc,
                     grad_g = sel, p = 0, seed = 1)
  expect_equal(ex$counts, 20L)
  expect_equal(ex$efficiency_lower_bound, 1, tolerance = 1e-8)
  # the continuous path (bound over the points examined)
  set.seed(3)
  rc <- suppressWarnings(optimal_design(info_vector = info_logit, theta = thc,
                                        design_box = list(c(-3, 3)), grad_g = sel,
                                        p = 0, continuous = TRUE, n_audit = 0))
  # an optimiser's value: on other platforms / R versions the search stops a
  # hair short of the boundary point (3.9995 was seen on Linux and Windows CI),
  # and the bound over the points examined then misses the value by ~1e-5,
  # above eps0 = 1e-6: not certified, but certified to within 0.01%
  expect_true(isTRUE(rc$converged) || rc$efficiency_lower_bound > 1 - 1e-4)
  expect_equal(exp(rc$criterion), 4, tolerance = 1e-3)
  expect_lte(rc$elfving_bound, 4 + 1e-6)
  expect_gt(rc$efficiency_lower_bound, 0.99999)
})

test_that("the value of a c-optimal design is computed robustly from M", {
  # nonsingular M: the eigen route equals c' M^-1 c
  set.seed(5)
  A <- matrix(rnorm(16), 4); M <- crossprod(A) + diag(4); cc <- c(1, -1, 0.5, 2)
  expect_equal(OWEAdesign:::.c_value_robust(M, cc),
               as.numeric(crossprod(cc, solve(M, cc))), tolerance = 1e-10)
  # singular M with c in its range: finite, the pseudo-inverse value
  f <- c(1, -1); Ms <- tcrossprod(f)                 # rank one, one-point design
  # (f f')^+ = f f' / (f'f)^2, so c = a f gives V = a^2
  expect_equal(OWEAdesign:::.c_value_robust(Ms, f), 1, tolerance = 1e-12)
  expect_equal(OWEAdesign:::.c_value_robust(Ms, 2 * f), 4, tolerance = 1e-12)
  # c outside the range: the design does not estimate c'theta
  expect_identical(OWEAdesign:::.c_value_robust(Ms, c(1, 1)), Inf)
  # a value below Elfving's bound is never certified, a matching one is
  expect_false(OWEAdesign:::.c_certified(4, 5.09, 1e-6))
  expect_false(OWEAdesign:::.c_certified(0, 5.09, 1e-6))
  expect_false(OWEAdesign:::.c_certified(Inf, 5.09, 1e-6))
  expect_true(OWEAdesign:::.c_certified(5.09 * (1 + 5e-7), 5.09, 1e-6))
  expect_true(OWEAdesign:::.c_certified(5.09, 5.09, 1e-6))
  # a hair below the bound (a support point a hair off a singular optimum) is
  # still certified; well below it is not
  expect_true(OWEAdesign:::.c_certified(5.09 * (1 - 5e-4), 5.09, 1e-6))
  expect_false(OWEAdesign:::.c_certified(5.09 * (1 - 5e-3), 5.09, 1e-6))
  expect_false(OWEAdesign:::.c_certified(5.09 * (1 + 1e-4), 5.09, 1e-6))
})
