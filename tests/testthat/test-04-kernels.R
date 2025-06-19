context("Kernel Functions")

# Helper function to create test data
create_test_data <- function() {
  set.seed(123)
  list(
    X_vec = rnorm(5),
    X_mat = matrix(rnorm(10), ncol = 2),
    X_mercer = seq(0, 1, length.out = 5)
  )
}

test_that("fbm kernel produces valid Gram matrix", {
  test_data <- create_test_data()

  # Test with vector input
  K_vec <- fbm(test_data$X_vec)
  expect_equal(dim(K_vec), c(5, 5))
  expect_true(isSymmetric(K_vec))

  # Test with matrix input
  K_mat <- fbm(test_data$X_mat)
  expect_equal(dim(K_mat), c(5, 5))

  # Test Hurst parameter
  K_hurst <- fbm(test_data$X_mat, Hurst = 0.7)
  expect_equal(attr(K_hurst, "parameters")$Hurst, 0.7)

  # Test attributes
  expect_equal(attr(K_vec, "kernel"), "fbm")
})

test_that("cfbm kernel produces valid centered Gram matrix", {
  test_data <- create_test_data()

  K <- cfbm(test_data$X_mat)
  expect_equal(dim(K), c(5, 5))
  expect_true(isSymmetric(K))

  # Check centering property (row means ~0)
  expect_lt(max(abs(rowMeans(K))), 1e-10)

  # Test attributes
  expect_equal(attr(K, "kernel"), "cfbm")
})

test_that("cfbm_sd produces standardized output", {
  test_data <- create_test_data()

  K <- cfbm_sd(test_data$X_mat)
  expect_equal(dim(K), c(5, 5))
  expect_true(isSymmetric(K))

  # Test Hurst parameter
  K_hurst <- cfbm_sd(test_data$X_mat, Hurst = 0.7)
  expect_equal(attr(K_hurst, "parameters")$Hurst, 0.7)

  # Test attributes
  expect_equal(attr(K_hurst, "kernel"), "cfbm (standardized)")
})

test_that("kronecker_delta kernel works correctly", {
  test_data <- create_test_data()

  # Test with identical inputs
  X <- c(1, 1, 2, 2)
  K <- kronecker_delta(X)
  expect_equal(diag(K), rep(1, length(X)))
  expect_equal(K[1,2], 1)
  expect_equal(K[1,3], 0)

  # Test centering
  K_centered <- kronecker_delta(X, center = TRUE)
  expect_equal(mean(K_centered), 0, tolerance = 1e-10)
})

test_that("rbf kernel produces valid results", {
  test_data <- create_test_data()

  K <- rbf(test_data$X_mat, lengthscale = 0.5)
  expect_equal(dim(K), c(5, 5))
  expect_true(all(K <= 1 & K >= 0))

  # Test lengthscale effect
  K_small <- rbf(test_data$X_mat, lengthscale = 0.1)
  K_large <- rbf(test_data$X_mat, lengthscale = 1)
  expect_true(mean(K_small) < mean(K_large))
})

test_that("polynomial kernel computes correctly", {
  test_data <- create_test_data()

  # Linear kernel
  K_lin <- polynomial(test_data$X_mat, d = 1)
  expect_equal(dim(K_lin), c(5, 5))

  # Quadratic kernel
  K_quad <- polynomial(test_data$X_mat, d = 2)
  expect_equal(dim(K_quad), c(5, 5))

  # Test offset
  K_offset <- polynomial(test_data$X_mat, d = 2, offset = 10)
  expect_true(all(K_offset > K_quad))
})

test_that("mercer kernel handles inputs correctly", {
  test_data <- create_test_data()

  K <- mercer(test_data$X_mercer, delta = 1)
  expect_equal(dim(K), c(5, 5))
  expect_true(isSymmetric(K))

  # Test delta parameter
  K_delta <- mercer(test_data$X_mercer, delta = 2)
  expect_true(mean(abs(K_delta)) < mean(abs(K)))

  # Test max_terms
  expect_silent(mercer(test_data$X_mercer, max_terms = 500))

  # Invalid delta parameter
  expect_error(
    mercer(test_data$X_mercer, delta = -0.1),
    "The delta parameter must be positive."
  )
})

test_that("All kernels handle edge cases", {
  # Empty input
  expect_error(fbm(matrix(numeric(0))), "X must contains at least one numeric value.")
  expect_error(cfbm(matrix(numeric(0))), "X must contains at least one numeric value.")
  expect_error(cfbm_sd(matrix(numeric(0))), "X must contains at least one numeric value.")
  expect_error(polynomial(matrix(numeric(0))), "X must contains at least one numeric value.")
  expect_error(rbf(matrix(numeric(0))), "X must contains at least one numeric value.")
  expect_error(mercer(matrix(numeric(0))), "X must contains at least one numeric value.")

  # Single observation
  expect_silent(kronecker_delta(1))
  expect_equal(dim(rbf(matrix(1))), c(1, 1))

  # Invalid parameters
  expect_error(fbm(matrix(1), Hurst = -1), "Hurst")
  expect_error(cfbm(matrix(1), Hurst = -1), "Hurst")
  expect_error(cfbm_sd(matrix(1), Hurst = -1), "Hurst")
  expect_error(mercer(1:5, delta = 0), "positive")
})

test_that("Kernel attributes are properly set", {
  test_data <- create_test_data()

  kernels <- list(
    fbm(test_data$X_mat),
    cfbm(test_data$X_mat),
    rbf(test_data$X_mat),
    polynomial(test_data$X_mat, d = 2),
    mercer(test_data$X_mercer)
  )

  for (k in kernels) {
    expect_true(!is.null(attr(k, "kernel")))
    expect_true(!is.null(attr(k, "parameters")))
  }
})
