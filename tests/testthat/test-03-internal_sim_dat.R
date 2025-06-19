context("Internal RKHS Utility Functions")

test_that("rkhs_norm works with matrices and lists", {
  X <- 1:4
  Y <- 10:13
  Ginv <- diag(1,4)
  Index <- expand.grid(1:2, 1:2)

  # Test with matrix input
  expect_silent(
    rkhs_norm(X, NULL, Ginv, Index)
  )
  expect_silent(
    rkhs_norm(X, Y, Ginv, Index)
  )

  # Test input validation
  expect_error(
    rkhs_norm(X, NULL, "not_a_matrix", Index),
    "Ginv must be a matrix"
  )
})

test_that("rkhs_norm_mat handles different input types", {
  X <- matrix(1:4, nrow = 2)
  Ginv <- diag(1,2)
  Index <- expand.grid(1:2, 1:2)

  # Test with matrix input
  result_list <- rkhs_norm_mat(X, Ginv, Index)
  expect_equal(dim(result_list), c(2, 2))

  # Test with list input
  X.list = list(X, matrix(2:5, nrow = 2))
  Ginv <- diag(1,4)
  result_mat <- rkhs_norm_mat(X.list, Ginv, Index)
  expect_equal(dim(result_mat), c(2, 2))

  # Test symmetry
  expect_equal(result_mat[1,2], result_mat[2,1])
})

test_that("cfbm_rkhs computes Gram matrix correctly", {
  X <- matrix(1:4, nrow = 2)
  Ginv <- diag(1,2)
  Index <- expand.grid(1:2, 1:2)

  # Test with default Hurst parameter
  result <- cfbm_rkhs(X, Ginv, Index)
  expect_equal(dim(result), c(2, 2))

  # Test with custom Hurst parameter
  result_hurst <- cfbm_rkhs(X, Ginv, Index, Hurst = 0.7)
  expect_equal(dim(result_hurst), c(2, 2))

  # Test input validation
  expect_error(
    cfbm_rkhs(X, Ginv, Index, Hurst = 1.5),
    "Hurst coefficient must be in \\(0,1\\]"
  )

  # Test with list input
  X.list <- list(X, matrix(2:5, nrow = 2))
  Ginv <- diag(1,4)
  result_mat <- cfbm_rkhs(X.list, Ginv, Index)
  expect_equal(dim(result_mat), c(2, 2))

  # Test symmetry
  expect_equal(result_mat[1,2], result_mat[2,1])
})

test_that("cfbm_rkhs_cross handles train-new comparisons", {
  X <- matrix(1:4, nrow = 2)
  Y <- matrix(2:3, nrow = 1)
  Ginv <- diag(1,2)
  Index <- expand.grid(1:2, 1:2)
  # Test matrix inputs
  result_mat <- cfbm_rkhs_cross(X, Y, Ginv, Index)
  expect_equal(dim(result_mat), c(1, 2))

  # Test list inputs
  X.list <- list(X, matrix(2:5, nrow = 2))
  Y.list <- list(matrix(3:6, nrow = 2), matrix(4:7, nrow = 2))
  Ginv <- diag(1,4)
  result <- cfbm_rkhs_cross(X.list, Y.list, Ginv, Index)
  expect_equal(dim(result), c(2, 2))

  # Test input validation
  expect_error(
    cfbm_rkhs_cross(X.list, matrix(3:6, nrow = 2), Ginv, Index),
    "must be of the same class"
  )
})

test_that("Edge cases are handled properly", {

  Ginv <- diag(1,4)
  Index <- expand.grid(1:2, 1:2)

  # Empty inputs
  expect_error(
    rkhs_norm(matrix(numeric(0)), NULL, Ginv, Index)
  )

  # Single element
  X_single <- matrix(1, ncol = 1)
  Index_single <- rep(1, 2)
  expect_silent(
    rkhs_norm_mat(X_single, diag(1), Index_single)
  )

  # Invalid indices
  X <- matrix(1:4, nrow = 2)
  Ginv <- diag(1,2)
  expect_error(
    rkhs_norm(X, NULL, Ginv, 'foo'),
    "Index matrix must be integer"
  )
})
