context("RKHS Kronecker Functions")

test_that("rkhs_norm_kron works for all cases", {

  # Test vector input
  G <- list(matrix(c(2,1,1,2), nrow = 2))
  decomp <- fire:::Kronecker_inv_helper(G, alpha = c(0.5))
  L_inv <- fire:::Lambda_inv(list(decomp$L1))
  x <- c(1,2)
  expect_silent(rkhs_norm_kron(x, G_list = decomp, L_inv = L_inv))
  expect_true(is.numeric(rkhs_norm_kron(x, G_list = decomp, L_inv = L_inv)))

  # Test matrix input
  G <- list(matrix(c(2,1,1,2), nrow=2), matrix(c(2,3,3,2), nrow=2))
  decomp <- fire:::Kronecker_inv_helper(G, alpha=c(0.5, 0.5))
  L_inv <- fire:::Lambda_inv(list(decomp$L1,decomp$L2))
  X <- matrix(rnorm(4), nrow=2)
  expect_silent(rkhs_norm_kron(X, G_list=decomp, L_inv=L_inv))

  # Test with Y input
  Y <- matrix(rnorm(4), nrow=2)
  expect_silent(rkhs_norm_kron(X, Y, G_list=decomp, L_inv=L_inv))

  # Test error conditions
  expect_error(rkhs_norm_kron("not numeric", G_list=decomp, L_inv=L_inv))
  expect_error(rkhs_norm_kron(X, G_list=list(), L_inv=L_inv))
})

test_that("Kronecker_norm_mat computes pairwise norms correctly", {
  # Setup
  G <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2))
  X <- list(matrix(rnorm(4),2), matrix(rnorm(4),2))

  # Test basic functionality
  result <- Kronecker_norm_mat(X, G, alpha = c(0.5,0.5))
  expect_equal(dim(result), c(2,2))
  expect_true(all(diag(result) >= 0)) # Diagonal should be non-negative

  # Test symmetry
  expect_equal(result[1,2], result[2,1])

  # Test parallel options
  expect_silent(Kronecker_norm_mat(X, G, alpha=c(0.5,0.5), os_type="Apple"))
  expect_silent(Kronecker_norm_mat(X, G, alpha=c(0.5,0.5), os_type="Windows"))
})

test_that("Kronecker_norm_cross computes cross norms correctly", {
  # Setup
  G <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2))
  Xtrain <- list(matrix(rnorm(4),2), matrix(rnorm(4),2))
  Xnew <- list(matrix(rnorm(4),2))

  # Test basic functionality
  result <- Kronecker_norm_cross(Xtrain, Xnew, G, alpha=c(0.5,0.5))
  expect_equal(dim(result), c(1,2))

  # Test with precomputed G_list
  decomp <- fire:::Kronecker_inv_helper(G, alpha=c(0.5,0.5))
  result2 <- Kronecker_norm_cross(Xtrain, Xnew, G, alpha=c(0.5,0.5), G_list=decomp)
  expect_equal(result, result2)
})

test_that("Special cases work correctly", {
  # Test m=1 case
  G1 <- list(matrix(c(2,1,1,2),2))
  x <- c(1,2)

  # With constant term
  decomp1 <- fire:::Kronecker_inv_helper(G1, alpha = 0.5, constant = TRUE)
  L_inv1 <- fire:::Lambda_inv(list(decomp1$L1))
  expect_silent(rkhs_norm_kron(x, G_list=decomp1, L_inv=L_inv1))

  # Without constant term
  decomp2 <- fire:::Kronecker_inv_helper(G1, alpha = 0.5, constant = FALSE)
  expect_silent(rkhs_norm_kron(x, G_list = decomp2, Ginv = MASS::ginv(decomp2$G1), constant = FALSE))
})

test_that("Error handling works", {
  # Dimension mismatches
  G <- list(matrix(c(2,1,1,2),2), matrix(c(2,3,3,2),2))
  X <- matrix(rnorm(6), 3) # Wrong dimensions
  Xnew <- matrix(rnorm(6), 3)
  expect_error(Kronecker_norm_mat(X, G, alpha = c(0.5, 0.5), validate = T),
               "Dimension mismatch")
  expect_error(Kronecker_norm_mat(X, Xnew, G, alpha = c(0.5, 0.5), validate = T),
               "Dimension mismatch")
  # Invalid os_type
  Xlist <- list(matrix(rnorm(4),2), matrix(rnorm(4),2))
  expect_error(Kronecker_norm_mat(Xlist, G, alpha=c(0.5,0.5), os_type="Linux", validate = T),
               "Invalid os_type")
})
