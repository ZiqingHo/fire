context("I-prior Gram Matrix Functions")

test_that("rbf_rkhs_kron produces valid RBF Gram matrices", {
  # Create test RKHS norm matrix
  nmat <- matrix(c(0, 1, 1, 1, 0, 0.5, 1, 0.5, 0), nrow = 3)

  # Test default lengthscale
  K1 <- rbf_rkhs_kron(nmat)
  expect_equal(dim(K1), c(3, 3))
  expect_true(all(diag(K1) == 1))  # Diagonal should be 1 (exp(0))
  expect_true(all(K1 >= 0 & K1 <= 1))  # All values between 0 and 1
  expect_true(isSymmetric(K1))

  # Test different lengthscale
  K2 <- rbf_rkhs_kron(nmat, lengthscale = 0.5)
  expect_true(all(K2 <= K1))  # Smaller lengthscale -> smaller values

  # Test invalid inputs
  expect_error(rbf_rkhs_kron("not a matrix"))
  expect_error(rbf_rhs_kron(nmat, lengthscale = -1))  # Negative lengthscale
})

test_that("rbf_rkhs_kron_cross produces valid cross RBF matrices", {
  nmat_cross <- matrix(c(0.5, 1.2, 0.8, 1.5), nrow = 2)

  K_cross <- rbf_rkhs_kron_cross(nmat_cross)
  expect_equal(dim(K_cross), c(2, 2))
  expect_true(all(K_cross > 0 & K_cross <= 1))

  # Test lengthscale effect
  K_cross_small <- rbf_rkhs_kron_cross(nmat_cross, lengthscale = 0.1)
  expect_true(all(K_cross_small < K_cross))
})

test_that("cfbm_rkhs_kron produces valid CFBM Gram matrices", {
  nmat <- matrix(c(0, 1, 1, 1, 0, 0.5, 1, 0.5, 0), nrow = 3)

  # Test default Hurst
  G1 <- cfbm_rkhs_kron(nmat)
  expect_equal(dim(G1), c(3, 3))
  expect_true(isSymmetric(G1))
  # Test different Hurst parameters
  G2 <- cfbm_rkhs_kron(nmat, Hurst = 0.8)
  expect_false(identical(G1, G2))

  # Test invalid inputs
  expect_error(cfbm_rkhs_kron(nmat, Hurst = 1.5))  # Hurst > 1
  expect_error(cfbm_rkhs_kron(nmat, Hurst = -0.1)) # Hurst < 0
})

test_that("cfbm_rkhs_kron_cross produces valid cross CFBM matrices", {
  nmat_train <- matrix(c(0, 1, 1, 0), nrow = 2)
  nmat_cross <- matrix(c(0.5, 1.2), nrow = 1)

  G_cross <- cfbm_rkhs_kron_cross(nmat_train, nmat_cross)
  expect_equal(dim(G_cross), c(1, 2))

  # Test Hurst parameter
  G_cross_hurst <- cfbm_rkhs_kron_cross(nmat_train, nmat_cross, Hurst = 0.3)
  expect_false(identical(G_cross, G_cross_hurst))
})

test_that("All functions maintain symmetry when appropriate", {
  nmat <- matrix(c(0, 1, 1, 0), nrow = 2)

  # RBF kernel should be symmetric
  K <- rbf_rkhs_kron(nmat)
  expect_equal(K, t(K))

  # CFBM kernel should be symmetric
  G <- cfbm_rkhs_kron(nmat)
  expect_equal(G, t(G))
})

test_that("Functions work together in pipeline", {
  # Simulate full pipeline
  G <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2))
  X_train <- list(matrix(rnorm(4),2), matrix(rnorm(4),2))
  X_new <- list(matrix(rnorm(4),2))

  # Compute RKHS norms
  nmat <- Kronecker_norm_mat(X_train, G, alpha = c(0.5,0.5))
  ncross <- Kronecker_norm_cross(X_train, X_new, G, alpha = c(0.5,0.5))

  # Test RBF pipeline
  K_train <- rbf_rkhs_kron(nmat)
  K_cross <- rbf_rkhs_kron_cross(ncross)
  expect_equal(nrow(K_cross), length(X_new))
  expect_equal(ncol(K_cross), length(X_train))

  # Test CFBM pipeline
  G_train <- cfbm_rkhs_kron(nmat)
  G_cross <- cfbm_rkhs_kron_cross(nmat, ncross)
  expect_equal(dim(G_train), c(2,2))
  expect_equal(dim(G_cross), c(1,2))
})
