context("Internal Kronecker Functions")

test_that("Kronecker_inv_helper works for all cases", {
  # Setup test data
  G1 <- list(matrix(c(2,1,1,2), nrow=2))
  G2 <- list(matrix(c(2,1,1,2), nrow=2), matrix(c(3,1,1,3), nrow=2))

  # Test m=1 without constant
  res1 <- Kronecker_inv_helper(G1, m=1, alpha=0.5, constant=FALSE)
  expect_true("G1" %in% names(res1))
  expect_null(res1$Q1)

  # Test m=1 with constant
  res2 <- Kronecker_inv_helper(G1, alpha=0.5, constant=TRUE) # test default m
  expect_true(all(c("Q1","L1") %in% names(res2)))
  expect_equal(ncol(res2$Q1), 2)

  # Test m=2
  res3 <- Kronecker_inv_helper(G2, m=2, alpha=c(0.5,0.5), constant=TRUE)
  expect_true(all(c("Q1","L1","Q2","L2") %in% names(res3)))
})

test_that("Lambda_inv computes correct inverses", {
    # 1D case
    L1 <- list(c(1,2))
    inv1 <- Lambda_inv(L1)
    expect_equal(inv1, c(1, 0.5))

    # 2D case
    L2 <- list(c(1,2), c(2,4))
    inv2 <- Lambda_inv(L2)
    expected <- c(0.5, 0.25, 0.25, 0.125)
    expect_equal(inv2, expected)

    # Edge cases
    expect_warning(Lambda_inv(list(c(1,0), c(1,1))), "At least one eigenvalue is 0.")
})

test_that("kron_mv performs correct multiplications", {
    Q <- matrix(c(1,0,0,1), nrow=2)
    X <- c(1,2,3,4)

    # Basic functionality
    res1 <- kron_mv(X, Q)
    expect_equal(res1, c(1,3,2,4))
})

test_that("Functions work together", {
    # Integrated test
    G <- list(matrix(c(2,1,1,2), nrow=2))
    decomp <- Kronecker_inv_helper(G, alpha=0.5, constant=TRUE)
    inv <- Lambda_inv(list(decomp$L1))
    res <- kron_mv(c(1,2), decomp$Q1)

    expect_equal(length(res), 2)
    expect_true(is.numeric(res))
})
