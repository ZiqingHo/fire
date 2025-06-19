context("pre_initial - Initial Parameter Estimation")

test_that("Basic functionality works for different kernel types", {
  # Setup test data
  X <- matrix(runif(10, 5, 10), ncol=2)
  Y <- runif(5, 5, 10)
  dat_T <- list(1:2)
  Index <- matrix(1:5, ncol = 1)

  # Test CFBM kernel
  expect_silent({
    init_cfbm <- pre_initial(X, Y, dat_T,
                             kernels = list(cfbm),
                             kernels_params = list(0.5),
                             Index = Index, constant = F,
                             kernel_iprior = "cfbm")
  })
  expect_length(init_cfbm$lambda, 2)
  expect_length(init_cfbm$sigma, 2)

  # Test RBF kernel
  expect_silent({
    init_rbf <- pre_initial(X, Y, dat_T,
                            kernels = list(rbf),
                            kernels_params = list(1),
                            Index = Index, constant = F,
                            kernel_iprior = "rbf")
  })
  expect_true(all(init_rbf$lambda > 0))
})


