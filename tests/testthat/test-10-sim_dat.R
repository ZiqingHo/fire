context("sim_dat function tests")

test_that("Basic functionality works for 1D case", {
  set.seed(1)
  result <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    N = 10,
    Ntrain = 5,
    constant = F
  )

  expect_type(result, "list")
  expect_named(result, c("X", "y"))
  expect_equal(dim(result$X), c(10, 3))
  expect_length(result$y, 10)
})

test_that("Basic functionality works for 2D case", {
  set.seed(1)
  result <- sim_dat(
    kernels = list(cfbm, rbf),
    kernels_param = list(0.5, 1),
    alpha = c(0.5, 0.5),
    dat_T = list(1:5, 1:2),
    N = 10,
    Ntrain = 5
  )

  expect_type(result, "list")
  expect_named(result, c("X", "y"))
  expect_equal(dim(result$X), c(10, 5, 2))
  expect_length(result$y, 10)
})

test_that("Basic functionality works for 3D case", {
  set.seed(123)
  result <- sim_dat(
    kernels = list(cfbm, cfbm, cfbm),
    kernels_param = list(0.5, 0.5, 0.5),
    alpha = c(0.5, 0.5, 0.5),
    dat_T = list(1:3, 1:2, 1:4),
    N = 10,
    Ntrain = 5
  )

  expect_type(result, "list")
  expect_named(result, c("X", "y"))
  expect_equal(dim(result$X), c(10, 3, 2, 4))
  expect_length(result$y, 10)
})

test_that("Different kernel types work", {
  set.seed(1)
  result_cfbm <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    kernel_iprior = "cfbm",
    constant = F,
    N = 5,
    Ntrain = 3
  )

  result_rbf <- sim_dat(
    kernels = list(rbf),
    kernels_param = list(1),
    alpha = c(0.5),
    dat_T = list(1:3),
    kernel_iprior = "rbf",
    constant = F,
    N = 5,
    Ntrain = 3
  )

  result_linear <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    kernel_iprior = "linear", constant = F,
    N = 5,
    Ntrain = 3
  )

  expect_true(all(result_cfbm$y != result_rbf$y))
  expect_true(all(result_cfbm$y != result_linear$y))
})

test_that("Parameter validation works", {
  # Mismatched lengths
  expect_error(sim_dat(
    kernels = list(cfbm, cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5, 0.5),
    dat_T = list(1:3, 1:2),
    N = 10,
    Ntrain = 5
  ), "kernels and kernels_param must have same length")

  # Ntrain >= N
  expect_error(sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    N = 5,
    Ntrain = 5
  ), "Ntrain must be smaller than N")
})

test_that("Non-integer time points work", {
  set.seed(123)
  result <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(seq(0, 1, length.out = 3)),
    constant = F,
    N = 5,
    Ntrain = 3
  )

  expect_equal(dim(result$X), c(5, 3))
})

test_that("Different parameter combinations work", {
  set.seed(123)
  # With intercept
  result_int <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    intercept_y = 1,
    intercept_x = 1,
    constant = F,
    N = 5,
    Ntrain = 3
  )

  # Without constant term
  result_noconst <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    constant = FALSE,
    N = 5,
    Ntrain = 3
  )

  expect_true(mean(result_int$y) > mean(result_noconst$y))
})

test_that("Different noise parameters affect results", {
  set.seed(1)
  result_lownoise <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    sigma = 0.1,
    constant = F,
    N = 5,
    Ntrain = 3
  )

  result_highnoise <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:3),
    sigma = 10,
    constant = F,
    N = 5,
    Ntrain = 3
  )

  expect_lt(sd(result_lownoise$y), sd(result_highnoise$y))
})

test_that("Output dimensions are correct for various cases", {
  set.seed(1)
  # 1D case
  result1d <- sim_dat(
    kernels = list(cfbm),
    kernels_param = list(0.5),
    alpha = c(0.5),
    dat_T = list(1:4),
    constant = F,
    N = 6,
    Ntrain = 4
  )
  expect_equal(dim(result1d$X), c(6, 4))

  # 2D case
  result2d <- sim_dat(
    kernels = list(cfbm, rbf),
    kernels_param = list(0.5, 1),
    alpha = c(0.5, 0.5),
    dat_T = list(1:3, 1:2),
    N = 6,
    Ntrain = 4
  )
  expect_equal(dim(result2d$X), c(6, 3, 2))

  # 4D case
  result4d <- sim_dat(
    kernels = list(cfbm, cfbm, cfbm, cfbm),
    kernels_param = list(0.5, 0.5, 0.5, 0.5),
    alpha = c(0.5, 0.5, 0.5, 0.5),
    dat_T = list(1:2, 1:2, 1:2, 1:2),
    N = 6,
    Ntrain = 4
  )
  expect_equal(dim(result4d$X), c(6, 2, 2, 2, 2))
})
