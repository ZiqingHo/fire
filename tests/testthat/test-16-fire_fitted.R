context("fitted.fire methods")

test_that("fitted.fire_matrix returns correct structure", {
  # Create mock fire_matrix object
  mock_model <- structure(
    list(
      lambda = c(1.1),
      noise = c(0.5),
      w = list(c(0.1, 0.2)),
      converged = TRUE
    ),
    class = c("fire_matrix"),
    training_data = matrix(rnorm(6), ncol=3),
    kernels = list(cfbm),
    kernel_iprior = "cfbm",
    kernels_params = list(0.5),
    intercept = 0.1,
    dat_T = list(1:3),
    dimensions = c(2, 3),
    constant = FALSE,
    os_type = 'Apple'
  )


  # Test fitted values
  fits <- fitted(mock_model)

  # Test output structure
  expect_s3_class(fits, "fire_fitted")
  expect_named(fits, c("yhat", "rmse", "residuals", "intercept", "model"))
  expect_length(fits$yhat, 2)
  expect_true(is.numeric(fits$rmse))
})

test_that("fitted.fire_tensor handles intercept correctly", {
  # Create mock fire_tensor object
  set.seed(1)
  mock_model <- structure(
    list(
      alpha = c(0.8),
      tau = c(1.5),
      noise = c(0.2),
      w = list(c(0.1, 0.2)),
      converged = TRUE
    ),
    class = c("fire_tensor"),
    original_response = 1:2,
    training_data = list(matrix(rnorm(4),2), matrix(rnorm(4),2)),
    center = F,
    constant = T,
    kernels = list(cfbm, cfbm),
    kernel_iprior = "cfbm",
    kernels_params = list(0.5, 0.5),
    intercept = 0.2,
    dat_T = list(T1 = 1:2, T2 = 1:2),
    dimensions = c(2, 2),
    sample_id = 1,
    os_type = 'Apple'
  )


  # Test fitted values
  fits <- fitted(mock_model)

  # Test output
  expect_equal(fits$intercept, 0.2)
  expect_length(fits$residuals, 2)
})

test_that("print.fire_fitted works", {
  # Create mock fitted object
  mock_fit <- structure(
    list(
      yhat = rnorm(50),
      rmse = 0.25,
      intercept = 0.1,
      residuals = rnorm(50)
    ),
    class = "fire_fitted"
  )

  # Test print output
  expect_output(print(mock_fit), "FIRe Model Fitted Results")
  expect_output(print(mock_fit), "Training RMSE")
  expect_output(print(mock_fit), "First 6 fitted values")
})

