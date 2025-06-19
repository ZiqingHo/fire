context("predict.fire methods")

test_that("predict.fire_matrix returns correct structure", {
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

  # Test prediction
  newdata <- matrix(rnorm(6), ncol=3)
  pred <- predict(mock_model, newdata)

  # Test output structure
  expect_s3_class(pred, "fire_prediction")
  expect_named(pred, c("yhat", "y.actual", "newdata", "model", "test_metrics"))
  expect_length(pred$yhat, 2)
})

test_that("predict.fire_tensor handles test data", {
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

  # Test prediction with test data
  newdata <- list(rnorm(4,2))
  pred <- predict(mock_model, newdata)

  # Test output
  expect_s3_class(pred, "fire_prediction")
  expect_equal(length(pred$yhat), 1)
  expect_true(is.null(pred$test_metrics))
})

test_that("print.fire_prediction works", {
  # Create mock prediction
  mock_pred <- structure(
    list(
      yhat = rnorm(50),
      y.actual = rnorm(50),
      test_metrics = list(rmse = 0.4)
    ),
    class = "fire_prediction"
  )

  # Test print output
  expect_output(print(mock_pred), "FIRe Model Predictions")
  expect_output(print(mock_pred), "Test RMSE")
})

test_that("predict handles errors", {
  mock_model <- structure(
    list(
      lambda = 1.0,
      noise = 0.5,
      w = list(c(0.1, 0.2)),
      converged = TRUE
    ),
    class = c("fire_model", "fire_matrix"),
    training_data = matrix(rnorm(20), ncol=2),
    dimensions = c(10, 2)
  )

  # Wrong dimensions
  expect_error(predict(mock_model, matrix(rnorm(9), ncol=3)),
               "columns matching training data")

  # Wrong Ynew length
  expect_error(predict(mock_model, matrix(rnorm(4), ncol=2), Ynew = 1:3),
               "must match number of rows")
})
