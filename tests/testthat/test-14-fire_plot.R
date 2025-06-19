context("plot.fire methods")

test_that("plot.fire_fitted produces ggplot output", {
  # Create mock fitted object
  mock_fit <- structure(
    list(
      yhat = rnorm(100),
      residuals = rnorm(100),
      rmse = 0.5
    ),
    class = "fire_fitted"
  )

  # Test plot output
  p <- plot(mock_fit)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Residuals vs Fitted")
  expect_true("Fitted" %in% p$labels$x)
})

test_that("plot.fire_prediction with test data produces grid output", {
  # Mock prediction with test data
  mock_pred <- structure(
    list(
      yhat = rnorm(50),
      y.actual = rnorm(50),
      test_metrics = list(rmse = 0.4)
    ),
    class = "fire_prediction"
  )

  # Test plot output
  p <- plot(mock_pred)
  expect_true(inherits(p, "gtable")) # gridExtra returns gtable
})

test_that("plot.fire_prediction without test data produces ggplot", {
  # Mock prediction without test data
  mock_pred <- structure(
    list(
      yhat = rnorm(50)
    ),
    class = "fire_prediction"
  )

  # Test plot output and warning
  expect_warning(p <- plot(mock_pred), "No test values provided")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Distribution of Predicted Values")
})
