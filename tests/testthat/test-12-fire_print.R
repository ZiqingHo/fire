context("print.fire methods")

test_that("print.fire_matrix generic exists", {
  expect_true(is.function(print.fire_matrix))
})

test_that("print.fire_matrix produces correct output", {
  # Create a mock fire_matrix object
  mock_model <- structure(
    list(
      training_data = matrix(1, ncol = 10, nrow = 4),
      converged = TRUE,
      niter = 10,
      mloglik = c(-2.5, -2.3, -2.1),
      lambda = c(1.1, 1.2, 1.3),
      noise = c(0.5, 0.4, 0.3)
    ),
    class = c("fire_matrix")
  )

  # Capture output
  output <- capture.output(print.fire_matrix(mock_model))

  # Check basic structure
  expect_true(any(grepl("FIRe Model \\(Matrix Input\\)", output)))
  expect_true(any(grepl("Converged in 10 iterations", output)))
  expect_true(any(grepl("Dimensions: 4 × 10", output)))
  expect_true(any(grepl("mLogLik: -2.1", output)))
  expect_true(any(grepl("Lambda: 1.300", output)))
  expect_true(any(grepl("Noise: 0.300", output)))
})

test_that("print.fire_tensor produces correct output", {
  # Create a mock fire_tensor object
  mock_model <- structure(
    list(
      converged = FALSE,
      niter = 15,
      mloglik = c(-5.1, -4.9, -4.8),
      alpha = c(0.8, 0.9, 1.0),
      tau = c(1.5, 1.4, 1.3),
      noise = c(0.2, 0.18, 0.15),
      dimensions = c(50, 4, 3)
    ),
    class = c("fire_tensor")
  )

  # Capture output
  output <- capture.output(print.fire_tensor(mock_model))

  # Check basic structure
  expect_true(any(grepl("FIRe Model \\(Tensor Input\\)", output)))
  expect_true(any(grepl("Stopped after 15 iterations", output)))
  expect_true(any(grepl("Dimensions: 50 × 4 × 3", output)))
  expect_true(any(grepl("mLogLik: -4.8", output)))
  expect_true(any(grepl("Alpha: 1.000", output)))
  expect_true(any(grepl("Tau: 1.300", output)))
  expect_true(any(grepl("Noise: 0.150", output)))
})

test_that("print methods handle edge cases", {
  # Test empty model
  empty_model <- structure(
    list(
      converged = FALSE,
      niter = 0,
      mloglik = numeric(0),
      lambda = numeric(0),
      noise = numeric(0),
      dimensions = c(0L, 0L)
    ),
    class = c("fire_matrix")
  )

  # Capture output
  output <- capture.output(print.fire_matrix(empty_model))

  # Check basic structure
  expect_true(any(grepl("FIRe Model \\(Matrix Input\\)", output)))
  expect_true(any(grepl("Stopped after 0 iterations", output)))
  expect_true(any(grepl("Dimensions: ", output)))
  expect_true(any(grepl("mLogLik: ", output)))
  expect_true(any(grepl("Estimated hyperparameters:", output)))
})

test_that("print methods return invisible x", {
  mock_model <- structure(
    list(
      converged = TRUE,
      niter = 5,
      mloglik = -1.2,
      lambda = 1.0,
      noise = 0.5,
      dimensions = c(10L, 2L)
    ),
    class = c("fire_matrix")
  )

  # Test return value
  out <- print(mock_model)
  expect_identical(out, mock_model)

})
