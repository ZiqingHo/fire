context("summary.fire methods")

test_that("summary.fire_matrix works with valid input", {
  # Create mock fire_matrix object
  mock_model <- structure(
    list(
      lambda = c(1.1, 1.2, 1.3),
      noise = c(0.5, 0.4, 0.3),
      converged = TRUE,
      niter = 10,
      mloglik = c(-2.5, -2.3, -2.1),
      duration = structure(5.2, units = "secs", class = "difftime"),
      training_rmse = 0.25
    ),
    class = c("fire_model", "fire_matrix"),
    dimensions = c(100L, 5L),
    kernels = list(cfbm),
    kernel_iprior = "cfbm",
    kernels_params = list(0.5),
    intercept = 0.1,
    convergence = list(
      maxiter = 200,
      tolerance = 1e-5,
      final_change = 1e-6
    ),
    training_data = matrix(rnorm(500), ncol=5),
    call = quote(fire(X = x, Y = y, kernels = cfbm))
  )

  # Generate summary
  s <- summary(mock_model)

  # Test class and structure
  expect_s3_class(s, "summary.fire_matrix")
  expect_named(s, c("coefficients", "convergence", "model_fit",
                    "timing", "dimensions", "kernel_info"))

  # Test print method
  expect_output(print(s), "FIRe Model Summary")
  expect_output(print(s), "Converged in 10 iterations")
})

test_that("summary.fire_tensor works with valid input", {
  # Create mock fire_tensor object
  mock_model <- structure(
    list(
      alpha = c(0.8, 0.9, 1.0),
      tau = c(1.5, 1.4, 1.3),
      noise = c(0.2, 0.18, 0.15),
      converged = FALSE,
      niter = 15,
      mloglik = c(-5.1, -4.9, -4.8),
      duration = structure(8.5, units = "secs", class = "difftime")
    ),
    class = c("fire_model", "fire_tensor"),
    dimensions = list(1:50, 1:4, 1:3),
    input_type = "array",
    kernels = list(cfbm, cfbm),
    kernel_iprior = "cfbm",
    kernels_params = list(0.5, 0.5),
    iprior_param = 0.5,
    intercept = 0.2,
    convergence = list(
      maxiter = 200,
      tolerance = 1e-5,
      final_change = 1e-4
    ),
    original_response = rnorm(50),
    call = quote(fire(X = x, Y = y, kernels = list(cfbm, cfbm)))
  )

    # Generate summary
    s <- summary(mock_model)

    # Test class and structure
    expect_s3_class(s, "summary.fire_tensor")
    expect_named(s, c("coefficients", "convergence", "model_fit",
                      "timing", "dimensions", "input_type",
                      "sample_id", "sample_size", "kernel_info"))

    # Test print method
    expect_output(print(s), "FIRe Tensor Model Summary")
    expect_output(print(s), "Failed to converge")
})

  test_that("summary handles edge cases", {
    # Empty model
    empty_model <- structure(
      list(
        lambda = numeric(0),
        noise = numeric(0),
        converged = FALSE,
        niter = 0,
        mloglik = numeric(0)
      ),
      class = c("fire_model", "fire_matrix"),
      dimensions = c(0L, 0L),
      call = quote(fire(X = x, Y = y))
    )

      expect_silent(summary(empty_model))
  })
