context("fire() function")

test_that("fire() generic works", {
  expect_true(is.function(fire))
  expect_true(inherits(fire, "function"))
})

test_that("fire.default() handles unsupported input types", {
  # Unsupported input types
  expect_error(fire(X = "character", Y = 1:3), "not supported")
  expect_error(fire(X = 1:10, Y = 1:10), "not supported")

  # Missing Y
  expect_error(fire(X = matrix(1:4, ncol=2)), "Y must be provided")

  # Invalid Y
  expect_error(fire(X = matrix(1:4, ncol=2), Y = letters[1:2]), "Y must be a numeric vector")
  expect_error(fire(X = matrix(1:4, ncol=2), Y = c(1, NA)), "NA values detected")
})

test_that("fire.data.frame() handles data frames correctly", {
  set.seed(1)
  df <- data.frame(a = runif(4,5,10), b = runif(4,5,10))
  y <- 1:4
  dat_T <- list(1:2)

  # Valid cases
  expect_message(fire(df, y, dat_T = dat_T))

  # Invalid cases
  expect_error(fire(df, y[-1]), "must match length")
  expect_error(fire(df[1:3,], y), "must match length")

  # NA handling
  df_na <- df
  df_na[1,1] <- NA
  expect_error(fire(df_na, y), "NA values detected")
})

test_that("fire.array() handles arrays correctly", {
  set.seed(1)
  arr <- array(runif(8, 5, 10), dim = c(4,2))
  y <- 1:4

  # Valid cases (first dimension matches)
  expect_message(fire(arr, y, dat_T = list(1:2)))

  # Invalid cases
  expect_error(fire(arr, 1:3), "must match length")

  # NA handling
  arr_na <- arr
  arr_na[1,1] <- NA
  expect_error(fire(arr_na, y), "NA values detected")
})

test_that("fire.list() handles lists correctly", {
  set.seed(1)
  lst <- list(matrix(runif(4,10,20), ncol = 2),
              matrix(runif(4,10,20), ncol=2),
              matrix(runif(4,10,20), ncol=2))
  y <- 1:3
  dat_T <- list(1:2, 1:2)
  fire(lst, y, dat_T = dat_T, kernels = list(cfbm, cfbm),
       kernels_params = list(0.5, 0.5),stop.eps = 1, maxiter = 3)

  # Invalid cases
  expect_error(fire(lst, 1:2), "must match length")

  # Dimension checking
  lst_bad <- list(matrix(1:4, ncol=2), matrix(1:6, ncol=3))
  expect_error(fire(lst_bad, 1:2), "identical dimensions")

  # NA handling
  lst_na <- lst
  lst_na[[1]][1,1] <- NA
  expect_error(fire(lst_na, y), "NA values detected")
})

test_that("Output structure is correct", {
  # Mock a simple case
  set.seed(1)
  df <- data.frame(a = runif(4,5,10), b = runif(4,5,10))
  y <- 1:4
  dat_T <- list(1:2)

  # Test the output class and structure
  result <- fire(df, y, dat_T = dat_T)
  expect_s3_class(result, "fire_model")
  expect_s3_class(result, "fire_matrix")
  expect_true("input_type" %in% names(attributes(result)))
  expect_true("call" %in% names(attributes(result)))

  # Check data frame conversion maintains original class info
  expect_equal(attr(result, "original_class"), "data.frame")
})

test_that("Method dispatch works correctly", {
  # Create test data
  set.seed(1)
  df <- data.frame(a = runif(4,5,10), b = runif(4,5,10))
  mat <- as.matrix(df)
  arr <- array(runif(16, 5, 10), dim = c(4,2,2))
  lst <- list(matrix(runif(4,5,10), ncol=2), matrix(runif(4,5,10), ncol=2))
  y <- 1:4
  dat_T <- list(1:2)
  # Test dispatch
  expect_s3_class(fire(df, y, dat_T), "fire_matrix")
  expect_s3_class(fire(mat, y, dat_T), "fire_matrix")
  expect_s3_class(fire(arr, y, dat_T = list(1:2, 1:2), kernels = list(cfbm,cfbm),
                       kernels_params = list(0.5, 0.5)), "fire_tensor")
  set.seed(1)
  lst <- list(matrix(runif(4,10,20), ncol = 2),
              matrix(runif(4,10,20), ncol=2),
              matrix(runif(4,10,20), ncol=2))
  expect_s3_class(fire(lst, y[1:3], dat_T = list(1:2, 1:2), kernels = list(cfbm,cfbm),
                      kernels_params = list(0.5, 0.5),stop.eps = 1,maxiter = 3), "fire_tensor")
})

test_that("Error messages are informative", {
  df <- data.frame(a = 1:4, b = 5:8)
  y <- 1:4

  # Test error messages contain relevant info
  expect_error(fire(df, y[-1]), "Number of rows in X")
  expect_error(fire(as.list(df), y), "must be arrays")
  expect_error(fire(unclass(df), y), "All elements of X must be arrays/matrices")
})
