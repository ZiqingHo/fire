context("Tensor Utility Functions")

test_that("tensor_sample() works correctly", {
  X <- array(1:24, dim = c(3, 2, 4))

  # Test sampling along mode 1
  samples_1 <- tensor_sample(X, 1)
  expect_type(samples_1, "list")
  expect_length(samples_1, 3)
  expect_named(samples_1, paste0("Sample-", 1:3))
  expect_equal(dim(samples_1[[1]]), c(2, 4))  # First sample should be 1x2x4

  # Test sampling along mode 3 (last mode)
  samples_3 <- tensor_sample(X, 3)
  expect_length(samples_3, 4)
  expect_equal(dim(samples_3[[1]]), c(3, 2))

  # Test error handling
  expect_error(tensor_sample(X, 2), "Sample mode must be either 1st mode or last mode")
  expect_error(tensor_sample("not a tensor", 1), "invalid 'length' argument")
})

test_that("vectorize_tensor() works correctly", {
  X <- array(1:8, dim = c(2, 2, 2))

  # Test with expand.grid indices
  idx_grid <- expand.grid(T1 = 1:2, T2 = 1:2, T3 = 1:2)
  vec_grid <- vectorize_tensor(X, idx_grid)
  expect_equal(vec_grid, as.vector(X))

  # Test with matrix indices
  idx_mat <- cbind(rep(1:2, each = 4), rep(c(1,1,2,2), 2), rep(1:2, 4))
  vec_mat <- vectorize_tensor(X, idx_mat)
  expect_length(vec_mat, 8)

  # Test error handling
  expect_error(vectorize_tensor(X, "invalid indices"))
})

test_that("unfolding() works correctly", {
  X <- array(1:8, dim = c(2, 2, 2))

  # Test mode-1 unfolding
  unfold1 <- unfolding(X, 1)
  expect_equal(dim(unfold1), c(2, 4))
  expect_equal(unfold1[1,], c(1, 3, 5, 7))

  # Test mode-2 unfolding
  unfold2 <- unfolding(X, 2)
  expect_equal(dim(unfold2), c(2, 4))

  # Test mode-3 unfolding
  unfold3 <- unfolding(X, 3)
  expect_equal(dim(unfold3), c(2, 4))

  # Test silent mode
  expect_silent(unfolding(X, 1, mode = FALSE))
  expect_message(unfolding(X, 1), "Mode-1 matricization")
})

test_that("dat_unfolding() works correctly", {
  X <- array(1:24, dim = c(3, 2, 4))

  # Test unfolding along first mode
  unf1 <- dat_unfolding(X, 1)
  expect_type(unf1, "list")
  expect_length(unf1, 2)  # Should unfold along remaining 2 modes
  expect_named(unf1, c("Mode-2", "Mode-3"))

  # Each mode unfolding should have same length as sample size
  expect_length(unf1[[1]], 3)
  expect_length(unf1[[2]], 3)

  # Test unfolding along last mode
  unf_last <- dat_unfolding(X, 3)
  expect_length(unf_last, 2)
  expect_named(unf_last, c("Mode-1", "Mode-2"))

  # Test message output
  expect_message(dat_unfolding(X, 1), "The 1st mode corresponds")
})

test_that("Edge cases are handled", {
  # Test empty inputs
  empty_arr <- array(numeric(0), dim = c(0, 2, 3))
  expect_error(tensor_sample(empty_arr, 1))
  expect_error(vectorize_tensor(empty_arr, matrix(1)), "cannot be empty")
})


test_that("vectorize_tensor() handles invalid indices correctly", {
  X <- array(1:8, dim = c(2, 2, 2))

  # Test single invalid index
  bad_idx <- matrix(c(1, 1, 3), nrow = 1)
  expect_error(
    vectorize_tensor(X, bad_idx),
    "Invalid indices found:\n• Dimension 3: index 3 is out of bounds (must be between 1 and 2)",
    fixed = TRUE
  )

  # Test multiple invalid indices in one row
  bad_idx <- matrix(c(0, 3, 2), nrow = 1)
  expect_error(
    vectorize_tensor(X, bad_idx),
    "Invalid indices found:\n• Dimension 1: index 0 is out of bounds (must be between 1 and 2)\n• Dimension 2: index 3 is out of bounds (must be between 1 and 2)",
    fixed = TRUE
  )

  # Test NA indices
  bad_idx <- matrix(c(1, NA, 1), nrow = 1)
  expect_error(
    vectorize_tensor(X, bad_idx),
    "• Dimension 2: index NA is out of bounds",
    fixed = TRUE
  )

  # Test invalid indices in multiple rows
  bad_idx <- matrix(c(1,1,3, 2,2,1), ncol = 3, byrow = TRUE)
  expect_error(
    vectorize_tensor(X, bad_idx),
    "• Dimension 3: index 3 is out of bounds",
    fixed = TRUE
  )
})

test_that("reshape_tensor works for 2D tensors", {
  Mat <- matrix(1:6, nrow = 3, ncol = 2)
  dim <- c(2)
  Index <- as.matrix(1:2, ncol = 1)
  N <- 3

  result <- reshape_tensor(Mat, dim, N, Index)

  # Check dimensions
  expect_equal(dim(result), c(N, dim))

  # Check values
  expect_equal(result[1,1], 1)
  expect_equal(result[1,2], 4)
  expect_equal(result[2,1], 2)
  expect_equal(result[2,2], 5)
  expect_equal(result[3,1], 3)
  expect_equal(result[3,2], 6)
})

test_that("reshape_tensor works for 3D tensors", {
  Mat <- matrix(1:8, nrow = 2)
  dim <- c(2, 2)
  Index <- expand.grid(1:2, 1:2)
  N <- 2

  result <- reshape_tensor(Mat, dim, N, Index)

  # Check dimensions
  expect_equal(dim(result), c(N, dim))

  # Check some values
  expect_equal(result[1,1,1], 1)
  expect_equal(result[1,2,2], 7)
  expect_equal(result[2,1,1], 2)
})

test_that("reshape_tensor works for higher-dimensional tensors", {
  Mat <- matrix(1:48, nrow = 3)
  dim <- c(2, 2, 2, 2)
  Index <- expand.grid(1:2, 1:2, 1:2, 1:2)
  N <- 3

  result <- reshape_tensor(Mat, dim, N, Index)

  # Check dimensions
  expect_equal(dim(result), c(N, dim))

  # Check some values
  expect_equal(result[1,1,1,1,1], 1)
  expect_equal(result[1,2,2,2,2], 46)
})

test_that("reshape_tensor handles edge cases", {
  # Empty matrix case
  Mat = matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(reshape_tensor(Mat,
                              dim = c(2,2), N = 0,
                              Index = matrix(numeric(0), ncol = 2)),
               "Mat is empty matrix")

  # Single element case
  Mat <- matrix(42, nrow = 1, ncol = 1)
  dim <- c(1)
  Index <- matrix(1, nrow = 1)
  N <- 1
  result <- reshape_tensor(Mat, dim, N, Index)
  expect_equal(result[1,1], 42)
})

test_that("reshape_tensor throws appropriate errors", {
  Mat <- matrix(1:6, nrow = 2)

  # Dimension mismatch
  expect_error(reshape_tensor(Mat, dim = c(3,3), N = 2,
                              Index = matrix(1, nrow = 2, ncol = 2)),
               "Number of columns in Mat must match number of rows in Index")

  # Index-dimension mismatch
  expect_error(reshape_tensor(Mat, dim = c(3,3,3), N = 2,
                              Index = matrix(1, nrow = 3, ncol = 2)),
               "Length of dim must match number of columns in Index")
})

test_that("valid Index matrices work", {
  # Test with correct expand.grid order
  Mat <- matrix(1:40, nrow = 2, ncol = 20)
  dim <- c(2, 2, 5)
  Index <- expand.grid(1:2, 1:2, 1:5)
  expect_silent(reshape_tensor(Mat, dim, N = 2, Index))

  # Test with different expand.grid order
  Index <- expand.grid(1:5, 1:2, 1:2)[, c(2, 3, 1)]
  expect_silent(reshape_tensor(Mat, dim, N = 2, Index))
})

test_that("invalid Index matrices are caught", {
  Mat <- matrix(1:24, nrow = 2)
  dim <- c(4, 3)

  # partially wrong indices
  Index <- expand.grid(1:4, 2:4)
  expect_error(
    reshape_tensor(Mat, dim, N = 2, Index),
    "Index matrix does not cover all combinations of dimensions:\n• Dimension 2: indices should be 1:3 but found 2:4\nExpected all combinations from expand.grid(1:4, 1:3)",
    fixed = TRUE
  )

  # Indices in wrong order
  Index <- expand.grid(1:3, 1:4)
  expect_error(
    reshape_tensor(Mat, dim, N = 2, Index),
    "Index matrix does not cover all combinations of dimensions:\n• Dimension 1: indices should be 1:4 but found 1:3\n• Dimension 2: indices should be 1:3 but found 1:4\nExpected all combinations from expand.grid(1:4, 1:3)",
    fixed = TRUE
  )
})

test_that("error messages are informative", {
  Mat <- matrix(1:24, nrow = 2)
  dim <- c(3, 2, 2)

  # Wrong order but complete coverage shouldn't error
  Index <- expand.grid(1:2, 1:2, 1:3)[, c(3,1,2)]
  expect_silent(reshape_tensor(Mat, dim, N = 2, Index))
})
