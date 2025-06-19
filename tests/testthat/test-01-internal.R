context("Internal utility functions")

test_that("suffix() generates correct ordinal suffixes", {
  expect_equal(suffix(1), "1st")
  expect_equal(suffix(2), "2nd")
  expect_equal(suffix(3), "3rd")
  expect_equal(suffix(4), "4th")
  expect_equal(suffix(11), "11th")  # Special case
  expect_equal(suffix(21), "21st")
  expect_equal(suffix(101), "101st")
})

test_that("suffix() handles edge cases", {
  expect_error(suffix("text"))  # Non-numeric input
  expect_error(suffix(NA))      # Missing value
})

test_that("parse_index() handles tensor subsetting", {
  X <- array(1:8, dim = c(2, 2, 2))

  expect_equal(parse_index(X, "1,1,1"), 1)
  expect_equal(parse_index(X, "2,1,2"), 6)
  expect_equal(dim(parse_index(X, ",,")), c(2, 2, 2))  # All empty indices
})

test_that("parse_index() throws errors for invalid input", {
  X <- array(1:8, dim = c(2, 2, 2))
  # Non-numeric indices
  expect_error(parse_index(X, "1,2,foo"), "Invalid index: 'foo'")
  expect_error(parse_index(X, "a,b,c"), "Invalid index: 'a'")

  # Mixed valid/invalid
  expect_error(parse_index(X, "1,2,3x"), "Invalid index: '3x'")
})

test_that("Basic functionality works", {
  # Test with two variables
  input1 <- list(c(1, 2), c(3, 4, 5))
  expected1 <- matrix(c(
    1, 1,
    1, 2,
    1, 3,
    2, 1,
    2, 2,
    2, 3
  ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("T1", "T2")))

  expect_equal(create_index_matrix(input1), expected1)

  # Test with non-integer time points
  input2 <- list(seq(0, 1, length.out = 2), seq(0, 1, length.out = 3))
  expected2 <- matrix(c(
    1, 1,
    1, 2,
    1, 3,
    2, 1,
    2, 2,
    2, 3
  ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("T1", "T2")))

  result2 <- create_index_matrix(input2)
  expect_equal(result2, expected2)
  expect_true(is.integer(result2))
})

test_that("Edge cases are handled", {
  # Empty input
  expect_error(create_index_matrix(list()))

  # Single variable with no time points
  expect_equal(nrow(create_index_matrix(list(numeric(0)))), 0)

  # Multiple variables with some having no time points
  expect_equal(nrow(create_index_matrix(list(c(1,2), numeric(0)))), 0)

  # Single time point
  expect_identical(
    create_index_matrix(list(1.5, 2.7)),
    matrix(c(1L, 1L), ncol = 2, dimnames = list(NULL, c("T1", "T2")))
  )

})

test_that("Column names are correct", {
  input <- list(c(1,2), c(3,4), c(5,6))
  result <- create_index_matrix(input)
  expect_equal(colnames(result), c("T1", "T2", "T3"))

  # Test with many columns
  many_cols <- replicate(10, c(1,2), simplify = FALSE)
  result_many <- create_index_matrix(many_cols)
  expect_equal(colnames(result_many), paste0("T", 1:10))
})

test_that("Order of combinations is correct", {
  input <- list(c(1,2), c(3,4))
  result <- create_index_matrix(input)
  # First column should vary fastest
  expect_equal(result[[1, "T1"]], 1)
  expect_equal(result[[2, "T1"]], 1)
  expect_equal(result[[3, "T1"]], 2)
  expect_equal(result[[4, "T1"]], 2)
})

test_that("Output is always integer", {
  # Non-integer inputs
  expect_true(is.integer(create_index_matrix(list(c(1.1, 2.2)))))
  expect_true(is.integer(create_index_matrix(list(c(0.1, 0.2, 0.3), c(0.4, 0.5)))))

  # Mixed integer and non-integer
  expect_true(is.integer(create_index_matrix(list(c(1, 2), c(3.1, 4.2, 5.3)))))
})
