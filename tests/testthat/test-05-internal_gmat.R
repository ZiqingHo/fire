context("gmat - Gram Matrix Construction")

test_that("Basic functionality works for supported kernels", {
  # Setup test data
  set.seed(123)
  dat1 <- matrix(rnorm(4), ncol=2)
  dat2 <- matrix(rnorm(6), ncol=3)

  # Test cfbm kernel
  expect_silent({
    G_cfbm <- gmat(kernels = list(cfbm),
                   kernels_params = list(0.5),
                   dat = list(dat1))
  })
  expect_type(G_cfbm, "list")
  expect_equal(length(G_cfbm), 1)
  expect_true(is.matrix(G_cfbm[[1]]))

  # Test rbf kernel
  expect_silent({
    G_rbf <- gmat(kernels = list(rbf),
                  kernels_params = list(1.0),
                  dat = list(dat1))
  })
  expect_true(all(diag(G_rbf[[1]]) > 0)) # Diagonal should be positive

  # Test polynomial kernel
  expect_silent({
    G_poly <- gmat(kernels = list(polynomial),
                   kernels_params = list(c(2, 1)), # d=2, offset=1
                   dat = list(dat1))
  })
})

test_that("Parameter validation works correctly", {
  dat <- list(matrix(rnorm(4), ncol=2))

  # Invalid kernel type
  expect_error(
    gmat(kernels = list(function(x) x),
         kernels_params = list(NA),
         dat = dat),
    "Unsupported kernel function"
  )

  # Missing parameters
  expect_error(
    gmat(kernels = list(rbf),
         kernels_params = list(NULL),
         dat = dat),
    "All elements of kernels_params must be non-NULL."
  )

  # Mismatched lengths
  expect_error(
    gmat(kernels = list(rbf, cfbm),
         kernels_params = list(1.0),
         dat = list(dat[[1]], dat[[1]])),
    "length"
  )
})

test_that("Centering works as expected", {
  dat <- list(matrix(rnorm(6), ncol=2))

  G_uncentered <- gmat(kernels = list(rbf),
                       kernels_params = list(1.0),
                       dat = dat,
                       center = FALSE)

  G_centered <- gmat(kernels = list(rbf),
                     kernels_params = list(1.0),
                     dat = dat,
                     center = TRUE)

  # Centered matrix should have smaller values
  expect_true(mean(abs(G_uncentered[[1]])) > mean(abs(G_centered[[1]])))

              # Centered matrix rows should sum to ~0
  expect_lt(max(abs(rowSums(G_centered[[1]]))), 1e-10)
})

test_that("All supported kernel types work", {
    dat <- list(1:4)

    kernels <- list(
                  cfbm,
                  fbm,
                  cfbm_sd,
                  rbf,
                  kronecker_delta,
                  polynomial,
                  mercer
                )

    params <- list(
                  0.5,                # cfbm Hurst
                  0.5,                # fbm Hurst
                  0.5,                # cfbm_sd Hurst
                  1.0,                # rbf lengthscale
                  NA,                 # kronecker_delta (no params)
                  c(2, 1),            # polynomial (d = 2, offset = 1)
                  0.1                 # mercer delta
                )

    for (i in seq_along(kernels)) {
      # For mercer specifically, don't test for silence
      if (identical(kernels[[i]], mercer)) {
        G <- gmat(kernels = list(kernels[[i]]),
                  kernels_params = list(params[[i]]),
                  dat = dat)
      } else {
        # For all other kernels, verify they run silently
        expect_silent({
          G <- gmat(kernels = list(kernels[[i]]),
                    kernels_params = list(params[[i]]),
                    dat = dat)
        })
      }
      expect_true(is.matrix(G[[1]]))
      expect_equal(dim(G[[1]]), rep(4, 2))
    }
})

test_that("Multiple kernels work together", {
      dat <- list(
        matrix(rnorm(6), ncol=2),  # Mode 1
        matrix(rnorm(9), ncol=3)   # Mode 2
      )

      G_multi <- gmat(
        kernels = list(rbf, cfbm),
        kernels_params = list(1.0, 0.5),
        dat = dat
      )

      expect_length(G_multi, 2)
      expect_equal(dim(G_multi[[1]]), c(3, 3))  # First kernel on 3 samples
      expect_equal(dim(G_multi[[2]]), c(3, 3))  # Second kernel on same 3 samples
})
