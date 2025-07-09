#' Fit a FIRE model
#'
#' A function to perform functional regression using I-priors and Reproducing Kernel Hilbert Space (RKHS).
#' The FIRE model parameters are estimated by using EM algorithm.
#'
#' @details
#' The \code{fire} function is able to take matrix, data.frame, array and list inputs,
#' the syntax is as per the default S3 method.
#'
#' @param sample_id The sample mode identifier, either the 1st or last mode
#' @param X A numeric inputs in matrix, data.frame, array or list form
#' @param ... Additional arguments passed to methods
#' @return An object of class \code{fire_matrix} or \code{fire_tensor}.
#' The \code{print()} and \code{summary()} methods show the corresponding model information.
#'
#' @section Methods:
#' This generic function has methods for different input types:
#' \itemize{
#'   \item \code{\link{fire.matrix}} for matrix/data.frame inputs
#'   \item \code{\link{fire.tensor}} for array/list inputs
#' }
#'
#' @seealso \code{\link{kernels_fire}}, \code{\link{Manure}}, \code{\link{Housing}}
#'
#' @examples
#' # Matrix input
#' data(Manure)
#' mod1 <- fire(X = Manure$absorp[1:5,], Y = Manure$y$DM[1:5],
#'  dat_T = list(1:700), stop.eps = 2, maxiter = 4)
#' summary(mod1)
#'
#' # Array input
#' data(Housing)
#' dat_T <- list(T1 = 1:4, T2 = 1:9)
#' mod2 <- fire(X = Housing$X[1:5,,], Y = Housing$y[1:5,2],
#' kernels = list(kronecker_delta, kronecker_delta),
#' kernels_params = list(NA, NA),
#' dat_T = dat_T, stop.eps = 2, maxiter = 4)
#' summary(mod2)
#' @export
fire <- function(X, Y, ...) {
  if (missing(Y)) stop("Response vector Y must be provided")
  if (!is.numeric(Y)) stop("Y must be a numeric vector")
  if (anyNA(Y)) stop("NA values detected in Y")

  UseMethod("fire")
}

#' @name fire
#' @export
fire.default <- function(X, Y,...,
                         dat_T, scale = TRUE,
                         kernels, kernels_params,
                         kernel_iprior = 'cfbm', iprior_param = NULL,
                         maxiter = 200, stop.eps = 1e-5, center = FALSE,
                         par_init = NULL, os_type = "Apple", asymptote = TRUE, sample_id = 1
                         ) {


  supported <- c("matrix", "array", "list", "data.frame")
  stop(
    "Input type '", class(X)[1], "' not supported. ",
    "Use one of: ", paste(supported, collapse = ", ")
  )
}

#' @rdname fire
#' @export
fire.data.frame <- function(X, Y, ...) {
  # Data frame validation
  if (!is.data.frame(X)) stop("X must be a data frame")
  if (nrow(X) != length(Y)) {
    stop("Number of rows in X (", nrow(X),
         ") must match length of Y (", length(Y), ")")
  }
  if (anyNA(X)) stop("NA values detected in X")

  # Convert and delegate to matrix method
  X_mat <- as.matrix(X)
  result <- fire.matrix(X_mat, Y, ...)

  structure(
    result,
    class = c("fire_matrix"),
    input_type = "data.frame",
    original_class = "data.frame",
    dimensions = dim(X_mat),
    call = match.call()
  )
}

#' @rdname fire
#' @export
fire.array <- function(X, Y, ...) {
  # Array validation
  if (!is.array(X)) stop("X must be an array")
  if (!(dim(X)[1] == length(Y) || tail(dim(X), 1) == length(Y))) {
    stop("Either first dimension (", dim(X)[1], ") or last dimension (",
         tail(dim(X), 1), ") of X must match length of Y (", length(Y), ")")
  }
  if (anyNA(X)) stop("NA values detected in X")

  # Convert to tensor and dispatch
  class(X) <- c("tensor", class(X))
  result <- fire.tensor(X, Y, ...)

  structure(
    result,
    class = c("fire_tensor"),
    input_type = "array",
    call = match.call()
  )
}

#' @rdname fire
#' @export
fire.list <- function(X, Y, ...) {
  # List validation
  if (!is.list(X)) stop("X must be a list")
  if (!all(sapply(X, is.array))) {
    stop("All elements of X must be arrays/matrices")
  }
  if (length(X) != length(Y)) {
    stop("Length of X (", length(X),
         ") must match length of Y (", length(Y), ")")
  }
  if (any(sapply(X, anyNA))) stop("NA values detected in list elements of X")

  dims <- sapply(X, dim)
  if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
    stop("All elements in X must have identical dimensions")
  }

  # Convert to tensor and dispatch
  class(X) <- c("tensor", class(X))
  result <- fire.tensor(X, Y, ...)

  structure(
    result,
    class = c("fire_tensor"),
    input_type = "list",
    call = match.call()
  )
}
