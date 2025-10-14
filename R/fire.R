#' Fit a FIRE Model
#'
#' A function to perform functional regression using I-priors with Reproducing Kernel Hilbert Space (RKHS) norm.
#' The FIRE model parameters are estimated using an EM algorithm.
#'
#' @param X A numeric input in \code{matrix}, \code{data.frame}, \code{array} or \code{list}
#' @param Y A numeric response vector
#' @param dat_T List of index sets for each mode
#' @param kernels List of kernel functions for each mode (see \code{\link{kernels_fire}})
#' @param kernels_params List of parameters for each kernel
#' @param kernel_iprior Type of I-prior kernel
#' @param iprior_param Parameter for I-prior kernel:
#' \itemize{
#'   \item \code{"cfbm"} - Hurst coefficient (default 0.5)
#'   \item \code{"rbf"} - lengthscale (default 1)
#'   \item \code{"linear"} - offset (default 0)
#'   \item \code{"poly"} - degree and offset (default c(2, mean(Y)))
#' }
#' @param control A list of control parameters (see Details)
#' @param ... Additional arguments passed to methods
#'
#' @details
#' The \code{control} argument can include the following parameters:
#' \itemize{
#'   \item{\code{scale}: Logical indicating whether to center the response (default TRUE)}
#'   \item{\code{maxiter}: Maximum number of EM iterations (default 200)}
#'   \item{\code{stop.eps}: Convergence tolerance (default 1e-3)}
#'   \item{\code{constant_g}: Logical indicating whether to include constant kernel term in g (default TRUE)}
#'   \item{\code{constant_h}: Logical indicating whether to include constant kernel term in h (default FALSE)}
#'   \item{\code{center}: Logical indicating whether to center the kernel g (default FALSE)}
#'   \item{\code{std}: Logical indicating whether to standardise the kernel g (default TRUE)}
#'   \item{\code{par_init}: Optional list of initial parameter values (lambda, noise)}
#'   \item{\code{os_type}: Operating system type for compatibility ("Apple" or "Windows", default "Apple")}
#'   \item{\code{cores}: Number of cores for parallel computation (default: detectCores() - 1)}
#'   \item{\code{asymptote}: Logical to use asymptotic initial values (default TRUE)}
#'   \item{\code{sample_id}: Which mode contains samples (default 1)}
#'   \item{\code{epsilon}: Small positive constant in initialisation of EM algorithm (default 1)}
#' }
#'
#' @section Methods:
#' This generic function has methods for different input types:
#' \itemize{
#'   \item \code{\link{fire.matrix}} for \code{matrix} or \code{data.frame} inputs
#'   \item \code{\link{fire.tensor}} for \code{array} or \code{list} inputs
#' }
#'
#' @return An \code{fire_matrix} or \code{fire_tensor} object. The \code{print()} and \code{summary()} methods display the model information.
#'
#' @seealso \code{\link{kernels_fire}}, \code{\link{Manure}}, \code{\link{Housing}}
#'
#' @examples
#' # Matrix input
#' data(Manure)
#' mod1 <- fire(X = Manure$absorp[1:5,], Y = Manure$y$DM[1:5],
#'  dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4))
#' summary(mod1)
#'
#' # Array input
#' data(Housing)
#' dat_T <- list(T1 = 1:4, T2 = 1:9)
#' mod2 <- fire(X = Housing$X[1:5,,], Y = Housing$y[1:5,2],
#' kernels = list(kronecker_delta, kronecker_delta),
#' kernels_params = list(NA, NA),
#' dat_T = dat_T, control = list(stop.eps = 2, maxiter = 4))
#' summary(mod2)
#' @export
fire <- function(X, Y, dat_T, kernels, kernels_params,
                 kernel_iprior = 'cfbm', iprior_param = NULL,
                 control = list(), ...) {
  if (missing(Y)) stop("Response vector Y must be provided")
  if (!is.numeric(Y)) stop("Y must be a numeric vector")
  if (anyNA(Y)) stop("NA values detected in Y")

  UseMethod("fire")
}

#' @name fire
#' @export
fire.default <- function(X, Y, dat_T, kernels, kernels_params,
                         kernel_iprior = 'cfbm', iprior_param = NULL,
                         control = list(), ...) {
  supported <- c("matrix", "array", "list", "data.frame")
  stop(
    "Input type '", class(X)[1], "' not supported.\n",
    "Use one of: ", paste(supported, collapse = ", ")
  )
}

#' @rdname fire
#' @export
fire.data.frame <- function(X, Y, dat_T, kernels = list(cfbm), kernels_params = list(0.5),
                            kernel_iprior = 'cfbm', iprior_param = NULL,
                            control = list(), ...) {
  # Data frame validation
  if (!is.data.frame(X)) stop("X must be a data frame")
  if (nrow(X) != length(Y)) {
    stop("Number of rows in X (", nrow(X),
         ") must match length of Y (", length(Y), ")")
  }
  if (anyNA(X)) stop("NA values detected in X")

  # Convert and delegate to matrix method
  X_mat <- as.matrix(X)
  result <- fire.matrix(X_mat, Y, dat_T, kernels, kernels_params,
                        kernel_iprior, iprior_param,
                        control, ...)

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
fire.array <- function(X, Y, dat_T, kernels, kernels_params,
                       kernel_iprior = 'cfbm', iprior_param = NULL,
                       control = list(), ...) {
  # Array validation
  if (!is.array(X)) stop("X must be an array")
  if (!(dim(X)[1] == length(Y) || tail(dim(X), 1) == length(Y))) {
    stop("Either first dimension (", dim(X)[1], ") or last dimension (",
         tail(dim(X), 1), ") of X must match length of Y (", length(Y), ")")
  }
  if (anyNA(X)) stop("NA values detected in X")

  # Convert to tensor and dispatch
  class(X) <- c("tensor", class(X))
  result <- fire.tensor(X, Y, dat_T, kernels, kernels_params,
                        kernel_iprior, iprior_param,
                        control, ...)

  structure(
    result,
    class = c("fire_tensor"),
    input_type = "array",
    call = match.call()
  )
}

#' @rdname fire
#' @export
fire.list <- function(X, Y, dat_T, kernels, kernels_params,
                      kernel_iprior = 'cfbm', iprior_param = NULL,
                      control = list(), ...) {
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
  result <- fire.tensor(X, Y, dat_T, kernels, kernels_params,
                        kernel_iprior, iprior_param,
                        control,...)

  structure(
    result,
    class = c("fire_tensor"),
    input_type = "list",
    call = match.call()
  )
}
