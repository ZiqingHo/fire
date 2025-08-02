#' @title Construct Kernel Gram Matrices
#' @name internal_gmat
#' @description Internal function to construct Gram matrices for different modes and kernels.
#'   Used internally by \code{\link{fire}}. Not exported or intended for direct use.
#' @keywords internal
NULL

#' @rdname internal_gmat
#'
#' @description Constructs Gram matrices for each mode using specified kernels and parameters.
#'   Handles multiple kernel types including fractional Brownian motion (fbm),
#'   centered fractional Brownian motion (cfbm), radial basis function (rbf),
#'   Kronecker delta, polynomial, and Mercer-like kernels.
#'
#' @param kernels List of kernel functions. The available kernels can be found in \code{\link{kernels_fire}}
#' @param kernels_params List of parameters for each kernel:
#' \itemize{
#'   \item{For \code{cfbm}/\code{fbm}: Hurst parameter (numeric)}
#'   \item{For \code{rbf}: lengthscale (numeric)}
#'   \item{For \code{polynomial}: degree `d` and offset (numeric vector length 2)}
#'   \item{For \code{kronecker_delta}: none (use NA)}
#'   \item{For \code{mercer}: delta parameter (numeric)}
#'   \item{For \code{matern}: nu, lengthsalce and sigma (numeric vector length 3)}
#'   }
#' @param dat List of data corresponding to each mode (one element per kernel)
#' @param center Logical indicating whether to center the kernel matrix
#' @param std Logical indicating whether to standardise the kernel matrix
#'
#' @return List of Gram matrices corresponding to each kernel
#'
#' @seealso \code{\link{kernels_fire}}
#'
#' @examples
#' kernels <- list(cfbm, rbf)
#' params <- list(0.5, 1.0)  # Hurst = 0.5, sigma = 1.0
#' data <- list(matrix(rnorm(4), ncol=2), matrix(rnorm(6), ncol=3))
#' G_list <- fire:::gmat(kernels, params, data)
gmat <- function(kernels, kernels_params, dat, center = FALSE, std = FALSE){
  # list of kernel functions g1, g2, g3
  # list of parameters for each kernel (use NA for kernels without parameters)
  # list of data corresponding to each mode

  Nkernels = length(kernels)

  if(length(kernels_params) != Nkernels){
    stop("Mismatched lengths: kernels_params has ", length(kernels_params),
         " elements but ", Nkernels, " kernels were specified")
    }

  if (any(sapply(kernels_params, is.null))) {
    stop("All elements of kernels_params must be non-NULL.")
  }

  # compute each matrix G1, G2, G3, save as a list without scale parameter alpha
  G_matrices <- lapply(1:Nkernels, function(i) {
    # Extract the kernel function and its corresponding parameter
    kernel_func = kernels[[i]]
    param = kernels_params[[i]]

    # Call the kernel function with the data and parameter (if applicable)
    if (identical(kernel_func, cfbm)) {
      # For cfbm kernel, pass Hurst parameter
      kernel_func(dat[[i]], Hurst = param, std = std)
    } else if (identical(kernel_func, fbm)) {
      # For fbm kernel, pass Hurst parameter
      kernel_func(dat[[i]], Hurst = param[[1]], std = std)
    }else if (identical(kernel_func, rbf)) {
      # For rbf kernel, pass sigma parameter
      kernel_func(dat[[i]], lengthscale = param[[1]], center = center, std = std)
    } else if (identical(kernel_func, kronecker_delta)) {
      # For kronecker_delta kernel, no parameter is needed
      kernel_func(dat[[i]], center = center, std = std)
    }else if (identical(kernel_func, polynomial)) {
      # For polynomial kernel, pass d parameter, offset parameter
      kernel_func(X = dat[[i]], d = param[[1]], offset = param[[2]], center = center, std = std)
    }else if (identical(kernel_func, mercer)) {
      # For mercer-like kernel, pass delta and max_term parameters
      kernel_func(dat[[i]], delta = param[[1]], std = std)
    }else if (identical(kernel_func, matern)) {
      # For matern kernel, pass nu, lengthscale, sigma parameters
      kernel_func(dat[[i]], nu = param[[1]], lengthscale = param[[2]], sigma = param[[3]], std = std)
    } else {
      stop("Unsupported kernel function")
    }
  })

  return(G_list = G_matrices)
}
