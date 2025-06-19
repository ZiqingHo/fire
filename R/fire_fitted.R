#' Obtain fitted values of FIRe model
#'
#' @param object A \code{fire_matrix} or \code{fire_tensor} object
#' @param ... Not used
#'
#' @return A list of class \code{fire_fitted} containing the fitted values,
#' training RMSE, residuals and intercept.
#'
#' The returned object has a \code{print} method that displays:
#' \itemize{
#'   \item Training RMSE
#'   \item Intercept value
#'   \item First few fitted values
#' }
#'
#' @examples
#' \dontrun{
#' data(Manure)
#' mod <- fire(X = Manure$absorp, Y = Manure$y$DM, dat_T = list(1:700))
#' fitted(mod)
#' }
#'
#' @seealso \code{\link{fire}}

#' @rdname fitted.fire
#' @export
fitted.fire_matrix <- function(object, ...) {
  # Extract components from model object
  X <- attr(object, "training_data")
  Y <- attr(object, "original_response")
  w <- object$w[[length(object$w)]]
  lambda <- tail(object$lambda, 1)

  # Get parameters from object attributes
  kernels <- attr(object, "kernels")
  kernels_params <- attr(object, "kernels_params")
  kernel_iprior <- attr(object, "kernel_iprior")
  iprior_param <- attr(object, "iprior_param")
  constant <- attr(object, "constant")
  dat_T <- attr(object, "dat_T")
  center <- attr(object, "center")
  intercept <- ifelse(is.null(attr(object, "intercept")), 0, attr(object, "intercept"))

  # Reconstruct dimensions and matrices
  d <- ncol(X)
  Index <- matrix(1:d, ncol = 1)

  G <- gmat(kernels = kernels,
            kernels_params = kernels_params,
            dat = dat_T,
            center = center)

  nmat <- Kronecker_norm_mat(X = X,
                             G = G,
                             alpha = c(1),
                             constant = constant,
                             Index = Index,
                             os_type = attr(object, "os_type"),
                             sample_id = 1)

  # Generate Gram matrix
  if (kernel_iprior == 'cfbm') {
    if (is.null(iprior_param)) iprior_param <- 0.5
    H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
  } else if (kernel_iprior == 'rbf') {
    if (is.null(iprior_param)) iprior_param <- 1
    H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
  } else if (kernel_iprior == 'linear') {
    H.tilde <- nmat + iprior_param
  } else if (kernel_iprior == 'poly'){
    H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
  }

  # Calculate fitted values and metrics
  Yfitted <- as.vector(intercept + lambda^2 * H.tilde %*% w)
  if (intercept != 0) {
    residuals <- Y - Yfitted + intercept  # Y stored is centered using mean(Y)
  } else {
    residuals <- Y - Yfitted
  }
  rmse <- sqrt(mean(residuals^2))

  # Create comprehensive output object
  result <- structure(
    list(
      yhat = Yfitted,
      rmse = rmse,
      residuals = residuals,
      intercept = intercept,
      model = object  # Include model reference
    ),
    class = c("fire_fitted", "list")
  )

  # Always return visibly for direct calls
  if (!is.object(result) || sys.nframe() == 1) {
    print(result)
    return(invisible(result))
  }

  return(result)
}

#' @rdname fitted.fire
#' @export
fitted.fire_tensor<- function(object, ...) {
  # Extract components from model object
  X <- attr(object, "training_data")
  Y <- attr(object, "original_response")
  w <- object$w[[length(object$w)]]
  alpha <- tail(object$alpha, 1)
  tau <- tail(object$tau, 1)

  # Get parameters from object attributes
  kernels <- attr(object, "kernels")
  kernels_params <- attr(object, "kernels_params")
  kernel_iprior <- attr(object, "kernel_iprior")
  iprior_param <- attr(object, "iprior_param")
  constant <- attr(object, "constant")
  dat_T <- attr(object, "dat_T")
  center <- attr(object, "center")
  os_type <- attr(object, "os_type")
  sample_id <- attr(object, "sample_id")
  intercept <- ifelse(is.null(attr(object, "intercept")), 0, attr(object, "intercept"))

  alpha_params = rep(alpha, length(kernels))

  Index <- create_index_matrix(dat_T)

  G <- gmat(kernels = kernels,
            kernels_params = kernels_params,
            dat = dat_T,
            center = center)

  nmat <- Kronecker_norm_mat(X = X,
                             G = G,
                             alpha = alpha_params,
                             constant = constant,
                             Index = Index,
                             os_type = attr(object, "os_type"),
                             sample_id = sample_id)

  # Generate Gram matrix
  if (kernel_iprior == 'cfbm') {
    if (is.null(iprior_param)) iprior_param <- 0.5
    H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
  } else if (kernel_iprior == 'rbf') {
    if (is.null(iprior_param)) iprior_param <- 1
    H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
  } else if (kernel_iprior == 'linear') {
    H.tilde <- nmat + iprior_param
  } else if (kernel_iprior == 'poly'){
    H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
  }

  # Calculate fitted values and metrics
  Yfitted <- as.vector(intercept + tau^2 * H.tilde %*% w)
  if (intercept != 0) {
    residuals <- Y - Yfitted + intercept  # Y stored is centered using mean(Y)
  } else {
    residuals <- Y - Yfitted
  }
  rmse <- sqrt(mean(residuals^2))

  # Create comprehensive output object
  result <- structure(
    list(
      yhat = Yfitted,
      rmse = rmse,
      residuals = residuals,
      intercept = intercept,
      model = object  # Include model reference
    ),
    class = c("fire_fitted", "list")
  )

  # Always return visibly for direct calls
  if (!is.object(result) || sys.nframe() == 1) {
    print(result)
    return(invisible(result))
  }

  return(result)
}

#' @rdname fitted.fire
#' @param x A \code{fire_fitted} object to print
#' @export
print.fire_fitted <- function(x, ...) {
  cat("FIRe Model Fitted Results\n")
  cat("-------------------------\n")
  cat(sprintf("Training RMSE: %.5f\n", x$rmse))
  cat(sprintf("Intercept: %.4f\n", x$intercept))
  cat("\nFirst 6 fitted values:\n")
  print(head(x$yhat))
  if (length(x$yhat) > 6) {
    cat(sprintf("[... %d more values not shown]\n", length(x$yhat) - 6))
  }
  invisible(x)
}

