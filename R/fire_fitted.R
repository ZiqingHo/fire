#' Fitted Values for FIRE Models
#'
#' @description
#' Computes fitted values for FIRE models with either matrix or tensor input data.
#'
#' @param object A model object of class \code{fire_matrix} or \code{fire_tensor}.
#' @param ... Not used.
#'
#' @return A list of class \code{fire_fitted} containing:
#' \itemize{
#'   \item \code{yhat}: Fitted values
#'   \item \code{rmse}: Training RMSE
#'   \item \code{residuals}: Model residuals
#'   \item \code{intercept}: Intercept value
#'   \item \code{model}: Reference to the original model object
#' }
#'
#' @seealso \code{\link{fire}}
#' @examples
#' # For matrix input
#' data(Manure)
#' mod <- fire(X = Manure$absorp[1:5,], Y = Manure$y$DM[1:5],
#'             dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4))
#' fitted(mod)
#'
#' # For tensor input
#' data(Housing)
#' dat_T <- list(T1 = 1:4, T2 = 1:9)
#' mod <- fire(X = Housing$X[1:5,,], Y = Housing$y[1:5,2],
#'             kernels = list(kronecker_delta, kronecker_delta),
#'             kernels_params = list(NA, NA),
#'             dat_T = dat_T, control = list(stop.eps = 2, maxiter = 4))
#' fitted(mod)
#' @name fitted.fire
NULL

#' @rdname fitted.fire
#' @method fitted fire_matrix
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
  constant_g <- attr(object, "constant_g")
  constant_h <- attr(object, "constant_h")
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
                             constant = constant_g,
                             Index = Index,
                             os_type = attr(object, "os_type"),
                             cores = attr(object, "cores"),
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

  if(constant_h){
    H.tilde <- 1 + H.tilde
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

  return(result)
}

#' @rdname fitted.fire
#' @method fitted fire_tensor
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
  constant_g <- attr(object, "constant_g")
  constant_h <- attr(object, "constant_h")
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
                             constant = constant_g,
                             Index = Index,
                             os_type = attr(object, "os_type"),
                             cores = attr(object, "cores"),
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

  if(constant_h){
    H.tilde <- 1 + H.tilde
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

  return(result)
}

#' Print FIRE Fitted Values
#'
#' @description
#' Prints summarized results for \code{fire_fitted} objects.
#'
#' @param x An object of class \code{fire_fitted}.
#' @param ... Not used.
#'
#' @return Invisibly returns the input object.
#' @method print fire_fitted
#' @seealso \code{\link{fitted.fire}}
#' @export
print.fire_fitted <- function(x, ...) {
  cat("FIRE Model Fitted Results\n")
  cat("-------------------------\n")
  cat(sprintf("Training RMSE: %.5f\n", x$rmse))
  cat(sprintf("Intercept: %.4f\n", x$intercept))

  if (length(x$yhat) > 6) {
    cat("\nFirst 6 fitted values:\n")
    print(head(x$yhat))
    cat(sprintf("[... %d more values not shown]\n", length(x$yhat) - 6))
  }else{
    cat("\nThe", length(x$yhat), "fitted values:\n")
    print(head(x$yhat))
  }
  invisible(x)
}

