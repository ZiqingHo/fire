#' Prediction methods for FIRe models
#'
#' Obtain predicted values from \code{fire_matrix} or \code{fire_tensor} object, optionally with test set evaluation.
#'
#' @param object A \code{fire_matrix} or \code{fire_tensor} object
#' @param newdata New data for prediction
#' @param Ynew Optional numeric vector of true response values for the \code{newdata}. Must be the same length as the number of rows in \code{newdata}.
#' @param ... Not used
#' @return A list of class \code{fire_prediction} containing the predicted values,
#' and including test RMSE, residuals if \code{Ynew} is provided.
#'
#' The returned object has a \code{print} method that displays:
#' \itemize{
#'   \item Number of predictions
#'   \item Test RMSE (if available)
#'   \item First few predicted values
#' }
#'
#' @examples
#' \dontrun{
#' data(Manure)
#' idx <- 1:200
#' mod <- fire(X = Manure$absorp[idx,], Y = Manure$y$DM[idx], dat_T = list(1:700))
#' predict(mod, newdata = Manure$absorp[-idx,], Ynew = Manure$y$DM[-idx])
#' }
#'
#' @seealso \code{\link{fire}}

#' @rdname predict.fire
#' @export
predict.fire_matrix <- function(object, newdata, Ynew = NULL, ...) {
  # Extract components from model object
  X_train <- attr(object, "training_data")
  Y_train <- attr(object, "original_response")
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
  os_type <- attr(object, "os_type")
  intercept <- ifelse(is.null(attr(object, "intercept")), 0, attr(object, "intercept"))

  # Validate newdata
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  if (ncol(newdata) != ncol(X_train)) {
    stop("New data must have ", ncol(X_train), " columns matching training data")
  }
  if (!is.null(Ynew) && length(Ynew) != nrow(newdata)) {
    stop("Ynew length must match number of rows in newdata")
  }

  # Reconstruct dimensions
  d <- ncol(X_train)
  Index <- matrix(1:d, ncol = 1)

  # Construct kernel matrices
  G <- gmat(kernels = kernels,
            kernels_params = kernels_params,
            dat = dat_T,
            center = center)

  nmat <- Kronecker_norm_mat(X = X_train,
                             G = G,
                             alpha = c(1),
                             constant = constant,
                             Index = Index,
                             os_type = os_type,
                             sample_id = 1)

  nmat.cross <- Kronecker_norm_cross(Xtrain = X_train,
                                     Xnew = newdata,
                                     G = G,
                                     alpha = c(1),
                                     constant = constant,
                                     Index = Index,
                                     os_type = os_type,
                                     sample_id = 1)

  # Generate cross kernel matrix
  if (kernel_iprior == 'cfbm') {
    if (is.null(iprior_param)) iprior_param <- 0.5
    Hcross.tilde <- cfbm_rkhs_kron_cross(nmat = nmat,
                                         nmat_cross = nmat.cross,
                                         Hurst = iprior_param)
  } else if (kernel_iprior == 'rbf') {
    if (is.null(iprior_param)) iprior_param <- 1
    Hcross.tilde <- rbf_rkhs_kron_cross(nmat_cross = nmat.cross,
                                        lengthscale = iprior_param)
  } else if (kernel_iprior == 'linear') {
    Hcross.tilde <- nmat.cross + iprior_param
  }else if (kernel_iprior == 'poly'){
    Hcross.tilde <- (nmat.cross  + iprior_param[2])^iprior_param[1]
    }


  # Calculate predictions
  Ypred <- as.vector(intercept + lambda^2 * Hcross.tilde %*% w)

  # Create comprehensive output object
  result <- structure(
    list(
      yhat = Ypred,
      y.actual = Ynew,
      newdata = newdata,
      model = object,
      test_metrics = if (!is.null(Ynew)) {
        test_residuals <- Ynew - Ypred
        list(
          rmse = sqrt(mean(test_residuals^2)),
          residuals = test_residuals
        )
      }
    ),
    class = c("fire_prediction", "list")
  )

  # Print output automatically when not assigned
  if (!is.object(result) || sys.nframe() == 1) {
    print(result)
    return(invisible(result))
  }

  return(result)
}

#' @rdname predict.fire
#' @export
predict.fire_tensor <- function(object, newdata, Ynew = NULL, ...) {
  # Extract components from model object
  X_train <- attr(object, "training_data")
  Y_train <- attr(object, "original_response")
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

  # construct H matrix
  G = gmat(kernels = kernels, kernels_params = kernels_params, dat = dat_T, center = center)
  nmat <- Kronecker_norm_mat(X = X_train,
                             G = G,
                             alpha = alpha_params,
                             constant = constant,
                             Index = Index,
                             os_type = os_type,
                             sample_id = sample_id)

  nmat.cross <- Kronecker_norm_cross(Xtrain = X_train,
                                     Xnew = newdata,
                                     G = G,
                                     alpha = alpha_params,
                                     constant = constant,
                                     Index = Index,
                                     os_type = os_type,
                                     sample_id = sample_id)

  # Generate cross kernel matrix
  if (kernel_iprior == 'cfbm') {
    if (is.null(iprior_param)) iprior_param <- 0.5
    Hcross.tilde <- cfbm_rkhs_kron_cross(nmat = nmat,
                                         nmat_cross = nmat.cross,
                                         Hurst = iprior_param)
  } else if (kernel_iprior == 'rbf') {
    if (is.null(iprior_param)) iprior_param <- 1
    Hcross.tilde <- rbf_rkhs_kron_cross(nmat_cross = nmat.cross,
                                        lengthscale = iprior_param)
  } else if (kernel_iprior == 'linear') {
    Hcross.tilde <- nmat.cross + iprior_param
  }else if (kernel_iprior == 'poly'){
    Hcross.tilde <- (nmat.cross  + iprior_param[2])^iprior_param[1]
  }

  # Calculate predictions
  Ypred <- as.vector(intercept + tau^2 * Hcross.tilde %*% w)

  # Create comprehensive output object
  result <- structure(
    list(
      yhat = Ypred,
      y.actual = Ynew,
      newdata = newdata,
      model = object,
      test_metrics = if (!is.null(Ynew)) {
        test_residuals <- Ynew - Ypred
        list(
          rmse = sqrt(mean(test_residuals^2)),
          residuals = test_residuals
        )
      }
    ),
    class = c("fire_prediction", "list")
  )

  # Print output automatically when not assigned
  if (!is.object(result) || sys.nframe() == 1) {
    print(result)
    return(invisible(result))
  }

  return(result)
}


#' @rdname predict.fire
#' @param x A \code{fire_prediction} object to print
#' @export
print.fire_prediction <- function(x, ...) {
  cat("FIRe Model Predictions\n")
  cat("----------------------\n")
  cat(sprintf("Number of predictions: %d\n", length(x$yhat)))
  if (!is.null(x$test_metrics)) {
    cat(sprintf("Test RMSE: %.5f\n", x$test_metrics$rmse))
  }

  cat("\nFirst 6 predicted values:\n")
  print(head(x$yhat))
  if (length(x$yhat) > 6) {
    cat(sprintf("[... %d more values not shown]\n", length(x$yhat) - 6))
  }
  invisible(x)
}
