#' Fitting and prediction methods for FIRE classification models
#'
#' @description
#' These functions provide fitted values, predictions, and print summaries
#' for \code{fire_class} objects.
#'
#' @param object A fitted \code{fire_class} model object.
#' @param newdata New input data (same structure as used in training).
#' @param Ynew Optional true class labels for \code{newdata}.
#' @param x An object returned by \code{fitted.fire_class} or
#'   \code{predict.fire_class}.
#' @param ... Not used.
#'
#' @return
#' \itemize{
#'   \item For \code{fitted.fire_class}: an object of class
#'   \code{"fire_class_fitted"} with elements \code{yhat}, \code{accuracy}, and
#'   \code{model}.
#'   \item For \code{predict.fire_class}: an object of class
#'   \code{"fire_class_prediction"} with elements \code{yhat}, \code{y.actual},
#'   \code{model}, and optional \code{accuracy}.
#' }
#'
#' @seealso
#' \code{\link{fire_class}}
#'
#' @examples
#' set.seed(42)
#' n_train <- 5; n_test  <- 2; n <- n_train + n_test
#' MatA <- matrix(c(1,1,1,0,0,0,0,0,0), nrow = 3, byrow = TRUE)
#' MatB <- matrix(c(0,0,0,0,0,0,1,1,1), nrow = 3, byrow = TRUE)
#' Y <- factor(c("A","A","B","B","A","A","B"), levels = c("A","B"))
#' # Generate matrices with small Gaussian noise
#' X <- lapply(seq_len(n), function(i) {base <- if (Y[i] == "A") MatA else MatB
#' noise <- matrix(rnorm(9, sd = 0.15), nrow = 3)
#' base + noise })
#' # Train/test split: first 5 train, last 2 test
#' X_train <- X[1:n_train]; X_test  <- X[(n_train + 1):n]
#' Y_train <- Y[1:n_train]; Y_test  <- Y[(n_train + 1):n]
#' dat_T <- list(1:3, 1:3)
#' mod <- fire_class(X = X_train, Y = Y_train, dat_T = dat_T,
#'  kernels = list(cfbm, cfbm), kernels_params = list(0.5, 0.5),
#'  class.labels = levels(Y_train), control = list(maxiter = 20, stop.eps = 1e-3))
#'
#' fit <- fitted(mod)
#' print(fit)
#'
#' pred <- predict(mod, newdata = X_test, Ynew = Y_test)
#' print(pred)
#'
#' @name fire_class-methods
#' @rdname fire_class-methods
#' @export
fitted.fire_class<- function(object, ...) {
  # Extract components from model object
  X <- attr(object, "training_data")
  Y <- attr(object, "original_response")
  true.class <- attr(object, "labels")
  w <- object$w[[length(object$w)]]
  alpha <- tail(object$alpha, 1)
  tau <- tail(object$tau, 1)

  # Get parameters from object attributes
  kernels <- attr(object, "kernels")
  kernels_params <- attr(object, "kernels_params")
  kernel_iprior <- attr(object, "kernel_iprior")
  iprior_param <- attr(object, "iprior_param")
  kernel_class <- attr(object, "kernel_class")
  num_class <- attr(object, "num_class")
  class.labels <- attr(object, "class.labels")
  constant_g <- attr(object, "constant_g")
  constant_h <- attr(object, "constant_h")
  dat_T <- attr(object, "dat_T")
  center <- attr(object, "center")
  std <- attr(object, "std")
  os_type <- attr(object, "os_type")
  sample_id <- attr(object, "sample_id")
  intercept <- ifelse(is.null(attr(object, "intercept")), 0, attr(object, "intercept"))

  alpha_params = rep(alpha, length(kernels))

  Index <- create_index_matrix(dat_T)

  G <- gmat(kernels = kernels,
            kernels_params = kernels_params,
            dat = dat_T,
            center = center, std = std)

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

  if(kernel_class == 'identity'){
    H.group <- diag(num_class)
  }else if(kernel_class == 'centred identity'){
    H.group <- diag(num_class) - matrix(1/num_class, num_class, num_class)
  }else if(is.matrix(kernel_class)){
    H.group <- kernel_class
  }

  H.tilde <- kronecker(H.tilde, H.group)
  # Calculate fitted values and pick the largest value
  Yfitted <- as.vector(intercept + tau^2 * H.tilde %*% w)
  fitted.class <- onehot_argmax(Yfitted, num_class) # dummy encoding
  fitted.class <- onehot_to_labels(fitted.class, class.labels)
  accuracy <- mean(true.class == fitted.class)
  # Create comprehensive output object
  result <- structure(
    list(
      yhat = fitted.class,
      accuracy = accuracy,
      model = object  # Include model reference
    ),
    class = c("fire_class_fitted", "list")
  )

  return(result)
}

#' @rdname fire_class-methods
#' @export
print.fire_class_fitted <- function(x, ...) {
  cat("FIRE Model Fitted Results (Classification)\n")
  cat("-------------------------\n")
  cat(sprintf("Accuracy: %.5f\n", x$accuracy))

  if (length(x$yhat) > 6) {
      cat("\nFirst 6 fitted category labels:\n")
      print(head(x$yhat))
      cat(sprintf("[... %d more not shown]\n", length(x$yhat) - 6))
    }else{
      cat("\nThe", length(x$yhat), "fitted category labels:\n")
      print(head(x$yhat))
    }

  invisible(x)
}

#' @rdname fire_class-methods
#' @export
predict.fire_class <- function(object, newdata, Ynew = NULL,
                                ...) {
  # Extract components from model object
  X_train <- attr(object, "training_data")
  Y_train <- attr(object, "original_response")
  labels <- levels(attr(object, "labels"))
  w <- object$w[[length(object$w)]]
  alpha <- tail(object$alpha, 1)
  tau <- tail(object$tau, 1)

  # Get parameters from object attributes
  kernels <- attr(object, "kernels")
  kernels_params <- attr(object, "kernels_params")
  kernel_iprior <- attr(object, "kernel_iprior")
  iprior_param <- attr(object, "iprior_param")
  kernel_class <- attr(object, "kernel_class")
  num_class <- attr(object, "num_class")
  class.labels <- attr(object, "class.labels")
  constant_g <- attr(object, "constant_g")
  constant_h <- attr(object, "constant_h")
  dat_T <- attr(object, "dat_T")
  center <- attr(object, "center")
  std <- attr(object, "std")
  os_type <- attr(object, "os_type")
  cores <- attr(object, "cores")
  sample_id <- attr(object, "sample_id")
  intercept <- ifelse(is.null(attr(object, "intercept")), 0, attr(object, "intercept"))

  alpha_params = rep(alpha, length(kernels))

  Index <- create_index_matrix(dat_T)

  # construct H matrix
  G = gmat(kernels = kernels, kernels_params = kernels_params, dat = dat_T,
           center = center, std = std)
  nmat <- Kronecker_norm_mat(X = X_train,
                             G = G,
                             alpha = alpha_params,
                             constant = constant_g,
                             Index = Index,
                             os_type = os_type, cores = cores,
                             sample_id = sample_id)

  nmat.cross <- Kronecker_norm_cross(Xtrain = X_train,
                                     Xnew = newdata,
                                     G = G,
                                     alpha = alpha_params,
                                     constant = constant_g,
                                     Index = Index,
                                     os_type = os_type, cores = cores,
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

  if(constant_h){
    Hcross.tilde <- 1 + Hcross.tilde
  }

  if(kernel_class == 'identity'){
    H.group <- diag(num_class)
  }else if(kernel_class == 'centred identity'){
    H.group <- diag(num_class) - matrix(1/num_class, num_class, num_class)
  }else if(is.matrix(kernel_class)){
    H.group <- kernel_class
  }
  Hcross.tilde <- kronecker(Hcross.tilde, H.group)

  # Calculate predictions and pick the largest
  Ypred <- as.vector(intercept + tau^2 * Hcross.tilde %*% w)
  pred.class <- onehot_argmax(Ypred, num_class) # dummy encoding, n x c
  pred.class <- onehot_to_labels(pred.class, class.labels) # length n

  # Create comprehensive output object
  result <- structure(
    list(
      yhat = pred.class,
      y.actual = Ynew,
      model = object,
      accuracy = if (!is.null(Ynew)) {
        mean(Ynew == pred.class)
      }
    ),
    class = c("fire_class_prediction", "list")
  )

  # Print output automatically when not assigned
  if (sys.nframe() == 1) {
    print(result)
    return(invisible(result))
  }

  return(result)
}

#' @rdname fire_class-methods
#' @export
print.fire_class_prediction <- function(x, ...) {
  cat("FIRE Model Predictions (Classification)\n")
  cat("----------------------\n")
  cat(sprintf("Number of predictions: %d\n", length(x$yhat)))
  if (!is.null(x$accuracy)) {
    cat(sprintf("Accuracy: %.5f\n", x$accuracy))
  }

  if (length(x$yhat) > 6) {
      cat("\nFirst 6 predicted labels:\n")
      print(head(x$yhat))
      cat(sprintf("[... %d more not shown]\n", length(x$yhat) - 6))
    }else{
      cat("\nThe",length(x$yhat) , "predicted labels:\n")
      print(head(x$yhat))
    }

  invisible(x)
}
