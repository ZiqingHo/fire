#' Fitted Values for FIRE Models
#'
#' @description
#' Computes fitted values for FIRE models with either matrix or tensor input data.
#'
#' @param object A model object of class \code{fire_matrix} or \code{fire_tensor}.
#' @param interval Logical indicating whether to compute credible band.
#' @param level Significance level for the intervals, between 0 and 1.
#' @param n.grid Number of grid points used to build credible band.
#' @param fixed A list of named vectors `c(index = i, value = v)` used to fix variables when building bands. Variables not specified are fixed at their column means.
#' @param ... Not used.
#' @return A list of class \code{fire_fitted} containing:
#' \itemize{
#'   \item \code{yhat}: Fitted values
#'   \item \code{rmse}: Training RMSE
#'   \item \code{residuals}: Model residuals
#'   \item \code{intercept}: Intercept value
#'   \item \code{model}: Reference to the original model object
#'   \item \code{CI}: A list with \code{upper}, \code{lower}, and \code{level} for training-point intervals.
#'   \item \code{Bands}: A list of length \eqn{d} (number of predictors). Each element \code{Bands[[j]]} is a list with:
#'     \describe{
#'       \item{\code{X.grid}}{Numeric vector of grid values for predictor \eqn{j}.}
#'       \item{\code{Ypred}}{Fitted mean on the grid (others fixed at column means).}
#'       \item{\code{upper}, \code{lower}}{Band limits on the grid.}
#'       \item{\code{level}}{Significance level.}
#'       \item{\code{fixed_at}}{Numeric vector of the column means used to fix the other predictors.}
#'       \item{\code{var_index}}{Integer index \eqn{j} of the varied predictor.}
#'     }
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
fitted.fire_matrix <- function(object, interval = FALSE, level = 0.05, n.grid = 100,
                               fixed = NULL, ...) {
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
  std <- attr(object, "std")
  intercept <- ifelse(is.null(attr(object, "intercept")), 0, attr(object, "intercept"))

  # Reconstruct dimensions and matrices
  d <- ncol(X)
  Index <- matrix(1:d, ncol = 1)

  G <- gmat(kernels = kernels,
            kernels_params = kernels_params,
            dat = dat_T,
            center = center, std = std)

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
  } else if (kernel_iprior == 'poly') {
    H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
  }

  if (constant_h) {
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
  Bands <- NULL

  # Construct confidence band
  if (interval) {

    noise <- tail(object$noise, 1)
    z <- qnorm(1 - level/2)
    sigma2 <- noise^2
    H <- lambda^2 * H.tilde
    n <- length(Y)
    Psi <- diag(n)/sigma2
    Vy <- H %*% H %*% Psi + Psi
    Vy.inv <- solve(Vy)
    Sigma.new <- (H) %*% (Psi - Psi %*% H %*% Vy.inv %*% t(H) %*% Psi) %*% t(H)
    se <- sqrt(diag(Sigma.new))
    upper <- Yfitted + z * se
    lower <- Yfitted - z * se

    # Baseline: fix others at column means, then override with user 'fixed'
    fix_at <- colMeans(as.matrix(X))

    if (!is.null(fixed)) {
      if (!is.list(fixed)) stop("'fixed' must be a list of named vectors like c(index=2, value=2).")

      for (el in fixed) {
        if (!is.atomic(el) || is.null(names(el))) {
          stop("Each element of 'fixed' must be a named vector like c(index=2, value=2).")
        }
        if (is.null(el[["index"]]) || is.null(el[["value"]])) {
          stop("Each element of 'fixed' must contain names 'index' and 'value'.")
        }

        idx <- as.integer(el[["index"]])
        val <- as.numeric(el[["value"]])

        if (is.na(idx) || idx < 1 || idx > d) {
          stop(sprintf("'fixed' index %s is out of bounds [1, %d].", as.character(el[["index"]]), d))
        }

        fix_at[idx] <- val
      }
    }

    # Per-variable credible bands (vary j; others fixed at fix_at)
    Bands <- vector("list", d)
    names(Bands) <- if (!is.null(colnames(X))) colnames(X) else paste0("X", 1:d)

    for (j in 1:d) {
      xj.min <- min(X[, j])
      xj.max <- max(X[, j])
      xj.grid <- seq(from = xj.min, to = xj.max, length.out = n.grid)

      X.grid <- matrix(rep(fix_at, each = n.grid), nrow = n.grid, byrow = FALSE)
      X.grid[, j] <- xj.grid

      nmat.cross <- Kronecker_norm_cross(Xtrain = X,
                                         Xnew = X.grid,
                                         G = G,
                                         alpha = c(1),
                                         constant = constant_g,
                                         Index = Index,
                                         os_type =  attr(object, "os_type"),
                                         cores = attr(object, "cores"),
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
      } else if (kernel_iprior == 'poly') {
        Hcross.tilde <- (nmat.cross + iprior_param[2])^iprior_param[1]
      }

      if (constant_h) {
        Hcross.tilde <- 1 + Hcross.tilde
      }

      # Calculate predictions
      Hcross <- lambda^2 * Hcross.tilde
      Ypred <- as.vector(intercept + Hcross %*% w)

      Sigma.new <- (Hcross) %*% (Psi - Psi %*% H %*% Vy.inv %*% t(H) %*% Psi) %*% t(Hcross) + sigma2
      se.new <- sqrt(diag(Sigma.new))

      Bands[[j]] <- list(
        X.grid   = xj.grid,
        Ypred    = Ypred,
        upper    = Ypred + z * se.new,
        lower    = Ypred - z * se.new,
        level    = level,
        fixed_at = fix_at,   # record actual fixed values used
        var_index = j
      )
    }
  }

  # Create comprehensive output object
  result <- structure(
    list(
      yhat = Yfitted,
      rmse = rmse,
      residuals = residuals,
      intercept = intercept,
      model = object,  # Include model reference
      CI = if (interval) {
        list(upper = upper, lower = lower, level = level)
      },
      Bands = Bands
    ),
    class = c("fire_fitted", "list")
  )

  return(result)
}


#' @rdname fitted.fire
#' @method fitted fire_tensor
#' @export
fitted.fire_tensor<- function(object,interval = FALSE, level = 0.05, ...) {
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

  # Calculate fitted values and metrics
  Yfitted <- as.vector(intercept + tau^2 * H.tilde %*% w)
  if (intercept != 0) {
    residuals <- Y + intercept - Yfitted  # Y stored is centered using mean(Y)
  } else {
    residuals <- Y - Yfitted
  }
  rmse <- sqrt(mean(residuals^2))

  if(interval){

    noise <- tail(object$noise, 1)
    z = qnorm(1 - level/2)
    sigma2 = noise^2
    H <- lambda^2 * H.tilde
    n <- length(Y)
    Psi <- diag(n)/sigma2
    Vy <- H %*% H %*% Psi + Psi
    Vy.inv <- solve(Vy)
    Sigma.new <- (H) %*% (Psi - Psi %*% H %*% Vy.inv %*% t(H) %*% Psi) %*% t(H)
    se <- sqrt(diag(Sigma.new))
    upper <- Yfitted + z * se
    lower <- Yfitted - z * se

  }
  # Create comprehensive output object
  result <- structure(
    list(
      yhat = Yfitted,
      rmse = rmse,
      residuals = residuals,
      intercept = intercept,
      model = object,  # Include model reference
      CI = if(interval){
        list(upper = upper, lower = lower, level = level)
      }
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

  if(!is.null(x$CI)){
    lower <- x$CI$lower
    upper <- x$CI$upper
    result <- data.frame(lower, x$yhat, upper)
    lower.name <- paste0(x$CI$level / 2 * 100, "%")
    upper.name <- paste0((1 - x$CI$level / 2) * 100, "%")
    names(result) <- c(lower.name, "Mean", upper.name)
    index <- min(length(x$yhat), 6)
    cat("\nFitted values:\n")
    print(result[seq_len(index),])
  }else{
    if (length(x$yhat) > 6) {
      cat("\nFirst 6 fitted values:\n")
      print(head(x$yhat))
      cat(sprintf("[... %d more values not shown]\n", length(x$yhat) - 6))
    }else{
      cat("\nThe", length(x$yhat), "fitted values:\n")
      print(head(x$yhat))
    }
  }

  invisible(x)
}

