#' FIRE Model with k-Fold Cross-Validation
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
#' @param nfolds Number of folds for cross-validation (default: 10)
#' @param seed Random seed for reproducibility (default: 1)
#' @param ... Additional arguments passed to methods
#'
#' @return A list containing:
#'   - cv_results: Data frame with results for each fold
#'   - mean_results: Vector of mean performance metrics across folds
#'   - fold_indices: List of indices used for each fold
#'   - final_model: Model trained on full data (if return_full_model = TRUE)
#'
#' @examples
#' data(Manure)
#' cv_result <- fire_cv(X = Manure$absorp[1:10,], Y = Manure$y$DM[1:10],
#'  kernels = list(cfbm), kernels_params = list(0.5),
#'  dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4), nfold = 2)
#'
#' @seealso \code{\link{kernels_fire}}
#'
#' @export
fire_cv <- function(X, Y, dat_T, kernels, kernels_params,
                    kernel_iprior = "cfbm", iprior_param = NULL,
                    control = list(), nfolds = 10, seed = 1, ...) {

  # Set seed for reproducibility
  set.seed(seed)

  # Initialize results data frame
  cv_results <- data.frame(
    fold = 1:nfolds,
    train_rmse = numeric(nfolds),
    test_rmse = numeric(nfolds),
    stringsAsFactors = FALSE
  )

  # Create folds
  folds <- caret::createFolds(1:length(Y), k = nfolds)

  # Perform k-fold CV
  for (i in 1:nfolds) {
    cat("\nProcessing fold", i, "of", nfolds, "\n")

    # Split into training and testing sets
    test_indices <- folds[[i]]
    train_indices <- setdiff(1:length(Y), test_indices)

    # Handle different input types
    if (is.list(X) && !is.data.frame(X)) {
      X.train <- X[train_indices]
      X.test <- X[test_indices]
    } else {
      # Convert to matrix if data.frame
      if (is.data.frame(X)) {
        X <- as.matrix(X)
      }

      if (is.matrix(X)) {

          X.train <- X[train_indices, ]
          X.test <- X[test_indices, ]

      } else if (is.array(X)) {
        ndim <- length(dim(X))
        idx_train <- rep(list(quote(expr=)), ndim)
        idx_train[[control$sample_id]] <- train_indices
        idx_test <- rep(list(quote(expr=)), ndim)
        idx_test[[control$sample_id]] <- test_indices
        X.train <- do.call(`[`, c(list(X), idx_train))
        X.test <- do.call(`[`, c(list(X), idx_test))
      } else {
        stop("X must be a list, matrix, array, or data.frame")
      }
    }


    Y.train <- Y[train_indices]
    Y.test <- Y[test_indices]

    # Fit FIRE model with error handling
    fold_result <- tryCatch({
      mod <- fire(
        X = X.train,
        Y = Y.train,
        dat_T = dat_T,
        kernels = kernels,
        kernels_params = kernels_params,
        kernel_iprior = kernel_iprior,
        iprior_param = iprior_param,
        control = control,
        ...
      )

      if(inherits(mod, "fire_tensor")){
        # Training performance

        train_rmse <- fitted.fire_tensor(mod)$rmse

        # Test performance
        pred <- predict.fire_tensor(mod, newdata = X.test, Ynew = Y.test)
        test_rmse <- pred$test_metrics$rmse
      }else{
        # Training performance

        train_rmse <- fitted.fire_matrix(mod)$rmse

        # Test performance
        pred <- predict.fire_matrix(mod, newdata = X.test, Ynew = Y.test)
        test_rmse <- pred$test_metrics$rmse

      }

      list(
        train_rmse = train_rmse,
        test_rmse = test_rmse
      )
    }, error = function(e) {
      cat("Error in fold", i, ":", e$message, "\n")
      list(
        train_rmse = NA,
        test_rmse = NA
      )
    })

    # Store results
    cv_results$train_rmse[i] <- fold_result$train_rmse
    cv_results$test_rmse[i] <- fold_result$test_rmse

    cat("Fold", i, "completed - Train RMSE:", fold_result$train_rmse,
        "Test RMSE:", fold_result$test_rmse, "\n")
  }

  # Calculate mean performance (excluding NA results)
  mean_train_rmse <- mean(cv_results$train_rmse, na.rm = TRUE)
  mean_test_rmse <- mean(cv_results$test_rmse, na.rm = TRUE)

  # Count successful folds
  successful_folds <- sum(!is.na(cv_results$test_rmse))

  # Return results
  list(
    cv_results = cv_results,
    mean_results = c(
      mean_train_rmse = mean_train_rmse,
      mean_test_rmse = mean_test_rmse,
      success_rate = successful_folds/nfolds
    ),
    fold_indices = folds
  )
}
