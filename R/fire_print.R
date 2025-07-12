#' Print Method for FIRE Models
#'
#' @description
#' Compact display of FIRE models for either matrix or tensor input data.
#'
#' @param x A model object of class \code{fire_matrix} or \code{fire_tensor}.
#' @param ... Not used.
#'
#' @return The input object (invisibly) for piping. Prints to console:
#' \itemize{
#'   \item Model type (matrix or tensor input)
#'   \item Convergence status
#'   \item Data dimensions
#'   \item Marginal log-likelihood
#'   \item Estimated hyperparameters
#' }
#'
#' @seealso \code{\link{fire}}
#' @examples
#' # For matrix input
#' data(Manure)
#' mod <- fire(X = Manure$absorp[1:10,], Y = Manure$y$DM[1:10],
#'             dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4))
#' print(mod)
#'
#' # For tensor input
#' data(Housing)
#' dat_T <- list(T1 = 1:4, T2 = 1:9)
#' mod <- fire(X = Housing$X[1:5,,], Y = Housing$y[1:5,2],
#'             kernels = list(kronecker_delta, kronecker_delta),
#'             kernels_params = list(NA, NA),
#'             dat_T = dat_T, control = list(stop.eps = 2, maxiter = 4))
#' print(mod)
#' @name print.fire
NULL

#' @rdname print.fire
#' @method print fire_matrix
#' @export
print.fire_matrix <- function(x, ...) {
  cat(cli::rule(left = "FIRE Model (Matrix Input)", col = "blue"), "\n")

  conv_msg <- if (x$converged) {
    cli::col_green(sprintf("Converged in %d iterations", x$niter))
  } else {
    cli::col_red(sprintf("Stopped after %d iterations", x$niter))
  }
  cat(conv_msg, "\n")

  cat(cli::style_bold("Dimension:"),
      paste(attr(x, 'dimensions')[2]), "\n")

  cat(cli::style_bold("mLogLik:"),
      format(tail(x$mloglik, 1), digits = 4), "\n")

  cat("Estimated hyperparameters:\n")
  cat(sprintf("  Lambda: %0.4f\n", tail(x$lambda, 1)))
  cat(sprintf("  Noise: %0.4f\n", tail(x$noise, 1)))

  invisible(x)
}

#' @rdname print.fire
#' @method print fire_tensor
#' @export
print.fire_tensor <- function(x, ...) {
  cat(cli::rule(left = "FIRE Model (Tensor Input)", col = "blue"), "\n")

  conv_msg <- if (x$converged) {
    cli::col_green(sprintf("Converged in %d iterations", x$niter))
  } else {
    cli::col_red(sprintf("Stopped after %d iterations", x$niter))
  }
  cat(conv_msg, "\n")

  cat(cli::style_bold("Dimensions:"),
      paste(dim(attr(x, 'training_data')[[1]]), collapse = " x "), "\n")

  cat(cli::style_bold("mLogLik:"),
      format(tail(x$mloglik, 1), digits = 4), "\n")

  cat("Estimated hyperparameters:\n")
  cat(sprintf("  Alpha: %0.4f\n", tail(x$alpha, 1)))
  cat(sprintf("  Tau: %0.4f\n", tail(x$tau, 1)))
  cat(sprintf("  Noise: %0.4f\n", tail(x$noise, 1)))

  invisible(x)
}
