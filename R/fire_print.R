#' Print method for FIRE models
#'
#' Compact display of \code{fire_matrix} or \code{fire_tensor} objects.
#'
#' @param x A \code{fire_matrix} or \code{fire_tensor} object
#' @param ... Not used
#'
#' @return The input object (invisibly) for piping. Prints to console:
#' \itemize{
#' \item{Basic model information}
#' \item{Convergence status}
#' \item{Dimensions}
#' \item{Estimated hyperparameters}
#' \item{Marginal log-likelihood}
#' }
#'
#' @examples
#' data(Manure)
#' mod <- fire(X = Manure$absorp[1:10,], Y = Manure$y$DM[1:10],
#'  dat_T = list(1:700), stop.eps = 2, maxiter = 4)
#' print(mod)
#'
#' @seealso \code{\link{fire}}
#' @name print.fire
NULL

#' @rdname print.fire
#' @export
print.fire_matrix <- function(x, ...) {
  cat(cli::rule(left = "FIRE Model (Matrix Input)", col = "blue"), "\n")

  conv_msg <- if (x$converged) {
    cli::col_green(sprintf("Converged in %d iterations", x$niter))
  } else {
    cli::col_red(sprintf("Stopped after %d iterations", x$niter))
  }
  cat(conv_msg, "\n")

  cat(cli::style_bold("Dimensions:"),
      paste(dim(x$training_data), collapse = " x "), "\n")

  cat(cli::style_bold("mLogLik:"),
      format(tail(x$mloglik, 1), digits = 4), "\n")

  cat("Estimated hyperparameters:\n")
  cat(sprintf("  Lambda: %0.4f\n", tail(x$lambda, 1)))
  cat(sprintf("  Noise: %0.4f\n", tail(x$noise, 1)))

  invisible(x)
}

#' @rdname print.fire
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
      paste(x$dimensions, collapse = " x "), "\n")

  cat(cli::style_bold("mLogLik:"),
      format(tail(x$mloglik, 1), digits = 4), "\n")

  cat("Estimated hyperparameters:\n")
  cat(sprintf("  Alpha: %0.4f\n", tail(x$alpha, 1)))
  cat(sprintf("  Tau: %0.4f\n", tail(x$tau, 1)))
  cat(sprintf("  Noise: %0.4f\n", tail(x$noise, 1)))

  invisible(x)
}


