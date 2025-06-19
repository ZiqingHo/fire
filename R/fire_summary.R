#' Summary method for FIRe models
#'
#' Provides a comprehensive summary of \code{fire_matrix} or \code{fire_tensor} object, including estimated parameters,
#' convergence information, model fit statistics, and kernel specifications.
#'
#' @param object A \code{fire_matrix} or \code{fire_tensor} object
#' @param ... Not used
#'
#' @return
#' For \code{summary.fire_matrix} and \code{summary.fire_tensor}:
#' Returns an object of class \code{summary.fire_matrix} or \code{summary.fire_tensor}
#' containing:
#' \itemize{
#'   \item Estimated parameters
#'   \item Convergence information
#'   \item Model fit statistics
#'   \item Kernel specifications
#'   \item Computation timing
#' }
#'
#' The \code{print} methods for these objects display nicely formatted output to the console.
#'
#' @examples
#' \dontrun{
#' data(Manure)
#' mod <- fire(X = Manure$absorp, Y = Manure$y$DM, dat_T = list(1:700))
#' summary(mod)
#' }

#' @rdname summary.fire
#' @export
summary.fire_matrix <- function(object, ...) {
  # Get kernel information from attributes
  kernels <- attr(object, "kernels")
  kernel_iprior <- attr(object, "kernel_iprior")
  kernel_params <- attr(object, "kernels_params")

  # Get the original call to extract kernel names as specified
  mc <- match.call()
  parent_call <- attr(object, "call")

  # Extract kernel names exactly as specified in the original call
  if (!is.null(parent_call$kernels)) {
    kernel_expr <- parent_call$kernels
    if (is.call(kernel_expr) && kernel_expr[[1]] == "list") {
      kernel_names <- sapply(as.list(kernel_expr)[-1], function(x) {
        if (is.call(x)) as.character(x[[1]]) else as.character(x)
      })
    } else {
      kernel_names <- as.character(kernel_expr)
    }
  } else {
    # Fallback (shouldn't happen since fire.matrix requires kernels)
    kernel_names <- sapply(kernels, function(k) as.character(substitute(k)))
  }

  structure(
    list(
      coefficients = c(
        lambda = tail(object$lambda, 1),
        noise = tail(object$noise, 1),
        intercept = ifelse(is.null(attr(object, "intercept")), 0,
                           attr(object, "intercept"))
      ),
      convergence = list(
        converged = object$converged,
        iterations = object$niter,
        maxiter = attr(object, "convergence")$maxiter,
        tolerance = attr(object, "convergence")$tolerance,
        final_change = attr(object, "convergence")$final_change
      ),
      model_fit = list(
        mloglik = tail(object$mloglik, 1),
        training_rmse = if(!is.null(object$training_rmse)) object$training_rmse else NA,
        n = nrow(attr(object, "training_data"))
      ),
      timing = object$duration,
      dimensions = attr(object, "dimensions"),
      kernel_info = list(
        kernels = kernel_names,
        kernel_params = kernel_params,
        kernel_iprior = kernel_iprior
      )
    ),
    class = "summary.fire_matrix"
  )
}

#' @rdname summary.fire
#' @param x A \code{summary.fire_matrix} or \code{summary.fire_tensor} object
#' @export
print.summary.fire_matrix <- function(x, ...) {
  # Header with colored title
  cat(cli::rule(left = crayon::blue("FIRe Model Summary"),
                width = getOption("width")), "\n")

  # Model information
  cat(crayon::bold("Data Dimensions:"),
      paste0(x$dimensions[1], " samples x ", x$dimensions[2], " features\n"))

  # Display kernel information with parameters
  kernel_text <- paste0("RKHS norm: ", x$kernel_info$kernels[1])

  # Add parameters if they exist
  if (!is.null(x$kernel_info$kernel_params[[1]])) {
    kernel_text <- paste0(kernel_text, " (",
                          paste(x$kernel_info$kernel_params[[1]]),
                          ")")
  }

  cat(crayon::bold("Kernels:\n"),
      paste0(kernel_text,
             " | I-prior: ", x$kernel_info$kernel_iprior, "\n\n"))

  # Coefficients section
  cat(crayon::bold("Estimated Parameters:"))
  coef_table <- data.frame(
    Parameter = c("Lambda", "Noise", "Intercept"),
    Value = format(unlist(x$coefficients), digits = 4, nsmall = 4),
    stringsAsFactors = FALSE
  )
  print(knitr::kable(coef_table, format = "simple", align = c('l', 'r')))
  cat("\n")

  # Convergence information
  cat(crayon::bold("Convergence:\n"))
  conv_details <- if(x$convergence$converged) {
    crayon::green(paste0("Converged in ", x$convergence$iterations, " iterations"))
  } else {
    crayon::red(paste0("Failed to converge after ", x$convergence$iterations,
                       "/", x$convergence$maxiter, " iterations"))
  }
  cat(paste(" Status:", conv_details), "\n")
  cat(paste(" Final change:", format(x$convergence$final_change, digits = 4)),
      paste0("(tolerance = ", format(x$convergence$tolerance, digits = 2), ")\n"))

  # Model fit statistics
  cat(crayon::bold("\nModel Fit:"))
  fit_table <- data.frame(
    Statistic = c("Marginal log-likelihood", "Sample size"),
    Value = c(format(x$model_fit$mloglik, digits = 4),
              x$model_fit$n),
    stringsAsFactors = FALSE
  )
  print(knitr::kable(fit_table, format = "simple", align = c('l', 'r')))

  # Timing information
  cat("\n")
  cat(crayon::bold("Computation time:"),
      format(x$timing, digits = 3), "\n")

  # Footer
  cat(cli::rule(width = getOption("width")))
}

#' @rdname summary.fire
#' @export
summary.fire_tensor <- function(object, ...) {
  # Get kernel information from attributes
  kernels <- attr(object, "kernels")
  kernel_iprior <- attr(object, "kernel_iprior")
  kernel_params <- attr(object, "kernels_params")
  dimensions <- attr(object, "dimensions")
  input_type <- attr(object, "input_type")
  sample_id <- attr(object, "sample_id")
  sample_size <- attr(object, "sample_size")
  # Get the original call to extract kernel names as specified
  parent_call <- attr(object, "call")

  # Extract kernel names exactly as specified in the original call
  if (!is.null(parent_call$kernels)) {
    kernel_expr <- parent_call$kernels
    if (is.call(kernel_expr) && kernel_expr[[1]] == "list") {
      kernel_names <- sapply(as.list(kernel_expr)[-1], function(x) {
        if (is.call(x)) as.character(x[[1]]) else as.character(x)
      })
    } else {
      kernel_names <- as.character(kernel_expr)
    }
  } else {
    # Fallback (shouldn't happen since fire_tensor requires kernels)
    kernel_names <- sapply(kernels, function(k) as.character(substitute(k)))
  }

  # Prepare dimension information
  dim_info <- if (input_type == "array") {
    paste0(dimensions, collapse = " x ")
  } else {
    paste0(length(dimensions[[1]]), " samples x ",
           length(dimensions) - 1, " modes")
  }

  structure(
    list(
      coefficients = c(
        alpha = tail(object$alpha, 1),
        tau = tail(object$tau, 1),
        noise = tail(object$noise, 1),
        intercept = ifelse(is.null(attr(object, "intercept")), 0,
                           attr(object, "intercept"))
      ),
      convergence = list(
        converged = object$converged,
        iterations = object$niter,
        maxiter = attr(object, "convergence")$maxiter,
        tolerance = attr(object, "convergence")$tolerance,
        final_change = attr(object, "convergence")$final_change
      ),
      model_fit = list(
        mloglik = tail(object$mloglik, 1),
        n = length(attr(object, "original_response"))
      ),
      timing = object$duration,
      dimensions = dim_info,
      input_type = input_type,
      sample_id = sample_id,
      sample_size = sample_size,
      kernel_info = list(
        kernels = kernel_names,
        kernel_params = kernel_params,
        kernel_iprior = kernel_iprior,
        iprior_param = attr(object, "iprior_param")
      )
    ),
    class = "summary.fire_tensor"
  )
}

#' @rdname summary.fire
#' @export
print.summary.fire_tensor <- function(x, ...) {
  # Header with colored title
  cat(cli::rule(left = crayon::blue("FIRe Tensor Model Summary"),
                width = getOption("width")), "\n")

  # Model information
  cat(crayon::bold("Data Dimensions:"), x$dimensions, "\n")
  cat(crayon::bold("Input Type:"), x$input_type, "\n")
  if (!is.null(x$sample_size)) {
    cat(crayon::bold("Sample size:"), x$sample_size, "\n")
  }

  # Display kernel information with parameters
  cat(crayon::bold("RKHS norm kernels:\n"))
  kernel_text <- mapply(function(kernel, params) {
    param_text <- if (!is.null(params)) {
      paste0(" (", paste(params, collapse = ", "), ")")
    } else ""
    paste0("- ", kernel, param_text)
  }, x$kernel_info$kernels, x$kernel_info$kernel_params)

  cat(paste(kernel_text, collapse = "\n"), "\n")
  cat("I-prior kernel:", x$kernel_info$kernel_iprior)
  if (!is.null(x$kernel_info$iprior_param)) {
    cat(" (parameter =", x$kernel_info$iprior_param, ")")
  }
  cat("\n")

  # Coefficients section
  cat(crayon::bold("Estimated Parameters:"), "\n")
  coef_table <- data.frame(
    Value = format(unlist(x$coefficients), digits = 4, nsmall = 4),
    stringsAsFactors = FALSE
  )
  print(knitr::kable(coef_table, format = "simple", align = c('l', 'r')))
  cat("\n")

  # Convergence information
  cat(crayon::bold("Convergence:\n"))
  conv_details <- if(x$convergence$converged) {
    crayon::green(paste0("Converged in ", x$convergence$iterations, " iterations"))
  } else {
    crayon::red(paste0("Failed to converge after ", x$convergence$iterations,
                       "/", x$convergence$maxiter, " iterations"))
  }
  cat(paste(" Status:", conv_details), "\n")
  cat(paste(" Final change:", format(x$convergence$final_change, digits = 4)),
      paste0("(tolerance = ", format(x$convergence$tolerance, digits = 2), ")\n"))

  # Model fit statistics
  cat(crayon::bold("\nModel Fit:"))
  fit_table <- data.frame(
    Statistic = c("Marginal log-likelihood", "Sample size"),
    Value = c(format(x$model_fit$mloglik, digits = 4),
              x$model_fit$n),
    stringsAsFactors = FALSE
  )
  print(knitr::kable(fit_table, format = "simple", align = c('l', 'r')))

  # Timing information
  cat("\n")
  cat(crayon::bold("Computation time:"),
      format(x$timing, digits = 3), "\n")

  # Footer
  cat(cli::rule(width = getOption("width")))
}
