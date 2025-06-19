#' Plot method for FIRe models
#'
#' The type of plot produced depends on
#' the class of the model object \code{fire_fitted} or \code{fire_prediction}
#'
#' @param x A \code{fire_fitted} or \code{fire_prediction} object
#' @param ... Not used
#' @return Returns a ggplot object or a gridExtra-arranged plot object.
#'   The exact return depends on the input:
#'   \itemize{
#'     \item For \code{fire_fitted} object: Returns a single ggplot object showing residuals vs fitted values
#'     \item For \code{fire_prediction} object:
#'       \itemize{
#'         \item When test data is available: Returns a gridExtra-arranged plot containing both
#'               residuals vs predicted and actual vs predicted plots
#'         \item When no test data is available: Returns a single ggplot object showing
#'               the distribution of predicted values with density overlay
#'       }
#'   }
#'
#' @examples
#' data(Manure)
#' idx <- 1:5
#' mod <- fire(X = Manure$absorp[idx,], Y = Manure$y$DM[idx],
#'  dat_T = list(1:700), stop.eps = 2, maxiter = 4)
#'
#' Yfitted = fitted(mod)
#' plot(Yfitted)
#'
#' Ypred = predict(mod, newdata = Manure$absorp[idx+10,],
#' Ynew = Manure$y$DM[idx+10])
#' plot(Ypred)

#' @rdname plot.fire
#' @export
plot.fire_fitted <- function(x, ...) {
  # Create residuals vs fitted plot
  p <- ggplot2::ggplot(data.frame(Fitted = x$yhat, Residuals = x$residuals),
                       ggplot2::aes(x = Fitted, y = Residuals)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Residuals vs Fitted",
                  subtitle = sprintf("Training RMSE: %.3f",
                                     x$rmse)) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}

#' @rdname plot.fire
#' @export
plot.fire_prediction <- function(x, ...) {
  if (!is.null(x$test_metrics)) {
    # If test values were provided, plot actual vs predicted with residuals
    df <- data.frame(
      Actual = x$y.actual,
      Predicted = x$yhat,
      Residuals = x$y.actual - x$yhat
    )

    # Common theme settings for both plots
    common_theme <- ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11),
        plot.subtitle = ggplot2::element_text(size = 9)
      )

    # Plot 1: Residuals vs Predictions
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = Predicted, y = Residuals)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "Residuals vs Predicted",
        subtitle = sprintf("Test RMSE: %.3f", x$test_metrics$rmse),
        x = "Predicted",
        y = "Residuals"
      ) +
      common_theme

    # Plot 2: Actual vs Predicted with blank subtitle
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = Predicted, y = Actual)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "Actual vs Predicted",
        subtitle = " ",  # Blank subtitle to maintain equal sizing
        x = "Predicted",
        y = "Actual"
      ) +
      common_theme

    # Combine plots with equal sizes
    p <- gridExtra::grid.arrange(
      p1, p2,
      ncol = 2,
      widths = c(1, 1)  # Equal widths for both plots
    )
    print(p)
    invisible(p)
  } else {
    # If no test values, show enhanced histogram
    warning("No test values provided - showing distribution of predictions")

    # Create data frame
    df_hist <- data.frame(Prediction = x$yhat)

    # Calculate automatic bins using Sturges' formula (default in ggplot2)
    n_bins <- nclass.Sturges(x$yhat)
    n_bins <- ifelse(n_bins > 10, n_bins, 10)

    # Calculate scaling factor for density curve
    bin_width <- diff(range(x$yhat))/n_bins
    scaling_factor <- length(x$yhat) * bin_width

    p <- ggplot2::ggplot(df_hist, ggplot2::aes(x = Prediction)) +
      ggplot2::geom_histogram(
        fill = "skyblue",
        color = "black",
        alpha = 0.8,
        bins = n_bins
      ) +
      ggplot2::stat_density(
        geom = "line",
        color = "darkblue",
        linewidth = 1,
        ggplot2::aes(y = ..density.. * scaling_factor)  # Scale density to counts
      ) +
      ggplot2::labs(
        title = "Distribution of Predicted Values",
        subtitle = sprintf("n = %d | bins = %d", length(x$yhat), n_bins),
        x = "Predicted",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9))

    print(p)
    invisible(p)
  }
}
