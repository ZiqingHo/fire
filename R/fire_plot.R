#' Plot for Fitted FIRE Models
#'
#' @description
#' Generates residuals vs fitted values plot for training data.
#'
#' @param x A \code{fire_fitted} object
#' @param ... Not used
#'
#' @return Returns a single \code{ggplot} object showing residuals vs fitted values
#'
#' @examples
#' data(Manure)
#' mod <- fire(X = Manure$absorp[1:5,], Y = Manure$y$DM[1:5],
#'             dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4))
#' fit <- fitted(mod)
#' plot(fit)
#' @method plot fire_fitted
#'
#' @seealso \code{\link{fire}}, \code{\link{fitted.fire}}
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

#' Plot for FIRE Predictions
#'
#' @description
#' Generates residuals vs predicted and actual vs predicted plots when test response data is provided.
#' Otherwise, generates the distribution of predicted values with density overlay.
#'
#' @param x A \code{fire_prediction} object
#' @param ... Not used
#'
#' @return
#'  \itemize{
#'    \item When test response data is available: Returns both
#'               residuals vs predicted and actual vs predicted plots
#'    \item When no test response data is available: Returns a single \code{ggplot} object showing
#'               the distribution of predicted values with density overlay
#'  }
#' @examples
#' data(Manure)
#' mod <- fire(X = Manure$absorp[1:5,], Y = Manure$y$DM[1:5],
#'             dat_T = list(1:700), control = list(stop.eps = 2, maxiter = 4))
#' pred <- predict(mod, newdata = Manure$absorp[6:10,],
#'                 Ynew = Manure$y$DM[6:10])
#' plot(pred)
#' @method plot fire_prediction
#'
#' @seealso \code{\link{fire}}, \code{\link{predict.fire}}
#' @export
plot.fire_prediction <- function(x, interval = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
  }


  if (!is.null(x$test_metrics)) {

    if(interval){ # plot prediction interval plot
      df <- data.frame(
        Actual = x$y.actual,
        Predicted = x$yhat,
        upper = x$CI$upper,
        lower = x$CI$lower
      )
      level <- x$CI$level
      df <- df[order(df$Predicted), ]
      df$Sample <- seq_len(nrow(df))

      p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample)) +
        # Interval with caps
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper, colour = "Interval"),
                      width = 0.4, size = 0.5) +
        # Actual value (solid circle)
        ggplot2::geom_point(ggplot2::aes(y = Actual, colour = "Actual"), shape = 16, size = 2) +
        ggplot2::scale_colour_manual(values = c("Interval" = "steelblue", "Actual" = "firebrick")) +
        ggplot2::scale_x_continuous(
          limits = c(0.5, max(df$Sample) + 0.5),          # start at 1, end at m
          breaks = function(lims) {               # let R pick, but ensure 1 and integers
            b <- scales::breaks_pretty()(lims)    # automatic nice spacing
            b <- unique(c(1, b))                  # force a tick at 1
            b <- as.integer(round(b))             # integer ticks only
            b[b >= 1 & b <= max(df$Sample)]
          },
          labels = scales::label_number(accuracy = 1)
        ) +
        ggplot2::labs(
          title = "Prediction Interval Plot",
          x = "Sample",
          y = "Value",
          colour = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
              legend.position = "bottom")
      print(p)
      invisible(p)
    }else{


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
    p <- p1 + p2
    print(p)
    invisible(p)
    }
  } else {
    # If no test values, show enhanced histogram
    warning("No test values provided - showing distribution of predictions")

    # Create data frame
    df_hist <- data.frame(Prediction = x$yhat)

    # Calculate automatic bins using Sturges' formula (default in ggplot2)
    n_bins <- max(nclass.Sturges(x$yhat), 10)

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
        ggplot2::aes(y = ggplot2::after_stat(density) * scaling_factor)  # Updated to use after_stat()
      ) +
      ggplot2::labs(
        title = "Distribution of Predicted Values",
        subtitle = sprintf("n = %d | %d bins", length(x$yhat), n_bins),
        x = "Predicted",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9))
    print(p)
    invisible(p)
    }

}
