#' Classification Utility Functions
#'
#' Collection of helper functions for classification problems
#' @name utils_classification
#' @aliases to_dummy_vector onehot_argmax
NULL

#' @rdname utils_classification
#' @param x Categorical vector (character, numeric, or factor)
#' @param levels Vector of category levels.
#' @return \code{to_dummy_vector}: A numeric vector of length \eqn{n \times} \code{length(levels)}
#' @examples
#' to_dummy_vector(c("A","B","A"), levels = c("A","B","C"))
#' to_dummy_vector(c(0,1,0,9), levels = 0:9)
#' @export
to_dummy_vector <- function(x, levels = NULL) {
  # Build factor with specified levels (or infer if not provided)
  f <- if (!is.null(levels)) {
    factor(x, levels = levels)
  } else {
    if (is.factor(x)) x else factor(x)
  }

  # If levels were provided, ensure no unseen labels slipped in
  if (!is.null(levels) && any(is.na(f))) {
    bad <- unique(x[is.na(f)])
    stop("Found values not in `levels`: ", paste(bad, collapse = ", "))
  }

  dummy_mat <- model.matrix(~ f - 1)   # n x length(levels)
  as.vector(t(dummy_mat))              # row-major flatten
}


#' @rdname utils_classification
#' @param vec Numeric vector of length \eqn{n \times c}, where \eqn{n} is the number
#' of observations and \eqn{c} is the number of categories.
#' @param c Integer giving the number of categories.
#' @return \code{onehot_argmax}: A numeric vector of length \eqn{n \times c} containing a one-hot encoding
#' where, within each block of size \eqn{c}, the maximum element is set to 1 and all others to 0.
#'
#' @export
onehot_argmax <- function(vec, c) {
  n <- length(vec) / c
  stopifnot(n == floor(n))

  mat <- matrix(vec, nrow = n, ncol = c, byrow = TRUE)
  max_idx <- max.col(mat, ties.method = "first")

  res <- matrix(0, n, c)
  res[cbind(seq_len(n), max_idx)] <- 1

  as.vector(t(res))
}

#' @rdname utils_classification
#' @param onehot A numeric vector of length \eqn{n \times c}, the output from
#'   \code{onehot_argmax}.
#' @return \code{onehot_to_labels}: A vector of length \eqn{n} containing the
#'   predicted category labels corresponding to the one-hot encoding.
#'
#' @examples
#' set.seed(1)
#' v <- runif(12) # n=4, c=3
#' r <- onehot_argmax(v, c = 3)
#' r
#' onehot_to_labels(r, levels = c("A", "B", "C"))
#' @export
onehot_to_labels <- function(onehot, levels) {
  c <- length(levels)
  n <- length(onehot) / c
  stopifnot(n == floor(n))

  mat <- matrix(onehot, nrow = n, ncol = c, byrow = TRUE)
  idx <- max.col(mat, ties.method = "first")

  levels[idx]
}

