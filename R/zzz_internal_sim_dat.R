#' @title Internal RKHS Norm Utilities
#' @name internal_sim_dat
#' @description Internal functions for computing RKHS norms and related operations.
#'   Used internally by \code{\link{sim_dat}}. Not exported or intended for direct use.
#' @keywords internal
NULL

#' @rdname internal_sim_dat
#'
#' @param X1 Sample matrix or vector
#' @param X2 Sample matrix or vector
#' @param Ginv (Generalized) Inverse of the Gram matrix
#' @param Index Matrix specifying tensor vectorization indices
#'
#' @return \code{rkhs_norm}: Numeric value representing the RKHS norm between X and Y
#'
#' @examples
#' G <- diag(1,4)
#' fire:::rkhs_norm(X1 = matrix(1:4, nrow=2),
#'           Ginv = solve(G),
#'           Index = expand.grid(1:2, 1:2))
rkhs_norm <- function(X1, X2 = NULL, Ginv, Index) {
  # Input validation
  if (!is.matrix(Ginv)) stop("Ginv must be a matrix")
  if (!is.matrix(Index)){
    Index <- as.matrix(Index)
  }

  if (is.matrix(X1)) {
    dim_list <- dim(X1)
    # vectorize tensor X1 and X2
    X1_vec <- vectorize_tensor(X1, Index)
    if (is.null(X2)) {
      X2_vec <- rep(0, prod(dim_list))
    } else {
      X2_vec <- vectorize_tensor(X2, Index)
    }
  } else {
    X1_vec <- X1
    if (is.null(X2)) {
      X2_vec <- rep(0, length(X1))
    } else {
      X2_vec <- X2
    }
  }

  # difference between X1 and X2
  tnsr <- X1_vec - X2_vec

  # compute norm
  r <- t(tnsr) %*% Ginv %*% tnsr

  return(r)
}

#' @rdname internal_sim_dat
#'
#' @param X Input data (list of samples or matrix where rows are samples)
#'
#' @return \code{rkhs_norm_mat}: Symmetric matrix of pairwise RKHS norms
#'
#' @examples
#' G <- diag(1,4)
#' fire:::rkhs_norm_mat(X = list(matrix(1:4,2), matrix(5:8,2)),
#'               Ginv = solve(G),
#'               Index = expand.grid(1:2, 1:2))
rkhs_norm_mat <- function(X, Ginv, Index) {
  # Input validation
  if (!(is.list(X) || is.matrix(X))) stop("X must be a list or matrix")

  # sample size
  if (is.list(X)) {
    N_sample <- length(X)
  } else {
    N_sample <- nrow(X)
  }

  # initialize a matrix to store the pairwise RKHS norm
  mat <- matrix(0, N_sample, N_sample)

  for (i in 1:N_sample) {
    for (j in i:N_sample) {
      if (is.list(X)) {
        mat[i,j] <- rkhs_norm(X[[i]], X[[j]], Ginv, Index)
      } else {
        mat[i,j] <- rkhs_norm(X[i,], X[j,], Ginv, Index)
      }
      mat[j,i] <- mat[i,j]
    }
  }

  return(mat)
}

#' @rdname internal_sim_dat
#'
#' @param Hurst Hurst coefficient in (0, 1]
#'
#' @return \code{cfbm_rkhs}: Centered Gram matrix for fractional Brownian motion in RKHS
#'
#' @examples
#' G <- diag(1,4)
#' fire:::cfbm_rkhs(X = list(matrix(1:4,2), matrix(5:8,2)),
#'           Ginv = solve(G),
#'           Index = expand.grid(1:2, 1:2),
#'           Hurst = 0.7)
cfbm_rkhs <- function(X, Ginv, Index, Hurst = 0.5) {
  # Input validation
  if (!is.numeric(Hurst) || Hurst <= 0 || Hurst > 1) {
    stop("Hurst coefficient must be in (0,1]")
  }

  if (is.list(X)) {
    N <- length(X)
  } else {
    N <- nrow(X)
  }

  # matrix that stores RKHS norm for each pair (X,X')
  mat <- rkhs_norm_mat(X, Ginv, Index)
  # power of Hurst
  mat <- mat^Hurst

  # Gram matrix
  G <- matrix(NA, N, N)

  fourth_term <- sum(mat)

  for (i in 1:N) {
    for (j in i:N) {
      G[i,j] <- -(N*N*(mat[i,j]) - N*sum(mat[i,]) - N*sum(mat[j,]) + fourth_term)/(2*N^2)
      G[j,i] <- G[i,j]
    }
  }

  return(G)
}

#' @rdname internal_sim_dat
#'
#' @param X_train Training data (list or matrix)
#' @param X_new New data (list or matrix, same type as X_train)
#'
#' @return \code{cfbm_rkhs_cross}: Cross-centered covariance matrix
#'
#' @examples
#' train <- list(matrix(1:4,2), matrix(5:8,2))
#' new <- list(matrix(9:12,2))
#' G <- diag(1,4)
#' fire:::cfbm_rkhs_cross(X_train = train,
#'                 X_new = new,
#'                 Ginv = solve(G),
#'                 Index = expand.grid(1:2, 1:2))
cfbm_rkhs_cross <- function(X_train, X_new, Ginv, Index, Hurst = 0.5) {

  if(is.list(X_train) != is.list(X_new)){
    stop("X_train and X_new must be of the same class, either list or matrix.")
  }
  # Determine if inputs are lists or matrices and get sizes
  if (is.list(X_train)) {
    N <- length(X_train)
    M <- length(X_new)
  } else {
    N <- nrow(X_train)
    M <- nrow(X_new)
  }

  # Initialize cross matrix
  cross_mat <- matrix(0, M, N)

  # Compute RKHS norms between new and training data
  for (i in 1:M) {
    for (j in 1:N) {
      if (is.list(X_train)) {
        norm_ij <- rkhs_norm(X_new[[i]], X_train[[j]], Ginv, Index)
      } else {
        norm_ij <- rkhs_norm(X_new[i,], X_train[j,], Ginv, Index)
      }
      cross_mat[i, j] <- norm_ij^Hurst
    }
  }

  # Compute row and column sums needed for centering
  train_mat <- rkhs_norm_mat(X_train, Ginv, Index)^Hurst

  train_row_sums <- rowSums(train_mat)
  train_total_sum <- sum(train_mat)

  # Compute the cross-centered matrix
  cross_centered <- matrix(0, M, N)

  for (i in 1:M) {
    new_row_sum <- sum(cross_mat[i,])
    for (j in 1:N) {
      cross_centered[i, j] <- -(N * cross_mat[i, j] - new_row_sum - train_row_sums[j] +
                                  train_total_sum) / (2 * N)
    }
  }

  return(cross_centered)
}

