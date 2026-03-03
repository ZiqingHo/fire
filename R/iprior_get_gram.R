#' Compute Gram Matrices for I-prior Models
#'
#' @description Functions to compute various Gram matrices for I-prior models using RKHS norms.
#'
#' @param nmat RKHS norm matrix of training data
#' @param nmat_cross Matrix of RKHS norm between training and test data
#' @param lengthscale Lengthscale parameter for RBF kernel
#' @param Hurst Hurst coefficient for cfbm kernel
#'
#' @return A Gram matrix of appropriate dimensions:
#' \itemize{
#'   \item For \code{rbf_rkhs_kron} and \code{cfbm_rkhs_kron}: Square training Gram matrix
#'   \item For \code{rbf_rkhs_kron_cross} and \code{cfbm_rkhs_kron_cross}: Cross Gram matrix between test and training data
#' }
#'
#' @examples
#' # Using internal functions to compute pairwise RKHS norm
#' G <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2))
#' X <- list(matrix(rnorm(4),2), matrix(rnorm(4),2))
#' nmat <- fire:::Kronecker_norm_mat(X, G, alpha=c(0.5,0.5))
#' ncross <- fire:::Kronecker_norm_cross(X, X[1], G, alpha=c(0.5,0.5))
#'
#' # Get Gram matrices
#' K_rbf <- rbf_rkhs_kron(nmat)
#' K_cfbm <- cfbm_rkhs_kron(nmat)
#' K_rbf_cross <- rbf_rkhs_kron_cross(ncross)
#' K_cfbm_cross <- cfbm_rkhs_kron_cross(nmat, ncross)
#'
#' @name iprior_get_gram
NULL

#' @rdname iprior_get_gram
#' @export
rbf_rkhs_kron <- function(nmat, lengthscale = 1) {
  lengthscale <- as.numeric(lengthscale)[1L]
  if (!is.finite(lengthscale) || lengthscale <= 0) {
    stop("lengthscale must be a positive finite numeric scalar.")
  }

  nmat <- as.matrix(nmat)
  if (!is.numeric(nmat)) {
    stop("nmat must contain numeric values.")
  }
  if (all(is.na(nmat))) {
    stop("nmat must contain numeric values (not all NA).")
  }

  exp(nmat / (-2 * lengthscale^2))
}

#' @rdname iprior_get_gram
#' @export
rbf_rkhs_kron_cross <- function(nmat_cross, lengthscale = 1) {
  lengthscale <- as.numeric(lengthscale)[1L]
  if (!is.finite(lengthscale) || lengthscale <= 0) {
    stop("lengthscale must be a positive finite numeric scalar.")
  }

  nmat_cross <- as.matrix(nmat_cross)
  if (!is.numeric(nmat_cross)) {
    stop("nmat_cross must contain numeric values.")
  }
  if (all(is.na(nmat_cross))) {
    stop("nmat_cross must contain numeric values (not all NA).")
  }

  exp(nmat_cross / (-2 * lengthscale^2))
}

#' @rdname iprior_get_gram
#' @export
cfbm_rkhs_kron <- function(nmat, Hurst = 0.5) {
  # Force Hurst to a numeric scalar
  Hurst <- as.numeric(Hurst)[1L]
  if (!is.finite(Hurst)) {
    stop("Hurst must be a finite numeric scalar.")
  }
  if (Hurst > 1 || Hurst <= 0) {
    stop("Hurst must be in (0, 1].")
  }

  # Force nmat to be numeric matrix
  nmat <- as.matrix(nmat)
  if (!is.numeric(nmat)) {
    stop("nmat must contain numeric values.")
  }
  if (all(is.na(nmat))) {
    stop("nmat must contain numeric values (not all NA).")
  }
  N <- nrow(nmat)
  mat <- nmat^Hurst
  I <- diag(N) - matrix(1/N, N, N)
  -0.5 * I %*% mat %*% I
}

#' @rdname iprior_get_gram
#' @export
cfbm_rkhs_kron_cross <- function(nmat, nmat_cross, Hurst = 0.5) {
  Hurst <- as.numeric(Hurst)[1L]
  if (!is.finite(Hurst)) {
    stop("Hurst must be a finite numeric scalar.")
  }
  if (Hurst > 1 || Hurst <= 0) {
    stop("Hurst must be in (0, 1].")
  }

  # Force nmat, nmat_cross to numeric matrices
  nmat <- as.matrix(nmat)
  nmat_cross <- as.matrix(nmat_cross)

  if (!is.numeric(nmat)) {
    stop("nmat must contain numeric values.")
  }
  if (!is.numeric(nmat_cross)) {
    stop("nmat_cross must contain numeric values.")
  }
  if (all(is.na(nmat))) {
    stop("nmat must contain numeric values (not all NA).")
  }
  if (all(is.na(nmat_cross))) {
    stop("nmat_cross must contain numeric values (not all NA).")
  }

  Ntrain <- nrow(nmat)
  Nnew   <- nrow(nmat_cross)

  norm_cross <- nmat_cross^Hurst
  norm_train <- nmat^Hurst

  rvec_new   <- rowSums(norm_cross)
  s_train    <- sum(norm_train)
  rvec_train <- colSums(norm_train)

  G_cross <- norm_cross -
    rvec_new / Ntrain -
    matrix(rvec_train / Ntrain, nrow = Nnew, ncol = Ntrain, byrow = TRUE) +
    s_train / (Ntrain^2)

  G_cross / (-2)
}
