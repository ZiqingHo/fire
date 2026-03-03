#' @title Internal Kronecker RKHS Norm Utilities
#' @name internal_kronecker
#' @description Internal functions for efficient RKHS norm computations using Kronecker structure:
#' \itemize{
#' \item{\code{Kronecker_inv_helper}: Precompute eigendecomposition for Kronecker-structured Gram matrices.
#'   Handles arbitrary number of modes with optional constant kernel terms.}
#' \item{\code{Lambda_inv}: Efficiently compute inverse of Kronecker-structured eigenvalue matrix.}
#' \item{\code{kron_mv}: Efficient Kronecker matrix-vector multiplication.}
#' }
#'
#' @return For each function:
#' \itemize{
#'   \item{\code{Kronecker_inv_helper}: Returns a list with:
#'     \itemize{
#'       \item{\code{G1} when m=1 and constant=FALSE (modified Gram matrix)}
#'       \item{\code{Q1,L1,Q2,L2,...} otherwise (eigenvectors and eigenvalues for each mode)}
#'     }
#'   }
#'   \item{\code{Lambda_inv}: Returns numeric vector of inverse eigenvalues}
#'   \item{\code{kron_mv}: Returns numeric vector of matrix-vector product}
#' }
#' @keywords internal
NULL

#' @rdname internal_kronecker
#'
#' @param G List of Gram matrices from \code{\link{gmat}}
#' @param alpha Numeric vector of scale parameters
#' @param constant Logical indicating whether to include constant kernel term
#'
#' @details
#' \code{Kronecker_inv_helper}: Special case avoids eigendecomposition when m = 1 and constant = FALSE for efficiency.
#'
#' @seealso \code{\link{gmat}}, \code{\link{fire}}
#' @examples
#' # Kronecker_inv_helper()
#' # 1D case without constant term
#' G1 <- list(matrix(c(2,1,1,2), 2))
#' res1 <- fire:::Kronecker_inv_helper(G1, alpha = 0.5, constant = FALSE)
#'
#' # 4D case with constant term
#' G4 <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2),
#'            matrix(c(4,1,1,4),2), matrix(c(5,1,1,5),2))
#' res4 <- fire:::Kronecker_inv_helper(G4, alpha = rep(0.5,4), constant = TRUE)
Kronecker_inv_helper <- function(G, alpha, constant = TRUE, tol = 1e-6) {
  m <- length(G)
  const <- as.integer(constant)
  result <- list()

  # Special case: m = 1 and no constant term
  if (m == 1 && !constant) {
    result$G1 <- alpha[1]^2 * G[[1]]
    return(result)
  }

  eig_with_fix <- function(G_i, alpha_i) {
    G_i <- (G_i + t(G_i)) / 2
    n <- nrow(G_i)

    G_tilde <- const * tcrossprod(rep(1, n)) + alpha_i^2 * G_i
    eg <- eigen(G_tilde)

    # force PSD
    if (any(eg$values < 0)) {
      neg_vals <- eg$values[eg$values < -tol]
      if (length(neg_vals) > 0) {
        warning("Some eigenvalues are significantly negative; kernel may not be PSD.")
      }
      eg$values[eg$values < tol] <- 0
    }
    eg
  }

  for (i in seq_len(m)) {
    eg <- eig_with_fix(G[[i]], alpha[i])
    result[[paste0("Q", i)]] <- eg$vectors
    result[[paste0("L", i)]] <- eg$values
  }

  return(result)
}

#' @rdname internal_kronecker
#'
#' @param L List of eigenvalue vectors from \code{Kronecker_inv_helper}
#'
#' @details  \code{Lambda_inv}:
#' Computes \eqn{\Lambda^{-1}} where \eqn{\Lambda = \Lambda_1 \otimes \Lambda_2 \otimes ... \otimes \Lambda_m}
#' @examples
#' # Lambda_inv()
#' L <- list(c(1,2), c(3,4))
#' inv <- fire:::Lambda_inv(L)
Lambda_inv <- function(L) {
  m <- length(L)

  if (m == 1) {
    lam <- L[[1]]
  } else {
    lam <- Reduce(function(x, y) as.vector(aperm(outer(x, y))), L)
  }
  # relative tolerance: only keep reasonably large eigenvalues
  rel_tol <- 1e-12
  tol <- rel_tol * max(abs(lam))
  nz <- abs(lam) > tol

  L_inv <- numeric(length(lam))
  L_inv[nz] <- 1 / lam[nz]

#  if (!all(nz)) warning("Some eigenvalues treated as zero; using pseudo-inverse.")
  L_inv
}
#' @rdname internal_kronecker
#'
#' @param Xvec Vectorized input tensor
#' @param Q Matrix from \code{Kronecker_inv_helper} output
#'
#' @examples
#' # kron_mv()
#' result <- fire:::kron_mv(rnorm(8), matrix(rnorm(4),2))
kron_mv <- function(Xvec, Q){

  nQ <- ncol(Q)
  nX <- length(Xvec)

  # Reshape and multiply
  dim(Xvec) <- c(nQ, nX/nQ)
  Xvec <- t(Xvec) %*% Q
  dim(Xvec) <- NULL

  return(Xvec)
}
