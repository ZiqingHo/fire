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
#' @seealso \code{\link{gmat}}
#' @examples
#' \dontrun{
#' # 1D case without constant term
#' G1 <- list(matrix(c(2,1,1,2), 2))
#' res1 <- fire:::Kronecker_inv_helper(G1, alpha = 0.5, constant = FALSE)
#'
#' # 4D case with constant term
#' G4 <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2),
#'            matrix(c(4,1,1,4),2), matrix(c(5,1,1,5),2))
#' res4 <- fire:::Kronecker_inv_helper(G4, alpha = rep(0.5,4), constant = TRUE)
#' }
Kronecker_inv_helper <- function(G, alpha, constant = TRUE){

  m <- length(G)

  # Convert logical constant to numeric (1/0)
  const <- as.integer(constant)
  result <- list()

  # Special case: m = 1 and no constant term
  if (m == 1 && !constant) {
    result$G1 <- alpha[1]^2 * G[[1]]
    return(result)
  }

  # Fast path for common case (m ≤ 3)
  if(m <= 3) {
    G1 <- (G[[1]] + t(G[[1]]))/2
    G1.tilde <- const * tcrossprod(rep(1,nrow(G1))) + alpha[1]^2 * G1
    G1.eigen = eigen(G1.tilde)
    result$Q1 = G1.eigen$vectors
    result$L1 = G1.eigen$values

    if(m >= 2){

      G2 <- (G[[2]] + t(G[[2]]))/2
      G2.tilde <- const * tcrossprod(rep(1,nrow(G2))) + alpha[2]^2 * G2
      G2.eigen = eigen(G2.tilde)
      result$Q2 = G2.eigen$vectors
      result$L2 = G2.eigen$values

      if(m == 3){

        G3 <- (G[[3]] + t(G[[3]]))/2
        G3.tilde <- const * tcrossprod(rep(1,nrow(G3))) + alpha[3]^2 * G3
        G3.eigen = eigen(G3.tilde)
        result$Q3 = G3.eigen$vectors
        result$L3 = G3.eigen$values
      }
    }
  } else {
    # General case for arbitrary m
    for (i in 1:m) {
      G_i <- G[[i]]

      # avoid numerical issues
      G_i <- (G_i + t(G_i)) / 2

      # Create modified Gram matrix
      I <- nrow(G_i)
      G_tilde <- const * tcrossprod(rep(1, I)) + alpha[i]^2 * G_i

      # Compute eigendecomposition
      eigen_decomp <- eigen(G_tilde)

      # Store results
      result[[paste0("Q", i)]] <- eigen_decomp$vectors
      result[[paste0("L", i)]] <- eigen_decomp$values
    }

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
#' \dontrun{
#' L <- list(c(1,2), c(3,4))
#' inv <- fire:::Lambda_inv(L)
#' }
Lambda_inv <- function(L){

  m <- length(L)

  if(m == 1){
    L.inv <- 1/L[[1]]
  }else{
    L.diag <- Reduce(function(x, y) as.vector(aperm(outer(x, y))), L)
    L.inv <- 1/L.diag
  }

  if(any(is.finite(L.inv) == FALSE)){
    warning('At least one eigenvalue is 0.')
  }

  return(L.inv)
}

#' @rdname internal_kronecker
#'
#' @param Xvec Vectorized input tensor
#' @param Q Matrix from \code{Kronecker_inv_helper} output
#'
#' @examples
#' \dontrun{
#' result <- fire:::kron_mv(rnorm(8), matrix(rnorm(4),2))
#' }
kron_mv <- function(Xvec, Q){

  nQ <- ncol(Q)
  nX <- length(Xvec)

  # Reshape and multiply
  dim(Xvec) <- c(nQ, nX/nQ)
  Xvec <- t(Xvec) %*% Q
  dim(Xvec) <- NULL

  return(Xvec)
}
