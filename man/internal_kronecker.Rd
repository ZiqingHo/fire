% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz_internal_kronecker_rkhs.R
\name{internal_kronecker}
\alias{internal_kronecker}
\alias{Kronecker_inv_helper}
\alias{Lambda_inv}
\alias{kron_mv}
\title{Internal Kronecker RKHS Norm Utilities}
\usage{
Kronecker_inv_helper(G, alpha, constant = TRUE)

Lambda_inv(L)

kron_mv(Xvec, Q)
}
\arguments{
\item{G}{List of Gram matrices from \code{\link{gmat}}}

\item{alpha}{Numeric vector of scale parameters}

\item{constant}{Logical indicating whether to include constant kernel term}

\item{L}{List of eigenvalue vectors from \code{Kronecker_inv_helper}}

\item{Xvec}{Vectorized input tensor}

\item{Q}{Matrix from \code{Kronecker_inv_helper} output}
}
\value{
For each function:
\itemize{
  \item{\code{Kronecker_inv_helper}: Returns a list with:
    \itemize{
      \item{\code{G1} when m=1 and constant=FALSE (modified Gram matrix)}
      \item{\code{Q1,L1,Q2,L2,...} otherwise (eigenvectors and eigenvalues for each mode)}
    }
  }
  \item{\code{Lambda_inv}: Returns numeric vector of inverse eigenvalues}
  \item{\code{kron_mv}: Returns numeric vector of matrix-vector product}
}
}
\description{
Internal functions for efficient RKHS norm computations using Kronecker structure:
\itemize{
\item{\code{Kronecker_inv_helper}: Precompute eigendecomposition for Kronecker-structured Gram matrices.
  Handles arbitrary number of modes with optional constant kernel terms.}
\item{\code{Lambda_inv}: Efficiently compute inverse of Kronecker-structured eigenvalue matrix.}
\item{\code{kron_mv}: Efficient Kronecker matrix-vector multiplication.}
}
}
\details{
\code{Kronecker_inv_helper}: Special case avoids eigendecomposition when m = 1 and constant = FALSE for efficiency.

\code{Lambda_inv}:
Computes \eqn{\Lambda^{-1}} where \eqn{\Lambda = \Lambda_1 \otimes \Lambda_2 \otimes ... \otimes \Lambda_m}
}
\examples{
# Kronecker_inv_helper()
# 1D case without constant term
G1 <- list(matrix(c(2,1,1,2), 2))
res1 <- fire:::Kronecker_inv_helper(G1, alpha = 0.5, constant = FALSE)

# 4D case with constant term
G4 <- list(matrix(c(2,1,1,2),2), matrix(c(3,1,1,3),2),
           matrix(c(4,1,1,4),2), matrix(c(5,1,1,5),2))
res4 <- fire:::Kronecker_inv_helper(G4, alpha = rep(0.5,4), constant = TRUE)
# Lambda_inv()
L <- list(c(1,2), c(3,4))
inv <- fire:::Lambda_inv(L)
# kron_mv()
result <- fire:::kron_mv(rnorm(8), matrix(rnorm(4),2))
}
\seealso{
\code{\link{gmat}}, \code{\link{fire}}
}
\keyword{internal}
