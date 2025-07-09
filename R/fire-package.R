#' FIRE: Functional I-prior Regression for Vectors and Tensors
#'
#' @description
#' Implements functional I-prior regression models for:
#' \itemize{
#'   \item Multivariate vector covariates (n samples × p features)
#'   \item Multi-dimensional tensor covariates
#' }
#' using reproducing kernel Hilbert spaces (RKHS) with EM algorithm estimation.
#'
#' @details
#' The FIRE package provides two specialized implementations:
#'
#' \strong{1. Vector Covariates} (via \code{\link{fire.matrix}}):
#' - For design matrices X of dimension n × p
#' - Each row represents one sample's p-dimensional covariate vector
#' - Kernel operates on the p-dimensional feature space
#'
#' \strong{2. Tensor Covariates} (via \code{\link{fire.tensor}}):
#' - For array/list inputs representing n samples of tensor data
#' - Supports arbitrary-order tensors (2D matrices, 3D arrays, etc.)
#' - Utilizes Kronecker product kernels
#'
#' @section S3 Methods:
#' The main \code{fire()} dispatches to:
#' \describe{
#'   \item{\code{fire.matrix}}{For n×p design matrices (vector covariates)}
#'   \item{\code{fire.array}, \code{fire.list}}{For tensor covariates}
#' }
#'
#' @author
#' \strong{Maintainer}: Ziqing Ho \email{hoziqing@hotmail.com}
#'
#' Contributors:
#' \itemize{
#'   \item Additional contributors (if any)
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/yourusername/FIRE}
#'   \item Report bugs at \url{https://github.com/yourusername/FIRE/issues}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
