#' Kernel Functions for FIRE
#'
#' A collection of kernel functions including:
#' \itemize{
#'   \item{\code{fbm}: Fractional Brownian motion kernel}
#'   \item{\code{cfbm}: Centered fractional Brownian motion kernel}
#'   \item{\code{cfbm_sd}: Standardized cfbm kernel}
#'   \item{\code{kronecker_delta}: Kronecker delta (identity) kernel}
#'   \item{\code{rbf}: Radial basis function (squared exponential) kernel}
#'   \item{\code{polynomial}: Polynomial kernel}
#'   \item{\code{mercer}: Mercer kernel with cosine basis}
#'   \item{\code{matern}: Matern kernel}
#' }
#'
#' @name kernels_fire
#' @aliases fbm cfbm cfbm_sd kronecker_delta rbf polynomial mercer
NULL

#' @rdname kernels_fire
#' @param X Input data (vector or matrix)
#' @param Hurst Hurst parameter for fbm/cfbm (between 0 and 1)
#' @param std Logical indicating whether to standardise the kernel matrix
#' @return A symmetric positive definite Gram matrix of size n x n where n is the
#'   number of observations in X. The matrix has attributes including:
#'   \itemize{
#'     \item{\code{kernel} - the type of kernel used}
#'     \item{\code{parameters} - the kernel parameters used}
#'   }
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(4), ncol=2)
#'
#' # Different kernels
#' fbm(X)
fbm <- function(X, Hurst = 0.5, std = FALSE) {
  # convert vector into matrix if necessary
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  if(Hurst >1 || Hurst <= 0){
    stop('Hurst must be (0,1]')
  }

  if(all(is.na(X))){
    stop('X must contains numeric values.')
  }

  # sample size
  N <- nrow(X)

  # Precompute the squared norms of each sample
  squared_norms <- rowSums(X^2)

  # Initialize the matrix A to store pairwise squared Euclidean distances ||xi - xj||^2*Hurst
  A <- matrix(0, N, N)

  # Get the indices of the upper triangular part of A
  index.mat <- upper.tri(A)
  index <- which(index.mat, arr.ind = TRUE)

  # Compute the cross-product matrix X %*% t(X)
  xcrossprod <- tcrossprod(X)

  # Extract the diagonal elements of xcrossprod for the row and column indices
  tmp1 <- diag(xcrossprod)[index[, 1]]  # ||x_i||^2
  tmp2 <- diag(xcrossprod)[index[, 2]]  # ||x_j||^2
  tmp3 <- xcrossprod[index]             # <x_i, x_j>

  # Compute the pairwise squared Euclidean distances
  A[index.mat] <- tmp1 + tmp2 - 2 * tmp3  # ||x_i - x_j||^2
  A <- A + t(A)  # Make A symmetric (diagonal is already zero)

  # Apply the Hurst exponent to the squared Euclidean distances
  A <- A^(Hurst)

  # each row is the squared norms raised to the Hurst exponent ||xi||^2*Hurst
  B <- matrix(rep(squared_norms^Hurst, N), nrow = N)

  # Gram matrix
  K <- -0.5 * (A - B - t(B))

  if(std){
    ones <- matrix(1, N, N)
    K_centered <- K - (ones %*% K)/N - (K %*% ones)/N + (ones %*% K %*% ones)/(N^2)
    K_std <- K_centered/(mean(diag(K_centered)))
    attr(K_std, "kernel") <- "fbm (standardised)"
    attr(K_std, "parameters") <- list(Hurst = Hurst)
    K_std
  }else{
    attr(K, "kernel") <- "fbm"
    attr(K, "parameters") <- list(Hurst = Hurst)
    K
  }

}


#' @rdname kernels_fire
#' @examples
#' cfbm(X)
#' @export
cfbm <- function(X, Hurst = 0.5, std = FALSE) {
  if (is.matrix(X)){
    N <- dim(X)[1]
  }else{
    N <- length(X)
  }
  if(Hurst >1 || Hurst <= 0){
    stop('Hurst must be (0,1]')
  }
  if(all(is.na(X))){
    stop('X must contains at least one numeric value.')
  }

  A <- matrix(0, N, N) # matrix of zeroes
  index.mat <- upper.tri(A) # entry above diagonal is TRUE o/w FALSE
  index <- which(index.mat, arr.ind = TRUE) # get the indices which entry is TRUE along column
  xcrossprod <- tcrossprod(X) # x %*% t(x)
  tmp1 <- diag(xcrossprod)[index[, 1]] # extract the diagonal elements following the row index
  tmp2 <- diag(xcrossprod)[index[, 2]] # along column index
  tmp3 <- xcrossprod[index] # extract x%*%t(x)[i,j]
  ## A is a matrix that each entry is the pairwise squared Euclidean distance
  # ||xi - xj||^2 = ||xi||^2 + ||xj||^2 - 2<xi,xj>
  # = Xcp[i,i] + Xcp[j,j] - 2Xcp[i,j]
  A[index.mat] <- tmp1 + tmp2 - 2 * tmp3 # the upper triangle is ||xi-xj||^2
  A <- A + t(A) # A is symmetric, diagonal is zero
  A <- abs(A) ^ Hurst # 1st term
  s <- sum(A) # 4th term; all the same regardless of i, j
  rvec <- colSums(A) # i-th entry is the sum_j ||xi-xj||^2
  rvec1 <- tcrossprod(rvec, rep(1, N)) # a matrix whose (i,k) entry is sum_j ||xi-xj||^2
  K <- (N^2 * A - N * rvec1 - N * t(rvec1) + s) / (-2* N^2)

  if(std){
    K_std <- K/(mean(diag(K)))
    attr(K_std, "kernel") <- "cfbm (standardised)"
    attr(K_std, "parameters") <- list(Hurst = Hurst)
    K_std
  }else{
    attr(K, "kernel") <- "cfbm"
    attr(K, "parameters") <- list(Hurst = Hurst)
    K
  }

}
#' @rdname kernels_fire
#' @param center Logical indicating whether to center the kernel matrix
#' @examples
#' kronecker_delta(X)
#' @export
kronecker_delta <- function(X, center = FALSE, std = FALSE){

  # sample size
  N <- if(is.matrix(X)) dim(X)[1] else length(X)

  # initialize Gram matrix
  K <- matrix(0, N, N)

  # kronecker delta kernel for a pair of (x,y)
  delta <- function(x,y){
    r <- as.numeric(identical(x,y))
    return(r)
  }

  for(i in 1:N){
    for(j in i:N){
      if(is.matrix(X)){
        K[i,j] <- delta(X[i,], X[j,])
        K[j,i] <- K[i,j]
      }else{
        K[i,j] <- delta(X[i], X[j])
        K[j,i] <- K[i,j]
      }
    }
  }

  # Center and/or standardize if requested
  if (center || std) {
    ones <- matrix(1, N, N)
    K_centered <- K - (ones %*% K)/N - (K %*% ones)/N + (ones %*% K %*% ones)/(N^2)

    if (std) {
      # Standardize by the average of the diagonal elements
      K_std <- K_centered / mean(diag(K_centered))
      attr(K_std, "kernel") <- "kronecker delta (standardised)"
      K_std
    } else {
      attr(K_centered, "kernel") <- "kronecker delta (centered)"
      K_centered
    }
  }else{
    attr(K, "kernel") <- "kronecker_delta"
    K
  }
}


#' @rdname kernels_fire
#' @param lengthscale Bandwidth parameter
#' @examples
#' rbf(X)
#' @export
rbf <- function(X, lengthscale = 1, center = FALSE, std = FALSE){
  if(all(is.na(X))){
    stop('X must contains numeric values.')
  }
  # sample size
  N <- if(is.matrix(X)) dim(X)[1] else length(X)

  # initialize Gram matrix
  K <- matrix(0, N, N)

  rbf_single <- function(X, Y, lengthscale = 1){
    r <- exp(- (sum((X-Y)^2))/(2 * lengthscale^2))
    return(r)
  }

  for(i in 1:N){
    for(j in i:N){
      if(is.matrix(X)){
        K[i,j] <- rbf_single(X[i,], X[j,], lengthscale)
        K[j,i] <- K[i,j]
      }else{
        K[i,j] <- rbf_single(X[i], X[j], lengthscale)
        K[j,i] <- K[i,j]
      }
    }
  }

  # Center and/or standardize if requested
  if (center || std) {
    ones <- matrix(1, N, N)
    K_centered <- K - (ones %*% K)/N - (K %*% ones)/N + (ones %*% K %*% ones)/(N^2)

    if (std) {
      # Standardize by the average of the diagonal elements
      K_std <- K_centered / mean(diag(K_centered))
      attr(K_std, "kernel") <- "rbf (standardised)"
      attr(K_std, "parameters") <- list(lengthscale = lengthscale)
      K_std
    } else {
      attr(K_centered, "kernel") <- "rbf (centered)"
      attr(K_centered, "parameters") <- list(lengthscale = lengthscale)
      K_centered
    }
  }else{
    attr(K, "kernel") <- "rbf"
    attr(K, "parameters") <- list(lengthscale = lengthscale)
    K
  }

}

#' @rdname kernels_fire
#' @param d Degree of polynomial
#' @param offset Constant offset in polynomial kernel
#' @examples
#' polynomial(X, d = 2, offset = 1)
#' @export
polynomial <- function(X, d = 1, offset = 0, center = FALSE, std = FALSE){
  if(all(is.na(X))){
    stop('X must contains at least one numeric value.')
  }
  # Sample size
  N <- if(is.matrix(X)) dim(X)[1] else length(X)

  # Initialize Gram matrix
  K <- matrix(0, N, N)

  polynomial_single <- function(X, Y, d = 1, offset = 0){
    # d: degree
    # offset: c
    return((as.numeric(X %*% Y) + offset)^d)
  }

  for(i in 1:N){
    for(j in i:N){
      if(is.matrix(X)){
        K[i,j] <- polynomial_single(X[i,], X[j,], d, offset)
        K[j,i] <- K[i,j]
      }else{
        K[i,j] <- polynomial_single(X[i], X[j], d, offset)
        K[j,i] <- K[i,j]
      }
    }
  }

  # Center and/or standardize if requested
  if (center || std) {
    ones <- matrix(1, N, N)
    K_centered <- K - (ones %*% K)/N - (K %*% ones)/N + (ones %*% K %*% ones)/(N^2)

    if (std) {
      # Standardize by the average of the diagonal elements
      K_std <- K_centered / mean(diag(K_centered))
      attr(K_std, "kernel") <- "polynomial (standardised)"
      attr(K_std, "parameters") <- list(d = d, offset = offset)
      K_std
    } else {
      attr(K_centered, "kernel") <- "polynomial (centered)"
      attr(K_centered, "parameters") <- list(d = d, offset = offset)
      K_centered
    }
  }else{
    attr(K, "kernel") <- "polynomial"
    attr(K, "parameters") <- list(d = d, offset = offset)
    K
  }
}

#' @rdname kernels_fire
#' @param delta Smoothness parameter for Mercer kernel
#' @param max_terms Maximum number of terms in Mercer series expansion
#' @examples
#' mercer(X)
#' @export
mercer <- function(X, delta = 1, max_terms = 1000, center = FALSE, std = FALSE) {
  if(all(is.na(X))){
    stop('X must contains at least one numeric value.')
  }

  N <- length(X)

  if (!is.numeric(X)) {
    stop("X must be a numeric vector with all elements between 0 and 1.")
  }
  if (any(X < 0 | X > 1)) {
    cat("Remark: X is converted into a numeric vector with all elements between 0 and 1.\n")
    X = seq(0, 1, length.out = N)
  }
  if(delta <= 0){
    stop('The delta parameter must be positive.')
  }


  K <- matrix(NA, N, N)

  # g(x,x') = sum_i^infinity i^-(1+delta) cos(pi i x) cos(pi i x')
  mercer_single <- function(X, Y, delta, max_terms = 1000) {
    if (delta <= 0) stop("delta must be positive")

    # Create vector of all terms to compute
    i <- 1:max_terms
    weights <- 1 / (i^(1 + delta))

    # Use trigonometric identity: cos(a)*cos(b) = 0.5[cos(a+b) + cos(a-b)]
    angle_sum <- pi * i %o% (X + Y)
    angle_diff <- pi * i %o% (X - Y)

    # Compute all terms at once
    terms <- 0.5 * weights * (cos(angle_sum) + cos(angle_diff))

    # Sum over all terms (rowSums if X/Y are vectors)
    if (length(X) == 1 && length(Y) == 1) {
      sum(terms)
    } else {
      rowSums(terms)
    }
  }

  for(i in 1:N){
    for(j in i:N){
      K[i,j] <- mercer_single(X[i],X[j],delta, max_terms)
      K[j,i] <- K[i,j]
    }
  }

  # Center and/or standardize if requested
  if (center || std) {
    ones <- matrix(1, N, N)
    K_centered <- K - (ones %*% K)/N - (K %*% ones)/N + (ones %*% K %*% ones)/(N^2)

    if (std) {
      # Standardize by the average of the diagonal elements
      K_std <- K_centered / mean(diag(K_centered))
      attr(K_std, "kernel") <- "mercer (standardised)"
      attr(K_std, "parameters") <- list(delta = delta)
      K_std
    } else {
      attr(K_centered, "kernel") <- "mercer (centered)"
      attr(K_centered, "parameters") <- list(delta = delta)
      K_centered
    }
  }else{
    attr(K, "kernel") <- "mercer"
    attr(K, "parameters") <- list(delta = delta)
    K
  }
}

#' @rdname kernels_fire
#' @param nu Smoothness parameter (typical values: 0.5, 1.5, 2.5)
#' @param sigma Scale parameter controlling the overall variance in Matern kernel
#' @examples
#' matern(X)
#' @export
matern<- function(X, nu = 1.5, lengthscale = 1.0, sigma = 1.0, center = FALSE, std = FALSE) {
  N <- if (is.matrix(X)) dim(X)[1] else length(X)

  # Compute pairwise Euclidean distances
  dists <- as.matrix(dist(X))
  r <- dists / lengthscale

  # Initialize kernel matrix
  K <- matrix(0, nrow = N, ncol = N)

  # Handle the r = 0 case separately to avoid NaNs in Bessel function
  r[r == 0] <- 1e-10


  # Handle special cases for common nu values (faster computation)
  if (nu == 0.5) {
    # Matern 1/2 kernel (equivalent to exponential kernel)
    K <- sigma^2 * exp(-r)
  } else if (nu == 1.5) {
    # Matern 3/2 kernel
    sqrt3_r <- sqrt(3) * r
    K <- sigma^2 * (1 + sqrt3_r) * exp(-sqrt3_r)
  } else if (nu == 2.5) {
    # Matern 5/2 kernel
    sqrt5_r <- sqrt(5) * r
    K <- sigma^2 * (1 + sqrt5_r + (5 * r^2)/3) * exp(-sqrt5_r)
  } else {
    # General case using Bessel function
    r[r == 0] <- 1e-10  # Avoid NaNs in Bessel function
    const <- sigma^2 * (2^(1 - nu)) / gamma(nu)
    scaled_r <- sqrt(2 * nu) * r
    K <- const * (scaled_r^nu) * besselK(scaled_r, nu)
  }
  # Set the diagonal (originally r == 0) to sigma^2
  diag(K) <- sigma^2

  # Center and/or standardize if requested
  if (center || std) {
    ones <- matrix(1, N, N)
    K_centered <- K - (ones %*% K)/N - (K %*% ones)/N + (ones %*% K %*% ones)/(N^2)

    if (std) {
      # Standardize by the average of the diagonal elements
      K_std <- K_centered / mean(diag(K_centered))
      attr(K_std, "kernel") <- "mercer (standardised)"
      attr(K_std, "parameters") <- list(nu = nu, lengthscale = lengthscale, sigma = sigma)
      K_std
       } else {
      attr(K_centered, "kernel") <- "mercer (centered)"
      attr(K_centered, "parameters") <- list(nu = nu, lengthscale = lengthscale, sigma = sigma)
      K_centered
      }
  } else {
    attr(K, "kernel") <- "matern"
    attr(K, "parameters") <- list(nu = nu, lengthscale = lengthscale, sigma = sigma)
    return(K)
  }
}
