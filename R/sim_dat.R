#' Simulate Multi-dimensional Functional Data with I-prior
#'
#' Generates simulated data for multi-dimensional functional regression models using
#' I-prior methodology.
#'
#' @param kernels List of kernel functions for each mode, may refer to \code{\link{kernels_fire}}
#' @param kernels_param List of parameters for each kernel function
#' @param alpha Vector of scale parameters for \code{kernels}
#' @param dat_T List of index sets for each mode
#' @param tau Scale parameter for the I-prior kernel
#' @param intercept_y Intercept term for response
#' @param intercept_x Intercept term for covariates
#' @param kernel_iprior Type of I-prior kernel ('cfbm', 'rbf', 'linear' or 'poly')
#' @param iprior_param Parameter for I-prior kernel (Hurst for cfbm, lengthscale for rbf)
#' @param sigma_v Standard deviation for random effects
#' @param sigma Noise standard deviation
#' @param sigma_w Standard deviation for weights
#' @param constant Logical indicating whether to include constant kernel term
#' @param center Whether to center the kernel matrices
#' @param N Total sample size
#' @param Ntrain Training sample size
#' @param os_type Operating system type ('Apple' or 'Windows')
#' @param cores Number of cores used in parallel computation
#'
#' @return A list containing:
#' \itemize{
#' \item X: Simulated covariate data
#' \item y: Simulated response vector
#' }
#'
#' @seealso \code{\link{kernels_fire}}
#'
#' @examples
#' # 2D example
#' set.seed(1)
#' dat_T <- list(1:3, 1:2)
#' sim_dat(kernels = list(cfbm, rbf), kernels_param = list(0.5, 1),
#' alpha = c(0.5, 0.5), dat_T = dat_T, N = 10, Ntrain = 5, os_type = 'Apple', cores = 1)
#'
#' # 3D example
#' set.seed(1)
#' dat_T <- list(1:3, 1:2, 1:4)
#' sim_dat(kernels = list(cfbm, cfbm, cfbm), kernels_param = list(0.5, 0.5,0.5),
#' alpha = c(0.5, 0.5, 0.5), dat_T = dat_T, N = 10, Ntrain = 5,
#' kernel_iprior = 'rbf', os_type = 'Apple', cores = 1)
#'
#' @export
sim_dat <- function(kernels, kernels_param, alpha,
                     dat_T, tau = 1, intercept_y = 0, intercept_x = 0, kernel_iprior = 'cfbm', iprior_param = NULL,
                     sigma_v = 1, sigma = 1, sigma_w = NULL, constant = TRUE, center = FALSE,
                     N, Ntrain, os_type = 'Apple', cores = NULL){
  # kernel: kernel function h1, h2
  # kernels_param: parameters of kernel;cfbm - Hurst; rbf - lengthscale
  # alpha: alpha1, alpha2, alpha3
  # dat_T: index of set T, eg: list(1:I1, 1:I2, 1:I3), or list(seq(0,1,length.out = I1))
  # tau: scale parameter in g(X,X')
  # intercept_y: add intercept on y
  # intercept_x: randomly generate the intercept using standard normal distribution for x(t)
  # sigma_w: sd of normal dist for w_j
  # kernel_iprior: kernel in iprior
  # iprior_param: parameter of kernel in iprior
  # sigma: noise parameter
  # sigma_v: sd of normal dist for v_j
  # constant: include the constant kernel (i.e. 1), default is yes
  # N: sample size
  # Ntrain: sample size of training set

  m <- length(kernels)
  if(m != length(kernels_param)) stop("kernels and kernels_param must have same length")
  if(m != length(alpha)) stop("kernels and alpha must have same length")
  if(m != length(dat_T)) stop("kernels and dat_T must have same length")
  if(Ntrain >= N) stop("Ntrain must be smaller than N")

  if(is.null(sigma_w)){
    sigma_w <- 1/sigma
  }

  # dimension of each mode
  d <- sapply(dat_T, length)
  dims <- prod(d)

  # create a matrix consists of all combinations of T1 x T2 x T3 ...
  if(m == 1){
    Index.T <- matrix(1:d, ncol = 1)
  }else{
    Index.T <- create_index_matrix(dat_T)
  }

  Gmat.list <- gmat(kernels = kernels, kernels_params = kernels_param, dat = dat_T, center = center)

  # Construct G matrix for arbitrary m dimensions
  if(m == 1){
  G <- if(constant) matrix(1, ncol = d[1], nrow = d[1]) + alpha[1]^2 * Gmat.list[[1]]
  else alpha[1]^2 * Gmat.list[[1]]
  }else{
    G <- if(constant) matrix(1, ncol = d[1], nrow = d[1]) + alpha[1]^2 * Gmat.list[[1]]
    else alpha[1]^2 * Gmat.list[[1]]
    for(j in 2:m) {
      if(constant) {
        G_j <- matrix(1, ncol = d[j], nrow = d[j]) + alpha[j]^2 * Gmat.list[[j]]
      } else {
        G_j <- alpha[j]^2 * Gmat.list[[j]]
      }
      G <- kronecker(X = G, Y = G_j)
    }
  }

  # initialize covariate of size N (sample size) * dims (dimension)
  X <- matrix(NA, nrow = N, ncol = dims)
  for(i in 1:N){
    ## randomly generate v_j
    v <- rnorm(dims, mean = 0, sd = sigma_v)

    ## generate X_i = Gv, i = 1,...,N

    X[i,] <- as.vector(intercept_x + G %*% v)
  }
  # reshape matrix X into 2nd- or 3rd-tensor
  if(m != 1){
    X <- reshape_tensor(Mat = X, dim = d, N = N, Index = Index.T)
  }

  if(m == 1){
    nmat <- Kronecker_norm_mat(X = X[1:Ntrain,], G = Gmat.list, alpha = alpha, Index = Index, constant = constant, os_type = os_type, sample_id = 1)
    nmat.cross <- Kronecker_norm_cross(Xtrain = X[1:Ntrain,], Xnew = X[-c(1:Ntrain),], G = Gmat.list, alpha = alpha, Index = Index, constant = constant, os_type = os_type, sample_id = 1)
  }else{
    X.list <- tensor_sample(X, sample_id = 1)
    X.train <- X.list[1:Ntrain]
    X.test <- X.list[-c(1:Ntrain)]
    nmat <- Kronecker_norm_mat(X = X.train, G = Gmat.list, alpha = alpha, Index = Index.T, constant = constant, os_type = os_type, sample_id = 1)
    nmat.cross <- Kronecker_norm_cross(Xtrain = X.train, Xnew = X.test, G = Gmat.list, alpha = alpha, Index = Index.T, constant = constant, os_type = os_type, sample_id = 1)
  }

  if (kernel_iprior == 'cfbm') {
    if (is.null(iprior_param)) {
      iprior_param <- 0.5
    }
    H <- tau^2 * cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
    Hcross <- tau^2 * cfbm_rkhs_kron_cross(nmat = nmat, nmat_cross = nmat.cross, Hurst = iprior_param)
  }else if(kernel_iprior == 'rbf') {
    if (is.null(iprior_param)) {
      iprior_param <- 1
    }
    H <- tau^2 * rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
    Hcross <- tau^2 * rbf_rkhs_kron_cross(nmat.cross, lengthscale = iprior_param)

  }else if(kernel_iprior == 'linear') {
    if (is.null(iprior_param)) {
      iprior_param <- 0
    }
    H <- tau^2 * (nmat + iprior_param)
    Hcross <- tau^2 * (nmat.cross + iprior_param)
  } else if(kernel_iprior == 'poly') {
    if (is.null(iprior_param)) {
      iprior_param <- c(2, 0)
    }
    H <- tau^2 * (nmat + iprior_param)
    Hcross <- tau^2 * ( (nmat.cross + iprior_param[2])^iprior_param[1])
  } else {
    stop(paste("Unsupported kernel_iprior:", kernel_iprior,
               "- must be 'cfbm', 'rbf', 'linear' or 'poly'"))
  }

  w <- rnorm(Ntrain, mean = 0, sd = sigma_w)
  ## generate f(X_n) = Hw, n = 1,...,N
  Y.train <- as.vector(H %*% w)
  Y.test <- as.vector(Hcross %*% w)
  Y <- c(Y.train, Y.test)
  Y <- intercept_y + Y + rnorm(N, sd = sigma)

  return(list(X = X, y = Y))

}
