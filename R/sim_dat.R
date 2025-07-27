#' Simulate Multi-dimensional Functional Data with I-prior
#'
#' @description
#' Generates simulated data for multi-dimensional functional regression models using
#' I-prior methodology.
#'
#' @param N Total sample size
#' @param Ntrain Training sample size.
#' @param kernels List of kernel functions for each mode
#' @param kernels_params List of parameters for each kernel function
#' @param alpha Vector of scale parameters for \code{kernels}
#' @param dat_T List of index sets for each mode
#' @param control A list of control parameters (see Details)
#'
#' @details The \code{control} argument can include the following parameters:
#' \itemize{
#'   \item \code{tau}: Scale parameter for the I-prior kernel (default: 1)
#'   \item \code{intercept_y}: Intercept term for response (default: 0)
#'   \item \code{intercept_x}: Intercept term for covariates (default: 0)
#'   \item \code{kernel_iprior}: Type of I-prior kernel ('cfbm', 'rbf', 'linear' or 'poly') (default: 'cfbm')
#'   \item \code{iprior_param}: Parameter for I-prior kernel (Hurst for cfbm, lengthscale for rbf) (default: NULL)
#'   \item \code{sigma_v}: Standard deviation for random effects (default: 1)
#'   \item \code{sigma}: Noise standard deviation (default: 1)
#'   \item \code{sigma_w}: Standard deviation for weights (default: NULL, which sets it to 1/sigma)
#'   \item \code{constant_g}: Logical indicating whether to include constant kernel term in g (default: TRUE)
#'   \item \code{constant_h}: Logical indicating whether to include constant kernel term in h (default: FALSE)
#'   \item \code{center}: Whether to center the kernel matrices (default: FALSE)
#'   \item \code{os_type}: Operating system type ('Apple' or 'Windows') (default: 'Apple')
#'   \item \code{cores}: Number of cores used in parallel computation (default: NULL)
#' }
#'
#' @return A list containing:
#' \itemize{
#' \item X: Simulated covariate data
#' \item y: Simulated response vector
#' }
#'
#' @examples
#' # 2D example
#' set.seed(1)
#' dat_T <- list(1:3, 1:2)
#' sim_dat(kernels = list(cfbm, rbf), kernels_params = list(0.5, 1),
#'         alpha = c(0.5, 0.5), dat_T = dat_T, N = 10, Ntrain = 5)
#'
#' # 3D example with control parameters
#' set.seed(1)
#' dat_T <- list(1:3, 1:2, 1:4)
#' sim_dat(kernels = list(cfbm, cfbm, cfbm), kernels_params = list(0.5, 0.5,0.5),
#'         alpha = c(0.5, 0.5, 0.5), dat_T = dat_T, N = 10, Ntrain = 5,
#'         control = list(kernel_iprior = 'rbf', sigma = 0.5))
#'
#' @export
sim_dat <- function(N, Ntrain,
                    kernels, kernels_params, alpha,
                    dat_T, control = list()){

  # Set default control parameters
  con <- list(
    tau = 1,
    intercept_y = 0,
    intercept_x = 0,
    kernel_iprior = 'cfbm',
    iprior_param = NULL,
    sigma_v = 1,
    sigma = 1,
    sigma_w = NULL,
    constant_g = TRUE,
    constant_h = FALSE,
    center = FALSE,
    os_type = 'Apple',
    cores = NULL
  )

  # Override defaults with user-supplied control parameters
  con[names(control)] <- control

  # Extract control parameters
  tau <- con$tau
  intercept_y <- con$intercept_y
  intercept_x <- con$intercept_x
  kernel_iprior <- con$kernel_iprior
  iprior_param <- con$iprior_param
  sigma_v <- con$sigma_v
  sigma <- con$sigma
  sigma_w <- con$sigma_w
  constant_g <- con$constant_g
  constant_h <- con$constant_h
  center <- con$center
  os_type <- con$os_type
  cores <- con$cores

  # Rest of the function remains the same...
  m <- length(kernels)
  if(m != length(kernels_params)) stop("kernels and kernels_params must have same length")
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

  Gmat.list <- gmat(kernels = kernels, kernels_params = kernels_params, dat = dat_T, center = center)

  # Construct G matrix for arbitrary m dimensions
  if(m == 1){
  G <- if(constant_g) matrix(1, ncol = d[1], nrow = d[1]) + alpha[1]^2 * Gmat.list[[1]]
  else alpha[1]^2 * Gmat.list[[1]]
  }else{
    G <- if(constant_g) matrix(1, ncol = d[1], nrow = d[1]) + alpha[1]^2 * Gmat.list[[1]]
    else alpha[1]^2 * Gmat.list[[1]]
    for(j in 2:m) {
      if(constant_g) {
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
    nmat <- Kronecker_norm_mat(X = X[1:Ntrain,], G = Gmat.list, alpha = alpha, Index = Index.T, constant = constant_g, os_type = os_type, sample_id = 1)
    nmat.cross <- Kronecker_norm_cross(Xtrain = X[1:Ntrain,], Xnew = X[-c(1:Ntrain),], G = Gmat.list, alpha = alpha, Index = Index.T, constant = constant_g, os_type = os_type, sample_id = 1)
  }else{
    X.list <- tensor_sample(X, sample_id = 1)
    X.train <- X.list[1:Ntrain]
    X.test <- X.list[-c(1:Ntrain)]
    nmat <- Kronecker_norm_mat(X = X.train, G = Gmat.list, alpha = alpha, Index = Index.T, constant = constant_g, os_type = os_type, sample_id = 1)
    nmat.cross <- Kronecker_norm_cross(Xtrain = X.train, Xnew = X.test, G = Gmat.list, alpha = alpha, Index = Index.T, constant = constant_g, os_type = os_type, sample_id = 1)
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

  if(constant_h){
    H <- 1 + H
    Hcross <- 1 + Hcross
  }
  w <- rnorm(Ntrain, mean = 0, sd = sigma_w)
  ## generate f(X_n) = Hw, n = 1,...,N
  Y.train <- as.vector(H %*% w)
  Y.test <- as.vector(Hcross %*% w)
  Y <- c(Y.train, Y.test)
  Y <- intercept_y + Y + rnorm(N, sd = sigma)

  return(list(X = X, y = Y))

}
