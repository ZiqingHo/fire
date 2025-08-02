#' @title Initial Hyperparameter Estimation for I-prior Models
#' @name iprior_initial
#' @description Internal function to estimate initial hyperparameters for I-prior EM algorithm.
#'   Computes asymptotically optimal starting values for lambda and sigma parameters.
#'   Used internally by \code{\link{fire}}. Not exported.
#'
#' @param X Covariate matrix/array
#' @param Y Response vector
#' @param dat_T List of index values for each mode
#' @param kernels List of kernel functions. Refer to \code{\link{kernels_fire}} for available kernel functions
#' @param kernels_params List of parameters for each kernel
#' @param center Logical indicating whether to center the kernel matrices in \code{gmat}
#' @param std Logical indicating whether to standardise the kernel matrices in \code{gmat}
#' @param G List of kernel matrices from \code{gmat}
#' @param Index Matrix of indices for tensor vectorization
#' @param kernel_iprior Kernel function for iprior model
#' @param iprior_param  Parameter of kernel in iprior
#' @param constant_g Logical indicating whether to include constant kernel term in g
#' @param constant_h Logical indicating whether to include constant kernel term in h
#' @param os_type Operating system for parallelization ("Apple" or "Windows")
#' @param cores Number of cores for parallel computation (default: detectCores() - 1)
#' @param sample_id Which mode represents samples, either the 1st or the last mode
#' @param epsilon Small positive constant in initialisation of EM algorithm
#'
#' @return List containing two vectors of the initial values for EM algorithm:
#' \itemize{
#'   \item \code{lambda}: Asymptote for lambda (scale parameter)
#'   \item \code{sigma}: Asymptote for noise parameter
#' }
#'
#' @examples
#' # vector covariates
#' set.seed(1)
#' X <- matrix(runif(20, 5, 10), ncol = 4) # 5 samples
#' Y <- runif(5, 5, 10)
#' dat_T <- list(1:4)
#' Index <- expand.grid(1:5)
#'
#' init_cfbm <- fire:::pre_initial(X, Y, dat_T,
#'                         kernels = list(cfbm),
#'                         kernels_params = list(0.5),
#'                         Index = Index, constant_g = FALSE, constant_h = FALSE,
#'                         kernel_iprior = "cfbm", cores = 1)
#'
#' # 3D covariates
#' X <- list(array(runif(20, 5, 10), dim = c(2,2,5)),array(runif(20, 5, 10), dim = c(2,2,5)))
#' Y <- runif(2, 5, 10)
#' dat_T <- list(1:2, 1:2, 1:5)
#' Index <- expand.grid(T3 = 1:5, T2 = 1:2, T1 = 1:2)
#' Index <- cbind(Index$T1, Index$T2, Index$T3)
#'
#' init_rbf <- fire:::pre_initial(X, Y, dat_T,
#'                        kernels = list(cfbm, rbf, cfbm),
#'                        kernels_params = list(0.5, 1, 0.5),
#'                        Index = Index, constant_g = TRUE, constant_h = FALSE,
#'                        kernel_iprior = "rbf", cores = 1)
#'
#' @seealso \code{\link{fire.matrix}}, \code{\link{fire.tensor}}, \code{\link{kernels_fire}}
#'
#' @keywords internal
pre_initial <- function(X, Y, dat_T, kernels, kernels_params, center = FALSE, std = TRUE,
                        G = NULL,
                        Index, kernel_iprior = 'cfbm', iprior_param = NULL,
                        constant_g = TRUE, constant_h = FALSE, os_type = "Apple", cores = NULL, sample_id = 1,
                        epsilon = 1e-6) {

  m <- length(kernels)
  N <- length(Y)

  # Convert tensor to list if needed
  if (m != 1) {
    if (!is.list(X)) {
      X <- tensor_sample(X, sample_id = sample_id)
    }
  }

  # Compute Gram matrix and norm matrix
  if(is.null(G)){
    G <- gmat(kernels = kernels, kernels_params = kernels_params,
              dat = dat_T, center = center, std = std)
  }

  nmat <- Kronecker_norm_mat(X = X, G = G, alpha = rep(1, length(kernels)),
                             constant = constant_g, Index = Index,
                             os_type = os_type, cores = cores, sample_id = sample_id)

  # Generate Gram matrix based on I-prior kernel choice
  if (kernel_iprior == 'cfbm') {
    if (is.null(iprior_param)) {
      iprior_param <- 0.5
    }
    H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
  } else if (kernel_iprior == 'rbf') {
    if (is.null(iprior_param)) {
      iprior_param <- 1
    }
    H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
  } else if (kernel_iprior == 'linear') {
    H.tilde <- nmat
  } else {
    stop(paste("Unsupported kernel_iprior:", kernel_iprior,
               "- must be 'cfbm', 'rbf', or 'linear'"))
  }

  if(constant_h){
    H.tilde = 1 + H.tilde
  }

  # Eigendecomposition
  eigen.Htilde <- eigen(H.tilde, symmetric = TRUE)
  U <- eigen.Htilde$values
  V <- eigen.Htilde$vectors

  # Compute initial parameter estimates
  vt.y <- crossprod(V, Y)
  lambda.tilde <- sqrt((1/N) * sum((vt.y^2) / U^2))
  sigma <- sqrt((1/N) * sum(vt.y^2))
  if(is.null(epsilon)){
    epsilon = min(lambda.tilde,sigma)/10
  }
  lambda <- sqrt(sigma * epsilon)


  return(list(lambda = c(sqrt(epsilon * lambda.tilde), epsilon),
              sigma = c(lambda, sigma)))
}
