#' Efficient RKHS Norm Computation Using Kronecker Structure
#'
#' @description Compute RKHS norms and related matrices using Kronecker product structure.
#'
#' @param X,Y Input data (vector/matrix/list)
#' @param Xtrain Training data
#' @param Xnew New data
#' @param G List of kernel matrices from \code{\link{gmat}}
#' @param G_list Precomputed eigendecomposition from internal \code{\link{Kronecker_inv_helper}}
#' @param L_inv Inverse eigenvalues from internal \code{\link{Lambda_inv}}
#' @param alpha Vector of scale parameters
#' @param Index Matrix of indices for tensor vectorization
#' @param constant Logical indicating whether to include constant kernel term in g
#' @param os_type Parallelization backend ("Windows" or "Apple")
#' @param cores Number of cores used for parallel computation
#' @param sample_id Mode representing samples, either 1st or last mode
#' @param validate Logical indicating whether to carry out input validation
#' @return
#' \describe{
#'   \item{\code{rkhs_norm_kron}}{Single norm value between X and Y}
#'   \item{\code{Kronecker_norm_mat}}{Matrix of pairwise norms between training samples}
#'   \item{\code{Kronecker_norm_cross}}{Matrix of norms between test and training samples}
#' }
#'
#' @details These internal functions provide efficient computation of RKHS norms by:
#' \itemize{
#'   \item Leveraging Kronecker product structure
#'   \item Supporting parallel computation
#'   \item Handling both constant and non-constant kernel terms
#' }
#'
#' @examples
#' # Using internal fire package functions
#' G <- list(matrix(runif(9),3,3), matrix(runif(4),2,2))
#' decomp <- fire:::Kronecker_inv_helper(G, alpha = c(0.5, 0.5))
#' L_inv <- fire:::Lambda_inv(list(decomp$L1, decomp$L2))
#'
#' # Compute norm between two matrices
#' X <- matrix(rnorm(6), nrow=3)
#' rkhs_norm_kron(X, G_list = decomp, L_inv = L_inv)
#'
#' # Compute pairwise norm matrix
#' X_list <- list(X, matrix(rnorm(6), nrow=3))
#' Kronecker_norm_mat(X_list, G, alpha = c(0.5, 0.5))
#'
#' # Compute cross norms
#' X_new <- list(matrix(rnorm(6), nrow=3))
#' Kronecker_norm_cross(X_list, X_new, G, alpha = c(0.5, 0.5))
#'
#' @name rkhs_norm_kronecker
#' @keywords internal
NULL

#' @importFrom parallel makeCluster stopCluster detectCores parLapply mclapply clusterExport
#' @importFrom MASS ginv
NULL

#' @rdname rkhs_norm_kronecker
#' @export
rkhs_norm_kron <- function(X, Y = NULL, G_list, L_inv = NULL, constant = TRUE, Index = NULL, Ginv = NULL,
                           validate = FALSE) {
  # number of modes
  if(is.vector(X)){
    m = 1
  }else{
    m = length(G_list)/2
  }

  if(validate){
    # Input validation
    if (!is.numeric(X)) {
      stop("X must be numeric.")
    }
    if (!is.null(Y)) {
      if (!is.numeric(Y)) {
        stop("Y must be numeric.")
      }
      if (length(Y) != length(X) || !all(dim(Y) == dim(X))) {
        stop("Y must have the same size or length as X.")
      }
    }

    if (is.vector(X)) {
      if (length(G_list) != 2 && constant) {
        stop("For vector inputs, G_list must contain exactly 2 components (Q1 and L1), but got ",
             length(G_list), " components instead.")
      }else if(length(G_list) != 1 && !constant){
        stop("For vector inputs without constant kernel term, G_list must contain exactly 1 component (G1), but got ",
             length(G_list), " components instead.")
      }
    } else {
      expected_dims <- length(G_list) / 2
      if (length(dim(X)) != expected_dims) {
        stop("Dimension mismatch: Input tensor has ", length(dim(X)),
             " modes but G_list has ", expected_dims,
             " components because G_list must have 2 components for each mode.")
      }
    }
    if (is.vector(X)) {
      if(constant){
        if (length(X) != length(L_inv)) {
          stop("Dimension mismatch: Vector input length (", length(X),
               ") must match L_inv length (", length(L_inv), ").")
        }
      }
    } else {
      total_elements <- prod(dim(X))
      if (total_elements != length(L_inv)) {
        stop("Size mismatch: Total elements in matrix/array (", total_elements,
             ") must match L_inv length (", length(L_inv), ").")
      }
    }
  }

  # Handle special case for m=1 without constant term
  if (!constant) {
    if(is.null(Y)){
      Y = rep(0,length(X))
    }
    tnsr = X - Y

    if (m != 1) stop("constant = FALSE only supported for m=1")
    if (is.null(Ginv)) Ginv <- MASS::ginv(G_list$Q1)
    return(as.numeric(crossprod(tnsr, Ginv) %*% tnsr))
  }

  if (is.null(L_inv) && constant) {
    stop("L_inv must be provided when constant = TRUE")
  }

  if(constant){
    if(m == 1){
      if(is.null(Y)){
        Y = rep(0,length(X))
      }
      tnsr = X - Y
      return(sum((crossprod(tnsr, G_list$Q1))^2 / G_list$L1))
    }else if (m == 2) {
      dim_list = dim(X)

      if(is.null(Y)){
        Y = matrix(0, nrow = dim_list[1], ncol = dim_list[2])
      }
      tnsr = X - Y

      X2 <- as.vector(tnsr %*% G_list$Q2)
      X1 <- as.vector(crossprod(
        matrix(X2, nrow = ncol(G_list$Q1), byrow = FALSE),
        G_list$Q1
      ))
      return(sum(L_inv * X1^2))
    }else if (m == 3) {
      dim_list = dim(X)

      # vectorize tensors X and Y
      X_vec = vectorize_tensor(X, Index)
      if(is.null(Y)){
        Y_vec = rep(0, prod(dim_list))
      }else{
        Y_vec = vectorize_tensor(Y, Index)
      }
      tnsr = X_vec - Y_vec

      X3 <- kron_mv(tnsr, G_list$Q3)
      X2 <- kron_mv(X3, G_list$Q2)
      X1 <- kron_mv(X2, G_list$Q1)
      return(sum(L_inv * X1^2))
    }else{
    # Apply Kronecker products sequentially for arbitrary m > 3
      dim_list = dim(X)

      # vectorize tensors X and Y
      X_vec = vectorize_tensor(X, Index)
      if(is.null(Y)){
        Y_vec = rep(0, prod(dim_list))
      }else{
        Y_vec = vectorize_tensor(Y, Index)
      }
      tnsr = X_vec - Y_vec

    result <- tnsr
    for (i in m:1) {
      Q <- G_list[[paste0("Q", i)]]
      if (is.null(Q)) stop(paste0("Q", i, " not found in G_list"))
      result <- kron_mv(result, Q)
    }
    return(sum(L_inv * result^2))
    }
  }

}
#' @rdname rkhs_norm_kronecker
#' @export
Kronecker_norm_mat <- function(X, G, alpha, constant = TRUE, Index = NULL, G_list = NULL,
                               os_type = "Apple", cores = NULL,
                               sample_id = 1, validate = FALSE) {

  if(validate){
  if (is.matrix(X) && length(G) != 1) {
    stop('Dimension mismatch: X is one-way tensor but there are', length(G), 'kernel matrices in G')
  }
  }

  # Set default number of cores (safe for CRAN)
  if(is.null(cores)) {
    cores <- if(interactive()) {
      max(1, parallel::detectCores() - 1)
    } else {
      1  # Default to 1 core during checking/examples
    }
  }
  cores <- min(cores, parallel::detectCores())


  # number of modes
  m <- length(G)

  # Convert to list format if needed
  if(m != 1 && !is.list(X)) X <- tensor_sample(X, sample_id = sample_id)

  # Compute eigendecomposition if not provided
  if(is.null(G_list)) G_list <- Kronecker_inv_helper(G = G, alpha = alpha, constant = constant)

  # Compute inverse eigenvalues - generalized for any m
  if(m == 1 && constant) {
    # Special case for m=1 with constant term
    L_inv <- diag(1/G_list$L1)
  }else if(m == 1) {
    # Special case for m=1 without constant term
    L_inv <- NULL  # Will use Ginv instead
  } else {
    # General case for m > 1
    L_list <- lapply(1:m, function(i) G_list[[paste0("L", i)]])
    L_inv <- Lambda_inv(L_list)
  }

  N_sample <- if(m == 1) nrow(X) else length(X)
  mat <- matrix(0, N_sample, N_sample)
  Ginv <- if(m == 1 && !constant) MASS::ginv(G_list$G1) else NULL

  # Parallel computation setup
  compute_pairwise_norms <- function(i) {
    norms <- numeric(N_sample)
    for(j in i:N_sample) {
      if(m == 1) {
        norms[j] <- rkhs_norm_kron(X = X[i,], Y = X[j,], Index = Index, G_list = G_list, L_inv = L_inv, constant = constant, Ginv = Ginv)
      } else {
        norms[j] <- rkhs_norm_kron(X = X[[i]], Y = X[[j]], Index = Index , G_list = G_list, L_inv = L_inv, constant = constant, Ginv = NULL)
      }
    }
    return(norms)
  }

  if(os_type == "Windows") {
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    parallel::clusterExport(cl, c("rkhs_norm_kron", "vectorize_tensor", "kron_mv"),
                            envir = environment())
    pairwise_norms <- parallel::parLapply(cl, 1:N_sample, compute_pairwise_norms)
    parallel::stopCluster(cl)
  } else if(os_type == "Apple") {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("parallel package required for mclapply")
    }
    pairwise_norms <- parallel::mclapply(1:N_sample, compute_pairwise_norms,
                                         mc.cores = cores)
  } else {
    stop("Invalid os_type. Use 'Windows' or 'Apple'")
  }

  # Fill matrix
  for(i in 1:N_sample) {
    mat[i, i:N_sample] <- pairwise_norms[[i]][i:N_sample]
    mat[i:N_sample, i] <- mat[i, i:N_sample]
  }

  return(mat)
}

#' @rdname rkhs_norm_kronecker
#' @export
Kronecker_norm_cross <- function(Xtrain, Xnew, G, alpha, constant = TRUE, Index = NULL,
                                 G_list = NULL, os_type = "Apple", cores = NULL,
                                 sample_id = 1, validate = FALSE) {

  if(validate){
  if (is.matrix(Xtrain) && length(G) != 1) {
    stop('Dimension mismatch: Xtrain and Xnew are one-way tensor but there are', length(G), 'kernel matrices in G')
  }
  }

  # Set default number of cores (safe for CRAN)
  if(is.null(cores)) {
    cores <- if(interactive()) {
      max(1, parallel::detectCores() - 1)
    } else {
      1  # Default to 1 core during checking/examples
    }
  }
  cores <- min(cores, parallel::detectCores())  # Never exceed available cores

  m <- length(G)

  # Convert to list format if needed
  if(m != 1) {
    if(!is.list(Xtrain)) Xtrain <- tensor_sample(Xtrain, sample_id = sample_id)
    if(!is.list(Xnew)) Xnew <- tensor_sample(Xnew, sample_id = sample_id)
  }

  # Compute eigendecomposition if not provided
  if(is.null(G_list)) G_list <- Kronecker_inv_helper(G, alpha, constant)

  # Compute inverse eigenvalues - generalized for any m
  if(m == 1 && constant) {
    # Special case for m=1 with constant term
    L_inv <- diag(1/G_list$L1)
  }else if(m == 1) {
    # Special case for m=1 without constant term
    L_inv <- NULL  # Will use Ginv instead
  } else {
    # General case for m > 1
    L_list <- lapply(1:m, function(i) G_list[[paste0("L", i)]])
    L_inv <- Lambda_inv(L_list)
  }

  Ntrain <- if(m == 1) nrow(Xtrain) else length(Xtrain)
  Nnew <- if(m == 1) nrow(Xnew) else length(Xnew)
  mat <- matrix(0, Nnew, Ntrain)
  Ginv <- if(m == 1 && !constant) MASS::ginv(G_list$G1) else NULL

  # Parallel computation setup
  compute_norms_cross <- function(i) {
    norms <- numeric(Ntrain)
    for(j in 1:Ntrain) {
      if(m == 1) {
        norms[j] <- rkhs_norm_kron(X = Xnew[i,], Y = Xtrain[j,], Index = Index, G_list = G_list, L_inv = L_inv, constant = constant, Ginv = Ginv)
      } else {
        norms[j] <- rkhs_norm_kron(X = Xnew[[i]], Y = Xtrain[[j]], Index = Index, G_list = G_list, L_inv = L_inv, constant = constant, Ginv = NULL)
      }
    }
    return(norms)
  }

  if(os_type == "Windows") {
    cl <- parallel::makeCluster(cores, type = "PSOCK")
    parallel::clusterExport(cl, c("rkhs_norm_kron", "vectorize_tensor", "kron_mv"),
                            envir = environment())
    norms_list <- parallel::parLapply(cl, 1:Nnew, compute_norms_cross)
    parallel::stopCluster(cl)
  } else if(os_type == "Apple") {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("parallel package required for mclapply")
    }
    norms_list <- parallel::mclapply(1:Nnew, compute_norms_cross,
                                     mc.cores = cores)
  } else {
    stop("Invalid os_type. Use 'Windows' or 'Apple'")
  }

  # Fill matrix
  for(i in 1:Nnew) mat[i,] <- norms_list[[i]]

  return(mat)
}
