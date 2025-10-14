#' Fit a FIRE Model for Classification
#'
#' Estimates a functional classification model using I-priors with Reproducing Kernel Hilbert Space (RKHS) norm.
#' The FIRE model parameters are estimated using an EM algorithm.
#'
#' @param X A numeric input in \code{array} or \code{list}
#' @param Y A categorical response vector
#' @param dat_T List of index sets for each mode
#' @param kernels List of kernel functions for each mode (see \code{\link{kernels_fire}})
#' @param kernels_params List of parameters for each kernel
#' @param kernel_iprior Type of I-prior kernel
#' @param iprior_param Parameter for I-prior kernel:
#' \itemize{
#'   \item \code{"cfbm"} - Hurst coefficient (default 0.5)
#'   \item \code{"rbf"} - lengthscale (default 1)
#'   \item \code{"linear"} - offset (default 0)
#'   \item \code{"poly"} - degree and offset (default c(2, mean(Y)))
#' }
#' @param kernel_class Class kernel. Either \code{"identity"}, \code{"centred identity"},
#'   or a user-supplied positive semi-definite \eqn{c \times c} matrix.
#' @param class.labels Vector of distinct class labels (length \eqn{c}).
#' @param control A list of control parameters (see Details)
#'
#' @return An \code{fire_class} object. The \code{print()} and \code{summary()} methods display the model information.
#'
#' @details
#' The \code{control} argument can include the following parameters:
#' \itemize{
#'   \item{\code{scale}: Logical indicating whether to center the response (default TRUE)}
#'   \item{\code{maxiter}: Maximum number of EM iterations (default 200)}
#'   \item{\code{stop.eps}: Convergence tolerance (default 1e-3)}
#'   \item{\code{constant_g}: Logical indicating whether to include constant kernel term in g (default TRUE)}
#'   \item{\code{constant_h}: Logical indicating whether to include constant kernel term in h (default FALSE)}
#'   \item{\code{center}: Logical indicating whether to center the kernel g (default FALSE)}
#'   \item{\code{std}: Logical indicating whether to standardise the kernel g (default TRUE)}
#'   \item{\code{par_init}: Optional list of initial parameter values (lambda, noise)}
#'   \item{\code{os_type}: Operating system type for compatibility ("Apple" or "Windows", default "Apple")}
#'   \item{\code{cores}: Number of cores for parallel computation (default: detectCores() - 1)}
#'   \item{\code{asymptote}: Logical to use asymptotic initial values (default TRUE)}
#'   \item{\code{sample_id}: Which mode contains samples (default 1)}
#'   \item{\code{epsilon}: Small positive constant in initialisation of EM algorithm (default 1)}
#' }
#'
#' @examples
#' set.seed(42)
#' n_train <- 5; n_test  <- 2; n <- n_train + n_test
#' MatA <- matrix(c(1,1,1,0,0,0,0,0,0), nrow = 3, byrow = TRUE)
#' MatB <- matrix(c(0,0,0,0,0,0,1,1,1), nrow = 3, byrow = TRUE)
#' Y <- factor(c("A","A","B","B","A","A","B"), levels = c("A","B"))
#' # Generate matrices with small Gaussian noise
#' X <- lapply(seq_len(n), function(i) {base <- if (Y[i] == "A") MatA else MatB
#' noise <- matrix(rnorm(9, sd = 0.15), nrow = 3)
#' base + noise })
#'
#' # Train/test split: first 5 train, last 2 test
#' X_train <- X[1:n_train]; X_test  <- X[(n_train + 1):n]
#' Y_train <- Y[1:n_train]; Y_test  <- Y[(n_train + 1):n]
#'
#' dat_T <- list(1:3, 1:3)
#' mod <- fire_class(X = X_train, Y = Y_train, dat_T = dat_T,
#'  kernels = list(cfbm, cfbm), kernels_params = list(0.5, 0.5),
#'  class.labels = levels(Y_train), control = list(maxiter = 20, stop.eps = 1e-3))
#'
#' @seealso \code{\link{kernels_fire}}, \code{\link{utils_classification}}
#'
#' @export
fire_class <- function(X, Y, dat_T,
                        kernels, kernels_params,
                        kernel_iprior = 'cfbm', iprior_param = NULL,
                        kernel_class = 'centred identity', class.labels,
                        control) {

  # Set default control parameters
  con <- list(
    scale = TRUE,
    maxiter = 200,
    stop.eps = 1e-3,
    constant_g = TRUE,
    constant_h = FALSE,
    center = FALSE,
    std = TRUE,
    par_init = NULL,
    os_type = "Apple",
    cores = NULL,
    asymptote = TRUE,
    sample_id = 1,
    epsilon = 1
  )
  # Override defaults with user-supplied control parameters
  con[names(control)] <- control

  scale = con$scale
  maxiter = con$maxiter
  stop.eps = con$stop.eps
  constant_g = con$constant_g
  constant_h = con$constant_h
  center = con$center
  std = con$std
  par_init = con$par_init
  os_type = con$os_type
  cores = con$cores
  asymptote = con$asymptote
  sample_id = con$sample_id
  epsilon = con$epsilon

  input_type <- ifelse(is.list(X), 'list', 'array')

  if (is.list(X)){
    if(sample_id > length(X)){
      stop('sample_id must be the 1st or last mode')
    }
  }else{
    if(sample_id > length(dim(X))){
      stop('sample_id must be the 1st or last mode')
    }
  }
  if(length(kernels) != length(kernels_params)) {
    stop("kernels and kernels_params must be lists of equal length")
  }

  if(is.null(cores)) {
    cores <- if(interactive()) {
      max(1, parallel::detectCores() - 1)
    } else {
      1  # Default to 1 core during checking/examples
    }
  }
  cores <- min(cores, parallel::detectCores())

  # convert the tensor from array into list
  if (!is.list(X)) {
    if (is.array(X)) {
      # For array input
      X <- tensor_sample(X, sample_id = sample_id)
    } else {
      stop("X must be either a list of arrays or a single array")
    }
  }
  original_dims <- dim(X[[1]])

  # Validate list elements if X is list
  if (is.list(X) && !all(sapply(X, is.array))) {
    stop("All elements of X must be arrays/matrices")
  }

  labels <- Y
  num_class <- length(class.labels)
  # convert Y into dummy variables
  Y <- to_dummy_vector(Y, levels = class.labels)

  # sample size
  N <- length(Y)

  intercept <- if(scale) mean(Y) else 0
  Y <- Y - intercept
  Ymean <- mean(Y)

  if(is.null(iprior_param)){
    if (kernel_iprior == 'cfbm') {
      iprior_param <- 0.5
    }else if(kernel_iprior == 'rbf'){
      iprior_param <- 1
    }else if(kernel_iprior == 'linear'){
      iprior_param <- 0
    }else if(kernel_iprior == 'poly'){
      iprior_param <- c(2, Ymean)
    }
  }

  Index <- create_index_matrix(dat_T)

  # precompute the Gram matrix for class
  if(kernel_class == 'identity'){
    H.group <- diag(num_class)
  }else if(kernel_class == 'centred identity'){
    H.group <- diag(num_class) - matrix(1/num_class, num_class, num_class)
  }else if(is.matrix(kernel_class)){
    H.group <- kernel_class
  }

  # precompute kernel matrices G_m's
  G <- gmat(kernels = kernels, kernels_params = kernels_params, dat = dat_T, center = center, std = std)

  if(is.null(par_init)){
    if(asymptote == F){
      noise = 1
      tau =  1
      alpha_init = 1 # set alpha1 = alpha2 = alpha3
      init_points = list(list(lambda = tau, noise = noise))
    }else{
      pre_init_result = pre_initial(X = X, Y = Y, dat_T = dat_T,
                                    kernels = kernels, kernels_params = kernels_params,
                                    center = center, std = std, G = G,
                                    Index = Index, kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                                    constant_g = constant_g, constant_h = constant_h,
                                    os_type = os_type, cores = cores,
                                    sample_id = sample_id, epsilon = epsilon,
                                    is.classification = TRUE, num_class = num_class, kernel_class = kernel_class)
      init_points = list(
        list(lambda = pre_init_result$lambda[1], noise = pre_init_result$lambda[2]),
        list(lambda = pre_init_result$sigma[1], noise = pre_init_result$sigma[2])
      )
      alpha_init = 1
    }
  }else{
    init_points = list(
      list(lambda = par_init[2], noise = par_init[3])
    )
    alpha_init = par_init[1]
  }

  run_EM = function(init){
    tau = init$lambda
    noise = init$noise
    w = rep(0,N)
    W = diag(N)
    loglik = noise_est = tau_est = alpha_est = c()
    w_est = vector('list')

    # Initialization
    niter = 0
    diff.loglik = 1

    # start iteration
    start.time = Sys.time()
    pb = txtProgressBar(min = 0, max = maxiter, style = 3)

    while (niter < maxiter){
      if(niter == 0){
        alpha = rep(alpha_init, length(kernels))

        nmat = Kronecker_norm_mat(X = X, G = G,
                                  alpha = alpha, constant = constant_g,
                                  Index = Index, os_type = os_type, cores = cores, sample_id = sample_id)

        # Generate Gram matrix based on I-prior kernel choice
        if (kernel_iprior == 'cfbm') {
          H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
        } else if (kernel_iprior == 'rbf') {
          H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
        } else if (kernel_iprior == 'linear') {
          H.tilde <- nmat + iprior_param
        }else if (kernel_iprior == 'poly'){
          H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
        } else {
          stop(paste("Unsupported kernel_iprior:", kernel_iprior,
                     "- must be 'cfbm', 'rbf', 'linear' or 'poly"))
        }

        if(constant_h){
          H.tilde <- 1 + H.tilde
        }

        H <- tau^2 * kronecker(H.tilde, H.group)

        eigen.H =  eigen(H, symmetric = TRUE)

        U = eigen.H$values
        V = eigen.H$vectors
        d = (U^2)/(noise^2) + noise^2 #eigenvalues of V_y
        Vy.inv = tcrossprod(V%*%diag(1/d),V)
        w = (H%*%Vy.inv%*%Y)/(noise^2)
        W = Vy.inv + tcrossprod(w)
      }

      ## optimization wrt alpha (alpha 1 = alpha2 = alpha3)
      res = optim(par = alpha_init, fn = Qfun_class,
                  X = X, Y = Y , dat_T = dat_T,
                  G = G, kernels = kernels, kernels_params = kernels_params,
                  kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                  kernel_class = kernel_class, num_class = num_class,
                  tau = tau, noise = noise,
                  Index = Index, W = W, w = w, constant_g = constant_g, constant_h = constant_h,
                  os_type = os_type, cores = cores, sample_id = sample_id,
                  method = 'L-BFGS-B', lower = 1e-3, upper = 1e4,
                  control = list(maxit = 2))

      alpha_init = alpha_est[niter + 1] = res$par

      alpha = rep(alpha_init, length(kernels))
      nmat = Kronecker_norm_mat(X = X, G = G,
                                alpha = alpha, constant = constant_g,
                                Index = Index, os_type = os_type, cores = cores, sample_id = sample_id)

      # Generate Gram matrix based on I-prior kernel choice
      if (kernel_iprior == 'cfbm') {
        H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
      } else if (kernel_iprior == 'rbf') {
        H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
      } else if (kernel_iprior == 'linear') {
        H.tilde <- nmat + iprior_param
      } else if(kernel_iprior == 'poly'){
        H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
      }

      if(constant_h){
        H.tilde <- 1 + H.tilde
      }
    #  H.tilde <- kronecker(H.tilde, H.group)
      # update tau
      res = optim(par = tau, fn = Qfun_class,
                  X = X, Y = Y , dat_T = dat_T,
                  G = G, kernels = kernels, kernels_params = kernels_params,
                  kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                  kernel_class = kernel_class, num_class = num_class,
                  noise = noise, alpha = alpha_init, H.tilde = H.tilde,
                  Index = Index, W = W, w = w, constant_g = constant_g, constant_h = constant_h,
                  os_type = os_type, cores = cores, sample_id = sample_id,
                  method = 'L-BFGS-B', lower = 1e-3, upper = 1e4,
                  control = list(maxit = 2))
      tau = tau_est[niter+1] = res$par
      H =  tau^2 * kronecker(H.tilde, H.group)

      eigen.H =  eigen(H, symmetric = TRUE)

      U = eigen.H$values
      V = eigen.H$vectors
      # update noise
      Hsq = V%*%(t(V) * U^2) #tcrossprod(V%*%diag(U^2),V)
      Nnoise = crossprod(Y) - 2 * crossprod(Y,H%*% w) + sum(Hsq * W)
      Dnoise = sum(diag(W))
      noise = noise_est[niter +1 ] = as.numeric(sqrt(sqrt(Nnoise / Dnoise)))

      # update w, W
      noise_sq = noise^2
      d = (U^2)/(noise_sq) + noise_sq #eigenvalues of V_y
      Vy.inv = tcrossprod(V%*%diag(1/d),V)
      w = w_est[[niter + 1]] = (H%*%Vy.inv%*%Y)/(noise_sq)
      W = Vy.inv + tcrossprod(w)

      loglik[niter + 1] = - N/2*log(2*pi) - 0.5 * crossprod(Y,Vy.inv)%*%Y - 0.5*sum(log(d))

      if(niter == 0){
        diff.loglik = loglik[niter+1]
      }else{
        diff.loglik = loglik[niter+1] - loglik[niter]
        if(loglik[niter+1] < loglik[niter]){
          break
        }
      }

      niter = niter + 1

      setTxtProgressBar(pb, niter)

      if (niter == maxiter | abs(diff.loglik) < stop.eps) {
        break
      }

    }

    close(pb)

    end.time = Sys.time()
    duration = end.time - start.time

    return(list(alpha = alpha_est, tau = tau_est, noise = noise_est,
                mloglik = loglik,
                niter = niter, duration = duration,
                w = w_est,
                init_tau = init[[1]], init_noise = init[[2]],
                diff.loglik = diff.loglik))

  }

  # Run EM for all initial points
  # results = lapply(init_points, run_EM)

  # Initialize results list
  results <- vector("list", length = length(init_points))

  # Run EM for initial points sequentially, stopping if convergence is achieved
  for (i in seq_along(init_points)) {
    results[[i]] <- run_EM(init_points[[i]])

    # Check convergence (just for reporting, not for stopping)
    current_converged <- results[[i]]$niter < maxiter && abs(results[[i]]$diff.loglik) < stop.eps

    if (current_converged) {
      if(is.null(par_init)) {
        if (i == 2) {
          message("Convergence achieved with the 2nd set of initial points.")
        } else {
          message("Convergence achieved with the 1st set of initial points.")
        }
      } else {
        message("Convergence achieved with the initial points that you set.")
      }
    }
  }

  # Clean up results in case we broke out early
  results <- results[!sapply(results, is.null)]

  # Select the best result based on convergence status and log-likelihood
  if (length(results) >= 2) {
    # Check which runs converged
    conv_status <- sapply(results, function(res) res$niter < maxiter && abs(res$diff.loglik) < stop.eps)

    if (any(conv_status)) {
      # Among converged runs, pick the one with the highest log-likelihood
      max_loglik_index <- which.max(sapply(results[conv_status], function(res) tail(res$mloglik, 1)))
      max_loglik_index <- which(conv_status)[max_loglik_index]  # Map back to original index
    } else {
      # If none converged, pick the best log-likelihood among all
      max_loglik_index <- which.max(sapply(results, function(res) tail(res$mloglik, 1)))
    }
  } else {
    # Only one result (either only one init point or all others NULL)
    max_loglik_index <- 1
  }

  best_result = results[[max_loglik_index]]

  # Add chosen initial values to the final result
  best_result$init_tau = init_points[[max_loglik_index]]$lambda
  best_result$init_noise = init_points[[max_loglik_index]]$noise

  # Add convergence status
  best_result$converged <- best_result$niter < maxiter && abs(best_result$diff.loglik) < stop.eps

  # Print convergence message
  if(best_result$converged) {
    message(sprintf(
      "Converged in %d iterations (mlogLik: %.3f, tol: %.1e)",
      best_result$niter,
      tail(best_result$mloglik, 1),
      stop.eps
    ))
  } else {
    if(best_result$niter < maxiter) {
      warning(sprintf(
        "Stopped early at iteration %d/%d despite non-decreasing mlogLik.
        Possible convergence to saddle point or very flat region.
        Final mlogLik: %.3f (rel. change: %.1e)",
        best_result$niter,
        maxiter,
        tail(best_result$mloglik, 1),
        abs(best_result$diff.loglik)
      ))
    } else {
      warning(sprintf(
        "Did not converge after %d iterations (final rel. change: %.1e)",
        maxiter,
        abs(best_result$diff.loglik)
      ))
    }
  }
  # Return structured result
  this_call <- match.call()

  output <- structure(
    best_result,
    class = c("fire_class"),
    input_type = input_type,
    dimensions = original_dims,
    call = this_call,
    timestamp = Sys.time(),
    training_data = X,
    original_response = Y,
    labels = labels,
    kernels = kernels,
    kernels_params = kernels_params,
    kernel_iprior = kernel_iprior,
    kernel_class = kernel_class,
    num_class = num_class,
    class.labels = class.labels,
    iprior_param = iprior_param,
    constant_g = constant_g,
    constant_h = constant_h,
    dat_T = dat_T,
    center = center,
    std = std,
    os_type = os_type,
    cores = cores,
    intercept = intercept,
    sample_id = sample_id,
    sample_size = N,
    convergence = list(
      converged = best_result$converged,
      niter = best_result$niter,
      maxiter = maxiter,
      tolerance = stop.eps,
      final_change = abs(best_result$diff.loglik)
    )
  )

  return(output)

}

# Keep Qfun_class as an unexported internal function in the same file
Qfun_class <- function(X, Y, dat_T,
                        G, kernels, kernels_params,
                        kernel_iprior, iprior_param,
                        kernel_class = 'centred identity', num_class,
                        tau, noise, alpha,
                        H.tilde = NULL, Index, W, w, constant_g = TRUE, constant_h = FALSE,
                        os_type = "Apple", cores = 1,
                        sample_id = 1){

  # set the overall scale parameter alpha_0 = 1 for identification
  alpha = rep(alpha, length(kernels))

  # use cfbm to construct gram matrix H without scale parameter tau
  if(is.null(H.tilde)){
    nmat = Kronecker_norm_mat(X = X, G = G,
                              alpha = alpha, constant = constant_g,
                              Index = Index, os_type = os_type, cores = cores, sample_id = sample_id)
    if (kernel_iprior == 'cfbm') {
      H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
    } else if (kernel_iprior == 'rbf') {
      H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
    } else if (kernel_iprior == 'linear') {
      H.tilde <- nmat
    } else if(kernel_iprior == 'poly'){
      H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
    } else {
      stop(paste("Unsupported kernel_iprior:", kernel_iprior,
                 "- must be 'cfbm', 'rbf', 'linear' or 'poly"))
    }}

  if(constant_h){
    H.tilde <- 1 + H.tilde
  }

  if(kernel_class == 'identity'){
    H.group <- diag(num_class)
  }else if(kernel_class == 'centred identity'){
    H.group <- diag(num_class) - matrix(1/num_class, num_class, num_class)
  }else if(is.matrix(kernel_class)){
    H.group <- kernel_class
  }

  H = tau^2 * kronecker(H.tilde, H.group)
  H.eigen = eigen(H)
  U = H.eigen$values
  V = H.eigen$vectors
  d = (U^2)/(noise^2) + noise^2 #eigenvalues of V_y
  Vy = tcrossprod(V%*%diag(d),V)
  Q = - 0.5 * sum(Y^2)/(noise^2) + crossprod(Y,H%*%w)/(noise^2) - 0.5*sum(Vy*W)

  return(-Q) # negative for maximization
}

