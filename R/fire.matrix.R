#' @title FIRE Model for Matrix Input
#' @rdname fire
#' @export
fire.matrix <- function(X, Y, dat_T,
                        kernels = list(cfbm), kernels_params = list(0.5),
                        kernel_iprior = 'cfbm', iprior_param = NULL,
                        control = list(),...) {

  # Set default control parameters
  con <- list(
    scale = TRUE,
    maxiter = 200,
    stop.eps = 1e-5,
    center = FALSE,
    constant_h = FALSE,
    par_init = NULL,
    os_type = "Apple",
    cores = NULL,
    asymptote = TRUE
  )
  # Override defaults with user-supplied control parameters
  con[names(control)] <- control

  scale = con$scale
  maxiter = con$maxiter
  stop.eps = con$stop.eps
  constant_h = con$constant_h
  center = con$center
  par_init = con$par_init
  os_type = con$os_type
  cores = con$cores
  asymptote = con$asymptote
  sample_id = con$sample_id

  # Input validation
  if (!is.matrix(X)) stop("X must be a matrix")
  if (nrow(X) != length(Y)) {
    stop("Number of rows in X (", nrow(X),
         ") must match length of Y (", length(Y), ")")
  }
  if (anyNA(X)) stop("NA values detected in X")
  if (anyNA(Y)) stop("NA values detected in Y")
  if (!kernel_iprior %in% c("cfbm", "rbf", "linear")) {
    stop("kernel_iprior must be one of: 'cfbm', 'rbf', 'linear'")
  }
  if(is.null(cores)) {
    cores <- if(interactive()) {
      max(1, parallel::detectCores() - 1)
    } else {
      1  # Default to 1 core during checking/examples
    }
  }
  cores <- min(cores, parallel::detectCores())
  # only work for without constant kernel term
  constant_g = FALSE

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

  N <- length(Y)
  d <- ncol(X)
  Index <- matrix(1:d, ncol = 1)

  # If asymptote is TRUE, use pre_initial to get initial points
  if (asymptote) {
    pre_init_result = pre_initial(X = X, Y = Y, dat_T = dat_T,
                                  kernels = kernels, kernels_params = kernels_params,
                                  Index = Index, kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                                  constant_g = constant_g, constant_h = constant_h,
                                  os_type = os_type, cores = cores,
                                  sample_id = 1)
    init_points = list(
      list(lambda = pre_init_result$lambda[1], noise = pre_init_result$lambda[2]),
      list(lambda = pre_init_result$sigma[1], noise = pre_init_result$sigma[2])
    )
  } else {
    # Use par_init if provided, otherwise default values
    if (is.null(par_init)) {
      init_points = list(list(lambda = 1, noise = 1))
    } else {
      # Check if par_init is a list of vectors
      if (is.list(par_init)) {
        # If par_init is a list of vectors, convert it to a list of lists
        init_points = lapply(par_init, function(x) list(lambda = x[1], noise = x[2]))
      } else {
        # If par_init is a vector, use it directly
        init_points = list(list(lambda = par_init[1], noise = par_init[2]))
      }
    }
  }
  # Function to run EM for a single initial point
  run_EM = function(init) {
    lambda = init$lambda
    noise = init$noise
    w = rep(0, N)
    W = diag(N)
    loglik = noise_est = lambda_est = c()
    w_est = list()

    # Precompute G and H.tilde
    G = gmat(kernels = kernels, kernels_params = kernels_params, dat = dat_T, center = center)
    nmat = Kronecker_norm_mat(X = X, G = G,
                              alpha = c(1), constant = constant_g,
                              Index = Index,  os_type = os_type, cores = cores, sample_id = 1)
    # Generate Gram matrix based on I-prior kernel choice
    if (kernel_iprior == 'cfbm') {
      H.tilde <- cfbm_rkhs_kron(nmat = nmat, Hurst = iprior_param)
    } else if (kernel_iprior == 'rbf') {
      H.tilde <- rbf_rkhs_kron(nmat = nmat, lengthscale = iprior_param)
    } else if (kernel_iprior == 'linear') {
      H.tilde <- nmat + iprior_param
    }else if (kernel_iprior == 'poly'){
      H.tilde <-  (nmat + iprior_param[2])^iprior_param[1]
    }

    if(constant_h){
      H.tilde <- 1 + H.tilde
    }

    # Start iteration
    niter = 0
    diff.loglik = 1
    start.time = Sys.time()

    pb <- txtProgressBar(min = 0, max = maxiter, style = 3)


    while (niter < maxiter) {
      ## Optimization wrt lambda
      res = optim(par = lambda, fn = Qfun_matrix,
                  X = X, Y = Y,
                  dat_T = dat_T, G = G,
                  kernels = kernels,kernels_params = kernels_params,
                  kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                  noise = noise,
                  H.tilde = H.tilde,
                  Index = Index,
                  W = W, w = w, constant_g = constant_g, constant_h = constant_h,
                  os_type = os_type, cores = cores,
                  method = 'L-BFGS-B',
                  control = list(maxit = 2))

      lambda_new = lambda_est[niter + 1] = res$par

      # Update lambda
      H = lambda_new^2 * H.tilde
      eigen.H = eigen(H, symmetric = TRUE)
      U = eigen.H$values
      V = eigen.H$vectors

      # Update noise
      Hsq = V %*% (t(V) * U^2)
      Nnoise = crossprod(Y) - 2 * crossprod(Y, H %*% w) + sum(Hsq * W)
      Dnoise = sum(diag(W))
      noise = noise_est[niter + 1] = as.numeric(sqrt(sqrt(Nnoise / Dnoise)))

      # Update w and W
      d = (U^2) / (noise^2) + noise^2
      Vy.inv = tcrossprod(V %*% diag(1 / d), V)
      w = w_est[[niter + 1]] = (H %*% Vy.inv %*% Y) / (noise^2)
      W = Vy.inv + tcrossprod(w)

      # Update log-likelihood
      loglik[niter + 1] = -N / 2 * log(2 * pi) - 0.5 * crossprod(Y, Vy.inv) %*% Y - 0.5 * sum(log(d))

      # Check convergence
      if (niter == 0) {
        diff.loglik = loglik[niter + 1]
      } else {
        diff.loglik = loglik[niter + 1] - loglik[niter]
        if (loglik[niter + 1] < loglik[niter] - 0.1) break # numerical optimization may have rounding error
      }

      niter = niter + 1
      setTxtProgressBar(pb, niter)
      if (niter == maxiter | abs(diff.loglik) < stop.eps) break
    }

    close(pb)
    end.time = Sys.time()
    duration = end.time - start.time

    return(list(lambda = lambda_est, noise = noise_est,
                mloglik = loglik,
                niter = niter, duration = duration,
                w = w_est,
                init_lambda = init$lambda, init_noise = init$noise,
                diff.loglik = diff.loglik))
  }

  # Run EM for all initial points
  results = lapply(init_points, run_EM)

  # Select the result with the maximum log-likelihood
  max_loglik_index = which.max(sapply(results, function(res) max(res$mloglik, na.rm = TRUE)))
  best_result = results[[max_loglik_index]]

  # Add chosen initial values to the final result
  best_result$init_lambda = init_points[[max_loglik_index]]$lambda
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
  output <- structure(
    best_result,
    class = c("fire_matrix"),
    input_type = "matrix",
    dimensions = dim(X),
    call = match.call(),
    timestamp = Sys.time(),
    training_data = X,
    original_response = Y,
    kernels = kernels,
    kernels_params = kernels_params,
    kernel_iprior = kernel_iprior,
    iprior_param = iprior_param,
    constant_g = constant_g,
    constant_h = constant_h,
    dat_T = dat_T,
    center = center,
    os_type = os_type,
    cores = cores,
    intercept = intercept,
    convergence = list(
      converged = best_result$converged,
      niter = best_result$niter,
      maxiter = maxiter,
      tolerance = stop.eps,
      final_change = abs(best_result$diff.loglik)
    )
  )
  output

}

# Keep Qfun_matrix as an unexported internal function in the same file
Qfun_matrix <- function(X, Y, dat_T, G, kernels, kernels_params,
                        kernel_iprior, iprior_param, lambda, noise,
                        H.tilde = NULL, Index, W, w, constant_g = TRUE, constant_h = FALSE,
                        os_type = "Apple", cores = NULL) {

  # sample size
  N = length(Y)

  # use cfbm to construct gram matrix H without scale parameter tau
  if(is.null(H.tilde)){
    nmat = Kronecker_norm_mat(X = X, G = G,
                              alpha = c(1), constant = constant_g,
                              Index = Index,  os_type = os_type, cores = cores, sample_id = 1)
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
      if (is.null(iprior_param)) {
        iprior_param <- 0
      }
      H.tilde <- nmat + iprior_param
    }else if (kernel_iprior == 'poly'){
      if (is.null(iprior_param)) {
        iprior_param <- c(2, mean(Y))
      }
      H.tilde <- (nmat + iprior_param[2])^iprior_param[1]
    } else {
      stop(paste("Unsupported kernel_iprior:", kernel_iprior,
                 "- must be 'cfbm', 'rbf', 'linear' or 'poly"))
    }
    if(constant_h){
      H.tilde <- 1 + H.tilde
    }
    }



  H = lambda^2 * H.tilde
  H.eigen = eigen(H)
  U = H.eigen$values
  V = H.eigen$vectors
  d = (U^2)/(noise^2) + noise^2 #eigenvalues of V_y
  Vy = tcrossprod(V%*%diag(d),V)
  Q = - 0.5 * sum(Y^2)/(noise^2) + crossprod(Y,H%*%w)/(noise^2) - 0.5*sum(Vy*W)

  return(-Q) # negative for maximization
}
