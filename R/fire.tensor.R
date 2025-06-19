#' @title FIRe Model for Tensor Input
#' @param Y A numeric response vector
#' @param scale Logical indicating whether to center the response by subtracting mean(Y)
#' @param dat_T List of index sets for each mode
#' @param kernels List of kernel functions for each mode, may refer \code{\link{kernels_fire}} for details
#' @param kernels_params List of parameters for each kernel
#' @param kernel_iprior Type of I-prior kernel
#' @param iprior_param Parameter for I-prior kernel:
#' \itemize{
#' \item{\code{"cfbm"}} - Hurst
#' \item{\code{"rbf"}} - lengthscale
#' \item{\code{"linear"}} - offset
#' \item{\code{"poly"}} - degree and offset
#' }
#' @param maxiter Maximum number of EM iterations
#' @param stop.eps Convergence tolerance
#' @param constant Logical indicating whether to include the constant kernel term
#' @param center Logical indicating whether to center the kernel matrix
#' @param par_init Optional list of initial parameter values (lambda, noise)
#' @param os_type Operating system type for compatibility ("Apple" or "Windows")
#' @param asymptote Logical to use asymptotic initial values
#' @rdname fire
#' @export
fire.tensor <- function(X, Y, ...,
                        dat_T, scale = TRUE,
                        kernels, kernels_params,
                        kernel_iprior = 'cfbm', iprior_param = NULL,
                        maxiter = 200, stop.eps = 1e-5, center = FALSE,
                        par_init = NULL, constant = TRUE, os_type = "Apple", asymptote = TRUE, sample_id = 1
){

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

  if(is.null(iprior_param)){
    if (kernel_iprior == 'cfbm') {
      iprior_param <- 0.5
    }else if(kernel_iprior == 'rbf'){
      iprior_param <- 1
    }else if(kernel_iprior == 'linear'){
      iprior_param <- 0
    }else if(kernel_iprior == 'poly'){
      iprior_param <- c(2, 0)
    }
  }

  # sample size
  N = length(Y)

  intercept <- if(scale) mean(Y) else 0
  Y <- Y - intercept

  Index <- create_index_matrix(dat_T)

  # precompute kernel matrices G_m's
  G = gmat(kernels = kernels, kernels_params = kernels_params, dat = dat_T, center = center)

  if(is.null(par_init)){
    if(asymptote == F){
      noise = 1
      tau =  1
      alpha_init = 1 # set alpha1 = alpha2 = alpha3
      init_points = list(list(lambda = tau, noise = noise))
    }else{
      pre_init_result = pre_initial(X = X, Y = Y, dat_T = dat_T,
                                    kernels = kernels, kernels_params = kernels_params,
                                    G = G,
                                    Index = Index, kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                                    constant = constant,
                                    os_type = os_type,
                                    sample_id = sample_id)
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
                                  alpha = alpha, constant = constant,
                                  Index = Index, os_type = os_type, sample_id = sample_id)

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

        H =  tau^2 * H.tilde

        eigen.H =  eigen(H, symmetric = TRUE)

        U = eigen.H$values
        V = eigen.H$vectors
        d = (U^2)/(noise^2) + noise^2 #eigenvalues of V_y
        Vy.inv = tcrossprod(V%*%diag(1/d),V)
        w = (H%*%Vy.inv%*%Y)/(noise^2)
        W = Vy.inv + tcrossprod(w)
      }

      ## optimization wrt alpha (alpha 1 = alpha2 = alpha3)
      res = optim(par = alpha_init, fn = Qfun_tensor,
                  X = X, Y = Y , dat_T = dat_T,
                  G = G, kernels = kernels, kernels_params = kernels_params,
                  kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                  tau = tau, noise = noise,
                  Index = Index, W = W, w = w, constant = constant,
                  os_type = os_type, sample_id = sample_id,
                  method = 'L-BFGS-B', lower = 1e-3, upper = 1e4,
                  control = list(maxit = 2))

      alpha_init = alpha_est[niter + 1] = res$par

      alpha = rep(alpha_init, length(kernels))
      nmat = Kronecker_norm_mat(X = X, G = G,
                                alpha = alpha, constant = constant,
                                Index = Index, os_type = os_type, sample_id = sample_id)

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

      # update tau
      res = optim(par = tau, fn = Qfun_tensor,
                  X = X, Y = Y , dat_T = dat_T,
                  G = G, kernels = kernels, kernels_params = kernels_params,
                  kernel_iprior = kernel_iprior, iprior_param = iprior_param,
                  noise = noise, alpha = alpha_init, H.tilde = H.tilde,
                  Index = Index, W = W, w = w, constant = constant,
                  os_type = os_type, sample_id = sample_id,
                  method = 'L-BFGS-B', lower = 1e-3, upper = 1e4,
                  control = list(maxit = 2))
      tau = tau_est[niter+1] = res$par
      H =  tau^2 * H.tilde

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

  # Initialize results list and convergence flag
  results <- vector("list", length = length(init_points))
  converged_flag <- FALSE

  # Run EM for initial points sequentially, stopping if convergence is achieved
  for (i in seq_along(init_points)) {
    if (!converged_flag) {
      results[[i]] <- run_EM(init_points[[i]])

      # Check if this run converged (using the same convergence criteria as later in the code)
      current_converged <- results[[i]]$niter < maxiter && abs(results[[i]]$diff.loglik) < stop.eps

      if (current_converged) {
        converged_flag <- TRUE
        message("Convergence achieved with initial point ", i, ". Skipping remaining initial points.")
        break
      }
    }
  }

  # Clean up results in case we broke out early
  results <- results[!sapply(results, is.null)]

  # Select the result with the maximum log-likelihood
  max_loglik_index = which.max(sapply(results, function(res) max(res$mloglik, na.rm = TRUE)))
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
    warning(sprintf(
      "Did not converge after %d iterations (final rel. change: %.1e)",
      maxiter,
      abs(best_result$diff.loglik)
    ))
  }

  # Return structured result
  structure(
    best_result,
    class = c("fire_tensor"),
    input_type = input_type,
    dimensions = original_dims,
    call = match.call(),
    timestamp = Sys.time(),
    training_data = X,
    original_response = Y,
    kernels = kernels,
    kernels_params = kernels_params,
    kernel_iprior = kernel_iprior,
    iprior_param = iprior_param,
    constant = constant,
    dat_T = dat_T,
    center = center,
    os_type = os_type,
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
}

# Keep Qfun_tensor as an unexported internal function in the same file
Qfun_tensor <- function(X, Y, dat_T,
                        G, kernels, kernels_params,
                        kernel_iprior, iprior_param,
                        tau, noise, alpha,
                        H.tilde = NULL, Index, W, w, constant = TRUE,
                        os_type = "Apple",
                        sample_id = 1){

  # sample size
  N = length(Y)
  # set the overall scale parameter alpha_0 = 1 for identification
  alpha = rep(alpha, length(kernels))

  # use cfbm to construct gram matrix H without scale parameter tau
  if(is.null(H.tilde)){
    nmat = Kronecker_norm_mat(X = X, G = G,
                              alpha = alpha, constant = constant,
                              Index = Index, os_type = os_type, sample_id = sample_id)
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

  H = tau^2 * H.tilde
  H.eigen = eigen(H)
  U = H.eigen$values
  V = H.eigen$vectors
  d = (U^2)/(noise^2) + noise^2 #eigenvalues of V_y
  Vy = tcrossprod(V%*%diag(d),V)
  Q = - 0.5 * sum(Y^2)/(noise^2) + crossprod(Y,H%*%w)/(noise^2) - 0.5*sum(Vy*W)

  return(-Q) # negative for maximization
}

