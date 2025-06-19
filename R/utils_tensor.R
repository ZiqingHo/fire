#' Tensor Utility Functions
#'
#' Collection of helper functions for tensor operations
#' @name utils_tensor
#' @aliases tensor_sample vectorize_tensor unfolding dat_unfolding
NULL

#' @rdname utils_tensor
#' @param X Input tensor
#' @param sample_id Sampling mode (1 or 4)
#' @return For \code{tensor_sample}: List of tensor samples
#' @export
#' @examples
#' # tensor_sample()
#' X <- array(1:24, dim = c(3,2,4))
#' tensor_sample(X, sample_id = 1)
tensor_sample <- function(X, sample_id = 1) {
  if (length(X) == 0) {
    stop("Input tensor X cannot be empty", call. = FALSE)
  }

  if (!sample_id %in% c(1, length(dim(X)))) {
    stop("Sample mode must be either 1st mode or last mode.")
  }

  num_mode <- length(dim(X)) - 1
  N <- dim(X)[sample_id]
  sample_list <- vector('list', N)
  names(sample_list) <- paste0('Sample-', 1:N)

  comma_str <- paste0(rep(',', num_mode), collapse = '')

  if (sample_id == 1) {
    for (i in 1:N) {
      index <- paste0(i, comma_str)
      sample_list[[i]] <- parse_index(X, index)
    }
  } else {
    for (i in 1:N) {
      index <- paste0(comma_str, i)
      sample_list[[i]] <- parse_index(X, index)
    }
  }

  sample_list
}

#' @rdname utils_tensor
#' @param X Input tensor
#' @param Index Matrix of indices
#' @return For `vectorize_tensor`: Vectorized tensor elements
#' @export
#' @examples
#' # vectorize_tensor()
#' X <- array(1:8, dim = c(2,2,2))
#' idx <- expand.grid(T1 = 1:2, T2 = 1:2, T3 = 1:2)
#' vectorize_tensor(X, idx)
#' idx <- cbind(rep(1:2,each = 4), rep(c(1,1,2,2),2), rep(1:2, 4))
#' vectorize_tensor(X, idx)
vectorize_tensor <- function(X, Index) {
  # Check for empty array
  if (length(X) == 0) {
    stop("Input tensor X cannot be empty", call. = FALSE)
  }

  # Convert Index to matrix if needed
  if (!is.matrix(Index)) {
    Index <- as.matrix(Index)
  }

  if(!is.numeric(Index)){
    stop("Index matrix must be integer")
  }

  # Validate index dimensions
  if (ncol(Index) != length(dim(X))) {
    stop(
      sprintf(
        "Index matrix must have %d columns (one for each tensor dimension), but has %d",
        length(dim(X)),
        ncol(Index)
      ),
      call. = FALSE
    )
  }

  # Check for invalid indices
  apply(Index, 1, function(idx) {
    invalid_dims <- which(idx < 1 | idx > dim(X) | is.na(idx))
    if (length(invalid_dims) > 0) {
      error_msg <- paste(
        "Invalid indices found:\n",
        paste(sprintf(
          "Dimension %d: index %d is out of bounds (must be between 1 and %d)",
          invalid_dims,
          idx[invalid_dims],
          dim(X)[invalid_dims]
        ), collapse = "\n"),
        sep = ""
      )
      stop(error_msg, call. = FALSE)
    }
    do.call(`[`, c(list(X), as.list(idx)))
  })
}

#' @rdname utils_tensor
#' @param Mat A matrix with N rows, each containing values to be assigned to the corresponding positions in the output tensor
#' @param dim Dimensions of output tensor
#' @param N Sample size
#' @param Index Matrix of indices
#' @return For `reshape_tensor`: A tensor with specified dimensions and the 1st mode corresponds to sampling mode
#' @export
#' @examples
#' # reshape_tensor()
#' N <- 2
#' Mat <- matrix(1:12, nrow = N)
#' dim <- c(2,3)
#' idx <- expand.grid(1:2, 1:3)
#' reshape_tensor(Mat, dim, N, idx)
reshape_tensor <- function(Mat, dim, N, Index) {
  if(all(is.na(Mat))){
    stop('Mat is empty matrix')
  }

  if (nrow(Index) != ncol(Mat)) {
    stop("Number of columns in Mat must match number of rows in Index")
  }
  if (length(dim) != ncol(Index)) {
    stop("Length of dim must match number of columns in Index")
  }

  if(dim(Mat)[1] == 1 && dim(Mat)[2] == 1){
    return(Mat)
  }
  # Validate Index against dim
  index_ranges <- lapply(seq_along(dim), function(i) {
    sort(unique(Index[,i]))
  })

  expected_ranges <- lapply(dim, function(d) 1:d)

  if (!identical(index_ranges, expected_ranges)) {
    # Find which dimensions have problems
    problems <- sapply(seq_along(dim), function(i) {
      !identical(index_ranges[[i]], expected_ranges[[i]])
    })

    error_msg <- paste(
      "Index matrix does not cover all combinations of dimensions:\n",
      paste(sprintf(
        "Dimension %d: indices should be 1:%d but found %s",
        which(problems),
        dim[problems],
        sapply(index_ranges[problems], function(x) paste(range(x), collapse = ":"))
      ), collapse = "\n"),
      "\nExpected all combinations from expand.grid(",
      paste0("1:", dim, collapse = ", "), ")",
      sep = ""
    )
    stop(error_msg, call. = FALSE)
  }

  # Create empty tensor with correct dimensions
  tensor_dims <- c(N, dim)
  tensor <- array(NA, dim = tensor_dims)

  # Create indexing list template
  idx_template <- vector("list", length(tensor_dims))
  idx_template[[1]] <- TRUE  # For the first dimension (N)

  # Fill the tensor
  for (n in 1:N) {
    idx_template[[1]] <- n  # Set the first index

    for (i in 1:nrow(Index)) {
      # Set the remaining indices from Index matrix
      idx_template[-1] <- as.list(Index[i, ])

      # Assign the value using do.call
      tensor <- do.call("[<-", c(list(tensor), idx_template, list(Mat[n, i])))
    }
  }

  return(tensor)
}

#' @rdname utils_tensor
#' @param n Unfolding mode
#' @param mode Logical to print status message
#' @return For `unfolding`: Unfolded matrix
#' @export
#' @examples
#' # unfolding()
#' X <- array(1:8, dim = c(2,2,2))
#' unfolding(X, 1)
unfolding <- function(X, n, mode = TRUE) {
  mode_dim <- dim(X)
  row_dim <- mode_dim[n]
  num_mode <- length(mode_dim)
  col_id <- (1:num_mode)[-n]
  col_dim <- prod(mode_dim[col_id])

  X_unfold <- matrix(aperm(X, perm = c(n, col_id)),
                     nrow = row_dim,
                     ncol = col_dim)

  if (mode) {
    message(paste('Mode-', n, ' matricization of tensor', sep = ''))
  }
  X_unfold
}

#' @rdname utils_tensor
#' @param X Input tensor
#' @param sample_id Sampling mode (1st or last mode)
#' @return For `dat_unfolding`: List of mode-n unfolded matrices
#' @export
#' @examples
#' # dat_unfolding()
#' X <- array(1:24, dim = c(3,2,4))
#' dat_unfolding(X, sample_id = 1)
dat_unfolding <- function(X, sample_id = 1, mode = TRUE) {
  mode_dim <- dim(X)
  sample_size <- mode_dim[sample_id]
  num_mode <- length(mode_dim)

  unfoldings <- vector('list', num_mode - 1)
  names(unfoldings) <- paste0('Mode-', (1:num_mode)[-sample_id])

  comma_str <- paste0(rep(',', num_mode - 1), collapse = '')

  if (sample_id == num_mode) {
    for (m in setdiff(1:num_mode, sample_id)) {
      unfoldings[[m]] <- lapply(1:sample_size, function(n) {
        index <- paste0(paste0(rep(',', sample_id - 1), collapse = ''), n)
        unfolding(parse_index(X, index), m, FALSE)
      })
    }
  } else {
    for (m in 1:num_mode) {
      if (m < sample_id) {
        unfoldings[[m]] <- lapply(1:sample_size, function(n) {
          index <- paste0(n, comma_str)
          unfolding(parse_index(X, index), m, FALSE)
        })
      } else if (m > sample_id) {
        unfoldings[[m - 1]] <- lapply(1:sample_size, function(n) {
          index <- paste0(n, comma_str)
          unfolding(parse_index(X, index), m - 1, FALSE)
        })
      }
    }
  }

  if (mode) {
    message(paste('The', suffix(sample_id), 'mode corresponds to the sample size.'))
  }
  unfoldings
}
