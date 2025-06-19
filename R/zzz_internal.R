#' @title Internal Utility Functions
#' @name internal_utils
#' @description Internal helper functions for string and index processing.
#'   Not intended for direct use by package users.
#' @keywords internal
NULL

#' @rdname internal_utils
#'
#' @param id An integer or numeric value
#' @return \code{suffix}: Character string with ordinal suffix appended (e.g. "1st", "2nd")
#'
#' @examples
#' \dontrun{
#' suffix(1)   # "1st"
#' suffix(22)  # "22nd"
#' suffix(13)  # "13th"
#' }
suffix <- function(id) {
  if(!is.numeric(id)){
    stop('Error: id must be numeric.')
  }
  if (id %% 100 %in% 11:13) return(paste0(id, 'th'))
  r <- switch(id %% 10 + 1, 'th', 'st', 'nd', 'rd','th','th','th','th','th','th','th')
  paste0(id, r)
}

#' @rdname internal_utils
#'
#' @param X A tensor or array to subset
#' @param string Character string of comma-separated indices (e.g. "1,2,")
#' @return \code{parse_index}: Subsetted tensor/array
#'
#' @examples
#' \dontrun{
#' X <- array(1:8, dim = c(2,2,2))
#' parse_index(X, "1,2,")  # Returns X[1,2,]
#' parse_index(X, ",1")    # Returns X[,1]
#' }
parse_index <- function(X, string) {
  str_comp <- strsplit(string, ",", fixed = TRUE)[[1]]

  index <- lapply(str_comp, function(x) {
    if (x == "") {
      TRUE
    } else {
      if (!grepl("^[0-9.]+$", x)) {
        stop("Invalid index: '", x, "'. Must be a numeric value or empty.", call. = FALSE)
      }
      as.numeric(x)
    }
  })

  if (substr(string, nchar(string), nchar(string)) == ",") {
    index <- c(index, TRUE)
  }
  do.call('[', c(list(X), index))
}

#' @rdname internal_utils
#'
#' @param dat_T A list of numeric vectors, where each vector represents the time points
#'        for a different variable.
#'
#' @return \code{create_index_matrix}: A numeric matrix where each row represents a unique combination of time points.
#'
#' @examples
#' \dontrun{
#' # Two variables
#' create_index_matrix(list(c(1, 2), c(3, 4, 5)))
#'
#' # Three variables
#' create_index_matrix(list(c(1, 2), c(3, 4), seq(0,1, length.out = 4)))
#' }
create_index_matrix <- function(dat_T) {
  m <- length(dat_T)

  # Create sequences of 1-based indices for each variable
  index_sequences <- lapply(dat_T, function(x) seq_along(x))

  # Create names in reverse order (T1 varies fastest)
  grid_args <- lapply(rev(seq_along(index_sequences)), function(i) index_sequences[[i]])
  names(grid_args) <- paste0("T", rev(seq_along(index_sequences)))

  # Generate all combinations
  index_grid <- do.call(expand.grid, grid_args)

  # Reorder columns to T1, T2, ..., Tm sequence and convert to integer matrix
  col_order <- rev(seq_along(index_sequences))
  index_matrix <- as.matrix(index_grid[, col_order, drop = FALSE])
  storage.mode(index_matrix) <- "integer"

  # Rename columns to match original format
  colnames(index_matrix) <- paste0("T", seq_along(dat_T))

  return(index_matrix)
}
