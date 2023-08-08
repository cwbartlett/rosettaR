# Check if correlation or covariance matrices are positive semi definite
is_pos_semidef <- function(matrix) {
  all(eigen(matrix, only.values = TRUE)$values >= 0)
}

# Converts all factor type variables in dataframe to character type.
fac2char <- function(df) {
  factorIndex <- sapply(df, is.factor)
  df[factorIndex] <- lapply(df[factorIndex], as.character)
  df
}

# Suppose you have many dataframes which you want to combine into a single dataframe.
# If there are variables unique to a single dataframe, then there is no way to
# calculate the pairwise covariance matrix. So, this function takes a list of
# datframes, determines which variables are shared between at least 2, then
# returns a single combined dataframe.
rosetta_bind <- function (x) {
  ## Find variables that are not shared with any other dataframe and those that
  ## are shared.
  col_names <- lapply(x, colnames)
  remove <- unlist(multi_diff(col_names))
  stay <- setdiff(unlist(col_names), remove)

  ## Remove variables that are not shared with any other dataframe.
  if (is.null(names(x))) { # If an unnamed list, then assign our own names
    names(x) <- 1:length(x)
  }
  removed_cols <- lapply(
    names(x),
    function(z) {x[[z]][, !(names(x[[z]]) %in% remove)]}
  )
  names(removed_cols) <- names(x)

  ## Add variables that are shared with other dataframes.
  for (i in 1:length(removed_cols)) {
    removed_cols[[i]][, setdiff(stay, colnames(removed_cols[[i]]))] <- NA
  }

  ## Bind dataframes that contain 2 or more shared variables.
  ret <- do.call("rbind", removed_cols)
  ret
}

# Performs a symmetric (both directions) setdiff() on a list of vectors.
multi_diff = function(x) {
  # Vector of all unique elements in x
  row_names <- sort(unique(unlist(x)))

  # Truth map. True if element i is in set j. Will apply matrix ops.
  map = as.matrix(as.data.frame(lapply(x, function(y){row_names %in% y})))
  rownames(map) <- row_names

  # Create non-conflicting column values
  col_vals = 2^seq(0, ncol(map) - 1)
  map = t(t(map) * col_vals)

  # Find row names that are not used anywhere else
  diff = lapply(col_vals, function(i) {
    names(which(rowSums(map) == i))
  })
  names(diff) = colnames(map)

  # return list of diffs
  diff
}
