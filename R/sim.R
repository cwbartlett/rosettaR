#' Simulate data for rosetta
#'
#' This function simulates 'complete' and 'missing' datasets. For each independent
#' dataset, one variable per domain will be set to missing.
#'
#' @param loading The measurement model for x. The '\code{fx}' argument as found in \code{psych::\link[psych]{sim.structure}}.
#' @param correlation The structure matrix of the latent variables. The '\code{Phi}' argument as found in \code{psych::\link[psych]{sim.structure}}.
#' @param factor_structure A named list. The list names are the factor names.
#' Each element is a character vector of feature names for the corresponding factor.
#' Should be ordered corresponding to the rows of the '\code{loading}' argument.
#' @param n_rows An integer for the number of rows for each independent dataset.
#' @param n_datasets An integer for the number of independent datasets.
#' @param seed An integer for the seed.
#'
#' @return Returns a list that contains the complete data and missing data.
#'
#' @importFrom psych sim.structure
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Data simulation example
#' #----------------------------------------------------------------------------
#' # By default, sim() will simulate 3 variables from 3 different domains.
#' d_sim <- sim()
#' str(d_sim)
#' complete_data <- d_sim$complete
#' missing_data <- d_sim$missing
#'
sim <- function(
  loading = matrix(
    c(.9, .8, .7, rep(0, 9),
      .6, .7, .8, rep(0, 9),
      .8, .9, .6),
    ncol = 3
  ),
  correlation = matrix(
    c(1, .2, .4,
      .2, 1, .3,
      .4, .3, 1),
    ncol = 3
  ),
  factor_structure = list(
    a = c(1, 2, 3),
    b = c(1, 2, 3),
    c = c(1, 2, 3)
  ),
  n_rows = 1000,
  n_datasets = 3,
  seed = NULL
) {
  complete_data <- sim_complete(
    loading = loading,
    correlation = correlation,
    factor_structure = factor_structure,
    n_rows = n_rows,
    n_datasets = n_datasets,
    seed = seed
  )

  missing_data <- sim_missing(
    complete = complete_data,
    factor_structure = factor_structure
  )

  ret <- list(complete = complete_data, missing = missing_data)

  factor_struc <- lapply(names(factor_structure), function(x) {paste(x, factor_structure[[x]], sep = "_")})
  names(factor_struc) <- names(factor_structure)
  attr(ret, "factor_structure") <- factor_struc

  ret
}

#-------------------------------------------------------------------------------
# helper functions
#-------------------------------------------------------------------------------
# simulate a complete dataset based on specified factor loadings and correlation structure.
sim_complete <- function(loading, correlation, factor_structure, n_rows, n_datasets, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }

  true_correlation <- psych::sim.structure(
    fx = loading,
    Phi = correlation
  )$model

  sim_data <- as.data.frame(
    MASS::mvrnorm(
      n = n_rows * n_datasets,
      mu = rep(0, nrow(true_correlation)),
      Sigma = true_correlation
    )
  )

  names(sim_data) <- unlist(lapply(names(factor_structure), function(x) {paste(x, factor_structure[[x]], sep = "_")}))

  split_grouping <- cut(seq(1, nrow(sim_data)), breaks = n_datasets, labels = FALSE)
  sim_data_list <- split(x = sim_data, f = split_grouping)
  sim_data_list
}

# For each independent dataset, set a variable within each domain to missing
sim_missing <- function(complete, factor_structure) {
  missing <- complete

  column_names <- lapply(names(factor_structure), function(x) {paste(x, factor_structure[[x]], sep = "_")})
  remove_var <- lapply(column_names, sample)

  # breaks if there are more datasets than domains
  for(i in seq_along(missing)) {
    missing[[i]][sapply(remove_var, "[[", i)] <- NULL
  }

  missing
}
