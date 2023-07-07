#' Combine imperfectly matched datasets
#'
#' Forms a complete dataset through latent class concatenation of imperfectly
#' matched dataset features.
#'
#' @param d A list of dataframes with imperfectly matched features.
#' @param factor_structure A named list. The list names are the factor names.
#' Each element is a character vector of feature names for the corresponding factor.
#'
#' @return List of dataframes which contain factor scores.
#'
#' @import lavaan
#' @importFrom DoE.wrapper lhs.design
#' @importFrom Matrix nearPD
#' @importFrom dplyr bind_rows
#'
#'
#' @export
#'
#' @examples
#' #----------------------------------------------------------------------------
#' # Rosetta example
#' #----------------------------------------------------------------------------
#' library(rosetta)
#'
#' # simulate data
#' d = sim()
#'
#' # check feature names
#' lapply(d$missing, names)
#'
#' # run rosetta
#' d_rosetta = rosetta(
#'    d = d$missing,
#'   factor_structure = list(
#'     a = c("a_1", "a_2", "a_3"),
#'     b = c("b_1", "b_2", "b_3"),
#'     c = c("c_1", "c_2", "c_3")
#'   )
#'  )
#'



rosetta = function(d,
                   factor_structure,
                   missing_corr='normal') {
  # Check arguments
  if(!all(unlist(lapply(d, is.data.frame)))) {
    stop("Check the 'd' argument in function rosetta::rosetta(). 'd' needs to be a list of dataframes.")
  } # TODO: maybe a better test here (e.g. "is.null", all so check whether || or | is more appropriate)
  if(length(names(factor_structure)) != length(factor_structure) || any(names(factor_structure) == "")) {
    stop("Check the 'factor_structure' argument in function rosetta::rosetta(). 'factor_structure' needs to be a named list.")
  }
  if(length(names(factor_structure)) != length(factor_structure) || any(names(factor_structure) == "")) {
    stop("Check the 'missing_corr' argument in function rosetta::rosetta(). 'missing_corr' can only take on values 'normal' or 'missing'.")
  }
  # check if there are any fully NA columns in each dataset and remove them
  for(i in seq_along(d)) {
    d[[i]] = Filter(function(x)!all(is.na(x)), d[[i]])
  }

  message(missing_corr)

  # step 1. unconstrained model
  ## lavaan RAM text
  lavaan_model =
    get_lavaan_model_text(factor_structure)

  ## if the dataset has a measure not shared by at least two sub-datasets, then the correlation matrix will have missing values
  ## while we could test the dataset for this, instead we force the user to specify if they want the missing correlations
  ## filled in or not.  If we did it automatically, then users might not have intended to have missing data.  This way
  ## thier eyes are open when they go into this procedure, since they opted in.
  if(all(missing_corr=='normal')){
    # combined data (now we want NAs in columns)
    d_bind = rosetta_bind(d)

    ## observed pairwise complete covariance matrix
    obs_cov = get_obs_cov(d_bind)
  } else if (all(missing_corr=='missing')){
    message("Using Steve's Algorithm")
    #===============================================================================
    # Algorithm for filling in missing correlations (S Buyske)
    #===============================================================================
    # 1. Define function to calculate frobenius norm of difference matrix.
    sm = function (mat, par) {
      # Store original correlation matrix and matrix that can be modified
      mat_par = mat
      # Get the missing value locations from the upper triangle
      mat_par[lower.tri(mat_par)] = 0
      index_na = which(is.na(mat_par), arr.ind = TRUE)
      # Restore modified matrix and assign values
      mat_par = mat
      for (i in 1:nrow(index_na)) {
        mat_par[index_na[i,1], index_na[i,2]] = par[i]
        mat_par[index_na[i,2], index_na[i,1]] = par[i]
      }
      # Difference between original matrix and nearest positive definite chosen matrix
      matt_diff = mat - Matrix::nearPD(mat_par, corr = TRUE, maxit = 500, conv.norm.type="F")[["mat"]]
      # Calculate Frobenius norm of difference matrix
      frob_norm = sum(matt_diff^2, na.rm = TRUE)^(1/2)
      frob_norm
    }

    data = dplyr::bind_rows(d)
    head(data)
    cov_mat = cov(data, use = "pairwise.complete.obs")
    cor_mat = cov2cor(cov_mat)

    # initial values
    n_initial = length(which(is.na(cov_mat)))/2
    par = DoE.wrapper::lhs.design(n_initial, nfactors = 1, default.levels = c(-1, 1))[[1]]

    # 2. Find values which minimize the frobenius norm
    val = optim(
      par = par,
      mat = cov_mat,
      fn = sm,
      lower = -1,
      upper = 1,
      method = "L-BFGS-B"
    )
    val[["par"]]

    # 3. Put the estimated values back in original matrix
    #    There should be a better way...
    mat_optim = cov_mat
    mat_optim[lower.tri(mat_optim)] = 0
    index_na = which(is.na(mat_optim), arr.ind = TRUE)
    mat_optim = cov_mat
    for (i in 1:nrow(index_na)) {
      mat_optim[index_na[i,1], index_na[i,2]] = val[["par"]][i]
      mat_optim[index_na[i,2], index_na[i,1]] = val[["par"]][i]
    }
    matrixcalc::is.positive.definite(mat_optim)
    obs_cov = mat_optim
  }

  # get number of rows per data set
  n_data_list = lapply(d,function(x)nrow(x)) |> unlist()

  # avg number of rows per dataset
  n_data_mean = mean(n_data_list) |> floor()

  ## the overall lavaan fit
  unconstrained_fit =
    lavaan::cfa(model = lavaan_model,
                sample.cov = obs_cov,
                sample.nobs = n_data_mean,
                std.lv = TRUE)


  ## factor covariance estimates
  unconstrained_factor_cov =
    get_fac_cov_estimates_lavaan(unconstrained_fit,factor_structure)

  # step 2. constrained model
  constrained_fit_list = lapply(d, function(x) {
    constrained_struc = lapply(factor_structure, function(y) {
      intersect(names(x), y)
    })

    # the constrained model text
    lavaan_model_constrained =
      get_lavaan_model_text(constrained_struc,
                            unconstrained_fit)

    # observed pairwise complete covariance matrix
    obs_cov = get_obs_cov(x)
    tmp_unlist = unlist(constrained_struc,recursive = TRUE)
    obs_cov_subset = obs_cov[tmp_unlist,tmp_unlist]


    ## the overall lavaan fit
    constrained_fit =
      lavaan::cfa(model = lavaan_model_constrained,
                  sample.cov = obs_cov,
                  sample.nobs = n_data_mean,
                  std.lv = TRUE)

    # model results
    constrained_factor_scores =
      lavaan::lavPredict(constrained_fit, newdata = x)

    constrained_factor_cov =
      get_fac_cov_estimates_lavaan(constrained_fit,
                                   constrained_struc)

    list(
      constrained_fit = constrained_fit,
      constrained_factor_scores = constrained_factor_scores,
      constrained_factor_cov = constrained_factor_cov
    )
  })

  out = lapply(constrained_fit_list, `[[`, "constrained_factor_scores")

  attr(out, "unconstrained_fit_lavaan_object") = unconstrained_fit
  attr(out, "factor_covariance") = unconstrained_factor_cov
  attr(out, "constrained_fit_factor_Scores") = lapply(constrained_fit_list, `[[`, "constrained_fit")

  return(out)
}

# Returns a character vector of the 'RAM' model for rosetta
get_lavaan_model_text = function(factor_structure,lavaan_obj=NULL) {

  ## check arguments
  if(length(names(factor_structure)) != length(factor_structure) || any(names(factor_structure) == "")) {
    stop("Check the 'factor_structure' argument in function rosetta:::get_lavaan_model_text(). 'factor_structure' needs to be a named list.")
  }

  x =  factor_structure
  factor_names = names(x)
  n_factors = length(factor_names)

  lavaan_text = NULL
  for(i in seq_along(factor_names)){
    factor_structure[[i]]
    fac_eq = paste0(factor_names[i]," =~ ",paste0(factor_structure[[i]],collapse = " + "))
    lavaan_text = c(lavaan_text, fac_eq)

    # constrain within-factor covariance to 1
    # 1* fixes a factor's covariance to 1
    factor_variance = paste0(factor_names[i]," ~~ 1*", factor_names[i])

    lavaan_text = c(lavaan_text, factor_variance)
  }

  # if lavaan obj is passed, we'll assume we're doing the constrained fit
  # for the constrained fit, between factor covariance is fixed to factor covariance estimated from the entire data set (unconstrained estimate)
  if(!is.null(lavaan_obj)){
    fac_cov_est = get_fac_cov_estimates_lavaan(lavaan_obj,factor_structure)
    for(i in seq_along(factor_names)){
      for(j in i:(n_factors)){
        if(!j==i){

          # x* fixes a factor's covariance to xs
          factor_covariance = paste0(factor_names[i]," ~~ ",
                                     fac_cov_est[(fac_cov_est$lhs==factor_names[i]) &
                                                   (fac_cov_est$rhs==factor_names[j]),"est"],"*", factor_names[j])

          lavaan_text = c(lavaan_text, factor_covariance)
        }
      }
    }
  }

  return(lavaan_text)
}

# Returns the observed pairwise complete covariance matrix.
get_obs_cov = function(d) {
  d = d[, colSums(is.na(d)) < nrow(d)] # Remove columns which only contain NA
  obs_cov = cov(d, method = "pearson", use = "pairwise.complete.obs")
  obs_cov
}

# Extract covariance estimates from lavaan model fit
get_fac_cov_estimates_lavaan = function(lavaan_model_obj,factor_structure) {
  pars_est = lavaan_model_obj |> lavaan::parameterestimates()
  factor_names = names(factor_structure)
  n_factors = length(factor_names)
  pars_est = pars_est[(pars_est$lhs %in% factor_names) &
                        (pars_est$rhs %in% factor_names) &
                        (pars_est$op %in% '~~'),
                      c('lhs','op','rhs','est')]
  return(pars_est)
}

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

d_sim = sim(
  loading = matrix(
    c(.6, .6, .6, rep(0, 9),
      .6, .6, .6, rep(0, 9),
      .6, .6, .6),
    ncol = 3
  ),
  correlation = matrix(
    c(1, .3, .3,
      .3, 1, .3,
      .3, .3, 1),
    ncol = 3
  ),
  factor_structure = list(
    a = c(1, 2, 3),
    b = c(1, 2, 3),
    c = c(1, 2, 3)
  ),
  n_rows = 1000,
  n_datasets = 3,
  seed = 123)

d_complete = rosetta(d = d_sim$complete,
                     factor_structure = list(
                       a = c("a_1", "a_2", "a_3"),
                       b = c("b_1", "b_2", "b_3"),
                       c = c("c_1", "c_2", "c_3")
                     ),
                     missing_corr = "normal")

d_missing = rosetta(d_sim$missing,
                    factor_structure = list(
                      a = c("a_1", "a_2", "a_3"),
                      b = c("b_1", "b_2", "b_3"),
                      c = c("c_1", "c_2", "c_3")
                    ))

require(tidyverse)

com_1 =  d_complete$`1` |> as.data.frame()
colnames(com_1) <- paste0(colnames(com_1),"_complete")
mis_1 = d_missing$`1`|> as.data.frame()
colnames(mis_1) <- paste0(colnames(mis_1),"_missing")
d1 =  bind_cols(mis_1,com_1)
d1$dataset = rep("dataset_1",nrow(com_1))


com_2 =  d_complete$`2` |> as.data.frame()
colnames(com_2) <- paste0(colnames(com_2),"_complete")
mis_2 = d_missing$`2`|> as.data.frame()
colnames(mis_2) <- paste0(colnames(mis_2),"_missing")
d2 =  bind_cols(mis_2,com_2)
d2$dataset = rep("dataset_2",nrow(com_2))

com_3 =  d_complete$`3` |> as.data.frame()
colnames(com_3) <- paste0(colnames(com_3),"_complete")
mis_3 = d_missing$`3`|> as.data.frame()
colnames(mis_3) <- paste0(colnames(mis_3),"_missing")
d3 =  bind_cols(mis_3,com_3)
d3$dataset = rep("dataset_3",nrow(com_3))

plot_df = bind_rows(d1,d2,d3)

ggplot(plot_df)+geom_point(aes(x = a_missing,
                               y = a_complete,color = dataset))+
  ggtitle("factor score recovery, loading =.1, factor correlation = .1")

ggplot(plot_df)+geom_point(aes(x = b_missing,
                               y = b_complete,color = dataset))+
  ggtitle("factor score recovery, loading =.1, factor correlation = .1")

ggplot(plot_df)+geom_point(aes(x = c_missing,
                               y = c_complete,color = dataset))+
  ggtitle("factor score recovery, loading =.6, factor correlation = .3")
