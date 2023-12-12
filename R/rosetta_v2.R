#' Combine imperfectly matched datasets
#'
#' Forms a complete dataset through latent class concatenation of imperfectly
#' matched dataset features.
#'
#' @param d A list of dataframes with imperfectly matched features.
#' @param factor_structure A named list. The list names are the factor names.
#' @param missing_corr 'normal'(default) or 'missing'. 'missing' tries to impute unobserved pairwise correlations
#' @param id_colnames optional names  of the columns that uniquely identify each row within a dataset (default is NULL)
#' Each element is a character vector of feature names for the corresponding factor.
#'
#' @return List of dataframes which contain factor scores.
#' @import optimParallel
#' @import parallel
#' @impoartFrom randtoolbox sobol
#' @import lavaan
#' @importFrom matrixcalc is.positive.definite
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
#'
#' # simulate data
#' d = sim()
#'
#' # check feature names
#' lapply(d$missing, names)
#'
#' # run rosetta
#' d_rosetta_missing = rosetta(
#'   d = d$missing,
#'   factor_structure = list(
#'     a = c("a_1", "a_2", "a_3"),
#'     b = c("b_1", "b_2", "b_3"),
#'     c = c("c_1", "c_2", "c_3")
#'   ),
#'   id_colnames = "ID"
#' )
#'
#' d_rosetta_complete = rosetta(
#'   d = d$complete,
#'   factor_structure = list(
#'     a = c("a_1", "a_2", "a_3"),
#'     b = c("b_1", "b_2", "b_3"),
#'     c = c("c_1", "c_2", "c_3")
#'   ),
#'   id_colnames = "ID"
#' )
#'


rosetta = function(d,
                   factor_structure,
                   missing_corr = 'normal',
                   id_colnames = NULL,
                   number_cores = 1) {
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

  # check all manifest vars in model are in the dataset
  vars_in_model = (factor_structure |> unlist())
  vars_in_model_test = rep(0,length(vars_in_model))
  for(i in seq_along(d)) {
    vars_in_model_test = vars_in_model_test +
      as.integer(vars_in_model %in% colnames(d[[i]]))
  }
  if(!all(vars_in_model_test > 0)){
    stop("manifest variables declared in factor structure but are not in data: ",
         paste(vars_in_model[vars_in_model_test==0], collapse = ", "))
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
    obs_cov = get_obs_cov(d_bind,
                          id_colnames)

  } else if (all(missing_corr=='missing')){
    message("Missing elements in covariance matrix. Using Steve's Matrix Imputation Algorithm")
    #===============================================================================
    # Algorithm for filling in missing correlations (S Buyske)
    #===============================================================================
    # 1. Define function to calculate frobenius norm of difference matrix.
    sm = function (par,mat) {
      # Store original correlation matrix and matrix that can be modified
      mat_par = mat
      # Get the missing value locations from the upper triangle
      mat_par[lower.tri(mat_par)] = 0
      index_na = which(is.na(mat_par), arr.ind = TRUE)
      # Restore modified matrix and assign values
      mat_par = mat

      # for (i in 1:nrow(index_na)) {
      #   mat_par[index_na[i,1], index_na[i,2]] = par[i]
      #   mat_par[index_na[i,2], index_na[i,1]] = par[i]
      # }

      # vectorizing so it's more slightly more efficient
      mat_par[index_na[,1], index_na[,2]] <- par
      mat_par[index_na[,2], index_na[,1]] <- par

      # Difference between original matrix and nearest positive definite chosen matrix
      matt_diff = mat - Matrix::nearPD(mat_par, corr = TRUE, maxit = 500, conv.norm.type="F")[["mat"]]
      # Calculate Frobenius norm of difference matrix
      frob_norm = sum(matt_diff^2, na.rm = TRUE)^(1/2)
      return(frob_norm)
    }

    data = dplyr::bind_rows(d)
    if(!is.null(id_colnames))data <- data[,-which(colnames(data) %in% id_colnames)]
    # head(data)
    #
    #     # # combined data (now we want NAs in columns)
    #     d_bind = rosetta_bind(d)

    ## observed pairwise complete covariance matrix
    cov_mat = get_obs_cov(data)

    cor_mat = stats::cov2cor(cov_mat)

    message("choosing initial values...")
    # initial values
    n_initial = length(which(is.na(cov_mat)))/2

    if(number_cores>1){
      par <- randtoolbox::sobol(n_initial)*2-1
      message("optimizing...")
      if(parallel::detectCores()<number_cores){
        stop("More cores requested for parallel processing than avaiable")
      }
      cl <- parallel::makeCluster(number_cores) # set the number of processor cores
      parallel::setDefaultCluster(cl=cl) # set 'cl' as default cluster
      val = optimParallel::optimParallel(
        par = par,
        mat = cor_mat,
        fn = sm,
        lower = -1,
        upper = 1,
        method = "L-BFGS-B")
      parallel::stopCluster(cl)
    } else {
      message("optimizing...")
      par = randtoolbox::sobol(n_initial)*2-1
      par = (randtoolbox::sobol(n_initial)*2-1)*.01
      # 2. Find values which minimize the frobenius norm
      val = stats::optim(
        par = par,
        mat = cor_mat,
        fn = sm,
        lower = -1,
        upper = 1,
        method = "L-BFGS-B"
      )
      if (val$convergence == 0) {
        message("Optimization for matrix imputation successful!")
      } else {
        warning("Optimization for matrix imputation failed to converge.")
      }

    }
    # 3. Put the estimated values back in original matrix
    #    There should be a better way...
    mat_optim = cor_mat
    mat_optim[lower.tri(mat_optim)] = 0
    index_na = which(is.na(mat_optim), arr.ind = TRUE)
    mat_optim = cor_mat
    for (i in 1:nrow(index_na)) {
      mat_optim[index_na[i,1], index_na[i,2]] = val[["par"]][i]
      mat_optim[index_na[i,2], index_na[i,1]] = val[["par"]][i]
    }

    mat_optim_cov = diag(sqrt(diag(cov_mat))) %*%
      mat_optim %*%
      diag(sqrt(diag(cov_mat)))

    colnames(mat_optim_cov) <- colnames(cor_mat)
    rownames(mat_optim_cov) <- rownames(cor_mat)
    if(! (all(eigen(mat_optim_cov)$values > 0) & isSymmetric(mat_optim_cov))){
      warning("after steve's matrix imputation algorithm, cov matrix is not positive semidefinite, attempting to coerce to positive semidefinite matrix")
      obs_cov = Matrix::nearPD(mat_optim_cov, corr = FALSE, maxit = 50000, conv.norm.type="F")$mat |> as.matrix()
    } else {
      obs_cov = mat_optim_cov
    }
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
    obs_cov = get_obs_cov(x,id_colnames)
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
      as.data.frame(lavaan::lavPredict(constrained_fit, newdata = x))

    # select ID rows
    if(!is.null(id_colnames)){

      id_df = as.data.frame(x[,id_colnames])
      colnames(id_df) <- id_colnames

      # connect IDs to factors scores
      constrained_factor_scores = dplyr::bind_cols(constrained_factor_scores,
                                                   id_df)
    }
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

    if(length(factor_structure[[i]])>0){
      factor_structure[[i]]
      fac_eq = paste0(factor_names[i]," =~ ",paste0(factor_structure[[i]],collapse = " + "))
      lavaan_text = c(lavaan_text, fac_eq)

      # constrain within-factor covariance to 1
      # 1* fixes a factor's covariance to 1
      factor_variance = paste0(factor_names[i]," ~~ 1*", factor_names[i])

      lavaan_text = c(lavaan_text, factor_variance)
    }
  }

  # if lavaan obj is passed, we'll assume we're doing the constrained fit
  # for the constrained fit, between factor covariance is fixed to factor covariance estimated from the entire data set (unconstrained estimate)
  if(!is.null(lavaan_obj)){
    fac_cov_est = get_fac_cov_estimates_lavaan(lavaan_obj,factor_structure)
    for(i in seq_along(factor_names)){
      for(j in i:(n_factors)){
        if(!j==i){
          if( (length(factor_structure[[i]])>0) &
              (length(factor_structure[[j]])>0) ){
            # x* fixes a factor's covariance to xs
            factor_covariance = paste0(factor_names[i]," ~~ ",
                                       fac_cov_est[(fac_cov_est$lhs==factor_names[i]) &
                                                     (fac_cov_est$rhs==factor_names[j]),"est"],"*", factor_names[j])

            lavaan_text = c(lavaan_text, factor_covariance)
          }
        }
      }
    }
  }


  return(lavaan_text)
}

# Returns the observed pairwise complete covariance matrix.
get_obs_cov = function(d, id_col = NULL) {
  d = d[, colSums(is.na(d)) < nrow(d)] # Remove columns which only contain NA

  if(!is.null(id_col))d = d[, -which(colnames(d) %in% id_col)]
  obs_cov = stats::cov(d, method = "pearson", use = "pairwise.complete.obs")
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
