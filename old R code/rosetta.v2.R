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
#' @importFrom sem sem fscores stdCoef
#' @importFrom DoE.wrapper lhs.design
#'
#' @importFrom Matrix nearPD
#' @importFrom easyr fac2char
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
#' d <- sim()
#'
#' # check feature names
#' lapply(d$missing, names)
#'
#' # run rosetta
#' d_rosetta <- rosetta(
#'    d = d$missing,
#'   factor_structure = list(
#'     a = c("a_1", "a_2", "a_3"),
#'     b = c("b_1", "b_2", "b_3"),
#'     c = c("c_1", "c_2", "c_3")
#'   )
#'  )
#'

rosetta.v3 <- function(d,
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
    d[[i]] <- Filter(function(x)!all(is.na(x)), d[[i]])
  }

  message(missing_corr)

  # step 1. unconstrained model
  ## sem RAM text
  lavaan_model <-
    get_lavaan_model_text(factor_structure)

  ## if the dataset has a measure not shared by at least two sub-datasets, then the correlation matrix will have missing values
  ## while we could test the dataset for this, instead we force the user to specify if they want the missing correlations
  ## filled in or not.  If we did it automatically, then users might not have intended to have missing data.  This way
  ## thier eyes are open when they go into this procedure, since they opted in.
  if(all(missing_corr=='normal')){
    # combined data (now we want NAs in columns)
    d_bind <- rosetta_bind(d)

    ## observed pairwise complete covariance matrix
    obs_cov <- get_obs_cov(d_bind)
  }
  else if (all(missing_corr=='missing')){
    message("Using Steve's Algorithm")
    #===============================================================================
    # Algorithm for filling in missing correlations (S Buyske)
    #===============================================================================
    # 1. Define function to calculate frobenius norm of difference matrix.
    sm <- function (mat, par) {
      # Store original correlation matrix and matrix that can be modified
      mat_par <- mat
      # Get the missing value locations from the upper triangle
      mat_par[lower.tri(mat_par)] <- 0
      index_na <- which(is.na(mat_par), arr.ind = TRUE)
      # Restore modified matrix and assign values
      mat_par <- mat
      for (i in 1:nrow(index_na)) {
        mat_par[index_na[i,1], index_na[i,2]] <- par[i]
        mat_par[index_na[i,2], index_na[i,1]] <- par[i]
      }
      # Difference between original matrix and nearest positive definite chosen matrix
      matt_diff <- mat - Matrix::nearPD(mat_par, corr = TRUE, maxit = 500, conv.norm.type="F")[["mat"]]
      # Calculate Frobenius norm of difference matrix
      frob_norm <- sum(matt_diff^2, na.rm = TRUE)^(1/2)
      frob_norm
    }

    data <- dplyr::bind_rows(d)
    head(data)
    cov_mat <- cov(data, use = "pairwise.complete.obs")
    cor_mat <- cov2cor(cov_mat)

    # initial values
    n_initial <- length(which(is.na(cov_mat)))/2
    par <- DoE.wrapper::lhs.design(n_initial, nfactors = 1, default.levels = c(-1, 1))[[1]]

    # 2. Find values which minimize the frobenius norm
    val <- optim(
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
    mat_optim <- cov_mat
    mat_optim[lower.tri(mat_optim)] <- 0
    index_na <- which(is.na(mat_optim), arr.ind = TRUE)
    mat_optim <- cov_mat
    for (i in 1:nrow(index_na)) {
      mat_optim[index_na[i,1], index_na[i,2]] <- val[["par"]][i]
      mat_optim[index_na[i,2], index_na[i,1]] <- val[["par"]][i]
    }
    matrixcalc::is.positive.definite(mat_optim)
    obs_cov <- mat_optim
  }

  ## the overall lavaan fit
  unconstrained_fit <- lavaan::cfa(model = lavaan_model,
                                   sample.cov = obs_cov,
                                   sample.nobs = nrow(data))

  ## factor covariance estimates
  unconstrained_factor_cov <-
    get_fac_cov_estimates_lavaan(unconstrained_fit,factor_structure)

  # step 2. constrained model
  constrained_fit_list <- lapply(d, function(x) {
    constrained_struc <- lapply(factor_structure, function(y) {
      intersect(names(x), y)
    })

    # the constrained sem RAM text
    lavaan_model_c <- get_lavaan_model_text(constrained_struc, unconstrained_fit)

    # observed pairwise complete covariance matrix
    obs_cov <- get_obs_cov(x)

    ## the overall lavaan fit
    constrained_fit <- lavaan::cfa(model = lavaan_model_c,
                                    sample.cov = obs_cov,
                                    sample.nobs = nrow(data))

    # model results
    constrained_factor_scores <- lavaan::lavPredict(constrained_fit, newdata = x)

    constrained_factor_cov <- get_fac_cov_estimates_lavaan(constrained_fit,
                                                           constrained_struc)

    list(
      constrained_fit = constrained_fit,
      constrained_factor_scores = constrained_factor_scores,
      constrained_factor_cov = constrained_factor_cov
    )
  })

  ret <- lapply(constrained_fit_list, `[[`, "constrained_factor_scores")
  #ret <- as.data.frame(do.call("rbind", ret))

  attr(ret, "unconstrained_fit") <- unconstrained_fit
  attr(ret, "factor_covariance") <- unconstrained_factor_cov
  attr(ret, "constrained_fit") <- lapply(constrained_fit_list, `[[`, "constrained_fit")

  ret
}

# Returns a character vector of the 'RAM' model for rosetta
sem_model <- function(factor_structure, sem_object = NULL) {
  ## check arguments
  if(length(names(factor_structure)) != length(factor_structure) || any(names(factor_structure) == "")) {
    stop("Check the 'factor_structure' argument in function rosetta:::sem_model(). 'factor_structure' needs to be a named list.")
  }
  ## Character values and inputs for SEM model
  vars <- stack(factor_structure)[["values"]]
  factors <- as.character(stack(factor_structure)[["ind"]])
  unique_factors <- unique(factors)
  n_factors <- length(factor_structure)
  n_vars <- length(vars)
  coefs <- paste("lam_", vars, sep = "")
  var_params <- paste("e_", vars, sep = "")
  fac_diag_values <- rep(1, n_factors)
  if (n_factors > 1) {
    factor_combs <- combn(unique_factors, 2)
  } else {
    factor_combs <- matrix(unique_factors)
  }
  n_combs <- ncol(factor_combs)
  fac_cor_params <- character(n_combs) # Initialize empty vector
  for (i in 1:n_combs) {
    fac_cor_params[i] <- paste("f_", paste(factor_combs[,i], collapse = "_"), sep = "")
  }
  if(is.null(sem_object)) { # Set to NA for unconstrained model
    fac_cor_values <- rep(NA, n_combs)
  } else { # Set to estimated values for constrained model
    if (n_factors > 1) {
      fac_cor_values <- fac_cov_estimates(sem_object) # Store dataframe of factor covariance estimates
      fac_cor_values <- fac_cor_values[match(fac_cor_params, fac_cor_values$covariance), ] # Make sure the covariance values are matched/ordered to the correct covariance
      fac_cor_values <- na.omit(fac_cor_values) # Remove any factor covariances that are not needed/matched
      fac_cor_values <- fac_cor_values$estimate
      fac_cor_params <- rep(NA, n_combs)
    }
  }

  ## Create the dataframe
  ### Columns
  if (n_factors > 1) {
    left <- c(factors, vars, unique_factors, factor_combs[1,])
    right <- c(vars, vars, unique_factors, factor_combs[2,])
    arrow <- c(rep("->", n_vars), rep("<->", n_vars), rep("<->", n_factors), rep("<->", n_combs))
    param <- c(coefs, var_params, rep(NA, n_factors), fac_cor_params)
    value <- c(rep(NA, n_vars), rep(NA, n_vars), fac_diag_values, fac_cor_values)
  } else {
    left <- c(factors, vars, unique_factors)
    right <- c(vars, vars, unique_factors)
    arrow <- c(rep("->", n_vars), rep("<->", n_vars), rep("<->", n_factors))
    param <- c(coefs, var_params, rep(NA, n_factors))
    value <- c(rep(NA, n_vars), rep(NA, n_vars), fac_diag_values)
  }

  ### Dataframe
  data <- data.frame(
    left = left,
    arrow = arrow,
    right = right,
    param = param,
    value = value,
    stringsAsFactors = FALSE
  )

  ### Create Empty character vector that will be filled. This avoids "writing" to disk.
  if (n_factors > 1) {
    text <- character(n_vars * 2 + n_factors + n_combs)
  } else {
    text <- character(n_vars * 2 + n_factors)
  }
  ### Fill in character vector
  for (i in 1:length(text)) { # Specify regression coefficients
    text[i] <- paste(
      paste(as.character(data[i,1:3]), collapse = " "),
      paste(as.character(data[i,4:5]), collapse = ", "),
      sep = ", "
    )
  }

  # Return SEM model
  model_text <- paste(text, collapse = "\n")
  specified_model <- sem::specifyModel(text = model_text, quiet = TRUE)
  specified_model
}

# Returns a character vector of the 'RAM' model for rosetta
get_lavaan_model_text <- function(factor_structure,lavaan_obj=NULL) {

  ## check arguments
  if(length(names(factor_structure)) != length(factor_structure) || any(names(factor_structure) == "")) {
    stop("Check the 'factor_structure' argument in function rosetta:::sem_model(). 'factor_structure' needs to be a named list.")
  }

  x =  factor_structure
  factor_names <- names(x)
  n_factors <- length(factor_names)
  lavaan_text = NULL
  for(i in seq_along(factor_names)){
    factor_structure[[i]]
    fac_eq = paste0(factor_names[i]," =~ ",paste0(factor_structure[[i]],collapse = " + "))
    lavaan_text = c(lavaan_text, fac_eq)

    # constrain within-factor covariance to 1
    factor_variance = paste0(factor_names[i]," ~~ 1*", factor_names[i])

    lavaan_text = c(lavaan_text, factor_variance)
  }

  if(!is.null(lavaan_obj)){
    fac_cov_est = get_fac_cov_estimates_lavaan(lavaan_obj,factor_structure)
    for(i in seq_along(factor_names)){
      for(j in i:(n_factors)){
        if(!j==i){
          # constrain factor variance to 1
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
get_obs_cov <- function(d) {
  d <- d[, colSums(is.na(d)) < nrow(d)] # Remove columns which only contain NA
  obs_cov <- cov(d, method = "pearson", use = "pairwise.complete.obs")
  obs_cov
}

# Returns factor covariance estimates from an sem object.
fac_cov_estimates <- function(sem) {
  n_factors <- sem$m - sem$n
  coef <- easyr::fac2char(sem::stdCoef(sem))
  n_coef <- nrow(sem::stdCoef(sem))
  n_combins <- ncol(combn(n_factors,2))
  factor_coef <- coef[(n_coef-n_combins+1):n_coef, 1:2]
  colnames(factor_coef) <- c("covariance", "estimate")
  rownames(factor_coef) <- NULL
  factor_coef
}


get_fac_cov_estimates_lavaan <- function(lavaan_model_obj,factor_structure) {
  pars_est = lavaan_model_obj |> lavaan::parameterestimates()
  factor_names = names(factor_structure)
  n_factors = length(factor_names)
  pars_est = pars_est[(pars_est$lhs %in% factor_names) &
                        (pars_est$rhs %in% factor_names) &
                        (pars_est$op %in% '~~'),
                      c('lhs','op','rhs','est')]
  return(pars_est)
}
