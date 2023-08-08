
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
  matt_diff <- mat - nearPD(mat_par, corr = TRUE, maxit = 500, conv.norm.type="F")[["mat"]]
  # Calculate Frobenius norm of difference matrix
  frob_norm <- sum(matt_diff^2, na.rm = TRUE)^(1/2)
  frob_norm
}

## 2. Find values which minimize the frobenius norm
val <- optim(
  par = par,
  mat = cov_mat,
  fn = sm,
  lower = -1,
  upper = 1,
  method = "L-BFGS-B"
)

val[["par"]]



######
## SEM function modified
######

sem_model <- function(factor_structure, sem_object = NULL) {
  message("starting sem_model")
  #factor_structure <- constrained_struc
  #sem_object <- unconstrained_fit
  ## check arguments
  if(length(names(factor_structure)) != length(factor_structure) || any(names(factor_structure) == "")) {
    stop("Check the 'factor_structure' argument in function rosetta:::sem_model(). 'factor_structure' needs to be a named list.")
  }
  ## Character values and inputs for SEM model
  vars <- stack(factor_structure)[["values"]]
  #message("vars")
  #message(vars)
  factors <- as.character(stack(factor_structure)[["ind"]])
  #message("factors")
  #message(factors)
  unique_factors <- unique(factors)
  #message("unique_factors")
  #message(unique_factors)
  n_factors <- length(factor_structure)
  #message("n_factors")
  #message(n_factors)
  n_vars <- length(vars)
  #message("n_vars")
  #message(n_vars)
  coefs <- paste("lam_", vars, sep = "")
  #message("coefs")
  #message(coefs)
  var_params <- paste("e_", vars, sep = "")
  #message("var_params")
  #message(var_params)
  fac_diag_values <- rep(1, n_factors)
  if (n_factors > 1) {
    factor_combs <- combn(unique_factors, 2)
  } else {
    factor_combs <- matrix(unique_factors)
  }
  #message("factor_combs")
  #message(factor_combs)
  n_combs <- ncol(factor_combs)
  #message("n_combs")
  #message(n_combs)
  fac_cor_params <- character(n_combs) # Initialize empty vector
  # message("fac_cor_params")
  #message(fac_cor_params)
  for (i in 1:n_combs) {
    #message("factor_combs[,i]")
    #message(factor_combs[,i])
    fac_cor_params[i] <- paste("f_", paste(factor_combs[,i], collapse = "_"), sep = "")
  }
  #message("fac_cor_params")
  #message(fac_cor_params)
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
  #message("fac_cor_values")
  #message(fac_cor_values)
  #message("fac_cor_params")
  #message(fac_cor_params)
  message("finished setting up character values and inputs")
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
  message("sorted out columns")
  message(left)
  message(arrow)
  message(right)
  message(param)
  message(value)
  ### Dataframe
  data <- data.frame(
    left = left,
    arrow = arrow,
    right = right,
    param = param,
    value = value,
    stringsAsFactors = FALSE
  )
  message("made dataframe")
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
  message("character vector created")
  # Return SEM model
  model_text <- paste(text, collapse = "\n")
  specified_model <- sem::specifyModel(text = model_text, quiet = TRUE)
  message("just ran sem-specify model")
  specified_model
}

