library(MASS)
library(matrixcalc)

#' DataGeneration
#'
#' Generates simulated data according to the manipulated conditions. 
#'
#' INPUT: Arguments required by the function
#' @param model: full model in lavaan syntax (string).
#' @param nclus: number of clusters of the data. 
#' @param ngroups: number of groups for the generated data.
#' @param N_g: number of observations per group.
#' @param reg_coeff: regression coefficients of cluster 1. Must a numeric vector of length 3.
#' @param balance: can take values "balanced" and "unbalanced". Determines if the clusters will have the same
#'        number of groups.
#' @param reliability: reliability per item (manipulates ratio between loading and unique variance).
#' @param NonInvSize: size of the loading non-invariance (for the simulation: 0.2 and 0.4).
#' @param NonInvItems: number of items affected by the non-invariance (for now, only 2).
#' @param NonInvG: proportion of groups affected by the non-invariance (e.g., 0.25, 0.5).
#' @param randomVarX: define the variance of the exogenous variable. If TRUE, then random values will be defined
#'        for all groups. If FALSE, the variance will be 1 for all groups.
#'
#' OUTPUT
#' @return SimData: generated data
#' @return exog_vars: generated variance of the variables
#' @return exog_cov: generated covariance between exogenous variables
#' @return endo_vars: generated variance of the endogenous variables
#' @return Theta: generated theta matrix
#' @return NonInvIdx: Index of the groups that presented the non-invariant loadings
#' 

# Note:
# In the code, we use exog and endog to name the latent variables, while in the paper we use F1, F2, ..., etc.
# The corresponding matches are:
# F1 = exog2
# F2 = exog1
# F3 = endog1
# F4 = endog2

# Note 2:
# To ease the reading of the code, the names of the regression parameters can be seen below:
#   # exog1 -> endog2 = B1
#   # exog1 -> endog1 = B2
#   # exog2 -> endog1 = B3
#   # endog1 -> endog2 = B4


DataGeneration <- function(model, nclus, ngroups, N_g,
                           reg_coeff, balance,
                           reliability, NonInvSize,
                           NonInvItems = 2, NonInvG,
                           randomVarX = T){
  
  # Get number of variables
  par_table <- lavaanify(model)
  lat_var <- lavNames(par_table, "lv")
  obs_var <- lavNames(par_table, "ov")
  m <- length(lat_var)
  p <- length(obs_var)
  
  # Identify type of latent variable
  endog1 <- lat_var[(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]
  
  # Reorder latent variables
  lat_var <- c(exog, endog)
  
  # Give names to the regression parameters (only for better understanding):
  B1 <- numeric(nclus)
  B2 <- numeric(nclus)
  B3 <- numeric(nclus)
  B4 <- numeric(nclus)
  
  # Get regression parameters for each cluster
  # Cluster 1
  B1[1] <- 0
  B2[1] <- reg_coeff
  B3[1] <- reg_coeff
  B4[1] <- reg_coeff
  
  # Cluster 2
  B1[2] <- reg_coeff
  B2[2] <- 0
  B3[2] <- reg_coeff
  B4[2] <- reg_coeff
  
  if (nclus == 4){
    # Cluster 3
    B1[3] <- reg_coeff
    B2[3] <- reg_coeff
    B3[3] <- 0
    B4[3] <- reg_coeff
    
    # Cluster 4
    B1[4] <- reg_coeff
    B2[4] <- reg_coeff
    B3[4] <- reg_coeff
    B4[4] <- 0
  }

  # Generate the structural parameters
  # BETA - Regression parameters
  # beta is cluster-specific. Only nclus matrices needed
  beta <- array(data = 0, dim = c(m, m, nclus), dimnames = list(lat_var, lat_var))
  colnames(beta) <- rownames(beta) <- lat_var
  
  # Initialize
  for(k in 1:nclus){
    # beta
    beta[endog2, exog[1], k] <- B1[k]
    beta[endog1, exog[1], k] <- B2[k]
    beta[endog1, exog[2], k] <- B3[k]
    beta[endog2, endog1, k] <- B4[k]
  }
  
  # psi is group- and cluster-specific. ngroups matrices are needed.
  # Generate enough psi matrices
  psi_g <- array(data = diag(m), dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
  
  # Generate group-specific values for psi sampling from uniform distributions
  if (randomVarX == T){
    exog_var1 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_var2 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_cov <- runif(n = ngroups, min = -0.30, max = 0.30)
    endo_var1 <- runif(n = ngroups, min = 0.75, max = 1.25) # Total endogenous variance
    endo_var2 <- runif(n = ngroups, min = 0.75, max = 1.25) # Total endogenous variance
  } else { # DEPRECATED
    exog_var1 <- rep(1, times = ngroups) 
    exog_var2 <- rep(1, times = ngroups) 
    exog_cov <- rep(1, times = ngroups) 
  }
  
  # Get a cluster label for each group (GperK)
  # Example: If balanced, G = 6, K = 2. GperK = 111222
  if (balance == "balanced"){
    GperK <- rep(x = 1:nclus, each = (ngroups/nclus))
  } else if (balance == "unbalanced"){
    largest <- ngroups*.75; smaller <- ngroups - largest # largest and smaller clusters
    GperK <- c(rep(x = 1, times = largest), rep(x = 2:nclus, each = smaller/(nclus - 1)))
  }
  
  # Insert the corresponding group- and cluster-specific parts of psi
  for(g in 1:ngroups){
    # Insert group-specific parts
    psi_g[exog[1], exog[1], g] <- exog_var1[g] # Group-specific
    psi_g[exog[2], exog[2], g] <- exog_var2[g] # Group-specific
    psi_g[exog[1], exog[2], g] <- exog_cov[g] # Group-specific
    psi_g[exog[2], exog[1], g] <- exog_cov[g] # Group-specific
    
    # Insert the group-and-cluster-specific parts (they use the cluster-specific betas)
    # For the endogenous part, start from the total var (endog_var) and subtract the explained by the regression
    psi_g[endog1, endog1, g] <- endo_var1[g] - ((B2[GperK[g]]^2 * exog_var1[g]) + 
                                                (B3[GperK[g]]^2 * exog_var2[g]) + 
                                                (2 * B2[GperK[g]] * B3[GperK[g]] * exog_cov[g])) # Group-specific
    
    psi_g[endog2, endog2, g] <- endo_var2[g] - ((B1[GperK[g]]^2 * exog_var1[g]) + 
                                                (B4[GperK[g]]^2 * endo_var1[g]) + 
                                                (2 * B1[GperK[g]] * B4[GperK[g]] * ((B2[GperK[g]] * exog_var1[g]) + (B3[GperK[g]] * exog_cov[g])))) # Group-specific
  }
  
  # Create the covariance matrix of the factors (phi in the paper)
  I <- diag(m) # Identity matrix
  cov_eta <- array(data = 0, dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
  # browser()
  for(g in 1:ngroups){
    cov_eta[, , g] <- solve(I - beta[, , GperK[g]]) %*% psi_g[, , g] %*% solve(t(I - beta[, , GperK[g]]))
  }
  
  # Make sure that the average covariance per cluster is 1
  # Add (or remove) the necessary values
  endo_variances <- as.data.frame(cbind(cov_eta[endog1, endog1, ], cov_eta[endog2, endog2, ], GperK))
  colnames(endo_variances) <- c("F3", "F4", "Cluster")
  
  for(k in 1:nclus){
    this_k <- endo_variances$Cluster == k
    rsd_F3 <- endo_variances$F3[this_k]
    rsd_F4 <- endo_variances$F4[this_k]
    
    needed_F3 <- mean(1 - rsd_F3)
    needed_F4 <- mean(1 - rsd_F4)
    
    endo_variances$F3[this_k] <- endo_variances$F3[this_k] + needed_F3
    endo_variances$F4[this_k] <- endo_variances$F4[this_k] + needed_F4
    
    cov_eta[endog1, endog1, this_k] <- endo_variances$F3[this_k]
    cov_eta[endog2, endog2, this_k] <- endo_variances$F4[this_k]
  }
  
  # browser()
  # Check that all cov matrices are positive definite
  # for(g in 1:ngroups){
  #   test <- F
  #   while(test == F){
  #     cov_eta[, , g] <- solve(I - beta[, , GperK[g]]) %*% psi_g[, , g] %*% solve(t(I - beta[, , GperK[g]]))
  #     test <- is.positive.definite(round(cov_eta[, , g], digits = 12))
  #   }
  #   # print(is.positive.definite(cov_eta[, , g]))
  # }
  
  # Generate the measurement model parameters
  # Lambda (depends on reliability lvl)
  if (reliability == "low"){load <- .4} else if (reliability == "high"){load <- .6}
  loadings <- sqrt(load)
  
  # Non-invariance
  # Create Invariant Lambda
  Lambda <- matrix(data = rep(x = c(1, rep(loadings, 4), rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m)
  
  # How many items are affected?
  # Size of the non-invariance (NonInvSize) - Number of items affected (NonInvItems)
  # Create non-invariant lambda
  LambdaNonInv <- matrix(data = rep(x = c(1, 
                                          (loadings - NonInvSize), 
                                          (loadings + NonInvSize), 
                                          rep(loadings, (4 - NonInvItems)), 
                                          rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m)
  
  # Generate Theta per group depending on reliability
  Theta <- array(data = 0, dim = c(p, p, ngroups))
  
  # Using Kim's papers approach
  for (g in 1:ngroups){
    Theta[, , g] <- diag(runif(n = p, min = ((1 - load) - .1), max = ((1 - load) + .1)))
  }
  
  # Generate sample covariance matrix (sigma) per groups
  Sigma <- array(data = 0, dim = c(p, p, ngroups))
  
  # Include non-invariance in the loadings according to NonInvG
  # Which groups are non-invariant? (random sampling)
  NonInvIdx <- sample(x = 1:ngroups, size = NonInvG*ngroups, replace = F)
  for(g in 1:ngroups){
    # Non-invariance - How many groups?
    if (g %in% NonInvIdx){
      Sigma[, , g] <- LambdaNonInv %*% cov_eta[, , g] %*% t(LambdaNonInv) + Theta[, , g]
    } else if (!c(g %in% NonInvIdx)){
      Sigma[, , g] <- Lambda %*% cov_eta[, , g] %*% t(Lambda) + Theta[, , g]
    }
  }

  
  # Data Generation final
  # For now, mu would be 0 as we are only interested in centered variables
  SimData <- c()
  for(g in 1:ngroups){
    tmp <- mvrnorm(n = N_g, mu = rep(0, p), Sigma = Sigma[, , g], empirical = T)
    SimData <- rbind(SimData, tmp)
  }
  
  # Add the final labels
  group <- rep(x = c(1:ngroups), each = N_g) # Group variable
  SimData <- cbind(SimData, group)
  colnames(SimData) <- c(obs_var, "group")
  
  # Return data
  return(list(SimData = SimData, 
              exog_vars = list(exog_var1, exog_var2), 
              exog_cov = exog_cov, 
              endo_vars = list(endo_var1, endo_var2),
              Theta = Theta,
              NonInvIdx = NonInvIdx))
}
