#' Evaluation
#'
#' Evaluation function to assess the performance of the model according to the simulation results. 
#'
#' INPUT: Arguments required by the function
#' @param beta 
#' @param z_gks
#' @param original
#' @param global_max
#' @param runs_LL
#' @param nclus
#' @param reg_coef
#' @param reg_diff

library("combinat")
library("FARI")

evaluation <- function(beta, z_gks, original, global_max, runs_LL, nclus, coeff, psi_gks){
  # Transform posterior classification probabilities to hard clustering for the sake of the evaluation
  # clusResult <- ifelse(test = z_gks > .5, yes = 1, no = 0)
  clusResult <- t(apply(z_gks, 1, function(x) as.numeric(x == max(x))))
  colnames(clusResult) <- paste("Cluster", seq_len(nclus))
  
  # Organize columns in the same way as the original order (NOT USED NOW)
  # i.e. first cluster has the first x groups, second cluster the second x groups, etc.
  # clusResult <- clusResult[, names(sort(apply(X = clusResult, MARGIN = 2, FUN = which.max)))]
  
  # Get "true" cluster labels
  true_or <- original
  for (i in 1:ncol(original)){
    original[original[, i] != 0, i] <- i
  }
  or_vec <- c(original)[c(original) != 0]
  or_vec <- factor(x = or_vec, levels = 1:nclus) # Just in case there is an empty cluster
  
  # Evaluation 1 - Cluster Recovery ----------------------------------------------------------------
  # Cluster Recovery 1. Misclassification Error Rate
  
  # Check all possible permutations of cluster results order
  perm <- permn(x = ncol(clusResult))
  n_perm <- length(perm)
  
  # Initialize necessary objects
  permutedClusters <- vector(mode = "list", length = n_perm)
  clusVec <- vector(mode = "list", length = n_perm)
  MisClassError <- numeric(n_perm)
  
  # Check every permutation individually
  for(i in 1:n_perm){
    permutedClusters[[i]] <- clusResult[, perm[[i]]]
    
    # Label the groups depending on the cluster 
    # (i.e., groups in cluster 1 will have a 1, groups in cluster 2 will have a 2, etc.)
    for (j in 1:ncol(clusResult)){
      permutedClusters[[i]][permutedClusters[[i]][, j] != 0, j] <- j
    }
    
    # Re order columns by order of groups (NOT USED NOW)
    # permutedClusters[[i]] <- permutedClusters[[i]][, names(sort(apply(X = permutedClusters[[i]], MARGIN = 2, FUN = which.max)))]
    
    # Get a vector of cluster labels
    clusVec[[i]] <- rowSums(permutedClusters[[i]])
    
    # Transform into a factor to add all the necessary levels (in case some cluster does not have any groups)
    clusVec[[i]] <- factor(x = clusVec[[i]], levels = 1:nclus)
    
    # Create a confusion matrix and calculate misclassification error rate
    conf_mat <- table(clusVec[[i]], or_vec)
    MisClassError[i] <- 1 - (sum(diag(conf_mat))/sum(conf_mat))
  }
  
  # Get the result from the best permutation
  MisClassErrorOut <- min(MisClassError)
  
  # Which permutation is the best one? (will be used later)
  bestPerm <- permutedClusters[[which.min(MisClassError)]]
  bestPermVec <- clusVec[[which.min(MisClassError)]]
  #browser()
  # Cluster Recovery 2. Adjusted Rand Index (ARI)
  ARI_res <- adjrandindex(part1 = bestPermVec, part2 = or_vec)
  
  # Cluster Recovery 3. Correct clustering?
  CorrectClus <- all(bestPermVec == or_vec)
  
  # Cluster Recovery 4. Frobenius Adjusted Rand Index (fARI) "fuzzy ARI"
  fARI <- fari(z_gks, true_or)$fari
  
  # Evaluation 2 - Regression Parameters recovery --------------------------------------------------
  # Re order beta parameters in the same order as the original clusters
  beta <- beta[colnames(bestPerm)]
  # Extract the regressions from the beta matrix
  betaVec <- lapply(X = 1:nclus, FUN = function(x){c(beta[[x]])[c(beta[[x]]) != 0]})
  
  # Root Mean Squared Error (RMSE)
  SquaredErrorB1 <- numeric(nclus)
  SquaredErrorB2 <- numeric(nclus)
  SquaredErrorB3 <- numeric(nclus)
  SquaredErrorB4 <- numeric(nclus)
  
  # Relative bias
  RelativeBiasB1 <- numeric(nclus)
  RelativeBiasB2 <- numeric(nclus)
  RelativeBiasB3 <- numeric(nclus)
  RelativeBiasB4 <- numeric(nclus)
  
  # Calculate corresponding measures
  # Cluster 1
  SquaredErrorB1[1] <- (0 - betaVec[[1]][2])^2
  SquaredErrorB2[1] <- (coeff - betaVec[[1]][1])^2
  SquaredErrorB3[1] <- (coeff - betaVec[[1]][3])^2
  SquaredErrorB4[1] <- (coeff - betaVec[[1]][4])^2
  
  RelativeBiasB1[1] <- NA #True is 0
  RelativeBiasB2[1] <- (betaVec[[1]][1] - (coeff))/(coeff)
  RelativeBiasB3[1] <- (betaVec[[1]][3] - (coeff))/(coeff)
  RelativeBiasB4[1] <- (betaVec[[1]][4] - (coeff))/(coeff)
  
  # Cluster 2
  SquaredErrorB1[2] <- (coeff - betaVec[[2]][2])^2
  SquaredErrorB2[2] <- (0 - betaVec[[2]][1])^2
  SquaredErrorB3[2] <- (coeff - betaVec[[2]][3])^2
  SquaredErrorB4[2] <- (coeff - betaVec[[2]][4])^2
  
  RelativeBiasB1[2] <- (betaVec[[2]][2] - (coeff))/(coeff)
  RelativeBiasB2[2] <- NA #True is 0
  RelativeBiasB3[2] <- (betaVec[[2]][3] - (coeff))/(coeff)
  RelativeBiasB4[2] <- (betaVec[[2]][4] - (coeff))/(coeff)
  
  if (nclus == 4){
    # Cluster 3
    SquaredErrorB1[3] <- (coeff - betaVec[[3]][2])^2
    SquaredErrorB2[3] <- (coeff - betaVec[[3]][1])^2
    SquaredErrorB3[3] <- (0 - betaVec[[3]][3])^2
    SquaredErrorB4[3] <- (coeff - betaVec[[3]][4])^2
    
    RelativeBiasB1[3] <- (betaVec[[3]][2] - (coeff))/(coeff)
    RelativeBiasB2[3] <- (betaVec[[3]][1] - (coeff))/(coeff)
    RelativeBiasB3[3] <- NA #True is 0
    RelativeBiasB4[3] <- (betaVec[[3]][4] - (coeff))/(coeff)
    
    # Cluster 4
    SquaredErrorB1[4] <- (coeff - betaVec[[4]][2])^2
    SquaredErrorB2[4] <- (coeff - betaVec[[4]][1])^2
    SquaredErrorB3[4] <- (coeff - betaVec[[4]][3])^2
    SquaredErrorB4[4] <- (0 - betaVec[[4]][4])^2
    
    RelativeBiasB1[4] <- (betaVec[[4]][2] - (coeff))/(coeff)
    RelativeBiasB2[4] <- (betaVec[[4]][1] - (coeff))/(coeff)
    RelativeBiasB3[4] <- (betaVec[[4]][3] - (coeff))/(coeff)
    RelativeBiasB4[4] <- NA #True is 0
    
  }
  
  # RMSE
  RMSE_B1 <- sqrt(mean(SquaredErrorB1))
  RMSE_B2 <- sqrt(mean(SquaredErrorB2))
  RMSE_B3 <- sqrt(mean(SquaredErrorB3))
  RMSE_B4 <- sqrt(mean(SquaredErrorB4))
  
  # Relative Bias
  RelBias_B1 <- mean(RelativeBiasB1, na.rm = T)
  RelBias_B2 <- mean(RelativeBiasB2, na.rm = T)
  RelBias_B3 <- mean(RelativeBiasB3, na.rm = T)
  RelBias_B4 <- mean(RelativeBiasB4, na.rm = T)
  
  # Evaluation 3 - Proportion of correct results ----------------------------------------------------
  # Proportion of global maxima
  # Round the loglikelihoods to avoid differences of floating points
  ProporGlobalMax <- mean(round(runs_LL, 2) == round(global_max, 2)) #/length(runs_LL)
  
  # EXTRA: Average exogenous variable variance ------ 
  exo_mean <- mean(unlist(lapply(X = psi_gks[, 1], FUN = '[[', 1)))
  cov_mean <- mean(unlist(lapply(X = psi_gks[, 1], FUN = '[[', 2)))
  
  # Return results ---------------------------------------------------------------------------------
  return(list(ClusRecovery = list(MisClassError = MisClassErrorOut,
                                  ARI = ARI_res, 
                                  CorrectClus = CorrectClus, 
                                  fARI = fARI),
              ParamRecovery = list(RelativeBias = list(RelBias_B1 = RelBias_B1,
                                                       RelBias_B2 = RelBias_B2,
                                                       RelBias_B3 = RelBias_B3,
                                                       RelBias_B4 = RelBias_B4),
                                   RMSE = list(RMSE_B1 = RMSE_B1,
                                               RMSE_B2 = RMSE_B2,
                                               RMSE_B3 = RMSE_B3, 
                                               RMSE_B4 = RMSE_B4)), 
              ProporGlobalMax = ProporGlobalMax, 
              exo_mean = exo_mean,
              cov_mean = cov_mean))
}

# Function to create the original cluster matrix
create_original <- function(balance, ngroups, nclus){
  if (balance == "unbalanced"){
    unb <- c(rep(0, ngroups), rep(1, (ngroups*.25)/(nclus - 1)))
    original <- matrix(data = c(rep(1, ngroups*.75), rep(unb, nclus - 1)), nrow = ngroups, ncol = nclus)
  } else {
    data <- c(rep(x = c(rep(1, (ngroups/nclus)), rep(0, ngroups)), times = nclus))
    data <- data[-c((length(data)-ngroups+1):length(data))]
    original <- matrix(data = data, nrow = ngroups, ncol = nclus)
  }
  return(original)
}



# computation of adjusted rand index
adjrandindex <- function(part1,part2){
  part1 <- as.numeric(part1)
  part2 <- as.numeric(part2)
  
  IM1 <- diag(max(part1))
  IM2 <- diag(max(part2))
  A <- IM1[part1,]
  B <- IM2[part2,]
  
  T <- t(A)%*%B
  N <- sum(T)
  Tc <- apply(T,2,sum)
  Tr <- apply(T,1,sum)
  a <- (sum(T^2) - N)/2
  b <- (sum(Tr^2) - sum(T^2))/2
  c <- (sum(Tc^2) - sum(T^2))/2
  d <- (sum(T^2) + N^2 - sum(Tr^2) - sum(Tc^2))/2
  ARI <- (choose(N,2)*(a + d) - ((a+b)*(a+c)+(c+d)*(b+d)))/(choose(N,2)^2 - ((a+b)*(a+c)+(c+d)*(b+d)))
  
  return(ARI)
}










