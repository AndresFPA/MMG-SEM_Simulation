library(lavaan)

# Set the working directory
setwd("C:/Users/User/OneDrive - Tilburg University/R/Simulation first paper/Final/Measurement model")
setwd("E:/OneDrive - Tilburg University/R/Simulation first paper/Final/Measurement model")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/Final/Measurement model")

# Source the relevant functions
source("DataGeneration.R")
source("evaluation.R")

# If we want to use the completely non-iterative model
# source("MMGSEM_noniter2.R")
# source("mgcfa_noniter_gls.R")

# Simulation Design
# Which factors are going to be tested? For now:
nclus <- c(2, 4) # Number of clusters
ngroups <- c(12, 24, 48) # Number of groups
coeff <- c(0.2, 0.3, 0.4) # Initial regression parameters
N_g <- c(50, 100, 200) # Sample size per groups
balance <- c("balanced", "unbalanced")
reliability <- c("high", "low")
NonInvSize <- c(0.2, 0.4, 0.6)
NonInvItems <- 2
NonInvG <- c(0.25, 0.50, 0.75)

model <- '
    # factor loadings
    F1 =~ x1 + x2 + x3 + x4 + x5
    F2 =~ z1 + z2 + z3 + z4 + z5
    F3 =~ m1 + m2 + m3 + m4 + m5
    F4 =~ y1 + y2 + y3 + y4 + y5
    
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'
Measur_model <- '
    # factor loadings
    F1 =~ x1 + x2 + x3 + x4 + x5
    F2 =~ z1 + z2 + z3 + z4 + z5
    F3 =~ m1 + m2 + m3 + m4 + m5
    F4 =~ y1 + y2 + y3 + y4 + y5
'

Struc_model <- '
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'

# Get design matrix
design <- expand.grid(nclus, ngroups, coeff, N_g, balance, reliability, model, NonInvSize, 
                      NonInvItems, NonInvG)
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "reliability", "model",
                      "NonInvSize", "NonInvItems", "NonInvG")

rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, ngroups, NonInvG, NonInvItems, NonInvSize, reliability)

# Function for the simulation
do_sim <- function(Design, RowDesign, K){
  # Create matrix to store results
  # 8 columns for: MisClassError, ARI, RMSE, and Relative bias
  # There are 3 columns for RMSE and Relative bias (one for each regression parameter)
  ResultsRow <- matrix(data = NA, nrow = (K), ncol = 2)
  colnames(ResultsRow) <- c("RMSE_lambda", "RMSE_theta")
  
  for(k in 1:K){
    #print(""); print(paste("Replication", k, "out of", K)); print("")
    print(paste("Replication", k, "out of", K))
    # Set seed per design condition (row) and replication (K)
    set.seed(RowDesign * k)
    
    # Generate data
    #SimData <- do.call(what = DataGeneration, args = Design[RowDesign, ])$SimData
    SimData <- DataGeneration(model = Design[RowDesign, "model"], 
                              nclus = Design[RowDesign, "nclus"], 
                              ngroups = Design[RowDesign, "ngroups"], 
                              reg_coeff = Design[RowDesign, "coeff"], 
                              N_g = Design[RowDesign, "N_g"], 
                              balance = Design[RowDesign, "balance"], 
                              reliability = Design[RowDesign, "reliability"], 
                              NonInvSize = Design[RowDesign, "NonInvSize"], 
                              NonInvItems = Design[RowDesign, "NonInvItems"], 
                              NonInvG = Design[RowDesign, "NonInvG"],
                              randomVarX = T)
    
    # Save data for replication purposes
    # save(SimData, file = paste0("Data/Cond", RowDesign, "Rep", k))
    
    #Non-Inv Included?
    NonInvIncl <- Design[RowDesign, "NonInvIncl"]
    NonInv <- c("F1 =~ x2", "F1 =~ x3",
                "F2 =~ z2", "F2 =~ z3",
                "F3 =~ m2", "F3 =~ m3",
                "F4 =~ y2", "F4 =~ y3")
    
    # Pre-center the data
    g_name <- as.character(unique(SimData$SimData[, "group"]))
    vars <- lavNames(lavaanify(Measur_model, auto = TRUE))
    lat_var <- lavNames(lavaanify(Measur_model, auto = TRUE), "lv")
    n_var <- length(vars)
    
    # # Center the data per group (so that the mean for all variables in each group is 0)
    # # Base R version for centering
    centered <- SimData$SimData
    group.idx <- match(SimData$SimData[,"group"], g_name)
    group.sizes <- tabulate(group.idx)
    group.means <- rowsum.default(as.matrix(SimData$SimData[,vars]),
                                  group = group.idx, reorder = FALSE,
                                  na.rm = FALSE)/group.sizes
    centered[,vars] <- SimData$SimData[,vars] - group.means[group.idx, drop = FALSE]
    
    # Add try-catch function in case any error appears. 
    # Non-invariances INCLUDED in the model
    # Run the model
    results <- lavaan::cfa(
      model = Measur_model, data = centered, group = "group",
      estimator = "ML", group.equal = "loadings",
      se = "none", test = "none", baseline = FALSE, h1 = FALSE,
      implied = FALSE, loglik = FALSE, 
      meanstructure = FALSE, group.partial = NonInv
    )
    
    EST <- lavInspect(results, "est") # Estimated measurement parameters
    lambda_est <- lapply(X = EST, "[[", "lambda")
    theta_est <- lapply(X = EST, "[[", "theta")
    
    # Evaluate the results
    # In case we got any error from running the model, fill the results with NA, and continue the sim
    if(class(results) == "try-error"){
      Evaluated <- NA
    } else {
      Evaluated <- evaluation(lambda_est = lambda_est, theta_est = theta_est, 
                              lambda = SimData$Lambda, theta = SimData$Theta, ngroups = Design[RowDesign, "ngroups"])
    }
    
    # Store the results
    ResultsRow[k, ] <- c(unlist(Evaluated)) 
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Result", "Row", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(ResultsRow)
}

# Set working directory for the results
# Post-IMPS
setwd("C:/Users/User/OneDrive - Tilburg University/R/Simulation first paper/Final/Measurement model/Results")
setwd("E:/OneDrive - Tilburg University/R/Simulation first paper/Final/Measurement model/Results")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/Final/Measurement model/Results")

# Create final results matrix 
# Everything is multiplied by 2 because we run the model twice (including and not including Non-Inv)
K <- 50 # Number of replications per condition
Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K, ncol = 2))
colnames(Results_final) <- c("RMSE_lambda", "RMSE_theta")
Results_final$Replication <- rep(x = 1:K, times = nrow(design))
Results_final$Condition <- rep(x = 1:nrow(design), each = K)

system.time(for(i in 60:250){
  cat("\n", "Condition", i, "out of", nrow(design), "\n")
  Results <- do_sim(Design = design, RowDesign = i, K = K)
  Results_final[(K*(i-1)+1):(i*K), 1:2] <- Results
})

save(Results_final, file = "FinalResults.Rdata")
save(design, file = "design.Rdata")
