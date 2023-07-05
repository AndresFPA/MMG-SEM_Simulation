library(lavaan)
library(semPlot)

# Set the working directory
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/Final")

# Source the relevant functions
source("MMGSEM.R")
source("E_Step.R")
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
  ResultsRow <- matrix(data = NA, nrow = (K * 2), ncol = 16)
  colnames(ResultsRow) <- c("MisClass", "ARI", "CorrectClus", "fARI",
                            "RelBias_B1", "RelBias_B2", "RelBias_B3", "RelBias_B4",
                            "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "Global Max %", 
                            "Exo_var", "Exo_cov", "Decreasing")
  
  # Create the original clustering matrix for comparison below
  original <- create_original(balance = Design[RowDesign, "balance"], 
                              ngroups = Design[RowDesign, "ngroups"], 
                              nclus = Design[RowDesign, "nclus"])
  
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
    
    # Add try-catch function in case any error appears. 
    # Non-invariances INCLUDED in the model
    # Run the model
    resultsK_F <- try(MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                             group = "group", nclus = Design[RowDesign, "nclus"], nstarts = 20, 
                             printing = F, partition = "hard", NonInv = NonInv, seed = (RowDesign * k), 
                             allG = T, fit = "factors"))
    
    # Run the model with original clustering to get proxy of the global maxima
    global_max_F <- try(MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                                 group = "group", nclus = Design[RowDesign, "nclus"], userStart = original, 
                                 printing = F, nstarts = 1, NonInv = NonInv, allG = T, fit = "factors")$loglikelihood)
    
    # Endogenous variance is group-specific
    resultsK_O <- NA
      
      # try(MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
      #                        group = "group", nclus = Design[RowDesign, "nclus"], nstarts = 20, 
      #                        printing = F, partition = "hard", NonInv = NonInv, seed = (RowDesign * k), 
      #                        allG = T, fit = "observed"))
    
    # Run the model with original clustering to get proxy of the global maxima
    global_max_O <- NA
      
      # try(MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
      #                          group = "group", nclus = Design[RowDesign, "nclus"], userStart = original, 
      #                          printing = F, nstarts = 1, NonInv = NonInv, allG = T, fit = "observed")$loglikelihood)
    
    # Evaluate the results
    # In case we got any error from running the model, fill the results with NA, and continue the sim
    if(class(resultsK_F) == "try-error"){
      Evaluated_F <- NA
    } else {
      Evaluated_F <- evaluation(beta = resultsK_F$beta_ks, z_gks = resultsK_F$z_gks, original = original, 
                                global_max = global_max_F, runs_LL = resultsK_F$runs_loglik,
                                nclus = Design[RowDesign, "nclus"], coeff = Design[RowDesign, "coeff"],
                                psi_gks = resultsK_F$psi_gks)
    }
    
    if(class(resultsK_O) == "try-error"){
      Evaluated_O <- NA
    } else {
      Evaluated_O <- NA 
      #evaluation(beta = resultsK_G$beta_ks, z_gks = resultsK_G$z_gks, original = original, 
      #                            global_max = global_max_G, runs_LL = resultsK_G$runs_loglik,
      #                            nclus = Design[RowDesign, "nclus"], coeff = Design[RowDesign, "coeff"],
      #                            psi_gks = resultsK_G$psi_gks)
    }
    
    # Store the results
    ResultsRow[k, ] <- c(unlist(Evaluated_F), resultsK_F$Decreasing) 
    ResultsRow[(K + k), ] <- NA
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Result", "Row", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(ResultsRow)
}

# Set working directory for the results
# Post-IMPS
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/Final/Results")

# Create final results matrix 
# Everything is multiplied by 2 because we run the model twice (including and not including Non-Inv)
K <- 50 # Number of replications per condition
Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K*2, ncol = 16))
colnames(Results_final) <- c("MisClass", "ARI", "CorrectClus", "fARI",
                             "RelBias_B1", "RelBias_B2", "RelBias_B3", "RelBias_B4",
                             "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "Global Max %", 
                             "Exo_var", "Exo_cov", "Decreasing")
Results_final$Replication <- rep(x = 1:K, times = nrow(design)*2)
Results_final$Condition <- rep(x = 1:nrow(design), each = K*2)
Results_final$NonInvIncl <- rep(rep(x = c("Factors", "Observed"), each = K), times = nrow(design))

system.time(for(i in 1:150){
  cat("\n", "Condition", i, "out of", nrow(design), "\n")
  Results <- do_sim(Design = design, RowDesign = i, K = K)
  Results_final[(K*2*(i-1)+1):(i*K*2), 1:15] <- Results
})

save(Results_final, file = "FinalResults.Rdata")
save(design, file = "design.Rdata")
