library(lavaan)

# Set the working directory
setwd("C:/Users/User/OneDrive - Tilburg University/R/Simulation first paper/FinalSim2")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/FinalSim2")

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
ngroups <- 24
coeff <- c(0.2, 0.3, 0.4) # Initial regression parameters
N_g <- c(50, 100, 200) # Sample size per groups
balance <- c("balanced", "unbalanced")
reliability <- c("low")
NonInvSize <- c(0.4, 0.6)
NonInvItems <- 2
ResRange <- c(0.2, 0.6)
NonInvG <- c(0.5, 1)
NonInvType <- c("fixed", "random")

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
                      NonInvItems, ResRange, NonInvG, NonInvType)
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "reliability", "model",
                      "NonInvSize", "NonInvItems", "ResRange", "NonInvG", "NonInvType")

rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, ngroups, NonInvG, NonInvItems, NonInvSize, reliability, ResRange)

# Function for the simulation
do_sim <- function(Design, RowDesign, K){
  # Create matrix to store results
  # 8 columns for: MisClassError, ARI, RMSE, and Relative bias
  # There are 3 columns for RMSE and Relative bias (one for each regression parameter)
  ResultsRow <- matrix(data = NA, nrow = (K * 6), ncol = 16)
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
                              ResRange = Design[RowDesign, "ResRange"],
                              NonInvG = Design[RowDesign, "NonInvG"],
                              NonInvType = Design[RowDesign, "NonInvType"],
                              randomVarX = T)
    
    #Non-Inv Included?
    NonInv <- c("F1 =~ x2", "F1 =~ x3",
                "F2 =~ z2", "F2 =~ z3",
                "F3 =~ m2", "F3 =~ m3",
                "F4 =~ y2", "F4 =~ y3")
    # browser()
    # Run the model 4 times depending on the constraints
    # 1. BOTH RES AND LOAD NON-INV ARE INCLUDED
    resultsG_both <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                       group = "group", nclus = Design[RowDesign, "nclus"], nstarts = 20, 
                       printing = F, partition = "hard", NonInv = NonInv, seed = (RowDesign * k), 
                       constraints = "loadings", allG = T, fit = "factors")
    
    # Get global maxima
    global_maxG_both <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                         group = "group", nclus = Design[RowDesign, "nclus"], userStart = original, 
                         printing = F, nstarts = 1, NonInv = NonInv, constraints = "loadings", 
                         allG = T, fit = "factors")$loglikelihood
    
    # 2. RESIDUALS ARE CLUSTER-SPECIFIC
    resultsK_both <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                            group = "group", nclus = Design[RowDesign, "nclus"], nstarts = 20, 
                            printing = F, partition = "hard", NonInv = NonInv, seed = (RowDesign * k), 
                            constraints = "loadings", allG = F, fit = "factors")
    
    # Get global maxima
    global_maxK_both <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                               group = "group", nclus = Design[RowDesign, "nclus"], userStart = original, 
                               printing = F, nstarts = 1, NonInv = NonInv, constraints = "loadings", 
                               allG = F, fit = "factors")$loglikelihood
    
    # ------------------------------------------------------------
    # 3. BOTH NON-INV ARE !IGNORED!
    resultsG_none <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                            group = "group", nclus = Design[RowDesign, "nclus"], nstarts = 20, 
                            printing = F, partition = "hard", seed = (RowDesign * k), NonInv = NULL,
                            constraints = c("loadings", "residuals"), allG = T, fit = "factors")
    
    # Get global maxima
    global_maxG_none <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                              group = "group", nclus = Design[RowDesign, "nclus"], userStart = original, 
                              printing = F, nstarts = 1, NonInv = NULL, 
                              constraints = c("loadings", "residuals"), 
                              allG = T, fit = "factors")$loglikelihood
    
    # 4. BOTH NON-INV ARE !IGNORED!
    resultsK_none <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                            group = "group", nclus = Design[RowDesign, "nclus"], nstarts = 20, 
                            printing = F, partition = "hard", seed = (RowDesign * k), NonInv = NULL,
                            constraints = c("loadings", "residuals"), allG = F, fit = "factors")
    
    # Get global maxima
    global_maxK_none <- MMGSEM(dat = SimData$SimData, step1model = Measur_model, step2model = Struc_model, 
                              group = "group", nclus = Design[RowDesign, "nclus"], userStart = original, 
                              printing = F, nstarts = 1, NonInv = NULL, 
                              constraints = c("loadings", "residuals"), 
                              allG = F, fit = "factors")$loglikelihood
    
    
    # ---------------------------------------------------------------
    # MG-SEM
    # Ignored
    mg_sem_fit <- sem(model = model, data = SimData$SimData, group = "group")
    betas_mg_sem <- lapply(X = lavInspect(object = mg_sem_fit, what = "est"), FUN = "[[", "beta")
    
    # Included
    mg_sem_fit_ign <- sem(model = model, data = SimData$SimData, group = "group", group.equal = c("loadings", "residuals"))
    betas_mg_sem_ign <- lapply(X = lavInspect(object = mg_sem_fit_ign, what = "est"), FUN = "[[", "beta")
    
    # ---------------------------------------------------------------
    # Evaluate the results
    EvaluatedG_both <- evaluation(beta = resultsG_both$beta_ks, z_gks = resultsG_both$z_gks, original = original, 
                                 global_max = global_maxG_both, runs_LL = resultsG_both$runs_loglik,
                                 nclus = Design[RowDesign, "nclus"], coeff = Design[RowDesign, "coeff"],
                                 psi_gks = resultsG_both$psi_gks)
    
    EvaluatedK_both <- evaluation(beta = resultsK_both$beta_ks, z_gks = resultsK_both$z_gks, original = original, 
                                  global_max = global_maxK_both, runs_LL = resultsK_both$runs_loglik,
                                  nclus = Design[RowDesign, "nclus"], coeff = Design[RowDesign, "coeff"],
                                  psi_gks = resultsK_both$psi_gks)
    
    
    EvaluatedG_none <- evaluation(beta = resultsG_none$beta_ks, z_gks = resultsG_none$z_gks, original = original, 
                                  global_max = global_maxG_none, runs_LL = resultsG_none$runs_loglik,
                                  nclus = Design[RowDesign, "nclus"], coeff = Design[RowDesign, "coeff"],
                                  psi_gks = resultsG_none$psi_gks)
    
    EvaluatedK_none <- evaluation(beta = resultsK_none$beta_ks, z_gks = resultsK_none$z_gks, original = original, 
                                 global_max = global_maxK_none, runs_LL = resultsK_none$runs_loglik,
                                 nclus = Design[RowDesign, "nclus"], coeff = Design[RowDesign, "coeff"],
                                 psi_gks = resultsK_none$psi_gks)
    
    Evaluated_MGSEM <- evaluation(beta = betas_mg_sem, original = original, mg_sem = T, nclus = Design[RowDesign, "nclus"],
                                  coeff = Design[RowDesign, "coeff"])
    
    Evaluated_MGSEM_ign <- evaluation(beta = betas_mg_sem_ign, original = original, mg_sem = T, nclus = Design[RowDesign, "nclus"],
                                      coeff = Design[RowDesign, "coeff"])
    
    # Store the results
    ResultsRow[k, ] <- unlist(c(EvaluatedG_both, resultsG_both$Decreasing))
    ResultsRow[(K + k), ] <- unlist(c(EvaluatedK_both, resultsK_both$Decreasing))
    ResultsRow[((K*2) + k), ] <- unlist(c(EvaluatedG_none, resultsG_none$Decreasing))
    ResultsRow[((K*3) + k), ] <- unlist(c(EvaluatedK_none, resultsK_none$Decreasing))
    ResultsRow[((K*4) + k), 9:12] <- unlist(c(Evaluated_MGSEM))
    ResultsRow[((K*5) + k), 9:12] <- unlist(c(Evaluated_MGSEM_ign))
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Result", "Row", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(ResultsRow)
}

# Set working directory for the results
# Post-IMPS
setwd("C:/Users/User/OneDrive - Tilburg University/R/Simulation first paper/FinalSim2/Results")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/FinalSim2/Results")

# Create final results matrix 
# Everything is multiplied by 2 because we run the model twice (including and not including Non-Inv)
K <- 50 # Number of replications per condition

Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K*6, ncol = 16))
colnames(Results_final) <- c("MisClass", "ARI", "CorrectClus", "fARI",
                             "RelBias_B1", "RelBias_B2", "RelBias_B3", "RelBias_B4",
                             "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "Global Max %", 
                             "Exo_var", "Exo_cov", "Decreasing")
Results_final$Replication <- rep(x = 1:K, times = nrow(design)*6)
Results_final$Condition <- rep(x = 1:nrow(design), each = K*6)
Results_final$NonInvIncl <- rep(rep(x = c("both", "bothK", "none", "noneK", "mg_sem", "mg_sem_ign"), each = K), times = nrow(design))

system.time(for(i in 1:1){
  cat("\n", "Condition", i, "out of", nrow(design), "\n")
  Results <- do_sim(Design = design, RowDesign = i, K = K)
  Results_final[(K*6*(i-1)+1):(i*K*6), 1:16] <- Results
})

save(Results_final, file = "FinalResults.Rdata")
save(design, file = "design.Rdata")
