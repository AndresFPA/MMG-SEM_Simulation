library(lavaan)
library(fpp3)
library(dplyr)
library(xtable)
library(ggpubr)
library(ggplot2)
library(qwraps2)
library(Cairo)

CairoWin()

# Set wd
setwd("C:/Users/User/OneDrive - Tilburg University/R/Simulation first paper/FinalSim2/Results")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/FinalSim2/Results")

# Load empty final results matrix
load("FinalResults.Rdata")
load("design.Rdata")

# Organize the data frame
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")
col_order <- c("Condition", "Replication", "NonInvIncl", "nclus", "ngroups", "coeff", "N_g",
               "balance", "reliability", "NonInvType", "NonInvSize", "ResRange", "NonInvItems", "NonInvG", "MisClass", 
               "ARI", "CorrectClus", "fARI", "RelBias_B1", "RelBias_B2", "RelBias_B3", "RelBias_B4", 
               "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "Global Max %", "Exo_var", "Exo_cov", "Decreasing")
Results_final <- Results_final[, col_order]
rm(col_order)

# Fill in the matrix with all results
ncond <- unique(Results_final$Condition) # How many conditions?
K <- length(unique(Results_final$Replication)) # How many replications?
for (i in ncond) {
  test <- NA
  test <- try(load(paste0("ResultRow", i, ".Rdata")))
  if(!c(class(test) == "try-error")){
    Results_final[(K*6*(i-1)+1):(i*K*6), 15:30] <- ResultsRow
  }
}

# Remove NA results
Results_final <- Results_final %>% filter(!is.na(RMSE_B1))
Results_final <- Results_final %>% filter(!is.na(MisClass)) # If I want to remove MG-SEM

# Remove unnecessary results
Results_final <- Results_final %>% filter(NonInvG == 0.5) 
# Results_final <- Results_final %>% filter(!c(NonInvG == 1 & NonInvType == "fixed")) 
# Results_final <- Results_final %>% filter(!c(NonInvG == .5 & NonInvType == "random")) 

# Save Results_final back-up for later
Results_final_back <- Results_final

# Re-order Non-Inv Incl
Results_final$NonInvIncl <- recode(Results_final$NonInvIncl, 
                                   both = "Model 1", 
                                   noneK = "Model 4",
                                   none = "Model 2",
                                   bothK = "Model 3")

Results_final$NonInvIncl <- factor(Results_final$NonInvIncl, levels = c("Model 1", 
                                                                        "Model 2", 
                                                                        "Model 3", 
                                                                        "Model 4",
                                                                        "mg_sem", "mg_sem_ign"))

# Check decreasing
sum(Results_final$Decreasing, na.rm = T) # No decreasing cases!

# Transform to numeric
# glimpse(Results_final)
# cl <- Results_final %>% select(MisClass:Exo_cov) %>% colnames()
# Results_final[, cl] <- Results_final %>% select(MisClass:Exo_cov) %>% mutate(across(where(is.character), ~ as.numeric(.x)))

# TRIAL GENERAL ------------------------------------------------------------------------------------
# Cluster Recovery
for_plots <- Results_final %>% group_by(N_g, coeff, NonInvIncl, NonInvType) %>% summarise(across(MisClass:Exo_var, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(ARI, fARI), names_to = "Measure", 
                               values_to = "ClusRecovery") %>% dplyr::select(Measure, ClusRecovery, N_g, coeff, NonInvIncl, NonInvType)

# for_plots_long <- for_plots_long %>% filter(coeff == 0.2)

CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ClusSim2.png", width = 650, height = 650)
CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ClusSim2.png", width = 650, height = 650)

ggplot(for_plots_long, aes(x = N_g, y = ClusRecovery, colour = Measure, linetype = NonInvType)) + 
  geom_line(size = 1) + 
  facet_grid(coeff ~ NonInvIncl) + labs(x = "Within-group Sample Size", y = "Cluster Recovery") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Model Specification", breaks = NULL, labels = NULL)) + 
  theme_bw() + theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_color_hue(labels = c("ARI", "fARI")) + scale_linetype_manual("Non-inv Type", values = c(1, 2), labels = c("Fixed", "Random")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Regression Parameters", breaks = NULL, labels = NULL), limits = c(0,1)) +
  ggtitle(label = "Cluster Recovery") + geom_point()

dev.off()

# # Parameter Recovery
# for_plots <- Results_final %>% group_by(N_g, coeff, NonInvIncl, NonInvType) %>% summarise(across(MisClass:Exo_var, mean))
# for_plots_long <- pivot_longer(data = for_plots, cols = c(RMSE_B1, RMSE_B2, RMSE_B3, RMSE_B4), names_to = "Parameter", 
#                                values_to = "RMSE") %>% dplyr::select(Parameter, RMSE, N_g, coeff, NonInvIncl, NonInvType)
# #for_plots_long <- for_plots_long %>% filter(coeff == 0.2)
# 
# CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)
# CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)
# 
# ggplot(for_plots_long, aes(x = N_g, y = RMSE, colour = Parameter, linetype = NonInvType)) + 
#   geom_line(size = 1) + 
#   facet_grid(coeff ~ NonInvIncl) + labs(x = "Within-group Sample Size", y = "Root Mean Squared Error (RMSE)") + 
#   scale_x_continuous(sec.axis = sec_axis(~ . , name = "Which non-invariance is included?", breaks = NULL, labels = NULL)) + 
#   theme_bw() + theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
#   scale_color_hue(labels = c(expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])), expression(paste(beta[4])))) +
#   scale_linetype_manual("Non-inv Type", values = c(1, 2), labels = c("Fixed", "Random")) +
#   scale_y_continuous(sec.axis = sec_axis(~ . , name = "Regression Parameters", breaks = NULL, labels = NULL)) +
#   ggtitle(label = "Parameter Recovery") + geom_point()
# 
# dev.off()

# Parameter Recovery (Only one RMSE)
Results_final$RMSE <- Results_final %>% select(RMSE_B1:RMSE_B4) %>% rowMeans()

for_plots <- Results_final %>% group_by(N_g, coeff, NonInvIncl, NonInvType) %>% summarise(across(MisClass:RMSE, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(RMSE), names_to = "Parameter", 
                               values_to = "RMSE") %>% dplyr::select(Parameter, RMSE, N_g, coeff, NonInvIncl, NonInvType)
#for_plots_long <- for_plots_long %>% filter(coeff == 0.2)

CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)
CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)

ggplot(for_plots_long, aes(x = N_g, y = RMSE, linetype = NonInvType)) + 
  geom_line(size = 1) + 
  facet_grid(coeff ~ NonInvIncl) + labs(x = "Within-group Sample Size", y = "Root Mean Squared Error (RMSE)") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Model Specification", breaks = NULL, labels = NULL)) + 
  theme_bw() + theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_linetype_manual("Non-inv Type", values = c(1, 2), labels = c("Fixed", "Random")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Regression Parameters", breaks = NULL, labels = NULL)) +
  ggtitle(label = "Parameter Recovery") + geom_point()

dev.off()

####################################################################################################
############################ TABLES - CLUSTER AND PARAMETER RECOVERY ###############################
####################################################################################################
Results_final <- Results_final_back
Results_final$RMSE <- Results_final %>% select(RMSE_B1:RMSE_B4) %>% rowMeans()

# Check mean results per simulation factor
a <- Results_final %>% group_by(NonInvIncl, nclus) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
b <- Results_final %>% group_by(NonInvIncl, ngroups) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
c <- Results_final %>% group_by(NonInvIncl, N_g) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
d <- Results_final %>% group_by(NonInvIncl, coeff) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
e <- Results_final %>% group_by(NonInvIncl, balance) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
f <- Results_final %>% group_by(NonInvIncl, reliability) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
g <- Results_final %>% group_by(NonInvIncl, NonInvSize) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
h <- Results_final %>% group_by(NonInvIncl, NonInvType) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
i <- Results_final %>% group_by(NonInvIncl, ResRange) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
j <- Results_final %>% group_by(NonInvIncl, NonInvG) %>% summarise(across(MisClass:RMSE, qwraps2::mean_sd, denote_sd = "paren", digits = 3))

list2 <- list(a, b, c, d, e, f, g, h, i, j)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[2]
  colnames(tmp)[1:2] <- c("Non-inv Included", "Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(1, ncol(tmp), c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g, h, i, j)
current <- current[, c("Factor", "Level", "Non-inv Included", "ARI", "CorrectClus", "fARI", 
                       "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "RMSE")]
# for_paper <- current %>% pivot_wider(names_from = `Non-inv Included`, values_from = c(ARI:RMSE_C)) %>% 
#   dplyr::select(Factor, Level, ARI_yes, CorrectClus_yes, fARI_yes, RMSE_A_yes, RMSE_B_yes, RMSE_C_yes)
for_paper <- current
colnames(for_paper) <- c("Factor", "Level", "Non-inv Included", "ARI", "CorrectClus", "fARI", "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "RMSE")

# Make wider
for_paper <- for_paper %>% filter(`Non-inv Included` %in% c("both", "bothK", "none", "noneK"))
for_paper <- pivot_wider(data = for_paper, names_from = `Non-inv Included`, values_from = ARI:RMSE)

# Re-organize table
# Get row indices of each factor
factors <- unique(for_paper$Factor)
for(i in 1:length(factors)){ 
  assign(x = paste0("rn_", factors[i]), value = which(for_paper$Factor == factors[i]))
}

for_paper <- for_paper[c(rn_coeff, rn_N_g, rn_nclus, rn_balance, rn_NonInvG,
                         rn_NonInvSize, rn_NonInvType, rn_ResRange), ]

rm(rn_coeff, rn_ngroups, rn_N_g, rn_nclus, rn_balance, rn_reliability, rn_NonInvG, rn_NonInvSize, rn_NonInvType, rn_ResRange)

# Add total - Cluster
yes_tot <- Results_final %>% filter(NonInvIncl %in% c("both", "bothK", "none", "noneK")) %>% group_by(NonInvIncl) %>%  
                             select(ARI, CorrectClus, fARI) %>% summarise(across(ARI:fARI, mean_sd, denote_sd = "paren")) %>% 
                             pivot_wider(names_from = "NonInvIncl", values_from = c("ARI", "CorrectClus", "fARI")) %>% 
                             select(ARI_both:fARI_noneK)

yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot %>% select(Factor, Level, ARI_both:fARI_noneK) %>% select(!contains("Correct"))

#Final - Cluster
for_paper_clus <- for_paper[, c("Factor", "Level", 
                                "ARI_both", "fARI_both",
                                "ARI_noneK", "fARI_noneK",
                                "ARI_bothK", "fARI_bothK",
                                "ARI_none", "fARI_none")]
for_paper_clus <- rbind(for_paper_clus, yes_tot)

# Final organization
for_paper_clus <- for_paper_clus[, c("Factor", "Level", 
                                "ARI_both", "fARI_both",
                                "ARI_none", "fARI_none",
                                "ARI_bothK", "fARI_bothK",
                                "ARI_noneK", "fARI_noneK")]

print(xtable(for_paper_clus, digits = 3), include.rownames = F)

rm(yes_tot, tmp)

# Add total - Parameter
yes_tot <- Results_final %>% filter(NonInvIncl %in% c("both", "bothK", "none", "noneK")) %>% group_by(NonInvIncl) %>%  
                             select(RMSE) %>% summarise(across(RMSE, mean_sd, denote_sd = "paren", digits = 3)) %>% 
                             pivot_wider(names_from = "NonInvIncl", values_from = c("RMSE"))
yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot %>% select(Factor, Level, both:noneK)
colnames(yes_tot) <- c("Factor", "Level", "RMSE_both", "RMSE_bothK", "RMSE_none", "RMSE_noneK")

#Final - Parameter
for_paper_par <- for_paper[, c("Factor", "Level", 
                               "RMSE_both", "RMSE_bothK", "RMSE_none", "RMSE_noneK"
                               )]
for_paper_par <- rbind(for_paper_par, yes_tot)

# Final organization
for_paper_par <- for_paper_par[, c("Factor", "Level", 
                                   "RMSE_both", "RMSE_none", "RMSE_bothK", "RMSE_noneK")]

# for_paper_par <- data.frame(lapply(for_paper_par, gsub, pattern = " ", replacement = "")) # Remove unnecessary spaces
print(xtable(for_paper_par, digits = 3), include.rownames = F)

# Extra - Average parameter recovery per regression parameter across all models
Results_final %>% ungroup() %>% select(RMSE_B1:RMSE_B4) %>% colMeans()


rm(yes_tot) 

Results_final %>% filter(NonInvG == 0.5) %>% group_by(NonInvIncl) %>%  
  select(RMSE_B1:RMSE_B4) %>% summarise(across(RMSE_B1:RMSE_B4, mean_sd, denote_sd = "paren", digits = 3))


