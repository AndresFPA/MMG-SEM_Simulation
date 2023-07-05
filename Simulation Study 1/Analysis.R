library(lavaan)
library(qwraps2)
library(fpp3)
library(dplyr)
library(xtable)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(Cairo)
CairoWin()

# Set wd
setwd("C:/Users/User/OneDrive - Tilburg University/R/Simulation first paper/Final/Results")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/R/Simulation first paper/Final/Results")

# Load empty final results matrix
load("FinalResults.Rdata")
load("design.Rdata")

# Organize the data frame
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")
col_order <- c("Condition", "Replication", "NonInvIncl", "nclus", "ngroups", "coeff", "N_g",
               "balance", "reliability", "NonInvSize", "NonInvItems", "NonInvG", "MisClass", 
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
    Results_final[(K*2*(i-1)+1):(i*K*2), 13:28] <- ResultsRow
  }
}

# Remove NA results (come from an unused manipulated condition)
Results_final <- Results_final %>% filter(!is.na(MisClass))

# Remove 0.2 loading non-invariance size (irrelevant for the results)
Results_final <- Results_final %>% filter(NonInvSize != 0.2)

# Check decreasing
sum(Results_final$Decreasing, na.rm = T) # No decreasing cases!

# Transform to numeric
# glimpse(Results_final)
# cl <- Results_final %>% select(MisClass:Exo_cov) %>% colnames()
# Results_final[, cl] <- Results_final %>% select(MisClass:Exo_cov) %>% mutate(across(where(is.character), ~ as.numeric(.x)))

####################################################################################################
##################################### PLOTS - CLUSTER RECOVERY #####################################
####################################################################################################

# Sample size, coefficient, etc
for_plots <- Results_final %>% group_by(N_g, coeff, nclus, NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(ARI, fARI), names_to = "Measure", 
                               values_to = "ClusRecovery") %>% dplyr::select(Measure, ClusRecovery, N_g, coeff, nclus, NonInvIncl, balance)
for_plots_long <- for_plots_long %>% filter(balance == "balanced")

# Balanced clusters
a <- ggplot(for_plots_long, aes(x = N_g, y = ClusRecovery, colour = Measure)) + 
  geom_line(linewidth = 1.2) + 
  facet_grid(nclus ~ coeff) + labs(x = "Within-group Sample Size", y = "Cluster Recovery") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Size of Regression Parameter", breaks = NULL, labels = NULL)) + 
  #theme_bw() + 
  theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_color_hue(labels = c("ARI", "fARI")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of Clusters", breaks = NULL, labels = NULL)) +
  geom_point() +
  ggtitle(label = "Balanced Clusters")

# Unbalanced clusters
for_plots <- Results_final %>% group_by(N_g, coeff, nclus, NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(ARI, fARI), names_to = "Measure", 
                               values_to = "ClusRecovery") %>% dplyr::select(Measure, ClusRecovery, N_g, coeff, nclus, NonInvIncl, balance)

for_plots_long <- for_plots_long %>% filter(balance == "unbalanced")

b <- ggplot(for_plots_long, aes(x = N_g, y = ClusRecovery, colour = Measure)) + 
  geom_line(linewidth = 1) + 
  facet_grid(nclus ~ coeff) + labs(x = "Within-group Sample Size", y = "Cluster Recovery") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Size of Regression Parameter", breaks = NULL, labels = NULL)) + 
  #theme_bw() + 
  theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_color_hue(labels = c("ARI", "fARI")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of Clusters", breaks = NULL, labels = NULL)) +
  geom_point() + 
  ggtitle(label = "Unbalanced Clusters")

CairoPNG(filename = "E:/OneDrive - Tilburg University/First paper/ClusRecov1.png", width = 650, height = 700)
ggarrange(plotlist = list(a, b), common.legend = T, legend = "bottom", ncol = 1)
dev.off()

####################################################################################################
#################################### PLOTS - PARAMETER RECOVERY ####################################
####################################################################################################

CairoPNG(filename = "E:/OneDrive - Tilburg University/First paper/ParRecov1.png", width = 650, height = 700)

# Sample size, coefficient, etc
for_plots <- Results_final %>% group_by(N_g, coeff, nclus, NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(RMSE_B1, RMSE_B2, RMSE_B3, RMSE_B4), names_to = "Parameter", 
                               values_to = "RMSE") %>% dplyr::select(Parameter, RMSE, N_g, coeff, nclus, NonInvIncl, balance)
for_plots_long <- for_plots_long %>% filter(balance == "balanced")

# Balanced clusters
d <- ggplot(for_plots_long, aes(x = N_g, y = RMSE, colour = Parameter)) + 
  geom_line(size = 1) + 
  facet_grid(nclus ~ coeff) + labs(x = "Within-group Sample Size", y = "Root Mean Squared Error (RMSE)") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Size of Regression Parameter", breaks = NULL, labels = NULL)) + 
  theme_bw() +
  theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  #scale_color_hue(labels = c("Beta 2", "Beta 1", "Beta 3", "Beta 4")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of Clusters", breaks = NULL, labels = NULL), limits = c(0, 0.16)) +
  geom_point() +
  ggtitle(label = "Balanced Clusters")

# Unbalanced
for_plots <- Results_final %>% group_by(N_g, coeff, nclus, NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(RMSE_B1, RMSE_B2, RMSE_B3, RMSE_B4), names_to = "Parameter", 
                               values_to = "RMSE") %>% dplyr::select(Parameter, RMSE, N_g, coeff, nclus, NonInvIncl, balance)
for_plots_long <- for_plots_long %>% filter(balance == "unbalanced")

e <- ggplot(for_plots_long, aes(x = N_g, y = RMSE, colour = Parameter)) + 
  geom_line(size = 1) + 
  facet_grid(nclus ~ coeff) + labs(x = "Within-group Sample Size", y = "Root Mean Squared Error (RMSE)") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Size of Regression Parameter", breaks = NULL, labels = NULL)) + 
  theme_bw() +
  theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  #scale_color_hue(labels = c("Beta 2", "Beta 1", "Beta 3", "Beta 4")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of Clusters", breaks = NULL, labels = NULL), limits = c(0, 0.16)) +
  geom_point() +
  ggtitle(label = "Unbalanced Clusters")

ggarrange(plotlist = list(d, e), common.legend = T, legend = "bottom", ncol = 1)
dev.off()

####################################################################################################
############################ TABLES - CLUSTER AND PARAMETER RECOVERY ###############################
####################################################################################################

# Check mean results per simulation factor
a <- Results_final %>% group_by(NonInvIncl, nclus) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
b <- Results_final %>% group_by(NonInvIncl, ngroups) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
c <- Results_final %>% group_by(NonInvIncl, N_g) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
d <- Results_final %>% group_by(NonInvIncl, coeff) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
e <- Results_final %>% group_by(NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
f <- Results_final %>% group_by(NonInvIncl, reliability) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
g <- Results_final %>% group_by(NonInvIncl, NonInvSize) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
h <- Results_final %>% group_by(NonInvIncl, NonInvG) %>% summarise(across(MisClass:Exo_var, qwraps2::mean_sd, denote_sd = "paren", digits = 3))

list2 <- list(a, b, c, d, e, f, g, h)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[2]
  colnames(tmp)[1:2] <- c("Non-inv Included", "Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(1, ncol(tmp), c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g, h)
current <- current[, c("Factor", "Level", "ARI", "CorrectClus", "fARI", 
                       "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4")]
# for_paper <- current %>% pivot_wider(names_from = `Non-inv Included`, values_from = c(ARI:RMSE_C)) %>% 
#   dplyr::select(Factor, Level, ARI_yes, CorrectClus_yes, fARI_yes, RMSE_A_yes, RMSE_B_yes, RMSE_C_yes)
for_paper <- current
colnames(for_paper) <- c("Factor", "Level", "ARI", "CorrectClus", "fARI", "beta1", "beta2", "beta3", "beta4")

# Re-organize table
# Get row indices of each factor
factors <- unique(for_paper$Factor)
for(i in 1:length(factors)){ 
  assign(x = paste0("rn_", factors[i]), value = which(for_paper$Factor == factors[i]))
}

for_paper <- for_paper[c(rn_coeff, rn_ngroups, rn_N_g, rn_nclus,
                         rn_balance, rn_reliability, rn_NonInvG,
                         rn_NonInvSize), ]

rm(rn_coeff, rn_ngroups, rn_N_g, rn_nclus, rn_balance, rn_reliability, rn_NonInvG, rn_NonInvSize)

# Add total - Cluster
yes_tot <- t(apply(Results_final[, c("ARI", "CorrectClus", "fARI")], 2, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot[, c("Factor", "Level", "ARI", "CorrectClus", "fARI")]

#Final - Cluster
for_paper_clus <- for_paper[, c("Factor", "Level", "ARI", "CorrectClus", "fARI")]
for_paper_clus <- rbind(for_paper_clus, yes_tot)
print(xtable(for_paper_clus, digits = 3), include.rownames = F)

rm(yes_tot, tmp)

# Add total - Parameter
yes_tot <- t(apply(Results_final[, c("RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4")], 2, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot[, c("Factor", "Level", "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4")]
colnames(yes_tot) <- c("Factor", "Level", "beta1", "beta2", "beta3", "beta4")

#Final - Parameter
for_paper_par <- for_paper[, c("Factor", "Level", "beta1", "beta2", "beta3", "beta4")]
for_paper_par <- rbind(for_paper_par, yes_tot)
# for_paper_par[, c("beta1", "beta2", "beta3", "beta4")] <- round(for_paper_par[, c("beta1", "beta2", "beta3", "beta4")], digits = 3)

rm(yes_tot)

print(xtable(for_paper_par, digits = 3), include.rownames = F)

####################################################################################################
##################################### TABLE 3 - GLOBAL MAXIMA ######################################
####################################################################################################
Global <- Results_final
Global$`Global Max %` <- ifelse(test = Global$`Global Max %` == 0, yes = 0, no = 1)

# Check mean results per simulation factor
a <- Global %>% group_by(NonInvIncl, nclus) %>% summarise(across(MisClass:Exo_var, mean))
b <- Global %>% group_by(NonInvIncl, ngroups) %>% summarise(across(MisClass:Exo_var, mean))
c <- Global %>% group_by(NonInvIncl, N_g) %>% summarise(across(MisClass:Exo_var, mean))
d <- Global %>% group_by(NonInvIncl, coeff) %>% summarise(across(MisClass:Exo_var, mean))
e <- Global %>% group_by(NonInvIncl, balance) %>% summarise(across(MisClass:Exo_var, mean))
f <- Global %>% group_by(NonInvIncl, reliability) %>% summarise(across(MisClass:Exo_var, mean))
g <- Global %>% group_by(NonInvIncl, NonInvSize) %>% summarise(across(MisClass:Exo_var, mean))
h <- Global %>% group_by(NonInvIncl, NonInvG) %>% summarise(across(MisClass:Exo_var, mean))

list2 <- list(a, b, c, d, e, f, g, h)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[2]
  colnames(tmp)[1:2] <- c("Non-inv Included", "Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(1, ncol(tmp), c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g, h)

Global2 <- current[, c("Factor", "Level", "Global Max %")]
1 - mean(Global2$`Global Max %`)
mean(Global$`Global Max %`)

Local <- Global %>% filter(`Global Max %` == 0) %>% select(Condition:NonInvG, `Global Max %`)
Local %>% count(nclus)
Local %>% count(ngroups)
Local %>% count(N_g)
Local %>% count(coeff)
Local %>% count(balance)
Local %>% count(reliability)
Local %>% count(NonInvSize)
Local %>% count(NonInvG)

####################################################################################################
######################################## COR fARI - RMSE ###########################################
####################################################################################################

# Total correlations
Results_final %>% select(fARI, RMSE_B1:RMSE_B4) %>% cor()

# Correlation depending on the number of clusters
Results_final %>% filter(nclus == 2) %>% select(fARI, RMSE_B1:RMSE_B4) %>% cor()
Results_final %>% filter(nclus == 4) %>% select(fARI, RMSE_B1:RMSE_B4) %>% cor()









