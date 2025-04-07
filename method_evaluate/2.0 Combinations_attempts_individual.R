library(purrr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)  # Load the gridExtra package
library(reshape2) # Needed for melting data frames for plotting
library(Matrix)
if (!require("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}
library(glmnet) # for LASSO


## Check pheno inculdes all 10K EUR sample (actually phenotype are often less than 10K)
# ID_10K <- read.table("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/EUR_10K.fam", header = F)[, 1]
# length(ID_10K, Trait_pheno)


# Step 0: Load PCs and covariates
cov_choice <- c("age_recruit", "sex", paste0("PC", 1:20))
total_covariates <- fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv") %>%
  select(all_of(c("eid", cov_choice)))
# Order covariates
total_covariates <- total_covariates[order(eid)]

## Set tuning dir
tune_dir <- "/gpfs/gibbs/pi/zhao/zb222/Data/ref_geno_data/EUR/4_fold_10K"

# Step 1: Gather necessary data
traits <- c("HDL", "LDL", "TC", "logTG")
n_traits <- length(traits)
pop <- "EUR"
epochs <- 1:4
annotations_info <- list(HDL = 83, LDL = 83, logTG = 83, TC = 83)
combined_ways <- c("JointPRS", "OLS", "LASSO", "Ridge")
# Initialize list to store coefficients
trait_coefficients <- list()
all_prs <- list()
ordered_annotations_enrich.p <- list()
ordered_annotations_prop.h <- list()
selected_annotations_prop.h <- list()
trait_coefficients <- list()
rank_test_result <- data.frame(
  trait = traits,
  tau_enrichmentP = numeric(n_traits),
  pvalue_enrichmentP = numeric(n_traits),
  tau_proph2 = numeric(n_traits),
  pvalue_proph2 = numeric(n_traits)
)
# Initialize combinations results dataframe
combined_ways_result <- data.frame(
  trait = rep(traits, each = length(epochs) * length(combined_ways)),
  epoch = rep(rep(epochs, each = length(combined_ways)), times = length(traits)),
  combined_way = factor(rep(combined_ways, times = length(traits) * length(epochs))),
  r_square = numeric(length(traits) * length(epochs) * length(combined_ways)),  # initialize with zeros
  stringsAsFactors = FALSE  # for other columns if needed
)
combined_ways_index <- 1
# Load binary annotation info
annot_info <- readLines("/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/annot_info/baselineLD_v2.2/baselineLD_v2.2_binary.txt")
annot_info <- paste0(annot_info, "L2_0")

for (trait in traits) {
  prs_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/original_sumstats/trait_", trait, "/tar_", pop)
  
  # Read fine-mapping and PRS files and assign column names
  K <- annotations_info[[trait]]
  
  # Arrange the order of annotations in ascending order of total heritability or enrichment_p
  heri_info_dir <- "/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/2.est_enrichment/result/"
  heri_info_raw <- read.table(paste0(heri_info_dir, trait, "_", pop, ".results"), header = T)
  heri_info_raw_binary <- heri_info_raw[heri_info_raw$Category %in% annot_info, ]
  
  ordered_annotations_enrich.p[[trait]] <- heri_info_raw_binary %>% 
    arrange(Enrichment_p) %>%
    pull(Category)
  
  ordered_annotations_prop.h[[trait]] <- heri_info_raw_binary %>% 
    arrange(Prop._h2) %>%
    pull(Category)
  
  selected_annotations_prop.h[[trait]] <- heri_info_raw_binary %>%
    arrange(desc(Prop._h2)) %>%
    filter(Prop._h2 > 0.5) %>%
    pull(Category)
  
  # Initialize a list to store data frames
  prs_data_list <- list()
  r_square_prs_results <- list()
  # Loop over the number of annotations
  for (k in 1:K) {
    # Construct the file name and path
    file_path <- paste0(prs_dir, "/baselineLD_v2.2_Annot_", k, "/", trait, "_", pop, "_PRScsx_baselineLD_v2.2_Annot_", k, ".sscore")
    
    # Read data using fread and select the required columns
    prs_data <- fread(file_path)[, .(IID = IID, PRS = SCORE1_AVG)]
    
    # Rename PRS column to reflect annotation name.
    setnames(prs_data, "PRS", ordered_annotations_enrich.p[[trait]][k])
    
    # Append to the list
    prs_data_list[[k]] <- prs_data
  }
  
  # Read in the Joint PRS data
  prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop, "/JointPRS/", trait, "_", pop, "_JointPRS.sscore"))[, .(IID = IID, PRS_Joint = SCORE1_AVG)]
  
  # Add Joint PRS to the list
  prs_data_list[[K + 1]] <- prs_JointPRS
  
  # Merge all data frames in the list by "IID"
  all_prs[[trait]] <- reduce(prs_data_list, full_join, by = "IID")
  
  # Loop to create a numeric vector for each annotation and one for JointPRS, one for linear combination
  for (k in 1:(K + 3)) {
    r_square_prs_results[[k]] <- numeric(4)  # 4 epochs
  }
  
  Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv"))
  setnames(Trait_pheno, c("eid", "pheno"))
  Trait_pheno[, pheno := scale(pheno)]
  
  trait_coefficients[[trait]] <- list()
  for (epoch in epochs) {
    Trait_tune_pheno_id <- fread(paste0(tune_dir, "/tune_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
    common_tune_id <- intersect(Trait_tune_pheno_id, Trait_pheno$eid)
    # Subset Covariates and phenos
    covariates_tune <- total_covariates[eid %in% common_tune_id]
    Trait_pheno_tune <- Trait_pheno[eid %in% common_tune_id]
    
    # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
    pheno_covariates_tune <- merge(Trait_pheno_tune, covariates_tune, by.x = "eid", by.y = "eid")
    df_prs_tune <- all_prs[[trait]][IID %in% common_tune_id]
    
    null_model_tune <- lm(pheno ~ ., data = pheno_covariates_tune[, -1])  # exclude eid
    adjusted_pheno_tune <- residuals(null_model_tune)
    # Scale all columns except the 'IID' column
    df_prs_tune_scaled <- df_prs_tune %>%
      mutate(across(-IID, scale))

    
    # # OLS approach
    # # Create the regression formula by collapsing the predictor names into a single string separated by '+'
    # formula_str_all <- paste("adjusted_pheno_tune", "~", paste(ordered_annotations_enrich.p[[trait]], collapse = " + "))
    # # Add the JointPRS result
    # formula_str_all <- paste(formula_str_all, "+ PRS_Joint")
    # # Fit the linear model using the dynamically created formula
    # full_model_tune_all <- lm(formula_str_all, data = df_prs_tune_scaled)
    # 
    # coeff_all <- coef(full_model_tune_all)
    # 
    # # Create formula for selected annotations and PRS_Joint
    # formula_str_selected <- paste("adjusted_pheno_tune", "~",
    #                               ifelse(length(selected_annotations_prop.h[[trait]]) > 0,
    #                                      paste(selected_annotations_prop.h[[trait]], collapse = " + "),
    #                                      "1"))  # Handle case with no selected annotations
    # formula_str_selected <- paste(formula_str_selected, "+ PRS_Joint")
    # 
    # # Fit the linear model using selected annotations
    # selected_model_tune <- lm(formula_str_selected, data = df_prs_tune_scaled)
    # 
    # # Extract coefficients (including intercept)
    # coeff_selected <- coef(selected_model_tune)
    
    # LASSO approach for full model (all annotations + PRS_Joint) alpha=1 corresponds to LASSO, alpha=0 corresponds to Ridge
    predictors_all <- c(ordered_annotations_enrich.p[[trait]], "PRS_Joint")
    formula_str_all <- paste("adjusted_pheno_tune ~ 0 +", paste(predictors_all, collapse = " + "))
    x_all <- model.matrix(as.formula(formula_str_all), data = df_prs_tune_scaled)
    y_all <- adjusted_pheno_tune

    # Fit LASSO with cross-validation
    cv_fit_all <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
    coeff_all <- coef(cv_fit_all, s = "lambda.min")[,1]
    names(coeff_all) <- rownames(coef(cv_fit_all, s = "lambda.min"))

    # LASSO approach for selected model (selected annotations + PRS_Joint)
    selected_predictors <- c(selected_annotations_prop.h[[trait]], "PRS_Joint")
    if (length(selected_predictors) == 0) selected_predictors <- "PRS_Joint"  # Ensure PRS is included

    formula_str_selected <- paste("adjusted_pheno_tune ~ 0 +", paste(selected_predictors, collapse = " + "))
    x_selected <- model.matrix(as.formula(formula_str_selected), data = df_prs_tune_scaled)
    y_selected <- adjusted_pheno_tune

    # Fit LASSO with cross-validation
    cv_fit_selected <- cv.glmnet(x = x_selected, y = y_selected, alpha = 0, standardize = FALSE)
    coeff_selected <- coef(cv_fit_selected, s = "lambda.min")[,1]
    names(coeff_selected) <- rownames(coef(cv_fit_selected, s = "lambda.min"))
    
    # Store coefficients
    trait_coefficients[[trait]][[epoch]] <- list()
    trait_coefficients[[trait]][[epoch]] <- coeff_selected
    
    # test set
    Trait_test_pheno_id <- fread(paste0(tune_dir, "/test_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
    common_test_id <- intersect(Trait_test_pheno_id, Trait_pheno$eid)
    # Subset Covariates and phenos
    covariates_test <- total_covariates[eid %in% common_test_id]
    Trait_pheno_test <- Trait_pheno[eid %in% common_test_id]
    
    # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
    pheno_covariates_test <- merge(Trait_pheno_test, covariates_test, by.x = "eid", by.y = "eid")
    df_prs_test <- all_prs[[trait]][IID %in% common_test_id]
    # Scale all columns except the 'IID' column
    df_prs_test_scaled <- df_prs_test %>%
      mutate(across(-IID, scale))
    null_model_test <- lm(pheno ~ ., data = pheno_covariates_test[, -1])  # exclude eid
    adjusted_pheno_test <- residuals(null_model_test)
    
    # Initialize a list to store linear models.
    linear_models <- list()
    
    # Loop through each PRS and fit the linear model
    for (k in 1:K) {
      formula <- as.formula(paste("adjusted_pheno_test ~", ordered_annotations_enrich.p[[trait]][k]))
      linear_models[[k]] <- lm(formula, data = df_prs_test_scaled)
      
      # Calculate R-squared and store it
      r_square_prs_results[[k]][epoch] <- summary(linear_models[[k]])$r.squared
    }
    # JointPRS model
    linear_models[[K+1]] <- lm(adjusted_pheno_test ~ PRS_Joint, data = df_prs_test_scaled)
    r_square_prs_results[[K+1]][epoch] <- summary(linear_models[[K+1]])$r.squared
    # Calculate residuals of linear combination model with all annotations
    ss_reg_all = coeff_all[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all[-1]
    residual_test_all = adjusted_pheno_test - ss_reg_all
    # Compute R-squared for linear combination
    r_square_prs_results[[K+2]][epoch] <- 1 - var(residual_test_all) / var(adjusted_pheno_test)
    
    # Calculate predicted values with selected annotations and PRS_Joint on test data
    ss_reg_selected <- coeff_selected[1] + 
      as.matrix(df_prs_test_scaled)[, c(selected_annotations_prop.h[[trait]], "PRS_Joint")] %*% 
      coeff_selected[-1]
    
    # Calculate residuals and R-squared
    residual_test_selected <- adjusted_pheno_test - ss_reg_selected
    r_square_prs_results[[K+3]][epoch] <- 1 - var(residual_test_selected) / var(adjusted_pheno_test)
    
    # Store combination results
    combined_ways_result$r_square[combined_ways_index] <- r_square_prs_results[[K+1]][epoch]  # JointPRS
    combined_ways_index <- combined_ways_index + 1
    
    # Calculate OLS result
    formula_str_all_intercept <- as.formula(paste("adjusted_pheno_tune ~ 1 +", paste(predictors_all, collapse = " + ")))
    full_model_tune_all <- lm(formula_str_all_intercept, data = df_prs_tune_scaled)
    coeff_OLS <- coef(full_model_tune_all)
    ss_reg_OLS = coeff_OLS[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_OLS[-1]
    residual_test_OLS = adjusted_pheno_test - ss_reg_OLS
    r_square_prs_OLS <- 1 - var(residual_test_OLS) / var(adjusted_pheno_test)
    
    combined_ways_result$r_square[combined_ways_index] <- r_square_prs_OLS 
    combined_ways_index <- combined_ways_index + 1
    
    ## Calculate predicted values with all annotations and PRS_Joint on test data with LASSO and Ridge approach
    # Fit LASSO with cross-validation
    cv_fit_LASSO <- cv.glmnet(x = x_all, y = y_all, alpha = 1, standardize = FALSE)
    coeff_LASSO <- coef(cv_fit_LASSO, s = "lambda.min")[,1]
    ss_reg_LASSO <- coeff_LASSO[1] + 
      as.matrix(df_prs_test_scaled)[,-1] %*% coeff_LASSO[-1]
    residual_test_LASSO = adjusted_pheno_test - ss_reg_LASSO
    r_square_prs_LASSO <- 1 - var(residual_test_LASSO) / var(adjusted_pheno_test)
    
    combined_ways_result$r_square[combined_ways_index] <- r_square_prs_LASSO 
    combined_ways_index <- combined_ways_index + 1
    
    # Fit Ridge with cross-validation
    cv_fit_Ridge <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
    coeff_Ridge <- coef(cv_fit_Ridge, s = "lambda.min")[,1]
    ss_reg_Ridge <- coeff_Ridge[1] + 
      as.matrix(df_prs_test_scaled)[,-1] %*% coeff_Ridge[-1]
    residual_test_Ridge = adjusted_pheno_test - ss_reg_Ridge
    r_square_prs_Ridge <- 1 - var(residual_test_Ridge) / var(adjusted_pheno_test)
    
    combined_ways_result$r_square[combined_ways_index] <- r_square_prs_Ridge 
    combined_ways_index <- combined_ways_index + 1
  }
  
  
  # Create an empty data frame for R-squared results
  r_square_results <- data.frame(Epoch = 1:4)
  
  
  # Dynamically add columns for each PRS result
  for (k in 1:(K + 3)) {
    prs_label <- if (k <= K) {
      sub("L2_0$", "",ordered_annotations_enrich.p[[trait]][k])
    } else if (k == K + 1) {
      "JointPRS"
    } else if (k == K + 2) {
      "weighted_all"
    } else {
      "weighted_selected"
    }
    r_square_results[[prs_label]] <- r_square_prs_results[[k]]
  }
  
  # Convert the data frame from wide to long format
  r_square_long <- melt(r_square_results, id.vars = "Epoch", variable.name = "PRS_Type", value.name = "R_Squared")
  # Calculate mean r2 across epochs
  r_square_mean <- r_square_long %>%
    group_by(PRS_Type) %>%
    summarise(Mean_R_Squared = mean(R_Squared, na.rm = TRUE)) %>%
    ungroup()
  
  ordered_r_square <- r_square_mean %>% 
    arrange(Mean_R_Squared) %>%
    pull(PRS_Type)
  
  ordered_r_square <- as.character(ordered_r_square)
  ordered_r_square <- ordered_r_square[!ordered_r_square %in% c("JointPRS", "weighted_all","weighted_selected")] # Keep only annotations r2 order
  
  # Perform Rank test with two ranks
  baseline_order <- 1:K
  enrichmentP_order <- match(ordered_r_square, rev(sub("L2_0$", "",ordered_annotations_enrich.p[[trait]])))
  proph2_order <- match(ordered_r_square, sub("L2_0$", "",ordered_annotations_prop.h[[trait]]))
  enrichmentP_rank_test <- cor.test(baseline_order, enrichmentP_order, method = "kendall")
  proph2_rank_test <- cor.test(baseline_order, proph2_order, method = "kendall")
  
  # Store the results
  idx_rank_test <- which(rank_test_result$trait == trait)
  rank_test_result$tau_enrichmentP[idx_rank_test] <- enrichmentP_rank_test$estimate
  rank_test_result$pvalue_enrichmentP[idx_rank_test] <- enrichmentP_rank_test$p.value
  
  rank_test_result$tau_proph2[idx_rank_test] <- proph2_rank_test$estimate
  rank_test_result$pvalue_proph2[idx_rank_test] <- proph2_rank_test$p.value
  # # Plot the results in order of enrich_p
  # p <- ggplot(r_square_long, aes(x = PRS_Type, y = R_Squared)) +
  #   geom_boxplot() +
  #   labs(x = "PRS Type", y = "R2") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #         plot.margin = margin(l = 3.5, unit = "cm"))
  # 
  # # # Print the plot for the current trait
  # # print(p)
  # ggsave(
  #   filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/R2_",trait,"_",pop,"_EnrichmentP.pdf"),
  #   plot = p,
  #   width = 20,
  #   height = 8,
  #   units = "in"
  # )
  
  
  # Reorder the PRS_Type in your plotting dataframe based on this new order
  prs_levels <- c(sub("L2_0$", "",ordered_annotations_prop.h[[trait]]), "JointPRS", "weighted_all", "weighted_selected")

  # Reorder the PRS_Type in your plotting dataframe based on the new combined order
  r_square_long$PRS_Type <- factor(r_square_long$PRS_Type, levels = prs_levels)
  
  
  # Create the boxplot with reordered annotations
  # p <- ggplot(r_square_long, aes(x = PRS_Type, y = R_Squared)) +
  #   geom_boxplot() +
  #   labs(x = "PRS Type", y = "R2") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #         plot.margin = margin(l = 3.5, unit = "cm"))
  # 
  # # Print the plot for the current trait
  # # print(p)
  # ggsave(
  #     filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/R2_",trait,"_",pop,"_Proph2.pdf"),
  #     plot = p,
  #     width = 20,
  #     height = 8,
  #     units = "in"
  #   )
}

rank_test_result

# 2.0 Loop over each trait in the coefficients list
# Convert coefficients to a tidy data frame
coeff_plots <- list()
for (trait in names(trait_coefficients)) {
  
  # Collect all unique coefficient names across all epochs of this trait
  coeffs_name <- c("(Intercept)",selected_annotations_prop.h[[trait]], "PRS_Joint")
  # Prepare data frame for the current trait
  trait_data <- data.frame()
  
  for (epoch in epochs) {
    coeff_values <- trait_coefficients[[trait]][[epoch]]
    
    # Create a data frame for the current epoch with all coefficients
    epoch_df <- data.frame(
      Coefficient = coeffs_name,
      Value = vapply(coeffs_name, function(c) {
        if (c %in% names(coeff_values)) coeff_values[c] else 0
      }, numeric(1)),
      Epoch = epoch,
      Trait = trait,
      stringsAsFactors = FALSE
    )
    
    epoch_df$Coefficient <- factor(epoch_df$Coefficient, levels = coeffs_name)
    trait_data <- rbind(trait_data, epoch_df)
  }
  
  # Ensure Epoch is ordered correctly
  trait_data$Epoch <- factor(trait_data$Epoch, levels = epochs)
  
  # Create the stacked bar plot
  coeff_plots[[trait]] <- ggplot(trait_data, aes(x = Epoch, y = Value, fill = Coefficient)) +
    geom_col(position = "stack") +
    labs(title = paste("Trait:", trait), x = "Epoch", y = "Coefficient Value") +
    theme_minimal()
  
}

# Arrange plots in a 2x2 grid
grid.arrange(
  grobs = coeff_plots,  # Use the list of plots
  nrow = 2,             # 2 rows
  ncol = 2              # 2 columns
)


##### Plot comparison between different combination ways of all annotations + JointPRS
combined_ways_result$combined_way <- factor(combined_ways_result$combined_way, levels = combined_ways)
combined_ways_plot <- ggplot(combined_ways_result, aes(x = trait, y = r_square, fill = combined_way)) +
  geom_boxplot() +
  labs(title = "Boxplot of r_square by Trait with different Combined Ways",
       x = "Traits",
       y = "R-Squared") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability


# Print the plot
print(combined_ways_plot)


##### 2.20 multi-population result#########
pops <- c("EAS", "SAS", "AMR", "AFR")
traits <- c("logTG") 

for (trait in traits){
  # Initialize dataframe for current trait
  combined_ways_result <- data.frame(
    pop = rep(pops, each = length(epochs) * length(combined_ways)),
    epoch = rep(rep(epochs, each = length(combined_ways)), times = length(pops)),
    combined_way = factor(rep(combined_ways, times = length(pops) * length(epochs))),
    r_square = numeric(length(pops) * length(epochs) * length(combined_ways)),  # initialize with zeros
    stringsAsFactors = FALSE  # for other columns if needed
    
  )
  combined_ways_index <- 1
  for (pop in pops){
    ## Set tuning dir
    tune_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/Data/ref_geno_data/", pop, "/4_fold")
    
    prs_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/original_sumstats/trait_", trait, "/tar_", pop)
    
    # Read fine-mapping and PRS files and assign column names
    K <- annotations_info[[trait]]
    
    # Initialize a list to store data frames
    prs_data_list <- list()
    r_square_prs_results <- list()
    # Loop over the number of annotations
    for (k in 1:K) {
      # Construct the file name and path
      file_path <- paste0(prs_dir, "/baselineLD_v2.2_Annot_", k, "/", trait, "_", pop, "_PRScsx_baselineLD_v2.2_Annot_", k, ".sscore")
      
      # Read data using fread and select the required columns
      prs_data <- fread(file_path)[, .(IID = IID, PRS = SCORE1_AVG)]
      
      # Rename PRS column to reflect annotation name.
      setnames(prs_data, "PRS", ordered_annotations_enrich.p[[trait]][k])
      
      # Append to the list
      prs_data_list[[k]] <- prs_data
    }
    
    # Read in the Joint PRS data
    prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop,"/JointPRS/UKB_", trait,"_JointPRS_EUR_", pop,"_prs_", pop, ".sscore"))[, .(IID = IID, PRS_Joint = SCORE1_AVG)]
    
    # Add Joint PRS to the list
    prs_data_list[[K + 1]] <- prs_JointPRS
    
    # Merge all data frames in the list by "IID"
    all_prs[[trait]] <- reduce(prs_data_list, full_join, by = "IID")
    
    # Loop to create a numeric vector for each annotation and one for JointPRS, one for linear combination
    for (k in 1:(K + 3)) {
      r_square_prs_results[[k]] <- numeric(4)  # 4 epochs
    }
    
    Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv"))
    setnames(Trait_pheno, c("eid", "pheno"))
    Trait_pheno[, pheno := scale(pheno)]
    
    trait_coefficients[[trait]] <- list()
    for (epoch in epochs) {
      Trait_tune_pheno_id <- fread(paste0(tune_dir, "/tune_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
      common_tune_id <- intersect(Trait_tune_pheno_id, Trait_pheno$eid)
      # Subset Covariates and phenos
      covariates_tune <- total_covariates[eid %in% common_tune_id]
      Trait_pheno_tune <- Trait_pheno[eid %in% common_tune_id]
      
      # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
      pheno_covariates_tune <- merge(Trait_pheno_tune, covariates_tune, by.x = "eid", by.y = "eid")
      df_prs_tune <- all_prs[[trait]][IID %in% common_tune_id]
      
      null_model_tune <- lm(pheno ~ ., data = pheno_covariates_tune[, -1])  # exclude eid
      adjusted_pheno_tune <- residuals(null_model_tune)
      # Scale all columns except the 'IID' column
      df_prs_tune_scaled <- df_prs_tune %>%
        mutate(across(-IID, scale))
      
      # LASSO approach for full model (all annotations + PRS_Joint) alpha=1 corresponds to LASSO, alpha=0 corresponds to Ridge
      predictors_all <- c(ordered_annotations_enrich.p[[trait]], "PRS_Joint")
      formula_str_all <- paste("adjusted_pheno_tune ~ 0 +", paste(predictors_all, collapse = " + "))
      x_all <- model.matrix(as.formula(formula_str_all), data = df_prs_tune_scaled)
      y_all <- adjusted_pheno_tune
      
      # Fit LASSO with cross-validation
      cv_fit_all <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
      coeff_all <- coef(cv_fit_all, s = "lambda.min")[,1]
      names(coeff_all) <- rownames(coef(cv_fit_all, s = "lambda.min"))
      
      # test set
      Trait_test_pheno_id <- fread(paste0(tune_dir, "/test_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
      common_test_id <- intersect(Trait_test_pheno_id, Trait_pheno$eid)
      # Subset Covariates and phenos
      covariates_test <- total_covariates[eid %in% common_test_id]
      Trait_pheno_test <- Trait_pheno[eid %in% common_test_id]
      
      # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
      pheno_covariates_test <- merge(Trait_pheno_test, covariates_test, by.x = "eid", by.y = "eid")
      df_prs_test <- all_prs[[trait]][IID %in% common_test_id]
      # Scale all columns except the 'IID' column
      df_prs_test_scaled <- df_prs_test %>%
        mutate(across(-IID, scale))
      null_model_test <- lm(pheno ~ ., data = pheno_covariates_test[, -1])  # exclude eid
      adjusted_pheno_test <- residuals(null_model_test)
      
      # Initialize a list to store linear models.
      linear_models <- list()
      
      # JointPRS model
      linear_models[[K+1]] <- lm(adjusted_pheno_test ~ PRS_Joint, data = df_prs_test_scaled)
      r_square_prs_results[[K+1]][epoch] <- summary(linear_models[[K+1]])$r.squared
      # Calculate residuals of linear combination model with all annotations
      ss_reg_all = coeff_all[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all[-1]
      residual_test_all = adjusted_pheno_test - ss_reg_all
      # Compute R-squared for linear combination
      r_square_prs_results[[K+2]][epoch] <- 1 - var(residual_test_all) / var(adjusted_pheno_test)
      
      
      # Store combination results
      combined_ways_result$r_square[combined_ways_index] <- r_square_prs_results[[K+1]][epoch]  # JointPRS
      combined_ways_index <- combined_ways_index + 1
      
      # Calculate OLS result
      formula_str_all_intercept <- as.formula(paste("adjusted_pheno_tune ~ 1 +", paste(predictors_all, collapse = " + ")))
      full_model_tune_all <- lm(formula_str_all_intercept, data = df_prs_tune_scaled)
      coeff_OLS <- coef(full_model_tune_all)
      ss_reg_OLS = coeff_OLS[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_OLS[-1]
      residual_test_OLS = adjusted_pheno_test - ss_reg_OLS
      r_square_prs_OLS <- 1 - var(residual_test_OLS) / var(adjusted_pheno_test)
      
      combined_ways_result$r_square[combined_ways_index] <- r_square_prs_OLS 
      combined_ways_index <- combined_ways_index + 1
      
      ## Calculate predicted values with all annotations and PRS_Joint on test data with LASSO and Ridge approach
      # Fit LASSO with cross-validation
      cv_fit_LASSO <- cv.glmnet(x = x_all, y = y_all, alpha = 1, standardize = FALSE)
      coeff_LASSO <- coef(cv_fit_LASSO, s = "lambda.min")[,1]
      ss_reg_LASSO <- coeff_LASSO[1] + 
        as.matrix(df_prs_test_scaled)[,-1] %*% coeff_LASSO[-1]
      residual_test_LASSO = adjusted_pheno_test - ss_reg_LASSO
      r_square_prs_LASSO <- 1 - var(residual_test_LASSO) / var(adjusted_pheno_test)
      
      combined_ways_result$r_square[combined_ways_index] <- r_square_prs_LASSO 
      combined_ways_index <- combined_ways_index + 1
      
      # Fit Ridge with cross-validation
      cv_fit_Ridge <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
      coeff_Ridge <- coef(cv_fit_Ridge, s = "lambda.min")[,1]
      ss_reg_Ridge <- coeff_Ridge[1] + 
        as.matrix(df_prs_test_scaled)[,-1] %*% coeff_Ridge[-1]
      residual_test_Ridge = adjusted_pheno_test - ss_reg_Ridge
      r_square_prs_Ridge <- 1 - var(residual_test_Ridge) / var(adjusted_pheno_test)
      
      combined_ways_result$r_square[combined_ways_index] <- r_square_prs_Ridge 
      combined_ways_index <- combined_ways_index + 1
    }
  }
  combined_ways_result$combined_way <- factor(combined_ways_result$combined_way, levels = combined_ways)
  combined_ways_plot <- ggplot(combined_ways_result, aes(x = pop, y = r_square, fill = combined_way)) +
    geom_boxplot() +
    labs(title = "Boxplot of r_square of logTG by Pops with different Combined Ways",
         x = "Pops",
         y = "R-Squared") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability
  
  
  # Print the plot
  print(combined_ways_plot)
}


###### 2.28 EUR results when number of combinations differ.#####
for (trait in traits) {
  prs_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/original_sumstats/trait_", trait, "/tar_", pop)
  
  # Read fine-mapping and PRS files and assign column names
  K <- annotations_info[[trait]]
  
  # Arrange the order of annotations in ascending order of total heritability or enrichment_p
  heri_info_dir <- "/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/2.est_enrichment/result/"
  heri_info_raw <- read.table(paste0(heri_info_dir, trait, "_", pop, ".results"), header = T)
  heri_info_raw_binary <- heri_info_raw[heri_info_raw$Category %in% annot_info, ]
  
  ordered_annotations_enrich.p[[trait]] <- heri_info_raw_binary %>% 
    arrange(Enrichment_p) %>%
    pull(Category)
  
  ordered_annotations_prop.h[[trait]] <- heri_info_raw_binary %>% 
    arrange(Prop._h2) %>%
    pull(Category)

  
  # Initialize a list to store data frames
  prs_data_list <- list()
  r_square_prs_results <- list()
  # Loop over the number of annotations
  for (k in 1:K) {
    # Construct the file name and path
    file_path <- paste0(prs_dir, "/baselineLD_v2.2_Annot_", k, "/", trait, "_", pop, "_PRScsx_baselineLD_v2.2_Annot_", k, ".sscore")
    
    # Read data using fread and select the required columns
    prs_data <- fread(file_path)[, .(IID = IID, PRS = SCORE1_AVG)]
    
    # Rename PRS column to reflect annotation name.
    setnames(prs_data, "PRS", ordered_annotations_enrich.p[[trait]][k])
    
    # Append to the list
    prs_data_list[[k]] <- prs_data
  }
  
  # Read in the Joint PRS data
  prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop, "/JointPRS/", trait, "_", pop, "_JointPRS.sscore"))[, .(IID = IID, PRS_Joint = SCORE1_AVG)]
  
  # Add Joint PRS to the list
  prs_data_list[[K + 1]] <- prs_JointPRS
  
  # Merge all data frames in the list by "IID"
  all_prs[[trait]] <- reduce(prs_data_list, full_join, by = "IID")
  
  # Loop to create a numeric vector for JointPRS + aggreated annotations.
  for (k in 1:(K + 1)) {
    r_square_prs_results[[k]] <- numeric(4)  # 4 epochs
  }
  
  Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv"))
  setnames(Trait_pheno, c("eid", "pheno"))
  Trait_pheno[, pheno := scale(pheno)]
  
  for (epoch in epochs) {
    Trait_tune_pheno_id <- fread(paste0(tune_dir, "/tune_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
    common_tune_id <- intersect(Trait_tune_pheno_id, Trait_pheno$eid)
    # Subset Covariates and phenos
    covariates_tune <- total_covariates[eid %in% common_tune_id]
    Trait_pheno_tune <- Trait_pheno[eid %in% common_tune_id]
    
    # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
    pheno_covariates_tune <- merge(Trait_pheno_tune, covariates_tune, by.x = "eid", by.y = "eid")
    df_prs_tune <- all_prs[[trait]][IID %in% common_tune_id]
    
    null_model_tune <- lm(pheno ~ ., data = pheno_covariates_tune[, -1])  # exclude eid
    adjusted_pheno_tune <- residuals(null_model_tune)
    # Scale all columns except the 'IID' column
    df_prs_tune_scaled <- df_prs_tune %>%
      mutate(across(-IID, scale))
    
    
    # # OLS approach
    # # Create the regression formula by collapsing the predictor names into a single string separated by '+'
    # formula_str_all <- paste("adjusted_pheno_tune", "~", paste(ordered_annotations_enrich.p[[trait]], collapse = " + "))
    # # Add the JointPRS result
    # formula_str_all <- paste(formula_str_all, "+ PRS_Joint")
    # # Fit the linear model using the dynamically created formula
    # full_model_tune_all <- lm(formula_str_all, data = df_prs_tune_scaled)
    # 
    # coeff_all <- coef(full_model_tune_all)
    # 
    # # Create formula for selected annotations and PRS_Joint
    # formula_str_selected <- paste("adjusted_pheno_tune", "~",
    #                               ifelse(length(selected_annotations_prop.h[[trait]]) > 0,
    #                                      paste(selected_annotations_prop.h[[trait]], collapse = " + "),
    #                                      "1"))  # Handle case with no selected annotations
    # formula_str_selected <- paste(formula_str_selected, "+ PRS_Joint")
    # 
    # # Fit the linear model using selected annotations
    # selected_model_tune <- lm(formula_str_selected, data = df_prs_tune_scaled)
    # 
    # # Extract coefficients (including intercept)
    # coeff_selected <- coef(selected_model_tune)
    
    # # LASSO approach for full model (all annotations + PRS_Joint) alpha=1 corresponds to LASSO, alpha=0 corresponds to Ridge
    # predictors_all <- c(ordered_annotations_enrich.p[[trait]], "PRS_Joint")
    # formula_str_all <- paste("adjusted_pheno_tune ~ 0 +", paste(predictors_all, collapse = " + "))
    # x_all <- model.matrix(as.formula(formula_str_all), data = df_prs_tune_scaled)
    # y_all <- adjusted_pheno_tune
    
    # # Fit LASSO with cross-validation
    # cv_fit_all <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
    # coeff_all <- coef(cv_fit_all, s = "lambda.min")[,1]
    # names(coeff_all) <- rownames(coef(cv_fit_all, s = "lambda.min"))
    # 
    # # LASSO approach for selected model (selected annotations + PRS_Joint)
    # selected_predictors <- c(selected_annotations_prop.h[[trait]], "PRS_Joint")
    # if (length(selected_predictors) == 0) selected_predictors <- "PRS_Joint"  # Ensure PRS is included
    # 
    # formula_str_selected <- paste("adjusted_pheno_tune ~ 0 +", paste(selected_predictors, collapse = " + "))
    # x_selected <- model.matrix(as.formula(formula_str_selected), data = df_prs_tune_scaled)
    # y_selected <- adjusted_pheno_tune
    # 
    # # Fit LASSO with cross-validation
    # cv_fit_selected <- cv.glmnet(x = x_selected, y = y_selected, alpha = 0, standardize = FALSE)
    # coeff_selected <- coef(cv_fit_selected, s = "lambda.min")[,1]
    # names(coeff_selected) <- rownames(coef(cv_fit_selected, s = "lambda.min"))
    
    # test set
    Trait_test_pheno_id <- fread(paste0(tune_dir, "/test_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
    common_test_id <- intersect(Trait_test_pheno_id, Trait_pheno$eid)
    # Subset Covariates and phenos
    covariates_test <- total_covariates[eid %in% common_test_id]
    Trait_pheno_test <- Trait_pheno[eid %in% common_test_id]
    
    # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
    pheno_covariates_test <- merge(Trait_pheno_test, covariates_test, by.x = "eid", by.y = "eid")
    df_prs_test <- all_prs[[trait]][IID %in% common_test_id]
    # Scale all columns except the 'IID' column
    df_prs_test_scaled <- df_prs_test %>%
      mutate(across(-IID, scale))
    null_model_test <- lm(pheno ~ ., data = pheno_covariates_test[, -1])  # exclude eid
    adjusted_pheno_test <- residuals(null_model_test)
    
    # JointPRS model
    JointPRS_model <- lm(adjusted_pheno_test ~ PRS_Joint, data = df_prs_test_scaled)
    r_square_prs_results[[1]][epoch] <- summary(JointPRS_model)$r.squared
    # Loop through each PRS and fit the linear model
    for (k in 2:(K+1)) {
      # Cumulatively select annotations
      formula_str_selected <- paste("adjusted_pheno_tune", "~", paste(ordered_annotations_prop.h[[trait]][1:(k-1)], collapse = " + "))
      formula_str_selected <- paste(formula_str_selected, "+ PRS_Joint")
      selected_model_tune <- lm(formula_str_selected, data = df_prs_tune_scaled)
      coeff_selected <- coef(selected_model_tune)
      
      
      # LASSO approach for full model (all annotations + PRS_Joint) alpha=1 corresponds to LASSO, alpha=0 corresponds to Ridge
      # formula_str_selected <- paste("adjusted_pheno_tune", "~ 0 + ", paste(ordered_annotations_prop.h[[trait]][1:(k-1)], collapse = " + "))
      # formula_str_selected <- paste(formula_str_selected, "+ PRS_Joint")
      # x_selected <- model.matrix(as.formula(formula_str_selected), data = df_prs_tune_scaled)
      # y_selected <- adjusted_pheno_tune
      # # # Fit LASSO with cross-validation
      # cv_fit_selected <- cv.glmnet(x = x_selected, y = y_selected, alpha = 0, standardize = FALSE)
      # coeff_selected <- coef(cv_fit_selected, s = "lambda.min")[,1]
      # names(coeff_selected) <- rownames(coef(cv_fit_selected, s = "lambda.min"))
      
      # Calculate predicted values with selected annotations and PRS_Joint on test data
      ss_reg_selected <- coeff_selected[1] + 
        as.matrix(df_prs_test_scaled)[, c(ordered_annotations_prop.h[[trait]][1:(k-1)], "PRS_Joint")] %*% 
        coeff_selected[-1]
      
      # Calculate residuals and R-squared
      residual_test_selected <- adjusted_pheno_test - ss_reg_selected
      r_square_prs_results[[k]][epoch] <- 1 - var(residual_test_selected) / var(adjusted_pheno_test)
    }
  }
  
  
  # Create an empty data frame for R-squared results
  r_square_results <- data.frame(Epoch = 1:4)
  
  
  # Dynamically add columns for each PRS result
  for (k in 1:(K + 1)) {
    prs_label <- if (k == 1) {
      "JointPRS"
    } else {
      as.character(k-1)
    }
    r_square_results[[prs_label]] <- r_square_prs_results[[k]]
  }
  
  # Convert the data frame from wide to long format
  r_square_long <- melt(r_square_results, id.vars = "Epoch", variable.name = "PRS_Type", value.name = "R_Squared")
  
  
  # Reorder the PRS_Type in your plotting dataframe based on this new order
  prs_levels <- c("JointPRS", as.character(1:K))
  
  # Reorder the PRS_Type in your plotting dataframe based on the new combined order
  r_square_long$PRS_Type <- factor(r_square_long$PRS_Type, levels = prs_levels)
  
  
  # Create the boxplot with reordered annotations
  p <- ggplot(r_square_long, aes(x = PRS_Type, y = R_Squared)) +
    geom_boxplot() +
    labs(title = paste("R-square Box Plot for", trait, "PRS in 4 fold cross-validation with hm3 panel, LASSO estimate"),
         subtitle = "Cumulative Annotations ordered by Prop._h2",
         x = "PRS Type", y = "R-square") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot for the current trait
  print(p)
}
  
  
