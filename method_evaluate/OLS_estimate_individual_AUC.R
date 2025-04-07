library(purrr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)  # Load the gridExtra package
library(reshape2) # Needed for melting data frames for plotting
library(Matrix)
library(pROC)
library(glmnet)


## EUR Population Results ###

# Step 0: Load PCs and covariates
cov_choice <- c("age_recruit", "sex", paste0("PC", 1:20))
total_covariates <- fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv") %>%
  select(all_of(c("eid", cov_choice)))
# Order covariates
total_covariates <- total_covariates[order(eid)]

## Set tuning dir
tune_dir <- "/gpfs/gibbs/pi/zhao/zb222/Data/ref_geno_data/EUR/4_fold_10K"

# Step 1: Gather necessary data
traits <- c("T2D", "BrC")  # Binary traits
n_traits <- length(traits)
pop <- "EUR"
epochs <- 1:4
annot_dir <- "/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/result/selected_annot/binary/baselineLD_v2.3_hm3_mega" 

# Initialize list to store coefficients
trait_coefficients <- list()
all_prs <- list()

# Initialize combinations results dataframe
combined_ways <- c("JointPRS", "hm3_mega", "JointPRS_annot")
combined_ways_result <- data.frame(
  trait = rep(traits, each = length(epochs)),
  Epochs = rep(epochs, times = n_traits),
  JointPRS = NA_real_,
  hm3_mega = NA_real_,
  JointPRS_annot = NA_real_  # Using underscore instead of hyphen
)


for (trait in traits) {
  prs_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_score/LD_v2.3_panel/trait_", trait, "/tar_", pop)
  
  # Read fine-mapping and PRS files and assign column names
  annotation_info <- readLines(paste0(annot_dir, "/selected_annot_", trait, "_", pop, ".txt"))
  K <- length(annotation_info)
  
  # Initialize a list to store data frames
  prs_data_list <- list()
  AUC_prs_results <- data.frame(Epochs = epochs)
  # Loop over the number of annotations
  for (k in 1:K) {
    # Construct the file name and path
    file_path <- paste0(prs_dir, "/baselineLD_v2.3_Annot_", k, "/", trait, "_", pop, "_JointPRS_baselineLD_v2.3_Annot_", k, ".sscore")
    
    # Read data using fread and select the required columns
    prs_data <- fread(file_path)[, .(IID = IID, PRS = SCORE1_AVG)]
    
    # Rename PRS column to reflect annotation name.
    setnames(prs_data, "PRS", annotation_info[k])
    
    # Append to the list
    prs_data_list[[k]] <- prs_data
  }
  
  # Read in the Joint PRS data
  prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop, "/JointPRS/", trait, "_", pop, "_JointPRS.sscore"))[, .(IID = IID, JointPRS = SCORE1_AVG)]
  
  # Add Joint PRS to the list
  prs_data_list[[K + 1]] <- prs_JointPRS
  
  # Merge all data frames in the list by "IID"
  all_prs[[trait]] <- reduce(prs_data_list, full_join, by = "IID")
  
  # Loop to create a numeric vector for each annotation and one for JointPRS, one for linear combination, one for AIC/BIC
  # for (type in c(annotation_info, "JointPRS", "OLS_all")) {
  #   AUC_prs_results[[type]] <- numeric(4)  # 4 epochs
  # }
  
  Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv"))
  setnames(Trait_pheno, c("eid", "pheno"))
  
  # trait_coefficients[[trait]] <- list()
  
  for (epoch in epochs) {
    Trait_tune_pheno_id <- fread(paste0(tune_dir, "/tune_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
    common_tune_id <- intersect(Trait_tune_pheno_id, Trait_pheno$eid)
    # Subset and phenos
    Trait_pheno_tune <- Trait_pheno[eid %in% common_tune_id][, -1] # match samples and exclude eid column
    Trait_pheno_tune <- as.numeric(unlist(Trait_pheno_tune))
    # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
    df_prs_tune <- all_prs[[trait]][IID %in% common_tune_id]
    
    # Scale all columns except the 'IID' column
    df_prs_tune_scaled <- df_prs_tune %>%
      mutate(across(-IID, scale))
    
    # Logistic Regression
    # Create the regression formula by collapsing the predictor names into a single string separated by '+'
    
    formula_str_all <- paste("Trait_pheno_tune", "~", paste(annotation_info, collapse = " + "))
    # Add the JointPRS result
    formula_str_all <- paste(formula_str_all, "+ JointPRS")
    # Fit the linear model using the dynamically created formula
    full_model_tune_all <- glm(formula_str_all, data = df_prs_tune_scaled, family=binomial(link="logit"))
    
    coeff_all <- coef(full_model_tune_all)
    
    # Logistic Regression with nonnegatice weights
    # formula_str_all <- paste("Trait_pheno_tune ~", 
    #                          paste(annotation_info, collapse = " + "), 
    #                          "+ JointPRS + 0")  # +0 removes intercept
    
    
    # Fit non-negative logistic regression
    full_model_logistic <- glmnet(
      x = df_prs_tune_scaled[,-1],
      y = Trait_pheno_tune,
      family = "binomial",
      lower.limits = 0,    # All coefficients >= 0
      upper.limits = Inf,  # No upper bounds
      lambda = 0,          # No regularization
      intercept = TRUE,    # Unconstrained intercept
      standardize = FALSE  # Keep original scale
    )
    
    # Extract coefficients (intercept first)
    coeff_all_nn_logistic <- as.numeric(coef(full_model_logistic, s = 0))
    names(coeff_all_nn_logistic) <- c("(Intercept)", colnames(df_prs_tune_scaled[,-1]))
    
    # test set
    Trait_test_pheno_id <- fread(paste0(tune_dir, "/test_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
    common_test_id <- intersect(Trait_test_pheno_id, Trait_pheno$eid)
    # Subset Covariates and phenos
    Trait_pheno_test <- Trait_pheno[eid %in% common_test_id][, -1]
    Trait_pheno_test <- as.numeric(unlist(Trait_pheno_test))
    # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
    df_prs_test <- all_prs[[trait]][IID %in% common_test_id]
    # Scale all columns except the 'IID' column
    df_prs_test_scaled <- df_prs_test %>%
      mutate(across(-IID, scale))

    
    # Initialize a list to store linear models.
    linear_models <- list()
    
    # Loop through each PRS and fit the linear model
    for (k in 1:K) {
      formula <- as.formula(paste("Trait_pheno_test ~", annotation_info[k]))
      linear_models[[k]] <- glm(formula, data=df_prs_test_scaled,family=binomial(link="logit"))
      glmfit_prob = predict(linear_models[[k]], type="response")
      glmfit_auc = roc(Trait_pheno_test, glmfit_prob, quiet=T, plot=F)$auc
      
      # # Calculate R-squared and store it
      AUC_prs_results[[annotation_info[k]]][epoch] <- glmfit_auc
      
      # Equivalently,
      # AUC_prs_results[[annotation_info[k]]][epoch] <- cor(adjusted_pheno_test, df_prs_test_scaled[[annotation_info[k]]])^2
    }
    # JointPRS model
    linear_models[[K+1]] <- glm(Trait_pheno_test ~ JointPRS, data = df_prs_test_scaled, family=binomial(link="logit"))
    glmfit_prob = predict(linear_models[[k+1]], type="response")
    glmfit_auc = roc(Trait_pheno_test, glmfit_prob, quiet=T, plot=F)$auc
    AUC_prs_results[["JointPRS"]][epoch] <- glmfit_auc
    # Calculate residuals of linear combination model with all annotations
    log_odds = coeff_all[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all[-1]
    manual_probs <- 1 / (1 + exp(-log_odds))
    # Compute AUC for combination
    AUC_prs_results[["Logistic_all"]][epoch] <- roc(Trait_pheno_test, as.numeric(manual_probs), quiet=T, plot=F)$auc
    
    
    # Calculate non-negative logistic regression predictions
    log_odds_nn_logistic <- coeff_all_nn_logistic[1] + as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all_nn_logistic[-1]
    manual_probs_nn_logistic <- 1 / (1 + exp(-log_odds_nn_logistic))
    AUC_prs_results[["NNLogistic_all"]][epoch] <- roc(Trait_pheno_test, as.numeric(manual_probs_nn_logistic), quiet=T, plot=F)$auc
    
    # # Calculate predicted values with selected annotations and PRS_Joint on test data
    # ss_reg_selected <- coeff_selected[1] + 
    #   as.matrix(df_prs_test_scaled)[, c(selected_annotations_prop.h[[trait]], "PRS_Joint")] %*% 
    #   coeff_selected[-1]
    # 
    # # Calculate residuals and R-squared
    # residual_test_selected <- adjusted_pheno_test - ss_reg_selected
    # AUC_prs_results[[K+3]][epoch] <- 1 - var(residual_test_selected) / var(adjusted_pheno_test)
    
    # Store combination results
    combined_ways_result_idx <- which(combined_ways_result$trait == trait & combined_ways_result$Epochs == epoch)
    combined_ways_result$JointPRS[combined_ways_result_idx] <- AUC_prs_results[["JointPRS"]][epoch]
    combined_ways_result$hm3_mega[combined_ways_result_idx] <- AUC_prs_results[["hm3_mega"]][epoch]
    combined_ways_result$JointPRS_annot[combined_ways_result_idx] <- AUC_prs_results[["NNLogistic_all"]][epoch]
  }
  
  
  # Create an empty data frame for R-squared results
  
  
  # Convert the data frame from wide to long format
  AUC_long <- melt(AUC_prs_results, id.vars = "Epochs", variable.name = "PRS_Type", value.name = "AUC")
  
  # # Reorder the PRS_Type in your plotting dataframe based on this new order
  # prs_levels <- c(ordered_annotations_prop.h[[trait]], "JointPRS", "weighted_all", "weighted_selected")
  
  # # Reorder the PRS_Type in your plotting dataframe based on the new combined order
  # AUC_long$PRS_Type <- factor(AUC_long$PRS_Type, levels = prs_levels)
  
  
  # # Create the boxplot with reordered annotations
  # p <- ggplot(AUC_long, aes(x = PRS_Type, y = AUC)) +
  #   geom_boxplot() +
  #   labs(title = paste("AUC Box Plot for Trait:", trait, " Pop: ",pop, "PRS in 4 fold cross-validation"),
  #        subtitle = "Annotations ordered by desc(Prop._h2)",
  #        x = "PRS Type", y = "AUC") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # # Print the plot for the current trait
  # print(p)
  
  # p <- ggplot(AUC_long, aes(x = PRS_Type, y = AUC)) +
  #   geom_boxplot() +
  #   labs(x = "PRS Type", y = "AUC") +
  #   theme_minimal() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1)  # Angled labels
  #     # plot.margin = margin(l = 2, unit = "cm")  # Left margin = 2cm (adjust as needed)
  #   )
  # 
  # 
  # # Save the results as pdf
  # n_categories <- length(unique(r_square_long$PRS_Type))
  # max_label_length <- max(nchar(as.character(r_square_long$PRS_Type)))
  # 
  # print(p)
  # 
  # adaptive_width <- 2 + n_categories * 0.8+ max_label_length * 0.1
  # adaptive_height <- 7 + (n_categories * 0.1)  # Taller if many categories
  # 
  # ggsave(
  #   filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/AUC_",trait,"_",pop,".pdf"),
  #   plot = p,
  #   width = adaptive_width,
  #   height = adaptive_height,
  #   units = "in"
  # )
  
}


##### Plot comparison between different combination ways of all annotations + JointPRS
##### Plot comparison between different combination ways of JointPRS, hm3, mega
combined_ways_result_long <- combined_ways_result %>%
  pivot_longer(
    cols = c("JointPRS", "hm3_mega", "JointPRS_annot"),
    names_to = "PRS_type",
    values_to = "R2"
  )
combined_ways_result_long$PRS_type <- factor(combined_ways_result_long$PRS_type, levels = combined_ways)
combined_ways_result_long$trait <- factor(combined_ways_result_long$trait, levels = traits)
combined_ways_plot <- ggplot(combined_ways_result_long, aes(x = trait, y = R2, fill = PRS_type)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7) +
  labs(
    x = "Trait",
    y = "R2",
    fill = "PRS Type"
  ) +
  scale_fill_manual(
    values = c("JointPRS" = "#1f77b4", 
               "hm3_mega" = "#ff7f0e", 
               "JointPRS_annot" = "#2ca02c"),
    labels = c("JointPRS", "hm3_mega", "JointPRS_annot")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Print the plot
print(combined_ways_plot)
ggsave(
  filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/AUC_EUR_JointPRS_annot.pdf"),
  plot = combined_ways_plot,
  width = 8,
  height = 5,
  units = "in"
)


##### Multi-Population Results #########

# Step 0: Load PCs and covariates
cov_choice <- c("age_recruit", "sex", paste0("PC", 1:20))
total_covariates <- fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/cov_data/adjust.csv") %>%
  select(all_of(c("eid", cov_choice)))
# Order covariates
total_covariates <- total_covariates[order(eid)]

## Set tuning dir
tune_dir <- "/gpfs/gibbs/pi/zhao/zb222/Data/ref_geno_data/EUR/4_fold"

# Step 1: Gather necessary data
traits <- c("T2D", "BrC")  # Binary traits
pops <- c("EAS", "AFR")
epochs <- 1:4
annot_dir <- "/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/result/selected_annot/binary/baselineLD_v2.3_hm3_mega" 

# Initialize list to store coefficients
trait_coefficients <- list()
all_prs <- list()

# Initialize combinations results dataframe
combined_ways <- c("JointPRS", "Logistic")
combined_ways_result <- data.frame(
  Epochs = epochs
)

for (trait in traits){
  for (pop in pops){
    ## Set tuning dir
    tune_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/Data/ref_geno_data/", pop, "/4_fold")
    
    prs_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_score/LD_v2.3_panel/trait_", trait, "/tar_", pop)
    # Read fine-mapping and PRS files and assign column names
    annotation_info <- readLines(paste0(annot_dir, "/selected_annot_", trait, "_", pop, ".txt"))
    K <- length(annotation_info)
    
    # Initialize a list to store data frames
    prs_data_list <- list()
    AUC_prs_results <- data.frame(Epochs = epochs)
    # Loop over the number of annotations
    for (k in 1:K) {
      # Construct the file name and path
      file_path <- paste0(prs_dir, "/baselineLD_v2.3_Annot_", k, "/", trait, "_", pop, "_JointPRS_baselineLD_v2.3_Annot_", k, ".sscore")
      
      # Read data using fread and select the required columns
      prs_data <- fread(file_path)[, .(IID = IID, PRS = SCORE1_AVG)]
      
      # Rename PRS column to reflect annotation name.
      setnames(prs_data, "PRS", annotation_info[k])
      
      # Append to the list
      prs_data_list[[k]] <- prs_data
    }
    
    # Read in the Joint PRS data
    prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop, "/JointPRS/", trait, "_", pop, "_JointPRS.sscore"))[, .(IID = IID, JointPRS = SCORE1_AVG)]
    
    # Add Joint PRS to the list
    prs_data_list[[K + 1]] <- prs_JointPRS
    
    # Merge all data frames in the list by "IID"
    all_prs[[trait]] <- reduce(prs_data_list, full_join, by = "IID")
    
    # Loop to create a numeric vector for each annotation and one for JointPRS, one for linear combination, one for AIC/BIC
    # for (type in c(annotation_info, "JointPRS", "OLS_all")) {
    #   AUC_prs_results[[type]] <- numeric(4)  # 4 epochs
    # }
    
    Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv"))
    setnames(Trait_pheno, c("eid", "pheno"))
    
    # trait_coefficients[[trait]] <- list()
    
    for (epoch in epochs) {
      Trait_tune_pheno_id <- fread(paste0(tune_dir, "/tune_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
      common_tune_id <- intersect(Trait_tune_pheno_id, Trait_pheno$eid)
      # Subset and phenos
      Trait_pheno_tune <- Trait_pheno[eid %in% common_tune_id][, -1] # match samples and exclude eid column
      Trait_pheno_tune <- as.numeric(unlist(Trait_pheno_tune))
      # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
      df_prs_tune <- all_prs[[trait]][IID %in% common_tune_id]
      
      # Scale all columns except the 'IID' column
      df_prs_tune_scaled <- df_prs_tune %>%
        mutate(across(-IID, scale))
      
      # Logistic Regression
      # Create the regression formula by collapsing the predictor names into a single string separated by '+'
      
      formula_str_all <- paste("Trait_pheno_tune", "~", paste(annotation_info, collapse = " + "))
      # Add the JointPRS result
      formula_str_all <- paste(formula_str_all, "+ JointPRS")
      # Fit the linear model using the dynamically created formula
      full_model_tune_all <- glm(formula_str_all, data = df_prs_tune_scaled, family=binomial(link="logit"))
      
      coeff_all <- coef(full_model_tune_all)
      
      # Fit non-negative logistic regression
      full_model_logistic <- glmnet(
        x = df_prs_tune_scaled[,-1],
        y = Trait_pheno_tune,
        family = "binomial",
        lower.limits = 0,    # All coefficients >= 0
        upper.limits = Inf,  # No upper bounds
        lambda = 0,          # No regularization
        intercept = TRUE,    # Unconstrained intercept
        standardize = FALSE  # Keep original scale
      )
      
      # Extract coefficients (intercept first)
      coeff_all_nn_logistic <- as.numeric(coef(full_model_logistic, s = 0))
      names(coeff_all_nn_logistic) <- c("(Intercept)", colnames(df_prs_tune_scaled[,-1]))
      
      # test set
      Trait_test_pheno_id <- fread(paste0(tune_dir, "/test_set_", epoch, ".ids"))[, V1]  # id from tune set, 8k
      common_test_id <- intersect(Trait_test_pheno_id, Trait_pheno$eid)
      # Subset Covariates and phenos
      Trait_pheno_test <- Trait_pheno[eid %in% common_test_id][, -1]
      Trait_pheno_test <- as.numeric(unlist(Trait_pheno_test))
      # Merge the data tables on 'eid' from Trait_pheno and 'IID' from covariates
      df_prs_test <- all_prs[[trait]][IID %in% common_test_id]
      # Scale all columns except the 'IID' column
      df_prs_test_scaled <- df_prs_test %>%
        mutate(across(-IID, scale))
      
      
      # Initialize a list to store linear models.
      linear_models <- list()
      
      # Loop through each PRS and fit the linear model
      for (k in 1:K) {
        formula <- as.formula(paste("Trait_pheno_test ~", annotation_info[k]))
        linear_models[[k]] <- glm(formula, data=df_prs_test_scaled,family=binomial(link="logit"))
        glmfit_prob = predict(linear_models[[k]], type="response")
        glmfit_auc = roc(Trait_pheno_test, glmfit_prob, quiet=T, plot=F)$auc
        
        # # Calculate R-squared and store it
        AUC_prs_results[[annotation_info[k]]][epoch] <- glmfit_auc
        
        # Equivalently,
        # AUC_prs_results[[annotation_info[k]]][epoch] <- cor(adjusted_pheno_test, df_prs_test_scaled[[annotation_info[k]]])^2
      }
      # JointPRS model
      linear_models[[K+1]] <- glm(Trait_pheno_test ~ JointPRS, data = df_prs_test_scaled, family=binomial(link="logit"))
      glmfit_prob = predict(linear_models[[k+1]], type="response")
      glmfit_auc = roc(Trait_pheno_test, glmfit_prob, quiet=T, plot=F)$auc
      AUC_prs_results[["JointPRS"]][epoch] <- glmfit_auc
      # Calculate residuals of linear combination model with all annotations
      log_odds = coeff_all[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all[-1]
      manual_probs <- 1 / (1 + exp(-log_odds))
      # Compute AUC for combination
      AUC_prs_results[["Logistic_all"]][epoch] <- roc(Trait_pheno_test, as.numeric(manual_probs), quiet=T, plot=F)$auc
      
      # Calculate non-negative logistic regression predictions
      log_odds_nn_logistic <- coeff_all_nn_logistic[1] + as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all_nn_logistic[-1]
      manual_probs_nn_logistic <- 1 / (1 + exp(-log_odds_nn_logistic))
      AUC_prs_results[["NNLogistic_all"]][epoch] <- roc(Trait_pheno_test, as.numeric(manual_probs_nn_logistic), quiet=T, plot=F)$auc
      
      # # Calculate predicted values with selected annotations and PRS_Joint on test data
      # ss_reg_selected <- coeff_selected[1] + 
      #   as.matrix(df_prs_test_scaled)[, c(selected_annotations_prop.h[[trait]], "PRS_Joint")] %*% 
      #   coeff_selected[-1]
      # 
      # # Calculate residuals and R-squared
      # residual_test_selected <- adjusted_pheno_test - ss_reg_selected
      # AUC_prs_results[[K+3]][epoch] <- 1 - var(residual_test_selected) / var(adjusted_pheno_test)
      
    
      
      # Store combination results
      combined_ways_result_idx <- which(combined_ways_result$trait == trait & combined_ways_result$Epochs == epoch)
      combined_ways_result$JointPRS[combined_ways_result_idx] <- AUC_prs_results[["JointPRS"]][epoch]
      combined_ways_result$hm3_mega[combined_ways_result_idx] <- AUC_prs_results[["hm3_mega"]][epoch]
      combined_ways_result$JointPRS_annot[combined_ways_result_idx] <- AUC_prs_results[["NNLogistic_all"]][epoch]
    }
    
    
    # Create an empty data frame for R-squared results
    
    
    # Convert the data frame from wide to long format
    AUC_long <- melt(AUC_prs_results, id.vars = "Epochs", variable.name = "PRS_Type", value.name = "AUC")
    
    # # Reorder the PRS_Type in your plotting dataframe based on this new order
    # prs_levels <- c(ordered_annotations_prop.h[[trait]], "JointPRS", "weighted_all", "weighted_selected")
    
    # # Reorder the PRS_Type in your plotting dataframe based on the new combined order
    # AUC_long$PRS_Type <- factor(AUC_long$PRS_Type, levels = prs_levels)
    
    
    # Create the boxplot with reordered annotations
    # p <- ggplot(AUC_long, aes(x = PRS_Type, y = AUC)) +
    #   geom_boxplot() +
    #   labs(title = paste("AUC Box Plot for Trait:", trait, " Pop: ",pop, "PRS in 4 fold cross-validation"),
    #        subtitle = "Annotations ordered by desc(Prop._h2)",
    #        x = "PRS Type", y = "AUC") +
    #   theme_minimal() +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # p <- ggplot(AUC_long, aes(x = PRS_Type, y = AUC)) +
    #   geom_boxplot() +
    #   labs(x = "PRS Type", y = "AUC") +
    #   theme_minimal() +
    #   theme(
    #     axis.text.x = element_text(angle = 45, hjust = 1),  # Angled labels
    #     plot.margin = margin(l = 2, unit = "cm")  # Left margin = 2cm (adjust as needed)
    #   )
    # 
    # 
    # # Save the results as pdf
    # n_categories <- length(unique(r_square_long$PRS_Type))
    # max_label_length <- max(nchar(as.character(r_square_long$PRS_Type)))
    # 
    # print(p)
    # 
    # adaptive_width <- 2 + n_categories * 0.8+ max_label_length * 0.1
    # adaptive_height <- 7 + (n_categories * 0.1)  # Taller if many categories
    # 
    # ggsave(
    #   filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/AUC_",trait,"_",pop,".pdf"),
    #   plot = p,
    #   width = adaptive_width,
    #   height = adaptive_height,
    #   units = "in"
    # )
  }
}