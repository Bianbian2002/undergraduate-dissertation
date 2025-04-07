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
if (!require("nnls", quietly = TRUE)) {
  install.packages("nnls")
}
library(nnls) # for non-negative OLS

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
traits <- c("HDL","LDL", "TC", "logTG", "Height", "BMI")  # Quantitative traits like "Height", "BMI", "HDL", "LDL", 'TC", "logTG"
pop <- "EUR"
n_traits <- length(traits)
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
  r_square_prs_results <- data.frame(Epochs = epochs)
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
  #   r_square_prs_results[[type]] <- numeric(4)  # 4 epochs
  # }
  
  Trait_pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/pheno_data/", trait, "/", trait, "_", pop, ".tsv"))
  setnames(Trait_pheno, c("eid", "pheno"))
  Trait_pheno[, pheno := scale(pheno)]
  
  # trait_coefficients[[trait]] <- list()
  
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
    
    
    # OLS approach
    # Create the regression formula by collapsing the predictor names into a single string separated by '+'
    formula_str_all <- paste("adjusted_pheno_tune", "~", paste(annotation_info, collapse = " + "))
    # Add the JointPRS result
    formula_str_all <- paste(formula_str_all, "+ JointPRS")
    # Fit the linear model using the dynamically created formula
    full_model_tune_all <- lm(formula_str_all, data = df_prs_tune_scaled)

    coeff_all <- coef(full_model_tune_all)
    
    # Non negative OLS
    full_model_tune_all_nnls <- nnls(as.matrix(df_prs_tune_scaled)[, -1], adjusted_pheno_tune)
    coeff_all_nnls <- coef(full_model_tune_all_nnls)
    
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
    # 
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
    # 
    # # Store coefficients
    # trait_coefficients[[trait]][[epoch]] <- list()
    # trait_coefficients[[trait]][[epoch]] <- coeff_selected
    
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
      formula <- as.formula(paste("adjusted_pheno_test ~", annotation_info[k]))
      linear_models[[k]] <- lm(formula, data = df_prs_test_scaled)

      # # Calculate R-squared and store it
      r_square_prs_results[[annotation_info[k]]][epoch] <- summary(linear_models[[k]])$r.squared
      
      # Equivalently,
      # r_square_prs_results[[annotation_info[k]]][epoch] <- cor(adjusted_pheno_test, df_prs_test_scaled[[annotation_info[k]]])^2
    }
    # JointPRS model
    linear_models[[K+1]] <- lm(adjusted_pheno_test ~ JointPRS, data = df_prs_test_scaled)
    r_square_prs_results[["JointPRS"]][epoch] <- summary(linear_models[[K+1]])$r.squared
    # Calculate residuals of linear combination model with all annotations
    ss_reg_all = coeff_all[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all[-1]
    residual_test_all = adjusted_pheno_test - ss_reg_all
    # Compute R-squared for linear combination
    r_square_prs_results[["OLS_all"]][epoch] <- 1 - var(residual_test_all) / var(adjusted_pheno_test)
    
    # Calculate nnls
    ss_reg_all_nnls = as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all_nnls
    residual_test_all_nnls = adjusted_pheno_test - ss_reg_all_nnls
    # Compute R-squared for linear combination
    r_square_prs_results[["NNOLS_all"]][epoch] <- 1 - var(residual_test_all_nnls) / var(adjusted_pheno_test)
    
    # # Calculate predicted values with selected annotations and PRS_Joint on test data
    # ss_reg_selected <- coeff_selected[1] + 
    #   as.matrix(df_prs_test_scaled)[, c(selected_annotations_prop.h[[trait]], "PRS_Joint")] %*% 
    #   coeff_selected[-1]
    # 
    # # Calculate residuals and R-squared
    # residual_test_selected <- adjusted_pheno_test - ss_reg_selected
    # r_square_prs_results[[K+3]][epoch] <- 1 - var(residual_test_selected) / var(adjusted_pheno_test)
    
    
    # Store combination results
    combined_ways_result_idx <- which(combined_ways_result$trait == trait & combined_ways_result$Epochs == epoch)
    combined_ways_result$JointPRS[combined_ways_result_idx] <- r_square_prs_results[["JointPRS"]][epoch]
    combined_ways_result$hm3_mega[combined_ways_result_idx] <- r_square_prs_results[["hm3_mega"]][epoch]
    combined_ways_result$JointPRS_annot[combined_ways_result_idx] <- r_square_prs_results[["NNOLS_all"]][epoch]
  }

  
  # Convert the data frame from wide to long format
  r_square_long <- melt(r_square_prs_results, id.vars = "Epochs", variable.name = "PRS_Type", value.name = "R_Squared")
  
  # # Reorder the PRS_Type in your plotting dataframe based on this new order
  # prs_levels <- c(ordered_annotations_prop.h[[trait]], "JointPRS", "weighted_all", "weighted_selected")
  
  # # Reorder the PRS_Type in your plotting dataframe based on the new combined order
  # r_square_long$PRS_Type <- factor(r_square_long$PRS_Type, levels = prs_levels)
  
  
  # # Create the boxplot with reordered annotations
  # p <- ggplot(r_square_long, aes(x = PRS_Type, y = R_Squared)) +
  #   geom_boxplot() +
  #   labs(title = paste("R-square Box Plot for Trait:", trait, " Pop: ",pop, "PRS in 4 fold cross-validation"),
  #        subtitle = "Annotations ordered by desc(Prop._h2)",
  #        x = "PRS Type", y = "R-square") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # # Print the plot for the current trait
  # print(p)

  # p <- ggplot(r_square_long, aes(x = PRS_Type, y = R_Squared)) +
  #   geom_boxplot() +
  #   labs(x = "PRS Type", y = "R-square") +
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
  #   filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/R2_",trait,"_",pop,".pdf"),
  #   plot = p,
  #   width = adaptive_width,
  #   height = adaptive_height,
  #   units = "in"
  # )
  
}


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
    filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/R2_EUR_JointPRS_annot.pdf"),
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


# Step 1: Gather necessary data
epochs <- 1:4
annot_dir <- "/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/result/selected_annot/binary/baselineLD_v2.3_hm3_mega" 

# Initialize list to store coefficients
trait_coefficients <- list()
all_prs <- list()

pops_group1 <- c("EAS","SAS","AFR","AMR")
pops_group2 <- c("EAS", "AFR")
n_pops_1 <- length(pops_group1)
n_pops_2 <- length(pops_group2)

traits_group1 <- c("HDL", "LDL", "logTG", "TC")  # "Height", "BMI", "T2D", "BrC"
traits_group2 <- c("Height", "BMI")
n_traits_1 <- length(traits_group1)
n_traits_2 <- length(traits_group2)

combined_ways <- c("JointPRS", "hm3_mega", "JointPRS_annot", "JointPRS_annot_pen_LASSO", "JointPRS_annot_pen_Ridge")
combined_ways_result_group1 <- data.frame(
  pop = rep(pops_group1, each = length(epochs) * n_traits_1),
  trait = rep(rep(traits_group1, each = length(epochs)), times = n_pops_1),
  Epochs = rep(epochs, times = n_traits_1 * n_pops_1),
  JointPRS = NA_real_,
  hm3_mega = NA_real_,
  JointPRS_annot = NA_real_,
  JointPRS_annot_pen_LASSO = NA_real_,
  JointPRS_annot_pen_Ridge = NA_real_
)
combined_ways_result_group2 <- data.frame(
  pop = rep(pops_group2, each = length(epochs) * n_traits_2),
  trait = rep(rep(traits_group2, each = length(epochs)), times = n_pops_2),
  Epochs = rep(epochs, times = n_traits_2 * n_pops_2),
  JointPRS = NA_real_,
  hm3_mega = NA_real_,
  JointPRS_annot = NA_real_,
  JointPRS_annot_pen_LASSO = NA_real_,
  JointPRS_annot_pen_Ridge = NA_real_
)
combined_ways_result <- rbind(combined_ways_result_group1, combined_ways_result_group2)

for (group_num in c(1,2)){
  for (trait in get(paste0("traits_group",group_num))){
    for (pop in get(paste0("pops_group",group_num))){
      ## Set tuning dir
      tune_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/Data/ref_geno_data/", pop, "/4_fold")
      
      prs_dir <- paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_score/LD_v2.3_panel/trait_", trait, "/tar_", pop)
      
      # Read fine-mapping and PRS files and assign column names
      annotation_info <- readLines(paste0(annot_dir, "/selected_annot_", trait, "_", pop, ".txt"))
      K <- length(annotation_info)
      
      # Initialize a list to store data frames
      prs_data_list <- list()
      r_square_prs_results <- data.frame(Epochs = epochs)
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
      if (trait %in% traits_group2){ 
        prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop, "/JointPRS/", trait, "_", pop, "_JointPRS.sscore"))[, .(IID = IID, JointPRS = SCORE1_AVG)]
      }else{
        prs_JointPRS <- fread(paste0("/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_scores_normal/hm3/trait_", trait, "/tar_", pop, "/JointPRS/UKB_",trait,"_JointPRS_EUR_",pop,"_prs_",pop,".sscore"))[, .(IID = IID, JointPRS = SCORE1_AVG)]
      }
      # some format as JointPRS/UKB_",trait,"_JointPRS_EUR_",pop,"_prs_",pop,".sscore" or /JointPRS/", trait, "_", pop, "_JointPRS.sscore"
      
      # Add Joint PRS to the list
      prs_data_list[[K + 1]] <- prs_JointPRS
      
      # Merge all data frames in the list by "IID"
      all_prs[[trait]] <- reduce(prs_data_list, full_join, by = "IID")
      
      # Loop to create a numeric vector for each annotation and one for JointPRS, one for linear combination, one for AIC/BIC
      # for (type in c(annotation_info, "JointPRS", "OLS_all")) {
      #   r_square_prs_results[[type]] <- numeric(4)  # 4 epochs
      # }
      
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
        
        # # LASSO approach for full model (all annotations + PRS_Joint) alpha=1 corresponds to LASSO, alpha=0 corresponds to Ridge
        formula_str_all <- paste("adjusted_pheno_tune", "~ 0 + ", paste(annotation_info, collapse = " + "))
        # Add the JointPRS result
        formula_str_all <- paste(formula_str_all, "+ JointPRS")
        
        x_all <- model.matrix(as.formula(formula_str_all), data = df_prs_tune_scaled)
        y_all <- adjusted_pheno_tune
  
        # Fit LASSO with cross-validation
        cv_fit_all <- cv.glmnet(x = x_all, y = y_all, alpha = 1, standardize = FALSE)
        coeff_all_LASSO <- coef(cv_fit_all, s = "lambda.min")[,1]
        names(coeff_all_LASSO) <- rownames(coef(cv_fit_all, s = "lambda.min"))
        
        # # Ridge approach for full model (all annotations + PRS_Joint) alpha=1 corresponds to LASSO, alpha=0 corresponds to Ridge
        
        # Fit Ridge with cross-validation
        cv_fit_all <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
        coeff_all_Ridge <- coef(cv_fit_all, s = "lambda.min")[,1]
        names(coeff_all_Ridge) <- rownames(coef(cv_fit_all, s = "lambda.min"))
        
        # OLS approach
        # Create the regression formula by collapsing the predictor names into a single string separated by '+'
        formula_str_all <- paste("adjusted_pheno_tune", "~", paste(annotation_info, collapse = " + "))
        # Add the JointPRS result
        formula_str_all <- paste(formula_str_all, "+ JointPRS")
        # Fit the linear model using the dynamically created formula
        full_model_tune_all <- lm(formula_str_all, data = df_prs_tune_scaled)
        coeff_all <- coef(full_model_tune_all)
        
        # Non negative OLS
        full_model_tune_all_nnls <- nnls(as.matrix(df_prs_tune_scaled)[, -1], adjusted_pheno_tune)
        coeff_all_nnls <- coef(full_model_tune_all_nnls)
        
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
        
        for (k in 1:K) {
          formula <- as.formula(paste("adjusted_pheno_test ~", annotation_info[k]))
          linear_models[[k]] <- lm(formula, data = df_prs_test_scaled)
          
          # # Calculate R-squared and store it
          r_square_prs_results[[annotation_info[k]]][epoch] <- summary(linear_models[[k]])$r.squared
          
          # Equivalently,
          # r_square_prs_results[[annotation_info[k]]][epoch] <- cor(adjusted_pheno_test, df_prs_test_scaled[[annotation_info[k]]])^2
        }
        # JointPRS model
        linear_models[[K+1]] <- lm(adjusted_pheno_test ~ JointPRS, data = df_prs_test_scaled)
        r_square_prs_results[["JointPRS"]][epoch] <- summary(linear_models[[K+1]])$r.squared
        # Calculate residuals of linear combination model with all annotations
        ss_reg_all = coeff_all[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all[-1]
        residual_test_all = adjusted_pheno_test - ss_reg_all
        # Compute R-squared for linear combination
        r_square_prs_results[["OLS_all"]][epoch] <- 1 - var(residual_test_all) / var(adjusted_pheno_test)
        
        # Calculate nnls
        ss_reg_all_nnls = as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all_nnls
        residual_test_all_nnls = adjusted_pheno_test - ss_reg_all_nnls
        # Compute R-squared for linear combination
        r_square_prs_results[["NNOLS_all"]][epoch] <- 1 - var(residual_test_all_nnls) / var(adjusted_pheno_test)
        
        ## Calculate LASSO
        ss_reg_all_LASSO = coeff_all_LASSO[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all_LASSO[-1]
        residual_test_all_LASSO = adjusted_pheno_test - ss_reg_all_LASSO
        # Compute R-squared for linear combination
        r_square_prs_results[["LASSO_all"]][epoch] <- 1 - var(residual_test_all_LASSO) / var(adjusted_pheno_test)
        
        ## Calculate Ridge
        ss_reg_all_Ridge = coeff_all_Ridge[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_all_Ridge[-1]
        residual_test_all_Ridge = adjusted_pheno_test - ss_reg_all_Ridge
        # Compute R-squared for linear combination
        r_square_prs_results[["Ridge_all"]][epoch] <- 1 - var(residual_test_all_Ridge) / var(adjusted_pheno_test)
        
        
        # 
        # # Calculate OLS result
        # formula_str_all_intercept <- as.formula(paste("adjusted_pheno_tune ~ 1 +", paste(predictors_all, collapse = " + ")))
        # full_model_tune_all <- lm(formula_str_all_intercept, data = df_prs_tune_scaled)
        # coeff_OLS <- coef(full_model_tune_all)
        # ss_reg_OLS = coeff_OLS[1] +  as.matrix(df_prs_test_scaled)[,-1] %*% coeff_OLS[-1]
        # residual_test_OLS = adjusted_pheno_test - ss_reg_OLS
        # r_square_prs_OLS <- 1 - var(residual_test_OLS) / var(adjusted_pheno_test)
        # 
        # combined_ways_result$r_square[combined_ways_index] <- r_square_prs_OLS 
        # combined_ways_index <- combined_ways_index + 1
        # 
        # ## Calculate predicted values with all annotations and PRS_Joint on test data with LASSO and Ridge approach
        # # Fit LASSO with cross-validation
        # cv_fit_LASSO <- cv.glmnet(x = x_all, y = y_all, alpha = 1, standardize = FALSE)
        # coeff_LASSO <- coef(cv_fit_LASSO, s = "lambda.min")[,1]
        # ss_reg_LASSO <- coeff_LASSO[1] + 
        #   as.matrix(df_prs_test_scaled)[,-1] %*% coeff_LASSO[-1]
        # residual_test_LASSO = adjusted_pheno_test - ss_reg_LASSO
        # r_square_prs_LASSO <- 1 - var(residual_test_LASSO) / var(adjusted_pheno_test)
        # 
        # combined_ways_result$r_square[combined_ways_index] <- r_square_prs_LASSO 
        # combined_ways_index <- combined_ways_index + 1
        # 
        # # Fit Ridge with cross-validation
        # cv_fit_Ridge <- cv.glmnet(x = x_all, y = y_all, alpha = 0, standardize = FALSE)
        # coeff_Ridge <- coef(cv_fit_Ridge, s = "lambda.min")[,1]
        # ss_reg_Ridge <- coeff_Ridge[1] + 
        #   as.matrix(df_prs_test_scaled)[,-1] %*% coeff_Ridge[-1]
        # residual_test_Ridge = adjusted_pheno_test - ss_reg_Ridge
        # r_square_prs_Ridge <- 1 - var(residual_test_Ridge) / var(adjusted_pheno_test)
        # 
        # combined_ways_result$r_square[combined_ways_index] <- r_square_prs_Ridge 
        # combined_ways_index <- combined_ways_index + 1
        # Store combination results
        combined_ways_result_idx <- which(combined_ways_result$pop == pop & combined_ways_result$trait == trait & combined_ways_result$Epochs == epoch)
        combined_ways_result$JointPRS[combined_ways_result_idx] <- r_square_prs_results[["JointPRS"]][epoch]
        combined_ways_result$hm3_mega[combined_ways_result_idx] <- r_square_prs_results[["hm3_mega"]][epoch]
        combined_ways_result$JointPRS_annot[combined_ways_result_idx] <- r_square_prs_results[["NNOLS_all"]][epoch]
        combined_ways_result$JointPRS_annot_pen_LASSO[combined_ways_result_idx] <- r_square_prs_results[["LASSO_all"]][epoch]
        combined_ways_result$JointPRS_annot_pen_Ridge[combined_ways_result_idx] <- r_square_prs_results[["Ridge_all"]][epoch]
      }
      r_square_long <- melt(r_square_prs_results, id.vars = "Epochs", variable.name = "PRS_Type", value.name = "R_Squared")
      
      # # Reorder the PRS_Type in your plotting dataframe based on this new order
      # prs_levels <- c(ordered_annotations_prop.h[[trait]], "JointPRS", "weighted_all", "weighted_selected")
      
      # # Reorder the PRS_Type in your plotting dataframe based on the new combined order
      # r_square_long$PRS_Type <- factor(r_square_long$PRS_Type, levels = prs_levels)
      
      
      # Create the boxplot with reordered annotations
      # p <- ggplot(r_square_long, aes(x = PRS_Type, y = R_Squared)) +
      #   geom_boxplot() +
      #   labs(x = "PRS Type", y = "R-square") +
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
      #   filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/R2_",trait,"_",pop,".pdf"),
      #   plot = p,
      #   width = adaptive_width,
      #   height = adaptive_height,
      #   units = "in"
      # )
      
    }
  }
}

##### Plot comparison between different combination ways of JointPRS, hm3_mega, LASSO, Ridge
combined_ways_result_long <- combined_ways_result %>%
  pivot_longer(
    cols = c("JointPRS", "hm3_mega", "JointPRS_annot","JointPRS_annot_pen_LASSO","JointPRS_annot_pen_Ridge"),
    names_to = "PRS_type",
    values_to = "R2"
  )
combined_ways_result_long$PRS_type <- factor(combined_ways_result_long$PRS_type, levels = combined_ways)
combined_ways_result_long$trait <- factor(combined_ways_result_long$trait, levels = traits)
prs_colors <- c(
  "JointPRS" = "#1f77b4", 
  "hm3_mega" = "#ff7f0e", 
  "JointPRS_annot" = "#2ca02c",
  "JointPRS_annot_pen_LASSO" = "#d62728",
  "JointPRS_annot_pen_Ridge" = "#9467bd"
)

# Create plot with facets
combined_ways_plot <- ggplot(
  combined_ways_result_long, 
  aes(x = trait, y = R2, fill = PRS_type)
) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  facet_wrap(~ pop, ncol = 4, scales = "free_x") + # Facet by population
  labs(
    x = "Trait",
    y = "R2",
    fill = "PRS Type"
  ) +
  scale_fill_manual(
    values = prs_colors,
    labels = c(
      "JointPRS", 
      "hm3_mega", 
      "JointPRS_annot",
      "JointPRS_annot(LASSO)",
      "JointPRS_annot(Ridge)"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"), # Facet header styling
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines") # Space between facets
  ) +
  guides(fill = guide_legend(nrow = 1)) # Horizontal legend
# Print the plot
print(combined_ways_plot)
ggsave(
  filename = paste0("/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/Figures/R2_minor_pop_JointPRS_annot.pdf"),
  plot = combined_ways_plot,
  width = 10,
  height = 5,
  units = "in"
)


