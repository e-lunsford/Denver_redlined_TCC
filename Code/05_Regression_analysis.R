#-------------------------------------------------------------------------------
# title: Regression analysis
# author: E Lunsford  
# date: 2025-09-10
#  
# This code is to apply regression models for the relationship:
# TCC ~ HOLC + SES variables. Models tested include OLS, WLS, SAR, and SEM
#
#
# Last Run: 09/10/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------



#################################################################################
#
# Load Libraries.
#
#################################################################################
library(AICcmodavg)
library(broom)
library(ggplot2)
library(ggspatial)
library(ggthemes)
library(knitr)
library(lmtest)
library(pacman)
library(sf)
library(spatialreg)
library(tidyverse)
library(dplyr)
library(lmtest)
library(spdep)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)

#################################################################################
#
# Load data and filter for analysis.
# For this analysis, we are focused on all cbg with at least 1 resident.
#
#################################################################################

load(file = "Data/combined_sf.RData")

# DIA CT are 83.88, 83.89, 9800.1, 158
# Largest DIA CT 9800.01
dia_id <- c("080319800011")

# Filter data for main analysis
df_main <- combined_sf %>%
  filter(
    # Remove DIA
    !GEOID %in% dia_id,
    # Filter for all cbg with at least 1 resident
    total_pop > 0)

# Quick look at tree canopy across holc grade
df_main %>%
  group_by(worst_holc) %>%
  summarise(
    mean_canopy = mean(nlcd_mean, na.rm = TRUE),
    median_canopy = median(nlcd_mean, na.rm = TRUE),
    n = n())

#################################################################################
#
# Test model assumptions
#
#################################################################################

# Test SLR with all variables
model_slr <- lm(formula = nlcd_mean ~ worst_holc + hisp_per + nhb_per + nhw_per + pov_per,
                data = df_main)

# Set graphics view
par(mfrow = c(2,2))

# view model plots
plot(model_slr)

# Reset graphics view
par(mfrow=c(1,1))

# View fitted vs residuals plot
plot(x = model_slr$fitted.values, 
     y = model_slr$residuals)
abline(h = 0, col = "red")

# Breusch-Pagan Test to test heteroskedasticity from lmtest package
bptest(model_slr)
# P < 0.05 --> reject H0 and conclude that there is evidence of heteroskedasticity

# Test for normality of residuals from stats package
res_slr <- residuals(model_slr)
shapiro.test(res_slr)
# P > 0.05 --> fail to reject H0 and conclude that there is no evidence against normality

# View QQ Plot of residuals
qqnorm(res_slr)
qqline(res_slr, col = "red")

#################################################################################
#
# Test for Spatial autocorrelation using spdep package
#
#################################################################################

# variable of interest
var_to_test <- "nlcd_mean"

# Remove NAs
sf_clean <- df_main[!is.na(df_main[[var_to_test]]), ]

# Create spatial weights matrix
neighbors <- spdep::poly2nb(sf_clean)
weights_list <- spdep::nb2listw(neighbours = neighbors, style = "W")

# Test Moran's I
moran.test(x = sf_clean$nlcd_mean, listw = weights_list)
# I stat = 0.609 --> closer to 1 similar values are clustered
# p-value < 2.2e-16 --> reject H0 and conclude that there is evidence of spatial autocorrelation

# Moran scatterplot
moran.plot(sf_clean[[var_to_test]], weights_list,
           labels = FALSE,
           xlab = "NLCD Tree Canopy",
           ylab = "Spatial lag",
           main = "Moran Scatterplot")

# Positive slope --> positive spatial correlation

# Quadrants:
# Upper right (High–High) & lower left (Low–Low): most points - clustering dominates.
# Upper left (Low–High) & lower right (High–Low): fewer points — are outliers.
# Spread: points fairly dispersed & still align with positive slope, suggesting moderate spatial autocorrelation.

#################################################################################
#
# Spatial Autocorrelation of full models for OLS, WLS, SEM, and SAR
#
#################################################################################

# Cap large weights
df_main$weight_invvar_capped <- pmin(1 / pmax(df_main$nlcd_variance, 1e-6), 100)

# Create spatial weights
neighbors <- poly2nb(df_main)
weights_list <- nb2listw(neighbours = neighbors, style = "W")

# OLS model
ols_model <- lm(formula = nlcd_mean ~ worst_holc + hisp_per + nhb_per + nhw_per + pov_per,
                data = df_main)

summary(ols_model)

# WLS model
wls_model <- lm(formula = nlcd_mean ~ worst_holc + hisp_per + nhb_per + nhw_per + pov_per,
                data = df_main,
                weights = df_main[["weight_invvar_capped"]])

summary(wls_model)

# SER model
sem_model <- errorsarlm(formula = nlcd_mean ~ worst_holc + hisp_per + nhb_per + nhw_per + pov_per,
                        data = df_main,
                        listw = weights_list)

summary(sem_model)

# SAR model
sar_model <- lagsarlm(formula = nlcd_mean ~ worst_holc + hisp_per + nhb_per + nhw_per + pov_per,
                      data = df_main,
                      listw = weights_list)

summary(sar_model)

# Moran's I
# OLS (I = 0.404)
moran.test(x = residuals(ols_model), listw = weights_list)

# WLS (I = 0.563)
moran.test(x = residuals(wls_model), listw = weights_list)

# SEM (I = -0.055)
moran.test(x = residuals(sem_model), listw = weights_list)

# SAR (I = -0.063)
moran.test(x = residuals(sar_model), listw = weights_list)

#################################################################################
#
# Model Set 1: OLS & WLS Loop of all combinations of + and *
# For this analysis, we are focused on all cbg with at least 1 resident
#
#################################################################################

# Define dependent variable
y_var <- "nlcd_mean"

# Define independent variables
x_vars <- c("worst_holc", "hisp_per", "nhb_per", "nhw_per", "pov_per")

# Variables allowed in interactions
interaction_vars <- c("hisp_per", "nhb_per", "nhw_per", "pov_per")

# Cap extremely large weights
df_main$weight_invvar_capped <- pmin(1 / pmax(df_main$nlcd_variance, 1e-6), 100)

# Generate all combinations for x variables
all_combinations <- unlist(lapply(1:length(x_vars), function(i) {
  combn(x_vars, i, simplify = FALSE)
}), recursive = FALSE)

# Create blank list to store results
# For summary
model_results_reg <- list()
model_results_ols <- list()
model_results_wls <- list()

# For model objects
results_model_reg <- list()
results_model_ols <- list()
results_model_wls <- list()

# Start loop for each combination
for(i in seq_along(all_combinations)) {
  
  # Main effects
  main_effects <- paste(all_combinations[[i]], collapse = " + ")
  
  # Interaction terms
  inter_candidates <- intersect(all_combinations[[i]], interaction_vars)
  
  if(length(all_combinations[[i]]) > 1) {
    if (length(inter_candidates) > 1) {
      interactions <- combn(inter_candidates, 2, function(x) paste(x[1], x[2], sep = "*"), simplify = TRUE)
      form_set <- list(
        paste(y_var, "~", main_effects),
        paste(y_var, "~", main_effects, "+", paste(interactions, collapse = "+"))
      )
    } else {
      form_set <- list(paste(y_var, "~", main_effects))
    }
    
  # Calculate interaction model
  for(formula_str in form_set) {
    
    # OLS
    model_ols <- lm(formula = as.formula(formula_str), 
                    data = df_main)
    ols_r2 <- summary(model_ols)$r.squared
    ols_aic <- AIC(model_ols)
    ols_bp <- bptest(model_ols)$p.value
    ols_moran <- moran.test(residuals(model_ols), listw = weights_list)$estimate[1]
    
    # WLS
    model_wls <- lm(formula = as.formula(formula_str), 
                    data = df_main,
                    weights = df_main[["weight_invvar_capped"]])
    wls_r2 <- summary(model_wls)$r.squared
    wls_aic <- AIC(model_wls)
    wls_bp <- bptest(model_wls)$p.value
    wls_moran <- moran.test(residuals(model_wls), listw = weights_list)$estimate[1]

    # Store comparisons
    model_results_reg[[formula_str]] <- data.frame(
        Model = c("OLS", "WLS"),
        R2 = c(ols_r2, wls_r2),
        AIC = c(ols_aic, wls_aic),
        BP_pval = c(ols_bp, wls_bp),
        moran_i = c(ols_moran, wls_moran)
      )
    
    # Store results
    results_model_ols[[formula_str]] <- model_ols
    results_model_wls[[formula_str]] <- model_wls
      
    # Track progress
    print(paste("Model", i, "Completed:", formula_str))
    }
  } else {
    
    # Calculate main effects model
    formula_str <- paste(y_var, "~", main_effects)
    
    # OLS
    model_ols <- lm(formula = as.formula(formula_str), 
                    data = df_main)
    
    # WLS
    model_wls <- lm(formula = as.formula(formula_str), 
                    data = df_main,
                    weights = df_main[["weight_invvar_capped"]])
    
    # Store results
    model_results_reg[[formula_str]] <- list(
      OLS = broom::tidy(model_ols, conf.int = TRUE), 
      WLS = broom::tidy(model_wls, conf.int = TRUE))
    
    results_model_reg[[formula_str]] <- list(
      OLS = model_ols,
      WLS = model_wls)
    
    results_model_ols[[formula_str]] <- model_ols
    results_model_wls[[formula_str]] <- model_wls
    
    # Track progress
    print(paste("Model Completed:", formula_str))
  }
}

# View model 1 for each
model_results_reg[[1]]
results_model_reg[[1]]

#################################################################################
#
# Model Set 2:  SEM and SAR Loop of all combinations of + and *
# For this analysis, we are focused on all cbg with at least 1 resident
#
#################################################################################

# Define dependent variable
y_var <- "nlcd_mean"

# Define independent variables
x_vars <- c("worst_holc", "hisp_per", "nhb_per", "nhw_per", "pov_per")

# Variables allowed in interactions
interaction_vars <- c("hisp_per", "nhb_per", "nhw_per", "pov_per")

# Cap large weights
df_main$weight_invvar_capped <- pmin(1 / pmax(df_main$nlcd_variance, 1e-6), 100)

# Generate combinations for x variables
all_combinations <- unlist(lapply(1:length(x_vars), function(i) {
  combn(x_vars, i, simplify = FALSE)
}), recursive = FALSE)

# Create spatial weights
sf_clean <- df_main[!is.na(df_main[[y_var]]), ]
neighbors <- poly2nb(sf_clean)
weights_list <- nb2listw(neighbors, style = "W")

# Create blank list to store results
# For summary
model_results_spat <- list()
model_results_sem <- list()
model_results_sar <- list()

# For model objects
results_model_spat <- list()
results_model_sem <- list()
results_model_sar <- list()

# Start loop for each combination
for (i in seq_along(all_combinations)) {
  
  # Main effects
  main_effects <- paste(all_combinations[[i]], collapse = " + ")
  
  # Interaction terms
  inter_candidates <- intersect(all_combinations[[i]], interaction_vars)
  
  if(length(all_combinations[[i]]) > 1) {
    if (length(inter_candidates) > 1) {
      interactions <- combn(inter_candidates, 2, function(x) paste(x[1], x[2], sep = "*"), simplify = TRUE)
      form_set <- list(
        paste(y_var, "~", main_effects),
        paste(y_var, "~", main_effects, "+", paste(interactions, collapse = "+"))
      )
    } else {
      form_set <- list(paste(y_var, "~", main_effects))
    }
    
    # Calculate interaction model
    for (formula_str in form_set) {
      
      # SEM
      model_sem <- errorsarlm(formula = as.formula(formula_str), 
                              data = sf_clean, 
                              listw = weights_list)
      sem_r2 <- model_sem$NK # pseudo-R2
      sem_aic <- AIC(model_sem)
      sem_moran <- moran.test(residuals(model_sem), listw = weights_list)$estimate[1]
      
      # SAR
      model_sar <- lagsarlm(formula = as.formula(formula_str), 
                            data = sf_clean, 
                            listw = weights_list)
      sar_r2 <- model_sar$NK # pseudo-R2
      sar_aic <- AIC(model_sar)
      sar_moran <- moran.test(residuals(model_sar), listw = weights_list)$estimate[1]
      
      # Store comparisons
      model_results_spat[[formula_str]] <- data.frame(
        Model = c("SEM", "SAR"),
        AIC = c(sem_aic, sar_aic),
        Moran_I = c(sem_moran, sar_moran)
      )
      
      # Store model results
      results_model_spat[[formula_str]] <- list(
        SEM = model_sem,
        SAR = model_sar)
      
      results_model_sem[[formula_str]] <- model_sem
      results_model_sar[[formula_str]] <- model_sar
      
    # Track progress
     print(paste("Model", i, "Completed:", formula_str))
    }
  } else {
    
    # Calculate main effect model
    formula_str <- paste(y_var, "~", main_effects)
    
    # SEM
    model_sem <- errorsarlm(formula = as.formula(formula_str), 
                            data = sf_clean, 
                            listw = weights_list)
    sem_aic <- AIC(model_sem)
    
    # SAR
    model_sar <- lagsarlm(formula = as.formula(formula_str), 
                          data = sf_clean, 
                          listw = weights_list)
    sar_aic <- AIC(model_sar)
    
    # Store results
    model_results_spat[[formula_str]] <- list(
      SEM = broom::tidy(model_sem, conf.int = TRUE), 
      SAR = broom::tidy(model_sar, conf.int = TRUE))
    
    results_model_spat[[formula_str]] <- list(
      SEM = model_sem,
      SAR = model_sar)
    
    results_model_sem[[formula_str]] <- model_sem
    results_model_sar[[formula_str]] <- model_sar

    # Track progress
    print(paste("Model Completed:", formula_str))
  }
}

# View model 1 for OLS & WLS (reg) and SAR & SEM (spat)
model_results_reg[[1]]
model_results_spat[[1]]

#################################################################################
#
# Model selection.
# Select best model based on AIC, include BP (for OLS & WLS), and Moran's I
#
#################################################################################

# Create function to extract model stats
get_model_stats <- function(model, weights_list) {
  
  # Set model type
  model_type <- class(model)[1]
  
  # Get AIC value
  aic_val <- tryCatch(AIC(object = model), error = function(e) NA)

  # Moran's I on residuals
  moran_val <- tryCatch({
    moran.test(x = residuals(model), listw = weights_list)$estimate[["Moran I statistic"]]
  }, error = function(e) NA)
  
  # BP test p-value (only works for OLS/WLS)
  bp_val <- tryCatch({
    if ("sarlm" %in% class(model)) {
      NA    # if SEM or SAR skip
    } else {
      bptest(model)$p.value # if OLS or WLS
    }
  }, error = function(e) NA)

  # Return 3-length vector
  return(c(AIC = aic_val, 
           Moran_I = moran_val, 
           BP_pval = bp_val))
}

# Create function to count significant terms
count_sig_terms <- function(model) {
  pvals <- tryCatch({
    if ("sarlm" %in% class(model)) {
      NA   # if SEM/SAR skip
    } else {
      summary(model)$coefficients[, "Pr(>|t|)"] # if OLS or WLS
    }
  }, error = function(e) NA)
  
  # Ignore intercept
  if("(Intercept)" %in% rownames(coef(model))) pvals <- pvals[-1]

  sum(pvals[-1] < 0.05, na.rm = TRUE)
}

# Create function to make a table for all models
make_model_table <- function(model_list, model_type, weights_list) {
  # Transpose and apply fxs
  stats_mat <- t(sapply(model_list, function(mod) {
    
    # Apply get_model_stats fx
    stats <- get_model_stats(model = mod, weights_list = weights_list)
    
    # Apply count_sig_terms fx
    sig_count <- count_sig_terms(model = mod)
    
    # Keep stats and sig_count
    c(stats, sig_count)
  }))

  # Create DF of applied functions
  model_stats_df <- as.data.frame(x = stats_mat, 
                                  stringsAsFactors = FALSE)
  
  # Make column type numeric
  numeric_cols <- c("AIC", "Moran_I", "BP_pval", "sig_count")
  
  # Add column for model formula
  model_stats_df$formula = names(model_list)
  
  # Add column for model number
  model_stats_df$model_num = seq_along(model_list)
  
  # Add column for OLS/WLS/SEM/SAR type
  model_stats_df$type = model_type
  
  # Add column name for significant term count
  colnames(model_stats_df)[4] <- "signif_terms"

  # Final DF of functions
  model_stats_df
}

# Build AIC tables
  # OLS
  aic_table_ols <- make_model_table(model_list = results_model_ols,
                                    model_type = "OLS",
                                    weights_list = weights_list)
  # WLS
  aic_table_wls <- make_model_table(model_list = results_model_wls,
                                    model_type = "WLS",
                                    weights_list = weights_list)
  # SEM
  aic_table_sem <- make_model_table(model_list = results_model_sem,
                                    model_type = "SEM",
                                    weights_list = weights_list)
  # SAR
  aic_table_sar <- make_model_table(model_list = results_model_sar,
                                    model_type = "SAR",
                                    weights_list = weights_list)

# Combine all AIC tables
aic_table_all <- bind_rows(aic_table_ols, 
                           aic_table_wls, 
                           aic_table_sem, 
                           aic_table_sar)

# Format AIC table
aic_table_grouped <- aic_table_all %>%
  # Ensure AIC is a numeric column
  dplyr::mutate(AIC = as.numeric(AIC)) %>%
  # Group by formula
  group_by(formula) %>%
  # Add delta, weight, and rank
  dplyr::mutate(
    delta_aic = AIC - min(AIC, na.rm = TRUE),
    weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic), na.rm = TRUE),
    rank = rank(AIC, ties.method = "min")) %>%
  # Order by AIC value
  arrange(AIC) %>%
  ungroup()

# View results
head(aic_table_grouped, n = 5)

#################################################################################
#
# Sensitivity Analyses: Use top 10 best AIC selected models
#
# Sensitivity A: Conduct SAR - Use majority HOLC grade instead of worst HOLC grade
# Sensitivity B: Conduct SAR - Use worst HOLC grade & filter spatial area to only include cbgs with an assigned HOLC grade
#
#################################################################################

# Get top 10 formulas from AIC table
top10_formulas <- aic_table_grouped %>%
  arrange(AIC) %>%
  slice(1:10) %>%
  pull(formula) %>%
  unique()

#################################################################################
#
# Model Set 3: Sensitivity Analysis A
# Conduct SAR but using majority HOLC grade instead of worst HOLC on top 10 AIC
#
#################################################################################

# Replace worst_holc with majority_holc
top10_majority_formulas <- stringr::str_replace(string = top10_formulas,
                                                pattern = "worst_holc",
                                                replacement = "majority_holc")

# Create blank list to store results
results_sar_senA <- list()

# Loop through formulas and run SAR only
for (i in seq_along(top10_majority_formulas)) {
  formula_str <- top10_majority_formulas[i]
  
  # Fit SAR model
  model_sar <- lagsarlm(formula = as.formula(formula_str),
                        data = sf_clean,
                        listw = weights_list)

  # Store model
  results_sar_senA[[formula_str]] <- model_sar
  
  # Track progress
  print(paste("Completed SAR model:", formula_str))
}

# View model estimates
# Best AIC[[1]]
sar_model46_senA <- as.data.frame(broom::tidy(results_sar_senA[[1]], conf.int = TRUE))

# Best AIC[[2]]
sar_model52_senA <- as.data.frame(broom::tidy(results_sar_senA[[2]], conf.int = TRUE))

# Run sensitivity on only HOLC as main effect
model1_sar_senA <- lagsarlm(formula = nlcd_mean ~ majority_holc,
                           data = sf_clean, 
                           listw = weights_list)

summary(model1_sar_senA)

#################################################################################
#
# Model Set 4: Sensitivity Analysis B
# Conduct SAR but remove all holc grade NA's
#
#################################################################################

# Sensitivity B data filter
worst_grade <- "no_holc"

# Filter for sensitivity B population of interest
df_senB <- df_main %>%
  filter(
    # Remove DIA
    !GEOID %in% dia_id,
    # Filter to at least 1 resident
    total_pop > 0,
    # Remove any cbg without an assigned HOLC grade
    !worst_holc %in% worst_grade)

# Create spatial weights for filtered data
sf_clean_senB <- df_senB[!is.na(df_senB[[y_var]]), ]
neighbors_senB <- poly2nb(sf_clean_senB)
weights_list_senB <- nb2listw(neighbors_senB, style = "W")

# Fit SAR models for top 10 formulas using sensitivity B population
# Create blank list to store results
results_sar_senB <- list()

# Loop through formulas and run SAR only
for (formula_str in top10_formulas) {
  
  # Fit SAR model
  model_sar <- lagsarlm(formula = as.formula(formula_str),
                        data = sf_clean_senB,
                        listw = weights_list_senB) 

  # Store model results
  results_sar_senB[[formula_str]] <- model_sar

  # Track progress
  print(paste("Completed SAR model (SenB):", formula_str))
}

# View model estimates
# Best AIC[[1]]
sar_model46_senB <- as.data.frame(broom::tidy(results_sar_senB[[1]], conf.int = TRUE))

# Best AIC[[2]]
sar_model52_senB <- as.data.frame(broom::tidy(results_sar_senB[[2]], conf.int = TRUE))

# Run sensitivity on only HOLC as main effect
model1_sar_senB <- lagsarlm(formula = nlcd_mean ~ worst_holc,
                            data = df_senB, listw = weights_list_senB)

summary(model1_sar_senB)

#################################################################################
#
# Calculate spatial impact estimates (direct/indirect/total) via spatialreg impacts fx
# Use tidy fxs to pull model estimates
# Save results to CSV
#
#################################################################################

# Define function for extraction
tidy_sar_impacts <- function(model_obj, impacts_summary) {
  
  # Function to tidy spatial impacts (direct/indirect/total)
    tidy_part <- function(part, prefix) {
      
    # Create starting data frame
    df <- as.data.frame(part[[1]]) %>% 
      rownames_to_column("term")
    
    # Add confidence intervals (CIs)
    ci <- as.data.frame(part[[2]]) %>%
      rownames_to_column("term") %>%
      select(term, `2.5%`, `97.5%`)
     
    # Change column names for consistency
    colnames(df)[-1] <- stringr::str_c(prefix, sep = "_", colnames(df)[-1])
    colnames(ci)[-1] <- stringr::str_c(prefix, sep = "_", c("ci_lower", "ci_upper"))

    # Join means, SEs, and CIs into one DF
    left_join(x = df,
              y = ci, 
              by = "term")
  }
  
  # Apply tidy_part fx
    # Direct
    direct_df <- tidy_part(part = impacts_summary$direct_sum, prefix = "direct")
    # Indirect
    indirect_df <- tidy_part(part = impacts_summary$indirect_sum, prefix = "indirect")
    # Total
    total_df <- tidy_part(part = impacts_summary$total_sum, prefix = "total")
    
  # Join impacts into one DF
  impacts_df <- direct_df %>%
    left_join(indirect_df, by = "term") %>%
    left_join(total_df, by = "term")
  
  # Set tidy model object to pull model estimates
  tidy_all_coefs <- broom::tidy(model_obj, conf.int = TRUE)
  
  # Pull rho term
  rho_df <- tidy_all_coefs %>%
    filter(term == "rho") %>%
    select(term, estimate, conf.low, conf.high) %>%
    rename(
      coef_estimate = estimate,
      coef_lower_ci = conf.low,
      coef_upper_ci = conf.high)
  
  # Pull estimates for variable estimates
  coefs_df <- tidy_all_coefs %>%
    filter(term != "rho") %>%
    rename(
      coef_estimate = estimate,
      coef_se = std.error,
      coef_z = statistic,
      coef_pvalue = p.value,
      coef_lower_ci = conf.low,
      coef_upper_ci = conf.high
    )
  
  # Combine main coefs and rhos
  all_coefs_df <- bind_rows(coefs_df, rho_df)

  # Join model estimates and impacts
  full_join(x = all_coefs_df, 
            y = impacts_df, 
            by = "term")
}

# Create list of models and weights objects
model_list <- list(
  model46        = list(model = results_model_sar[[46]], impacts_w = weights_list),
  model52        = list(model = results_model_sar[[52]], impacts_w = weights_list),
  model1         = list(model = results_model_sar[[1]],  impacts_w = weights_list),
  model46_SensA  = list(model = results_sar_senA[[1]],   impacts_w = weights_list),
  model52_SensA  = list(model = results_sar_senA[[2]],   impacts_w = weights_list),
  model1_SensA   = list(model = model1_sar_senA,         impacts_w = weights_list),
  model46_SensB  = list(model = results_sar_senB[[1]],   impacts_w = weights_list_senB),
  model52_SensB  = list(model = results_sar_senB[[2]],   impacts_w = weights_list_senB),
  model1_SensB   = list(model = model1_sar_senB,         impacts_w = weights_list_senB)
)

# Loop over model list
all_impacts_df <- bind_rows(
  # Apply functions
  imap(model_list, function(x, model_name) {
    # Get spatial impacts
    imp <- impacts(obj = x$model,       # spatial regression object
                   listw = x$impacts_w, # spatial weights
                   R = 1000)            # number of simulations
    # Get tidy results 
    tidy_sar_impacts(model_obj = x$model, 
                     impacts_summary = summary(object = imp,
                                               zstats = TRUE, 
                                               p.vals = TRUE)) %>%
    # Add model name
    mutate(model = model_name)
  })
) %>%
  # Select columns for DF
  select(model, term, everything())

# View results
head(x = all_impacts_df, n = 10)

# Save as Rdata
save(x = all_impacts_df,
     file = "Outputs/all_impacts_df.RData")
  
# Save as csv
readr::write_csv(x = all_impacts_df,
                 file = "Outputs/all_impacts_df.csv")


#################################################################################
