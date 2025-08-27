#-------------------------------------------------------------------------------
# title: Weighted least squares regression analysis
# author: E Lunsford  
# date: 2025-08-27
#  
# This code is conducted a weighted least squares regression analysis
#
#
# Last Run: 08/27/2025 and code was in working order 
# using R 4.5.1 and RStudio 2025.05.1+513 
#
#-------------------------------------------------------------------------------


#################################################################################
# Load Libraries.
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
library(spdep)
library(stringr)
library(tidyverse)
library(tidycensus)

#################################################################################
# Load data and filter for analysis.
# For this analysis, we are focused on all cbg with at least 1 resident
#################################################################################

load(file = "Data/acs_tcc_holc_final.RData")

# DIA CT are 83.88, 83.89, 9800.1, 158
# Largest DIA CT 9800.01
dia_id <- c("080319800011")

df_main <- acs_tcc_holc_final %>%
  filter(
    # Remove DIA
    !GEOID %in% dia_id,
    # Filter for all cbg with at least 1 resident
         total_pop > 0)

#################################################################################
# Categorical HOLC with quantitative TCC - comparison of means test
#################################################################################

# ANOVA
one.way <- aov(nlcd_mean ~ worst_holc,
               data = df_main)

summary(one.way)

# Check homoscedasity
par(mfrow = c(2,2))
plot(one.way)

two.way <- aov(nlcd_mean ~ worst_holc + hisp_per,
               data = df_main)

summary(two.way)
plot(two.way)

# Test SLR with all variables
model_slr <- lm(nlcd_mean ~ worst_holc + hisp_per + nhb_per + nhw_per + pov_per,
                data = df_main)

plot(model_slr)
par(mfrow=c(1,1))

plot(model_slr$fitted.values, model_slr$residuals)
abline(h = 0, col = "red")

#################################################################################
# Breusch-Pagan Test to test heteroskedasticity
#################################################################################

bptest(model_slr)
# P < 0.05 --> reject H0 and conclude that there is evidence of heteroskedasticity

summary(model_slr)

#################################################################################
# Test for normality of residuals
#################################################################################

res_slr <- residuals(model_slr)
shapiro.test(res_slr)
# P > 0.05 --> fail to reject H0 and conclude that there is no evidence against normality

qqnorm(res_slr)
qqline(res_slr, col = "red")

#################################################################################
# Test for Spatial autocorrelation using spdep package
#################################################################################

# variable of interest
var_to_test <- "nlcd_mean"

# Remove NAs
sf_clean <- df_main[!is.na(df_main[[var_to_test]]), ]

# Create spatial weights matrix
neighbors <- poly2nb(sf_clean)
weights_list <- nb2listw(neighbours = neighbors, style = "W")

# Test Moran's I
moran.test(x = sf_clean$nlcd_mean, listw = weights_list)
# I stat = 0.609 --> closer to 1 similar values (high–high or low–low) are clustered
# p-value < 2.2e-16 --> reject H0 and conclude that there is evidence of spatial autocorrelation

# Create Moran scatterplot
moran.plot(sf_clean[[var_to_test]], weights_list,
           labels = FALSE,
           xlab = "NLCD Tree Canopy",
           ylab = "Spatial lag",
           main = "Moran Scatterplot")

# Positive slope --> positive spatial correlation, cbg with higher TCC values tend to be 
# surrounded by other high TCC cbg, and low TCC areas tend to be near other low TCC areas.

# Quadrants:
# Upper right (High–High) and lower left (Low–Low) have a large share of points, clustering of similar values dominates the pattern.
# Upper left (Low–High) and lower right (High–Low) quadrants have fewer points — are spatial outliers where a cbg differs from its surroundings.
# Spread: points fairly dispersed but still align with positive slope, suggesting moderate spatial autocorrelation.


#################################################################################
# Model Set 1: SLR Loop of all combinations of + and *
# For this analysis, we are focused on all cbg with at least 1 resident
#################################################################################

# Define dependent variable
y_var <- "nlcd_mean"

# Define independent variables
x_vars <- c("worst_holc", "hisp_per", "nhb_per", "nhw_per", "pov_per")

# Generate all combinations for x variables
all_combinations <- unlist(lapply(1:length(x_vars), function(i) {
  combn(x_vars, i, simplify = FALSE)
}), recursive = FALSE)

# Create blank list to store results
  # For summary
  model_results_slr <- list()
  # For model objects
  results_model_slr <- list()

# Start loop for each combination
for(i in seq_along(all_combinations)) {
  
  # Main effects
  main_effects <- paste(all_combinations[[i]], collapse = " + ")
  
  # Interaction Terms
  if(length(all_combinations[[i]]) > 1) {
    interactions <- combn(all_combinations[[i]], 2, function(x) paste(x[1], x[2], sep = "*"),
                          simplify = TRUE)
    
    # Loop with and without interactions
    for(use_interactions in c(TRUE, FALSE)) {
      if(use_interactions) {
        interaction_terms <- paste(interactions, collapse = " + ")
        formula_str <- paste(y_var, "~", main_effects, "+", interaction_terms)
      } else {
        formula_str <- paste(y_var, "~", main_effects)
      }
      
      # Fit model
      model <- lm(as.formula(formula_str), data = df_main)
      
      # Store results
      model_results_slr[[formula_str]] <- summary(model)
      results_model_slr[[formula_str]] <- model
      
      # Track Progress
      message(paste("Model", i, "Completed:", formula_str))
    } 
  } else {
    
    # If only one variable, no interaction terms are possible
    formula_str <- paste(y_var, "~", main_effects)
    model <- lm(as.formula(formula_str), data = df_main)
    
    # Store results (both summary and model, same as above)
    model_results_slr[[formula_str]] <- summary(model)
    results_model_slr[[formula_str]] <- model
    
    message(paste("Model Completed:", formula_str))
  }
}

# Get names of all models
names(model_results_slr)

# View specific model results
knitr::kable(broom::tidy(model_results_slr[[1]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_slr[[2]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_slr[[3]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_slr[[4]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_slr[[5]], conf.int = TRUE))

################################################################################
# Select best models based on AIC using AICcmodavg package
################################################################################

# Get AIC table
aic_table_slr <- aictab(cand.set = results_model_slr)

# Get top 10 models by AIC
top10_names <- aic_table_slr$Modnames[1:10]
top10_models <- results_model_slr[top10_names]

# loop through bptest with top 10 models using lmtest package
bptest_results <- list()

for(name in top10_names) {
  model <- results_model_slr[[name]]
  test_result <- bptest(model)
  bptest_results[[name]] <- test_result
}

bptest_summary_slr <- data.frame(
  Modnames = names(bptest_results),
  BP_statistic = sapply(bptest_results, function(x) x$statistic),
  p_value = sapply(bptest_results, function(x) x$p.value)) %>%
  left_join(aic_table_slr, by = "Modnames")

#################################################################################
# Model Set 2:  WLS Loop of all combinations of + and *
# For this analysis, we are focused on all cbg with at least 1 resident
#################################################################################

# Define dependent variable
y_var <- "nlcd_mean"

# Define independent variables
x_vars <- c("worst_holc", "hisp_per", "nhb_per", "nhw_per", "pov_per")

# Variables allowed in interactions
interaction_vars <- c("hisp_per", "nhb_per", "nhw_per", "pov_per")

# Define weight
df_main$weight_invvar <- 1 / pmax(df_main$nlcd_variance, 1e-6)

# Generate all combinations for x variables
all_combinations <- unlist(lapply(1:length(x_vars), function(i) {
  combn(x_vars, i, simplify = FALSE)
}), recursive = FALSE)

# Create blank lists to store results
  # For summary
  model_results_wls <- list()
  # For model
  results_model_wls <- list()

# Loop for each combination
for (i in seq_along(all_combinations)) {
  
  # Main effects
  main_effects <- paste(all_combinations[[i]], collapse = " + ")
  
  # Interaction Terms: only between variables in `interaction_vars`
  inter_candidates <- intersect(all_combinations[[i]], interaction_vars)
  
  # Loop interactions
  if (length(inter_candidates) > 1) {
    interactions <- combn(inter_candidates, 2, function(x) paste(x[1], x[2], sep = "*"), simplify = TRUE)
    
    # Loop with and without interactions
        for (use_interactions in c(TRUE, FALSE)) {
      if (use_interactions) {
        interaction_terms <- paste(interactions, collapse = " + ")
        formula_str <- paste(y_var, "~", main_effects, "+", interaction_terms)
      } else {
        formula_str <- paste(y_var, "~", main_effects)
      }
      
      # Fit model
      model <- lm(as.formula(formula_str),
                  data = df_main,
                  weights = df_main[["weight_invvar"]])
      
      # Store results
      # As summary
      model_results_wls[[formula_str]] <- summary(model)
      # As model
      results_model_wls[[formula_str]] <- model
      
      # Track progress
      print(paste("Model", i, "Completed:", formula_str))
    }
  } else {
    
    # Only main effects possible
    formula_str <- paste(y_var, "~", main_effects)
    model <- lm(as.formula(formula_str),
                data = df_main,
                weights = df_main[["weight_invvar"]])
    
    # Store results
    # As summary
    model_results_wls[[formula_str]] <- summary(model)
    # As model
    results_model_wls[[formula_str]] <- model
    
    # Track progress
    print(paste("Model Completed:", formula_str))
  }
}

# Get names of all models
names(model_results_wls)

# View specific model results
knitr::kable(broom::tidy(model_results_wls[[1]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_wls[[6]], conf.int = TRUE))

#individual SES
knitr::kable(broom::tidy(model_results_wls[[2]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_wls[[3]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_wls[[4]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_wls[[5]], conf.int = TRUE))

################################################################################
# Select best models based on AIC using AICcmodavg package
################################################################################

# Get AIC table
aic_table_wls <- aictab(cand.set = results_model_wls)

# Get top 10 models by AIC
top10_names <- aic_table_wls$Modnames[1:10]
top10_models <- results_model_wls[top10_names]

# Loop through Breusch–Pagan test with top 10 models
bptest_results_wls <- list()
for (name in top10_names) {
  model <- top10_models[[name]]
  test_result <- bptest(model)
  bptest_results_wls[[name]] <- test_result
}

# Create summary table
bptest_summary_wls <- data.frame(
  Modnames = names(bptest_results_wls),
  BP_statistic = sapply(bptest_results_wls, function(x) x$statistic),
  p_value = sapply(bptest_results_wls, function(x) x$p.value),
  signif_counts = sapply(top10_models, function(mod) {
    sum(coef(summary(mod))[, "Pr(>|t|)"] < 0.05, na.rm = TRUE)
  })
) %>%
  left_join(aic_table_wls, by = "Modnames")


# View results
tidy(results_model_wls[[6]], conf.int = TRUE)

#model_6 <- as.data.frame(tidy(model_results_op1_wls[[6]], conf.int = TRUE))

#write_csv(x = model_6, file = "model_6.csv")


#################################################################################
# Model Set 3: Sensitivity Analysis A
# WLS Loop of all combinations of + and *
# Conduct WLS but using majority holc grade instead of worst holc
#################################################################################

# Define dependent variable
y_var <- "nlcd_mean"

# Define independent variables
x_vars <- c("majority_holc", "hisp_per", "nhb_per", "nhw_per", "pov_per")

# Variables allowed in interactions
interaction_vars <- c("hisp_per", "nhb_per", "nhw_per", "pov_per")

# Define weight
df_main$weight_invvar <- 1 / pmax(df_main$nlcd_variance, 1e-6)

# Generate all combinations for x variables
all_combinations <- unlist(lapply(1:length(x_vars), function(i) {
  combn(x_vars, i, simplify = FALSE)
}), recursive = FALSE)

# Create blank list to store results
  # For summary
  model_results_senA <- list()
  # For model
  results_model_senA <- list()

# Start loop for each combination
for(i in seq_along(all_combinations)) {
  
  # Main effects
  main_effects <- paste(all_combinations[[i]], collapse = " + ")
  
  # Interaction Terms: only between variables in `interaction_vars`
  inter_candidates <- intersect(all_combinations[[i]], interaction_vars)
  
  # Loop interactions
  if(length(all_combinations[[i]]) > 1) {
    interactions <- combn(all_combinations[[i]], 2, function(x) paste(x[1], x[2], sep = "*"),
                          simplify = TRUE)
    
    # Loop with and without interactions
    for(use_interactions in c(TRUE, FALSE)) {
      if(use_interactions) {
        interaction_terms <- paste(interactions, collapse = " + ")
        formula_str <- paste(y_var, "~", main_effects, "+", interaction_terms)
      } else {
        formula_str <- paste(y_var, "~", main_effects)
      }
      
      # Fit model
      model <- lm(as.formula(formula_str), 
                  data = df_main,
                  weights = df_main[["weight_invvar"]])
      
      # Store results
      # As Summary
      model_results_senA[[formula_str]] <- summary(model)
      # As models
      results_model_senA[[formula_str]] <- model
      
      # Track Progress
      print(paste("Model", i, "Completed:", formula_str))
    } 
  } else {
    
    # Only main effects possible
    formula_str <- paste(y_var, "~", main_effects)
    model <- lm(as.formula(formula_str), 
                data = df_main,
                weights = df_main[["weight_invvar"]])
    
    # Store results
    model_results_senA[[formula_str]] <- summary(model)
    
    # Track progress
    print(paste("Model Completed:", formula_str))
  }
}

# Get names of all models
names(model_results_senA)

# Select best models based on AIC using AICcmodavg package
# Get AIC table
aic_table_senA <- aictab(cand.set = results_model_senA)

# Get top 10 models by AIC
top10_names <- aic_table_senA$Modnames[1:10]
top10_models <- results_model_senA[top10_names]

# Loop through Breusch–Pagan test with top 10 models
bptest_results_senA <- list()
for (name in top10_names) {
  model <- top10_models[[name]]
  test_result <- bptest(model)
  bptest_results_senA[[name]] <- test_result
}

# Create summary table
bptest_summary_senA <- data.frame(
  Modnames = names(bptest_results_senA),
  BP_statistic = sapply(bptest_results_senA, function(x) x$statistic),
  p_value = sapply(bptest_results_senA, function(x) x$p.value),
  signif_counts = sapply(top10_models, function(mod) {
    sum(coef(summary(mod))[, "Pr(>|t|)"] < 0.05, na.rm = TRUE)
  })
) %>%
  left_join(aic_table_senA, by = "Modnames")

# View specific model results
knitr::kable(broom::tidy(model_results_senA[[1]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_wls[[1]], conf.int = TRUE))

#################################################################################
# Model Set 4: Sensitivity Analysis B
# WLS Loop of all combinations of + and *
# Conduct WLS but using remove all holc grade NA's
#################################################################################

worst_grade <- "NA"

df_sen <- acs_tcc_holc_final %>%
  filter(!GEOID %in% dia_id,
         total_pop > 0,
         !worst_holc %in% worst_grade,
         !majority_holc %in% worst_grade)

# Define dependent variable
y_var <- "nlcd_mean"

# Define independent variables
x_vars <- c("worst_holc", "hisp_per", "nhb_per", "nhw_per", "pov_per")

# Variables allowed in interactions
interaction_vars <- c("hisp_per", "nhb_per", "nhw_per", "pov_per")

# Define weight
df_sen$weight_invvar <- 1 / pmax(df_sen$nlcd_variance, 1e-6)

# Generate all combinations for x variables
all_combinations <- unlist(lapply(1:length(x_vars), function(i) {
  combn(x_vars, i, simplify = FALSE)
}), recursive = FALSE)

# Create blank list to store results
# For summary
  model_results_senB <- list()
# For model
  results_model_senB <- list()

# Start loop for each combination
for(i in seq_along(all_combinations)) {
  
  # Main effects
  main_effects <- paste(all_combinations[[i]], collapse = " + ")
  
  # Interaction Terms: only between variables in `interaction_vars`
  inter_candidates <- intersect(all_combinations[[i]], interaction_vars)
  
  # Interaction terms only if there are ≥ 2 allowed interaction variables
  if (length(all_combinations[[i]]) > 1) {
    if (length(inter_candidates) > 1) {
      interactions <- combn(inter_candidates, 2, function(x) paste(x[1], x[2], sep = "*"), simplify = TRUE)
    } else {
      interactions <- character(0)  # No interactions possible
    }
    
    # Loop with and without interactions
    for (use_interactions in c(TRUE, FALSE)) {
      if (use_interactions && length(interactions) > 0) {
        interaction_terms <- paste(interactions, collapse = " + ")
        formula_str <- paste(y_var, "~", main_effects, "+", interaction_terms)
      } else {
        formula_str <- paste(y_var, "~", main_effects)
      }
      
      # Fit model
      model <- lm(as.formula(formula_str), 
                  data = df_sen,
                  weights = df_sen[["weight_invvar"]])
      
      # Store results
      # As Summary
      model_results_senB[[formula_str]] <- summary(model)
      # As models
      results_model_senB[[formula_str]] <- model
      
      # Track Progress
      print(paste("Model", i, "Completed:", formula_str))
    } 
  } else {
    
    # Only main effects
    formula_str <- paste(y_var, "~", main_effects)
    model <- lm(as.formula(formula_str), 
                data = df_sen,
                weights = df_sen[["weight_invvar"]])
    
    # Store results
    model_results_senB[[formula_str]] <- summary(model)
    results_model_senB[[formula_str]] <- model
    
    print(paste("Model Completed:", formula_str))
  }
}

# Get names of all models
names(model_results_senB)

# Select best models based on AIC using AICcmodavg package
# Get AIC table
aic_table_senB <- aictab(cand.set = results_model_senB)

# Get top 10 models by AIC
top10_names <- aic_table_senB$Modnames[1:10]
top10_models <- results_model_senB[top10_names]

# Loop through Breusch–Pagan test with top 10 models
bptest_results_senB <- list()
for (name in top10_names) {
  model <- top10_models[[name]]
  test_result <- bptest(model)
  bptest_results_senB[[name]] <- test_result
}

# Create summary table
bptest_summary_senB <- data.frame(
  Modnames = names(bptest_results_senB),
  BP_statistic = sapply(bptest_results_senB, function(x) x$statistic),
  p_value = sapply(bptest_results_senB, function(x) x$p.value),
  signif_counts = sapply(top10_models, function(mod) {
    sum(coef(summary(mod))[, "Pr(>|t|)"] < 0.05, na.rm = TRUE)
  })
) %>%
  left_join(aic_table_senB, by = "Modnames")

# View specific model results
knitr::kable(broom::tidy(model_results_senB[[1]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_wls[[1]], conf.int = TRUE))
knitr::kable(broom::tidy(model_results_senA[[1]], conf.int = TRUE))
