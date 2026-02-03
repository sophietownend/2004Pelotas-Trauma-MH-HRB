#-------------------------------------------------------------------------------
# PELOTAS ANALYSES - TRAUMA-MH-CMD
#-------------------------------------------------------------------------------
# Complete case analysis
#-------------------------------------------------------------------------------
# Sophie Townend; 26/01/2026
#-------------------------------------------------------------------------------

library(dplyr)
library(sandwich)
library(lmtest)
library(ggplot2)
library(broom)
library(car)

###### Variable Dictionary ###### 

### Exposure - Trauma ###
# y11_cte_recode: cumulative trauma load up to age 6 (ordinal none/1/2/3+)
# y15_cte_recode: cumulative trauma load up to age 15 for additional analyses

### Mediators - MH ###
# y15_externalising: externalising problems at age 15 (continuous)
# y15_internalising: internalising problems at age 15 (continuous)

### Outcomes - Health Risk Behaviours ###
# y18_alcohol: current problematic alcohol use at age 18 (binary)
# y18_smoking: current smoking at age 18 (binary)
# y18_druguse_recode: current illicit drug use at age 18 (binary) 
# joverall_pa: total physical activity at 18 (Total acceleration Euclidean Norm Minus One; continuous)
# jpercultraprocessadoscal: Proportion from ultra-processed foods consumption relative to the total energy intake at 18 (in calories; continuous)
# jmvpab5: Moderate to vigorous physical activity at 18 (5 min bouted; continous) in sensivity analyses

### Confounders ###
# asexo: adolescent sex (binary, 0=female; 1=male)
# aethnicity: adolescent ethnicity (binary)
# ae83: maternal alcohol consumption during pregnancy (binary)
# afumott: maternal smoking during pregnancy (binary)
# aescmae: maternal education years at birth (continuous)
# arendtot: monthly family income at birth (continous)
# birthday: day of the year that the adolescent was born (from 1 to 365; proxy for cohort birth order; continuous)

### Auxiliary variables for multiple imputation ###
# y6_cte_recode: cumulative trauma load up to age 6 (ordinal none/1/2/3+)
# y15_cte_recode: cumulative trauma load up to age 15 (ordinal none/1/2/3+)
# y11_internalising: internalising problems at age 11 (continuous)
# y11_externalising: externalising problems at age 11 (continuous)
# y18_internalising: internalising problems at age 18 (continuous)
# y18_externalising: externalising problems at age 18 (continuous)
# goverall_pa: total physical activity at 11 (Total acceleration Euclidean Norm Minus One; continuous)
# gmvpab5: Moderate to vigorous physical activity at 11 (5 min bouted, time spent in minutes; continous)
# gpercultraprocessados: Proportion from ultra-processed foods consumption relative to the total energy intake at 11 (in calories; continuous)

###### 1. Complete case analyses: Preparation ###### 
dat <- read.csv("trauma-MH-CMD_260126.csv")

# set variable types
# text covariates
character_variables <- c('idalt20231018')
dat[character_variables] <- lapply(dat[character_variables],as.character)
str(dat[character_variables])

# factor covariates
factor_variables <- c('asexo','aethnicity','ae83','afumott','y11_alltrauma','y6_alltrauma','y15_alltrauma',
                      'y18_alcohol','y18_smoking','y18_druguse_recode','y11_cte','y6_cte','y15_cte')
dat[factor_variables] <- lapply(dat[factor_variables],as.factor)
str(dat[factor_variables])

# numeric covariates - note, including cumulative trauma here, linear assumption holds
numeric_variables <- c('aescmae','arendtot','edi3m','birthday','y15_internalising','y15_externalising',
                       'y11_internalising','y11_externalising','y18_internalising','y18_externalising',
                       'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                       'gpercultraprocessados','y11_cte_recode','y15_cte_recode')
dat[numeric_variables] <- sapply(dat[numeric_variables],as.numeric)
str(dat[numeric_variables])

# recode binary outcome levels (0/1) 
# smoking
dat$y18_smoking <- recode_factor(dat$y18_smoking,
                                   "Not a current smoker" = 0,
                                   "Current smoker" = 1)
table(dat$y18_smoking)
# alcohol
dat$y18_alcohol <- recode_factor(dat$y18_alcohol,
                                   "No problematic alcohol use" = 0,
                                   "Problematic alcohol use" = 1)
table(dat$y18_alcohol)

# recode other categorical variables with labels
# sex
dat$asexo <- recode_factor(dat$asexo,
                          "female" = 0,
                          "male" = 1)
table(dat$asexo)
# maternal alcohol
dat$ae83 <- recode_factor(dat$ae83,
                           "n„o" = 0,
                           "sim" = 1)
table(dat$ae83)
# maternal smoking
dat$afumott <- recode_factor(dat$afumott,
                          "No" = 0,
                          "Yes" = 1)
table(dat$afumott)

### create a transformed externalising variable for the ext-alcohol models
dat$y15_externalising_transformed <- ((dat$y15_externalising + 1) / 10)^-1
psych::describe(dat$y15_externalising)
psych::describe(dat$y15_externalising_transformed)


### define complete cases:
# base on exposure, mediators, confounders, all outcomes (note, excluding bouted PA - sensitivity)
vars_to_check <- c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
dat$complete_case <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")

complete_dat <- dat[complete.cases(dat[, vars_to_check]), ]
# N = 1174

# save a recoded dat file
write.csv(dat, "trauma-MH-CMD_recode_260126.csv", row.names = F)
# save a complete dat file
write.csv(complete_dat, "trauma-MH-CMD_complete_dat_260126.csv", row.names = F)



#### 1.1 Get summary info for Table 1 descriptives in complete case #######

#### Descriptive statistics (Table 1)
continuous_vars = c('y11_cte_recode','y15_internalising','y15_externalising',
              'y11_internalising','aescmae','arendtot','edi12m','birthday',
              'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
              'gpercultraprocessados','y11_externalising','y18_internalising','y18_externalising','y15_cte_recode','y6_cte_recode')
binary_vars = c('asexo','aethnicity','ae83','afumott',
             'y18_alcohol','y18_smoking','y18_druguse_recode','y15_alcohol','y15_smoking',
             'y15_druguse','y11_alcohol','y11_smoking')

# Function for binary variables
describe_binary <- function(x, var_name) {
  n_total <- length(x)
  n_missing <- sum(is.na(x))
  n_present <- sum(!is.na(x))
  missing_pct <- 100 * n_missing / n_total
  
  x_non_missing <- x[!is.na(x)]
  n_endorsed <- sum(x_non_missing == 1)
  pct_endorsed <- 100 * n_endorsed / n_present
  ci <- prop.test(n_endorsed, n_present)$conf.int * 100
  
  data.frame(
    variable = var_name,
    type = "binary",
    n_present = n_present,
    n_missing = n_missing,
    pct_missing = round(missing_pct, 2),
    n_endorsed = n_endorsed,
    pct_endorsed = round(pct_endorsed, 2),
    ci_lower = round(ci[1], 2),
    ci_upper = round(ci[2], 2),
    mean = NA, sd = NA, se = NA, mean_ci_lower = NA, mean_ci_upper = NA
  )
}

# Function for continuous variables
describe_continuous <- function(x, var_name) {
  n_total <- length(x)
  n_missing <- sum(is.na(x))
  n_present <- sum(!is.na(x))
  missing_pct <- 100 * n_missing / n_total
  
  x_non_missing <- x[!is.na(x)]
  m <- mean(x_non_missing)
  s <- sd(x_non_missing)
  se <- s / sqrt(n_present)
  ci <- m + c(-1, 1) * 1.96 * se
  
  data.frame(
    variable = var_name,
    type = "continuous",
    n_present = n_present,
    n_missing = n_missing,
    pct_missing = round(missing_pct, 2),
    n_endorsed = NA,
    pct_endorsed = NA,
    ci_lower = NA,
    ci_upper = NA,
    mean = round(m, 2),
    sd = round(s, 2),
    se = round(se, 2),
    mean_ci_lower = round(ci[1], 2),
    mean_ci_upper = round(ci[2], 2)
  )
}

# Main function to apply summaries
summarise_selected_vars <- function(complete_dat, binary_vars, continuous_vars) {
  binary_summaries <- lapply(binary_vars, function(var) {
    describe_binary(complete_dat[[var]], var)
  }) %>% bind_rows()
  
  continuous_summaries <- lapply(continuous_vars, function(var) {
    describe_continuous(complete_dat[[var]], var)
  }) %>% bind_rows()
  
  bind_rows(binary_summaries, continuous_summaries)
}

results <- summarise_selected_vars(complete_dat, binary_vars, continuous_vars)

# Export to CSV
write.csv(results, "descriptive_statistics_complete_case.csv", row.names = FALSE)


###### 2. Complete case analyses: Regressions XM ###### 

#### 2.1. Linear regression XM (continuous outcomes) - adjusted ####
# Define outcome and predictor variables
outcomes <- c("y15_internalising", "y15_externalising")
predictors <- c("y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_XM_continuous_adjusted.csv")


#plotting diagnostics
OutcomeList <- data.frame(rbind("y15_internalising", "y15_externalising"))

pdf("diagnostic-plots/lm_diagnostics_XM_continuous_adjusted.pdf")
par(mfrow=c(2,2))

for(i in 1:nrow(OutcomeList)){
  Y <- complete_dat[,which(colnames(complete_dat) %in% as.character(OutcomeList[i,1]))]
  plot(lm(Y ~ as.numeric(complete_dat[,"y11_cte_recode"]) + complete_dat[,"asexo"] + complete_dat[,"aethnicity"] + complete_dat[,"ae83"] + complete_dat[,"afumott"] + complete_dat[,"aescmae"] + complete_dat[,"arendtot"] + complete_dat[,"birthday"]),main=OutcomeList[i,1])
}

dev.off()


#### 2.2. Linear regression XM (continuous outcomes) - unadjusted ####
# Define outcome and predictor variables
outcomes <- c("y15_internalising", "y15_externalising")
predictors <- c("y11_cte_recode")

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_XM_continuous_unadjusted.csv")


#plotting diagnostics 
OutcomeList <- data.frame(rbind("y15_internalising", "y15_externalising"))

pdf("diagnostic-plots/lm_diagnostics_XM_continuous_unadjusted.pdf")
par(mfrow=c(2,2))

for(i in 1:nrow(OutcomeList)){
  Y <- complete_dat[,which(colnames(complete_dat) %in% as.character(OutcomeList[i,1]))]
  plot(lm(Y ~ as.numeric(complete_dat[,"y11_cte_recode"])),main=OutcomeList[i,1])
}

dev.off()


###### 3. Complete case analyses: Regressions MY ###### 
## first include only int (or only ext) as a predictor (unadjusted model)
## THEN correct for baseline confounds and exposure (adjusted for baseline confounds), 
## THEN additionally adjust for the other mediator (fully adjusted)

#### 3.1. Linear regression MY (continuous outcomes) internalising - unadjusted ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_internalising")

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_continuous_unadjusted.csv")


#### 3.2. Linear regression MY (continuous outcomes) externalising - unadjusted ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_externalising")

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_continuous_unadjusted.csv")



#### 3.3. Linear regression MY (continuous outcomes) internalising - adjusted confounds ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_continuous_adjustedconf.csv")


#### 3.4. Linear regression MY (continuous outcomes) externalising - adjusted confounds ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_externalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_continuous_adjustedconf.csv")


#### 3.5. Linear regression MY (continuous outcomes) internalising - adjusted confounds + ext ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_internalising", "y15_externalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_continuous_adjustedconf_andExt.csv")


#### 3.6. Linear regression MY (continuous outcomes) externalising - adjusted confounds + int ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_externalising", "y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_continuous_adjustedconf_andInt.csv")



#### 3.7. Logistic regression MY (binary outcomes) internalising - unadjusted ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")  
predictors <- c("y15_internalising")

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_binary_unadjusted.csv")


#### 3.8. Logistic regression MY (binary outcomes) externalising - unadjusted ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_smoking", "y18_druguse_recode")  ### NOTE - remove alcohol, non-linear externalising
predictors <- c("y15_externalising") 

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_binary_unadjusted.csv")


#### 3.9. Logistic regression MY (binary outcomes) internalising - adjusted confounds ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")  
predictors <- c("y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_binary_adjusted_adjustedconf.csv")


#### 3.10. Logistic regression MY (binary outcomes) externalising - adjusted confounds ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_smoking", "y18_druguse_recode")   ## NOTE remove alcohol - non-linear
predictors <- c("y15_externalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_binary_adjusted_adjustedconf.csv")


#### 3.11. Logistic regression MY (binary outcomes) internalising - adjusted confounds + ext ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")  
predictors <- c("y15_internalising", "y15_externalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_binary_adjusted_adjustedconf_andExt.csv")


#### 3.12. Logistic regression MY (binary outcomes) externalising - adjusted confounds + int ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_smoking", "y18_druguse_recode")   ## NOTE remove alcohol - non-linear
predictors <- c("y15_externalising", "y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_binary_adjusted_adjustedconf_andInt.csv")



#### 3.13. Logistic regression MY: transformed ext / alcohol - unadjusted ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol")  
predictors <- c("y15_externalising_transformed")

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_transformedExt_alc_unadjusted.csv")


#### 3.14. Logistic regression MY: transformed ext / alcohol - adjusted confounds ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol")  
predictors <- c("y15_externalising_transformed", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_transformedExt_alc_adjustedconf.csv")


#### 3.15. Logistic regression MY: transformed ext / alcohol - adjusted confounds + int ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol")  
predictors <- c("y15_externalising_transformed","y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_transformedExt_alc_adjustedconf_andInt.csv")





###### 4. Complete case analyses: Sex interaction regressions XY ###### 

#### 4.1. Sex interaction XY (continuous outcomes) - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y11_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_continuous_adjusted.csv")


#### 4.2. Sex interaction XY (continuous outcomes) - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y11_cte_recode", "asexo")                          # Main effects
interaction_term <- c("y11_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo")                                           # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_continuous_unadjusted.csv")


#### 4.3. Sex interaction XY (binary outcomes) - adjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y11_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_binary_adjusted.csv")


#### 4.4. Sex interaction XY (binary outcomes) - unadjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y11_cte_recode", "asexo")                          # Main effects
interaction_term <- c("y11_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo")                                           # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_binary_unadjusted.csv")


###### 5. Complete case analyses: Sex interaction regressions XM ###### 

#### 5.1. Sex interaction XM (continuous outcomes) - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("y15_internalising", "y15_externalising")     
predictors <- c("y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y11_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XM_continuous_adjusted.csv")


#### 5.2. Sex interaction XM (continuous outcomes) - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("y15_internalising", "y15_externalising")     
predictors <- c("y11_cte_recode", "asexo")                          # Main effects
interaction_term <- c("y11_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo")                                           # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XM_continuous_unadjusted.csv")


###### 6. Complete case analyses: Sex interaction regressions MY ###### 

#### 6.1. Sex interaction MY - internalising interaction (continuous outcomes) - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_internalising", "y15_externalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_continuous_adjusted.csv")


#### 6.2. Sex interaction MY - externalising interaction (continuous outcomes) - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_externalising", "y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_continuous_adjusted.csv")


#### 6.3. Sex interaction MY - internalising interaction (continuous outcomes) - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_internalising", "asexo")                          # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_continuous_unadjusted.csv")


#### 6.4. Sex interaction MY - externalising interaction (continuous outcomes) - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_externalising", "asexo")                          # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_continuous_unadjusted.csv")


#### 6.5. Sex interaction MY - internalising interaction (binary outcomes) - adjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_internalising", "y15_externalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_binary_adjusted.csv")


#### 6.6. Sex interaction MY - externalising interaction (binary outcomes) - adjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_externalising", "y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_binary_adjusted.csv")


#### 6.7. Sex interaction MY - internalising interaction (binary outcomes) - unadjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_internalising", "asexo")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_binary_unadjusted.csv")


#### 6.8. Sex interaction MY - externalising interaction (binary outcomes) - unadjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_externalising", "asexo")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_binary_unadjusted.csv")


#### 6.9. Sex interaction MY - transformed externalising / alcohol - adjusted ####
outcomes <- c("y18_alcohol")     
predictors <- c("y15_externalising_transformed", "y15_internalising", "y11_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_externalising_transformed", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_transformedexternalising_alc_adjusted.csv")


#### 6.10. Sex interaction MY - transformed externalising / alcohol - unadjusted ####
outcomes <- c("y18_alcohol")     
predictors <- c("y15_externalising_transformed", "asexo")    # Main effects
interaction_term <- c("y15_externalising_transformed", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_transformedexternalising_alc_unadjusted.csv")



#### 7.1. (post hoc) Simple effects of sex interaction MY - internalising*sex for drug use - adjusted ####
# Define variables
outcome <- "y18_druguse_recode"                      # Binary outcome (0/1)
continuous_var <- "y15_internalising"            # Continuous predictor
factor_var <- "asexo"                         # Grouping factor
covariates <- c("y15_externalising","aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")           # Additional covariates
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott") # Convert to factor

# Ensure factor variables are coded correctly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Get levels of the factor
levels_group <- levels(complete_dat[[factor_var]])

# Storage for results
results_list <- list()

for (level in levels_group) {
  # Subset data to one level of the factor
  subset_df <- complete_dat %>% filter(!!as.name(factor_var) == level)
  
  # Fit logistic regression model
  predictors <- c(continuous_var, covariates)
  formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  model <- glm(formula, data = subset_df, family = binomial)
  
  # Get robust SEs
  vcov_robust <- vcovHC(model, type = "HC1")
  coefs <- coef(model)
  se_robust <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  # Sample size
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  # Loop through coefficients
  for (i in 2:length(coefs)) {  # Skip intercept
    coef_name <- names(coefs)[i]
    beta <- coefs[i]
    robust_se <- se_robust[i]
    z_value <- beta / robust_se
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds Ratio
    OR <- exp(beta)
    OR_CI_Lower <- exp(ci_lower)
    OR_CI_Upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Group_Level = level,
      Predictor = coef_name,
      Beta = beta,
      RobustSE = robust_se,
      zValue = z_value,
      pValue = p_value,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_CI_Lower,
      OR_CI_Upper = OR_CI_Upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "simple_effects_int_by_sex_y18_druguse.csv", row.names = FALSE)


#### 7.2. (post hoc) Simple effects of sex interaction MY - ext*sex for overall PA - adjusted ####
# Define variables
outcome <- "joverall_pa"                      # Continuous outcome
continuous_var <- "y15_externalising"         # Continuous predictor
factor_var <- "asexo"                         # Grouping variable (moderator)
covariates <- c("y15_internalising", "y11_cte_recode", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")           # Additional covariates
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott") # Convert to factor

# Ensure factors
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Levels of the factor
levels_group <- levels(complete_dat[[factor_var]])

# Store results
results_list <- list()

# Loop over each level of the factor
for (level in levels_group) {
  # Subset data to current group level
  subset_df <- complete_dat %>% filter(!!as.name(factor_var) == level)
  
  # Define formula: outcome ~ continuous_var + covariates
  formula <- as.formula(paste(outcome, "~", paste(c(continuous_var, covariates), collapse = " + ")))
  
  # Fit model
  model <- lm(formula, data = subset_df)
  summary_model <- summary(model)
  
  # Robust SE
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Group_Level = level,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "simple_effects_externalising_by_sex_joverall_pa.csv", row.names = FALSE)

## plot
# Define variables
outcome <- "joverall_pa"
continuous_var <- "y15_externalising"
group_var <- "asexo"  # Factor with 2+ levels

# Ensure group is a factor
complete_dat[[group_var]] <- as.factor(complete_dat[[group_var]])

pdf("results/sex_interactions/simple_effects_externalising_by_sex_joverall_pa.pdf")
ggplot(complete_dat, aes_string(x = continuous_var, y = outcome, color = group_var)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  labs(
    title = paste("Simple Slopes of", continuous_var, "by", group_var),
    x = continuous_var,
    y = outcome,
    color = group_var
  ) +
  theme_minimal(base_size = 14)
dev.off()

#### 7.3. (post hoc) Simple effects of sex interaction MY - ext*sex for PA bouted - adjusted ####
# Define variables
outcome <- "jmvpab5"                          # Continuous outcome
continuous_var <- "y15_externalising"         # Continuous predictor
factor_var <- "asexo"                         # Grouping variable (moderator)
covariates <- c("y15_internalising", "y11_cte_recode", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")           # Additional covariates
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott") # Convert to factor

# Ensure factors
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Levels of the factor
levels_group <- levels(complete_dat[[factor_var]])

# Store results
results_list <- list()

# Loop over each level of the factor
for (level in levels_group) {
  # Subset data to current group level
  subset_df <- complete_dat %>% filter(!!as.name(factor_var) == level)
  
  # Define formula: outcome ~ continuous_var + covariates
  formula <- as.formula(paste(outcome, "~", paste(c(continuous_var, covariates), collapse = " + ")))
  
  # Fit model
  model <- lm(formula, data = subset_df)
  summary_model <- summary(model)
  
  # Robust SE
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Group_Level = level,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "simple_effects_externalising_by_sex_jmvpab5.csv", row.names = FALSE)

## plot
# Define variables
outcome <- "jmvpab5"
continuous_var <- "y15_externalising"
group_var <- "asexo"  # Factor with 2+ levels

# Ensure group is a factor
complete_dat[[group_var]] <- as.factor(complete_dat[[group_var]])

pdf("results/sex_interactions/simple_effects_externalising_by_sex_jmvpab5.pdf")
ggplot(complete_dat, aes_string(x = continuous_var, y = outcome, color = group_var)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  labs(
    title = paste("Simple Slopes of", continuous_var, "by", group_var),
    x = continuous_var,
    y = outcome,
    color = group_var
  ) +
  theme_minimal(base_size = 14)
dev.off()



# MEDIATION ANALYSES -> ALL RUN IN STATA  --------------------------------
###### 8. Complete case analyses: Preparation of data for Stata ###### 
dat <- read.csv("trauma-MH-CMD_complete_dat_260126.csv")

# set variable types
# text covariates
character_variables <- c('idalt20231018')
dat[character_variables] <- lapply(dat[character_variables],as.character)
str(dat[character_variables])

# factor covariates
factor_variables <- c('asexo','aethnicity','ae83','afumott','y11_alltrauma','y6_alltrauma','y15_alltrauma',
                      'y18_alcohol','y18_smoking','y18_druguse_recode','y11_cte','y6_cte','y15_cte')
dat[factor_variables] <- lapply(dat[factor_variables],as.factor)
str(dat[factor_variables])

# numeric covariates - note, including cumulative trauma here, linear assumption holds
numeric_variables <- c('aescmae','arendtot','edi3m','birthday','y15_internalising','y15_externalising',
                       'y11_internalising','y11_externalising','y18_internalising','y18_externalising',
                       'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                       'gpercultraprocessados','y11_cte_recode','y15_cte_recode','y15_externalising_transformed')
dat[numeric_variables] <- sapply(dat[numeric_variables],as.numeric)
str(dat[numeric_variables])

# define complete cases:
# base on exposure, mediators, confounders, all outcomes (note, excluding bouted PA - sensitivity - and trauma up to 15 - just mediation)
vars_to_check <- c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
dat$complete_case <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")

# save a csv file for stata
write.csv(dat, "trauma-MH-CMD_complete_dat_stata_260126.csv", row.names = F)





# Complete case analyses for trauma up to age 15  --------------------------------

###### 1. Complete case analyses: Preparation ###### 
dat <- read.csv("trauma-MH-CMD_040425_UTF8.csv")

# set variable types
# text covariates
character_variables <- c('idalt20231018')
dat[character_variables] <- lapply(dat[character_variables],as.character)
str(dat[character_variables])

# factor covariates
factor_variables <- c('asexo','aethnicity','ae83','afumott','y11_alltrauma','y6_alltrauma','y15_alltrauma',
                      'y18_alcohol','y18_smoking','y18_druguse_recode','y11_cte','y6_cte','y15_cte')
dat[factor_variables] <- lapply(dat[factor_variables],as.factor)
str(dat[factor_variables])

# numeric covariates - note, including cumulative trauma here, linear assumption holds
numeric_variables <- c('aescmae','arendtot','edi3m','birthday','y15_internalising','y15_externalising',
                       'y11_internalising','y11_externalising','y18_internalising','y18_externalising',
                       'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                       'gpercultraprocessados','y11_cte_recode','y15_cte_recode')
dat[numeric_variables] <- sapply(dat[numeric_variables],as.numeric)
str(dat[numeric_variables])

# recode binary outcome levels (0/1) 
# drug use
dat$y18_druguse <- recode_factor(dat$y18_druguse_recode,
                                 "Not a current drug user" = 0,
                                 "Current drug user" = 1)
table(dat$y18_druguse)
# smoking
dat$y18_smoking <- recode_factor(dat$y18_smoking,
                                 "Not a current smoker" = 0,
                                 "Current smoker" = 1)
table(dat$y18_smoking)
# alcohol
dat$y18_alcohol <- recode_factor(dat$y18_alcohol,
                                 "No problematic alcohol use" = 0,
                                 "Problematic alcohol use" = 1)
table(dat$y18_alcohol)

# recode other variables (categorical variables with labels)
# sex
dat$asexo <- recode_factor(dat$asexo,
                           "female" = 0,
                           "male" = 1)
table(dat$asexo)
# maternal alcohol
dat$ae83 <- recode_factor(dat$ae83,
                          "n„o" = 0,
                          "sim" = 1)
table(dat$ae83)
# maternal smoking
dat$afumott <- recode_factor(dat$afumott,
                             "No" = 0,
                             "Yes" = 1)
table(dat$afumott)

##### create a transformed externalising variable for the ext-alcohol models
dat$y15_externalising_transformed <- ((dat$y15_externalising + 1) / 10)^-1
psych::describe(dat$y15_externalising)
psych::describe(dat$y15_externalising_transformed)

# define complete cases:
# base on exposure, mediators, confounders, all outcomes (note, excluding bouted PA - sensitivity - and trauma up to 15 - just mediation)
vars_to_check <- c("y15_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
dat$complete_case <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")

table(dat$complete_case)

complete_dat <- dat[complete.cases(dat[, vars_to_check]), ]

#### 1.1 Get summary info for Table 1 descriptives in complete case #######

#### Descriptive statistics (Table 1)
continuous_vars = c('y15_cte_recode','y15_internalising','y15_externalising',
                    'y11_internalising','aescmae','arendtot','edi12m','birthday',
                    'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                    'gpercultraprocessados','y11_externalising','y18_internalising','y18_externalising','y11_cte_recode','y6_cte_recode')
binary_vars = c('asexo','aethnicity','ae83','afumott',
                'y18_alcohol','y18_smoking','y18_druguse_recode','y15_alcohol','y15_smoking',
                'y15_druguse','y11_alcohol','y11_smoking','y15_alltrauma_recode')

# Function for binary variables (assumed coded as 0/1)
describe_binary <- function(x, var_name) {
  n_total <- length(x)
  n_missing <- sum(is.na(x))
  n_present <- sum(!is.na(x))
  missing_pct <- 100 * n_missing / n_total
  
  x_non_missing <- x[!is.na(x)]
  n_endorsed <- sum(x_non_missing == 1)
  pct_endorsed <- 100 * n_endorsed / n_present
  ci <- prop.test(n_endorsed, n_present)$conf.int * 100
  
  data.frame(
    variable = var_name,
    type = "binary",
    n_present = n_present,
    n_missing = n_missing,
    pct_missing = round(missing_pct, 2),
    n_endorsed = n_endorsed,
    pct_endorsed = round(pct_endorsed, 2),
    ci_lower = round(ci[1], 2),
    ci_upper = round(ci[2], 2),
    mean = NA, sd = NA, se = NA, mean_ci_lower = NA, mean_ci_upper = NA
  )
}

# Function for continuous variables
describe_continuous <- function(x, var_name) {
  n_total <- length(x)
  n_missing <- sum(is.na(x))
  n_present <- sum(!is.na(x))
  missing_pct <- 100 * n_missing / n_total
  
  x_non_missing <- x[!is.na(x)]
  m <- mean(x_non_missing)
  s <- sd(x_non_missing)
  se <- s / sqrt(n_present)
  ci <- m + c(-1, 1) * 1.96 * se
  
  data.frame(
    variable = var_name,
    type = "continuous",
    n_present = n_present,
    n_missing = n_missing,
    pct_missing = round(missing_pct, 2),
    n_endorsed = NA,
    pct_endorsed = NA,
    ci_lower = NA,
    ci_upper = NA,
    mean = round(m, 2),
    sd = round(s, 2),
    se = round(se, 2),
    mean_ci_lower = round(ci[1], 2),
    mean_ci_upper = round(ci[2], 2)
  )
}

# Main function to apply summaries
summarise_selected_vars <- function(complete_dat, binary_vars, continuous_vars) {
  binary_summaries <- lapply(binary_vars, function(var) {
    describe_binary(complete_dat[[var]], var)
  }) %>% bind_rows()
  
  continuous_summaries <- lapply(continuous_vars, function(var) {
    describe_continuous(complete_dat[[var]], var)
  }) %>% bind_rows()
  
  bind_rows(binary_summaries, continuous_summaries)
}

# Run the function on your data frame
# Replace your_data with your actual data frame name
results <- summarise_selected_vars(complete_dat, binary_vars, continuous_vars)

# Export to CSV
write.csv(results, "descriptive_statistics_complete_case_y15_cte.csv", row.names = FALSE)



###### 2. Complete case analyses: Regressions y15_cte XM ###### 

#### 2.1. Linear regression y15_cte XM (continuous outcomes) - adjusted ####
# Define outcome and predictor variables
outcomes <- c("y15_internalising", "y15_externalising")
predictors <- c("y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_XM_continuous_adjusted_y15_cte.csv")

#plotting diagnostics
OutcomeList <- data.frame(rbind("y15_internalising", "y15_externalising"))

pdf("lm_diagnostics_XM_continuous_adjusted_y15_cte.pdf")
par(mfrow=c(2,2))

for(i in 1:nrow(OutcomeList)){
  Y <- complete_dat[,which(colnames(complete_dat) %in% as.character(OutcomeList[i,1]))]
  plot(lm(Y ~ as.numeric(complete_dat[,"y15_cte_recode"]) + complete_dat[,"asexo"] + complete_dat[,"aethnicity"] + complete_dat[,"ae83"] + complete_dat[,"afumott"] + complete_dat[,"aescmae"] + complete_dat[,"arendtot"] + complete_dat[,"birthday"]),main=OutcomeList[i,1])
}

dev.off()


#### 2.2. Linear regression y15_cte XY (continuous outcomes) - unadjusted ####
# Define outcome and predictor variables
outcomes <- c("y15_internalising", "y15_externalising")
predictors <- c("y15_cte_recode")

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_XM_continuous_unadjusted_y15_cte.csv")


#plotting diagnostics 
OutcomeList <- data.frame(rbind("y15_internalising", "y15_externalising"))

pdf("lm_diagnostics_XM_continuous_unadjusted_y15_cte.pdf")
par(mfrow=c(2,2))

for(i in 1:nrow(OutcomeList)){
  Y <- complete_dat[,which(colnames(complete_dat) %in% as.character(OutcomeList[i,1]))]
  plot(lm(Y ~ as.numeric(complete_dat[,"y15_cte_recode"])),main=OutcomeList[i,1])
}

dev.off()



###### 3. Complete case analyses: Regressions y15_cte MY ###### 
#### 3.1. Linear regression MY (continuous outcomes) internalising - unadjusted ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_internalising")

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_continuous_unadjusted_confirmation.csv")


#### 3.2. Linear regression MY (continuous outcomes) externalising - unadjusted ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_externalising")

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_continuous_unadjusted_confirmation.csv")



#### 3.3. Linear regression MY (continuous outcomes) internalising - adjusted confounds ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_continuous_adjustedconf_y15_cte.csv")


#### 3.4. Linear regression MY (continuous outcomes) externalising - adjusted confounds ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_externalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_continuous_adjustedconf_y15_cte.csv")


#### 3.5. Linear regression MY (continuous outcomes) internalising - adjusted confounds + ext ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_internalising", "y15_externalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_continuous_adjustedconf_andExt_y15_cte.csv")


#### 3.6. Linear regression MY (continuous outcomes) externalising - adjusted confounds + int ####
# Define outcome and predictor variables
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")
predictors <- c("y15_externalising", "y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")

# Specify factor predictors (ensure proper variable type)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Convert to factor
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Build formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results list
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_continuous_adjustedconf_andInt_y15_cte.csv")



#### 3.7. Logistic regression MY (binary outcomes) internalising - unadjusted ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")  
predictors <- c("y15_internalising")

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_binary_unadjusted_confirmation.csv")


#### 3.8. Logistic regression MY (binary outcomes) externalising - unadjusted ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_smoking", "y18_druguse_recode")  ### NOTE - remove alcohol, non-linear
predictors <- c("y15_externalising") 

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_binary_unadjusted.csv")


#### 3.9. Logistic regression MY (binary outcomes) internalising - adjusted confounds ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")  
predictors <- c("y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_binary_adjusted_adjustedconf_y15_cte.csv")


#### 3.10. Logistic regression MY (binary outcomes) externalising - adjusted confounds ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_smoking", "y18_druguse_recode")   ## NOTE remove alcohol - non-linear
predictors <- c("y15_externalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_binary_adjusted_adjustedconf_y15_cte.csv")


#### 3.11. Logistic regression MY (binary outcomes) internalising - adjusted confounds + ext ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")  
predictors <- c("y15_internalising", "y15_externalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_INT_binary_adjusted_adjustedconf_andExt_y15_cte.csv")


#### 3.12. Logistic regression MY (binary outcomes) externalising - adjusted confounds + int ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_smoking", "y18_druguse_recode")   ## NOTE remove alcohol - non-linear
predictors <- c("y15_externalising", "y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_EXT_binary_adjusted_adjustedconf_andInt_y15_cte.csv")


#### 3.13. Logistic regression MY: transformed ext / alcohol - unadjusted ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol")  
predictors <- c("y15_externalising_transformed")

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_binary_transformedExt_alc_unadjusted_confirmation.csv")


#### 3.14. Logistic regression MY: transformed ext / alcohol - adjusted confounds ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol")  
predictors <- c("y15_externalising_transformed", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_binary_transformedExt_alc_adjustedconf_y15_cte.csv")


#### 3.15. Logistic regression MY: transformed ext / alcohol - adjusted confounds + int ####
# Define binary outcomes and predictors
binary_outcomes <- c("y18_alcohol")  
predictors <- c("y15_externalising_transformed","y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor formula string
predictor_formula <- paste(predictors, collapse = " + ")

# Initialize results container
results_list <- list()

for (outcome in binary_outcomes) {
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust covariance matrix
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))  # 95% CI
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and their 95% CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_MY_binary_transformedExt_alc_adjustedconf_andInt_y15_cte.csv")


###### 4. Complete case analyses: Sex interaction regressions XY y15_cte ###### 

#### 4.1. Sex interaction XY (continuous outcomes) y15_cte - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_continuous_adjusted_y15_cte.csv")


#### 4.2. Sex interaction XY (continuous outcomes) y15_cte - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_cte_recode", "asexo")                          # Main effects
interaction_term <- c("y15_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo")                                           # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_continuous_unadjusted_y15_cte.csv")


#### 4.3. Sex interaction XY (binary outcomes) y15_cte - adjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_binary_adjusted_y15_cte.csv")


#### 4.4. Sex interaction XY (binary outcomes) y15_cte - unadjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_cte_recode", "asexo")                          # Main effects
interaction_term <- c("y15_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo")                                           # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XY_binary_unadjusted_y15_cte.csv")


###### 5. Complete case analyses: Sex interaction regressions XM ###### 

#### 5.1. Sex interaction XM (continuous outcomes) y15_cte - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("y15_internalising", "y15_externalising")     
predictors <- c("y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XM_continuous_adjusted_y15_cte.csv")


#### 5.2. Sex interaction XM (continuous outcomes) y15_cte - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("y15_internalising", "y15_externalising")     
predictors <- c("y15_cte_recode", "asexo")                          # Main effects
interaction_term <- c("y15_cte_recode", "asexo")                    # Interaction (sex * trauma)
factor_vars <- c("asexo")                                           # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_XM_continuous_unadjusted_y15_cte.csv")


###### 6. Complete case analyses: Sex interaction regressions MY y15_cte ###### 

#### 6.1. Sex interaction MY - internalising interaction (continuous outcomes) y15_cte - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_internalising", "y15_externalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_continuous_adjusted_y15_cte.csv")


#### 6.2. Sex interaction MY - externalising interaction (continuous outcomes) y15_cte - adjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_externalising", "y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_continuous_adjusted_y15_cte.csv")


#### 6.3. Sex interaction MY - internalising interaction (continuous outcomes) (CONFIRMATION) - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_internalising", "asexo")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_continuous_unadjusted_confirmation.csv")


#### 6.4. Sex interaction MY - externalising interaction (continuous outcomes) CONFIRMATION - unadjusted ####
# Define continuous outcomes and predictors
outcomes <- c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5")     
predictors <- c("y15_externalising", "asexo")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- lm(formula, data = complete_dat)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_continuous_unadjusted_confirmation.csv")


#### 6.5. Sex interaction MY - internalising interaction (binary outcomes) y15_cte - adjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_internalising", "y15_externalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_binary_adjusted_y15_cte.csv")


#### 6.6. Sex interaction MY - externalising interaction (binary outcomes) y15_cte - adjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_externalising", "y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_binary_adjusted_y15_cte.csv")


#### 6.7. Sex interaction MY - internalising interaction (binary outcomes) - unadjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_internalising", "asexo")    # Main effects
interaction_term <- c("y15_internalising", "asexo")                    # Interaction (sex * internalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_internalising_MY_binary_unadjusted_confirmation.csv")


#### 6.8. Sex interaction MY - externalising interaction (binary outcomes) - unadjusted ####
outcomes <- c("y18_alcohol", "y18_smoking", "y18_druguse_recode")     
predictors <- c("y15_externalising", "asexo")    # Main effects
interaction_term <- c("y15_externalising", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_externalising_MY_binary_unadjusted_confirmation.csv")


#### 6.9. Sex interaction MY - transformed externalising / alcohol y15_cte - adjusted ####
outcomes <- c("y18_alcohol")     
predictors <- c("y15_externalising_transformed", "y15_internalising", "y15_cte_recode", "asexo", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")    # Main effects
interaction_term <- c("y15_externalising_transformed", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_transformedexternalising_alc_adjusted_y15_cte.csv")


#### 6.10. Sex interaction MY - transformed externalising / alcohol - unadjusted ####
outcomes <- c("y18_alcohol")     
predictors <- c("y15_externalising_transformed", "asexo")    # Main effects
interaction_term <- c("y15_externalising_transformed", "asexo")                    # Interaction (sex * externalising)
factor_vars <- c("asexo")          # Categorical predictors

# Ensure factors are treated properly
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Construct predictor and interaction strings
main_effects <- paste(predictors, collapse = " + ")
interaction_str <- paste(interaction_term, collapse = " * ")
full_formula <- paste(main_effects, "+", interaction_str)

# Initialize results container
results_list <- list()

# Loop over each outcome variable
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~", full_formula))
  model <- glm(formula, data = complete_dat, family = binomial)
  summary_model <- summary(model)
  
  # Robust SE and 95% CI
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- suppressMessages(coefci(model, vcov. = vcov_robust, level = 0.95))
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    zval <- summary_model$coefficients[i, "z value"]
    pval <- summary_model$coefficients[i, "Pr(>|z|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    # Odds ratios and CI
    OR <- exp(beta)
    OR_lower <- exp(ci_lower)
    OR_upper <- exp(ci_upper)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Outcome = outcome,
      Predictor = coef_name,
      LogOdds = beta,
      StdError = se,
      RobustSE = robust,
      zValue = zval,
      pValue = pval,
      LogOdds_CI_Lower = ci_lower,
      LogOdds_CI_Upper = ci_upper,
      OddsRatio = OR,
      OR_CI_Lower = OR_lower,
      OR_CI_Upper = OR_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "Results_sexinteraction_transformedexternalising_alc_unadjusted_confirmation.csv")




#### 7.1. (post hoc) Simple effects of sex interaction MY - ext*sex for overall PA - adjusted ####
# Define variables
outcome <- "joverall_pa"                      # Continuous outcome
continuous_var <- "y15_externalising"         # Continuous predictor
factor_var <- "asexo"                         # Grouping variable (moderator)
covariates <- c("y15_internalising", "y15_cte_recode", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")           # Additional covariates
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott") # Convert to factor

# Ensure factors
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Levels of the factor
levels_group <- levels(complete_dat[[factor_var]])

# Store results
results_list <- list()

# Loop over each level of the factor
for (level in levels_group) {
  # Subset data to current group level
  subset_df <- complete_dat %>% filter(!!as.name(factor_var) == level)
  
  # Define formula: outcome ~ continuous_var + covariates
  formula <- as.formula(paste(outcome, "~", paste(c(continuous_var, covariates), collapse = " + ")))
  
  # Fit model
  model <- lm(formula, data = subset_df)
  summary_model <- summary(model)
  
  # Robust SE
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Group_Level = level,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "simple_effects_externalising_by_sex_joverall_pa.csv", row.names = FALSE)

## plot
# Define variables
outcome <- "joverall_pa"
continuous_var <- "y15_externalising"
group_var <- "asexo"  # Factor with 2+ levels

# Ensure group is a factor
complete_dat[[group_var]] <- as.factor(complete_dat[[group_var]])

pdf("simple_effects_externalising_by_sex_joverall_pa.pdf")
ggplot(complete_dat, aes_string(x = continuous_var, y = outcome, color = group_var)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  labs(
    title = paste("Simple Slopes of", continuous_var, "by", group_var),
    x = continuous_var,
    y = outcome,
    color = group_var
  ) +
  theme_minimal(base_size = 14)
dev.off()

#### 7.2. (post hoc) Simple effects of sex interaction MY - ext*sex for PA bouted - adjusted ####
# Define variables
outcome <- "jmvpab5"                          # Continuous outcome
continuous_var <- "y15_externalising"         # Continuous predictor
factor_var <- "asexo"                         # Grouping variable (moderator)
covariates <- c("y15_internalising", "y15_cte_recode", "aethnicity", "ae83", "afumott", "aescmae", "arendtot", "birthday")           # Additional covariates
factor_vars <- c("asexo", "aethnicity", "ae83", "afumott") # Convert to factor

# Ensure factors
complete_dat[factor_vars] <- lapply(complete_dat[factor_vars], as.factor)

# Levels of the factor
levels_group <- levels(complete_dat[[factor_var]])

# Store results
results_list <- list()

# Loop over each level of the factor
for (level in levels_group) {
  # Subset data to current group level
  subset_df <- complete_dat %>% filter(!!as.name(factor_var) == level)
  
  # Define formula: outcome ~ continuous_var + covariates
  formula <- as.formula(paste(outcome, "~", paste(c(continuous_var, covariates), collapse = " + ")))
  
  # Fit model
  model <- lm(formula, data = subset_df)
  summary_model <- summary(model)
  
  # Robust SE
  vcov_robust <- vcovHC(model, type = "HC1")
  robust_se <- sqrt(diag(vcov_robust))
  ci_robust <- coefci(model, vcov. = vcov_robust, level = 0.95)
  
  n_obs <- nobs(model)
  df_resid <- df.residual(model)
  
  for (i in 2:length(coef(model))) {  # Skip intercept
    coef_name <- names(coef(model))[i]
    beta <- coef(model)[i]
    se <- summary_model$coefficients[i, "Std. Error"]
    tval <- summary_model$coefficients[i, "t value"]
    pval <- summary_model$coefficients[i, "Pr(>|t|)"]
    robust <- robust_se[i]
    ci_lower <- ci_robust[i, 1]
    ci_upper <- ci_robust[i, 2]
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Group_Level = level,
      Predictor = coef_name,
      Beta = beta,
      StdError = se,
      RobustSE = robust,
      tValue = tval,
      pValue = pval,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      DF = df_resid,
      N = n_obs
    )
  }
}

# Combine and export
results_df <- do.call(rbind, results_list)
write.csv(results_df, "simple_effects_externalising_by_sex_jmvpab5.csv", row.names = FALSE)

## plot
# Define variables
outcome <- "jmvpab5"
continuous_var <- "y15_externalising"
group_var <- "asexo"  # Factor with 2+ levels

# Ensure group is a factor
complete_dat[[group_var]] <- as.factor(complete_dat[[group_var]])

pdf("simple_effects_externalising_by_sex_jmvpab5.pdf")
ggplot(complete_dat, aes_string(x = continuous_var, y = outcome, color = group_var)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  labs(
    title = paste("Simple Slopes of", continuous_var, "by", group_var),
    x = continuous_var,
    y = outcome,
    color = group_var
  ) +
  theme_minimal(base_size = 14)
dev.off()


# MEDIATION ANALYSES FOR AGE 15 TRAUMA MODELS -> ALL RUN IN STATA  --------------------------------
###### 9. Complete case analyses: Preparation of data for Stata ###### 
dat <- read.csv("trauma-MH-CMD_complete_dat_y15_cte.csv")

# set variable types
# text covariates
character_variables <- c('idalt20231018')
dat[character_variables] <- lapply(dat[character_variables],as.character)
str(dat[character_variables])

# factor covariates
factor_variables <- c('asexo','aethnicity','ae83','afumott','y11_alltrauma','y6_alltrauma','y15_alltrauma',
                      'y18_alcohol','y18_smoking','y18_druguse_recode','y11_cte','y6_cte','y15_cte')
dat[factor_variables] <- lapply(dat[factor_variables],as.factor)
str(dat[factor_variables])

# numeric covariates - note, including cumulative trauma here, linear assumption holds
numeric_variables <- c('aescmae','arendtot','edi3m','birthday','y15_internalising','y15_externalising',
                       'y11_internalising','y11_externalising','y18_internalising','y18_externalising',
                       'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                       'gpercultraprocessados','y11_cte_recode','y15_cte_recode','y15_externalising_transformed')
dat[numeric_variables] <- sapply(dat[numeric_variables],as.numeric)
str(dat[numeric_variables])

# define complete cases:
# base on exposure, mediators, confounders, all outcomes (note, excluding bouted PA - sensitivity)
vars_to_check <- c("y15_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
dat$complete_case <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")

complete_dat <- dat[complete.cases(dat[, vars_to_check]), ]

# save a csv file for stata
write.csv(complete_dat, "trauma-MH-CMD_complete_dat_stata_y15_cte.csv", row.names = F)






