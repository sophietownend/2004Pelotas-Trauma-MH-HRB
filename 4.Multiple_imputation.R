#-------------------------------------------------------------------------------
# PELOTAS ANALYSES - TRAUMA-MH-CMD
#-------------------------------------------------------------------------------
# Multiple Imputation (redo)
#-------------------------------------------------------------------------------
# Sophie Townend; 26/01/26
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(sandwich)
library(lmtest)
library(ggplot2)
library(broom)
library(purrr)
library(readr)
library(dagitty)
library(midoc)
library(mice)
library(VIM)
library(readstata13)
library(DiagrammeR)
library(knitr)
library(smcfcs)
library(mitools)

###### Variable Dictionary ###### 
### Exposure - Trauma ###
# y11_cte_recode: cumulative trauma load up to age 6 (ordinal none/1/2/3+)
# y15_cte_recode: cumulative trauma load up to age 15 for additional analyses

### Mediators - MH ###
# y15_internalising: internalising problems at age 15 (continuous)
# y15_externalising: externalising problems at age 15 (continuous)

### Outcomes - Health Risk Behaviours ###
# y18_alcohol: current problematic alcohol use at age 18 (binary)
# y18_smoking: current smoking at age 18 (binary)
# y18_druguse_recode: current illicit drug use at age 18 (binary) 
# joverall_pa: total physical activity at 18 (Total acceleration Euclidean Norm Minus One; continuous)
# jpercultraprocessadoscal: Proportion from ultra-processed foods consumption relative to the total energy intake at 18 (in calories; continuous)
# jmvpab5: Moderate to vigorous physical activity at 18 (5 min bouted; continous) in sensivity analyses

### Confounders ###
# asexo: adolescent sex (binary, 0=female; 1=male)
# aethnicity: adolescent ethnicity (binary, white/other)
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
# gmvpab5: Moderate to vigorous physical activity at 18 (5 min bouted, time spent in minutes; continous)
# gpercultraprocessados: Proportion from ultra-processed foods consumption relative to the total energy intake at 11 (in calories; continuous)
# y15_alcohol: LIFETIME alcohol use, rather than current problematic at 18
# y15_smoking: LIFETIME smoking, rather than current at 18
# y15_druguse: LIFETIME drug use, rather than current at 18
# y11_alcohol: LIFETIME alcohol use, rather than current problematic at 18
# y11_smoking: LIFETIME smoking, rather than current at 18


###### 1. Data Preparation ###### 
dat <- read.csv("trauma-MH-CMD_recode_260126.csv")

# set variable types
# text covariates
character_variables <- c('idalt20231018')
dat[character_variables] <- lapply(dat[character_variables],as.character)
str(dat[character_variables])

# factor covariates
factor_variables <- c('asexo','aethnicity','ae83','afumott','y11_alltrauma','y6_alltrauma','y15_alltrauma',
                      'y18_alcohol','y18_smoking','y18_druguse_recode','y11_cte','y6_cte','y15_cte','y15_alcohol','y15_smoking',
                      'y15_druguse','y11_alcohol','y11_smoking')
dat[factor_variables] <- lapply(dat[factor_variables],as.factor)
str(dat[factor_variables])

# numeric covariates - note, including cumulative trauma here, linear assumption holds
numeric_variables <- c('aescmae','arendtot','edi3m','birthday','y15_internalising','y15_externalising',
                       'y11_internalising','y11_externalising','y18_internalising','y18_externalising',
                       'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                       'gpercultraprocessados','y11_cte_recode','y15_cte_recode','y6_cte_recode',
                       'aa07','y6_ctspc','y11_ctspc','y15_ctspc')
dat[numeric_variables] <- sapply(dat[numeric_variables],as.numeric)
str(dat[numeric_variables])

# check recoded binary outcome levels (0/1) 
table(dat$y18_druguse_recode)
table(dat$y18_smoking)
table(dat$y18_alcohol)
table(dat$asexo)
table(dat$ae83)
table(dat$afumott)
table(dat$aethnicity)

# recode binary auxiliary variables:
# y15_druguse
dat$y15_druguse <- recode_factor(dat$y15_druguse,
                                 "Never used drugs" = 0,
                                 "Used drugs at least once" = 1)
table(dat$y15_druguse)
# y15_smoking
dat$y15_smoking <- recode_factor(dat$y15_smoking,
                                 "Never smoked" = 0,
                                 "Smoked" = 1)
table(dat$y15_smoking)
# y15_alcohol
dat$y15_alcohol <- recode_factor(dat$y15_alcohol,
                                 "Never been drunk" = 0,
                                 "Has been drunk" = 1)
table(dat$y15_alcohol)
# y11_smoking
dat$y11_smoking <- recode_factor(dat$y11_smoking,
                                 "Never smoked" = 0,
                                 "Smoked" = 1)
table(dat$y11_smoking)
# y11_alcohol
dat$y11_alcohol <- recode_factor(dat$y11_alcohol,
                                 "Never had an alcoholic drink" = 0,
                                 "Drank alcohol" = 1)
table(dat$y11_alcohol)



###### 2. Data Exploration: missingness ###### 

# (a) Explore the missingness patterns
# explore missingness in exposure, outcome(s), and covariates.
md.pattern(dat[, c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")],rotate.names = TRUE)

pattern <- as.data.frame(md.pattern(dat[, c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                                            "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                                            "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")],rotate.names = TRUE))

write.csv(pattern,"missingness_pattern_analysis.csv")

# explore missingness in exposure, outcome(s), covariates, and auxiliary
md.pattern(dat[, c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal",
                   "y6_cte_recode",
                   "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                   "y11_alcohol","y11_smoking","y15_alcohol","y15_smoking","y15_druguse","goverall_pa","gpercultraprocessados",
                   "gmvpab5","aa07","y6_ctspc","y11_ctspc","y15_ctspc")],rotate.names = TRUE)

pattern <- as.data.frame(md.pattern(dat[, c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                                            "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                                            "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal",
                                            "y6_cte_recode",
                                            "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                                            "y11_alcohol","y11_smoking","y15_alcohol","y15_smoking","y15_druguse","goverall_pa","gpercultraprocessados",
                                            "gmvpab5","aa07","y6_ctspc","y11_ctspc","y15_ctspc")],rotate.names = TRUE))

write.csv(pattern,"missingness_pattern_auxiliary.csv")

# (b) Using the countNA function
# counts amount of missingness for each individual across all variables
misscount <- numeric(nrow(dat))
for (i in 1:nrow(dat)) {
  misscount[i] <- countNA(dat[i, c(
    "y11_cte_recode", "y15_externalising", "y15_internalising", 
    "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
    "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal"
  )])
}
# Display the counts in a table
table(misscount)

# percentages
round(table(misscount) / sum(table(misscount)) * 100, 2)
# even if CCA was valid, we'd only be basing it on ~28% of our sample - not efficient


########### 2.1. Summary of missing data ##############
# Load libraries
library(dplyr)
library(tidyr)

summarize_dataset <- function(data, cont_vars = NULL, bin_vars = NULL) {
  results <- list()
  total_n <- nrow(data)
  
  # Continuous variable summary
  if (!is.null(cont_vars)) {
    cont_summary <- data %>%
      summarise(across(all_of(cont_vars), list(
        n = ~sum(!is.na(.)),
        missing = ~sum(is.na(.)),
        missing_pct = ~sum(is.na(.)) / total_n * 100,
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        skewness = ~e1071::skewness(., na.rm = TRUE, type = 2),
        kurtosis = ~e1071::kurtosis(., na.rm = TRUE, type = 2)
      ), .names = "{.col}_{.fn}"))
    
    results$continuous <- cont_summary
  }
  
  # Binary variable summary
  if (!is.null(bin_vars)) {
    bin_summary <- data %>%
      summarise(across(all_of(bin_vars), list(
        n = ~sum(!is.na(.)),
        missing = ~sum(is.na(.)),
        missing_pct = ~sum(is.na(.)) / total_n * 100,
        endorsement_n = ~sum(. == 1, na.rm = TRUE),
        endorsement_pct = ~mean(. == 1, na.rm = TRUE) * 100
      ), .names = "{.col}_{.fn}"))
    
    results$binary <- bin_summary
  }
  
  return(results)
}

# Call the function
summary_whole <- summarize_dataset(dat,
                                   cont_vars = c('y11_cte_recode','y15_internalising','y15_externalising',
                                                 'y11_internalising','aescmae','arendtot','edi3m','birthday',
                                                 'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                                                 'gpercultraprocessados','y11_externalising','y18_internalising','y18_externalising','y15_cte_recode','y6_cte_recode'),
                                   bin_vars = c('asexo','aethnicity','ae83','afumott',
                                                'y18_alcohol','y18_smoking','y18_druguse_recode','y15_alcohol','y15_smoking',
                                                'y15_druguse','y11_alcohol','y11_smoking'))

# View results
summary_whole$continuous
summary_whole$binary

# convert to long format
cont_long <- tidyr::pivot_longer(summary_whole$continuous, 
                                 cols = everything(), 
                                 names_to = "variable_stat", 
                                 values_to = "value")

bin_long <- tidyr::pivot_longer(summary_whole$binary, 
                                cols = everything(), 
                                names_to = "variable_stat", 
                                values_to = "value")

# Write to CSV
write.csv(cont_long, "continuous_summary.csv", row.names = FALSE)
write.csv(bin_long, "binary_summary.csv", row.names = FALSE)

### 2.2. Summary of missing data with skewness/kurtosis p-values for continuous variables ####
# for skewness/kurtosis:
install.packages("e1071")
library(e1071)
# Function to compute skewness, kurtosis, and their p-values
compute_skew_kurt <- function(data, vars) {
  results <- lapply(vars, function(var) {
    x <- data[[var]]
    x_clean <- x[!is.na(x)]
    n <- length(x_clean)
    
    if (n < 8) {
      # Not enough data for stable estimates
      return(data.frame(
        variable = var,
        n = n,
        skewness = NA,
        skew_p = NA,
        kurtosis = NA,
        kurt_p = NA
      ))
    }
    
    # Compute statistics
    skew_val <- skewness(x_clean, type = 2)
    kurt_val <- kurtosis(x_clean, type = 2)
    
    # Standard errors under normality
    skew_se <- sqrt(6 / n)
    kurt_se <- sqrt(24 / n)
    
    # Z-scores and two-tailed p-values
    skew_z <- skew_val / skew_se
    kurt_z <- kurt_val / kurt_se
    
    skew_p <- 2 * (1 - pnorm(abs(skew_z)))
    kurt_p <- 2 * (1 - pnorm(abs(kurt_z)))
    
    data.frame(
      variable = var,
      n = n,
      skewness = skew_val,
      skew_p = skew_p,
      kurtosis = kurt_val,
      kurt_p = kurt_p
    )
  })
  
  # Combine into one data frame
  do.call(rbind, results)
}

result_df <- compute_skew_kurt(
  data = dat,
  vars = c('y11_cte_recode','y15_internalising','y15_externalising',
           'y11_internalising','aescmae','arendtot','edi3m','birthday',
           'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
           'gpercultraprocessados','y11_internalising','y11_externalising','y18_internalising','y18_externalising',
           'y15_cte_recode','y6_cte_recode')
)

# Save to CSV
write.csv(result_df, "skew_kurt_summary.csv", row.names = FALSE)


#### 2.3. Summary stats for each group (complete/incomplete) ######
dat$ethnicity_2 <- recode_factor(dat$aethnicity,
                                 "Other" = 0,
                                 "White" = 1)
table(dat$ethnicity_2)

summarize_by_group_simple <- function(data, group_var, cont_vars = NULL, bin_vars = NULL) {
  group_var <- rlang::sym(group_var)
  results <- list()
  
  # Compute group sizes to use for missing data stats
  group_sizes <- data %>%
    group_by(!!group_var) %>%
    summarise(group_n = n(), .groups = "drop")
  
  # Continuous summary
  if (!is.null(cont_vars)) {
    cont_summary <- data %>%
      group_by(!!group_var) %>%
      summarise(across(all_of(cont_vars), list(
        n = ~sum(!is.na(.)),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE)
      ), .names = "{.col}_{.fn}"), .groups = "drop") %>%
      left_join(group_sizes, by = rlang::as_string(group_var)) %>%
      mutate(across(ends_with("_n"), ~ group_n - ., .names = "{.col}_missing")) %>%
      mutate(across(ends_with("_missing"), ~ (. / group_n) * 100, .names = "{.col}_pct")) %>%
      select(-group_n)
    
    results$continuous <- cont_summary
  }
  
  # Binary summary
  if (!is.null(bin_vars)) {
    bin_summary <- data %>%
      group_by(!!group_var) %>%
      summarise(across(all_of(bin_vars), list(
        n = ~sum(!is.na(.)),
        endorsement_n = ~sum(. == 1, na.rm = TRUE),
        endorsement_pct = ~mean(. == 1, na.rm = TRUE) * 100
      ), .names = "{.col}_{.fn}"), .groups = "drop") %>%
      left_join(group_sizes, by = rlang::as_string(group_var)) %>%
      mutate(across(ends_with("_n"), ~ group_n - ., .names = "{.col}_missing")) %>%
      mutate(across(ends_with("_missing"), ~ (. / group_n) * 100, .names = "{.col}_pct")) %>%
      select(-group_n)
    
    results$binary <- bin_summary
  }
  
  # Merge summaries
  if (!is.null(cont_vars) & !is.null(bin_vars)) {
    merged <- left_join(results$continuous, results$binary, by = rlang::as_string(group_var))
    return(merged)
  } else if (!is.null(cont_vars)) {
    return(results$continuous)
  } else if (!is.null(bin_vars)) {
    return(results$binary)
  } else {
    stop("No variables specified for summarization.")
  }
}

# Run summary
summary_by_group <- summarize_by_group_simple(dat,
                                              group_var = "complete_case",
                                              cont_vars = c('y11_cte_recode','y15_internalising','y15_externalising',
                                                            'y11_internalising','aescmae','arendtot','edi3m','birthday',
                                                            'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                                                            'gpercultraprocessados','y11_externalising','y18_internalising','y18_externalising','y18_externalising','y15_cte_recode','y6_cte_recode'),
                                              bin_vars = c('asexo','aethnicity','ethnicity_2','ae83','afumott',
                                                           'y18_alcohol','y18_smoking','y18_druguse_recode', 'y15_alcohol','y15_smoking',
                                                           'y15_druguse','y11_alcohol','y11_smoking'))

summary_long <- summary_by_group %>%
  pivot_longer(
    cols = -complete_case, 
    names_to = c("variable", "stat"),
    names_pattern = "^(.*)_(n|mean|sd|endorsement_n|endorsement_pct|n_missing|n_missing_pct)$",
    values_to = "value"
  )

# Export to CSV
write.csv(summary_long, "summary_by_group.csv", row.names = FALSE)



#### 3. Associations with missingness ####
library(dplyr)
library(purrr)
library(broom)
library(readr)

# Make a complete case indicator 
# base on exposure, mediators, confounders, all outcomes 
vars_to_check <- c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
dat$complete_case <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")

# recoding a complete case variable so it's binary 0/1
dat$complete_case2 <- recode_factor(dat$complete_case,
                                    "complete" = 0,
                                    "incomplete" = 1)
table(dat$complete_case2)

###### complete case indicator for analysis variables excluding age 15 MH
vars_to_check <- c("y11_cte_recode", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
dat$complete_case_no15 <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")

# recoding a complete case variable so it's binary
dat$complete_case2_no15 <- recode_factor(dat$complete_case_no15,
                                         "complete" = 0,
                                         "incomplete" = 1)
table(dat$complete_case2_no15)


# Create MI data frame -- note, adding auxiliary variables (birth weight [aa07], total score on the Conflict Tactics Scale Parent-Child Version age 6 and 11)
MI_df <- dat[,c("idalt20231018","y11_cte_recode","asexo","aethnicity","ae83","afumott","aescmae",
                "arendtot","birthday","y15_internalising","y15_externalising","y18_alcohol",
                "y18_smoking","y18_druguse_recode","joverall_pa","jpercultraprocessadoscal","jmvpab5","y6_cte_recode",
                "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                "y11_alcohol","y11_smoking","y15_alcohol","y15_smoking","y15_druguse","goverall_pa","gpercultraprocessados",
                "gmvpab5","aa07","y6_ctspc","y11_ctspc","y15_ctspc")]
str(MI_df)
head(MI_df)

MI_df$completecase <- ifelse(apply(MI_df,1,anyNA)==F,1,0)
table(MI_df$completecase)

### complete case in MI dataframe
vars_to_check <- c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")

# Create new variable to flag complete cases
MI_df$complete_case <- ifelse(complete.cases(MI_df[, vars_to_check]), "complete", "incomplete")

# recoding a complete case variable
MI_df$complete_case2 <- recode_factor(MI_df$complete_case,
                                      "complete" = 0,
                                      "incomplete" = 1)
table(MI_df$complete_case2)


### Univariate association of each variable with missingness as an outcome 
# Define predictor variables
predictors <- c('y11_cte_recode','y15_internalising','y15_externalising',
                'aescmae','arendtot','birthday',
                'joverall_pa','jmvpab5', 'jpercultraprocessadoscal','goverall_pa','gmvpab5',
                'gpercultraprocessados','y11_internalising','y11_externalising','y18_internalising',
                'y18_externalising','y15_cte_recode','y6_cte_recode','asexo','aethnicity','ae83','afumott',
                'y18_alcohol','y18_smoking','y18_druguse_recode', 'y15_alcohol','y15_smoking',
                'y15_druguse','y11_alcohol','y11_smoking',"aa07","y6_ctspc","y11_ctspc","y15_ctspc") 

# Create an empty data frame to store results
results <- data.frame(
  Predictor = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  Significant = character(),
  Sample_Size = integer(),
  stringsAsFactors = FALSE
)

# Loop through each predictor and run logistic regression
for (var in predictors) {
  formula <- as.formula(paste("complete_case2 ~", var))
  model <- glm(formula, data = dat, family = binomial)
  
  # Extract OR, 95% CI, and p-value
  coef_summary <- summary(model)$coefficients
  est <- coef(model)[2]
  se <- coef_summary[2, 2]
  p_val <- coef_summary[2, 4]
  OR <- exp(est)
  CI <- exp(c(est - 1.96 * se, est + 1.96 * se))
  n <- nobs(model)  # Sample size
  
  # Determine significance
  sig_flag <- ifelse(p_val < 0.05, "Yes", "No")
  
  # Append to results
  results <- rbind(results, data.frame(
    Predictor = var,
    OR = round(OR, 3),
    CI_lower = round(CI[1], 3),
    CI_upper = round(CI[2], 3),
    p_value = round(p_val, 4),
    Significant = sig_flag,
    Sample_Size = n,
    stringsAsFactors = FALSE
  ))
}

write.csv(results, "univariate/univariate_regression_results_with_missingness_outcome.csv", row.names = FALSE) 



#### 4. Association of each confounder and auxiliary variable with exposure, mediators, and outcomes ######
#### trauma at age 11
predictors <- c('asexo','aethnicity','ae83','afumott','y15_alcohol','y15_smoking',
                'y15_druguse','y11_alcohol','y11_smoking','aescmae','arendtot','birthday','goverall_pa','gmvpab5',
                'gpercultraprocessados','y11_internalising','y11_externalising','y18_internalising','y18_externalising',
                'y15_cte_recode','y6_cte_recode',"aa07","y6_ctspc","y11_ctspc","y15_ctspc")
outcome_var <- "y11_cte_recode"

run_single_regression <- function(predictor, data, outcome) {
  formula <- as.formula(paste(outcome, "~", predictor))
  
  # Subset to complete cases for this model
  model_data <- data %>% select(all_of(c(outcome, predictor))) %>% na.omit()
  n_obs <- nrow(model_data)
  
  model <- lm(formula, data = model_data)
  
  tidy_out <- tidy(model, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = predictor,
      significant = ifelse(p.value < 0.05, TRUE, FALSE),
      n = n_obs
    )
  
  return(tidy_out)
}

# Run regressions and collect results
results <- lapply(predictors, run_single_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
dir.create("Associations_with_analysis_variables")
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_y11_cte_recode.csv")

# And for complete cases only
vars_to_check <- c("y11_cte_recode", "y15_externalising", "y15_internalising", 
                   "asexo","aethnicity","ae83","afumott","aescmae","arendtot","birthday", 
                   "y18_alcohol", "y18_smoking", "y18_druguse_recode", "joverall_pa", "jpercultraprocessadoscal")
dat$complete_case <- ifelse(complete.cases(dat[, vars_to_check]), "complete", "incomplete")
complete_dat <- dat[complete.cases(dat[, vars_to_check]), ]

results <- lapply(predictors, run_single_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_y11_cte_recode.csv")


#### internalising at age 15
outcome_var <- "y15_internalising"
# Run regressions and collect results
results <- lapply(predictors, run_single_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_y15_internalising.csv")

# And complete_dat
results <- lapply(predictors, run_single_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_y15_internalising.csv")


#### externalising at age 15
outcome_var <- "y15_externalising"
# Run regressions and collect results
results <- lapply(predictors, run_single_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_y15_externalising.csv")

# and complete dat
results <- lapply(predictors, run_single_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_y15_externalising.csv")


#### Overall PA at 18
outcome_var <- "joverall_pa"
# Run regressions and collect results
results <- lapply(predictors, run_single_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_joverall_pa.csv")

# complete dat
results <- lapply(predictors, run_single_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_joverall_pa.csv")


#### Bouted PA at 18
outcome_var <- "jmvpab5"
# Run regressions and collect results
results <- lapply(predictors, run_single_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_jmvpab5.csv")

# complete dat
results <- lapply(predictors, run_single_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_jmvpab5.csv")


#### UPF diet at 18
outcome_var <- "jpercultraprocessadoscal"
# Run regressions and collect results
results <- lapply(predictors, run_single_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_jpercultraprocessadoscal.csv")

# complete dat
results <- lapply(predictors, run_single_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_jpercultraprocessadoscal.csv")


#### alcohol use at 18
outcome_var <- "y18_alcohol"

# Function to run logistic regression and extract results
run_logistic_regression <- function(predictor, data, outcome) {
  formula <- as.formula(paste(outcome, "~", predictor))
  
  model_data <- data %>% select(all_of(c(outcome, predictor))) %>% na.omit()
  n_obs <- nrow(model_data)
  
  model <- glm(formula, data = model_data, family = binomial)
  
  tidy_out <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = predictor,
      significant = ifelse(p.value < 0.05, TRUE, FALSE),
      n = n_obs
    )
  
  return(tidy_out)
}

# Run regressions and combine results
results <- lapply(predictors, run_logistic_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_y18_alcohol.csv")

# complete dat
results <- lapply(predictors, run_logistic_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_y18_alcohol.csv")


#### smoking at 18
outcome_var <- "y18_smoking"
# Run regressions and combine results
results <- lapply(predictors, run_logistic_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_y18_smoking.csv")

# complete dat
results <- lapply(predictors, run_logistic_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_y18_smoking.csv")


#### drug use 18
outcome_var <- "y18_druguse_recode"
# Run regressions and combine results
results <- lapply(predictors, run_logistic_regression, data = dat, outcome = outcome_var)
results_df <- bind_rows(results)

# Reorder columns
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/regression_summary_y18_druguse.csv")

# complete dat
results <- lapply(predictors, run_logistic_regression, data = complete_dat, outcome = outcome_var)
results_df <- bind_rows(results)
results_df <- results_df %>%
  select(predictor, term, estimate, conf.low, conf.high, std.error, statistic, p.value, significant, n)

# Write to CSV
write_csv(results_df, "Associations_with_analysis_variables/completedat_regression_summary_y18_druguse.csv")





#### 5. Missing data patterns ####
# First, check: shouldn't include age 6 or age 15 CTE to impute age 11 CTE, too colinear (cumulative trauma)
# scatter of CTE at 6 v 11
ggplot(MI_df, aes(x = y11_cte_recode, y = y6_cte_recode)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of y11 vs y6",
       x = "x values",
       y = "y values")

ggplot(MI_df, aes(x = y11_cte_recode, y = y6_cte_recode)) +
  geom_jitter(width = 0.2, height = 0.2, color = "blue") +
  labs(title = "Jittered Scatter Plot")
# perfect triangle - don't include as auxiliary

# Second: Examine data added by auxiliary variables in missing data patterns
# Partially observed variables in my analysis model:
#  ○ Exposure (y11_cte_recode)
#  ○ Mediators (y15_internalising / externalising)
#  ○ Outcomes (drug use, alcohol, smoking, PA, UPF)
#  ○ Baseline confounders (aethnicity, aescmae)
#  ○ + all auxiliary have some missingness
# Fully observed variables in my analysis model:
#  ○ Confounders (asexo, ae83, afumott, arendtot, birthday)

### Exposure: trauma
md.pattern(dat[, c("y11_cte_recode", "y6_ctspc", "y11_ctspc")],rotate.names = TRUE)
md.pattern(dat[, c("y11_cte_recode", "y11_internalising", "y11_externalising")],rotate.names = TRUE)

### Mediators: y15_internalising y15_externalising
md.pattern(dat[, c("y15_internalising", "y15_externalising", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("y15_internalising", "y15_externalising", "y11_internalising", "y11_externalising", "y6_ctspc", "y11_ctspc")],rotate.names = TRUE)
md.pattern(dat[, c("y15_internalising", "y15_externalising", "y11_internalising", "y11_externalising", "y18_internalising", "y18_externalising")],rotate.names = TRUE)

### Outcomes: drug use, alcohol, smoking, PA, UPF
md.pattern(dat[, c("y18_druguse_recode", "y18_smoking","y18_alcohol","joverall_pa","jpercultraprocessadoscal")],rotate.names = TRUE)
### y18_druguse_recode
md.pattern(dat[, c("y18_druguse_recode", "y11_alcohol", "y11_smoking", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("y18_druguse_recode", "y11_alcohol", "y11_smoking", "y15_druguse", "y15_smoking", "y15_alcohol")],rotate.names = TRUE)
### y18_smoking
md.pattern(dat[, c("y18_smoking", "y11_alcohol", "y11_smoking", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("y18_smoking", "y11_alcohol", "y11_smoking", "y15_druguse", "y15_smoking", "y15_alcohol")],rotate.names = TRUE)
### y18_alcohol
md.pattern(dat[, c("y18_alcohol", "y11_alcohol", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("y18_alcohol", "y11_alcohol", "y15_druguse", "y15_smoking", "y15_alcohol")],rotate.names = TRUE)
### joverall_pa
md.pattern(dat[, c("joverall_pa", "goverall_pa", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("joverall_pa", "goverall_pa", "y6_ctspc", "y11_ctspc", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
### jmvpab5
md.pattern(dat[, c("jmvpab5", "gmvpab5", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("jmvpab5", "gmvpab5", "y6_ctspc", "y11_ctspc", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
### jpercultraprocessadoscal
md.pattern(dat[, c("jpercultraprocessadoscal", "gpercultraprocessados", "y11_internalising", "y11_externalising")],rotate.names = TRUE)
md.pattern(dat[, c("jpercultraprocessadoscal", "gpercultraprocessados", "y6_ctspc", "y11_ctspc", "y11_internalising", "y11_externalising")],rotate.names = TRUE)



# 6 Mice Imputation model - 75  --------------------------------
# Create MI data frame (exclude pp_ID, exclude y15_substance use, y15_ctspc, y6_cte_recode)
MI_df <- dat[,c("y11_cte_recode","asexo","aethnicity","ae83","afumott","aescmae",
                "arendtot","birthday","y15_internalising","y15_externalising","y18_alcohol",
                "y18_smoking","y18_druguse_recode","joverall_pa","jpercultraprocessadoscal","jmvpab5",
                "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                "y11_alcohol","y11_smoking","goverall_pa","gpercultraprocessados",
                "gmvpab5","aa07","y6_ctspc","y11_ctspc")]
str(MI_df)
head(MI_df)

#install.packages("mice")
library(mice)
str(MI_df)
options(max.print = 1000000)
getOption('max.print')

#Perform a "dryrun" 
testimp <- mice(data=MI_df, maxit = 0,
                printFlag=FALSE)
testimp

### define prediction matrix ####
predMatrix_1<-matrix(NA, nrow=29 ,ncol=29)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df)
predMatrix_1
predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aethnicity",]        <-c(1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aescmae",]           <-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["y15_internalising",] <-c(1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0) 
predMatrix_1["y15_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0) 
predMatrix_1["y18_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0) 
predMatrix_1["y18_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0) 
predMatrix_1["y18_druguse_recode",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0) 
predMatrix_1["joverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
predMatrix_1["jpercultraprocessadoscal",]         <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0) 
predMatrix_1["jmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0) 
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["y11_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1) 
predMatrix_1["y18_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y18_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["y11_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1) 
predMatrix_1["goverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["gpercultraprocessados",]             <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["gmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["aa07",]              <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0)   
predMatrix_1["y6_ctspc",]          <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1)
predMatrix_1["y11_ctspc",]         <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0)
predMatrix_1    

testimp2 <- mice(data=MI_df, predictorMatrix = predMatrix_1,
                 maxit = 0, printFlag=FALSE)
# check 
testimp2$predictorMatrix 
testimp2$formulas 
testimp2 # check methods 


### Imputation using stratified imputation ####
#Create separate male and female datasets
MI_df_female <- subset(MI_df, asexo==0)
MI_df_male <- subset(MI_df, asexo==1)

str(MI_df_female)
str(MI_df_male)
table(MI_df_female$asexo)
table(MI_df_male$asexo)

# dry run for each group (check prediction matrix and method)
testimp_female <- mice(data=MI_df_female, predictorMatrix = predMatrix_1,
                       maxit=0, printFlag=FALSE, seed=456)
testimp_male <- mice(data=MI_df_male, predictorMatrix = predMatrix_1,
                       maxit=0, printFlag=FALSE, seed=456)
testimp_female 
testimp_male

#Perform the imputation - maxit=25, m=75
bygroup_female_75<-mice(MI_df_female, maxit=25, m=75, seed=123, printFlag=FALSE, 
                        predictorMatrix = predMatrix_1)
bygroup_male_75<-mice(MI_df_male, maxit=25, m=75, seed=123, printFlag=FALSE, 
                      predictorMatrix = predMatrix_1)

# Combine the two sets of imputations
bygroup_75 <- rbind(bygroup_female_75, bygroup_male_75)

# extract first imputation (to check)
bygroup_75_1 <- complete(bygroup_75,1)
head(bygroup_75_1) 

dir.create("../mice_260126")
setwd("../mice_260126")
save(bygroup_75, file = "imp_75.rda")
#load("imp_75.rda")

### Analyse imputed dataset with test lm models ####
bygroup_est_joverall_pa <- with(bygroup_75, lm(joverall_pa ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jpercultraprocessadoscal <- with(bygroup_75, lm(jpercultraprocessadoscal ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jmvpab5 <- with(bygroup_75, lm(jmvpab5 ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y18_alcohol <- with(bygroup_75, glm(y18_alcohol ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_smoking <- with(bygroup_75, glm(y18_smoking ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_druguse <- with(bygroup_75, glm(y18_druguse_recode ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y15_externalising <- with(bygroup_75, lm(y15_externalising ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y15_internalising <- with(bygroup_75, lm(y15_internalising ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)
bygroup_est_y15_externalising_pool <- pool(bygroup_est_y15_externalising)
bygroup_est_y15_internalising_pool <- pool(bygroup_est_y15_internalising)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE)
summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE)
summary(bygroup_est_jmvpab5_pool, conf.int = TRUE)
summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE)
summary(bygroup_est_y18_smoking_pool, conf.int = TRUE)
summary(bygroup_est_y18_druguse_pool, conf.int = TRUE)
summary(bygroup_est_y15_externalising_pool, conf.int = TRUE)
summary(bygroup_est_y15_internalising_pool, conf.int = TRUE)

# Int/Ext 
# unadjusted - Int
bygroup_est_joverall_pa <- with(bygroup_75, lm(joverall_pa ~ y15_internalising))
bygroup_est_jpercultraprocessadoscal <- with(bygroup_75, lm(jpercultraprocessadoscal ~ y15_internalising))
bygroup_est_jmvpab5 <- with(bygroup_75, lm(jmvpab5 ~ y15_internalising))
bygroup_est_y18_alcohol <- with(bygroup_75, glm(y18_alcohol ~ y15_internalising, family = binomial))
bygroup_est_y18_smoking <- with(bygroup_75, glm(y18_smoking ~ y15_internalising, family = binomial))
bygroup_est_y18_druguse <- with(bygroup_75, glm(y18_druguse_recode ~ y15_internalising, family = binomial))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE)
summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE)
summary(bygroup_est_jmvpab5_pool, conf.int = TRUE)
summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE)
summary(bygroup_est_y18_smoking_pool, conf.int = TRUE)
summary(bygroup_est_y18_druguse_pool, conf.int = TRUE)

# unadjusted - Ext
bygroup_est_joverall_pa <- with(bygroup_75, lm(joverall_pa ~ y15_externalising))
bygroup_est_jpercultraprocessadoscal <- with(bygroup_75, lm(jpercultraprocessadoscal ~ y15_externalising))
bygroup_est_jmvpab5 <- with(bygroup_75, lm(jmvpab5 ~ y15_externalising))
bygroup_est_y18_alcohol <- with(bygroup_75, glm(y18_alcohol ~ y15_externalising, family = binomial))
bygroup_est_y18_smoking <- with(bygroup_75, glm(y18_smoking ~ y15_externalising, family = binomial))
bygroup_est_y18_druguse <- with(bygroup_75, glm(y18_druguse_recode ~ y15_externalising, family = binomial))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE)
summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE)
summary(bygroup_est_jmvpab5_pool, conf.int = TRUE)
summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE)
summary(bygroup_est_y18_smoking_pool, conf.int = TRUE)
summary(bygroup_est_y18_druguse_pool, conf.int = TRUE)

# adjusted for baseline confounders + trauma + one another
bygroup_est_joverall_pa <- with(bygroup_75, lm(joverall_pa ~ y15_internalising + y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jpercultraprocessadoscal <- with(bygroup_75, lm(jpercultraprocessadoscal ~ y15_internalising + y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jmvpab5 <- with(bygroup_75, lm(jmvpab5 ~ y15_internalising + y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y18_alcohol <- with(bygroup_75, glm(y18_alcohol ~ y15_internalising + y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_smoking <- with(bygroup_75, glm(y18_smoking ~ y15_internalising + y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_druguse <- with(bygroup_75, glm(y18_druguse_recode ~ y15_internalising + y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE)
summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE)
summary(bygroup_est_jmvpab5_pool, conf.int = TRUE)
summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE)
summary(bygroup_est_y18_smoking_pool, conf.int = TRUE)
summary(bygroup_est_y18_druguse_pool, conf.int = TRUE)


### The maximum number of iterations / traceplots ####
# how many iterations might be needed before the trend starts to stabilise
plot(bygroup_75, c("y11_cte_recode", "aethnicity", "aescmae"), layout = c(2, 4))
plot(bygroup_75, c("y15_internalising", "y15_externalising"), layout = c(2, 2))
plot(bygroup_75, c("y18_alcohol", "y18_smoking", "y18_druguse_recode"), layout = c(2, 3))
plot(bygroup_75, c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5"), layout = c(2, 3))


### Distribution of observed and imputed values - histograms ####
par(mfrow=c(1,3))
hist(bygroup_75$data$y11_cte_recode, xlab="Cumulative trauma at 11 years",
     main="Observed values")
abline(v=mean(bygroup_75$data$y11_cte_recode, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$y11_cte_recode[[1]], xlab="Cumulative trauma at 11 years",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$y11_cte_recode[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$y11_cte_recode), xlab="Cumulative trauma at 11 years",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$y11_cte_recode)),lty=2, col="red")

# internalising
par(mfrow=c(1,3))
hist(bygroup_75$data$y15_internalising, xlab="Internalising at 15 years",
     main="Observed values")
abline(v=mean(bygroup_75$data$y15_internalising, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$y15_internalising[[1]], xlab="Internalising at 15 years",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$y15_internalising[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$y15_internalising), xlab="Internalising at 15 years",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$y15_internalising)),lty=2, col="red")

# externalising
par(mfrow=c(1,3))
hist(bygroup_75$data$y15_externalising, xlab="Externalising at 15 years",
     main="Observed values")
abline(v=mean(bygroup_75$data$y15_externalising, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$y15_externalising[[1]], xlab="Externalising at 15 years",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$y15_externalising[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$y15_externalising), xlab="Externalising at 15 years",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$y15_externalising)),lty=2, col="red")

# aescmae
par(mfrow=c(1,3))
hist(bygroup_75$data$aescmae, xlab="Maternal education",
     main="Observed values")
abline(v=mean(bygroup_75$data$aescmae, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$aescmae[[1]], xlab="Maternal education",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$aescmae[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$aescmae), xlab="Maternal education",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$aescmae)),lty=2, col="red")

# joverall_pa
par(mfrow=c(1,3))
hist(bygroup_75$data$joverall_pa, xlab="Overall PA",
     main="Observed values")
abline(v=mean(bygroup_75$data$joverall_pa, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$joverall_pa[[1]], xlab="Overall PA",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$joverall_pa[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$joverall_pa), xlab="Overall PA",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$joverall_pa)),lty=2, col="red")

# jpercultraprocessadoscal
par(mfrow=c(1,3))
hist(bygroup_75$data$jpercultraprocessadoscal, xlab="UPF intake",
     main="Observed values")
abline(v=mean(bygroup_75$data$jpercultraprocessadoscal, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$jpercultraprocessadoscal[[1]], xlab="UPF intake",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$jpercultraprocessadoscal[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$jpercultraprocessadoscal), xlab="UPF intake",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$jpercultraprocessadoscal)),lty=2, col="red")

# jmvpab5
par(mfrow=c(1,3))
hist(bygroup_75$data$jmvpab5, xlab="bouted PA",
     main="Observed values")
abline(v=mean(bygroup_75$data$jmvpab5, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_75$imp$jmvpab5[[1]], xlab="bouted PA",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_75$imp$jmvpab5[[1]]),lty=2, col="red")

hist(unlist(bygroup_75$imp$jmvpab5), xlab="bouted PA",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_75$imp$jmvpab5)),lty=2, col="red")


### Strip plots ####
# Compare observed (blue) and imputed (red) data
stripplot(bygroup_75, y11_cte_recode + aethnicity + aescmae ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup_75, y15_internalising + y15_externalising ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup_75, y18_alcohol + y18_smoking + y18_druguse_recode + joverall_pa + jpercultraprocessadoscal + jmvpab5 ~ .imp,
          pch = 20, cex = 1.2)


### Density plots ####
densityplot(bygroup_75, data = ~ y11_cte_recode) 
densityplot(bygroup_75, data = ~ aethnicity) 
densityplot(bygroup_75, data = ~ aescmae) # See below for extra checks
densityplot(bygroup_75, data = ~ y15_internalising + y15_externalising) 
densityplot(bygroup_75, data = ~ y18_alcohol + y18_smoking + y18_druguse_recode) 
densityplot(bygroup_75, data = ~ joverall_pa + jpercultraprocessadoscal + jmvpab5) 



### The number of imputations (compare to m=50, m=100) ####
#Perform the imputations with different imputation (m) numbers
# 50
bygroup_female_50<-mice(MI_df_female, maxit=25, m=50, seed=123, printFlag=FALSE, 
                     predictorMatrix = predMatrix_1)
bygroup_male_50<-mice(MI_df_male, maxit=25, m=50, seed=123, printFlag=FALSE, 
                   predictorMatrix = predMatrix_1)

# Combine the two sets of imputations 
bygroup_50 <- rbind(bygroup_female_50, bygroup_male_50)
# extract first imputation (to check)
bygroup_50_1 <- complete(bygroup_50,1)
head(bygroup_50_1) 
save(bygroup_50, file = "imp_50.rda")

# 100
bygroup_female_100<-mice(MI_df_female, maxit=25, m=100, seed=123, printFlag=FALSE, 
                        predictorMatrix = predMatrix_1)
bygroup_male_100<-mice(MI_df_male, maxit=25, m=100, seed=123, printFlag=FALSE, 
                      predictorMatrix = predMatrix_1)

# Combine the two sets of imputations 
bygroup_100 <- rbind(bygroup_female_100, bygroup_male_100)
# extract first imputation (to check)
bygroup_100_1 <- complete(bygroup_100,1)
head(bygroup_100_1) 
save(bygroup_100, file = "imp_100.rda")

# compare differences in lm results across imputation numbers
# test - relationship between CTE and Ext
imp_m1_m75 <- with(bygroup_75,
                             lm(y15_externalising ~ y11_cte_recode + 
                                  asexo + aethnicity + ae83 + afumott + 
                                  aescmae + arendtot + birthday))
imp_m1_m75_res <- summary(pool(imp_m1_m75),
                                    conf.int=TRUE)[2,c(2,3,7,8)]

imp_m1_m50 <- with(bygroup_50,
                   lm(y15_externalising ~ y11_cte_recode + 
                        asexo + aethnicity + ae83 + afumott + 
                        aescmae + arendtot + birthday))
imp_m1_m50_res <- summary(pool(imp_m1_m50),
                          conf.int=TRUE)[2,c(2,3,7,8)]

imp_m1_m100 <- with(bygroup_100,
                             lm(y15_externalising ~ y11_cte_recode + 
                                  asexo + aethnicity + ae83 + afumott + 
                                  aescmae + arendtot + birthday))
imp_m1_m100_res <- summary(pool(imp_m1_m100),
                                    conf.int=TRUE)[2,c(2,3,7,8)]

# obtain the estimate in the complete dat data, when no values are missing. 
fulldatamod <- lm(y15_externalising ~ y11_cte_recode + 
                    asexo + aethnicity + ae83 + afumott + 
                    aescmae + arendtot + birthday, 
                  data = complete_dat)
fulldatamod_res <- round(c(
  Estimate = summary(fulldatamod)$coefficients[2, 1],
  StdError = summary(fulldatamod)$coefficients[2, 2],
  confint(fulldatamod)[2, ]
), 3)

# combine results into table
tab_bfdur_diff_m <- rbind(
  fulldatamod_res,
  imp_m1_m50_res,
  imp_m1_m75_res,
  imp_m1_m100_res)
tab_bfdur_diff_m <- round(tab_bfdur_diff_m,3)
tab_bfdur_diff_m <- cbind(c("Complete data", 50, 75, 100),
                          tab_bfdur_diff_m)
colnames(tab_bfdur_diff_m) <- c("m", "Est", "SE", "95%CI_L", "95%CI_U")
rownames(tab_bfdur_diff_m) <- NULL
kable(tab_bfdur_diff_m)

# test - relationship between externalising and smoking (controlling for int)
imp_m1_m50 <- with(bygroup_50,
                   glm(y18_smoking ~ y15_externalising + y15_internalising + y11_cte_recode + 
                        asexo + aethnicity + ae83 + afumott + 
                        aescmae + arendtot + birthday, family=binomial))
imp_m1_m50_res <- summary(pool(imp_m1_m50),
                          conf.int=TRUE)[2,c(2,3,7,8)] 

imp_m1_m75 <- with(bygroup_75,
                   glm(y18_smoking ~ y15_externalising + y15_internalising + y11_cte_recode + 
                        asexo + aethnicity + ae83 + afumott + 
                        aescmae + arendtot + birthday, family=binomial))
imp_m1_m75_res <- summary(pool(imp_m1_m75),
                          conf.int=TRUE)[2,c(2,3,7,8)]

imp_m1_m100 <- with(bygroup_100,
                    glm(y18_smoking ~ y15_externalising + y15_internalising + y11_cte_recode + 
                         asexo + aethnicity + ae83 + afumott + 
                         aescmae + arendtot + birthday, family=binomial))
imp_m1_m100_res <- summary(pool(imp_m1_m100),
                           conf.int=TRUE)[2,c(2,3,7,8)]

# obtain the estimate in the complete dat data, when no values are missing. 
fulldatamod <- glm(y18_smoking ~ y15_externalising + y15_internalising + y11_cte_recode + 
                    asexo + aethnicity + ae83 + afumott + 
                    aescmae + arendtot + birthday, 
                   data = complete_dat, family=binomial)
fulldatamod_res <- round(c(
  Estimate = summary(fulldatamod)$coefficients[2, 1],
  StdError = summary(fulldatamod)$coefficients[2, 2],
  confint(fulldatamod)[2, ]
), 3)

# combine results into table
tab_bfdur_diff_m <- rbind(
  fulldatamod_res,
  imp_m1_m50_res,
  imp_m1_m75_res,
  imp_m1_m100_res)
tab_bfdur_diff_m <- round(tab_bfdur_diff_m,3)
tab_bfdur_diff_m <- cbind(c("Complete data", 50, 75, 100),
                          tab_bfdur_diff_m)
colnames(tab_bfdur_diff_m) <- c("m", "Est","SE", "95%CI_L", "95%CI_U")
rownames(tab_bfdur_diff_m) <- NULL
kable(tab_bfdur_diff_m)

# test - relationship between internalising and drug use (controlling for ext)
imp_m1_m50 <- with(bygroup_50,
                   glm(y18_druguse_recode ~ y15_internalising + y15_externalising + y11_cte_recode + 
                         asexo + aethnicity + ae83 + afumott + 
                         aescmae + arendtot + birthday, family=binomial))
imp_m1_m50_res <- summary(pool(imp_m1_m50),
                          conf.int=TRUE)[2,c(2,3,7,8)]

imp_m1_m75 <- with(bygroup_75,
                   glm(y18_druguse_recode ~ y15_internalising + y15_externalising + y11_cte_recode + 
                         asexo + aethnicity + ae83 + afumott + 
                         aescmae + arendtot + birthday, family=binomial))
imp_m1_m75_res <- summary(pool(imp_m1_m75),
                          conf.int=TRUE)[2,c(2,3,7,8)]

imp_m1_m100 <- with(bygroup_100,
                    glm(y18_druguse_recode ~ y15_internalising + y15_externalising + y11_cte_recode + 
                          asexo + aethnicity + ae83 + afumott + 
                          aescmae + arendtot + birthday, family=binomial))
imp_m1_m100_res <- summary(pool(imp_m1_m100),
                           conf.int=TRUE)[2,c(2,3,7,8)]

# obtain the estimate in the complete dat data, when no values are missing. 
fulldatamod <- glm(y18_druguse_recode ~ y15_internalising + y15_externalising + y11_cte_recode + 
                     asexo + aethnicity + ae83 + afumott + 
                     aescmae + arendtot + birthday, 
                   data = complete_dat, family=binomial)
fulldatamod_res <- round(c(
  Estimate = summary(fulldatamod)$coefficients[2, 1],
  StdError = summary(fulldatamod)$coefficients[2, 2],
  confint(fulldatamod)[2, ]
), 3)

# combine results into table
tab_bfdur_diff_m <- rbind(
  fulldatamod_res,
  imp_m1_m50_res,
  imp_m1_m75_res,
  imp_m1_m100_res)
tab_bfdur_diff_m <- round(tab_bfdur_diff_m,3)
tab_bfdur_diff_m <- cbind(c("Complete data", 50, 75, 100),
                          tab_bfdur_diff_m)
colnames(tab_bfdur_diff_m) <- c("m", "Est", "SE", "95%CI_L", "95%CI_U")
rownames(tab_bfdur_diff_m) <- NULL
kable(tab_bfdur_diff_m)


### cross-tabs ####
# density plot for aescmae (maternal education) was worrying - check imputed v original data
# Loop through each imputed dataset
for (i in 1:bygroup_75$m) {
  cat("\nImputed dataset:", i, "\n")
  dat_i <- complete(bygroup_75, i)  # extract the ith dataset
  print(table(dat_i$aescmae, useNA = "ifany"))
}

# Collect proportions across imputations
prop_list <- lapply(1:bygroup_75$m, function(i) {
  dat_i <- complete(bygroup_75, i)
  prop.table(table(dat_i$aescmae))
})
# Combine into one data frame
props_combined <- do.call(rbind, prop_list)
mean_props <- colMeans(props_combined)
print(mean_props)
# Create a long-format data frame with the imputation number included
long_data <- complete(bygroup_75, action = "long", include = TRUE)
# Check 
head(long_data)
# Tabulate aescmae by .imp 
table(long_data$aescmae, long_data$.imp, useNA = "ifany")

### y15_internalising --- no values of 17 in orig or imputation - check the rest of the 
# imputation here (note, also compared to a norm imputed version instead of pmm - no difference in results)
# Loop through each imputed dataset
for (i in 1:bygroup_75$m) {
  cat("\nImputed dataset:", i, "\n")
  dat_i <- complete(bygroup_75, i)  # extract the ith dataset
  print(table(dat_i$y15_internalising, useNA = "ifany"))
}
# Collect proportions across imputations
prop_list <- lapply(1:bygroup_75$m, function(i) {
  dat_i <- complete(bygroup_75, i)
  prop.table(table(dat_i$y15_internalising))
})
# Combine into one data frame
props_combined <- do.call(rbind, prop_list)
mean_props <- colMeans(props_combined)
print(mean_props)
# Create a long-format data frame with the imputation number included
long_data <- complete(bygroup_75, action = "long", include = TRUE)
head(long_data)
# Tabulate y15_internalising by .imp (which is like Stata's _mj)
table(long_data$y15_internalising, long_data$.imp, useNA = "ifany")


## Examine prevalence/mean differences in observed v imputed ####
# determine the difference in log-odds resulting from the current MAR-based imputation of (y18_smoking/y18_druguse)
# between observed cases and those which are imputed within each dataset

# number of imputed datasets
numimps <- 75

## smoking
MI_imp75 <- complete(bygroup_75, "broad", inc=TRUE)
MI_imp75$incomp <- as.integer(ici(MI_imp75$y18_smoking.0)) 
MI_imp75_difflogodds <- 1:numimps
for (i in 1:numimps) {
  tempvar <- paste0("y18_smoking.", i)
  x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = MI_imp75)
  print(x$coefficients)
  MI_imp75_difflogodds[i] <- x$coefficients[2]
}

mean(MI_imp75_difflogodds)
exp(mean(MI_imp75_difflogodds))


## drug use
MI_imp75$incomp <- as.integer(ici(MI_imp75$y18_druguse_recode.0)) 
MI_imp75_difflogodds <- 1:numimps
for (i in 1:numimps) {
  tempvar <- paste0("y18_druguse_recode.", i)
  x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = MI_imp75)
  print(x$coefficients)
  MI_imp75_difflogodds[i] <- x$coefficients[2]
}

mean(MI_imp75_difflogodds)
exp(mean(MI_imp75_difflogodds))


## test alcohol
MI_imp75$incomp <- as.integer(ici(MI_imp75$y18_alcohol.0))   
MI_imp75_difflogodds <- 1:numimps
for (i in 1:numimps) {
  tempvar <- paste0("y18_alcohol.", i)
  x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = MI_imp75)
  print(x$coefficients)
  MI_imp75_difflogodds[i] <- x$coefficients[2]
}

mean(MI_imp75_difflogodds)
exp(mean(MI_imp75_difflogodds))


### Tipping Point Analysis ####
# further checks on smoking / drug use assumptions

# Impute a variable and add on (e.g.) 2 points/units to the imputed data (add 5, 
# add 10), and see how sensitive the conclusions are to the imputation
# if we were off by 2/5/10, would this impact our conclusions?
# In other words, a range of values from one extreme where responders have 4 times 
# the odds of the outcome compared to non-responders, to the other extreme where 
# non-responders have 4 times the odds of the outcome compared to responders.
# Here we are looping through a range of delta in increments of 0.25 units. The code carries out 
# the imputation and fits the regression model of interests before saving the specific parameter 
# relating to the confounder-adjusted association. Further lines fit additional models 
# which derive some Marginal Sensitivity Parameters to aid the interpretation of the Tipping Point 
# results. 

# NOTE, MSP: marginal sensitivity parameter (mean difference between responders and non-responders)
# NOTE, IMOR: the ignorable missingness odds-ratio which describes the odd of (the outcome) for non-responders compared to responders.

# prep
install.packages("devtools", dependencies = TRUE)
library(devtools)
install_github("moreno-betancur/mice")
library(mice)

###### Smoking -----------------------------------------------------------------
# Create MI data frame -- removing most age 15 (missingness)
MI_df <- dat[,c("y11_cte_recode","asexo","aethnicity","ae83","afumott","aescmae",
                "arendtot","birthday","y15_internalising","y15_externalising","y18_alcohol",
                "y18_smoking","y18_druguse_recode","joverall_pa","jpercultraprocessadoscal","jmvpab5",
                "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                "y11_alcohol","y11_smoking","goverall_pa","gpercultraprocessados",
                "gmvpab5","aa07","y6_ctspc","y11_ctspc")]
str(MI_df)
head(MI_df)

MI_df_mnar <- MI_df 
MI_df_mnar$M_y18_smoking <- as.integer(as.character(MI_df_mnar$y18_smoking))
str(MI_df_mnar$M_y18_smoking)
names(MI_df_mnar)

# set up a new prediction matrix for the new imputation
# adding missing indicator for smoking, and letting this new indicator predict any 
# partially observed variables in the analysis model for smoking
predMatrix_1<-matrix(NA, nrow=30 ,ncol=30)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df_mnar)
predMatrix_1

predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1) # missing smoking as predictor
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aethnicity",]        <-c(1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aescmae",]           <-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["y15_internalising",] <-c(1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["y15_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["y18_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0)
predMatrix_1["y18_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0)
predMatrix_1["y18_druguse_recode",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0) 
predMatrix_1["joverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)  
predMatrix_1["jpercultraprocessadoscal",]         <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 
predMatrix_1["jmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0) 
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0) 
predMatrix_1["goverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0)
predMatrix_1["gpercultraprocessados",]             <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["gmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["aa07",]              <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)   
predMatrix_1["y6_ctspc",]          <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
predMatrix_1["y11_ctspc",]         <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
predMatrix_1["M_y18_smoking",]     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predMatrix_1    

# Set-up predictor matrix for unidentifiable part 
# in this case the whole matrix is zeroes because the unidentifiable part 
# of the imputation model contains a single constant (delta*M) rather than 
# additional contributions from the other variables in the dataset
predSens<-matrix(NA, nrow=30 ,ncol=30)
colnames(predSens)<-paste(":",names(MI_df_mnar[1:30]),sep="")
rownames(predSens)<-names(MI_df_mnar[1:30])
predSens
predSens["y11_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["asexo",]                       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aethnicity",]                  <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["ae83",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["afumott",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aescmae",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["arendtot",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["birthday",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_druguse_recode",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["joverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jpercultraprocessadoscal",]    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["goverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gpercultraprocessados",]       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aa07",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y6_ctspc",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_ctspc",]                   <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["M_y18_smoking",]               <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens

#Set-up list with sensitivity parameter values
pSens<-rep(list(list("")), ncol(MI_df_mnar[1:30]))
names(pSens)<-names(MI_df_mnar[1:30])
pSens

# set up vector describing manner of imputation for each variable
narfcs_method <- c("pmm","","logreg","","",
                   "pmm","","","pmm","pmm",
                   "logreg","logregSens","logreg","pmm","pmm",
                   "pmm","pmm","pmm","pmm","pmm",
                   "pmm","logreg","logreg","pmm","pmm",
                   "pmm","pmm","pmm","pmm","")

# choose a smaller number of datasets given computational burden
narfcs_numimps <- 25

# used to collect the parameters of interest
k <-0
tipping <- as.data.frame(array(dim=c(dim=17,5)))
colnames(tipping) <- c("csp", "msp", "beta", "se_beta","sampprev")

# Looping over delta values (j)
for (j in seq.int(-2, 2, by = 0.25)) {
  k <- k+1 
  cat("Running delta j =", j, " (iteration", k, ")\n") ### trying to debug
  
  # specify a delta value for the prediction equation for (smoking)
  pSens[["y18_smoking"]]<-list(c(j))
  
  # NARFCS imputation
  y18_smoking_NARFCS<-mice(MI_df_mnar[1:30],
                           m=narfcs_numimps, method=narfcs_method, 
                           predictorMatrix=predMatrix_1, predictorSens=predSens, 
                           parmSens=pSens, seed=260873, print=F, maxit=25)
  
  # create a dataset which can be used to determine the MSP 
  # the imputed data that MICE creates is stored much more efficiently than -ice- so requires some processing
  smoking_NARFCS <- complete(y18_smoking_NARFCS, "broad", inc=TRUE)
  
  # create a flag to indicate when (smoking) was missing 
  smoking_NARFCS$incomp <- as.integer(ici(smoking_NARFCS$y18_smoking.0))
  
  # scroll through imputed variables and derive a succession of MSPs
  smoking_NARFCS_difflogodds <- 1:narfcs_numimps
  for (i in 1:narfcs_numimps) {
    tempvar <- paste0("y18_smoking.", i)
    x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = smoking_NARFCS)
    smoking_NARFCS_difflogodds[i] <- x$coefficients[2]
  }
  
  # derive the average MSP across imputed datasets for the given delta value
  # and output along with delta and the newly estimated regression parameter
  newest <- round(summary(pool(with(y18_smoking_NARFCS, glm(formula = y18_smoking ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial(link = "logit"))))), 3)
  
  # derive the prevalence of Y in the whole population
  wholesampprev <- round(summary(pool(with(y18_smoking_NARFCS, glm(formula = y18_smoking ~ 1, family = binomial(link = "identity"))))), 3)
  
  tipping[k,1] <- j
  tipping[k,2] <- mean(smoking_NARFCS_difflogodds)
  tipping[k,3] <- newest[2,1]
  tipping[k,4] <- newest[2,2]
  tipping[k,5] <- wholesampprev[1,1]  
  
  tipping$lower <- tipping$beta - 1.96*tipping$se_beta
  tipping$upper <- tipping$beta + 1.96*tipping$se_beta
  tipping$imor <- exp(tipping$msp)
  
  print(tipping[k,])
}

names(tipping)
head(tipping)

# Convert all beta coefficients and their CIs to odds ratios
tipping$or_beta <- exp(tipping$beta)
tipping$or_lower <- exp(tipping$lower)
tipping$or_upper <- exp(tipping$upper)
print(tipping)

getwd()
dir.create("TippingPoint")
write.csv(tipping, "mice_090725/TippingPoint/tipping_results_CTE_smoking.csv", row.names = FALSE)
tipping <- read.csv("mice_090725/TippingPoint/tipping_results_CTE_smoking.csv")

### plotting
library(ggplot2)

# OR / delta
p <- ggplot(tipping, aes(x = csp, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_smoking",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_OR_delta_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# OR / imor
p <- ggplot(tipping, aes(x = imor, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_smoking",
    x = expression("imor"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_OR_imor_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# OR / sampprev
p <- ggplot(tipping, aes(x = sampprev, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_smoking",
    x = expression("sampprev"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_OR_sampprev_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)

#### log the y-axis (and x-axis, for imor) 
tipping$log_or <- log(tipping$or_beta)
tipping$log_or_lower <- log(tipping$or_lower)
tipping$log_or_upper <- log(tipping$or_upper)
tipping$log_imor <- log(tipping$imor) 
print(tipping)

# log(OR) / delta
p <- ggplot(tipping, aes(x = csp, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_smoking",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_log_OR_delta_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# log(OR) / log(imor)
p <- ggplot(tipping, aes(x = log_imor, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_smoking",
    x = expression("log_imor"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_log_OR_imor_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)


##### Drug use -----------------------------------------------------------------
MI_df_mnar <- MI_df 
MI_df_mnar$M_y18_druguse <- as.integer(ici(MI_df_mnar$y18_druguse_recode))
str(MI_df_mnar$M_y18_druguse)
names(MI_df_mnar)

# set up a new prediction matrix for the new imputation
predMatrix_1<-matrix(NA, nrow=30 ,ncol=30)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df_mnar)
predMatrix_1
predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1) # missing drug use as predictor
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aethnicity",]        <-c(1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aescmae",]           <-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["y15_internalising",] <-c(1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["y15_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["y18_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0)
predMatrix_1["y18_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0)
predMatrix_1["y18_druguse_recode",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0) 
predMatrix_1["joverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)  
predMatrix_1["jpercultraprocessadoscal",]         <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 
predMatrix_1["jmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0) 
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0) 
predMatrix_1["goverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0)
predMatrix_1["gpercultraprocessados",]             <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["gmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["aa07",]              <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)   
predMatrix_1["y6_ctspc",]          <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
predMatrix_1["y11_ctspc",]         <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
predMatrix_1["M_y18_druguse",]     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predMatrix_1    

# Set-up predictor matrix for unidentifiable part 
predSens<-matrix(NA, nrow=30 ,ncol=30)
colnames(predSens)<-paste(":",names(MI_df_mnar[1:30]),sep="")
rownames(predSens)<-names(MI_df_mnar[1:30])
predSens
predSens["y11_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["asexo",]                       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aethnicity",]                  <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["ae83",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["afumott",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aescmae",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["arendtot",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["birthday",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_druguse_recode",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["joverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jpercultraprocessadoscal",]    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["goverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gpercultraprocessados",]       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aa07",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y6_ctspc",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_ctspc",]                   <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["M_y18_druguse",]               <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens

#Set-up list with sensitivity parameter values
pSens<-rep(list(list("")), ncol(MI_df_mnar[1:30]))
names(pSens)<-names(MI_df_mnar[1:30])
pSens

# set up vector describing manner of imputation for each variable
narfcs_method <- c("pmm","","logreg","","",
                   "pmm","","","pmm","pmm",
                   "logreg","logreg","logregSens","pmm","pmm",
                   "pmm","pmm","pmm","pmm","pmm",
                   "pmm","logreg","logreg","pmm","pmm",
                   "pmm","pmm","pmm","pmm","")

# choose a smaller number of datasets given computational burden
narfcs_numimps <- 25

# used to collect the parameters of interest
k <-0
tipping <- as.data.frame(array(dim=c(dim=17,5)))
colnames(tipping) <- c("csp", "msp", "beta", "se_beta","sampprev")

# Looping over delta values (j)
for (j in seq.int(-2, 2, by = 0.25)) {
  k <- k+1 
  cat("Running delta j =", j, " (iteration", k, ")\n") ### trying to debug
  
  # specify a delta value for the prediction equation for (drug use)
  pSens[["y18_druguse_recode"]]<-list(c(j))
  
  # NARFCS imputation
  y18_druguse_NARFCS<-mice(MI_df_mnar[1:30],
                           m=narfcs_numimps, method=narfcs_method, 
                           predictorMatrix=predMatrix_1, predictorSens=predSens, 
                           parmSens=pSens, seed=260873, print=F, maxit=25)
  
  # create a dataset which can be used to determine the MSP 
  # the imputed data that MICE creates is stored much more efficiently than -ice- so requires some processing
  druguse_NARFCS <- complete(y18_druguse_NARFCS, "broad", inc=TRUE)
  
  # create a flag to indicate when (drug use) was missing 
  druguse_NARFCS$incomp <- as.integer(ici(druguse_NARFCS$y18_druguse_recode.0))
  
  # scroll through imputed variables and derive a succession of MSPs
  druguse_NARFCS_difflogodds <- 1:narfcs_numimps
  for (i in 1:narfcs_numimps) {
    tempvar <- paste0("y18_druguse_recode.", i)
    x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = druguse_NARFCS)
    druguse_NARFCS_difflogodds[i] <- x$coefficients[2]
  }
  
  # derive the average MSP across imputed datasets for the given delta value
  # and output along with delta and the newly estimated regression parameter
  newest <- round(summary(pool(with(y18_druguse_NARFCS, glm(formula = y18_druguse_recode ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial(link = "logit"))))), 3)
  
  # derive the prevalence of Y in the whole population
  wholesampprev <- round(summary(pool(with(y18_druguse_NARFCS, glm(formula = y18_druguse_recode ~ 1, family = binomial(link = "identity"))))), 3)
  
  tipping[k,1] <- j
  tipping[k,2] <- mean(druguse_NARFCS_difflogodds)
  tipping[k,3] <- newest[2,1]
  tipping[k,4] <- newest[2,2]
  tipping[k,5] <- wholesampprev[1,1]  
  
  tipping$lower <- tipping$beta - 1.96*tipping$se_beta
  tipping$upper <- tipping$beta + 1.96*tipping$se_beta
  tipping$imor <- exp(tipping$msp)
  
  print(tipping[k,])
}

names(tipping)

# Convert all beta coefficients and their CIs to odds ratios
tipping$or_beta <- exp(tipping$beta)
tipping$or_lower <- exp(tipping$lower)
tipping$or_upper <- exp(tipping$upper)
print(tipping)
getwd()
write.csv(tipping, "TippingPoint/tipping_results_CTE_druguse.csv", row.names = FALSE)

#### log the y-axis (and x-axis, for imor) 
tipping$log_or <- log(tipping$or_beta)
tipping$log_or_lower <- log(tipping$or_lower)
tipping$log_or_upper <- log(tipping$or_upper)
tipping$log_imor <- log(tipping$imor) 
print(tipping)

# log(OR) / delta
p <- ggplot(tipping, aes(x = csp, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_druguse",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_log_OR_delta_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# log(OR) / log(imor)
p <- ggplot(tipping, aes(x = log_imor, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y11_cte_recode Effect on y18_druguse",
    x = expression("log_imor"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_log_OR_imor_CTE.png",
  plot = p,width = 7,height = 5,dpi = 300
)



###### Smoking and externalising -----------------------------------------------
MI_df_mnar <- MI_df 
MI_df_mnar$M_y18_smoking <- as.integer(ici(MI_df_mnar$y18_smoking))
names(MI_df_mnar)

# set up a new prediction matrix for the new imputation
predMatrix_1<-matrix(NA, nrow=30 ,ncol=30)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df_mnar)
predMatrix_1
predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1) # missing smoking as predictor
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aethnicity",]        <-c(1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aescmae",]           <-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["y15_internalising",] <-c(1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["y15_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing smoking as predictor
predMatrix_1["y18_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0)
predMatrix_1["y18_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0)
predMatrix_1["y18_druguse_recode",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0) 
predMatrix_1["joverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)  
predMatrix_1["jpercultraprocessadoscal",]         <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 
predMatrix_1["jmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0) 
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0) 
predMatrix_1["goverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0)
predMatrix_1["gpercultraprocessados",]             <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["gmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["aa07",]              <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)   
predMatrix_1["y6_ctspc",]          <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
predMatrix_1["y11_ctspc",]         <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
predMatrix_1["M_y18_smoking",]     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predMatrix_1    

# Set-up predictor matrix for unidentifiable part 
predSens<-matrix(NA, nrow=30 ,ncol=30)
colnames(predSens)<-paste(":",names(MI_df_mnar[1:30]),sep="")
rownames(predSens)<-names(MI_df_mnar[1:30])
predSens
predSens["y11_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["asexo",]                       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aethnicity",]                  <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["ae83",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["afumott",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aescmae",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["arendtot",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["birthday",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_druguse_recode",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["joverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jpercultraprocessadoscal",]    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["goverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gpercultraprocessados",]       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aa07",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y6_ctspc",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_ctspc",]                   <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["M_y18_smoking",]               <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens

#Set-up list with sensitivity parameter values
pSens<-rep(list(list("")), ncol(MI_df_mnar[1:30]))
names(pSens)<-names(MI_df_mnar[1:30])
pSens

# set up vector describing manner of imputation for each variable
narfcs_method <- c("pmm","","logreg","","",
                   "pmm","","","pmm","pmm",
                   "logreg","logregSens","logreg","pmm","pmm",
                   "pmm","pmm","pmm","pmm","pmm",
                   "pmm","logreg","logreg","pmm","pmm",
                   "pmm","pmm","pmm","pmm","")

# choose a smaller number of datasets given computational burden
narfcs_numimps <- 25

# used to collect the parameters of interest
k <-0
tipping <- as.data.frame(array(dim=c(dim=17,5)))
colnames(tipping) <- c("csp", "msp", "beta", "se_beta","sampprev")

# Looping over delta values (j)
for (j in seq.int(-2, 2, by = 0.25)) {
  k <- k+1 
  cat("Running delta j =", j, " (iteration", k, ")\n") ### trying to debug
  
  # specify a delta value for the prediction equation for (smoking)
  pSens[["y18_smoking"]]<-list(c(j))
  
  # NARFCS imputation
  y18_smoking_NARFCS<-mice(MI_df_mnar[1:30],
                           m=narfcs_numimps, method=narfcs_method, 
                           predictorMatrix=predMatrix_1, predictorSens=predSens, 
                           parmSens=pSens, seed=260873, print=F, maxit=25)
  
  # create a dataset which can be used to determine the MSP 
  # the imputed data that MICE creates is stored much more efficiently than -ice- so requires some processing
  smoking_NARFCS <- complete(y18_smoking_NARFCS, "broad", inc=TRUE)
  
  # create a flag to indicate when (smoking) was missing 
  smoking_NARFCS$incomp <- as.integer(ici(smoking_NARFCS$y18_smoking.0))
  
  # scroll through imputed variables and derive a succession of MSPs
  smoking_NARFCS_difflogodds <- 1:narfcs_numimps
  for (i in 1:narfcs_numimps) {
    tempvar <- paste0("y18_smoking.", i)
    x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = smoking_NARFCS)
    smoking_NARFCS_difflogodds[i] <- x$coefficients[2]
  }
  
  # derive the average MSP across imputed datasets for the given delta value
  # and output along with delta and the newly estimated regression parameter
  newest <- round(summary(pool(with(y18_smoking_NARFCS, glm(formula = y18_smoking ~ y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial(link = "logit"))))), 3)
  
  # derive the prevalence of Y in the whole population
  wholesampprev <- round(summary(pool(with(y18_smoking_NARFCS, glm(formula = y18_smoking ~ 1, family = binomial(link = "identity"))))), 3)
  
  tipping[k,1] <- j
  tipping[k,2] <- mean(smoking_NARFCS_difflogodds)
  tipping[k,3] <- newest[2,1]
  tipping[k,4] <- newest[2,2]
  tipping[k,5] <- wholesampprev[1,1]  
  
  tipping$lower <- tipping$beta - 1.96*tipping$se_beta
  tipping$upper <- tipping$beta + 1.96*tipping$se_beta
  tipping$imor <- exp(tipping$msp)
  
  print(tipping[k,])
}
names(tipping)
write.csv(tipping, "mice_260126/TippingPoint/tipping_results_externalising_smoking.csv", row.names = FALSE)
tipping <- read.csv("mice_260126/TippingPoint/tipping_results_externalising_smoking.csv")


### plot 
library(ggplot2)
tipping$or_beta <- exp(tipping$beta)
tipping$or_lower <- exp(tipping$lower)
tipping$or_upper <- exp(tipping$upper)
print(tipping)

# OR / delta
p <- ggplot(tipping, aes(x = csp, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_smoking",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_OR_delta_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# OR / imor
p <- ggplot(tipping, aes(x = imor, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_smoking",
    x = expression("imor"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_OR_imor_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# OR / sampprev
p <- ggplot(tipping, aes(x = sampprev, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_smoking",
    x = expression("sampprev"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_OR_sampprev_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

#### log the y-axis (and x-axis, for imor) 
tipping$log_or <- log(tipping$or_beta)
tipping$log_or_lower <- log(tipping$or_lower)
tipping$log_or_upper <- log(tipping$or_upper)
tipping$log_imor <- log(tipping$imor) 
print(tipping)

# log(OR) / delta
p <- ggplot(tipping, aes(x = csp, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_smoking",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_log_OR_delta_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# log(OR) / log(imor)
p <- ggplot(tipping, aes(x = log_imor, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_smoking",
    x = expression("log_imor"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/smoking_log_OR_imor_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)


###### Drug use and externalising ----------------------------------------------
MI_df_mnar <- MI_df
MI_df_mnar$M_y18_druguse <- as.integer(ici(MI_df_mnar$y18_druguse_recode))
names(MI_df_mnar)

# set up a new prediction matrix for the new imputation
predMatrix_1<-matrix(NA, nrow=30 ,ncol=30)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df_mnar)
predMatrix_1
predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1) # missing drug use as predictor
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aethnicity",]        <-c(1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["aescmae",]           <-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["y15_internalising",] <-c(1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["y15_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1) # missing drug use as predictor
predMatrix_1["y18_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0)
predMatrix_1["y18_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0)
predMatrix_1["y18_druguse_recode",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0) 
predMatrix_1["joverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0)  
predMatrix_1["jpercultraprocessadoscal",]         <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 
predMatrix_1["jmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0) 
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_internalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y18_externalising",] <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["y11_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,0) 
predMatrix_1["y11_smoking",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0) 
predMatrix_1["goverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0)
predMatrix_1["gpercultraprocessados",]             <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["gmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0) 
predMatrix_1["aa07",]              <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)   
predMatrix_1["y6_ctspc",]          <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0)
predMatrix_1["y11_ctspc",]         <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
predMatrix_1["M_y18_druguse",]     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predMatrix_1    

# Set-up predictor matrix for unidentifiable part 
predSens<-matrix(NA, nrow=30 ,ncol=30)
colnames(predSens)<-paste(":",names(MI_df_mnar[1:30]),sep="")
rownames(predSens)<-names(MI_df_mnar[1:30])
predSens
predSens["y11_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["asexo",]                       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aethnicity",]                  <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["ae83",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["afumott",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aescmae",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["arendtot",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["birthday",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_druguse_recode",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["joverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jpercultraprocessadoscal",]    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["jmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y15_cte_recode",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_internalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y18_externalising",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_alcohol",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_smoking",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["goverall_pa",]                 <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gpercultraprocessados",]       <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["gmvpab5",]                     <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["aa07",]                        <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y6_ctspc",]                    <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["y11_ctspc",]                   <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens["M_y18_druguse",]               <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
predSens

#Set-up list with sensitivity parameter values
pSens<-rep(list(list("")), ncol(MI_df_mnar[1:30]))
names(pSens)<-names(MI_df_mnar[1:30])
pSens

# set up vector describing manner of imputation for each variable
narfcs_method <- c("pmm","","logreg","","",
                   "pmm","","","pmm","pmm",
                   "logreg","logreg","logregSens","pmm","pmm",
                   "pmm","pmm","pmm","pmm","pmm",
                   "pmm","logreg","logreg","pmm","pmm",
                   "pmm","pmm","pmm","pmm","")

# choose a smaller number of datasets given computational burden
narfcs_numimps <- 25

# used to collect the parameters of interest
k <-0
tipping <- as.data.frame(array(dim=c(dim=17,5)))
colnames(tipping) <- c("csp", "msp", "beta", "se_beta","sampprev")

# Looping over delta values (j)
for (j in seq.int(-2, 2, by = 0.25)) {
  k <- k+1 
  cat("Running delta j =", j, " (iteration", k, ")\n") ### trying to debug
  
  # specify a delta value for the prediction equation for (drug use)
  pSens[["y18_druguse_recode"]]<-list(c(j))
  
  # NARFCS imputation
  y18_druguse_NARFCS<-mice(MI_df_mnar[1:30],
                           m=narfcs_numimps, method=narfcs_method, 
                           predictorMatrix=predMatrix_1, predictorSens=predSens, 
                           parmSens=pSens, seed=260873, print=F, maxit=25)
  
  # create a dataset which can be used to determine the MSP 
  # the imputed data that MICE creates is stored much more efficiently than -ice- so requires some processing
  druguse_NARFCS <- complete(y18_druguse_NARFCS, "broad", inc=TRUE)
  
  # create a flag to indicate when (drug use) was missing 
  druguse_NARFCS$incomp <- as.integer(ici(druguse_NARFCS$y18_druguse_recode.0))
  
  # scroll through imputed variables and derive a succession of MSPs
  druguse_NARFCS_difflogodds <- 1:narfcs_numimps
  for (i in 1:narfcs_numimps) {
    tempvar <- paste0("y18_druguse_recode.", i)
    x <- glm(formula =get(tempvar) ~ incomp, family = "binomial", data = druguse_NARFCS)
    druguse_NARFCS_difflogodds[i] <- x$coefficients[2]
  }
  
  # derive the average MSP across imputed datasets for the given delta value
  # and output along with delta and the newly estimated regression parameter
  newest <- round(summary(pool(with(y18_druguse_NARFCS, glm(formula = y18_druguse_recode ~ y15_externalising + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial(link = "logit"))))), 3)
  
  # derive the prevalence of Y in the whole population
  wholesampprev <- round(summary(pool(with(y18_druguse_NARFCS, glm(formula = y18_druguse_recode ~ 1, family = binomial(link = "identity"))))), 3)
  
  tipping[k,1] <- j
  tipping[k,2] <- mean(druguse_NARFCS_difflogodds)
  tipping[k,3] <- newest[2,1]
  tipping[k,4] <- newest[2,2]
  tipping[k,5] <- wholesampprev[1,1]  
  
  tipping$lower <- tipping$beta - 1.96*tipping$se_beta
  tipping$upper <- tipping$beta + 1.96*tipping$se_beta
  tipping$imor <- exp(tipping$msp)
  
  print(tipping[k,])
}
names(tipping)
write.csv(tipping, "mice_260126/TippingPoint/tipping_results_externalising_druguse.csv", row.names = FALSE)

tipping$or_beta <- exp(tipping$beta)
tipping$or_lower <- exp(tipping$lower)
tipping$or_upper <- exp(tipping$upper)
print(tipping)

# OR / delta
p <- ggplot(tipping, aes(x = csp, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_druguse",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_OR_delta_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# OR / imor
p <- ggplot(tipping, aes(x = imor, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_druguse",
    x = expression("imor"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_OR_imor_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# OR / sampprev
p <- ggplot(tipping, aes(x = sampprev, y = or_beta)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = or_lower, ymax = or_upper), width = 0.1, color = "steelblue") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_druguse",
    x = expression("sampprev"),
    y = expression("Estimated Coefficient (OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_OR_sampprev_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

#### log the y-axis (and x-axis, for imor) 
tipping$log_or <- log(tipping$or_beta)
tipping$log_or_lower <- log(tipping$or_lower)
tipping$log_or_upper <- log(tipping$or_upper)
tipping$log_imor <- log(tipping$imor)
print(tipping)

# log(OR) / delta
p <- ggplot(tipping, aes(x = csp, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_druguse",
    x = expression(Delta ~ "(Sensitivity Parameter)"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_log_OR_delta_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)

# log(OR) / log(imor)
p <- ggplot(tipping, aes(x = log_imor, y = log_or)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = log_or_lower, ymax = log_or_upper), width = 0.1, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(
    title = "Tipping Point Analysis for y15_externalising Effect on y18_druguse",
    x = expression("log_imor"),
    y = expression("Estimated Coefficient log(OR)")
  ) +
  theme_minimal(base_size = 14)
ggsave(
  filename = "TippingPoint/druguse_log_OR_imor_ext.png",
  plot = p,width = 7,height = 5,dpi = 300
)



## Monte Carlo error ####
# function from https://thestatsgeek.com/2022/02/23/how-many-imputations-with-mice-assessing-monte-carlo-error-after-multiple-imputation-in-r/
# Check whether the Monte Carlo error of B is approximately 10 per cent of its standard error, and whether 75 imputed datasets are sufficient

# function:
miceMCError <- function(pooledRes) {
  monteCarloSE <- sqrt(pooledRes$pooled$b/pooledRes$m)
  ciLower <- pooledRes$pooled$estimate - qt(0.975,df=pooledRes$m-1)*monteCarloSE
  ciUpper <- pooledRes$pooled$estimate + qt(0.975,df=pooledRes$m-1)*monteCarloSE
  mcTable <- cbind(pooledRes$pooled$estimate, monteCarloSE, ciLower, ciUpper)
  colnames(mcTable) <- c("Estimate", "Monte Carlo SE", "95% CI lower limit", "95% CI upper limit")
  print(mcTable)
  print("Warning: 95% CI only quantifies Monte-Carlo uncertainty!")
}

## unadjusted models
# smoking
bygroup_est_y18_smoking <- with(bygroup_75, glm(y18_smoking ~ y11_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_smoking))
summary(pool(bygroup_est_y18_smoking), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_smoking)
miceMCError(pooled)

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup_75, glm(y18_alcohol ~ y11_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
summary(pool(bygroup_est_y18_alcohol), conf.int = TRUE, exponentiate = TRUE) # this provides ORs 
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)

# drug use
bygroup_est_y18_druguse <- with(bygroup_75, glm(y18_druguse_recode ~ y11_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_druguse))
summary(pool(bygroup_est_y18_druguse), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_druguse)
miceMCError(pooled)

# overall PA
bygroup_est_joverall_pa <- with(bygroup_75, lm(joverall_pa ~ y11_cte_recode))
summary(pool(bygroup_est_joverall_pa))
summary(pool(bygroup_est_joverall_pa), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_joverall_pa)
miceMCError(pooled)

# bouted PA
bygroup_est_jmvpab5 <- with(bygroup_75, lm(jmvpab5 ~ y11_cte_recode))
summary(pool(bygroup_est_jmvpab5))
summary(pool(bygroup_est_jmvpab5), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_jmvpab5)
miceMCError(pooled)

# UPF
bygroup_est_jpercultraprocessadoscal <- with(bygroup_75, lm(jpercultraprocessadoscal ~ y11_cte_recode))
summary(pool(bygroup_est_jpercultraprocessadoscal))
summary(pool(bygroup_est_jpercultraprocessadoscal), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_jpercultraprocessadoscal)
miceMCError(pooled)
#### all monte carlo errors are less than 10% of the standard errors

## adjusted
# smoking
bygroup_est_y18_smoking <- with(bygroup_75, glm(y18_smoking ~ y11_cte_recode + asexo + 
                                                  aethnicity + ae83 + afumott + aescmae + 
                                                  arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_smoking))
summary(pool(bygroup_est_y18_smoking), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_smoking)
miceMCError(pooled)
rm(bygroup_est_y18_smoking)

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup_75, glm(y18_alcohol ~ y11_cte_recode + asexo + 
                                                  aethnicity + ae83 + afumott + aescmae + 
                                                  arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
summary(pool(bygroup_est_y18_alcohol), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)
rm(bygroup_est_y18_alcohol)

# drug use
bygroup_est_y18_druguse <- with(bygroup_75, glm(y18_druguse_recode ~ y11_cte_recode + asexo + 
                                                  aethnicity + ae83 + afumott + aescmae + 
                                                  arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_druguse))
summary(pool(bygroup_est_y18_druguse), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_druguse)
miceMCError(pooled)
rm(bygroup_est_y18_druguse)

# overall PA
bygroup_est_joverall_pa <- with(bygroup_75, lm(joverall_pa ~ y11_cte_recode + asexo + 
                                                 aethnicity + ae83 + afumott + aescmae + 
                                                 arendtot + birthday))
summary(pool(bygroup_est_joverall_pa))
summary(pool(bygroup_est_joverall_pa), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_joverall_pa)
miceMCError(pooled)
rm(bygroup_est_joverall_pa)

# bouted PA
bygroup_est_jmvpab5 <- with(bygroup_75, lm(jmvpab5 ~ y11_cte_recode + asexo + 
                                             aethnicity + ae83 + afumott + aescmae + 
                                             arendtot + birthday))
summary(pool(bygroup_est_jmvpab5))
summary(pool(bygroup_est_jmvpab5), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_jmvpab5)
miceMCError(pooled)
rm(bygroup_est_jmvpab5)

# UPF
bygroup_est_jpercultraprocessadoscal <- with(bygroup_75, lm(jpercultraprocessadoscal ~ y11_cte_recode + asexo + 
                                                              aethnicity + ae83 + afumott + aescmae + 
                                                              arendtot + birthday))
summary(pool(bygroup_est_jpercultraprocessadoscal))
summary(pool(bygroup_est_jpercultraprocessadoscal), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_jpercultraprocessadoscal)
miceMCError(pooled)
rm(bygroup_est_jpercultraprocessadoscal)
#### all monte carlo errors are less than 10% of the standard errors


#### export to Stata (y11 CTE model, 75 imps) ####
# convert to long format 
completed_data_75 <- complete(bygroup_75, action = "long", include = TRUE)
# rename the .imp and .id columns
completed_data_75_clean <- completed_data_75 %>%
  rename(imp = .imp, id = .id)

# change factor to numeric 
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(asexo = as.numeric(as.character(asexo)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(ae83 = as.numeric(as.character(ae83)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(afumott = as.numeric(as.character(afumott)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_alcohol = as.numeric(as.character(y18_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_smoking = as.numeric(as.character(y18_smoking)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_druguse = as.numeric(as.character(y18_druguse_recode)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_alcohol = as.numeric(as.character(y11_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_smoking = as.numeric(as.character(y11_smoking)))

# Check
table(completed_data_75_clean$asexo, useNA = "ifany")
table(completed_data_75_clean$aethnicity, useNA = "ifany")
table(completed_data_75_clean$ae83, useNA = "ifany")
table(completed_data_75_clean$afumott, useNA = "ifany")
table(completed_data_75_clean$y18_alcohol, useNA = "ifany")
table(completed_data_75_clean$y18_smoking, useNA = "ifany")
table(completed_data_75_clean$y18_druguse_recode, useNA = "ifany")
table(completed_data_75_clean$y11_alcohol, useNA = "ifany")
table(completed_data_75_clean$y11_smoking, useNA = "ifany")

# Check 2
str(completed_data_75_clean$asexo)
str(completed_data_75_clean$aethnicity)
str(completed_data_75_clean$ae83)
str(completed_data_75_clean$afumott)
str(completed_data_75_clean$y18_alcohol)
str(completed_data_75_clean$y18_smoking)
str(completed_data_75_clean$y18_druguse_recode)
str(completed_data_75_clean$y11_alcohol)
str(completed_data_75_clean$y11_smoking)

# write to a .dta file
library(haven)
write_dta(completed_data_75_clean, "imputed_data_long.dta")






# 7. Transformed externalising imputation model: Passive imputation, y11_cte --------------------------------
# As a non-linear relationship between externalising problems and problematic alcohol 
# use were found, we passively impute a transformed externalising variable in separate models 

# Create MI data frame (exclude pp_ID, exclude y15_substance use, y15_ctspc, y6_cte_recode)
MI_df <- dat[,c("y11_cte_recode","asexo","aethnicity","ae83","afumott","aescmae",
                "arendtot","birthday","y15_internalising","y15_externalising","y18_alcohol",
                "y18_smoking","y18_druguse_recode","joverall_pa","jpercultraprocessadoscal","jmvpab5",
                "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                "y11_alcohol","y11_smoking","goverall_pa","gpercultraprocessados",
                "gmvpab5","aa07","y6_ctspc","y11_ctspc")]
str(MI_df)
head(MI_df)

install.packages("mice")

library(mice)
str(MI_df)
options(max.print = 1000000)
getOption('max.print')


##### create a transformed externalising variable for the ext-alcohol models
MI_df$y15_externalising_transformed <- (((MI_df$y15_externalising + 1) / 10)^-1)
psych::describe(MI_df$y15_externalising)
psych::describe(MI_df$y15_externalising_transformed)
hist(MI_df$y15_externalising)
hist(MI_df$y15_externalising_transformed)

#Perform a "dryrun" 
testimp <- mice(data=MI_df, maxit = 0,
                printFlag=FALSE)
testimp
testimp$formulas

### define prediction matrix ####
predMatrix_1<-matrix(NA, nrow=30 ,ncol=30)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df)
predMatrix_1

predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aethnicity",]        <-c(1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) 
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aescmae",]           <-c(1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1) 
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["y15_internalising",] <-c(1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1)
predMatrix_1["y15_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0) 
predMatrix_1["y18_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1) 
predMatrix_1["y18_smoking",]       <-c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1) 
predMatrix_1["y18_druguse_recode",]       <-c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1) 
predMatrix_1["joverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1) 
predMatrix_1["jpercultraprocessadoscal",]         <-c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1) 
predMatrix_1["jmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1) 
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_internalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1)
predMatrix_1["y11_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y18_internalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y18_externalising",] <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y11_alcohol",]       <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_smoking",]       <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["goverall_pa",]       <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["gpercultraprocessados",]             <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["gmvpab5",]           <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["aa07",]              <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1)   
predMatrix_1["y6_ctspc",]          <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1)
predMatrix_1["y11_ctspc",]         <-c(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1)
predMatrix_1["y15_externalising_transformed",] <-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1    

# define passive imputation method
meth <- testimp$method
meth["y15_externalising_transformed"] <- "~I(((y15_externalising + 1)/10)^-1)"
meth

# if we want to run mice with this new matrix, we pull this matrix into the command
testimp2 <- mice(data=MI_df, predictorMatrix = predMatrix_1, method = meth,
                 maxit = 0, printFlag=FALSE)

# check it's worked
testimp2$predictorMatrix 
testimp2$formulas 
testimp2$method 

### Imputation using stratified imputation ####
#Create separate male and female datasets
MI_df_female <- subset(MI_df, asexo==0)
MI_df_male <- subset(MI_df, asexo==1)

str(MI_df_female)
str(MI_df_male)
table(MI_df_female$asexo)
table(MI_df_male$asexo)

# NOTE - unclear whether I should remove asexo as a predictor in my matrices
# as it will now have the same value for all individuals in each subgroup... test

# dry run for each group (check prediction matrix and method).   
testimp_female <- mice(data=MI_df_female, predictorMatrix = predMatrix_1, method = meth,
                       maxit=0, printFlag=FALSE, seed=456)
testimp_male <- mice(data=MI_df_male, predictorMatrix = predMatrix_1, method = meth,
                     maxit=0, printFlag=FALSE, seed=456)
testimp_female 
testimp_male

#Perform the imputation - maxit=25, m=75
bygroup_female<-mice(MI_df_female, maxit=25, m=75, seed=123, printFlag=FALSE, method = meth,
                     predictorMatrix = predMatrix_1)
bygroup_male<-mice(MI_df_male, maxit=25, m=75, seed=123, printFlag=FALSE, method = meth,
                   predictorMatrix = predMatrix_1)
# note - get warning 'Number of logged events: 1'

# Combine the two sets of imputations using `rbind`
bygroup <- rbind(bygroup_female, bygroup_male)
bygroup

setwd("../mice_090725")
save(bygroup, file = "imp_75_transformed_ext.rda")

#To access the first imputed dataset 
bygroup1 <- complete(bygroup,1)
head(bygroup1)


### The maximum number of iterations / traceplots ####
plot(bygroup, c("y11_cte_recode", "aethnicity", "aescmae"), layout = c(2, 4))
plot(bygroup, c("y15_internalising", "y15_externalising", "y15_externalising_transformed"), layout = c(2, 3))
plot(bygroup, c("y18_alcohol", "y18_smoking", "y18_druguse"), layout = c(2, 3))
plot(bygroup, c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5"), layout = c(2, 3))

### Strip plots ####
# Compare observed (blue) and imputed (red) data (continuous variables; m = 50)
stripplot(bygroup, y11_cte_recode + aethnicity + aescmae ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup, y15_internalising + y15_externalising + y15_externalising_transformed ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup, y18_alcohol + y18_smoking + y18_druguse + joverall_pa + jpercultraprocessadoscal + jmvpab5 ~ .imp,
          pch = 20, cex = 1.2)

### Density plots ####
densityplot(bygroup, data = ~ y11_cte_recode) 
densityplot(bygroup, data = ~ aethnicity) 
densityplot(bygroup, data = ~ aescmae) 
densityplot(bygroup, data = ~ y15_internalising + y15_externalising + y15_externalising_transformed) 
densityplot(bygroup, data = ~ y18_alcohol + y18_smoking + y18_druguse_recode) 
densityplot(bygroup, data = ~ joverall_pa + jpercultraprocessadoscal + jmvpab5) 


### Distribution of observed and imputed values - histograms ####
par(mfrow=c(1,3))
hist(bygroup$data$y15_externalising_transformed, xlab="y15_externalising_transformed",
     main="Observed values")
abline(v=mean(bygroup$data$y15_externalising_transformed, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup$imp$y15_externalising_transformed[[1]], xlab="y15_externalising_transformed",
     main="Imputed values (first imp)")
abline(v=mean(bygroup$imp$y15_externalising_transformed[[1]]),lty=2, col="red")

hist(unlist(bygroup$imp$y15_externalising_transformed), xlab="y15_externalising_transformed",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup$imp$y15_externalising_transformed)),lty=2, col="red")


### Analyse imputed dataset with practice lm models ####
bygroup_est_joverall_pa <- with(bygroup, lm(joverall_pa ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jpercultraprocessadoscal <- with(bygroup, lm(jpercultraprocessadoscal ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jmvpab5 <- with(bygroup, lm(jmvpab5 ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y18_alcohol <- with(bygroup, glm(y18_alcohol ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_smoking <- with(bygroup, glm(y18_smoking ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_druguse <- with(bygroup, glm(y18_druguse_recode ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y15_externalising <- with(bygroup, lm(y15_externalising ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y15_internalising <- with(bygroup, lm(y15_internalising ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y15_externalising_trans <- with(bygroup, lm(y15_externalising_transformed ~ y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)
bygroup_est_y15_externalising_pool <- pool(bygroup_est_y15_externalising)
bygroup_est_y15_internalising_pool <- pool(bygroup_est_y15_internalising)
bygroup_est_y15_externalising_trans_pool <- pool(bygroup_est_y15_externalising_trans)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE) 
summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE) 
summary(bygroup_est_jmvpab5_pool, conf.int = TRUE) 
summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE) 
summary(bygroup_est_y18_smoking_pool, conf.int = TRUE) 
summary(bygroup_est_y18_druguse_pool, conf.int = TRUE) 
summary(bygroup_est_y15_externalising_pool, conf.int = TRUE) 
summary(bygroup_est_y15_internalising_pool, conf.int = TRUE) 
summary(bygroup_est_y15_externalising_trans_pool, conf.int = TRUE)

# now for Ext transformed
bygroup_est_y18_alcohol_trans <- with(bygroup, glm(y18_alcohol ~ y15_internalising + y15_externalising_transformed + y11_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
#Pool the results
bygroup_est_y18_alcohol_trans_pool <- pool(bygroup_est_y18_alcohol_trans)
# summarise results
summary(bygroup_est_y18_alcohol_trans_pool, conf.int = TRUE) 

## Monte Carlo error ####
# function from https://thestatsgeek.com/2022/02/23/how-many-imputations-with-mice-assessing-monte-carlo-error-after-multiple-imputation-in-r/

# function:
miceMCError <- function(pooledRes) {
  monteCarloSE <- sqrt(pooledRes$pooled$b/pooledRes$m)
  ciLower <- pooledRes$pooled$estimate - qt(0.975,df=pooledRes$m-1)*monteCarloSE
  ciUpper <- pooledRes$pooled$estimate + qt(0.975,df=pooledRes$m-1)*monteCarloSE
  mcTable <- cbind(pooledRes$pooled$estimate, monteCarloSE, ciLower, ciUpper)
  colnames(mcTable) <- c("Estimate", "Monte Carlo SE", "95% CI lower limit", "95% CI upper limit")
  print(mcTable)
  print("Warning: 95% CI only quantifies Monte-Carlo uncertainty!")
}

## unadjusted models to inspect MC error
# smoking
bygroup_est_y18_smoking <- with(bygroup, glm(y18_smoking ~ y11_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_smoking))
summary(pool(bygroup_est_y18_smoking), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_smoking)
miceMCError(pooled)

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup, glm(y18_alcohol ~ y11_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
summary(pool(bygroup_est_y18_alcohol), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)

# drug use
bygroup_est_y18_druguse <- with(bygroup, glm(y18_druguse_recode ~ y11_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_druguse))
summary(pool(bygroup_est_y18_druguse), conf.int = TRUE, exponentiate = TRUE) # this provides ORs
# pool for mc error
pooled <- pool(bygroup_est_y18_druguse)
miceMCError(pooled)

# overall PA
bygroup_est_joverall_pa <- with(bygroup, lm(joverall_pa ~ y11_cte_recode))
summary(pool(bygroup_est_joverall_pa))
summary(pool(bygroup_est_joverall_pa), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_joverall_pa)
miceMCError(pooled)

# bouted PA
bygroup_est_jmvpab5 <- with(bygroup, lm(jmvpab5 ~ y11_cte_recode))
summary(pool(bygroup_est_jmvpab5))
summary(pool(bygroup_est_jmvpab5), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_jmvpab5)
miceMCError(pooled)

# UPF
bygroup_est_jpercultraprocessadoscal <- with(bygroup, lm(jpercultraprocessadoscal ~ y11_cte_recode))
summary(pool(bygroup_est_jpercultraprocessadoscal))
summary(pool(bygroup_est_jpercultraprocessadoscal), conf.int = TRUE) 
# pool for mc error
pooled <- pool(bygroup_est_jpercultraprocessadoscal)
miceMCError(pooled)
#### all monte carlo errors are less than 10% of the standard errors

## adjusted
# smoking
bygroup_est_y18_smoking <- with(bygroup, glm(y18_smoking ~ y11_cte_recode + asexo + 
                                                  aethnicity + ae83 + afumott + aescmae + 
                                                  arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_smoking))
# pool for mc error
pooled <- pool(bygroup_est_y18_smoking)
miceMCError(pooled)
rm(bygroup_est_y18_smoking)

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup, glm(y18_alcohol ~ y11_cte_recode + asexo + 
                                                  aethnicity + ae83 + afumott + aescmae + 
                                                  arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)
rm(bygroup_est_y18_alcohol)

# drug use
bygroup_est_y18_druguse <- with(bygroup, glm(y18_druguse_recode ~ y11_cte_recode + asexo + 
                                                  aethnicity + ae83 + afumott + aescmae + 
                                                  arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_druguse))
# pool for mc error
pooled <- pool(bygroup_est_y18_druguse)
miceMCError(pooled)
rm(bygroup_est_y18_druguse)

# overall PA
bygroup_est_joverall_pa <- with(bygroup, lm(joverall_pa ~ y11_cte_recode + asexo + 
                                                 aethnicity + ae83 + afumott + aescmae + 
                                                 arendtot + birthday))
summary(pool(bygroup_est_joverall_pa))
# pool for mc error
pooled <- pool(bygroup_est_joverall_pa)
miceMCError(pooled)
rm(bygroup_est_joverall_pa)

# bouted PA
bygroup_est_jmvpab5 <- with(bygroup, lm(jmvpab5 ~ y11_cte_recode + asexo + 
                                             aethnicity + ae83 + afumott + aescmae + 
                                             arendtot + birthday))
summary(pool(bygroup_est_jmvpab5))
# pool for mc error
pooled <- pool(bygroup_est_jmvpab5)
miceMCError(pooled)
rm(bygroup_est_jmvpab5)

# UPF
bygroup_est_jpercultraprocessadoscal <- with(bygroup, lm(jpercultraprocessadoscal ~ y11_cte_recode + asexo + 
                                                              aethnicity + ae83 + afumott + aescmae + 
                                                              arendtot + birthday))
summary(pool(bygroup_est_jpercultraprocessadoscal))
# pool for mc error
pooled <- pool(bygroup_est_jpercultraprocessadoscal)
miceMCError(pooled)
rm(bygroup_est_jpercultraprocessadoscal)
#### all monte carlo errors are less than 10% of the standard errors


#### export to Stata (y11, transformed externalising, 75 imps) #####
# first, convert to long format 
completed_data_75 <- complete(bygroup, action = "long", include = TRUE)

# rename the .imp and .id columns
completed_data_75_clean <- completed_data_75 %>%
  rename(imp = .imp, id = .id)

str(completed_data_75_clean)

# change to numeric 
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(asexo = as.numeric(as.character(asexo)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(ae83 = as.numeric(as.character(ae83)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(afumott = as.numeric(as.character(afumott)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_alcohol = as.numeric(as.character(y18_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_smoking = as.numeric(as.character(y18_smoking)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_druguse = as.numeric(as.character(y18_druguse_recode)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_alcohol = as.numeric(as.character(y11_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_smoking = as.numeric(as.character(y11_smoking)))

# Check
table(completed_data_75_clean$asexo, useNA = "ifany")
table(completed_data_75_clean$aethnicity, useNA = "ifany")
table(completed_data_75_clean$ae83, useNA = "ifany")
table(completed_data_75_clean$afumott, useNA = "ifany")
table(completed_data_75_clean$y18_alcohol, useNA = "ifany")
table(completed_data_75_clean$y18_smoking, useNA = "ifany")
table(completed_data_75_clean$y18_druguse_recode, useNA = "ifany")
table(completed_data_75_clean$y11_alcohol, useNA = "ifany")
table(completed_data_75_clean$y11_smoking, useNA = "ifany")

# Check 2
str(completed_data_75_clean$asexo)
str(completed_data_75_clean$aethnicity)
str(completed_data_75_clean$ae83)
str(completed_data_75_clean$afumott)
str(completed_data_75_clean$y18_alcohol)
str(completed_data_75_clean$y18_smoking)
str(completed_data_75_clean$y18_druguse_recode)
str(completed_data_75_clean$y11_alcohol)
str(completed_data_75_clean$y11_smoking)

# write to a .dta file
library(haven)
write_dta(completed_data_75_clean, "imputed_data_long_transExt.dta")





# 8. Age 15 cumulative trauma imputation model with 75 imputed datasets   --------------------------------
# Create MI data frame
MI_df <- dat[,c("y11_cte_recode","asexo","aethnicity","ae83","afumott","aescmae",
                "arendtot","birthday","y15_internalising","y15_externalising","y18_alcohol",
                "y18_smoking","y18_druguse_recode","joverall_pa","jpercultraprocessadoscal","jmvpab5",
                "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                "y11_alcohol","y11_smoking","goverall_pa","gpercultraprocessados",
                "gmvpab5","aa07","y6_ctspc","y11_ctspc")]
str(MI_df)
head(MI_df)

library(mice)
str(MI_df)
options(max.print = 1000000)
getOption('max.print')

#Perform a "dryrun" 
testimp <- mice(data=MI_df, maxit = 0,
                printFlag=FALSE)
testimp
testimp$formulas

### define prediction matrix ####
predMatrix_1<-matrix(NA, nrow=29 ,ncol=29)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df)
predMatrix_1

# exclude y11_cte from all predictions
# y15_cte_recode is included as a predictor for all analysis variables, instead of y11
predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aethnicity",]        <-c(0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aescmae",]           <-c(0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["y15_internalising",] <-c(0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0) 
predMatrix_1["y15_externalising",] <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0) 
predMatrix_1["y18_alcohol",]       <-c(0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0) 
predMatrix_1["y18_smoking",]       <-c(0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0)
predMatrix_1["y18_druguse_recode",]       <-c(0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0) 
predMatrix_1["joverall_pa",]       <-c(0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0)  
predMatrix_1["jpercultraprocessadoscal",]         <-c(0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0) 
predMatrix_1["jmvpab5",]           <-c(0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0)
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["y11_internalising",] <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_externalising",] <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1) 
predMatrix_1["y18_internalising",] <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y18_externalising",] <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_alcohol",]       <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["y11_smoking",]       <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1) 
predMatrix_1["goverall_pa",]       <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["gpercultraprocessados",]             <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["gmvpab5",]           <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["aa07",]              <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)   
predMatrix_1["y6_ctspc",]          <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1)
predMatrix_1["y11_ctspc",]         <-c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0)
predMatrix_1    

testimp2 <- mice(data=MI_df, predictorMatrix = predMatrix_1,
                 maxit = 0, printFlag=FALSE)
testimp2$predictorMatrix 
testimp2$formulas 
testimp2 

### Imputation using stratified imputation: y15_cte_recode ####
MI_df_female <- subset(MI_df, asexo==0)
MI_df_male <- subset(MI_df, asexo==1)

str(MI_df_female)
str(MI_df_male)
table(MI_df_female$asexo)
table(MI_df_male$asexo)

# dry run for each group (check prediction matrix and method)
testimp_female <- mice(data=MI_df_female, predictorMatrix = predMatrix_1,
                       maxit=0, printFlag=FALSE, seed=456)
testimp_male <- mice(data=MI_df_male, predictorMatrix = predMatrix_1,
                     maxit=0, printFlag=FALSE, seed=456)
testimp_female 
testimp_male
testimp_female$method
testimp_male$method

#Perform the imputation - maxit=25, m=75
bygroup_female_15cte_75<-mice(MI_df_female, maxit=25, m=75, seed=123, printFlag=FALSE, 
                           predictorMatrix = predMatrix_1)
bygroup_male_15cte_75<-mice(MI_df_male, maxit=25, m=75, seed=123, printFlag=FALSE, 
                         predictorMatrix = predMatrix_1)
bygroup_female_15cte_75
bygroup_male_15cte_75

# Combine the two sets of imputations 
bygroup_15cte_75 <- rbind(bygroup_female_15cte_75, bygroup_male_15cte_75)
bygroup_15cte_75

# to check
bygroup_15cte_1 <- complete(bygroup_15cte_75,1)
head(bygroup_15cte_1) 

setwd("../mice_260126")
save(bygroup_15cte_75, file = "imp_75_15cte.rda")
load(file="imp_75_15cte.rda")

### Analyse imputed dataset with practice lm models ####
bygroup_est_joverall_pa <- with(bygroup_15cte_75, lm(joverall_pa ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jpercultraprocessadoscal <- with(bygroup_15cte_75, lm(jpercultraprocessadoscal ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jmvpab5 <- with(bygroup_15cte_75, lm(jmvpab5 ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y18_alcohol <- with(bygroup_15cte_75, glm(y18_alcohol ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_smoking <- with(bygroup_15cte_75, glm(y18_smoking ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_druguse <- with(bygroup_15cte_75, glm(y18_druguse_recode ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y15_externalising <- with(bygroup_15cte_75, lm(y15_externalising ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y15_internalising <- with(bygroup_15cte_75, lm(y15_internalising ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)
bygroup_est_y15_externalising_pool <- pool(bygroup_est_y15_externalising)
bygroup_est_y15_internalising_pool <- pool(bygroup_est_y15_internalising)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE)
rm(bygroup_est_joverall_pa)
rm(bygroup_est_joverall_pa_pool)

summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE)
rm(bygroup_est_jpercultraprocessadoscal)
rm(bygroup_est_jpercultraprocessadoscal_pool)

summary(bygroup_est_jmvpab5_pool, conf.int = TRUE)
rm(bygroup_est_jmvpab5)
rm(bygroup_est_jmvpab5_pool)

summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE)
rm(bygroup_est_y18_alcohol)
rm(bygroup_est_y18_alcohol_pool)

summary(bygroup_est_y18_smoking_pool, conf.int = TRUE)
rm(bygroup_est_y18_smoking)
rm(bygroup_est_y18_smoking_pool)

summary(bygroup_est_y18_druguse_pool, conf.int = TRUE)
rm(bygroup_est_y18_druguse)
rm(bygroup_est_y18_druguse_pool)

summary(bygroup_est_y15_externalising_pool, conf.int = TRUE)
rm(bygroup_est_y15_externalising)
rm(bygroup_est_y15_externalising_pool)

summary(bygroup_est_y15_internalising_pool, conf.int = TRUE)
rm(bygroup_est_y15_internalising)
rm(bygroup_est_y15_internalising_pool)

# now for Int/Ext 
bygroup_est_joverall_pa <- with(bygroup_15cte_75, lm(joverall_pa ~ y15_internalising + y15_externalising + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jpercultraprocessadoscal <- with(bygroup_15cte_75, lm(jpercultraprocessadoscal ~ y15_internalising + y15_externalising + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_jmvpab5 <- with(bygroup_15cte_75, lm(jmvpab5 ~ y15_internalising + y15_externalising + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
bygroup_est_y18_alcohol <- with(bygroup_15cte_75, glm(y18_alcohol ~ y15_internalising + y15_externalising + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_smoking <- with(bygroup_15cte_75, glm(y18_smoking ~ y15_internalising + y15_externalising + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
bygroup_est_y18_druguse <- with(bygroup_15cte_75, glm(y18_druguse_recode ~ y15_internalising + y15_externalising + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))

#Pool the results
bygroup_est_joverall_pa_pool <- pool(bygroup_est_joverall_pa)
bygroup_est_jpercultraprocessadoscal_pool <- pool(bygroup_est_jpercultraprocessadoscal)
bygroup_est_jmvpab5_pool <- pool(bygroup_est_jmvpab5)
bygroup_est_y18_alcohol_pool <- pool(bygroup_est_y18_alcohol)
bygroup_est_y18_smoking_pool <- pool(bygroup_est_y18_smoking)
bygroup_est_y18_druguse_pool <- pool(bygroup_est_y18_druguse)

# summarise results
summary(bygroup_est_joverall_pa_pool, conf.int = TRUE)
rm(bygroup_est_joverall_pa)
rm(bygroup_est_joverall_pa_pool)

summary(bygroup_est_jpercultraprocessadoscal_pool, conf.int = TRUE)
rm(bygroup_est_jpercultraprocessadoscal)
rm(bygroup_est_jpercultraprocessadoscal_pool)

summary(bygroup_est_jmvpab5_pool, conf.int = TRUE)
rm(bygroup_est_jmvpab5)
rm(bygroup_est_jmvpab5_pool)

summary(bygroup_est_y18_alcohol_pool, conf.int = TRUE)
rm(bygroup_est_y18_alcohol)
rm(bygroup_est_y18_alcohol_pool)

summary(bygroup_est_y18_smoking_pool, conf.int = TRUE)
rm(bygroup_est_y18_smoking)
rm(bygroup_est_y18_smoking_pool)

summary(bygroup_est_y18_druguse_pool, conf.int = TRUE)
rm(bygroup_est_y18_druguse)
rm(bygroup_est_y18_druguse_pool)

### The maximum number of iterations / traceplots (int/norm) ####
plot(bygroup_15cte_75, c("y15_cte_recode","y11_cte_recode", "aethnicity", "aescmae"), layout = c(2, 4))
plot(bygroup_15cte_75, c("y15_internalising", "y15_externalising"), layout = c(2, 2))
plot(bygroup_15cte_75, c("y18_alcohol", "y18_smoking", "y18_druguse_recode"), layout = c(2, 3))
plot(bygroup_15cte_75, c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5"), layout = c(2, 3))


### Distribution of observed and imputed values - histograms ####
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$y15_cte_recode, xlab="Cumulative trauma at 15 years",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$y15_cte_recode, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$y15_cte_recode[[1]], xlab="Cumulative trauma at 15 years",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$y15_cte_recode[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$y15_cte_recode), xlab="Cumulative trauma at 15 years",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$y15_cte_recode)),lty=2, col="red")

# internalising
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$y15_internalising, xlab="Internalising at 15 years",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$y15_internalising, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$y15_internalising[[1]], xlab="Internalising at 15 years",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$y15_internalising[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$y15_internalising), xlab="Internalising at 15 years",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$y15_internalising)),lty=2, col="red")


#externalising
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$y15_externalising, xlab="Externalising at 15 years",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$y15_externalising, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$y15_externalising[[1]], xlab="Externalising at 15 years",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$y15_externalising[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$y15_externalising), xlab="Externalising at 15 years",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$y15_externalising)),lty=2, col="red")


#aescmae
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$aescmae, xlab="Maternal education",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$aescmae, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$aescmae[[1]], xlab="Maternal education",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$aescmae[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$aescmae), xlab="Maternal education",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$aescmae)),lty=2, col="red")


# joverall_pa
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$joverall_pa, xlab="Overall PA",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$joverall_pa, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$joverall_pa[[1]], xlab="Overall PA",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$joverall_pa[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$joverall_pa), xlab="Overall PA",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$joverall_pa)),lty=2, col="red")


#jpercultraprocessadoscal
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$jpercultraprocessadoscal, xlab="UPF intake",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$jpercultraprocessadoscal, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$jpercultraprocessadoscal[[1]], xlab="UPF intake",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$jpercultraprocessadoscal[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$jpercultraprocessadoscal), xlab="UPF intake",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$jpercultraprocessadoscal)),lty=2, col="red")


#jmvpab5
par(mfrow=c(1,3))
hist(bygroup_15cte_75$data$jmvpab5, xlab="bouted PA",
     main="Observed values")
abline(v=mean(bygroup_15cte_75$data$jmvpab5, na.rm=TRUE),
       lty=2, col="red")

hist(bygroup_15cte_75$imp$jmvpab5[[1]], xlab="bouted PA",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_75$imp$jmvpab5[[1]]),lty=2, col="red")

hist(unlist(bygroup_15cte_75$imp$jmvpab5), xlab="bouted PA",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_75$imp$jmvpab5)),lty=2, col="red")


### Strip plots ####
# Compare observed (blue) and imputed (red) data (continuous variables; m = 50)
stripplot(bygroup_15cte_75, y15_cte_recode + y11_cte_recode + aethnicity + aescmae ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup_15cte_75, y15_internalising + y15_externalising ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup_15cte_75, y18_alcohol + y18_smoking + y18_druguse_recode + joverall_pa + jpercultraprocessadoscal + jmvpab5 ~ .imp,
          pch = 20, cex = 1.2)


### Density plots ####
densityplot(bygroup_15cte_75, data = ~ y15_cte_recode) 
densityplot(bygroup_15cte_75, data = ~ y11_cte_recode) 
densityplot(bygroup_15cte_75, data = ~ aethnicity) 
densityplot(bygroup_15cte_75, data = ~ aescmae) 
densityplot(bygroup_15cte_75, data = ~ y15_internalising + y15_externalising) 
densityplot(bygroup_15cte_75, data = ~ y18_alcohol + y18_smoking + y18_druguse_recode) 
densityplot(bygroup_15cte_75, data = ~ joverall_pa + jpercultraprocessadoscal + jmvpab5) 


## Monte Carlo error ####
# function from https://thestatsgeek.com/2022/02/23/how-many-imputations-with-mice-assessing-monte-carlo-error-after-multiple-imputation-in-r/
miceMCError <- function(pooledRes) {
  monteCarloSE <- sqrt(pooledRes$pooled$b/pooledRes$m)
  ciLower <- pooledRes$pooled$estimate - qt(0.975,df=pooledRes$m-1)*monteCarloSE
  ciUpper <- pooledRes$pooled$estimate + qt(0.975,df=pooledRes$m-1)*monteCarloSE
  mcTable <- cbind(pooledRes$pooled$estimate, monteCarloSE, ciLower, ciUpper)
  colnames(mcTable) <- c("Estimate", "Monte Carlo SE", "95% CI lower limit", "95% CI upper limit")
  print(mcTable)
  print("Warning: 95% CI only quantifies Monte-Carlo uncertainty!")
}


## unadjusted models 
# smoking
bygroup_est_y18_smoking <- with(bygroup_15cte_75, glm(y18_smoking ~ y15_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_smoking))
# pool for mc error
pooled <- pool(bygroup_est_y18_smoking)
miceMCError(pooled)

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup_15cte_75, glm(y18_alcohol ~ y15_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)

# drug use
bygroup_est_y18_druguse <- with(bygroup_15cte_75, glm(y18_druguse_recode ~ y15_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_druguse))
# pool for mc error
pooled <- pool(bygroup_est_y18_druguse)
miceMCError(pooled)

# overall PA
bygroup_est_joverall_pa <- with(bygroup_15cte_75, lm(joverall_pa ~ y15_cte_recode))
summary(pool(bygroup_est_joverall_pa))
# pool for mc error
pooled <- pool(bygroup_est_joverall_pa)
miceMCError(pooled)

# bouted PA
bygroup_est_jmvpab5 <- with(bygroup_15cte_75, lm(jmvpab5 ~ y15_cte_recode))
summary(pool(bygroup_est_jmvpab5))
# pool for mc error
pooled <- pool(bygroup_est_jmvpab5)
miceMCError(pooled)

# UPF
bygroup_est_jpercultraprocessadoscal <- with(bygroup_15cte_75, lm(jpercultraprocessadoscal ~ y15_cte_recode))
summary(pool(bygroup_est_jpercultraprocessadoscal))
# pool for mc error
pooled <- pool(bygroup_est_jpercultraprocessadoscal)
miceMCError(pooled)
#### all monte carlo errors are less than 10% of the standard errors

## adjusted
# smoking
bygroup_est_y18_smoking <- with(bygroup_15cte_75, glm(y18_smoking ~ y15_cte_recode + asexo + 
                                               aethnicity + ae83 + afumott + aescmae + 
                                               arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_smoking))
# pool for mc error
pooled <- pool(bygroup_est_y18_smoking)
miceMCError(pooled)
rm(bygroup_est_y18_smoking)

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup_15cte_75, glm(y18_alcohol ~ y15_cte_recode + asexo + 
                                               aethnicity + ae83 + afumott + aescmae + 
                                               arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)
rm(bygroup_est_y18_alcohol)

# drug use
bygroup_est_y18_druguse <- with(bygroup_15cte_75, glm(y18_druguse_recode ~ y15_cte_recode + asexo + 
                                               aethnicity + ae83 + afumott + aescmae + 
                                               arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_druguse))
# pool for mc error
pooled <- pool(bygroup_est_y18_druguse)
miceMCError(pooled)
rm(bygroup_est_y18_druguse)

# overall PA
bygroup_est_joverall_pa <- with(bygroup_15cte_75, lm(joverall_pa ~ y15_cte_recode + asexo + 
                                              aethnicity + ae83 + afumott + aescmae + 
                                              arendtot + birthday))
summary(pool(bygroup_est_joverall_pa))
# pool for mc error
pooled <- pool(bygroup_est_joverall_pa)
miceMCError(pooled)
rm(bygroup_est_joverall_pa)

# bouted PA
bygroup_est_jmvpab5 <- with(bygroup_15cte_75, lm(jmvpab5 ~ y15_cte_recode + asexo + 
                                          aethnicity + ae83 + afumott + aescmae + 
                                          arendtot + birthday))
summary(pool(bygroup_est_jmvpab5))
# pool for mc error
pooled <- pool(bygroup_est_jmvpab5)
miceMCError(pooled)
rm(bygroup_est_jmvpab5)

# UPF
bygroup_est_jpercultraprocessadoscal <- with(bygroup_15cte_75, lm(jpercultraprocessadoscal ~ y15_cte_recode + asexo + 
                                                           aethnicity + ae83 + afumott + aescmae + 
                                                           arendtot + birthday))
summary(pool(bygroup_est_jpercultraprocessadoscal))
# pool for mc error
pooled <- pool(bygroup_est_jpercultraprocessadoscal)
miceMCError(pooled)
rm(bygroup_est_jpercultraprocessadoscal)
#### all monte carlo errors are less than 10% of the standard errors



#### export to Stata (y15, 75 imps) #####
# first, convert to long format 
completed_data_75 <- complete(bygroup_15cte_75, action = "long", include = TRUE)

# rename the .imp and .id columns
completed_data_75_clean <- completed_data_75 %>%
  rename(imp = .imp, id = .id)
str(completed_data_75_clean)

# change to numeric 
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(asexo = as.numeric(as.character(asexo)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(ae83 = as.numeric(as.character(ae83)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(afumott = as.numeric(as.character(afumott)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_alcohol = as.numeric(as.character(y18_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_smoking = as.numeric(as.character(y18_smoking)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_druguse = as.numeric(as.character(y18_druguse_recode)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_alcohol = as.numeric(as.character(y11_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_smoking = as.numeric(as.character(y11_smoking)))

# Check
table(completed_data_75_clean$asexo, useNA = "ifany")
table(completed_data_75_clean$aethnicity, useNA = "ifany")
table(completed_data_75_clean$ae83, useNA = "ifany")
table(completed_data_75_clean$afumott, useNA = "ifany")
table(completed_data_75_clean$y18_alcohol, useNA = "ifany")
table(completed_data_75_clean$y18_smoking, useNA = "ifany")
table(completed_data_75_clean$y18_druguse_recode, useNA = "ifany")
table(completed_data_75_clean$y11_alcohol, useNA = "ifany")
table(completed_data_75_clean$y11_smoking, useNA = "ifany")

# Check 2
str(completed_data_75_clean$asexo)
str(completed_data_75_clean$aethnicity)
str(completed_data_75_clean$ae83)
str(completed_data_75_clean$afumott)
str(completed_data_75_clean$y18_alcohol)
str(completed_data_75_clean$y18_smoking)
str(completed_data_75_clean$y18_druguse_recode)
str(completed_data_75_clean$y11_alcohol)
str(completed_data_75_clean$y11_smoking)

# write to a .dta file
library(haven)
write_dta(completed_data_75_clean, "imputed_data_long_y15.dta")



# 9. Transformed externalising imputation model: Passive imputation, y15_cte --------------------------------
# Create MI data frame (exclude pp_ID, exclude y15_substance use, y15_ctspc, y6_cte_recode)
MI_df <- dat[,c("y11_cte_recode","asexo","aethnicity","ae83","afumott","aescmae",
                "arendtot","birthday","y15_internalising","y15_externalising","y18_alcohol",
                "y18_smoking","y18_druguse_recode","joverall_pa","jpercultraprocessadoscal","jmvpab5",
                "y15_cte_recode","y11_internalising","y11_externalising","y18_internalising","y18_externalising",
                "y11_alcohol","y11_smoking","goverall_pa","gpercultraprocessados",
                "gmvpab5","aa07","y6_ctspc","y11_ctspc")]
str(MI_df)
head(MI_df)

library(mice)
str(MI_df)
options(max.print = 1000000)
getOption('max.print')

#### Passive imputation method
##### create a transformed externalising variable for the ext-alcohol models
MI_df$y15_externalising_transformed <- (((MI_df$y15_externalising + 1) / 10)^-1)
psych::describe(MI_df$y15_externalising)
psych::describe(MI_df$y15_externalising_transformed)
hist(MI_df$y15_externalising)
hist(MI_df$y15_externalising_transformed)

#Perform a "dryrun" 
testimp <- mice(data=MI_df, maxit = 0,
                printFlag=FALSE)
testimp
testimp$formulas

### define prediction matrix ####
predMatrix_1<-matrix(NA, nrow=30 ,ncol=30)
colnames(predMatrix_1)<-rownames(predMatrix_1)<-names(MI_df)
predMatrix_1

# exclude y11_cte from all predictions
# y15_cte_recode is included as a predictor for all analysis variables, instead of y11
predMatrix_1["y11_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["asexo",]             <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aethnicity",]        <-c(0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1) 
predMatrix_1["ae83",]              <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["afumott",]           <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["aescmae",]           <-c(0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1) 
predMatrix_1["arendtot",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors  
predMatrix_1["birthday",]          <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # fully observed, no predictors
predMatrix_1["y15_internalising",] <-c(0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1) 
predMatrix_1["y15_externalising",] <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0) 
predMatrix_1["y18_alcohol",]       <-c(0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1) 
predMatrix_1["y18_smoking",]       <-c(0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,1)
predMatrix_1["y18_druguse_recode",]       <-c(0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,1) 
predMatrix_1["joverall_pa",]       <-c(0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1)  
predMatrix_1["jpercultraprocessadoscal",]         <-c(0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1) 
predMatrix_1["jmvpab5",]           <-c(0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1)
predMatrix_1["y15_cte_recode",]    <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_internalising",] <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y11_externalising",] <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y18_internalising",] <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y18_externalising",] <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["y11_alcohol",]       <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1) 
predMatrix_1["y11_smoking",]       <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1) 
predMatrix_1["goverall_pa",]       <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["gpercultraprocessados",]             <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["gmvpab5",]           <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1) 
predMatrix_1["aa07",]              <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1)   
predMatrix_1["y6_ctspc",]          <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1)
predMatrix_1["y11_ctspc",]         <-c(0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1)
predMatrix_1["y15_externalising_transformed",] <-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
predMatrix_1    

# define passive imputation method
meth <- testimp$method
meth["y15_externalising_transformed"] <- "~I(((y15_externalising + 1)/10)^-1)"
meth

# if we want to run mice with this new matrix, we pull this matrix into the command
testimp2 <- mice(data=MI_df, predictorMatrix = predMatrix_1, method = meth,
                 maxit = 0, printFlag=FALSE)

# check it's worked
testimp2$predictorMatrix 
testimp2$formulas 
testimp2$method 

### Imputation using stratified imputation: y15_cte_recode, passive imp for ext ####
#Create separate male and female datasets
MI_df_female <- subset(MI_df, asexo==0)
MI_df_male <- subset(MI_df, asexo==1)

str(MI_df_female)
str(MI_df_male)
table(MI_df_female$asexo)
table(MI_df_male$asexo)

# dry run for each group (check prediction matrix and method).   
testimp_female <- mice(data=MI_df_female, predictorMatrix = predMatrix_1, method = meth,
                       maxit=0, printFlag=FALSE, seed=456)
testimp_male <- mice(data=MI_df_male, predictorMatrix = predMatrix_1, method = meth,
                     maxit=0, printFlag=FALSE, seed=456)
testimp_female 
testimp_male

#Perform the imputation - maxit=25, m=75
bygroup_female<-mice(MI_df_female, maxit=25, m=75, seed=123, printFlag=FALSE, method = meth,
                     predictorMatrix = predMatrix_1)
bygroup_male<-mice(MI_df_male, maxit=25, m=75, seed=123, printFlag=FALSE, method = meth,
                   predictorMatrix = predMatrix_1)

# Combine the two sets of imputations
bygroup_15cte_ext <- rbind(bygroup_female, bygroup_male)
bygroup_15cte_ext

setwd("mice_260126")
save(bygroup_15cte_ext, file = "imp_75_15cte_transformed_ext.rda")
load(file="imp_75_15cte_transformed_ext.rda")


### The maximum number of iterations / traceplots ####
plot(bygroup_15cte_ext, c("y15_cte_recode", "y11_cte_recode", "aethnicity", "aescmae"), layout = c(2, 4))
plot(bygroup_15cte_ext, c("y15_internalising", "y15_externalising", "y15_externalising_transformed"), layout = c(2, 3))
plot(bygroup_15cte_ext, c("y18_alcohol", "y18_smoking", "y18_druguse_recode"), layout = c(2, 3))
plot(bygroup_15cte_ext, c("joverall_pa", "jpercultraprocessadoscal", "jmvpab5"), layout = c(2, 3))


### Strip plots ####
stripplot(bygroup_15cte_ext, y15_cte_recode + y11_cte_recode + aethnicity + aescmae ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup_15cte_ext, y15_internalising + y15_externalising + y15_externalising_transformed ~ .imp,
          pch = 20, cex = 1.2)
stripplot(bygroup_15cte_ext, y18_alcohol + y18_smoking + y18_druguse_recode + joverall_pa + jpercultraprocessadoscal + jmvpab5 ~ .imp,
          pch = 20, cex = 1.2)


### Density plots ####
densityplot(bygroup_15cte_ext, data = ~ y15_cte_recode) 
densityplot(bygroup_15cte_ext, data = ~ y11_cte_recode) 
densityplot(bygroup_15cte_ext, data = ~ aethnicity) 
densityplot(bygroup_15cte_ext, data = ~ aescmae)
densityplot(bygroup_15cte_ext, data = ~ y15_internalising + y15_externalising + y15_externalising_transformed) 
densityplot(bygroup_15cte_ext, data = ~ y18_alcohol + y18_smoking + y18_druguse_recode) 
densityplot(bygroup_15cte_ext, data = ~ joverall_pa + jpercultraprocessadoscal + jmvpab5) 


### Distribution of observed and imputed values - histograms ####
par(mfrow=c(1,3))
hist(bygroup_15cte_ext$data$y15_externalising_transformed, xlab="y15_externalising_transformed",
     main="Observed values")
abline(v=mean(bygroup_15cte_ext$data$y15_externalising_transformed, na.rm=TRUE),
       lty=2, col="red")
hist(bygroup_15cte_ext$imp$y15_externalising_transformed[[1]], xlab="y15_externalising_transformed",
     main="Imputed values (first imp)")
abline(v=mean(bygroup_15cte_ext$imp$y15_externalising_transformed[[1]]),lty=2, col="red")
hist(unlist(bygroup_15cte_ext$imp$y15_externalising_transformed), xlab="y15_externalising_transformed",
     main="Imputed values (pooled)")
abline(v=mean(unlist(bygroup_15cte_ext$imp$y15_externalising_transformed)),lty=2, col="red")


### Analyse imputed dataset with practice lm models ####
bygroup_est_y15_externalising_trans <- with(bygroup_15cte_ext, lm(y15_externalising_transformed ~ y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday))
#Pool the results
bygroup_est_y15_externalising_trans_pool <- pool(bygroup_est_y15_externalising_trans)
# summarise results
summary(bygroup_est_y15_externalising_trans_pool, conf.int = TRUE)

# now for Ext transformed
bygroup_est_y18_alcohol_trans <- with(bygroup_15cte_ext, glm(y18_alcohol ~ y15_internalising + y15_externalising_transformed + y15_cte_recode + asexo + aethnicity + ae83 + afumott + aescmae + arendtot + birthday, family = binomial))
#Pool the results
bygroup_est_y18_alcohol_trans_pool <- pool(bygroup_est_y18_alcohol_trans)
# summarise results
summary(bygroup_est_y18_alcohol_trans_pool, conf.int = TRUE) # B=-0.08 (p<.001) for transformed ext here


## Monte Carlo error ####
# function from https://thestatsgeek.com/2022/02/23/how-many-imputations-with-mice-assessing-monte-carlo-error-after-multiple-imputation-in-r/
miceMCError <- function(pooledRes) {
  monteCarloSE <- sqrt(pooledRes$pooled$b/pooledRes$m)
  ciLower <- pooledRes$pooled$estimate - qt(0.975,df=pooledRes$m-1)*monteCarloSE
  ciUpper <- pooledRes$pooled$estimate + qt(0.975,df=pooledRes$m-1)*monteCarloSE
  mcTable <- cbind(pooledRes$pooled$estimate, monteCarloSE, ciLower, ciUpper)
  colnames(mcTable) <- c("Estimate", "Monte Carlo SE", "95% CI lower limit", "95% CI upper limit")
  print(mcTable)
  print("Warning: 95% CI only quantifies Monte-Carlo uncertainty!")
}

# alcohol 
bygroup_est_y18_alcohol <- with(bygroup_15cte_ext, glm(y18_alcohol ~ y15_cte_recode, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)

## adjusted alcohol 
bygroup_est_y18_alcohol <- with(bygroup_15cte_ext, glm(y18_alcohol ~ y15_cte_recode + asexo + 
                                               aethnicity + ae83 + afumott + aescmae + 
                                               arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)
rm(bygroup_est_y18_alcohol)
#### all monte carlo errors are less than 10% of the standard errors

# transformed ext
bygroup_est_ext <- with(bygroup_15cte_ext, lm(y15_externalising_transformed ~ y15_cte_recode + asexo + 
                                             aethnicity + ae83 + afumott + aescmae + 
                                             arendtot + birthday))
summary(pool(bygroup_est_ext))
# pool for mc error
pooled <- pool(bygroup_est_ext)
miceMCError(pooled)
rm(bygroup_est_ext)

# alcohol - ext 
bygroup_est_y18_alcohol <- with(bygroup_15cte_ext, glm(y18_alcohol ~ y15_externalising_transformed + y15_internalising + y15_cte_recode + asexo + 
                                                      aethnicity + ae83 + afumott + aescmae + 
                                                      arendtot + birthday, family = binomial))
summary(pool(bygroup_est_y18_alcohol))
# pool for mc error
pooled <- pool(bygroup_est_y18_alcohol)
miceMCError(pooled)
rm(bygroup_est_y18_alcohol)
#### monte carlo errors are less than 10% of the standard errors


#### export to Stata (y15, transformed externalising, 75 imps) #####
# first, convert to long format 
completed_data_75 <- complete(bygroup_15cte_ext, action = "long", include = TRUE)

# rename the .imp and .id columns
completed_data_75_clean <- completed_data_75 %>%
  rename(imp = .imp, id = .id)

str(completed_data_75_clean)

# change to numeric 
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(asexo = as.numeric(as.character(asexo)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(ae83 = as.numeric(as.character(ae83)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(afumott = as.numeric(as.character(afumott)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_alcohol = as.numeric(as.character(y18_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_smoking = as.numeric(as.character(y18_smoking)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y18_druguse_recode = as.numeric(as.character(y18_druguse_recode)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_alcohol = as.numeric(as.character(y11_alcohol)))
completed_data_75_clean <- completed_data_75_clean %>%
  mutate(y11_smoking = as.numeric(as.character(y11_smoking)))

# Check
table(completed_data_75_clean$asexo, useNA = "ifany")
table(completed_data_75_clean$aethnicity, useNA = "ifany")
table(completed_data_75_clean$ae83, useNA = "ifany")
table(completed_data_75_clean$afumott, useNA = "ifany")
table(completed_data_75_clean$y18_alcohol, useNA = "ifany")
table(completed_data_75_clean$y18_smoking, useNA = "ifany")
table(completed_data_75_clean$y18_druguse_recode, useNA = "ifany")
table(completed_data_75_clean$y11_alcohol, useNA = "ifany")
table(completed_data_75_clean$y11_smoking, useNA = "ifany")

# Check 2
str(completed_data_75_clean$asexo)
str(completed_data_75_clean$aethnicity)
str(completed_data_75_clean$ae83)
str(completed_data_75_clean$afumott)
str(completed_data_75_clean$y18_alcohol)
str(completed_data_75_clean$y18_smoking)
str(completed_data_75_clean$y18_druguse_recode)
str(completed_data_75_clean$y11_alcohol)
str(completed_data_75_clean$y11_smoking)

# write to a .dta file
library(haven)
write_dta(completed_data_75_clean, "imputed_data_long_y15_transExt.dta")

