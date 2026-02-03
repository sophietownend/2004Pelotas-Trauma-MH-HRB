# 2004Pelotas-Trauma-MH-CMD
This repository contains the analysis code used for data preparation, imputation, and statistical analyses for the manuscript:

> *Mental health problems as mediators between childhood trauma and cardiometabolic risk behaviours in the 2004 Pelotas Birth Cohort*

The study examines whether internalising and externalising mental health problems mediate longitudinal associations between childhood trauma exposure and health risk behaviours (HRBs) in a population-based cohort of Brazilian adolescents.

---

## Overview of analyses

The main analyses include:

1. **Exposure–mediator associations**  
   Linear regression models examining associations between cumulative trauma exposure up to age 11 and mental health problems at age 15.

2. **Mediator–outcome associations**  
   Linear and logistic regression models examining associations between mental health problems at age 15 and health risk behaviours at age 18.

3. **Sex interactions**  
   Sensitivity analyses to assess potential sex differences.

4. **Counterfactual mediation analyses**  
   - Exposure: Cumulative trauma up to age 11  
   - Mediators: Internalising and externalising problems at age 15  
   - Outcomes: Smoking, problematic alcohol use, illicit drug use, overall and bouted physical activity, and ultra-processed food intake at age 18
   - Sensitivity analyses assessing whether indirect effects were driven by a single mediator.

5. **Additional sensitivity analyses**  
   Repetition of all models using cumulative trauma exposure up to age 15 as the exposure.

---

The following confounders were adjusted for in all analyses: adolescent sex, adolescent ethnicity, maternal alcohol consumption during pregnancy, maternal smoking during pregnancy, maternal years of education at birth, monthly family income at birth, and day of the year that the adolescent was born (proxy for cohort birth order). 

---

## Repository contents

The following scripts are included:

1. **`1.Data_preparation.R`**  
   R code for data preparation and variable derivation.

2. **`2.Complete_case_analyses.R`**  
   R code for regression analyses using complete-case data.

3. **`3.Complete_case_mediation.do`**  
   Stata code for counterfactual mediation analyses using complete-case data.

4. **`Multiple_imputation.R`**  
   R code for multivariate imputation by chained equations. Four imputation models were run as described in the Supplementary Materials.

5. **`Imputed_analyses.do`**  
   Stata code for counterfactual mediation analyses using imputed data (N = 4,229).

---

## Data availability

The data used in this study are not publicly available. Applications to use the data can be made by contacting the 2004 Pelotas Birth Cohort coordinators (https://epidemio-ufpel.org.br/corpo-docente-ppgep/). 

---

## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/ 
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png 
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg  
<img width="378" height="146" alt="image" src="https://github.com/user-attachments/assets/27e7d514-c0f4-497d-8be2-5d5b2bed5553" />
