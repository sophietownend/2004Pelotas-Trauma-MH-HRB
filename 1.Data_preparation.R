#-------------------------------------------------------------------------------
# PELOTAS ANALYSES - TRAUMA-MH-CMD
#-------------------------------------------------------------------------------
# Dataset preparation 
#-------------------------------------------------------------------------------
# Sophie Townend; 02/04/2025
#-------------------------------------------------------------------------------

library(dplyr)

###### Variable Dictionary ###### 

### Exposure - Trauma ###
# y11_cte_recode: cumulative trauma load up to age 6 (ordinal none/1/2/3+)
# y15_cte_recode: cumulative trauma load up to age 15 for additional analyses

### Mediators - MH ###
# hp1emotion: emotion problems at 15 (continuous)
# hp1conduct: conduct problems at 15 (continuous)
# hp1hyper: hyperactivity at 15 (continuous)
# hp1peer: peer problems at 15 (continuous)
## conduct + hyper = y15_externalising (continuous)
## emotion + peer = y15_internalising (continuous)

### Outcomes - Health Risk Behaviours ###
# y18_alcohol: current problematic alcohol use at age 18 (binary)
# y18_smoking: current smoking at age 18 (binary)
# y18_druguse_recode: current illicit drug use at age 18 (binary) 
# joverall_pa: total physical activity at 18 (Total acceleration Euclidean Norm Minus One; continuous)
# jmvpab5: Moderate to vigorous physical activity at 18 (5 min bouted; continous)
# jpercultraprocessadoscal: Proportion from ultra-processed foods consumption relative to the total energy intake at 18 (in calories; continuous)

### Confounders ###
# asexo: adolescent sex (binary)
# aethnicity: adolescent ethnicity (binary)
# ae83: maternal alcohol consumption during pregnancy (binary)
# afumott: maternal smoking during pregnancy (binary)
# aescmae: maternal education years at birth (continuous)
# arendtot: monthly family income at birth (continous)
# birthday: day of the year that the adolescent was born (from 1 to 365; proxy for cohort birth order; continuous)

### Auxiliary variables for multiple imputation ###
# y6_alltrauma: exposure to any trauma up to age 6 (binary yes/no)
# y6_cte: cumulative trauma load up to age 6 (ordinal none/1/2/3+)
# gp1emotion: emotion problems at 11 (continuous)
# gp1conduct: conduct problems at 11 (continuous)
# gp1hyper: hyperactivity at 11 (continuous)
# gp1peer: peer problems at 11 (continuous)
# pemotion: emotion problems at 18 (continuous)
# pconduct: conduct problems at 18 (continuous)
# phyper: hyperactivity at 18 (continuous)
# ppeer: peer problems at 18 (continuous)
# goverall_pa: total physical activity at 11 (Total acceleration Euclidean Norm Minus One; continuous)
# gmvpab5: Moderate to vigorous physical activity at 18 (5 min bouted; continous)
# gpercultraprocessadosg: Proportion from ultra-processed foods consumption relative to the total energy intake at 11 (in calories; continuous)

###### Dataset Preparation ###### 
dat <- read.csv("trauma-MH-CMD.csv", na.strings = c(''))

# Exclude 2 pps with no data
# from N=4231 to 4229
dat <- dat[-c(which(dat$merge == 'using only (2)')),]

### recode ordinal trauma variables --------------------------------------------
# y6_cte, y11_cte and y15_cte
str(dat$y6_cte)
str(dat$y11_cte)
str(dat$y15_cte)

table(dat$y11_cte)
dat$y11_cte_recode <- recode_factor(dat$y11_cte, 'No trauma exposure'='0', '1 trauma exposure'='1', '2 trauma exposures'='2', '3 or more trauma exposures'='3')
table(dat$y11_cte_recode)

table(dat$y15_cte) # auxiliary
dat$y15_cte_recode <- recode_factor(dat$y15_cte, 'No trauma exposure'='0', '1 trauma exposure'='1', '2 trauma exposures'='2', '3 or more trauma exposures'='3')
table(dat$y15_cte_recode)

table(dat$y6_cte) # auxiliary
dat$y6_cte_recode <- recode_factor(dat$y6_cte, 'No trauma exposure'='0', '1 trauma exposure'='1', '2 trauma exposures'='2', '3 or more trauma exposures'='3')
table(dat$y6_cte_recode)

### create SDQ internalising/externalising scores - y15 ------------------------
psych::describe(dat$hp1emotion)
psych::describe(dat$hp1peer)
psych::describe(dat$hp1conduct)
psych::describe(dat$hp1hyper)

# y15_internalising
dat <- dat %>%
  mutate(y15_internalising = hp1emotion + hp1peer)
psych::describe(dat$y15_internalising)

# y15_externalising
dat <- dat %>%
  mutate(y15_externalising = hp1conduct + hp1hyper)
psych::describe(dat$y15_externalising)

### create SDQ internalising/externalising scores - y11 (auxiliary)  -----------
psych::describe(dat$gp1emotion)
psych::describe(dat$gp1peer)
psych::describe(dat$gp1conduct)
psych::describe(dat$gp1hyper)

# y11_internalising
dat <- dat %>%
  mutate(y11_internalising = gp1emotion + gp1peer)
psych::describe(dat$y11_internalising)

# y11_externalising
dat <- dat %>%
  mutate(y11_externalising = gp1conduct + gp1hyper)
psych::describe(dat$y11_externalising)

### create SDQ internalising/externalising scores - y18 (auxiliary) ------------
psych::describe(dat$pemotion)
psych::describe(dat$ppeer)
psych::describe(dat$pconduct)
psych::describe(dat$phyper)

# y18_internalising
dat <- dat %>%
  mutate(y18_internalising = pemotion + ppeer)
psych::describe(dat$y18_internalising)

# y18_externalising
dat <- dat %>%
  mutate(y18_externalising = pconduct + phyper)
psych::describe(dat$y18_externalising)


### Recode drug use - y18 ------------------------------------------------------ 
# y18_jca35 (marijuana) 
unique(dat$y18_jca35)
dat <- dat %>%
  mutate(
    y18_jca35_ord = case_when(
      y18_jca35 == "Nunca usei" ~ 0,
      y18_jca35 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca35 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca35 == "Uso de vez em quando" ~ 3,
      y18_jca35 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca35 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca35)
table(dat$y18_jca35, dat$y18_jca35_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_marijuana = case_when(
      is.na(y18_jca35_ord) ~ NA_real_,
      y18_jca35_ord < 3    ~ 0,
      y18_jca35_ord > 2    ~ 1
    )
  )
table(dat$y18_jca35_ord)
table(dat$y18_marijuana)

# y18_jca36 (injected cocaine) 
unique(dat$y18_jca36)
dat <- dat %>%
  mutate(
    y18_jca36_ord = case_when(
      is.na(y18_jca36) ~ NA_real_,
      y18_jca36 == "Nunca usei" ~ 0,
      y18_jca36 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca36 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca36 == "Uso de vez em quando" ~ 3,
      y18_jca36 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca36 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca36)
table(dat$y18_jca36, dat$y18_jca36_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_cocaine_injected = case_when(
      is.na(y18_jca36_ord) ~ NA_real_,
      y18_jca36_ord < 3    ~ 0,
      y18_jca36_ord > 2    ~ 1
    )
  )
table(dat$y18_jca36_ord)
table(dat$y18_cocaine_injected)

# y18_jca38 (pills to get high) 
unique(dat$y18_jca38)
dat <- dat %>%
  mutate(
    y18_jca38_ord = case_when(
      is.na(y18_jca38) ~ NA_real_,
      y18_jca38 == "Nunca usei" ~ 0,
      y18_jca38 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca38 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca38 == "Uso de vez em quando" ~ 3,
      y18_jca38 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca38 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca38)
table(dat$y18_jca38, dat$y18_jca38_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_pills = case_when(
      is.na(y18_jca38_ord) ~ NA_real_,
      y18_jca38_ord < 3    ~ 0,
      y18_jca38_ord > 2    ~ 1
    )
  )
table(dat$y18_jca38_ord)
table(dat$y18_pills)

# y18_jca39 (sniffed cocaine) 
unique(dat$y18_jca39)
dat <- dat %>%
  mutate(
    y18_jca39_ord = case_when(
      is.na(y18_jca39) ~ NA_real_,
      y18_jca39 == "Nunca usei" ~ 0,
      y18_jca39 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca39 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca39 == "Uso de vez em quando" ~ 3,
      y18_jca39 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca39 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca39)
table(dat$y18_jca39, dat$y18_jca39_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_cocaine_sniffed = case_when(
      is.na(y18_jca39_ord) ~ NA_real_,
      y18_jca39_ord < 3    ~ 0,
      y18_jca39_ord > 2    ~ 1
    )
  )
table(dat$y18_jca39_ord)
table(dat$y18_cocaine_sniffed)

# y18_jca40 (inhalants) 
unique(dat$y18_jca40)
dat <- dat %>%
  mutate(
    y18_jca40_ord = case_when(
      is.na(y18_jca40) ~ NA_real_,
      y18_jca40 == "Nunca usei" ~ 0,
      y18_jca40 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca40 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca40 == "Uso de vez em quando" ~ 3,
      y18_jca40 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca40 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca40)
table(dat$y18_jca40, dat$y18_jca40_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_inhalants = case_when(
      is.na(y18_jca40_ord) ~ NA_real_,
      y18_jca40_ord < 3    ~ 0,
      y18_jca40_ord > 2    ~ 1
    )
  )
table(dat$y18_jca40_ord)
table(dat$y18_inhalants)

# y18_jca41 (heroin) 
unique(dat$y18_jca41)
dat <- dat %>%
  mutate(
    y18_jca41_ord = case_when(
      is.na(y18_jca41) ~ NA_real_,
      y18_jca41 == "Nunca usei" ~ 0,
      y18_jca41 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca41 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca41 == "Uso de vez em quando" ~ 3,
      y18_jca41 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca41 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca41)
table(dat$y18_jca41, dat$y18_jca41_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_heroin = case_when(
      is.na(y18_jca41_ord) ~ NA_real_,
      y18_jca41_ord < 3    ~ 0,
      y18_jca41_ord > 2    ~ 1
    )
  )
table(dat$y18_jca41_ord)
table(dat$y18_heroin)

# y18_jca42 (ecstasy) 
unique(dat$y18_jca42)
dat <- dat %>%
  mutate(
    y18_jca42_ord = case_when(
      is.na(y18_jca42) ~ NA_real_,
      y18_jca42 == "Nunca usei" ~ 0,
      y18_jca42 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca42 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca42 == "Uso de vez em quando" ~ 3,
      y18_jca42 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca42 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca42)
table(dat$y18_jca42, dat$y18_jca42_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_ecstasy = case_when(
      is.na(y18_jca42_ord) ~ NA_real_,
      y18_jca42_ord < 3    ~ 0,
      y18_jca42_ord > 2    ~ 1
    )
  )
table(dat$y18_jca42_ord)
table(dat$y18_ecstasy)

# y18_jca43 (marijuana with crack) 
unique(dat$y18_jca43)
dat <- dat %>%
  mutate(
    y18_jca43_ord = case_when(
      is.na(y18_jca43) ~ NA_real_,
      y18_jca43 == "Nunca usei" ~ 0,
      y18_jca43 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca43 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca43 == "Uso de vez em quando" ~ 3,
      y18_jca43 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca43 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca43)
table(dat$y18_jca43, dat$y18_jca43_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_marijuana_crack = case_when(
      is.na(y18_jca43_ord) ~ NA_real_,
      y18_jca43_ord < 3    ~ 0,
      y18_jca43_ord > 2    ~ 1
    )
  )
table(dat$y18_jca43_ord)
table(dat$y18_marijuana_crack)


# y18_jca44 (crack) 
unique(dat$y18_jca44)
dat <- dat %>%
  mutate(
    y18_jca44_ord = case_when(
      is.na(y18_jca44) ~ NA_real_,
      y18_jca44 == "Nunca usei" ~ 0,
      y18_jca44 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca44 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca44 == "Uso de vez em quando" ~ 3,
      y18_jca44 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca44 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca44)
table(dat$y18_jca44, dat$y18_jca44_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_crack = case_when(
      is.na(y18_jca44_ord) ~ NA_real_,
      y18_jca44_ord < 3    ~ 0,
      y18_jca44_ord > 2    ~ 1
    )
  )
table(dat$y18_jca44_ord)
table(dat$y18_crack)

# y18_jca45 (LCD/acid - hallucinogens)
unique(dat$y18_jca45)
dat <- dat %>%
  mutate(
    y18_jca45_ord = case_when(
      is.na(y18_jca45) ~ NA_real_,
      y18_jca45 == "Nunca usei" ~ 0,
      y18_jca45 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca45 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca45 == "Uso de vez em quando" ~ 3,
      y18_jca45 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca45 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca45)
table(dat$y18_jca45, dat$y18_jca45_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_hallucinogens = case_when(
      is.na(y18_jca45_ord) ~ NA_real_,
      y18_jca45_ord < 3    ~ 0,
      y18_jca45_ord > 2    ~ 1
    )
  )
table(dat$y18_jca45_ord)
table(dat$y18_hallucinogens)

# y18_jca46 (glue) 
unique(dat$y18_jca46)
dat <- dat %>%
  mutate(
    y18_jca46_ord = case_when(
      is.na(y18_jca46) ~ NA_real_,
      y18_jca46 == "Nunca usei" ~ 0,
      y18_jca46 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca46 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca46 == "Uso de vez em quando" ~ 3,
      y18_jca46 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca46 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca46)
table(dat$y18_jca46, dat$y18_jca46_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_glue = case_when(
      is.na(y18_jca46_ord) ~ NA_real_,
      y18_jca46_ord < 3    ~ 0,
      y18_jca46_ord > 2    ~ 1
    )
  )
table(dat$y18_jca46_ord)
table(dat$y18_glue)

# y18_jca47 (other)
unique(dat$y18_jca47)
dat <- dat %>%
  mutate(
    y18_jca47_ord = case_when(
      is.na(y18_jca47) ~ NA_real_,
      y18_jca47 == "Nunca usei" ~ 0,
      y18_jca47 %in% c("Só experimentei", "S√≥ experimentei") ~ 1,
      y18_jca47 %in% c("Já usei, mas não uso mais", "J√° usei, mas n√£o uso mais") ~ 2,
      y18_jca47 == "Uso de vez em quando" ~ 3,
      y18_jca47 %in% c("Uso só nos finais de semana", "Uso s√≥ nos finais de semana") ~ 4,
      y18_jca47 == "Uso todos os dias, ou quase todos os dias" ~ 5
    )
  )
table(dat$y18_jca47)
table(dat$y18_jca47, dat$y18_jca47_ord, useNA = "ifany")

dat <- dat %>%
  mutate(
    y18_other = case_when(
      is.na(y18_jca47_ord) ~ NA_real_,
      y18_jca47_ord < 3    ~ 0,
      y18_jca47_ord > 2    ~ 1
    )
  )
table(dat$y18_jca47_ord)
table(dat$y18_other)

### all drug use
dat <- dat %>%
  mutate(
    y18_druguse_recode = case_when(
      y18_marijuana == 1 | y18_cocaine_injected == 1 | y18_pills == 1 | y18_cocaine_sniffed == 1 | y18_inhalants == 1 | y18_heroin == 1 | y18_ecstasy == 1 | y18_marijuana_crack == 1 | y18_crack == 1 | y18_hallucinogens == 1 | y18_glue == 1 | y18_other == 1 ~ 1,
      y18_marijuana == 0 & y18_cocaine_injected == 0 & y18_pills == 0 & y18_cocaine_sniffed == 0 & y18_inhalants == 0 & y18_heroin == 0 & y18_ecstasy == 0 & y18_marijuana_crack == 0 & y18_crack == 0 & y18_hallucinogens == 0 & y18_glue == 0 & y18_other == 0 ~ 0,
      is.na(y18_marijuana) | is.na(y18_cocaine_injected) | is.na(y18_pills) | is.na(y18_cocaine_sniffed) | is.na(y18_inhalants) | is.na(y18_heroin) | is.na(y18_ecstasy) | is.na(y18_marijuana_crack) | is.na(y18_crack) | is.na(y18_hallucinogens) | is.na(y18_glue) | is.na(y18_other) ~ NA_real_
    )
  )
table(dat$y18_druguse_recode, useNA = "ifany")


# save this file
write.csv(dat, "trauma-MH-CMD_260126.csv", row.names = F)

