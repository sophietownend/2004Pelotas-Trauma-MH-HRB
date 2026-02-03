******************************************************************
*       Childhood trauma, adolescent mental health, &            *
*   CMD health risk behaviours in the 2004 Pelotas Birth Cohort  *
******************************************************************

**************************
* Complete Case Analyses *
**************************

//AUTHOR: Sophie Townend (st2325@bath.ac.uk)
//about
Stata/BE 18.0 for Mac (Apple Silicon)
Revision 13 Jul 2023
Copyright 1985-2023 StataCorp LLC

// check info on dataset
describe

// VARIABLE DICTIONARY
// Exposure - Trauma ###
* y11_cte_recode: cumulative trauma load up to age 11 (ordinal none/1/2/3+)
* y15_cte_recode: cumulative trauma load up to age 15 (additional analyses)

// Mediators - MH ###
* y15_internalising: internalising problems at age 15 (continuous)
* y15_externalising: externalising problems at age 15 (continuous)

// Outcomes - Health Risk Behaviours ###
* y18_alcohol: current problematic alcohol use at age 18 (binary)
* y18_smoking: current smoking at age 18 (binary)
* y18_druguse_recode: current illicit drug use at age 18 (binary)
* joverall_pa: total physical activity at 18 (Total acceleration Euclidean Norm Minus One; continuous)
* jpercultraprocessadoscal: Proportion from ultra-processed foods consumption relative to the total energy intake at 18 (in calories; continuous)
* jmvpab5: Moderate to vigorous physical activity at 18 (5 min bouted; continous) in sensivity analyses

// Confounders ###
* asexo: adolescent sex (binary, 0=female; 1=male)
* aethnicity: adolescent ethnicity (binary)
* ae83: maternal alcohol consumption during pregnancy (binary)
* afumott: maternal smoking during pregnancy (binary)
* aescmae: maternal education years at birth (continuous)
* arendtot: monthly family income at birth (continous)
* birthday: day of the year that the adolescent was born (from 1 to 365; proxy for cohort birth order; continuous)

// numeric variables
describe y11_cte_recode y15_cte_recode y15_internalising y15_externalising joverall_pa jpercultraprocessadoscal jmvpab5 aescmae arendtot birthday y15_externalising_transformed

// categorical/binary variables
describe y18_alcohol y18_smoking y18_druguse_recode asexo aethnicity ae83 afumott

mean y11_cte_recode
proportion y11_cte_recode
mean y15_cte_recode
proportion y15_cte_recode
proportion asexo
proportion aethnicity
proportion afumott
proportion ae83
mean arendtot
mean aescmae
mean y15_externalising
mean y15_internalising
proportion y18_alcohol
proportion y18_smoking
proportion y18_druguse_recode
mean joverall_pa
mean jmvpab5
mean jpercultraprocessadoscal


// test a few regression analyses to ensure consistent to R output
regress y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
regress y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
logistic y18_druguse_recode y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 


**********************************************************

* MEDIATION ANALYSES - GFORMULA *

* check up to date gformula package
which gformula

* continuous outcome, multiple mediators, no confounders, no XM interaction *

* overall PA - unadjusted
log using "gformula_results_joverall_pa_unadjusted.txt", text replace
gformula joverall_pa y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - unadjusted 
log using "gformula_results_jpercultraprocessadoscal_unadjusted.txt", text replace
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - unadjusted
log using "gformula_results_jmvpab5_unadjusted.txt", text replace
gformula jmvpab5 y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*continuous outcome, multiple mediators, baseline confounders, no XM interaction 

* overall PA - adjusted
log using "gformula_results_joverall_pa_adjusted.txt", text replace
gformula joverall_pa y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - adjusted
log using "gformula_results_jpercultraprocessadoscal_adjusted.txt", text replace
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - adjusted
log using "gformula_results_jmvpab5_adjusted.txt", text replace
gformula jmvpab5 y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*binary outcome, multiple mediators, no baseline confounders, no XM interaction 

* alcohol use at age 18 - unadjusted --- transformed externalising
log using "gformula_results_y18_alcohol_unadjusted_transformedExt.txt", text replace
gformula y18_alcohol y11_cte_recode y15_internalising y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising y15_externalising_transformed, y15_internalising: y11_cte_recode, y15_externalising_transformed: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - unadjusted
log using "gformula_results_y18_smoking_unadjusted.txt", text replace
gformula y18_smoking y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - unadjusted
log using "gformula_results_y18_druguse_unadjusted.txt", text replace
gformula y18_druguse_recode y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close


*binary outcome, multiple mediators, baseline confounders, no XM interaction 

* alcohol use at age 18 - adjusted (using transformed externalising)
log using "gformula_results_y18_alcohol_adjusted_transformedExt.txt", text replace
gformula y18_alcohol y11_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - adjusted
log using "gformula_results_y18_smoking_adjusted.txt", text replace
gformula y18_smoking y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - adjusted
log using "gformula_results_y18_druguse_adjusted.txt", text replace
gformula y18_druguse_recode y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close


* continuous outcome, ONE mediator, no confounders, no XM interaction *

* overall PA - internalising unadjusted
log using "gformula_results_joverall_pa_unadjusted_int.txt", text replace
gformula joverall_pa y11_cte_recode y15_internalising, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator (y15_internalising) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising , y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* overall PA - externalising unadjusted
log using "gformula_results_joverall_pa_unadjusted_ext.txt", text replace
gformula joverall_pa y11_cte_recode y15_externalising, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator (y15_externalising) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - internalising unadjusted 
log using "gformula_results_jpercultraprocessadoscal_unadjusted_int.txt", text replace
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator (y15_internalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - externalising unadjusted 
log using "gformula_results_jpercultraprocessadoscal_unadjusted_ext.txt", text replace
gformula jpercultraprocessadoscal y11_cte_recode y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator (y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - internalising unadjusted
log using "gformula_results_jmvpab5_unadjusted_int.txt", text replace
gformula jmvpab5 y11_cte_recode y15_internalising, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator (y15_internalising) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - externalising unadjusted
log using "gformula_results_jmvpab5_unadjusted_ext.txt", text replace
gformula jmvpab5 y11_cte_recode y15_externalising, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator (y15_externalising) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*continuous outcome, ONE mediator, baseline confounders, no XM interaction 

* overall PA - internalising adjusted
log using "gformula_results_joverall_pa_adjusted_int.txt", text replace
gformula joverall_pa y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising  asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* overall PA - externalising adjusted
log using "gformula_results_joverall_pa_adjusted_ext.txt", text replace
gformula joverall_pa y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - internalising adjusted
log using "gformula_results_jpercultraprocessadoscal_adjusted_int.txt", text replace
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - externalising adjusted
log using "gformula_results_jpercultraprocessadoscal_adjusted_ext.txt", text replace
gformula jpercultraprocessadoscal y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - internalising adjusted
log using "gformula_results_jmvpab5_adjusted_int.txt", text replace
gformula jmvpab5 y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - externalising adjusted
log using "gformula_results_jmvpab5_adjusted_ext.txt", text replace
gformula jmvpab5 y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*binary outcome, ONE mediator, no baseline confounders, no XM interaction 

* alcohol use at age 18 - internalising unadjusted
log using "gformula_results_y18_alcohol_unadjusted_int.txt", text replace
gformula y18_alcohol y11_cte_recode y15_internalising, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* alcohol use at age 18 - externalising unadjusted - TRANSFORMED
log using "gformula_results_y18_alcohol_unadjusted_Ext.txt", text replace
gformula y18_alcohol y11_cte_recode y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_externalising:regress) ///
equations(y18_alcohol: y11_cte_recode y15_externalising_transformed, y15_externalising_transformed: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - internalising unadjusted
log using "gformula_results_y18_smoking_unadjusted_int.txt", text replace
gformula y18_smoking y11_cte_recode y15_internalising, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - externalising unadjusted
log using "gformula_results_y18_smoking_unadjusted_ext.txt", text replace
gformula y18_smoking y11_cte_recode y15_externalising, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_externalising) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - internalising unadjusted
log using "gformula_results_y18_druguse_unadjusted_int.txt", text replace
gformula y18_druguse_recode y11_cte_recode y15_internalising, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - externalising unadjusted
log using "gformula_results_y18_druguse_unadjusted_ext.txt", text replace
gformula y18_druguse_recode y11_cte_recode y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_externalising) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close


*binary outcome, ONE mediator, baseline confounders, no XM interaction 

* alcohol use at age 18 - internalising adjusted
log using "gformula_results_y18_alcohol_adjusted_int.txt", text replace
gformula y18_alcohol y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* alcohol use at age 18 - externalising adjusted - TRANSFORMED
log using "gformula_results_y18_alcohol_adjusted_Ext.txt", text replace
gformula y18_alcohol y11_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - internalising adjusted
log using "gformula_results_y18_smoking_adjusted_int.txt", text replace
gformula y18_smoking y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - externalising adjusted
log using "gformula_results_y18_smoking_adjusted_ext.txt", text replace
gformula y18_smoking y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - internalising adjusted
log using "gformula_results_y18_druguse_adjusted_int.txt", text replace
gformula y18_druguse_recode y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - externalising adjusted
log using "gformula_results_y18_druguse_adjusted_ext.txt", text replace
gformula y18_druguse_recode y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close



***************************************************************
* Rerun mediation analyses for cumulative trauma up to age 15 *

// check info on dataset
describe

// numeric variables
describe y15_cte_recode y15_internalising y15_externalising joverall_pa jpercultraprocessadoscal jmvpab5 aescmae arendtot birthday y15_externalising_transformed

// categorical/binary variables
describe y18_alcohol y18_smoking y18_druguse asexo aethnicity ae83 afumott

//check a few regressions to ensure consistency with R output
logistic y18_alcohol y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
logistic y18_smoking y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
logistic y18_druguse y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
regress joverall_pa y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
regress jpercultraprocessadoscal y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
regress jmvpab5 y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 

**********************************************************

* MEDIATION ANALYSES - GFORMULA *

* check up to date gformula package
which gformula

* 50 bootstrap, 100,000 sim

* continuous outcome, multiple mediators, no confounders, no XM interaction *

* overall PA - unadjusted
log using "y15_cte_gformula_results_joverall_pa_unadjusted.txt", text replace
gformula joverall_pa y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - unadjusted 
log using "y15_cte_gformula_results_jpercultraprocessadoscal_unadjusted.txt", text replace
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - unadjusted
log using "y15_cte_gformula_results_jmvpab5_unadjusted.txt", text replace
gformula jmvpab5 y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*continuous outcome, multiple mediators, baseline confounders, no XM interaction 

* overall PA - adjusted
log using "y15_cte_gformula_results_joverall_pa_adjusted.txt", text replace
gformula joverall_pa y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - adjusted
log using "y15_cte_gformula_results_jpercultraprocessadoscal_adjusted.txt", text replace
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - adjusted
log using "y15_cte_gformula_results_jmvpab5_adjusted.txt", text replace
gformula jmvpab5 y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*binary outcome, multiple mediators, no baseline confounders, no XM interaction 

* alcohol use at age 18 - unadjusted --- transformed ext
log using "y15_cte_gformula_results_y18_alcohol_unadjusted_transformedExt.txt", text replace
gformula y18_alcohol y15_cte_recode y15_internalising y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising y15_externalising_transformed, y15_internalising: y15_cte_recode, y15_externalising_transformed: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - unadjusted
log using "y15_cte_gformula_results_y18_smoking_unadjusted.txt", text replace
gformula y18_smoking y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - unadjusted
log using "y15_cte_gformula_results_y18_druguse_unadjusted.txt", text replace
gformula y18_druguse_recode y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

*binary outcome, multiple mediators, baseline confounders, no XM interaction 

* alcohol use at age 18 - adjusted (transformed externalising)
log using "y15_cte_gformula_results_y18_alcohol_adjusted_transformedExt.txt", text replace
gformula y18_alcohol y15_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - adjusted
log using "y15_cte_gformula_results_y18_smoking_adjusted.txt", text replace
gformula y18_smoking y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - adjusted
log using "y15_cte_gformula_results_y18_druguse_adjusted.txt", text replace
gformula y18_druguse_recode y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close


* continuous outcome, ONE mediator, no confounders, no XM interaction *

* overall PA - internalising unadjusted
log using "y15_cte_gformula_results_joverall_pa_unadjusted_int.txt", text replace
gformula joverall_pa y15_cte_recode y15_internalising, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator (y15_internalising) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising , y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* overall PA - externalising unadjusted
log using "y15_cte_gformula_results_joverall_pa_unadjusted_ext.txt", text replace
gformula joverall_pa y15_cte_recode y15_externalising, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator (y15_externalising) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - internalising unadjusted 
log using "y15_cte_gformula_results_jpercultraprocessadoscal_unadjusted_int.txt", text replace
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator (y15_internalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - externalising unadjusted 
log using "y15_cte_gformula_results_jpercultraprocessadoscal_unadjusted_ext.txt", text replace
gformula jpercultraprocessadoscal y15_cte_recode y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator (y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - internalising unadjusted
log using "y15_cte_gformula_results_jmvpab5_unadjusted_int.txt", text replace
gformula jmvpab5 y15_cte_recode y15_internalising, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator (y15_internalising) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - externalising unadjusted
log using "y15_cte_gformula_results_jmvpab5_unadjusted_ext.txt", text replace
gformula jmvpab5 y15_cte_recode y15_externalising, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator (y15_externalising) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*continuous outcome, ONE mediator, baseline confounders, no XM interaction 

* overall PA - internalising adjusted
log using "y15_cte_gformula_results_joverall_pa_adjusted_int.txt", text replace
gformula joverall_pa y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising  asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* overall PA - externalising adjusted
log using "y15_cte_gformula_results_joverall_pa_adjusted_ext.txt", text replace
gformula joverall_pa y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - internalising adjusted
log using "y15_cte_gformula_results_jpercultraprocessadoscal_adjusted_int.txt", text replace
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* UPF - externalising adjusted
log using "y15_cte_gformula_results_jpercultraprocessadoscal_adjusted_ext.txt", text replace
gformula jpercultraprocessadoscal y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - internalising adjusted
log using "y15_cte_gformula_results_jmvpab5_adjusted_int.txt", text replace
gformula jmvpab5 y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close

* bouted PA (sensitivity) - externalising adjusted
log using "y15_cte_gformula_results_jmvpab5_adjusted_ext.txt", text replace
gformula jmvpab5 y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim
log close


*binary outcome, ONE mediator, no baseline confounders, no XM interaction 

* alcohol use at age 18 - internalising unadjusted
log using "y15_cte_gformula_results_y18_alcohol_unadjusted_int.txt", text replace
gformula y18_alcohol y15_cte_recode y15_internalising, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* alcohol use at age 18 - externalising unadjusted - TRANSFORMED
log using "y15_cte_gformula_results_y18_alcohol_unadjusted_Ext.txt", text replace
gformula y18_alcohol y15_cte_recode y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_externalising:regress) ///
equations(y18_alcohol: y15_cte_recode y15_externalising_transformed, y15_externalising_transformed: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - internalising unadjusted
log using "y15_cte_gformula_results_y18_smoking_unadjusted_int.txt", text replace
gformula y18_smoking y15_cte_recode y15_internalising, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - externalising unadjusted
log using "y15_cte_gformula_results_y18_smoking_unadjusted_ext.txt", text replace
gformula y18_smoking y15_cte_recode y15_externalising, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_externalising) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - internalising unadjusted
log using "y15_cte_gformula_results_y18_druguse_unadjusted_int.txt", text replace
gformula y18_druguse_recode y15_cte_recode y15_internalising, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - externalising unadjusted
log using "y15_cte_gformula_results_y18_druguse_unadjusted_ext.txt", text replace
gformula y18_druguse_recode y15_cte_recode y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_externalising) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close


*binary outcome, ONE mediator, baseline confounders, no XM interaction 

* alcohol use at age 18 - internalising adjusted
log using "y15_cte_gformula_results_y18_alcohol_adjusted_int.txt", text replace
gformula y18_alcohol y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* alcohol use at age 18 - externalising adjusted - TRANSFORMED
log using "y15_cte_gformula_results_y18_alcohol_adjusted_Ext.txt", text replace
gformula y18_alcohol y15_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - internalising adjusted
log using "y15_cte_gformula_results_y18_smoking_adjusted_int.txt", text replace
gformula y18_smoking y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* smoking at age 18 - externalising adjusted
log using "y15_cte_gformula_results_y18_smoking_adjusted_ext.txt", text replace
gformula y18_smoking y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - internalising adjusted
log using "y15_cte_gformula_results_y18_druguse_adjusted_int.txt", text replace
gformula y18_druguse_recode y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close

* drug use at age 18 - externalising adjusted
log using "y15_cte_gformula_results_y18_druguse_adjusted_ext.txt", text replace
gformula y18_druguse_recode y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(100000) minsim logOR
log close




