******************************************************************
*       Childhood trauma, adolescent mental health, &            *
*   CMD health risk behaviours in the 2004 Pelotas Birth Cohort  *
******************************************************************

********************
* IMPUTED ANALYSES *
********************

//AUTHOR: Sophie Townend (st2325@bath.ac.uk)
about
Stata/BE 18.0 for Mac (Apple Silicon)
Revision 13 Jul 2023
Copyright 1985-2023 StataCorp LLC


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


**********************************************************************
* Mediation analyses for main model: cumulative trauma up to age 11  *
**********************************************************************

*** set the data as mi:
mi import flong, id(id) m(imp) imputed(y11_cte_recode aethnicity aescmae y15_internalising y15_externalising y18_alcohol y18_smoking y18_druguse_recode joverall_pa jpercultraprocessadoscal jmvpab5 y15_cte_recode y11_internalising y11_externalising y18_internalising y18_externalising y11_alcohol y11_smoking goverall_pa gpercultraprocessados gmvpab5 y6_ctspc y11_ctspc)

*** check Stata interpreted correctly:
mi describe
mi varying

*Check Data*
// Issues with binary coding - check correct
tab asexo
tab y18_alcohol
tab afumott
tab aethnicity
tab y18_druguse_recode

* add labels for interpretation
label define sex_lbl 0 "female" 1 "male"
label values asexo sex_lbl
* check it's worked:
codebook asexo
* add label for interpretation
label define malcohol_lbl 0 "nâo" 1 "sim"
label values ae83 malcohol_lbl
* check it's worked:
codebook ae83
* add label for interpretation
label define msmoke_lbl 0 "No" 1 "Yes"
label values afumott msmoke_lbl
* check it's worked:
codebook afumott
* add label for interpretation
label define alcohol_lbl 0 "No problematic alcohol use" 1 "Problematic alcohol use"
label values y18_alcohol alcohol_lbl
* check it's worked:
codebook y18_alcohol
* add label for interpretation
label define smoking_lbl 0 "Not a current smoker" 1 "Current smoker"
label values y18_smoking smoking_lbl
* check it's worked:
codebook y18_smoking
* recode values 
codebook y18_druguse_recode
recode y18_druguse_recode (1=0) (2=1)
* add label for interpretation
label define drug_lbl 0 "Not a current drug user" 1 "Current drug user"
label values y18_druguse_recode drug_lbl
* check it's worked:
codebook y18_druguse_recode

*** compare results to those obtained in R **** 
*Descriptive statistics in whole sample*
//total sample 
mi estimate: mean y11_cte_recode
mi estimate: proportion y11_cte_recode
mi estimate: proportion asexo
mi estimate: mean aescmae
mi estimate: proportion aethnicity
mi estimate: proportion afumott
mi estimate: proportion ae83
mi estimate: mean arendtot
mi estimate: mean birthday
mi estimate: mean y15_internalising
mi estimate: mean y15_externalising
mi estimate: proportion y18_alcohol
mi estimate: proportion y18_smoking
mi estimate: proportion y18_druguse_recode 
mi estimate: mean joverall_pa
mi estimate: mean jmvpab5
mi estimate: mean jpercultraprocessadoscal

summarize _mi_id _mi_miss _mi_m
summarize _mi_m
tab _mi_m

*Monte Carlo Error - check whether the Monte Carlo error of B is approximately 10 per cent of its standard error (so one value below coefficient vs SE of the coefficient)
mi estimate, mcerror: logistic y18_alcohol y11_cte_recode
mi estimate, mcerror: logistic y18_smoking y11_cte_recode
mi estimate, mcerror: logistic y18_druguse_recode y11_cte_recode
mi estimate, mcerror: regress joverall_pa y11_cte_recode
mi estimate, mcerror: regress jmvpab5 y11_cte_recode
mi estimate, mcerror: regress jpercultraprocessadoscal y11_cte_recode
mi estimate, mcerror: logistic y18_alcohol y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: logistic y18_smoking y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: logistic y18_druguse_recode y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: regress joverall_pa y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: regress jmvpab5 y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: regress jpercultraprocessadoscal y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
*** all consistent with R, and less than 10% of the SE

*** test internalising/externalising correlation ***
mi passive: egen z_y15_internalising = std(y15_internalising)
mi passive: egen z_y15_externalising = std(y15_externalising)
mi estimate: regress z_y15_internalising z_y15_externalising

*** test substance use HRB correlations ***
mi estimate: logit y18_alcohol y18_smoking
mi estimate: logit y18_alcohol y18_druguse_recode
mi estimate: logit y18_smoking y18_druguse_recode


**********************************************************

* REGRESSION ANALYSES *

//unadjusted X-M
mi estimate: regress y15_internalising y11_cte_recode 
mi estimate: regress y15_externalising y11_cte_recode 
//adjusted X-M
mi estimate: regress y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

//unadjusted M-Y - internalising
mi estimate, or: logistic y18_alcohol y15_internalising
mi estimate, or: logistic y18_smoking y15_internalising
mi estimate, or: logistic y18_druguse_recode y15_internalising
mi estimate: regress joverall_pa y15_internalising 
mi estimate: regress jmvpab5 y15_internalising 
mi estimate: regress jpercultraprocessadoscal y15_internalising 
//adjusted M-Y confounds - internalising
mi estimate, or: logistic y18_alcohol y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_smoking y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
//adjusted M-Y confounds + other mediator - internalising
mi estimate, or: logistic y18_alcohol y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_smoking y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  

//unadjusted M-Y - externalising
** no alcohol **
mi estimate, or: logistic y18_smoking y15_externalising
mi estimate, or: logistic y18_druguse_recode y15_externalising
mi estimate: regress joverall_pa y15_externalising 
mi estimate: regress jmvpab5 y15_externalising 
mi estimate: regress jpercultraprocessadoscal y15_externalising 
//adjusted M-Y confounds - externalising
** no alcohol **
mi estimate, or: logistic y18_smoking y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
//adjusted M-Y confounds + other mediator - externalising
** no alcohol **
mi estimate, or: logistic y18_smoking y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday


**********************************************************

* IMPUTED SEX DIFFERENCES *

* X-Y sex interaction - unadjusted *
mi estimate, or: logistic y18_alcohol c.y11_cte_recode##i.asexo 
mi estimate, or: logistic y18_smoking c.y11_cte_recode##i.asexo 
mi estimate, or: logistic y18_druguse_recode c.y11_cte_recode##i.asexo 
mi estimate: regress joverall_pa c.y11_cte_recode##i.asexo
mi estimate: regress jmvpab5 c.y11_cte_recode##i.asexo
mi estimate: regress jpercultraprocessadoscal c.y11_cte_recode##i.asexo 

* X-Y sex interaction - adjusted *
mi estimate, or: logistic y18_alcohol c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* X-Y sex interaction - unadjusted *
mi estimate: regress y15_internalising c.y11_cte_recode##i.asexo
mi estimate: regress y15_externalising c.y11_cte_recode##i.asexo

* X-Y sex interaction - adjusted *
mi estimate: regress y15_internalising c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress y15_externalising c.y11_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* M-Y sex interaction - unadjusted *
* internalising:
mi estimate, or: logistic y18_alcohol c.y15_internalising##i.asexo 
mi estimate, or: logistic y18_smoking c.y15_internalising##i.asexo 
mi estimate, or: logistic y18_druguse_recode c.y15_internalising##i.asexo 
mi estimate: regress joverall_pa c.y15_internalising##i.asexo
mi estimate: regress jmvpab5 c.y15_internalising##i.asexo
mi estimate: regress jpercultraprocessadoscal c.y15_internalising##i.asexo 
* externalising:
mi estimate, or: logistic y18_smoking c.y15_externalising##i.asexo 
mi estimate, or: logistic y18_druguse_recode c.y15_externalising##i.asexo 
mi estimate: regress joverall_pa c.y15_externalising##i.asexo
mi estimate: regress jmvpab5 c.y15_externalising##i.asexo
mi estimate: regress jpercultraprocessadoscal c.y15_externalising##i.asexo 

* M-Y sex interaction - adjusted with confounds *
* internalising:
mi estimate, or: logistic y18_alcohol c.y15_internalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_internalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_internalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_internalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_internalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_internalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
* externalising:
mi estimate, or: logistic y18_smoking c.y15_externalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_externalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_externalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_externalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_externalising##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* M-Y sex interaction - additionally adjusted for other mediator *
* internalising:
mi estimate, or: logistic y18_alcohol c.y15_internalising##i.asexo y15_externalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_internalising##i.asexo y15_externalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_internalising##i.asexo y15_externalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_internalising##i.asexo y15_externalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_internalising##i.asexo y15_externalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_internalising##i.asexo y15_externalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
* externalising:
mi estimate, or: logistic y18_smoking c.y15_externalising##i.asexo y15_internalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_externalising##i.asexo y15_internalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_externalising##i.asexo y15_internalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_externalising##i.asexo y15_internalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_externalising##i.asexo y15_internalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

*significant interaction effects for overall PA
mi estimate: regress joverall_pa y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==0
mi estimate: regress joverall_pa y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==1
*significant interaction effects for bouted PA
mi estimate: regress jmvpab5 y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==0
mi estimate: regress jmvpab5 y15_externalising y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==1

*significant interaction effects for drug use
mi estimate, or: logistic y18_druguse_recode y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==0
mi estimate, or: logistic y18_druguse_recode y15_internalising y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==1


**********************************************************

* MEDIATION ANALYSES *
* final: 50 bootstrap, 100,000 sim

//Use github R code to extract: https://github.com/gemmahammerton/gformula_1993_Pelotas/blob/main/1b%20Stocker%20et%20al%20R%20code.R
** Applies Rubin's rules, creates figures and table


* continuous outcome, multiple mediators, no confounders, no XM interaction *

* overall PA - unadjusted
log using "gformula_imputed_results_joverall_pa_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - unadjusted 
log using "gformula_imputed_results_jpercultraprocessadoscal_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensivity) - unadjusted 
log using "gformula_imputed_results_jmvpab5_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close



*continuous outcome, multiple mediators, baseline confounders, no XM interaction 

* overall PA - adjusted with sex interaction 
log using "gformula_imputed_results_joverall_pa_adjusted_sexint.log", text replace
forvalues impdata = 1/75 { 
				preserve
				di "Imp number is:  " `impdata'
				di "`impdata'"
				keep if _mi_m==`impdata'
* Create the interaction term
    gen sexXext = asexo * y15_externalising
				
gformula joverall_pa y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
derived(sexXext) derrules(sexXext:asexo*y15_externalising) ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - adjusted 
log using "gformula_imputed_results_jpercultraprocessadoscal_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensivity) - adjusted with sex interaction 
log using "gformula_imputed_results_jmvpab5_adjusted_sexint.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
* Create the interaction term
    gen sexXext = asexo * y15_externalising
	
gformula jmvpab5 y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
derived(sexXext) derrules(sexXext:asexo*y15_externalising) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


*binary outcome, multiple mediators, no baseline confounders, no XM interaction 
/// Note - no alcohol model here - need other imputation model

* smoking at age 18 - unadjusted 
log using "gformula_imputed_results_smoking_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - unadjusted 
log using "gformula_imputed_results_druguse_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y11_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising y15_externalising, y15_internalising: y11_cte_recode, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, multiple mediators, baseline confounders, no XM interaction 
/// Note - no alcohol model here - need other imputation model

* smoking at age 18 - adjusted 
log using "gformula_imputed_results_smoking_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - adjusted with sex interaction *** note- internalising, not trauma
log using "gformula_imputed_results_druguse_adjusted_sexint.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'

	* Create the interaction term
    gen sexXint = asexo * y15_internalising
	
gformula y18_druguse_recode y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXint, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXint, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
derived(sexXint) derrules(sexXint:asexo*y15_internalising) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close



* continuous outcome, ONE mediator, no confounders, no XM interaction *

* overall PA - internalising unadjusted 
log using "gformula_imputed_results_joverall_pa_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y11_cte_recode y15_internalising, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator (y15_internalising) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising , y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* overall PA - externalising unadjusted 
log using "gformula_imputed_results_joverall_pa_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y11_cte_recode y15_externalising, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator (y15_externalising) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - internalising unadjusted 
log using "gformula_imputed_results_jpercultraprocessadoscal_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator (y15_internalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - externalising unadjusted 
log using "gformula_imputed_results_jpercultraprocessadoscal_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y11_cte_recode y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator (y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - internalising unadjusted 
log using "gformula_imputed_results_jmvpab5_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y11_cte_recode y15_internalising, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator (y15_internalising) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - externalising unadjusted 
log using "gformula_imputed_results_jmvpab5_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y11_cte_recode y15_externalising, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator (y15_externalising) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


*continuous outcome, ONE mediator, baseline confounders, no XM interaction 

* overall PA - internalising adjusted
log using "gformula_imputed_results_joverall_pa_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_internalising  asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close

* overall PA - externalising adjusted with sex interaction
log using "gformula_imputed_results_joverall_pa_adjusted_ext_sexint.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
	* Create the interaction term
    gen sexXext = asexo * y15_externalising
	
gformula joverall_pa y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, ///
mediation outcome(joverall_pa) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
derived(sexXext) derrules(sexXext:asexo*y15_externalising) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - internalising adjusted
log using "gformula_imputed_results_jpercultraprocessadoscal_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - externalising adjusted
log using "gformula_imputed_results_jpercultraprocessadoscal_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - internalising adjusted
log using "gformula_imputed_results_jmvpab5_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - externalising adjusted with sex interaction
log using "gformula_imputed_results_jmvpab5_adjusted_ext_sexint.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
	* Create the interaction term
    gen sexXext = asexo * y15_externalising
	
gformula jmvpab5 y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, ///
mediation outcome(jmvpab5) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXext, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
derived(sexXext) derrules(sexXext:asexo*y15_externalising) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close



*binary outcome, ONE mediator, no baseline confounders, no XM interaction 
/// Note - no alcohol model here for externalising - need other imputation model

* alcohol at age 18 - internalising unadjusted
log using "gformula_imputed_results_alcohol_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y11_cte_recode y15_internalising, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - internalising unadjusted
log using "gformula_imputed_results_smoking_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y11_cte_recode y15_internalising, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - externalising unadjusted
log using "gformula_imputed_results_smoking_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y11_cte_recode y15_externalising, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_externalising) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - internalising unadjusted
log using "gformula_imputed_results_druguse_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y11_cte_recode y15_internalising, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising, y15_internalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close

* drug use at age 18 - externalising unadjusted
log using "gformula_imputed_results_druguse_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y11_cte_recode y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_externalising) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_externalising, y15_externalising: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, ONE mediator, baseline confounders, no XM interaction 
/// Note - no alcohol model here for externalising - need other imputation model

* alcohol at age 18 - internalising adjusted
log using "gformula_imputed_results_alcohol_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - internalising adjusted
log using "gformula_imputed_results_smoking_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - externalising adjusted
log using "gformula_imputed_results_smoking_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - internalising adjusted with sex interaction
log using "gformula_imputed_results_druguse_adjusted_int_sexint.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'

	* Create the interaction term
    gen sexXint = asexo * y15_internalising
	
gformula y18_druguse_recode y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXint, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday sexXint, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
derived(sexXint) derrules(sexXint:asexo*y15_internalising) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - externalising adjusted
log using "gformula_imputed_results_druguse_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y11_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y11_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close



****************************************************************
*    Mediation analyses for transformed ext - alcohol model    *
****************************************************************

*** set the data as mi:
mi import flong, id(id) m(imp) imputed(y11_cte_recode aethnicity aescmae y15_internalising y15_externalising y18_alcohol y18_smoking y18_druguse_recode joverall_pa jpercultraprocessadoscal jmvpab5 y15_cte_recode y11_internalising y11_externalising y18_internalising y18_externalising y11_alcohol y11_smoking goverall_pa gpercultraprocessados gmvpab5 y6_ctspc y11_ctspc y15_externalising_transformed)

*** check Stata interpreted correctly:
mi describe
mi varying

*Check Data*
describe y11_cte_recode y15_cte_recode y15_internalising y15_externalising joverall_pa jpercultraprocessadoscal jmvpab5 aescmae arendtot birthday y15_externalising_transformed

// Issues with binary coding - check correct
tab asexo
tab y18_alcohol
tab afumott
tab aethnicity
tab y18_druguse_recode

* add labels for interpretation
label define sex_lbl 0 "female" 1 "male"
label values asexo sex_lbl
* check it's worked:
codebook asexo
* add label for interpretation
label define malcohol_lbl 0 "nâo" 1 "sim"
label values ae83 malcohol_lbl
* check it's worked:
codebook ae83
* add label for interpretation
label define msmoke_lbl 0 "No" 1 "Yes"
label values afumott msmoke_lbl
* check it's worked:
codebook afumott
* add label for interpretation
label define alcohol_lbl 0 "No problematic alcohol use" 1 "Problematic alcohol use"
label values y18_alcohol alcohol_lbl
* check it's worked:
codebook y18_alcohol
* add label for interpretation
label define smoking_lbl 0 "Not a current smoker" 1 "Current smoker"
label values y18_smoking smoking_lbl
* check it's worked:
codebook y18_smoking
* recode values 
codebook y18_druguse_recode
recode y18_druguse_recode (1=0) (2=1)
* add label for interpretation
label define drug_lbl 0 "Not a current drug user" 1 "Current drug user"
label values y18_druguse_recode drug_lbl
* check it's worked:
codebook y18_druguse_recode

*** compare results to those obtained in R **** 
*Descriptive statistics in whole sample*
//total sample 
mi estimate: mean y11_cte_recode
mi estimate: proportion asexo
mi estimate: mean aescmae
mi estimate: proportion aethnicity
mi estimate: proportion afumott
mi estimate: proportion ae83
mi estimate: mean arendtot
mi estimate: mean birthday
mi estimate: mean y15_internalising
mi estimate: mean y15_externalising
mi estimate: proportion y18_alcohol
mi estimate: proportion y18_smoking
mi estimate: proportion y18_druguse
mi estimate: mean joverall_pa
mi estimate: mean jmvpab5
mi estimate: mean jpercultraprocessadoscal
mi estimate: mean y15_externalising_transformed


*Monte Carlo Error - check whether the Monte Carlo error of B is approximately 10 per cent of its standard error (so one value below coefficient vs SE of the coefficient)
* X-M
mi estimate, mcerror: regress y15_externalising_transformed y11_cte_recode
* M-Y - unadjusted
mi estimate, mcerror: logistic y18_alcohol y15_externalising_transformed
* M-Y - adjusted for baseline conf
mi estimate, mcerror: logistic y18_alcohol y15_externalising_transformed y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
* M-Y - additionally adjusted for int 
mi estimate, mcerror: logistic y18_alcohol y15_externalising_transformed y15_internalising y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
*** all consistent with R, and less than 10% of the SE


**********************************************************

* REGRESSION ANALYSES *

//unadjusted X-M (as test)
mi estimate: regress y15_internalising y11_cte_recode 
mi estimate: regress y15_externalising y11_cte_recode 
mi estimate: regress y15_externalising_transformed y11_cte_recode 

//adjusted X-M (as test)
mi estimate: regress y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress y15_externalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress y15_externalising_transformed y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

//unadjusted M-Y - externalising TRANSFORMED
mi estimate, or: logistic y18_alcohol y15_externalising_transformed
//adjusted M-Y confounds 
mi estimate, or: logistic y18_alcohol y15_externalising_transformed y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
//adjusted M-Y confounds + mediator 
mi estimate, or: logistic y18_alcohol y15_externalising_transformed y15_internalising y11_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 


**********************************************************

* IMPUTED SEX DIFFERENCES *

* M-Y sex interaction - unadjusted *
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising_transformed##i.asexo 

* M-Y sex interaction - adjusted with confounds *
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising_transformed##i.asexo y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* M-Y sex interaction - additionally adjusted for other mediator *
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising_transformed##i.asexo y15_internalising y11_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday


**********************************************************

* MEDIATION ANALYSES - transformed externalising - alcohol *
* final: 50 bootstrap, 100,000 sim

*binary outcome, multiple mediators, no baseline confounders, no XM interaction 

* alcohol at age 18 - unadjusted 
log using "gformula_imputed_results_alcohol_transformedExt_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y11_cte_recode y15_internalising y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising y15_externalising_transformed, y15_internalising: y11_cte_recode, y15_externalising_transformed: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, multiple mediators, baseline confounders, no XM interaction 

* alcohol at age 18 - adjusted 
log using "gformula_imputed_results_alcohol_transformedExt_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y11_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_internalising y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close



*binary outcome, ONE mediator, no baseline confounders, no XM interaction 

* alcohol at age 18 - transformed externalising unadjusted
log using "gformula_imputed_results_alcohol_unadjusted_transformedExt.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y11_cte_recode y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_externalising_transformed, y15_externalising_transformed: y11_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, ONE mediator, baseline confounders, no XM interaction 

* alcohol at age 18 - transformed externalising adjusted
log using "gformula_imputed_results_alcohol_adjusted_transformedExt.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y11_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y11_cte_recode) mediator(y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y11_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y11_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


****************************************************************
*    Mediation analyses for cumulative trauma up to age 15     *
****************************************************************

*** set the data as mi:
mi import flong, id(id) m(imp) imputed(y11_cte_recode aethnicity aescmae y15_internalising y15_externalising y18_alcohol y18_smoking y18_druguse_recode joverall_pa jpercultraprocessadoscal jmvpab5 y15_cte_recode y11_internalising y11_externalising y18_internalising y18_externalising y11_alcohol y11_smoking goverall_pa gpercultraprocessados gmvpab5 y6_ctspc y11_ctspc)

*** check Stata interpreted correctly:
mi describe
mi varying

*Check Data*
// Issues with binary coding - check correct
tab asexo
tab y18_alcohol
tab afumott
tab aethnicity
tab y18_druguse_recode

* add labels for interpretation
label define sex_lbl 0 "female" 1 "male"
label values asexo sex_lbl
* check it's worked:
codebook asexo
* add label for interpretation
label define malcohol_lbl 0 "nâo" 1 "sim"
label values ae83 malcohol_lbl
* check it's worked:
codebook ae83
* add label for interpretation
label define msmoke_lbl 0 "No" 1 "Yes"
label values afumott msmoke_lbl
* check it's worked:
codebook afumott
* add label for interpretation
label define alcohol_lbl 0 "No problematic alcohol use" 1 "Problematic alcohol use"
label values y18_alcohol alcohol_lbl
* check it's worked:
codebook y18_alcohol
* add label for interpretation
label define smoking_lbl 0 "Not a current smoker" 1 "Current smoker"
label values y18_smoking smoking_lbl
* check it's worked:
codebook y18_smoking
* recode values 
codebook y18_druguse_recode
recode y18_druguse_recode (1=0) (2=1)
* add label for interpretation
label define drug_lbl 0 "Not a current drug user" 1 "Current drug user"
label values y18_druguse_recode drug_lbl
* check it's worked:
codebook y18_druguse_recode

*** compare results to those obtained in R **** 
*Descriptive statistics in whole sample*
//total sample 
mi estimate: mean y15_cte_recode
mi estimate: proportion y15_cte_recode
mi estimate: proportion asexo
mi estimate: mean aescmae
mi estimate: proportion aethnicity
mi estimate: proportion afumott
mi estimate: proportion ae83
mi estimate: mean arendtot
mi estimate: mean birthday
mi estimate: mean y15_internalising
mi estimate: mean y15_externalising
mi estimate: proportion y18_alcohol
mi estimate: proportion y18_smoking
mi estimate: proportion y18_druguse_recode
mi estimate: mean joverall_pa
mi estimate: mean jmvpab5
mi estimate: mean jpercultraprocessadoscal

*Monte Carlo Error - check whether the Monte Carlo error of B is approximately 10 per cent of its standard error (so one value below coefficient vs SE of the coefficient)
mi estimate, mcerror: logistic y18_alcohol y15_cte_recode
mi estimate, mcerror: logistic y18_smoking y15_cte_recode
mi estimate, mcerror: logistic y18_druguse_recode y15_cte_recode
mi estimate, mcerror: regress joverall_pa y15_cte_recode
mi estimate, mcerror: regress jmvpab5 y15_cte_recode
mi estimate, mcerror: regress jpercultraprocessadoscal y15_cte_recode
mi estimate, mcerror: logistic y18_alcohol y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: logistic y18_smoking y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: logistic y18_druguse_recode y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: regress joverall_pa y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: regress jmvpab5 y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
mi estimate, mcerror: regress jpercultraprocessadoscal y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
*** all less than 10% of the SE


**********************************************************

* REGRESSION ANALYSES - age 15 *

//unadjusted X-Y (test)
mi estimate, or: logistic y18_alcohol y15_cte_recode
mi estimate, or: logistic y18_smoking y15_cte_recode
mi estimate, or: logistic y18_druguse_recode y15_cte_recode
mi estimate: regress joverall_pa y15_cte_recode 
mi estimate: regress jmvpab5 y15_cte_recode 
mi estimate: regress jpercultraprocessadoscal y15_cte_recode 
//adjusted X-Y (test)
mi estimate, or: logistic y18_alcohol y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

//unadjusted X-M
mi estimate: regress y15_internalising y15_cte_recode 
mi estimate: regress y15_externalising y15_cte_recode 
//adjusted X-M
mi estimate: regress y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

//unadjusted M-Y - internalising
mi estimate, or: logistic y18_alcohol y15_internalising
mi estimate, or: logistic y18_smoking y15_internalising
mi estimate, or: logistic y18_druguse_recode y15_internalising
mi estimate: regress joverall_pa y15_internalising 
mi estimate: regress jmvpab5 y15_internalising 
mi estimate: regress jpercultraprocessadoscal y15_internalising 
//adjusted M-Y confounds - internalising
mi estimate, or: logistic y18_alcohol y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_smoking y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
//adjusted M-Y confounds + other mediator - internalising
mi estimate, or: logistic y18_alcohol y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_smoking y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  

//unadjusted M-Y - externalising
** no alcohol **
mi estimate, or: logistic y18_smoking y15_externalising
mi estimate, or: logistic y18_druguse_recode y15_externalising
mi estimate: regress joverall_pa y15_externalising 
mi estimate: regress jmvpab5 y15_externalising 
mi estimate: regress jpercultraprocessadoscal y15_externalising 
//adjusted M-Y confounds - externalising
** no alcohol **
mi estimate, or: logistic y18_smoking y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
//adjusted M-Y confounds + other mediator - externalising
** no alcohol **
mi estimate, or: logistic y18_smoking y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate, or: logistic y18_druguse_recode y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress joverall_pa y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jmvpab5 y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday  
mi estimate: regress jpercultraprocessadoscal y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday


**********************************************************

* IMPUTED SEX DIFFERENCES - age 15 *

* X-Y sex interaction - unadjusted *
mi estimate, or: logistic y18_alcohol c.y15_cte_recode##i.asexo 
mi estimate, or: logistic y18_smoking c.y15_cte_recode##i.asexo 
mi estimate, or: logistic y18_druguse_recode c.y15_cte_recode##i.asexo 
mi estimate: regress joverall_pa c.y15_cte_recode##i.asexo
mi estimate: regress jmvpab5 c.y15_cte_recode##i.asexo
mi estimate: regress jpercultraprocessadoscal c.y15_cte_recode##i.asexo 

* X-Y sex interaction - adjusted *
mi estimate, or: logistic y18_alcohol c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* X-M sex interaction - unadjusted *
mi estimate: regress y15_internalising c.y15_cte_recode##i.asexo
mi estimate: regress y15_externalising c.y15_cte_recode##i.asexo

* X-M sex interaction - adjusted *
mi estimate: regress y15_internalising c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress y15_externalising c.y15_cte_recode##i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* M-Y sex interaction - unadjusted *
* internalising:
mi estimate, or: logistic y18_alcohol c.y15_internalising##i.asexo 
mi estimate, or: logistic y18_smoking c.y15_internalising##i.asexo 
mi estimate, or: logistic y18_druguse_recode c.y15_internalising##i.asexo 
mi estimate: regress joverall_pa c.y15_internalising##i.asexo
mi estimate: regress jmvpab5 c.y15_internalising##i.asexo
mi estimate: regress jpercultraprocessadoscal c.y15_internalising##i.asexo 
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising##i.asexo 
mi estimate, or: logistic y18_smoking c.y15_externalising##i.asexo 
mi estimate, or: logistic y18_druguse_recode c.y15_externalising##i.asexo 
mi estimate: regress joverall_pa c.y15_externalising##i.asexo
mi estimate: regress jmvpab5 c.y15_externalising##i.asexo
mi estimate: regress jpercultraprocessadoscal c.y15_externalising##i.asexo 

* M-Y sex interaction - adjusted with confounds *
* internalising:
mi estimate, or: logistic y18_alcohol c.y15_internalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_internalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_internalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_internalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_internalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_internalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_externalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_externalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_externalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_externalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_externalising##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* M-Y sex interaction - additionally adjusted for other mediator *
* internalising:
mi estimate, or: logistic y18_alcohol c.y15_internalising##i.asexo y15_externalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_internalising##i.asexo y15_externalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_internalising##i.asexo y15_externalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_internalising##i.asexo y15_externalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_internalising##i.asexo y15_externalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_internalising##i.asexo y15_externalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_smoking c.y15_externalising##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate, or: logistic y18_druguse_recode c.y15_externalising##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress joverall_pa c.y15_externalising##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jmvpab5 c.y15_externalising##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress jpercultraprocessadoscal c.y15_externalising##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

*significant interaction effects for overall PA
mi estimate: regress joverall_pa y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==0
mi estimate: regress joverall_pa y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==1
*significant interaction effects for bouted PA
mi estimate: regress jmvpab5 y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==0
mi estimate: regress jmvpab5 y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==1

*significant interaction effects for smoking
mi estimate, or: logistic y18_smoking y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==0
mi estimate, or: logistic y18_smoking y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday if asexo==1



**********************************************************

* MEDIATION ANALYSES - age 15 *

* continuous outcome, multiple mediators, no confounders, no XM interaction *

* overall PA - unadjusted
log using "gformula_imputed_results15_joverall_pa_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - unadjusted 
log using "gformula_imputed_results15_jpercultraprocessadoscal_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensivity) - unadjusted 
log using "gformula_imputed_results15_jmvpab5_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator (y15_internalising y15_externalising) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close



*continuous outcome, multiple mediators, baseline confounders, no XM interaction 

* overall PA - adjusted 
log using "gformula_imputed_results15_joverall_pa_adjusted.log", text replace
			forvalues impdata = 1/75 { 
				preserve
				di "Imp number is:  " `impdata'
				di "`impdata'"
				keep if _mi_m==`impdata'
				
			gformula joverall_pa y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
			mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
			commands(joverall_pa:regress, y15_internalising:regress, y15_externalising:regress) ///
			equations(joverall_pa: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
			linexp  ///
			samples(50) seed(79) moreMC sim(10000) minsim

			return list 
			restore
			}
			log close


* UPF - adjusted 
log using "gformula_imputed_results15_jpercultraprocessadoscal_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensivity) - adjusted 
log using "gformula_imputed_results15_jmvpab5_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


*binary outcome, multiple mediators, no baseline confounders, no XM interaction 

/// Note - no alcohol model here - need other imputation model

* smoking at age 18 - unadjusted 
log using "gformula_imputed_results15_smoking_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - unadjusted 
log using "gformula_imputed_results15_druguse_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y15_cte_recode y15_internalising y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising y15_externalising, y15_internalising: y15_cte_recode, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, multiple mediators, baseline confounders, no XM interaction 
/// Note - no alcohol model here - need other imputation model

* smoking at age 18 - adjusted 
log using "gformula_imputed_results15_smoking_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - adjusted 
log using "gformula_imputed_results15_druguse_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* continuous outcome, ONE mediator, no confounders, no XM interaction *

* overall PA - internalising unadjusted 
log using "gformula_imputed_results15_joverall_pa_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y15_cte_recode y15_internalising, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator (y15_internalising) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising , y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* overall PA - externalising unadjusted 
log using "gformula_imputed_results15_joverall_pa_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y15_cte_recode y15_externalising, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator (y15_externalising) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - internalising unadjusted
log using "gformula_imputed_results15_jpercultraprocessadoscal_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator (y15_internalising) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - externalising unadjusted 
log using "gformula_imputed_results15_jpercultraprocessadoscal_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y15_cte_recode y15_externalising, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator (y15_externalising) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - internalising unadjusted 
log using "gformula_imputed_results15_jmvpab5_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y15_cte_recode y15_internalising, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator (y15_internalising) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - externalising unadjusted 
log using "gformula_imputed_results15_jmvpab5_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y15_cte_recode y15_externalising, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator (y15_externalising) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


*continuous outcome, ONE mediator, baseline confounders, no XM interaction 

* overall PA - internalising adjusted
log using "gformula_imputed_results15_joverall_pa_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_internalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_internalising  asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* overall PA - externalising adjusted
log using "gformula_imputed_results15_joverall_pa_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula joverall_pa y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(joverall_pa) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(joverall_pa:regress, y15_externalising:regress) ///
equations(joverall_pa: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - internalising adjusted
log using "gformula_imputed_results15_jpercultraprocessadoscal_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_internalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* UPF - externalising adjusted
log using "gformula_imputed_results15_jpercultraprocessadoscal_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jpercultraprocessadoscal y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jpercultraprocessadoscal) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jpercultraprocessadoscal:regress, y15_externalising:regress) ///
equations(jpercultraprocessadoscal: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - internalising adjusted
log using "gformula_imputed_results15_jmvpab5_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_internalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


* bouted PA (sensitivity) - externalising adjusted
log using "gformula_imputed_results15_jmvpab5_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula jmvpab5 y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(jmvpab5) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(jmvpab5:regress, y15_externalising:regress) ///
equations(jmvpab5: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim

return list 
restore
}
log close


*binary outcome, ONE mediator, no baseline confounders, no XM interaction 
/// Note - no alcohol model here for externalising - need other imputation model

* alcohol at age 18 - internalising unadjusted
log using "gformula_imputed_results15_alcohol_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y15_cte_recode y15_internalising, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - internalising unadjusted
log using "gformula_imputed_results15_smoking_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y15_cte_recode y15_internalising, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - externalising unadjusted
log using "gformula_imputed_results15_smoking_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y15_cte_recode y15_externalising, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_externalising) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - internalising unadjusted
log using "gformula_imputed_results15_druguse_unadjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y15_cte_recode y15_internalising, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising, y15_internalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - externalising unadjusted
log using "gformula_imputed_results15_druguse_unadjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y15_cte_recode y15_externalising, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_externalising) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_externalising, y15_externalising: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, ONE mediator, baseline confounders, no XM interaction 
/// Note - no alcohol model here for externalising - need other imputation model

* alcohol at age 18 - internalising adjusted
log using "gformula_imputed_results15_alcohol_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - internalising adjusted
log using "gformula_imputed_results15_smoking_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_internalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* smoking at age 18 - externalising adjusted
log using "gformula_imputed_results15_smoking_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_smoking y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_smoking) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_smoking:logit, y15_externalising:regress) ///
equations(y18_smoking: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - internalising adjusted
log using "gformula_imputed_results15_druguse_adjusted_int.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_internalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_internalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_internalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


* drug use at age 18 - externalising adjusted
log using "gformula_imputed_results15_druguse_adjusted_ext.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_druguse_recode y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_druguse_recode) exposure(y15_cte_recode) mediator(y15_externalising) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_druguse_recode:logit, y15_externalising:regress) ///
equations(y18_druguse_recode: y15_cte_recode y15_externalising asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close





****************************************************************
*    Mediation analyses for cumulative trauma up to age 15     *
*    transformed ext - alcohol model                           *
****************************************************************

*** set the data as mi:
mi import flong, id(id) m(imp) imputed(y11_cte_recode aethnicity aescmae y15_internalising y15_externalising y18_alcohol y18_smoking y18_druguse_recode joverall_pa jpercultraprocessadoscal jmvpab5 y15_cte_recode y11_internalising y11_externalising y18_internalising y18_externalising y11_alcohol y11_smoking goverall_pa gpercultraprocessados gmvpab5 y6_ctspc y11_ctspc y15_externalising_transformed)

*** check Stata interpreted correctly:
mi describe
mi varying

*Check Data*
// Issues with binary coding - check correct
tab asexo
tab y18_alcohol
tab afumott
tab aethnicity
tab y18_druguse_recode

* add labels for interpretation
label define sex_lbl 0 "female" 1 "male"
label values asexo sex_lbl
* check it's worked:
codebook asexo
* add label for interpretation
label define malcohol_lbl 0 "nâo" 1 "sim"
label values ae83 malcohol_lbl
* check it's worked:
codebook ae83
* add label for interpretation
label define msmoke_lbl 0 "No" 1 "Yes"
label values afumott msmoke_lbl
* check it's worked:
codebook afumott
* add label for interpretation
label define alcohol_lbl 0 "No problematic alcohol use" 1 "Problematic alcohol use"
label values y18_alcohol alcohol_lbl
* check it's worked:
codebook y18_alcohol
* add label for interpretation
label define smoking_lbl 0 "Not a current smoker" 1 "Current smoker"
label values y18_smoking smoking_lbl
* check it's worked:
codebook y18_smoking
* add label for interpretation
label define drug_lbl 0 "Not a current drug user" 1 "Current drug user"
label values y18_druguse_recode drug_lbl
* check it's worked:
codebook y18_druguse_recode

*** compare results to those obtained in R **** 
*Descriptive statistics in whole sample*
//total sample 
mi estimate: mean y11_cte_recode
mi estimate: mean y15_cte_recode
mi estimate: proportion asexo
mi estimate: mean aescmae
mi estimate: proportion aethnicity
mi estimate: proportion afumott
mi estimate: proportion ae83
mi estimate: mean arendtot
mi estimate: mean birthday
mi estimate: mean y15_internalising
mi estimate: mean y15_externalising
mi estimate: proportion y18_alcohol
mi estimate: proportion y18_smoking
mi estimate: proportion y18_druguse_recode
mi estimate: mean joverall_pa
mi estimate: mean jmvpab5
mi estimate: mean jpercultraprocessadoscal

mi estimate: mean y15_externalising_transformed


*Monte Carlo Error - check whether the Monte Carlo error of B is approximately 10 per cent of its standard error (so one value below coefficient vs SE of the coefficient)
* X-M
mi estimate, mcerror: regress y15_externalising_transformed y15_cte_recode
* M-Y - unadjusted
mi estimate, mcerror: logistic y18_alcohol y15_externalising_transformed
* M-Y - adjusted for baseline conf
mi estimate, mcerror: logistic y18_alcohol y15_externalising_transformed y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
* M-Y - additionally adjusted for int 
mi estimate, mcerror: logistic y18_alcohol y15_externalising_transformed y15_internalising y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday
*** all less than 10% of the SE


**********************************************************

* REGRESSION ANALYSES *

//unadjusted X-Y (test)
mi estimate, or: logistic y18_alcohol y15_cte_recode

//adjusted X-Y (test)
mi estimate, or: logistic y18_alcohol y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday


//unadjusted X-M
mi estimate: regress y15_internalising y15_cte_recode 
mi estimate: regress y15_externalising y15_cte_recode 
mi estimate: regress y15_externalising_transformed y15_cte_recode 

//adjusted X-M
mi estimate: regress y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
mi estimate: regress y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday
mi estimate: regress y15_externalising_transformed y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

//unadjusted M-Y - internalising (test)
mi estimate, or: logistic y18_alcohol y15_internalising
//adjusted M-Y confounds 
mi estimate, or: logistic y18_alcohol y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
//adjusted M-Y confounds + mediator 
mi estimate, or: logistic y18_alcohol y15_internalising y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 


//unadjusted M-Y - externalising (test)
mi estimate, or: logistic y18_alcohol y15_externalising
//adjusted M-Y confounds - externalising
mi estimate, or: logistic y18_alcohol y15_externalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
//adjusted M-Y confounds + mediator - externalising
mi estimate, or: logistic y18_alcohol y15_externalising y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 


//unadjusted M-Y - externalising TRANSFORMED
mi estimate, or: logistic y18_alcohol y15_externalising_transformed
//adjusted M-Y confounds - externalising
mi estimate, or: logistic y18_alcohol y15_externalising_transformed y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 
//adjusted M-Y confounds + mediator - externalising
mi estimate, or: logistic y18_alcohol y15_externalising_transformed y15_internalising y15_cte_recode i.asexo i.aethnicity i.afumott i.ae83 aescmae arendtot birthday 


**********************************************************

* IMPUTED SEX DIFFERENCES *

* M-Y sex interaction - unadjusted *
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising_transformed##i.asexo 

* M-Y sex interaction - adjusted with confounds *
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising_transformed##i.asexo y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday

* M-Y sex interaction - additionally adjusted for other mediator *
* externalising:
mi estimate, or: logistic y18_alcohol c.y15_externalising_transformed##i.asexo y15_internalising y15_cte_recode i.aethnicity i.afumott i.ae83 aescmae arendtot birthday




**********************************************************

* MEDIATION ANALYSES - transformed externalising - alcohol *
* final: 50 bootstrap, 100,000 sim

*binary outcome, multiple mediators, no baseline confounders, no XM interaction 

* alcohol at age 18 - unadjusted 
log using "gformula_imputed_results15_alcohol_transformedExt_unadjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y15_cte_recode y15_internalising y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising y15_externalising_transformed, y15_internalising: y15_cte_recode, y15_externalising_transformed: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, multiple mediators, baseline confounders, no XM interaction 

* alcohol at age 18 - adjusted 
log using "gformula_imputed_results15_alcohol_transformedExt_adjusted.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y15_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_internalising y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_internalising:regress, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_internalising y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_internalising: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close



*binary outcome, ONE mediator, no baseline confounders, no XM interaction 

* alcohol at age 18 - transformed externalising unadjusted
log using "gformula_imputed_results15_alcohol_unadjusted_transformedExt.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y15_cte_recode y15_externalising_transformed, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_externalising_transformed) ///
commands(y18_alcohol:logit, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_externalising_transformed, y15_externalising_transformed: y15_cte_recode) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


*binary outcome, ONE mediator, baseline confounders, no XM interaction 

* alcohol at age 18 - transformed externalising adjusted
log using "gformula_imputed_results15_alcohol_adjusted_transformedExt.log", text replace
forvalues impdata = 1/75 { 
	preserve
	di "Imp number is:  " `impdata'
	di "`impdata'"
	keep if _mi_m==`impdata'
	
gformula y18_alcohol y15_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, ///
mediation outcome(y18_alcohol) exposure(y15_cte_recode) mediator(y15_externalising_transformed) base_confs(asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
commands(y18_alcohol:logit, y15_externalising_transformed:regress) ///
equations(y18_alcohol: y15_cte_recode y15_externalising_transformed asexo aethnicity afumott ae83 aescmae arendtot birthday, y15_externalising_transformed: y15_cte_recode asexo aethnicity afumott ae83 aescmae arendtot birthday) ///
linexp  ///
samples(50) seed(79) moreMC sim(10000) minsim logOR

return list 
restore
}
log close


