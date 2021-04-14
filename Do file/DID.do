use https://github.com/scunning1975/mixtape/raw/master/nsw_mixtape.dta, clear
su re78 if treat
gen y1 = r(mean)
su re78 if treat==0
gen y0 = r(mean)
gen ate = y1-y0
su ate
di 6349.144 - 4554.801
* ATE is 1794.34 
drop if treat==0
drop y1 y0 ate
compress


* Reload experimental group data
use https://github.com/scunning1975/mixtape/raw/master/nsw_mixtape.dta, clear
drop if treat==0

* Now merge in the CPS controls from footnote 2 of Table 2 (Dehejia and Wahba 2002)
append using https://github.com/scunning1975/mixtape/raw/master/cps_mixtape.dta
gen agesq=age*age
gen agecube=age*age*age
gen edusq=educ*educ
gen educube=educ*educ*educ
gen u74 = 0 if re74!=.
replace u74 = 1 if re74==0
gen u75 = 0 if re75!=.
replace u75 = 1 if re75==0
gen interaction1 = educ*re74
gen re74sq=re74^2
gen re74cube=re74^3
gen re75sq=re75^2
gen re75cube=re75^3
gen interaction2 = u74*hisp

* Now estimate the propensity score
logit treat age agesq agecube educ edusq marr nodegree black hisp re74 re75 u74 u75 interaction1 
predict pscore

* Checking mean propensity scores for treatment and control groups
su pscore if treat==1, detail
su pscore if treat==0, detail

* Now look at the propensity score distribution for treatment and control groups
histogram pscore, by(treat) binrescale


***Replication w/ OLS
**ols up to quaratic
clear
use "/Users/justinkao/Desktop/The University of Texas at Austin/Causal Inference/Replication 2/GIthub/Data/nsw_mixtape.dta", clear
reg treat age agesq educ edusq marr nodegree black hisp re74 re74sq re75 re75sq u74 u75, robust
predict pscore_ols_1
su pscore_ols_1 if treat==1, detail
su pscore_ols_1 if treat==0, detail
histogram pscore_ols_1, by(treat) binrescale
graph export pscore_ols_1_1.png, replace
* Trimming the propensity score
drop if pscore_ols_1 <= 0.1 
drop if pscore_ols_1 >= 0.9
histogram pscore_ols_1, by(treat) binrescale
graph export pscore_ols_1_2.png, replace

* Manual with non-normalized weights using all the data
drop d1 d0 s1 s0 y1 y0 ht norm
gen d1=treat/pscore_ols_1
gen d0=(1-treat)/(1-pscore_ols_1)
egen s1=sum(d1)
egen s0=sum(d0)
gen y1=treat*re78/pscore_ols_1
gen y0=(1-treat)*re78/(1-pscore_ols_1)
gen ht=y1-y0
* Manual with normalized weights
replace y1=(treat*re78/pscore_ols_1)/(s1/_N)
replace y0=((1-treat)*re78/(1-pscore_ols_1))/(s0/_N)
gen norm=y1-y0
su ht norm
* ATT under non-normalized weights is -$4,047.76
* ATT under normalized weights is -$4967.571


**ols up to cube
clear
use "/Users/justinkao/Desktop/The University of Texas at Austin/Causal Inference/Replication 2/GIthub/Data/nsw_mixtape.dta", clear
reg treat age agesq agecube educ edusq educube marr nodegree black hisp re74 re74sq re74cube re75 re75sq re75cube u74 u75, robust  
predict pscore_ols_2
su pscore_ols_2 if treat==1, detail
su pscore_ols_2 if treat==0, detail
histogram pscore_ols_2, by(treat) binrescale
graph export pscore_ols_2_1.png, replace
* Trimming the propensity score
drop if pscore_ols_2 <= 0.1 
drop if pscore_ols_2 >= 0.9
histogram pscore_ols_2, by(treat) binrescale
graph export pscore_ols_2_2.png, replace

* Manual with non-normalized weights using all the data
drop d1 d0 s1 s0 y1 y0 ht norm
gen d1=treat/pscore_ols_2
gen d0=(1-treat)/(1-pscore_ols_2)
egen s1=sum(d1)
egen s0=sum(d0)
gen y1=treat*re78/pscore_ols_2
gen y0=(1-treat)*re78/(1-pscore_ols_2)
gen ht=y1-y0
* Manual with normalized weights
replace y1=(treat*re78/pscore_ols_2)/(s1/_N)
replace y0=((1-treat)*re78/(1-pscore_ols_2))/(s0/_N)
gen norm=y1-y0
su ht norm
* ATT under non-normalized weights is -$2254.879
* ATT under normalized weights is -$3918.454

*logit up to quaratic
clear
use "/Users/justinkao/Desktop/The University of Texas at Austin/Causal Inference/Replication 2/GIthub/Data/nsw_mixtape.dta", clear
logit treat age agesq educ edusq marr nodegree black hisp re74 re74sq re75 re75sq u74 u75, robust
predict pscore_logit_1
su pscore_logit_1 if treat==1, detail
su pscore_logit_1 if treat==0, detail
histogram pscore_logit_1, by(treat) binrescale
graph export pscore_logit_1_1.png, replace
* Trimming the propensity score
drop if pscore_logit_1 <= 0.1 
drop if pscore_logit_1 >= 0.9
histogram pscore_logit_1, by(treat) binrescale
graph export pscore_logit_1_2.png, replace

* Manual with non-normalized weights using all the data
drop d1 d0 s1 s0 y1 y0 ht norm
gen d1=treat/pscore_logit_1
gen d0=(1-treat)/(1-pscore_logit_1)
egen s1=sum(d1)
egen s0=sum(d0)
gen y1=treat*re78/pscore_logit_1
gen y0=(1-treat)*re78/(1-pscore_logit_1)
gen ht=y1-y0
* Manual with normalized weights
replace y1=(treat*re78/pscore_logit_1)/(s1/_N)
replace y0=((1-treat)*re78/(1-pscore_logit_1))/(s0/_N)
gen norm=y1-y0
su ht norm
* ATT under non-normalized weights is $1101.713
* ATT under normalized weights is $1246.698

*logit up to cube
clear
use "/Users/justinkao/Desktop/The University of Texas at Austin/Causal Inference/Replication 2/GIthub/Data/nsw_mixtape.dta", clear
logit treat age agesq agecube educ edusq educube marr nodegree black hisp re74 re74sq re74cube re75 re75sq re75cube u74 u75, robust  
predict pscore_logit_2
su pscore_logit_2 if treat==1, detail
su pscore_logit_2 if treat==0, detail
histogram pscore_logit_2, by(treat) binrescale
graph export pscore_logit_2_1.png, replace
* Trimming the propensity score
drop if pscore_logit_2 <= 0.1 
drop if pscore_logit_2 >= 0.9
histogram pscore_logit_2, by(treat) binrescale
graph export pscore_logit_2_2.png, replace

* Manual with non-normalized weights using trimmed data
gen d1=treat/pscore_logit_2
gen d0=(1-treat)/(1-pscore_logit_2)
egen s1=sum(d1)
egen s0=sum(d0)
gen y1=treat*re78/pscore_logit_2
gen y0=(1-treat)*re78/(1-pscore_logit_2)
gen ht=y1-y0
* Manual with normalized weights using trimmed data
replace y1=(treat*re78/pscore_logit_2)/(s1/_N)
replace y0=((1-treat)*re78/(1-pscore_logit_2))/(s0/_N)
gen norm=y1-y0
su ht norm
* ATT under non-normalized weights is $1,744.356
* ATT under normalized weights is $1,635.114




























