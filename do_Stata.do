//////time-varying cox
///all cause
capture postclose allcause
postfile allcause str30 var hr ll ul p using "F:\蛋白质组与死亡率\数据集和代码\allcause", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdloga0a0c4dh25-stdlogq9y4h2 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post allcause ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcause
use "F:\蛋白质组与死亡率\数据集和代码\allcause.dta", clear

///cancer
capture postclose cancer
postfile cancer str30 var hr ll ul p using "F:\蛋白质组与死亡率\数据集和代码\cancer", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failurecancer)
foreach var of varlist stdloga0a0c4dh25-stdlogq9y4h2 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post cancer ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose cancer
use "F:\蛋白质组与死亡率\数据集和代码\cancer.dta", clear

///cvd
capture postclose cvd
postfile cvd str30 var hr ll ul p using "F:\蛋白质组与死亡率\数据集和代码\cvd", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failurecvd)
foreach var of varlist stdloga0a0c4dh25-stdlogq9y4h2 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post cvd ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose cvd
use "F:\蛋白质组与死亡率\数据集和代码\cvd.dta", clear


////// prospective cox
///all cause
capture postclose allcauseb
postfile allcauseb str30 var hr ll ul p using "F:\蛋白质组与死亡率\数据集和代码\allcauseb", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_b.dta", clear
stset exit, enter(entry) failure(failureall)
foreach var of varlist stdloga0a0c4dh25-stdlogq9y4h2 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post allcauseb ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcauseb
use "F:\蛋白质组与死亡率\数据集和代码\allcauseb.dta", clear

///cancer
capture postclose cancerb
postfile cancerb str30 var hr ll ul p using "F:\蛋白质组与死亡率\数据集和代码\cancerb", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_b.dta", clear
stset exit, enter(entry) failure(failurecancer)
foreach var of varlist stdloga0a0c4dh25-stdlogq9y4h2 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post cancerb ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose cancerb
use "F:\蛋白质组与死亡率\数据集和代码\cancerb.dta", clear

///cvd
capture postclose cvdb
postfile cvdb str30 var hr ll ul p using "F:\蛋白质组与死亡率\数据集和代码\cvdb", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_b.dta", clear
stset exit, enter(entry) failure(failurecvd)
foreach var of varlist stdloga0a0c4dh25-stdlogq9y4h2 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post cvdb ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose cvdb
use "F:\蛋白质组与死亡率\数据集和代码\cvdb.dta", clear



////// sensitivity analysis
/// sensia full adjusted model
capture postclose allcause_sensia
postfile allcause_sensia str30 var hra lla ula pa using "F:\蛋白质组与死亡率\数据集和代码\allcause_sensia", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex bmi i.smoke i.alc met i.hyper i.t2d i.cvd i.cancer, vce(cluster id)
    matrix M = r(table)
    post allcause_sensia ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcause_sensia
use "F:\蛋白质组与死亡率\数据集和代码\allcause_sensia.dta", clear

/// sensib exclude in three years
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
bysort id (entry): gen first_entry = entry[1]
bysort id (entry): gen last_exit  = exit[_N]
gen timedeath = last_exit - first_entry
bys id: egen died = max(failureall)
gen drop_id = (died==1 & timedeath <= 3)
drop if drop_id==1
save "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensib.dta"

capture postclose allcause_sensib
postfile allcause_sensib str30 var hrb llb ulb pb using "F:\蛋白质组与死亡率\数据集和代码\allcause_sensib", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensib.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post allcause_sensib ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcause_sensib
use "F:\蛋白质组与死亡率\数据集和代码\allcause_sensib.dta", clear

/// sensic exclude bmi>30
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
bysort id (time): gen baseline_bmi = bmi[1]
drop if baseline_bmi >= 30
save "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensic.dta"

capture postclose allcause_sensic
postfile allcause_sensic str30 var hrc llc ulc pc using "F:\蛋白质组与死亡率\数据集和代码\allcause_sensic", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensic.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post allcause_sensic ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcause_sensic
use "F:\蛋白质组与死亡率\数据集和代码\allcause_sensic.dta", clear

/// sensid exclude baseline disease
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
bysort id (time): gen base_t2d = t2d[1]
bysort id (time): gen base_cvd  = cvd[1]
bysort id (time): gen base_cancer = cancer[1]
drop if base_t2d==1 | base_cvd==1 | base_cancer==1
save "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensid.dta"

capture postclose allcause_sensid
postfile allcause_sensid str30 var hrd lld uld pd using "F:\蛋白质组与死亡率\数据集和代码\allcause_sensid", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensid.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post allcause_sensid ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcause_sensid
use "F:\蛋白质组与死亡率\数据集和代码\allcause_sensid.dta", clear


///sensie exclude age>70 years
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
bysort id (time): gen base_age = ageadj[1]
drop if base_age >=70
save "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensie.dta"

capture postclose allcause_sensie
postfile allcause_sensie str30 var hre lle ule pe using "F:\蛋白质组与死亡率\数据集和代码\allcause_sensie", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality_sensie.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex, vce(cluster id)
    matrix M = r(table)
    post allcause_sensie ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose allcause_sensie
use "F:\蛋白质组与死亡率\数据集和代码\allcause_sensie.dta", clear


////// sungroup analysis
/// by sex
capture postclose subgroup_male
postfile subgroup_male str30 var hrb llb ulb pb using "F:\蛋白质组与死亡率\数据集和代码\subgroup_male", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex if sex==1, vce(cluster id)
    matrix M = r(table)
    post subgroup_male ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose subgroup_male
use "F:\蛋白质组与死亡率\数据集和代码\subgroup_male.dta", clear

capture postclose subgroup_female
postfile subgroup_female str30 var hrb llb ulb pb using "F:\蛋白质组与死亡率\数据集和代码\subgroup_female", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex if sex==0, vce(cluster id)
    matrix M = r(table)
    post subgroup_female ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose subgroup_female
use "F:\蛋白质组与死亡率\数据集和代码\subgroup_female.dta", clear


use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox c.`var'##i.sexadj ageadj, vce(cluster id)
}


/// by age
capture postclose subgroup_age1
postfile subgroup_age1 str30 var hrb llb ulb pb using "F:\蛋白质组与死亡率\数据集和代码\subgroup_age1", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex if age_baseline<60, vce(cluster id)
    matrix M = r(table)
    post subgroup_age1 ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose subgroup_age1
use "F:\蛋白质组与死亡率\数据集和代码\subgroup_age1.dta", clear

capture postclose subgroup_age2
postfile subgroup_age2 str30 var hrb llb ulb pb using "F:\蛋白质组与死亡率\数据集和代码\subgroup_age2", replace
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox `var' ageadj i.sex if age_baseline>=60, vce(cluster id)
    matrix M = r(table)
    post subgroup_age2 ("`var'") (M[1,1]) (M[5,1]) (M[6,1]) (M[4,1])
}
postclose subgroup_age2
use "F:\蛋白质组与死亡率\数据集和代码\subgroup_age2.dta", clear

use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, id(id) enter(entry) failure(failureall)
foreach var of varlist stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909 {
    stcox c.`var'##c.ageadj i.sexadj, vce(cluster id)
}


////// restricted cubic spline
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
stset exit, enter(time entry) failure(failureall == 1) id(code_id)

local proteins stdlogp10643 stdlogp01700 stdlogp01701 stdlogp01834 stdlogp01714 stdlogp02747 stdlogp06312 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp19652 stdlogq9hc29 stdlogp35579 stdlogq92835 stdlogp04004 stdlogp10909

local outdir "F:\蛋白质组与死亡率\数据集和代码\RCS_results_3knots"
capture mkdir "`outdir'"
tempname posth
tempfile rcs_summary
postfile `posth' str20 protein double k1 k2 k3 xmin xmax p_overall p_nonlinear using `rcs_summary', replace

foreach p of local proteins {
    di as text "=================================================="
    di as result "Running RCS for: `p'"
    di as text "=================================================="
    capture confirm variable `p'
    if _rc {
        di as error "Variable `p' not found. Skipped."
        continue
    }
    capture drop __touse
    capture drop __rcs*
    gen byte __touse = (_st == 1) & !missing(`p', ageadj, sexadj)
    quietly count if __touse
    if r(N) == 0 {
        di as error "No eligible observations for `p'. Skipped."
        capture drop __touse
        capture drop __rcs*
        continue
    }
    di as text "Eligible observations = " r(N)
    quietly _pctile `p' if __touse, p(10 50 90)
    local k1 = r(r1)
    local k2 = r(r2)
    local k3 = r(r3)
    di as text "RCS knots: `k1' `k2' `k3'"
    quietly _pctile `p' if __touse, p(1 99)
    local xmin = r(r1)
    local xmax = r(r2)
    di as text "Prediction range: `xmin' to `xmax'"
    if (`k1' >= `k2' | `k2' >= `k3') {
        di as error "Duplicate or non-increasing knots for `p'. RCS skipped."
        post `posth' ("`p'") (`k1') (`k2') (`k3') ///
            (`xmin') (`xmax') (.) (.)
        capture drop __touse
        capture drop __rcs*
        continue
    }
    mkspline __rcs = `p' if __touse, cubic knots(`k1' `k2' `k3')
    capture noisily stcox __rcs1 __rcs2 ageadj i.sexadj ///
        if __touse, vce(cluster code_id)
    if _rc {
        di as error "Cox model failed for `p'. Skipped."
        capture drop __touse
        capture drop __rcs*
        continue
    }
    quietly test __rcs1 __rcs2
    local p_overall = r(p)
    quietly test __rcs2
    local p_nonlinear = r(p)
    di as text "P overall = " as result %8.6f `p_overall'
    di as text "P nonlinear = " as result %8.6f `p_nonlinear'
    capture estimates drop __rcsmodel
    estimates store __rcsmodel
    preserve
        clear
        set obs 102
        gen double `p' = `xmin' + (`xmax' - `xmin') * (_n - 1) / 100 in 1/101
        replace `p' = 0 in 102
        gen byte isref = (_n == 102)
        mkspline __grcs = `p', cubic knots(`k1' `k2' `k3')
        quietly summarize __grcs1 if isref, meanonly
        local ref1 = r(mean)
        quietly summarize __grcs2 if isref, meanonly
        local ref2 = r(mean)
        drop if isref
        gen double d1 = __grcs1 - `ref1'
        gen double d2 = __grcs2 - `ref2'
        estimates restore __rcsmodel
        scalar __b1 = _b[__rcs1]
        scalar __b2 = _b[__rcs2]
        matrix __V = e(V)
        scalar __v11 = __V["__rcs1","__rcs1"]
        scalar __v22 = __V["__rcs2","__rcs2"]
        scalar __v12 = __V["__rcs1","__rcs2"]
        gen double loghr = d1 * __b1 + d2 * __b2
        gen double se = sqrt(d1^2 * __v11 + d2^2 * __v22 + 2 * d1 * d2 * __v12)
        gen double hr = exp(loghr)
        gen double lb = exp(loghr - 1.96 * se)
        gen double ub = exp(loghr + 1.96 * se)
        gen double p_overall = `p_overall'
        gen double p_nonlinear = `p_nonlinear'
        sort `p'
        local po_text = cond(`p_overall' < 0.001, "<0.001", string(`p_overall', "%5.3f"))
        local pn_text = cond(`p_nonlinear' < 0.001, "<0.001", string(`p_nonlinear', "%5.3f"))
        twoway ///
            (rarea ub lb `p', color(gs12%45) lcolor(none)) ///
            (line hr `p', lcolor(navy) lwidth(medthick)), ///
            yline(1, lpattern(dash) lcolor(gs8)) ///
            xline(0, lpattern(dot) lcolor(gs8)) ///
            ytitle("HR for all-cause mortality") ///
            xtitle("`p'") ///
            note("P overall = `po_text'; P nonlinear = `pn_text'", size(small)) ///
            legend(off) ///
            graphregion(color(white)) ///
            plotregion(color(white)) ///
            name(G_`p', replace)
        graph export "`outdir'\RCS_`p'.pdf", as(pdf) replace
        save "`outdir'\RCS_prediction_`p'.dta", replace
    restore
    post `posth' ("`p'") (`k1') (`k2') (`k3') ///
        (`xmin') (`xmax') (`p_overall') (`p_nonlinear')
    capture estimates drop __rcsmodel
    capture drop __touse
    capture drop __rcs*
}


postclose `posth'
preserve
    use `rcs_summary', clear
    format k1 k2 k3 xmin xmax %9.4f
    format p_overall p_nonlinear %9.6f
    sort p_overall
    list, noobs sep(0)
    save "`outdir'\RCS_summary.dta", replace
    export excel using "`outdir'\RCS_summary.xlsx", firstrow(variables) replace
restore


//////score and mortality
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
gen score_acti=0.463*stdlogp10643+0.43*stdlogp01700+0.423*stdlogp01701+0.376*stdlogp01834+0.359*stdlogp01714+0.354*stdlogp02747+0.348*stdlogp06312-0.172*stdlogp27918
gen score_regu=0.381*stdlogp10909+0.269*stdlogp04004+0.147*stdlogp08603+0.12*stdlogp04003
egen stdscore_acti=std(score_acti)
egen stdscore_regu=std(score_regu)
gen score_imba=stdscore_acti-stdscore_regu
egen stdscore_imba=std(score_imba)
save "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", replace

stset exit, id(id) enter(entry) failure(failureall)
stcox stdscore_imba ageadj i.sex, vce(cluster id)
matrix list r(table)
stcox stdscore_acti ageadj i.sex, vce(cluster id)
matrix list r(table)
stcox stdscore_regu ageadj i.sex, vce(cluster id)
matrix list r(table)

use "F:\蛋白质组与死亡率\数据集和代码\3.pro_mortality_score.dta", clear
replace death=0 if death==.
mixed stdscore_imba i.death##i.time ageadj i.sexadj || id:
testparm i.death
testparm i.time
testparm i.death##i.time

mixed stdscore_acti i.death##i.time ageadj i.sexadj || id:
testparm i.death
testparm i.time
testparm i.death##i.time

mixed stdscore_regu i.death##i.time ageadj i.sexadj || id:
testparm i.death
testparm i.time
testparm i.death##i.time


//////external validation
use "F:\蛋白质组与死亡率\数据集和代码\3.CHNS_pro_mortality.dta", clear

gen score_acti=0.463*stdlogp10643+0.43*stdlogp01700+0.423*stdlogp01701+0.376*stdlogp01834+0.359*stdlogp01714+0.354*stdlogp02747+0.348*stdlogp06312
gen score_regu=0.381*stdlogp10909+0.269*stdlogp04004+0.147*stdlogp08603+0.12*stdlogp04003
egen stdscore_acti=std(score_acti)
egen stdscore_regu=std(score_regu)
gen score_imba=stdscore_acti-stdscore_regu
egen stdscore_imba=std(score_imba)

logit death stdscore_imba age i.sex,or
logit death stdscore_acti age i.sex,or
logit death stdscore_regu age i.sex,or


//////PCA of pathway intensity
use "F:\蛋白质组与死亡率\数据集和代码\2.pro_mortality.dta", clear
pca stdlogp01701 stdlogp01700 stdlogp01714 stdlogp01834 stdlogp06312 stdlogp10643 stdlogp02747 stdlogp04003 stdlogp08603 stdlogp27918 stdlogp04004 stdlogp10909
estat loadings
predict pathway_score1, score

pca stdlogp01701 stdlogp01700 stdlogp01714 stdlogp01834 stdlogp06312 
estat loadings
predict pathway_score2, score

mixed pathway_score1 i.death##i.time || id:
testparm i.death
testparm i.time
testparm i.death##i.time

mixed pathway_score2 i.death##i.time || id:
testparm i.death
testparm i.time
testparm i.death##i.time


//////////// dietary modification
///screening
capture postclose result_dietscore
postfile result_dietscore str30 depvar str30 var coef p using"F:\蛋白质组与死亡率\5.膳食可调节性\result_dietscore", replace
use "F:\蛋白质组与死亡率\5.膳食可调节性\GNHS_score_diet.dta",clear
foreach depvar of varlist stdscore_imba stdscore_acti stdscore_regu {
	foreach var of varlist stdeacarbo-serving_coffee{
	mixed `depvar' `var' ageadj i.sexadj energy ||code_id:
	matrix M = r(table)
    post result_dietscore ("`depvar'") ("`var'") (M[1,1]) (M[4,1])
	}
}
postclose result_dietscore
use "F:\蛋白质组与死亡率\5.膳食可调节性\result_dietscore.dta",clear
sort p

/// score and energy-adjusted carbohydrate
use "F:\蛋白质组与死亡率\5.膳食可调节性\GNHS_score_diet.dta", clear
mixed stdscore_imba stdeacarbo ageadj i.sexadj energy ||code_id:
mixed stdscore_acti stdeacarbo ageadj i.sexadj energy ||code_id: 
mixed stdscore_regu stdeacarbo ageadj i.sexadj energy ||code_id:

/// energysubstitution model
mixed stdscore_imba epcarbo ageadj i.sexadj energy ||code_id:
mixed stdscore_imba epcarbo epprotein ageadj i.sexadj energy ||code_id:
mixed stdscore_imba epcarbo epfat ageadj i.sexadj energy ||code_id:

mixed stdscore_regu epcarbo ageadj i.sexadj energy ||code_id:
mixed stdscore_regu epcarbo epprotein ageadj i.sexadj energy ||code_id:
mixed stdscore_regu epcarbo epfat ageadj i.sexadj energy ||code_id:


///regulatory proteins
mixed stdlogp08603 stdeacarbo ageadj i.sexadj energy ||code_id:
mixed stdlogp04003 stdeacarbo ageadj i.sexadj energy ||code_id:
mixed stdlogp04004 stdeacarbo ageadj i.sexadj energy ||code_id:
mixed stdlogp10909 stdeacarbo ageadj i.sexadj energy ||code_id:


/// risk modification 
use "F:\蛋白质组与死亡率\5.膳食可调节性\GNHS_score_diet.dta",clear
keep idtime stdeacarbo eacarbo epcarbo energy
merge 1:1 idtime using "F:\蛋白质组与死亡率\数据集和代码\3.pro_mortality_score.dta"
keep if _merge==3

bysort id: egen mean_eacarbo = mean(eacarbo)
bysort id: egen mean_epcarbo = mean(epcarbo)

xtile cea=mean_eacarbo,n(3)

gen cep=1 if mean_epcarbo<50
replace cep=3 if mean_epcarbo>65
replace cep=2 if cep==.

gen stdscore_regu_rev = -stdscore_regu
gen stdlogp04004_rev = -stdlogp04004
gen stdlogp10909_rev = -stdlogp10909

stset exit, id(id) enter(entry) failure(failureall)
stcox stdscore_imba ageadj i.sex if cea==1, vce(cluster id)
stcox stdscore_imba ageadj i.sex if cea==2, vce(cluster id)
stcox stdscore_imba ageadj i.sex if cea==3, vce(cluster id)
stcox c.stdscore_imba##c.stdeacarbo ageadj i.sex, vce(cluster id)

stset exit, id(id) enter(entry) failure(failureall)
stcox stdscore_regu_rev ageadj i.sex if cea==1, vce(cluster id)
stcox stdscore_regu_rev ageadj i.sex if cea==2, vce(cluster id)
stcox stdscore_regu_rev ageadj i.sex if cea==3, vce(cluster id)
stcox c.stdscore_regu_rev##c.stdeacarbo ageadj i.sex, vce(cluster id)

stcox stdlogp04004_rev ageadj i.sex if cea==1, vce(cluster id)
stcox stdlogp04004_rev ageadj i.sex if cea==2, vce(cluster id)
stcox stdlogp04004_rev ageadj i.sex if cea==3, vce(cluster id)
stcox c.stdlogp04004_rev##c.stdeacarbo ageadj i.sex, vce(cluster id)

stcox stdlogp10909_rev ageadj i.sex if cea==1, vce(cluster id)
stcox stdlogp10909_rev ageadj i.sex if cea==2, vce(cluster id)
stcox stdlogp10909_rev ageadj i.sex if cea==3, vce(cluster id)
stcox c.stdlogp10909_rev##c.stdeacarbo ageadj i.sex, vce(cluster id)

stset exit, id(id) enter(entry) failure(failureall)
stcox stdscore_imba ageadj i.sex if cep==1, vce(cluster id)
stcox stdscore_imba ageadj i.sex if cep==2, vce(cluster id)
stcox stdscore_imba ageadj i.sex if cep==3, vce(cluster id)
stcox c.stdscore_imba##c.epcarbo ageadj i.sex, vce(cluster id)

stset exit, id(id) enter(entry) failure(failureall)
stcox stdscore_regu_rev ageadj i.sex if cep==1, vce(cluster id)
stcox stdscore_regu_rev ageadj i.sex if cep==2, vce(cluster id)
stcox stdscore_regu_rev ageadj i.sex if cep==3, vce(cluster id)
stcox c.stdscore_regu_rev##c.epcarbo ageadj i.sex, vce(cluster id)

stcox stdlogp04004_rev ageadj i.sex if cep==1, vce(cluster id)
stcox stdlogp04004_rev ageadj i.sex if cep==2, vce(cluster id)
stcox stdlogp04004_rev ageadj i.sex if cep==3, vce(cluster id)
stcox c.stdlogp04004_rev##c.epcarbo ageadj i.sex, vce(cluster id)

stcox stdlogp10909_rev ageadj i.sex if cep==1, vce(cluster id)
stcox stdlogp10909_rev ageadj i.sex if cep==2, vce(cluster id)
stcox stdlogp10909_rev ageadj i.sex if cep==3, vce(cluster id)
stcox c.stdlogp10909_rev##c.epcarbo ageadj i.sex, vce(cluster id)


////// analysis of trial data
use "F:\蛋白质组与死亡率\6. RCT概念性验证\trial_use.dta", clear

/// within group difference
mixed stdscore_imba i.time || code_id: if intervention == "HFLC"
mixed stdscore_acti i.time || code_id: if intervention == "HFLC"
mixed stdscore_regu i.time || code_id: if intervention == "HFLC"

mixed stdscore_imba i.time || code_id: if intervention == "LFHC"
mixed stdscore_acti i.time || code_id: if intervention == "LFHC"
mixed stdscore_regu i.time || code_id: if intervention == "LFHC"

///between group difference
mixed stdscore_imba i.time##i.group || code_id:
mixed stdscore_acti i.time##i.group || code_id:
mixed stdscore_regu i.time##i.group || code_id:

///within and between group difference for CRegS proteins
mixed stdlogp08603 i.time || code_id: if intervention == "HFLC"
mixed stdlogp08603 i.time || code_id: if intervention == "LFHC"

mixed stdlogp04003 i.time || code_id: if intervention == "HFLC"
mixed stdlogp04003 i.time || code_id: if intervention == "LFHC"

mixed stdlogp04004 i.time || code_id: if intervention == "HFLC"
mixed stdlogp04004 i.time || code_id: if intervention == "LFHC"

mixed stdlogp10909 i.time || code_id: if intervention == "HFLC"
mixed stdlogp10909 i.time || code_id: if intervention == "LFHC"

mixed stdlogp08603 i.time##i.group || code_id:
mixed stdlogp04003 i.time##i.group || code_id:
mixed stdlogp04004 i.time##i.group || code_id:
mixed stdlogp10909 i.time##i.group || code_id:






















