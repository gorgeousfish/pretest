/*==============================================================================
    California Proposition 99 Tobacco Control Analysis
==============================================================================*/

clear all
set more off

* Load data: prefer a local shipped/downloaded copy and fall back to GitHub only
* when no local file is available.
capture noisily findfile prop99_smoking.dta
if _rc == 0 {
    use "`r(fn)'", clear
}
else {
    webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"
    webuse prop99_smoking.dta, clear
}

* Set panel structure
xtset state year

/*==============================================================================
    Analysis 1: Full sample (1970-2000), M = 5
==============================================================================*/


pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)

/*==============================================================================
    Analysis 2: Full sample, smaller threshold M = 1 (stricter)
==============================================================================*/


pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(1) nograph


/*==============================================================================
    Analysis 3: Shorter window (1985-1995), M = 5
==============================================================================*/



keep if year >= 1985 & year <= 1995
pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5)

/*==============================================================================
    Analysis 4: Overall Mode (Cumulative-Violation Path) -- full sample
    Reload full sample before running overall so this is not a continuation of
    the 1985-1995 windowed sample from Analysis 3.
    Use overall when cumulative violations are the empirical question.
    This replaces the iterative violation path and sets kappa = 1.
==============================================================================*/

* Reload full dataset to restore the full 1970-2000 sample
capture noisily findfile prop99_smoking.dta
if _rc == 0 {
    use "`r(fn)'", clear
}
else {
    webuse set "https://raw.githubusercontent.com/gorgeousfish/pretest/main/"
    webuse prop99_smoking.dta, clear
}
xtset state year

pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5) overall
