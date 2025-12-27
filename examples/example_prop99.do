/*==============================================================================
    California Proposition 99 Tobacco Control Analysis
==============================================================================*/

clear all
set more off

* Load data
findfile prop99_smoking.dta
use "`r(fn)'", clear

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
    Analysis 4: Overall Mode (Less Conservative)
    Use overall violations mode when iterative mode is too strict.
    This assumes cumulative errors are bounded, so kappa = 1.
==============================================================================*/

pretest cigsale, treatment(treated) time(year) treat_time(1989) threshold(5) overall

