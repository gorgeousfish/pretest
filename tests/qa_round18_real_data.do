/*******************************************************************************
* Quality Assurance Test - Round 18
* Testing pretest-stata with REAL datasets
* 
* Author: Master Test Architect (Murat)
* Date: 2025-12-21
* 
* Test Strategy:
*   1. Stata built-in nlswork data
*   2. Mathematical correctness verification
*   3. Boundary condition testing
*   4. Reasonableness checks
*******************************************************************************/

clear all
set more off
capture log close
log using "/Users/cxy/Desktop/2026project/pretest/pretest-stata/tests/qa_round18_real_data.log", replace text

* Set up paths (absolute paths as required)
local base_path "/Users/cxy/Desktop/2026project/pretest/pretest-stata"
adopath + "`base_path'"
adopath + "`base_path'/ado"

* Load Mata functions
capture do "`base_path'/_pretest_mata.do"

di as text ""
di as text "=============================================================="
di as text "  QUALITY ASSURANCE TEST - ROUND 18"
di as text "  Real Data Testing with Reasonableness Verification"
di as text "=============================================================="
di as text ""

* Initialize test counters
local n_tests = 0
local n_pass = 0
local n_fail = 0
local n_warning = 0

********************************************************************************
* TEST 1: Stata Built-in nlswork Data (REAL DATA)
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 1: Stata nlswork Data (Union Wage Premium) - REAL DATA"
di as text "=============================================================="

* Load nlswork data
webuse nlswork, clear

* Data description
di as text ""
di as text "nlswork data: National Longitudinal Survey of Young Women (1968-1988)"
di as text "Outcome: ln_wage (log hourly wage)"
di as text "Treatment: union membership (as quasi-treatment)"

* Keep complete cases for key variables
keep if !missing(ln_wage, union, year, idcode)

* Create continuous time index (year has gaps)
egen time = group(year)
xtset idcode time

* Summary
di as text ""
di as text "Data Summary:"
tabstat ln_wage, by(time) stat(n mean sd)
tab union time if time <= 5

* Test with treat_time = 6
local ++n_tests
di as text ""
di as text "Testing with treat_time = 6..."

capture noisily pretest ln_wage, treatment(union) time(time) threshold(0.5) treat_time(6) nograph

if _rc == 0 {
    local s_pre = e(S_pre)
    local kappa = e(kappa)
    local delta = e(delta_bar)
    local ci_lo = e(ci_lower)
    local ci_hi = e(ci_upper)
    local ci_width = `ci_hi' - `ci_lo'
    local pretest_pass = e(pretest_pass)
    
    di as text ""
    di as text "Results:"
    di as text "  S_pre = " %9.4f `s_pre'
    di as text "  kappa = " %9.4f `kappa'
    di as text "  delta_bar (ATT estimate) = " %9.4f `delta'
    di as text "  CI = [" %9.4f `ci_lo' ", " %9.4f `ci_hi' "]"
    di as text "  CI width = " %9.4f `ci_width'
    di as text "  Pretest pass = " `pretest_pass'
    
    * Reasonableness checks for log wage
    * Union wage premium typically 10-20% (0.1-0.2 in log terms)
    local reasonable = 1
    
    di as text ""
    di as text "Reasonableness Assessment:"
    
    if `s_pre' > 1 {
        di as error "  [?] S_pre = " %9.4f `s_pre' " seems high for wage data"
        local reasonable = 0
    }
    else {
        di as text "  [OK] S_pre magnitude is plausible"
    }
    
    if abs(`delta') > 1 {
        di as error "  [?] delta_bar > 1 seems unreasonable for log wage"
        local reasonable = 0
    }
    else {
        di as text "  [OK] delta_bar magnitude is plausible for log wage"
    }
    
    if `ci_width' > 5 {
        di as error "  [?] CI width > 5 seems too wide"
        local reasonable = 0
    }
    else {
        di as text "  [OK] CI width is reasonable"
    }
    
    if `reasonable' == 1 {
        di as result "  [PASS] All results appear reasonable"
        local ++n_pass
    }
    else {
        di as text "  [WARNING] Some results need investigation"
        local ++n_warning
    }
}
else {
    di as error "  [FAIL] Command failed with error " _rc
    local ++n_fail
}

********************************************************************************
* TEST 2: Mathematical Formula Verification
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 2: Mathematical Formula Verification"
di as text "=============================================================="

* Create synthetic data with KNOWN parameters
clear
set seed 54321
local n_units = 200
local T_total = 6
local t0 = 4

* Known violation pattern: nu_t = 0.05 for all t
local true_nu = 0.05

set obs `=`n_units' * `T_total''
gen id = ceil(_n / `T_total')
bysort id: gen t = _n
gen D = (id <= `n_units'/2)

* Generate Y with known structure
gen Y = 10 + 0.5 * t + rnormal(0, 0.5)

* Add cumulative violation for treated group
forvalues s = 2/`T_total' {
    replace Y = Y + `true_nu' if D == 1 & t >= `s'
}

* Add treatment effect for post-treatment
replace Y = Y + 2 if D == 1 & t >= `t0'

xtset id t

* Manual calculation of nu_t
di as text ""
di as text "Manual calculation of iterative violations nu_t:"
forvalues time = 2/`=`t0'-1' {
    local prev_time = `time' - 1
    
    * Mean for treated
    qui sum Y if D == 1 & t == `time'
    local Y_t_D1 = r(mean)
    qui sum Y if D == 1 & t == `prev_time'
    local Y_tm1_D1 = r(mean)
    
    * Mean for control
    qui sum Y if D == 0 & t == `time'
    local Y_t_D0 = r(mean)
    qui sum Y if D == 0 & t == `prev_time'
    local Y_tm1_D0 = r(mean)
    
    * nu_t = (Y_t - Y_{t-1})_D=1 - (Y_t - Y_{t-1})_D=0
    local nu_`time' = (`Y_t_D1' - `Y_tm1_D1') - (`Y_t_D0' - `Y_tm1_D0')
    di as text "  nu_`time' = " %9.4f `nu_`time'' " (expected: ~" %9.4f `true_nu' ")"
}

* Manual calculation of S_pre (p=2)
local sum_nu_sq = 0
local T_pre_minus_1 = `t0' - 2
forvalues time = 2/`=`t0'-1' {
    local sum_nu_sq = `sum_nu_sq' + (`nu_`time'')^2
}
local manual_S_pre = sqrt(`sum_nu_sq' / `T_pre_minus_1')
di as text ""
di as text "Manual S_pre (p=2) = " %9.4f `manual_S_pre'

* Run pretest
local ++n_tests
di as text ""
di as text "Running pretest..."
capture noisily pretest Y, treatment(D) time(t) threshold(0.5) treat_time(`t0') p(2) nograph

if _rc == 0 {
    local pkg_S_pre = e(S_pre)
    local pkg_kappa = e(kappa)
    local pkg_delta = e(delta_bar)
    
    di as text ""
    di as text "Comparison:"
    di as text "  Manual S_pre  = " %9.6f `manual_S_pre'
    di as text "  Package S_pre = " %9.6f `pkg_S_pre'
    local diff_S = abs(`manual_S_pre' - `pkg_S_pre')
    di as text "  Difference    = " %9.6f `diff_S'
    
    * Verify kappa formula
    local T_post = `T_total' - `t0' + 1
    local sum_sq = 0
    forvalues s = 1/`T_post' {
        local sum_sq = `sum_sq' + `s'^2
    }
    local manual_kappa = sqrt(`sum_sq' / `T_post')
    di as text ""
    di as text "Kappa verification (p=2, T_post=`T_post'):"
    di as text "  Manual kappa  = " %9.6f `manual_kappa'
    di as text "  Package kappa = " %9.6f `pkg_kappa'
    local diff_kappa = abs(`manual_kappa' - `pkg_kappa')
    di as text "  Difference    = " %9.6e `diff_kappa'
    
    if `diff_S' < 0.01 & `diff_kappa' < 1e-10 {
        di as result "  [PASS] Mathematical formulas verified"
        local ++n_pass
    }
    else {
        di as error "  [FAIL] Formula discrepancy detected"
        local ++n_fail
    }
}
else {
    di as error "  [FAIL] Command failed"
    local ++n_fail
}

********************************************************************************
* TEST 3: T_post = 1 Boundary (BUG-001 verification)
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 3: Boundary Condition - T_post = 1 (BUG-001 verification)"
di as text "=============================================================="

clear
set seed 99999
set obs 300
gen id = ceil(_n/3)
gen time = mod(_n-1, 3) + 1
gen treat = (id <= 50)
gen y = 10 + 0.5*time + rnormal(0, 1)
replace y = y + 2 if treat == 1 & time >= 3
xtset id time

local ++n_tests
di as text "Testing with T=3, t0=3 (T_post=1)..."

capture noisily pretest y, treatment(treat) time(time) threshold(0.5) treat_time(3) nograph

if _rc == 0 {
    local se_val = e(se_delta_bar)
    local ci_conv_lo = e(ci_conv_lower)
    local ci_conv_hi = e(ci_conv_upper)
    local ci_conv_width = `ci_conv_hi' - `ci_conv_lo'
    
    di as text "  se_delta_bar = " %12.6e `se_val'
    di as text "  Conventional CI width = " %12.6e `ci_conv_width'
    
    if `se_val' < 0.01 {
        di as error "  [FAIL] BUG-001 CONFIRMED: se_delta_bar abnormally small"
        local ++n_fail
    }
    else {
        di as result "  [PASS] se_delta_bar appears reasonable"
        local ++n_pass
    }
}
else {
    di as error "  [FAIL] Command failed"
    local ++n_fail
}

********************************************************************************
* TEST 4: Cluster Standard Errors (BUG-008 verification)
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 4: Cluster Standard Errors (BUG-008 verification)"
di as text "=============================================================="

clear
set seed 88888
set obs 2500
gen id = ceil(_n / 5)
bysort id: gen t = _n
gen D = (id <= 250)
gen cluster_id = ceil(id / 10)
gen Y = 10 + 0.5*t + 2*D + rnormal(0, 1)
xtset id t

local ++n_tests

* Test without clustering
di as text ""
di as text "Without clustering:"
capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) nograph
if _rc == 0 {
    local ci_width_no_cluster = e(ci_upper) - e(ci_lower)
    local f_alpha_no_cluster = e(f_alpha)
    local se_no_cluster = e(se_delta_bar)
    di as text "  CI width = " %9.4f `ci_width_no_cluster'
    di as text "  f_alpha = " %9.4f `f_alpha_no_cluster'
    di as text "  se_delta_bar = " %9.4f `se_no_cluster'
}

* Test with clustering
di as text ""
di as text "With clustering:"
capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) cluster(cluster_id) nograph
if _rc == 0 {
    local ci_width_cluster = e(ci_upper) - e(ci_lower)
    local f_alpha_cluster = e(f_alpha)
    local se_cluster = e(se_delta_bar)
    di as text "  CI width = " %9.4f `ci_width_cluster'
    di as text "  f_alpha = " %9.4f `f_alpha_cluster'
    di as text "  se_delta_bar = " %9.4f `se_cluster'
    
    local ci_ratio = `ci_width_cluster' / `ci_width_no_cluster'
    local f_ratio = `f_alpha_cluster' / `f_alpha_no_cluster'
    local se_ratio = `se_cluster' / `se_no_cluster'
    
    di as text ""
    di as text "Ratios (cluster/no-cluster):"
    di as text "  CI width ratio = " %9.2f `ci_ratio' "x"
    di as text "  f_alpha ratio = " %9.2f `f_ratio' "x"
    di as text "  se ratio = " %9.2f `se_ratio' "x"
    
    if `ci_ratio' > 10 {
        di as error "  [FAIL] BUG-008 CONFIRMED: Cluster SE ratio = " %9.1f `ci_ratio' "x (abnormal)"
        local ++n_fail
    }
    else {
        di as result "  [PASS] Cluster SE ratio appears reasonable"
        local ++n_pass
    }
}
else {
    di as error "  [FAIL] Cluster command failed"
    local ++n_fail
}

********************************************************************************
* TEST 5: Large Violation Rejection
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 5: Large Violation Should Be Rejected"
di as text "=============================================================="

clear
set seed 77777
set obs 500
gen id = ceil(_n / 5)
bysort id: gen t = _n
gen D = (id <= 50)

gen Y = 10 + rnormal(0, 0.5)
forvalues s = 2/5 {
    replace Y = Y + 0.8 if D == 1 & t >= `s'
}
replace Y = Y + 2 if D == 1 & t >= 4

xtset id t

local ++n_tests
di as text "Testing with large pre-treatment violation (S_pre should be ~0.8)..."
di as text "Using threshold M = 0.3 (should reject)..."

capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.3) nograph

if _rc == 0 {
    local S_pre = e(S_pre)
    local phi = e(phi)
    local pretest_pass = e(pretest_pass)
    
    di as text "  S_pre = " %9.4f `S_pre'
    di as text "  phi = " `phi'
    di as text "  pretest_pass = " `pretest_pass'
    
    if `pretest_pass' == 0 & `phi' == 1 {
        di as result "  [PASS] Correctly rejected (S_pre > M)"
        local ++n_pass
    }
    else {
        di as error "  [FAIL] Should have rejected but passed"
        local ++n_fail
    }
}
else {
    di as error "  [FAIL] Command failed"
    local ++n_fail
}

********************************************************************************
* TEST 6: Perfect Parallel Trends
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 6: Perfect Parallel Trends (S_pre should be ~0)"
di as text "=============================================================="

clear
set seed 66666
set obs 1000
gen id = ceil(_n / 5)
bysort id: gen t = _n
gen D = (id <= 100)

gen Y = 10 + 0.5 * t + 2 * D + rnormal(0, 0.3)
replace Y = Y + 3 if D == 1 & t >= 4

xtset id t

local ++n_tests
di as text "Testing with perfect parallel trends..."

capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) nograph

if _rc == 0 {
    local S_pre = e(S_pre)
    local pretest_pass = e(pretest_pass)
    local delta = e(delta_bar)
    
    di as text "  S_pre = " %9.4f `S_pre' " (should be ~0 due to sampling noise)"
    di as text "  pretest_pass = " `pretest_pass'
    di as text "  delta_bar = " %9.4f `delta' " (true effect = 3)"
    
    if `S_pre' < 0.15 & `pretest_pass' == 1 {
        di as result "  [PASS] S_pre is small and test passed"
        local ++n_pass
    }
    else if `S_pre' >= 0.15 {
        di as error "  [FAIL] S_pre unexpectedly large for parallel trends data"
        local ++n_fail
    }
    else {
        di as error "  [FAIL] Unexpected pretest result"
        local ++n_fail
    }
    
    if abs(`delta' - 3) < 0.5 {
        di as result "  [PASS] delta_bar is close to true effect (3)"
    }
    else {
        di as text "  [WARNING] delta_bar differs from true effect"
    }
}
else {
    di as error "  [FAIL] Command failed"
    local ++n_fail
}

********************************************************************************
* TEST 7: Overall Mode kappa = 1
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 7: Overall Mode - kappa should equal 1"
di as text "=============================================================="

local ++n_tests

di as text ""
di as text "Overall mode:"
capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) overall nograph
if _rc == 0 {
    local kappa_overall = e(kappa)
    di as text "  kappa = " %9.6f `kappa_overall'
    
    if abs(`kappa_overall' - 1) < 0.001 {
        di as result "  [PASS] Overall mode kappa = 1 as expected"
        local ++n_pass
    }
    else {
        di as error "  [FAIL] Overall mode kappa should be 1, got " %9.4f `kappa_overall'
        local ++n_fail
    }
}
else {
    di as error "  [FAIL] Overall mode command failed"
    local ++n_fail
}

********************************************************************************
* TEST 8: Different p values
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST 8: Different p values (monotonicity check)"
di as text "=============================================================="

local ++n_tests

* Create data with known violations
clear
set seed 55555
set obs 600
gen id = ceil(_n / 6)
bysort id: gen t = _n
gen D = (id <= 50)
gen Y = 10 + 0.3*t + rnormal(0, 0.3)
forvalues s = 2/6 {
    replace Y = Y + 0.1 if D == 1 & t >= `s'
}
replace Y = Y + 2 if D == 1 & t >= 4
xtset id t

di as text "Testing S_pre for different p values..."
di as text "(S_pre should increase as p increases)"

local S_p1 = .
local S_p2 = .
local S_p10 = .

capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) p(1) nograph
if _rc == 0 {
    local S_p1 = e(S_pre)
    di as text "  p=1:  S_pre = " %9.4f `S_p1'
}

capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) p(2) nograph
if _rc == 0 {
    local S_p2 = e(S_pre)
    di as text "  p=2:  S_pre = " %9.4f `S_p2'
}

capture noisily pretest Y, treatment(D) time(t) treat_time(4) threshold(0.5) p(10) nograph
if _rc == 0 {
    local S_p10 = e(S_pre)
    di as text "  p=10: S_pre = " %9.4f `S_p10'
}

* Check monotonicity: S_pre(p=1) <= S_pre(p=2) <= S_pre(p=10)
if `S_p1' <= `S_p2' + 0.001 & `S_p2' <= `S_p10' + 0.001 {
    di as result "  [PASS] S_pre is monotonic in p"
    local ++n_pass
}
else {
    di as error "  [FAIL] S_pre is not monotonic in p"
    local ++n_fail
}

********************************************************************************
* Summary
********************************************************************************
di as text ""
di as text "=============================================================="
di as text "TEST SUMMARY - ROUND 18"
di as text "=============================================================="
di as text ""
di as text "Total tests: " `n_tests'
di as result "  Passed: " `n_pass'
if `n_fail' > 0 {
    di as error "  Failed: " `n_fail'
}
else {
    di as text "  Failed: " `n_fail'
}
di as text "  Warnings: " `n_warning'
di as text ""
local pass_rate = round(`n_pass' / `n_tests' * 100, 0.1)
di as text "Pass rate: " %5.1f `pass_rate' "%"

if `n_fail' > 0 {
    di as text ""
    di as error "CONFIRMED BUGS:"
    di as error "  - BUG-001: T_post=1 se_delta_bar issue (if failed)"
    di as error "  - BUG-008: Cluster SE abnormally large (if failed)"
}

di as text ""
di as text "=============================================================="
di as text "End of Round 18 Quality Assurance Test"
di as text "=============================================================="

log close
