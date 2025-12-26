*! _pretest_check_coverage.ado
*! Version 0.1.0
*!
*! Purpose:
*!   Validates that all (time, treatment) cell combinations contain observations.
*!   This is a prerequisite for the pretest method, which requires the complete
*!   parameter vector theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)
*!   to compute pre-treatment severity S_pre (Mikhaeil & Harshaw, 2025, Section 3.1).
*!
*! Syntax:
*!   _pretest_check_coverage depvar, treatment(varname) time(varname) [options]
*!
*! Required:
*!   depvar           - Outcome variable
*!   treatment(var)   - Binary treatment indicator (0 = control, 1 = treated)
*!   time(var)        - Time period variable
*!
*! Options:
*!   strict           - Return error if any cell is missing (default: warn only)
*!   silent           - Suppress diagnostic output
*!   minobs(integer)  - Minimum observations per cell (default: 1)
*!
*! Returned values:
*!   r(coverage_ok)      - 1 if all cells have >= minobs observations, 0 otherwise
*!   r(n_missing_cells)  - Count of cells with zero observations
*!   r(n_low_cells)      - Count of cells with 0 < n < minobs observations
*!   r(n_valid_obs)      - Total non-missing observations
*!   r(n_time_periods)   - Number of distinct time periods
*!   r(missing_cells)    - String listing cells with zero observations
*!   r(low_cells)        - String listing cells with low observations
*!   r(coverage)         - Matrix: (time, n_treated, n_control)
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

program define _pretest_check_coverage, rclass
    version 17.0
    
    syntax varlist(min=1 max=1 numeric) , ///
        TREATment(varname numeric) ///
        TIME(varname numeric) ///
        [STRICT]  ///
        [SILENT]  ///
        [MINOBS(integer 1)]
    
    local outcome `varlist'
    
    // Extract unique time values from non-missing observations
    qui levelsof `time' if !missing(`outcome') & !missing(`treatment'), local(time_vals)
    local n_times : word count `time_vals'
    
    // Count total valid (non-missing) observations
    qui count if !missing(`outcome') & !missing(`treatment') & !missing(`time')
    local n_valid_obs = r(N)
    
    // Initialize tracking variables
    local n_missing_cells = 0
    local n_low_cells = 0
    local missing_cells ""
    local low_cells ""
    local coverage_ok = 1
    
    // Construct coverage matrix: rows = time periods, cols = (time, n_treated, n_control)
    tempname coverage_mat
    matrix `coverage_mat' = J(`n_times', 3, .)
    matrix colnames `coverage_mat' = time n_treated n_control
    local row = 1
    
    // Iterate over each time period to check cell coverage
    foreach t of local time_vals {
        // Count observations in treated group (D=1) at time t
        qui count if `time' == `t' & `treatment' == 1 & !missing(`outcome')
        local n_t1 = r(N)
        
        // Count observations in control group (D=0) at time t
        qui count if `time' == `t' & `treatment' == 0 & !missing(`outcome')
        local n_t0 = r(N)
        
        // Record counts in coverage matrix
        matrix `coverage_mat'[`row', 1] = `t'
        matrix `coverage_mat'[`row', 2] = `n_t1'
        matrix `coverage_mat'[`row', 3] = `n_t0'
        
        // Flag cells with insufficient observations
        if `n_t1' < `minobs' {
            if `n_t1' == 0 {
                local n_missing_cells = `n_missing_cells' + 1
                local missing_cells "`missing_cells' (`time'=`t', D=1)"
            }
            else {
                local n_low_cells = `n_low_cells' + 1
                local low_cells "`low_cells' (`time'=`t', D=1, n=`n_t1')"
            }
            local coverage_ok = 0
        }
        
        if `n_t0' < `minobs' {
            if `n_t0' == 0 {
                local n_missing_cells = `n_missing_cells' + 1
                local missing_cells "`missing_cells' (`time'=`t', D=0)"
            }
            else {
                local n_low_cells = `n_low_cells' + 1
                local low_cells "`low_cells' (`time'=`t', D=0, n=`n_t0')"
            }
            local coverage_ok = 0
        }
        
        local ++row
    }
    
    // Display diagnostic output unless silent mode is enabled
    if "`silent'" == "" {
        if `coverage_ok' == 0 {
            di as text _n "Time-Treatment Coverage Diagnostic"
            di as text "{hline 55}"
            di as text "  `time'" _col(18) "Treated" _col(33) "Control" _col(48) "Status"
            di as text "{hline 55}"
            
            local row = 1
            foreach t of local time_vals {
                local n_t1 = `coverage_mat'[`row', 2]
                local n_t0 = `coverage_mat'[`row', 3]
                
                // Classify cell status based on observation counts
                if `n_t1' == 0 | `n_t0' == 0 {
                    local status "MISSING"
                    di as text "  `t'" _col(18) as result %6.0f `n_t1' ///
                        _col(33) as result %6.0f `n_t0' ///
                        _col(48) as error "`status'"
                }
                else if `n_t1' < 10 | `n_t0' < 10 {
                    local status "LOW"
                    di as text "  `t'" _col(18) as result %6.0f `n_t1' ///
                        _col(33) as result %6.0f `n_t0' ///
                        _col(48) as text "`status'"
                }
                else {
                    local status "OK"
                    di as text "  `t'" _col(18) as result %6.0f `n_t1' ///
                        _col(33) as result %6.0f `n_t0' ///
                        _col(48) as result "`status'"
                }
                
                local ++row
            }
            
            di as text "{hline 55}"
            
            // Display warnings for problematic cells
            if `n_missing_cells' > 0 {
                di as error "Warning: `n_missing_cells' cell(s) have zero observations"
                di as text "  Missing cells:`missing_cells'"
                di as text "  Estimates cannot be computed for these time periods."
            }
            if `n_low_cells' > 0 {
                di as text "Note: `n_low_cells' cell(s) have fewer than `minobs' observations"
            }
            
            di as text _n "Recommendations:"
            di as text "  1. Restrict sample to time periods with complete data"
            di as text "  2. Verify the treatment time specification"
            di as text "  3. Check for data quality issues"
            
            // In strict mode, terminate with error if coverage is incomplete
            // Rationale: The pretest method requires the complete parameter vector
            //   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)
            // to compute S_pre = ((1/(T_pre-1)) * sum|nu_t|^p)^{1/p}.
            // Missing time periods make S_pre undefined (Mikhaeil & Harshaw, 2025, Section 3.1).
            if "`strict'" != "" {
                di as error _n "Error: Incomplete (time, treatment) coverage"
                di as error "The pretest method requires observations at ALL time periods"
                di as error "for both treatment and control groups. The severity measure"
                di as error "S_pre cannot be computed with missing time periods."
                di as text _n "Possible solutions:"
                di as text "  - Restrict analysis to time periods with complete data"
                di as text "  - Consider alternative methods for incomplete panels"
                exit 120
            }
        }
        else {
            di as text "Coverage check passed: all cells have sufficient observations"
        }
    }
    
    // Return results
    return scalar coverage_ok = `coverage_ok'
    return scalar n_missing_cells = `n_missing_cells'
    return scalar n_low_cells = `n_low_cells'
    return scalar n_valid_obs = `n_valid_obs'
    return scalar n_time_periods = `n_times'
    return local missing_cells "`missing_cells'"
    return local low_cells "`low_cells'"
    return matrix coverage = `coverage_mat'
end
