*! _pretest_validate v0.1.0 - Data Validation Module for pretest
*!
*! Purpose:
*!   Validates input data structure and variable requirements for
*!   difference-in-differences analysis under the conditional extrapolation
*!   framework.
*!
*! Syntax:
*!   _pretest_validate outcome, treatment(varname) time(varname)
*!
*! Validation Rules:
*!   - All variables must exist in the dataset
*!   - All variables must be numeric (required for DID estimation)
*!   - Treatment indicator D must be binary {0, 1} per Section 2.1
*!   - Missing values are detected and reported
*!
*! Error Codes:
*!   101 - Variable not found in dataset
*!   102 - Variable is not numeric type
*!   103 - Treatment variable contains values other than 0 and 1
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test:
*!   Valid Inference when Testing Violations of Parallel Trends for
*!   Difference-in-Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470
*!   See Section 2.1 for treatment indicator specification.

program define _pretest_validate, rclass
    version 17.0
    
    // Parse syntax with string types to provide custom error messages
    syntax anything(name=outcome) , TREATment(string) TIME(string)
    
    local treatment `treatment'
    local time `time'
    
    // =========================================================================
    // Section 1: Outcome Variable Validation
    // =========================================================================
    
    // Check existence
    capture confirm variable `outcome'
    if _rc {
        di as error "Error 101: Outcome variable '`outcome'' not found"
        exit 101
    }
    
    // Check numeric type (required for computing group means)
    capture confirm numeric variable `outcome'
    if _rc {
        di as error "Error 102: Outcome variable '`outcome'' must be numeric"
        exit 102
    }
    
    // =========================================================================
    // Section 2: Treatment Variable Validation
    // Per Section 2.1 of Mikhaeil & Harshaw (2025), D must be binary {0, 1}
    // where D=0 denotes control group and D=1 denotes treatment group.
    // =========================================================================
    
    // Check existence
    capture confirm variable `treatment'
    if _rc {
        di as error "Error 101: Treatment variable '`treatment'' not found"
        exit 101
    }
    
    // Check numeric type
    capture confirm numeric variable `treatment'
    if _rc {
        di as error "Error 102: Treatment variable '`treatment'' must be numeric"
        exit 102
    }
    
    // Verify binary structure: only {0, 1} values permitted (missing allowed)
    qui count if `treatment' != 0 & `treatment' != 1 & !missing(`treatment')
    if r(N) > 0 {
        di as error "Error 103: Treatment variable '`treatment'' must contain only 0 and 1"
        di as error "  Found " r(N) " observations with invalid values"
        
        // Display invalid values for diagnostics
        qui tab `treatment' if `treatment' != 0 & `treatment' != 1 & !missing(`treatment'), matrow(badvals)
        di as error "  Invalid values:"
        matrix list badvals, noheader
        
        exit 103
    }
    
    // =========================================================================
    // Section 3: Time Variable Validation
    // Time variable indexes periods t = 1, ..., T with treatment at t0.
    // =========================================================================
    
    // Check existence
    capture confirm variable `time'
    if _rc {
        di as error "Error 101: Time variable '`time'' not found"
        exit 101
    }
    
    // Check numeric type (required for period ordering)
    capture confirm numeric variable `time'
    if _rc {
        di as error "Error 102: Time variable '`time'' must be numeric"
        exit 102
    }
    
    // =========================================================================
    // Section 4: Missing Value Detection
    // Missing observations are excluded from DID estimation. Users are warned
    // when missing values are present to ensure transparency in sample size.
    // =========================================================================
    
    local has_warning = 0
    qui count
    local n_total = r(N)
    
    // Outcome variable
    qui count if missing(`outcome')
    if r(N) > 0 {
        local n_miss_y = r(N)
        local pct_miss = 100 * `n_miss_y' / `n_total'
        di as text "Warning: Outcome variable '`outcome'' has " ///
            as result `n_miss_y' as text " missing values (" ///
            as result %5.2f `pct_miss' as text "%)"
        local has_warning = 1
    }
    
    // Treatment variable
    qui count if missing(`treatment')
    if r(N) > 0 {
        local n_miss_d = r(N)
        local pct_miss = 100 * `n_miss_d' / `n_total'
        di as text "Warning: Treatment variable '`treatment'' has " ///
            as result `n_miss_d' as text " missing values (" ///
            as result %5.2f `pct_miss' as text "%)"
        local has_warning = 1
    }
    
    // Time variable
    qui count if missing(`time')
    if r(N) > 0 {
        local n_miss_t = r(N)
        local pct_miss = 100 * `n_miss_t' / `n_total'
        di as text "Warning: Time variable '`time'' has " ///
            as result `n_miss_t' as text " missing values (" ///
            as result %5.2f `pct_miss' as text "%)"
        local has_warning = 1
    }
    
    // =========================================================================
    // Section 5: Summary Statistics and Return Values
    // Report sample composition for treated (D=1) and control (D=0) groups.
    // =========================================================================
    
    // Count valid observations (complete cases)
    qui count if !missing(`outcome') & !missing(`treatment') & !missing(`time')
    local n_valid = r(N)
    
    // Count by treatment status
    qui count if `treatment' == 1 & !missing(`outcome') & !missing(`time')
    local n_treated = r(N)
    
    qui count if `treatment' == 0 & !missing(`outcome') & !missing(`time')
    local n_control = r(N)
    
    // Display validation summary
    di as text _n "Data Validation Summary:"
    di as text "  Total observations:       " as result %10.0fc `n_total'
    di as text "  Valid observations:       " as result %10.0fc `n_valid'
    di as text "  Treated (D=1):            " as result %10.0fc `n_treated'
    di as text "  Control (D=0):            " as result %10.0fc `n_control'
    
    if `has_warning' {
        di as text _n "  Note: Missing values will be excluded from analysis"
    }
    
    // Return scalar results
    return scalar n_total = `n_total'
    return scalar n_valid = `n_valid'
    return scalar n_treated = `n_treated'
    return scalar n_control = `n_control'
    return scalar has_warning = `has_warning'
    
    // Return variable names
    return local outcome "`outcome'"
    return local treatment "`treatment'"
    return local time "`time'"
end
