*! _pretest_validate v0.1.0 - Data Validation Module
*!
*! Description:
*!   Validates input data for the pretest command.
*!   - Checks variable existence (error 101)
*!   - Checks variable type is numeric (error 102)
*!   - Checks treatment variable is binary 0/1 (error 103)
*!   - Detects and warns about missing values
*!
*! Syntax:
*!   _pretest_validate outcome, treatment(varname) time(varname)
*!
*! Error codes:
*!   101 - Variable not found
*!   102 - Variable type error (non-numeric)
*!   103 - Treatment variable not binary 0/1
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), Section 2.1
*!   Treatment indicator D ∈ {0, 1} binary requirement

program define _pretest_validate, rclass
    version 17.0
    
    // Parse syntax - use string to accept variable names for custom error handling
    // Not using varlist/varname to avoid Stata's native errors at parse stage
    syntax anything(name=outcome) , TREATment(string) TIME(string)
    
    local treatment `treatment'
    local time `time'
    
    // ========================================
    // 1. Check outcome variable existence and type
    // ========================================
    
    // 1.1 Existence check
    capture confirm variable `outcome'
    if _rc {
        di as error "Error 101: Outcome variable '`outcome'' not found"
        exit 101
    }
    
    // 1.2 Numeric type check
    capture confirm numeric variable `outcome'
    if _rc {
        di as error "Error 102: Outcome variable '`outcome'' must be numeric"
        exit 102
    }
    
    // ========================================
    // 2. Check treatment variable existence, type, and binary
    // ========================================
    
    // 2.1 Existence check
    capture confirm variable `treatment'
    if _rc {
        di as error "Error 101: Treatment variable '`treatment'' not found"
        exit 101
    }
    
    // 2.2 Numeric type check
    capture confirm numeric variable `treatment'
    if _rc {
        di as error "Error 102: Treatment variable '`treatment'' must be numeric"
        exit 102
    }
    
    // 2.3 Binary check (only 0, 1, and missing allowed)
    // Paper requirement: D ∈ {0, 1} for control and treatment groups
    qui count if `treatment' != 0 & `treatment' != 1 & !missing(`treatment')
    if r(N) > 0 {
        di as error "Error 103: Treatment variable '`treatment'' must contain only 0 and 1"
        di as error "  Found " r(N) " observations with values other than 0, 1, or missing"
        
        // Display summary of invalid values
        qui tab `treatment' if `treatment' != 0 & `treatment' != 1 & !missing(`treatment'), matrow(badvals)
        di as error "  Invalid values found:"
        matrix list badvals, noheader
        
        exit 103
    }
    
    // ========================================
    // 3. Check time variable existence and type
    // ========================================
    
    // 3.1 Existence check
    capture confirm variable `time'
    if _rc {
        di as error "Error 101: Time variable '`time'' not found"
        exit 101
    }
    
    // 3.2 Numeric type check
    capture confirm numeric variable `time'
    if _rc {
        di as error "Error 102: Time variable '`time'' must be numeric"
        exit 102
    }
    
    // ========================================
    // 4. Missing value detection and warnings
    // ========================================
    
    local has_warning = 0
    
    // 4.1 Outcome variable missing values
    qui count if missing(`outcome')
    if r(N) > 0 {
        local n_miss_y = r(N)
        qui count
        local n_total = r(N)
        local pct_miss = 100 * `n_miss_y' / `n_total'
        di as text "Warning: Outcome variable '`outcome'' has " ///
            as result `n_miss_y' as text " missing values (" ///
            as result %5.2f `pct_miss' as text "%)"
        local has_warning = 1
    }
    
    // 4.2 Treatment variable missing values
    qui count if missing(`treatment')
    if r(N) > 0 {
        local n_miss_d = r(N)
        qui count
        local n_total = r(N)
        local pct_miss = 100 * `n_miss_d' / `n_total'
        di as text "Warning: Treatment variable '`treatment'' has " ///
            as result `n_miss_d' as text " missing values (" ///
            as result %5.2f `pct_miss' as text "%)"
        local has_warning = 1
    }
    
    // 4.3 Time variable missing values
    qui count if missing(`time')
    if r(N) > 0 {
        local n_miss_t = r(N)
        qui count
        local n_total = r(N)
        local pct_miss = 100 * `n_miss_t' / `n_total'
        di as text "Warning: Time variable '`time'' has " ///
            as result `n_miss_t' as text " missing values (" ///
            as result %5.2f `pct_miss' as text "%)"
        local has_warning = 1
    }
    
    // ========================================
    // 5. Return validation results
    // ========================================
    
    // Compute basic statistics
    qui count if !missing(`outcome') & !missing(`treatment') & !missing(`time')
    local n_valid = r(N)
    
    qui count
    local n_total = r(N)
    
    // Count treated and control groups
    qui count if `treatment' == 1 & !missing(`outcome') & !missing(`time')
    local n_treated = r(N)
    
    qui count if `treatment' == 0 & !missing(`outcome') & !missing(`time')
    local n_control = r(N)
    
    // Display validation summary
    di as text _n "Data Validation Summary:"
    di as text "  Total observations:       " as result %10.0fc `n_total'
    di as text "  Valid observations:       " as result %10.0fc `n_valid'
    di as text "  Treated observations:     " as result %10.0fc `n_treated'
    di as text "  Control observations:     " as result %10.0fc `n_control'
    
    if `has_warning' {
        di as text _n "  Note: Missing values will be excluded from analysis"
    }
    
    // Return results
    return scalar n_total = `n_total'
    return scalar n_valid = `n_valid'
    return scalar n_treated = `n_treated'
    return scalar n_control = `n_control'
    return scalar has_warning = `has_warning'
    
    return local outcome "`outcome'"
    return local treatment "`treatment'"
    return local time "`time'"
end
