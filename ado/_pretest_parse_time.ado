*! _pretest_parse_time.ado
*! version 0.1.0
*!
*! Time structure parsing for the pretest package.
*!
*! This module parses the time variable and computes key temporal parameters
*! required for the pretest methodology:
*!
*!   T      - Total number of time periods (unique values in timevar)
*!   t0     - Treatment time: first post-treatment period
*!   T_pre  - Number of pre-treatment periods: t0 - 1
*!   T_post - Number of post-treatment periods: T - T_pre
*!
*! The time index is normalized to consecutive integers 1, 2, ..., T internally,
*! regardless of the original time values in the data.
*!
*! Syntax:
*!   _pretest_parse_time timevar [, treatment(varname) treat_time(#) is_panel(#)]
*!
*! Options:
*!   treatment()  - Treatment indicator variable (for inferring t0)
*!   treat_time() - User-specified treatment time (takes precedence)
*!   is_panel()   - Panel data indicator (0 = repeated cross-sections, 1 = panel)
*!
*! Exit codes:
*!   100 - Treatment time not specified and cannot be inferred
*!   104 - Insufficient pre-treatment periods (T_pre < 2)
*!   109 - Specified treatment time not found in data
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470
*!
*!   Section 2.1: The analyst works with T time periods indexed t = 1, ..., T.
*!   Treatment is administered at t0, with T_pre = t0 - 1 pre-treatment periods.
*!   The methodology requires T_pre >= 2 for the iterative violations to be
*!   well-defined.

program define _pretest_parse_time, rclass
    version 17.0
    
    syntax varlist(min=1 max=1 numeric) ///
        [, TREATment(varname numeric) ///
           TREat_time(real -999) ///
           IS_panel(integer 0)]
    
    local timevar `varlist'
    
    // -------------------------------------------------------------------------
    // Step 1: Extract unique time values and compute T
    // -------------------------------------------------------------------------
    // Per Section 2.1, we work with T time periods. The original time values
    // are mapped to consecutive integers 1, ..., T for internal computations.
    
    qui levelsof `timevar', local(time_levels)
    local T : word count `time_levels'
    
    // Store original time values in a matrix for mapping back to user units
    tempname time_vals
    matrix `time_vals' = J(`T', 1, .)
    local i = 1
    foreach t of local time_levels {
        matrix `time_vals'[`i', 1] = `t'
        local ++i
    }
    
    // Extract time range bounds
    local t_min : word 1 of `time_levels'
    local t_max : word `T' of `time_levels'
    
    // -------------------------------------------------------------------------
    // Step 2: Detect non-consecutive time values
    // -------------------------------------------------------------------------
    // The methodology uses consecutive period indices internally. If the
    // original time variable has gaps (e.g., 2006, 2008, 2010), we note this
    // and proceed with the mapping to 1, 2, 3, etc.
    
    local has_gap = 0
    if `t_max' - `t_min' + 1 != `T' {
        local has_gap = 1
        di as text "Note: Time variable has non-consecutive values" ///
            " (`T' periods spanning `t_min' to `t_max')."
        di as text "      Internal consecutive indexing applied."
    }
    
    // -------------------------------------------------------------------------
    // Step 3: Determine treatment time t0
    // -------------------------------------------------------------------------
    // t0 is the first post-treatment period. Priority order:
    //   1. User-specified via treat_time() option
    //   2. Inferred from time-varying treatment indicator (panel data only)
    //   3. Error if neither available
    
    local t0 = .
    local t0_method = ""
    
    // --- Method 1: User-specified treatment time ---
    if `treat_time' != -999 {
        local t0_orig = `treat_time'
        
        // Validate that the specified time exists in the data
        local t0_found = 0
        local t0_index = 0
        local i = 1
        foreach t of local time_levels {
            if `t' == `t0_orig' {
                local t0_found = 1
                local t0_index = `i'
            }
            local ++i
        }
        
        // Handle invalid treatment time
        if `t0_found' == 0 {
            if `t0_orig' < `t_min' | `t0_orig' > `t_max' {
                di as error "treat_time(`t0_orig') is outside the data" ///
                    " time range [`t_min', `t_max']"
                exit 109
            }
            else {
                // Time is within range but not observed (gap)
                di as error "treat_time(`t0_orig') is not an observed" ///
                    " time value in the data"
                di as error "  Observed time values: `time_levels'"
                exit 109
            }
        }
        
        // Map to consecutive index
        local t0 = `t0_index'
        local t0_method "user-specified"
    }
    // --- Method 2: Infer from time-varying treatment in panel data ---
    else if `is_panel' == 1 & "`treatment'" != "" {
        // For panel data, attempt to infer t0 from treatment timing
        capture xtset
        if _rc == 0 {
            local panelvar = r(panelvar)
            
            // Check if treatment varies within units over time
            tempvar treat_sd
            qui bysort `panelvar': egen `treat_sd' = sd(`treatment')
            
            qui sum `treat_sd'
            if r(max) > 0 {
                // Treatment is time-varying; find first treated period
                qui levelsof `timevar' if `treatment' == 1, local(treat_times)
                
                if "`treat_times'" != "" {
                    local t0_orig : word 1 of `treat_times'
                    
                    // Map to consecutive index
                    local i = 1
                    foreach t of local time_levels {
                        if `t' == `t0_orig' {
                            local t0 = `i'
                        }
                        local ++i
                    }
                    local t0_method "inferred from time-varying treatment"
                    di as text "Note: Treatment time t0=`t0_orig'" ///
                        " (index `t0') inferred from treatment indicator."
                }
                else {
                    di as error "No observations with treatment=1 found."
                    di as error "  Please specify treat_time(#)."
                    exit 100
                }
            }
            else {
                // Treatment is time-invariant: block-adoption design
                // This is the standard case per Section 1.2 of the paper
                di as text ""
                di as text "{bf:Block-adoption design detected}"
                di as text "  Treatment is time-invariant within units, which is"
                di as text "  typical for block-adoption DID designs where all"
                di as text "  treated units receive treatment at the same time."
                di as text ""
                di as error "  Please specify {bf:treat_time(#)} to indicate t0,"
                di as text "  the first post-treatment period."
                di as text "  Example: treat_time(2010) if treatment began in 2010."
                exit 100
            }
        }
        else {
            di as error "Panel structure not set. Please specify treat_time(#)."
            exit 100
        }
    }
    // --- Method 3: Repeated cross-sections require explicit specification ---
    else {
        di as error "Treatment time cannot be inferred from the data."
        di as error "  Please specify treat_time(#), where # is the first"
        di as error "  post-treatment period."
        di as error "  Example: treat_time(2010) if treatment began in 2010."
        exit 100
    }
    
    // -------------------------------------------------------------------------
    // Step 4: Compute derived time parameters
    // -------------------------------------------------------------------------
    // Per Section 2.1:
    //   T_pre  = t0 - 1  (number of pre-treatment periods)
    //   T_post = T - T_pre (number of post-treatment periods)
    
    local T_pre = `t0' - 1
    local T_post = `T' - `T_pre'
    
    // -------------------------------------------------------------------------
    // Step 5: Validate T_pre >= 2
    // -------------------------------------------------------------------------
    // Per Section 2.1: "We presume that the pre-treatment period is of size
    // at least T_pre >= 2, so that at least one pre-treatment violation of
    // parallel trends is well-defined."
    //
    // This is required because the iterative violation nu_t is defined for
    // t = 2, ..., t0-1, which requires at least two pre-treatment periods.
    
    if `T_pre' < 2 {
        di as error "Insufficient pre-treatment periods."
        di as error "  The pretest methodology requires at least 2 pre-treatment"
        di as error "  periods (T_pre >= 2) for the iterative violations to be"
        di as error "  well-defined."
        di as error "  Current: T=`T', t0=`t0', T_pre=`T_pre'."
        di as error "  Reference: Mikhaeil & Harshaw (2025), Section 2.1."
        exit 104
    }
    
    // -------------------------------------------------------------------------
    // Step 6: Return results
    // -------------------------------------------------------------------------
    
    // Retrieve original t0 value for user-facing output
    local t0_orig_val : word `t0' of `time_levels'
    
    // Scalar returns
    return scalar T = `T'           // Total number of periods
    return scalar t0 = `t0'         // Treatment time (consecutive index)
    return scalar t0_orig = `t0_orig_val'  // Treatment time (original units)
    return scalar T_pre = `T_pre'   // Number of pre-treatment periods
    return scalar T_post = `T_post' // Number of post-treatment periods
    return scalar has_gap = `has_gap'      // Indicator for non-consecutive time
    return scalar t_min = `t_min'   // Minimum observed time value
    return scalar t_max = `t_max'   // Maximum observed time value
    
    // Local returns
    return local time_levels "`time_levels'"  // All unique time values
    return local t0_method "`t0_method'"      // Method used to determine t0
    
    // Matrix return: mapping from consecutive index to original time values
    return matrix time_vals = `time_vals'
    
end
