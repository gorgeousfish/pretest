*! _pretest_parse_time v0.1.0 - Time Structure Parsing Module
*!
*! Description:
*!   Parses time variable structure and computes key time parameters.
*!   - T (total periods): Number of unique time values in data
*!   - t0 (treatment time): First post-treatment period
*!   - T_pre (pre-treatment periods): t0 - 1, periods t ∈ {1, ..., t0-1}
*!   - T_post (post-treatment periods): T - T_pre, periods t ∈ {t0, ..., T}
*!
*! Syntax:
*!   _pretest_parse_time timevar [, treatment(varname) treat_time(#) is_panel(#)]
*!
*! Options:
*!   treatment()  - Treatment variable (for inferring t0 from time-varying treatment)
*!   treat_time() - User-specified treatment time (highest priority)
*!   is_panel()   - Panel data indicator (0/1)
*!
*! Error codes:
*!   100 - Required parameter missing (treat_time not specified for repeated cross-section)
*!   104 - T_pre < 2 (at least 2 pre-treatment periods required)
*!   109 - treat_time outside data time range
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), Section 2.1
*!   - T periods indexed t = 1, ..., T
*!   - t0 is the time at which treatment was administered
*!   - T_pre = t0 - 1 is the number of pre-treatment periods
*!   - T_post = T - T_pre is the number of post-treatment periods
*!   - Requires T_pre >= 2 (iterative violation nu_t only defined for t >= 2)

program define _pretest_parse_time, rclass
    version 17.0
    
    // Parse syntax
    syntax varlist(min=1 max=1 numeric) ///
        [, TREATment(varname numeric) ///
           TREat_time(real -999) ///
           IS_panel(integer 0)]
    
    local timevar `varlist'
    
    // ========================================
    // 1. Get unique time values and compute T
    // ========================================
    
    // Get all unique time values (sorted ascending)
    qui levelsof `timevar', local(time_levels)
    local T : word count `time_levels'
    
    // Store time values in Stata matrix (for later access)
    tempname time_vals
    matrix `time_vals' = J(`T', 1, .)
    local i = 1
    foreach t of local time_levels {
        matrix `time_vals'[`i', 1] = `t'
        local ++i
    }
    
    // Get time range
    local t_min : word 1 of `time_levels'
    local t_max : word `T' of `time_levels'
    
    // ========================================
    // 2. Check for gaps in time variable
    // ========================================
    
    local has_gap = 0
    if `t_max' - `t_min' + 1 != `T' {
        local has_gap = 1
        di as text "Note: Time variable has non-consecutive values (`T' periods from `t_min' to `t_max')."
        di as text "      Internal mapping applied automatically."
    }
    
    // ========================================
    // 3. Determine treatment time t0
    // ========================================
    // Priority: user-specified > data inference > error
    
    local t0 = .
    local t0_method = ""
    
    // Mode 1: User explicitly specifies treat_time
    if `treat_time' != -999 {
        local t0_orig = `treat_time'
        
        // Validate treat_time is in data range
        // Note: treat_time is original value, need to check if exists in data
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
        
        // If treat_time not in data, check if outside range
        if `t0_found' == 0 {
            // treat_time is not a valid time value in data
            if `t0_orig' < `t_min' | `t0_orig' > `t_max' {
                di as error "Error 109: treat_time(`t0_orig') is outside data time range [`t_min', `t_max']"
                exit 109
            }
            else {
                // treat_time in range but not a valid time point (gap case)
                di as error "Error 109: treat_time(`t0_orig') is not a valid time value in the data"
                di as error "  Available time values: `time_levels'"
                exit 109
            }
        }
        
        // Use consecutive index as internal t0
        local t0 = `t0_index'
        local t0_method "user-specified"
    }
    // Mode 2: Infer from panel data with time-varying treatment
    else if `is_panel' == 1 & "`treatment'" != "" {
        // Check if treatment varies within individuals
        capture xtset
        if _rc == 0 {
            local panelvar = r(panelvar)
            
            // Compute treatment standard deviation within each individual
            tempvar treat_sd
            qui bysort `panelvar': egen `treat_sd' = sd(`treatment')
            
            // Check if time-varying treatment exists
            qui sum `treat_sd'
            if r(max) > 0 {
                // Treatment varies within some individuals, can infer t0
                // Find minimum time when treatment first becomes 1
                qui levelsof `timevar' if `treatment' == 1, local(treat_times)
                
                if "`treat_times'" != "" {
                    local t0_orig : word 1 of `treat_times'
                    
                    // Find corresponding consecutive index
                    local i = 1
                    foreach t of local time_levels {
                        if `t' == `t0_orig' {
                            local t0 = `i'
                        }
                        local ++i
                    }
                    local t0_method "inferred from time-varying treatment"
                    di as text "Note: t0=`t0_orig' (index `t0') inferred from time-varying treatment"
                }
                else {
                    // No observations with treatment==1
                    di as error "Error 100: No observations with treatment=1 found"
                    di as error "  Cannot infer treatment time. Please specify treat_time(#)"
                    exit 100
                }
            }
            else {
                // Treatment is time-invariant - typical for block-adoption design
                // Use friendly message instead of error
                di as text ""
                di as text "{bf:Note: Block-adoption design detected}"
                di as text "  Treatment indicator is time-invariant within individuals."
                di as text "  This is typical for block-adoption DID designs where:"
                di as text "    - Treatment group (D=1) is always treated after t₀"
                di as text "    - Control group (D=0) is never treated"
                di as text ""
                di as error "  → Please specify {bf:treat_time(#)} to indicate t₀"
                di as text "    t₀ = first post-treatment period (when treatment takes effect)"
                di as text "    Example: treat_time(2010) if policy was implemented in 2010"
                exit 100
            }
        }
        else {
            // xtset failed (should not happen since is_panel=1)
            di as error "Error 100: Panel structure not set. Please specify treat_time(#)"
            exit 100
        }
    }
    // Mode 3: Repeated cross-section, must specify
    else {
        di as error "Error 100: For repeated cross-sections, treat_time(#) is required"
        di as error "  Hint: t0 is the first post-treatment period"
        di as error "  Example: treat_time(2010) if treatment occurred in 2010"
        exit 100
    }
    
    // ========================================
    // 4. Compute time parameters
    // ========================================
    
    // T_pre = t0 - 1 (using consecutive index)
    local T_pre = `t0' - 1
    
    // T_post = T - T_pre
    local T_post = `T' - `T_pre'
    
    // ========================================
    // 5. Validate T_pre >= 2
    // ========================================
    // Paper Section 2.1: "We presume that the pre-treatment period is of size at least T_pre >= 2"
    // Reason: Iterative violation nu_t only defined for t >= 2
    
    if `T_pre' < 2 {
        di as error "Error 104: At least 2 pre-treatment periods required (T_pre >= 2)"
        di as error "  Current configuration: T=`T', t0=`t0', T_pre=`T_pre'"
        di as error "  Paper reference: Section 2.1"
        di as error "  Reason: Iterative violation nu_t is only defined for t >= 2"
        exit 104
    }
    
    // ========================================
    // 6. Return results
    // ========================================
    
    // Get original t0 value (for display)
    local t0_orig_val : word `t0' of `time_levels'
    
    return scalar T = `T'
    return scalar t0 = `t0'
    return scalar t0_orig = `t0_orig_val'
    return scalar T_pre = `T_pre'
    return scalar T_post = `T_post'
    return scalar has_gap = `has_gap'
    return scalar t_min = `t_min'
    return scalar t_max = `t_max'
    
    return local time_levels "`time_levels'"
    return local t0_method "`t0_method'"
    
    // Return time mapping matrix
    return matrix time_vals = `time_vals'
end
