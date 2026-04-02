*! _pretest_detect_panel.ado
*! Panel structure detection for difference-in-differences analysis
*!
*! Part of the pretest package for valid inference when testing violations
*! of parallel trends in difference-in-differences designs.
*!
*! Description:
*!   Detects whether the current dataset has panel structure via xtset.
*!   The pretest methodology applies to both:
*!     - Panel data: n subjects tracked over T time periods
*!     - Repeated cross-sections: independent samples at each time period
*!   This module identifies which structure is present to inform subsequent
*!   estimation procedures.
*!
*! Syntax:
*!   _pretest_detect_panel
*!
*! Returns:
*!   r(is_panel)    - 1 if panel data, 0 if repeated cross-section
*!   r(panelvar)    - Panel identifier variable (if is_panel=1)
*!   r(timevar)     - Time variable (if is_panel=1)
*!   r(n_units)     - Number of unique panel units (if is_panel=1)
*!   r(n_periods)   - Number of distinct time periods (if is_panel=1)
*!   r(is_balanced) - 1 if balanced panel, 0 otherwise (if is_panel=1)
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470
*!   See Section 2.1 for data structure requirements.

program define _pretest_detect_panel, rclass
    version 17.0
    
    // Query xtset configuration
    capture xtset
    local xtset_rc = _rc
    
    if `xtset_rc' == 0 {
        // Panel structure detected via xtset
        local panelvar = r(panelvar)
        local timevar = r(timevar)
        
        // Validate panel variable is properly specified
        if "`panelvar'" != "" & "`panelvar'" != "." {
            return scalar is_panel = 1
            return local panelvar "`panelvar'"
            return local timevar "`timevar'"
            
            // Count unique panel units (n in the notation of the paper)
            tempvar unit_tag
            qui egen `unit_tag' = tag(`panelvar')
            qui count if `unit_tag' == 1
            local n_units = r(N)
            return scalar n_units = `n_units'
            
            // Determine panel balance and period count
            // Use levelsof to handle potential gaps in time variable
            qui levelsof `timevar', local(time_levels)
            local n_periods : word count `time_levels'
            qui count
            local n_obs = r(N)
            local is_balanced = (`n_obs' == `n_units' * `n_periods')
            return scalar is_balanced = `is_balanced'
            return scalar n_periods = `n_periods'
        }
        else {
            // xtset configured without panel variable (time-series only)
            // Treat as repeated cross-section for DID purposes
            return scalar is_panel = 0
            return local panelvar ""
            return local timevar ""
        }
    }
    else {
        // No xtset configuration: repeated cross-sectional data
        return scalar is_panel = 0
        return local panelvar ""
        return local timevar ""
    }
end
